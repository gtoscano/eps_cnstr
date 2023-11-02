#include "nlp.hpp"
#include <sys/stat.h>
#include <cassert>
#include <iostream>
#include <vector>
#include <float.h>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <unistd.h>
#include <fmt/core.h>
#include <filesystem>
#include <iomanip>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

int nvars;
int ncons;

namespace rnd
{
    static auto dev = std::random_device();
    static auto gen = std::mt19937{dev()};
    static auto dist_0_1 = std::uniform_real_distribution<double>(0, 1);

    bool flip(const double p)
    {
        return (dist_0_1(gen) < p);
    }

    int n_to_m_int(const int n, const int m)
    {
        auto dist_n_m_int = std::uniform_int_distribution<int>(n, m);
        return (dist_n_m_int(gen));
    }

    double n_to_m_real(const double n, const double m)
    {
        auto dist_n_m_real = std::uniform_real_distribution<double>(n, m);
        return (dist_n_m_real(gen));
    }
};


struct parcel_t {
    std::string name;
    std::vector<int> location;
    double amount;
    std::vector<std::vector<int> > bmps;
}; // { name: "":string, location: [s: int,h: int,u: int], bmps: [[bmp:integer, ...,], [...], ...], amount: :double}

struct var_t {
    std::string name;
    std::vector<int> location;
    double amount;
    int bmp;
};

void to_json(json &j, const var_t &p) {
    j = json{{"name",  p.name},
             {"location", p.location},
             {"amount",   p.amount},
             {"bmp", p.bmp}};
}

void from_json(const json &j, var_t &p) {
    j.at("name").get_to(p.name);
    j.at("location").get_to(p.location);
    j.at("amount").get_to(p.amount);
    j.at("bmp").get_to(p.bmp);
}

// constructor 63 0.9 0 1 0
EPA_NLP::EPA_NLP(std::string prefix_file, double max_constr, int pollutant_id, double limit_alpha, int cnstr_eval) {
    this->prefix = prefix_file;
    this->limit_alpha = limit_alpha;
    this->cnstr_eval = cnstr_eval;
    pollutant_idx = pollutant_id;
    this->total_cost = 1.0;
    this->total_acres = 1.0;
    //read_global_files(prefix_file, pollutant_id);
    this->max_constr = max_constr * sum_load_valid[pollutant_idx];
    has_content = false;
    msu_cbpo_path = getEnvVar("MSU_CBPO_PATH", "/opt/opt4cast");
}

void EPA_NLP::update_reduction(double max_constr) {
    this->max_constr = (1.0 - max_constr) * sum_load_valid[pollutant_idx];
}

EPA_NLP::~EPA_NLP() {}

bool EPA_NLP::get_nlp_info(
        Index &n,
        Index &m,
        Index &nnz_jac_g,
        Index &nnz_h_lag,
        IndexStyleEnum &index_style
) {
    int nvariables;
    int nconstraints; 

    get_info(nvariables, nconstraints, fmt::format("{}/output/nsga3/{}/config/", msu_cbpo_path, prefix));
/*
    initial_x.resize(nvariables);
    for (auto &x: initial_x) {
        x = 0.01;//rnd::n_to_m_real(0.0, 1.0);
    }
    */
    n = nvariables;
    m = nconstraints;
    
    nvars = n;
    ncons = m;

    int tmp_limit_bmp_counter = 0;
    for (auto const&[key, val]: limit_bmps_dict)
        tmp_limit_bmp_counter += limit_vars_dict[key].size();

    nnz_jac_g = nvars * 2 + tmp_limit_bmp_counter;
    nnz_h_lag = (int) (nvars * (nvars + 1)) / 2.0;
    index_style = TNLP::C_STYLE;

    return true;
}



// EFICCIENCY BEGIN
void EPA_NLP::compute_ef_keys() {

    for (const auto& pair : efficiency_) {
        ef_keys_.push_back(pair.first);
    }

    // Sort the vector of keys
    std::sort(ef_keys_.begin(), ef_keys_.end());
}

void EPA_NLP::compute_ef_size() {

    int counter = 0;

    for (const auto &key: ef_keys_) {
        auto bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                ++counter;
            }
        }
    }

    ef_size_ = counter;
    ef_begin_ = 0;
    ef_end_ = counter;
}

void EPA_NLP::normalize_efficiency() {
    int counter = 0;
    for (const auto &key: ef_keys_) {
        auto bmp_groups =  efficiency_[key];
        std::vector<std::vector<double>> grps_tmp;
        for (const auto& bmp_group : bmp_groups) {
            std::vector<double> grp_tmp;
            //The representation has a slack variable to represent the unused space
            double sum = x[counter];
            ++counter;

            for (const auto& bmp : bmp_group) {
                grp_tmp.push_back(x[counter]);
                sum += x[counter];
                ++counter;
            }

            for (auto &bmp: grp_tmp) {
                bmp = bmp / sum;
            }
            grps_tmp.push_back(grp_tmp);
        }
        ef_x_[key] = grps_tmp;
    }
}

int EPA_NLP::compute_ef() {
    double total_cost = 0;

    std::vector<double> fx(2, 0.0);
    std::vector<double> g(2, 0.0);
    int nconstraints = 2;

    compute_eta();

    std::vector<double> pt_load(3, 0.0);
    for (const auto &key: ef_keys_) {
        auto bmp_groups =  efficiency_[key];
        std::vector<double> prod(3, 1);
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto s = out[0];
        auto u = out[2];
        auto state_id = lrseg_.at(s)[1];
        auto alpha = amount_[key];
        auto bmp_group_idx = 0;
        for (const auto &bmp_group: bmp_groups) {
            std::vector<double> sigma(3, 0.0);
            auto bmp_idx = 0;
            auto sigma_cnstr = 0.0;

            for (const auto &bmp: bmp_group) {
                auto pct = ef_x_[key][bmp_group_idx][bmp_idx];
                sigma_cnstr += pct;
                double amount = pct * alpha;
                auto bmp_cost_idx = fmt::format("{}_{}", state_id, bmp);
                double cost = amount * bmp_cost_[bmp_cost_idx];

                std::string s_tmp = fmt::format("{}_{}_{}", bmp, s, u);
                auto eta = eta_dict_[s_tmp];
                sigma[0] += eta[0] * pct;
                sigma[1] += eta[1] * pct;
                sigma[2] += eta[2] * pct;
                sigma_cnstr += pct;
                ++bmp_idx;
            }
            prod[0] *= (1.0 - sigma[0]);
            prod[1] *= (1.0 - sigma[1]);
            prod[2] *= (1.0 - sigma[2]);
            //g[nconstraints] = (1.0 - sigma_cnstr);

            ++bmp_group_idx;
            ++nconstraints;
        }
        if(!phi_dict_.contains(key)) {
            std::cout<<"nooooooo.... no key here\n"<<key<<"\n";
        }
        auto phi = phi_dict_[key];

        pt_load[0] += phi[0] * alpha * prod[0];
        pt_load[1] += phi[1] * alpha * prod[1];
        pt_load[2] += phi[2] * alpha * prod[2];

    }

    //g[0] =  (sum_load_invalid[0] + pt_load[0]) - 0.8*(sum_load_invalid[0] + sum_load_valid[0]) ;
    double load_ub = sum_load_invalid_[load_to_opt_] + sum_load_valid_[load_to_opt_];
    double load_obtained = sum_load_invalid_[load_to_opt_] + pt_load[load_to_opt_];

    double load_lb = 0.8 * load_ub; //20% reduction
    //g[0] = (load_obtained < load_lb)?(load_obtained-load_lb)/load_ub:0;

    double cost_ub = 1000000.0;
    double cost_steps = 1000000.0;
    g[1] = (total_cost > cost_ub)? (cost_ub-total_cost)/cost_steps:0;

    fx[0] = total_cost;
    fx[1] = sum_load_invalid_[load_to_opt_] + pt_load[load_to_opt_];
    fmt::print("Total Cost: {}. Total Load: {}.\n", fx[0], fx[1]);

    //write_files_csv2 (selected_bmps, emo_uuid, exec_uuid);
    //write_files_animal_csv2 (selected_bmps, emo_uuid, exec_uuid);

    //int nrows = write_and_send_parquet_file2(selected_bmps, emo_uuid, exec_uuid);
    return 0;
}

void EPA_NLP::compute_eta() {
    auto redis = sw::redis::Redis(REDIS_URL);

    for (const auto &[key, bmp_groups]: efficiency_) {
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto s = out[0];
        auto u = out[2];
        for (const auto &bmp_group: bmp_groups) {
            for (const auto &bmp: bmp_group) {
                std::string s_tmp = fmt::format("{}_{}_{}", bmp, s, u);
                std::vector <std::string> eta_tmp;
                redis.lrange(s_tmp, 0, -1, std::back_inserter(eta_tmp));

                if (!eta_tmp.empty()) {
                    eta_dict_[s_tmp] = {stof(eta_tmp[0]), stof(eta_tmp[1]), stof(eta_tmp[2])};
                } else {
                    eta_dict_[s_tmp] = {0.0, 0.0, 0.0};
                    //std::cout << "ETA not found"<<key<<"\n";
                }
            }
        }
    }
}

void EPA_NLP::load(const std::string& filename) {

    // Open the JSON file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        exit(-1);
        return;
    }

    // Parse the JSON file directly into a nlohmann::json object
    json json_obj = json::parse(file);
    std::vector<std::string> keys_to_check = {"efficiency", "amount", "bmp_cost", "animal_unit", "lrseg", "scenario_data_str", "u_u_group", "counties" };
    for (const auto& key : keys_to_check) {
        if (!json_obj.contains(key)) {
            std::cout << "The JSON object does not contain the key '" << key << "'\n";
            exit(-1);
        }
    }

    // Access the JSON data
    amount_ = json_obj["amount"].get<std::unordered_map<std::string, double>>();

    auto efficiency_tmp = json_obj["efficiency"].get<std::unordered_map<std::string, std::vector<std::string>>>();
    for (const auto& [key, value] : efficiency_tmp) {
        std::vector<std::string> bmp_group;
        for (const auto& bmp_load_src : value) {
            std::vector<std::string> bmp_split;
            misc_utilities::split_str(bmp_load_src, '_', bmp_split);
            // if bmp_split[0] is in valid_lc_bmps, then add it to the bmp_group
            if(std::ranges::find(valid_lc_bmps_, bmp_split[0]) != valid_lc_bmps_.end()) {
                bmp_group.push_back(bmp_load_src);
            }
        }
        if (bmp_group.size() > 0) {
            efficiency[key] = bmp_group;
        }

    }
    
    bmp_cost_ = json_obj["bmp_cost"].get<std::unordered_map<std::string, double>>();
    auto lrseg_tmp = json_obj["lrseg"].get<std::unordered_map<std::string, std::vector<int>>>();
    scenario_data_str_ = "empty"+json_obj["scenario_data_str"].get<std::string>();

    scenario_data_str_ = "empty_38_6611_256_6_8_59_1_6608_158_2_31_8_410";
    //fmt::print("scenario_data_str_: {}\n", scenario_data_str_);
    for (const auto& [key, vec] : lrseg_tmp) {
        if (vec.size() == 4) {
            lrseg_dict_[std::stoi(key)] = std::make_tuple(vec[0], vec[1], vec[2], vec[3]);
        } else {
            std::cerr << "The vector size is not 4.\n";
        }
    }

    auto u_u_group_tmp = json_obj["u_u_group"].get<std::unordered_map<std::string, int>>();
    for (const auto& [key, value] : u_u_group_tmp) {
        u_u_group_dict[std::stoi(key)] = value;
    }


    auto geography_county_tmp = json_obj["counties"].get<std::unordered_map<std::string, std::tuple<int, int, std::string, std::string, std::string>>>();
    for (const auto& [key, val] : geography_county_tmp){
        geography_county_[std::stoi(key)] = val;
    }

    sum_load_valid_ = json_obj["sum_load_valid"].get<std::vector<double>>();
    sum_load_invalid_ = json_obj["sum_load_invalid"].get<std::vector<double>>();

}

// EFICCIENCY END 

bool EPA_NLP::get_n_and_m(
        Index &n,
        Index &m
) {
    int nconstraints = 1;
    int nvariables = 0;

    for (const auto &key: ef_keys_) {
        auto bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
                ++nconstraints;
            for (const auto &bmp : bmp_group) {
                ++nvariables; 
            }
        }
    }

    nconstraints += 1;
    std::cout << "Number of parcels " << ef_keys.size() << std::endl;
    std::cout << "Number of variables " << nvariables << std::endl;
    std::cout << "Number of constraints " << nconstraints << std::endl;
    nvars = nvariables;
    n = nvars;

    m = nconstraints;
    ncons = m;
    return true;
}

bool EPA_NLP::random_x(
        Index &n,
        Index &m,
        Number *x
) {
    int nvariables = 0;

    std::vector<double> my_initial_x;
    for (int i = 0; i < n; ++i) {
        x[i] = 0.0;
    }

    for (const auto &key: ef_keys_) {
        auto bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                x[nvariables] = rnd::n_to_m_real(0.0, 1.0);
                ++nvariables;
            }
        }
    }


    initial_x.resize(my_initial_x.size());
	std::copy(my_initial_x.begin(), my_initial_x.end(), initial_x.begin());

    return true;
}


bool EPA_NLP::get_bounds_info(
        Index n,
        Number *x_l,
        Number *x_u,
        Index m,
        Number *g_l,
        Number *g_u
) {
    for (Index i = 0; i < nvars; i++) {
        x_l[i] = 0.0;
        x_u[i] = 1.0; 
    }

    g_l[0] = 0.0;
    g_u[0] = this->max_constr;
    for (Index i = 1; i < m; ++i) {
        //for (Index i = 1; i < m; ++i) {
        g_l[i] = 0.0;
        g_u[i] = 1.0;
    }
    /*
    int cnstr_counter = m - limit_bmps_dict.size();
    for (auto const&[key, val]: limit_bmps_dict) {
        g_l[cnstr_counter] = val[0];

        g_u[cnstr_counter] = val[1];
        std::cout << val[0] << " " << val[1] << std::endl;
        ++cnstr_counter;
    }
    */

    return true;
}

bool EPA_NLP::my_get_starting_point(
        Number *x
) {
    // initialize to the given starting point
    for (int i = 0; i < nvars; ++i) {
        x[i] = initial_x[i];
    }

    return true;
}

bool EPA_NLP::set_starting_point(
        std::vector<double> my_x
) {
    // initialize to the given starting point
    initial_x.resize(my_x.size());
    std::copy(my_x.begin(), my_x.end(), initial_x.begin());

    return true;
}
bool EPA_NLP::get_starting_point(
        Index n,
        bool init_x,
        Number *x,
        bool init_z,
        Number *z_L,
        Number *z_U,
        Index m,
        bool init_lambda,
        Number *lambda
) {
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    // initialize to the given starting point
    for (int i = 0; i < nvars; ++i) {
        x[i] = initial_x[i];
        //      x[i] = rnd::n_to_m_real(0.0, 1.0);
    }

    return true;
}

bool EPA_NLP::eval_f(
        Index n,
        const Number *x,
        bool new_x,
        Number &obj_value
) {
    assert(n == nvars);
    double fitness = 0.0;

    int idx = 0;

    for (const auto &key: ef_keys_) {

        auto bmp_groups =  efficiency_[key];
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto load_src = out[2];
        auto state_id = lrseg_.at(lrseg)[1];
        auto alpha = amount_[key];
        auto bmp_group_idx = 0;
        auto bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                double amount = pct * alpha;
                auto bmp_cost_idx = fmt::format("{}_{}", state_id, bmp);
                double cost = x[idx] * amount * bmp_cost_[bmp_cost_idx];
                fitness += cost;
                ++idx;
            }
        }
    }


    obj_value = fitness;
    //assert(obj_value>0);
    return true;
}

bool EPA_NLP::eval_grad_f(
        Index n,
        const Number *x,
        bool new_x,
        Number *grad_f
) {
    int idx = 0;

    for (const auto &key: ef_keys_) {

        auto bmp_groups =  efficiency_[key];
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto load_src = out[2];
        auto state_id = lrseg_.at(lrseg)[1];
        auto alpha = amount_[key];
        auto bmp_group_idx = 0;
        auto bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                double amount = pct * alpha;
                auto bmp_cost_idx = fmt::format("{}_{}", state_id, bmp);
                grad_f[idx] =  amount * bmp_cost_[bmp_cost_idx];
                fitness += cost;
                ++idx;
            }
        }
    }


    return true;
}

bool EPA_NLP::my_eval_g(
        int n,
        const Number *x,
        std::vector<double> &g
) {

    assert(n == nvars);
        
    double total_cost = 0;
    std::vector<double> fx(2, 0.0);
    std::vector<double> g(2, 0.0);
    int nconstraints = 2;

    std::vector<double> pt_load(3, 0.0);
    for (const auto &key: ef_keys_) {
        auto bmp_groups =  efficiency_[key];
        std::vector<double> prod(3, 1);
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto load_src = out[2];
        auto state_id = lrseg_.at(lrseg)[1];
        auto alpha = amount_[key];
        auto bmp_group_idx = 0;
        for (const auto &bmp_group: bmp_groups) {
            std::vector<double> sigma(3, 0.0);
            auto bmp_idx = 0;
            auto sigma_cnstr = 0.0;

            for (const auto &bmp: bmp_group) {
                auto pct = x[idx];
                sigma_cnstr += pct;
                double amount = pct * alpha;

                std::string s_tmp = fmt::format("{}_{}_{}", bmp, s, u);
                auto eta = eta_dict_[s_tmp];
                sigma[0] += eta[0] * pct;
                sigma[1] += eta[1] * pct;
                sigma[2] += eta[2] * pct;
                sigma_cnstr += pct;
                ++bmp_idx;
            }
            prod[0] *= (1.0 - sigma[0]);
            prod[1] *= (1.0 - sigma[1]);
            prod[2] *= (1.0 - sigma[2]);
            //g[nconstraints] = (1.0 - sigma_cnstr);

            ++bmp_group_idx;
            ++nconstraints;
        }
        if(!phi_dict_.contains(key)) {
            std::cout<<"nooooooo.... no key here\n"<<key<<"\n";
        }
        auto phi = phi_dict_[key];

        pt_load[0] += phi[0] * alpha * prod[0];
        pt_load[1] += phi[1] * alpha * prod[1];
        pt_load[2] += phi[2] * alpha * prod[2];

    }

    //g[0] =  (sum_load_invalid[0] + pt_load[0]) - 0.8*(sum_load_invalid[0] + sum_load_valid[0]) ;
    double load_ub = sum_load_invalid_[load_to_opt_] + sum_load_valid_[load_to_opt_];
    double load_obtained = sum_load_invalid_[load_to_opt_] + pt_load[load_to_opt_];

    double load_lb = 0.8 * load_ub; //20% reduction
    //g[0] = (load_obtained < load_lb)?(load_obtained-load_lb)/load_ub:0;

    double cost_ub = 1000000.0;
    double cost_steps = 1000000.0;
    g[1] = (total_cost > cost_ub)? (cost_ub-total_cost)/cost_steps:0;

    fx[0] = total_cost;
    fx[1] = sum_load_invalid_[load_to_opt_] + pt_load[load_to_opt_];
    fmt::print("Total Cost: {}. Total Load: {}.\n", fx[0], fx[1]);
    return true;
}


bool EPA_NLP::eval_g(
        Index n,
        const Number *x,
        bool new_x,
        Index m,
        Number *g
) {
    assert(n == nvars);
    int parcel_idx = 0;
    int nconstraints = 1;
    double post_treatment_load = 0.0;

    std::unordered_map<int, double> sum_limit;
    double total_cost = 0.0;
    double total_acres = 0.0;
    for (auto const&[key, val]: limit_bmps_dict) {
        sum_limit[key] = 0.0;
    }

    for (int j(0); j < n;) {
        double alpha = alpha_vec[parcel_idx];
        double phi = phi_vec[parcel_idx][pollutant_idx];
        //[0] = s, [1] = h, [2] = u
        int s = s_h_u_vec[parcel_idx][0];
        int u = s_h_u_vec[parcel_idx][2];
        int acc_idx = 0;
        auto bmps = bmp_grp_src_links_vec[parcel_idx];
        double prod = 1.0;
        for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
            double sigma = 0.0;
            double sigma_cnstr = 0.0;
            for (int bc = 0; bc < bmp_counter; ++bc) {
                int bmp_idx = bmps[acc_idx + bc];

                double pct = x[j];
                //total_cost += x[j] * alpha * tau_dict[bmp_idx];
                //total_acres += x[j] * alpha;
                ++j;
                //b_s_u
                std::string s_tmp = fmt::format("{}_{}_{}", bmp_idx, s, u);

                double eta_value = 0.0;
                if (eta_dict.find(s_tmp) != eta_dict.end()) {
                    eta_value = eta_dict[s_tmp][pollutant_idx];
                }
                if (sum_limit.find(bmp_idx) != sum_limit.end()) {
                    sum_limit[bmp_idx] = sum_limit[bmp_idx] + alpha * pct;
                }

                sigma += eta_value * pct;
                sigma_cnstr += pct;
            }
            g[nconstraints] = sigma_cnstr;
            ++nconstraints;
            acc_idx += bmp_counter;
            prod *= (1.0 - sigma);
        }
        post_treatment_load += phi * alpha * prod;
        ++parcel_idx;
    }
    g[0] = post_treatment_load;
    for (auto const&[key, val]: limit_bmps_dict) {
        auto vars = limit_vars_dict[key];
        auto alphas = limit_alpha_dict[key];
        double tmp_sum = 0.0;
        for (int i(0); i < vars.size(); ++i) {
            tmp_sum += x[vars[i]] * alphas[i];
        }
        g[nconstraints] = tmp_sum;
        ++nconstraints;
    }
    /*
    g[nconstraints] = 0.5;//total_cost;
    g[nconstraints+1] = 0.5; //total_acres;

    */
    return true;
}

// return the structure or values of the Jacobian
bool EPA_NLP::eval_jac_g(
        Index n,
        const Number *x,
        bool new_x,
        Index m,
        Index nele_jac,
        Index *iRow,
        Index *jCol,
        Number *values
) {
    assert(n == nvars);
    if (values == NULL) {
        // return the structure of the Jacobian

        // this particular Jacobian is not dense
        for (int i = 0; i < nvars; ++i) {
            iRow[i] = 0;
            jCol[i] = i;
        }
        int jac_index = nvars;
        int jac_row = 1;
        int nparcels = bmp_grp_src_counter_vec.size();

        for (int parcel_idx(0); parcel_idx < nparcels; ++parcel_idx) {
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    iRow[jac_index] = jac_row;
                    jCol[jac_index] = jac_index - nvars;
                    ++jac_index;
                }
                ++jac_row;
            }
        }

        for (auto const&[key, val]: limit_bmps_dict) {
            for (auto &var_idx: limit_vars_dict[key]) {
                iRow[jac_index] = jac_row;
                jCol[jac_index] = var_idx;
                ++jac_index;
            }
            ++jac_row;

        }
        /*
        for (int i = 0; i < nvars; ++i) {
            iRow[jac_index] = jac_row;
            jCol[jac_index] = i;
            ++jac_index;
        }
        ++jac_row;
        for (int i = 0; i < nvars; ++i) {
            iRow[jac_index] = jac_row;
            jCol[jac_index] = i;
            ++jac_index;
        }
        */

    } else {
        // return the values of the Jacobian of the constraints
        int parcel_idx = 0;
        for (int j(0); j < nvars;) {
            double alpha = alpha_vec[parcel_idx];
            double phi = phi_vec[parcel_idx][pollutant_idx];
            int s = s_h_u_vec[parcel_idx][0];
            int u = s_h_u_vec[parcel_idx][2];
            int acc_idx = 0;
            auto bmps = bmp_grp_src_links_vec[parcel_idx];
            double prod = 1.0;
            int ngroups = (int) bmp_grp_src_counter_vec[parcel_idx].size();
            std::vector<double> sigma_less(bmps.size(), 0.0);
            std::vector<double> sigma_full(ngroups, 0.0);
            int j_saved = j;
            int grp_counter = 0;
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
                double sigma = 0.0;
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    int bmp_idx = bmps[acc_idx + bc];
                    double pct = x[j];
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp_idx, s, u);

                    if (eta_dict.find(s_tmp) != eta_dict.end()) {
                        sigma += eta_dict[s_tmp][pollutant_idx] * pct;
                    }
                    ++j;
                }

                sigma_full[grp_counter] = sigma;
                j -= bmp_counter;

                for (int bc = 0; bc < bmp_counter; ++bc) {
                    int bmp_idx = bmps[acc_idx + bc];
                    double pct = x[j];
                    ++j;
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp_idx, s, u);
                    double eta_value = 0.0;
                    if (eta_dict.find(s_tmp) != eta_dict.end()) {
                        eta_value = eta_dict[s_tmp][pollutant_idx];
                    }

                    sigma_less[acc_idx + bc] = sigma - eta_value * pct;
                    //values[nvars+cnstr_counter] = sigma_less[acc_idx + bc];
                }
                acc_idx += bmp_counter;
                prod *= 1.0 - sigma;
                ++grp_counter;
            }
            grp_counter = 0;
            j = j_saved;
            acc_idx = 0;
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    int bmp_idx = bmps[acc_idx + bc];
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp_idx, s, u);
                    double eta_value = 0.0;
                    if (eta_dict.find(s_tmp) != eta_dict.end()) {
                        eta_value = eta_dict[s_tmp][pollutant_idx];
                    }
                    values[j] = prod / (1.0 - sigma_full[grp_counter]);
                    values[j] *= (1.0 - sigma_less[acc_idx + bc]);
                    values[j] *= alpha * phi * (-eta_value);
                    ++j;
                }
                acc_idx += bmp_counter;
                ++grp_counter;
            }

            ++parcel_idx;
        }
        int nparcels = bmp_grp_src_counter_vec.size();
        int var_counter = 0;
        for (int parcel_idx(0); parcel_idx < nparcels; ++parcel_idx) {
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    values[nvars + var_counter] = 1;
                    ++var_counter;
                }
            }
        }
        int limit_cnstr_counter = nvars * 2;

        for (auto const&[key, val]: limit_bmps_dict) {
            for (auto &my_alpha: limit_alpha_dict[key]) {
                values[limit_cnstr_counter] = my_alpha;
                ++limit_cnstr_counter;
            }
        }


    }
    return true;
}
//return the structure or values of the Hessian
bool EPA_NLP::eval_h(
        Index n,
        const Number *x,
        bool new_x,
        Number obj_factor,
        Index m,
        const Number *lambda,
        bool new_lambda,
        Index nele_hess,
        Index *iRow,
        Index *jCol,
        Number *values
) {
    assert(n == nvars);

    if (values == NULL) {
        Index idx = 0;
        for (Index row = 0; row < nvars; row++) {
            for (Index col = 0; col <= row; col++) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }

        assert(idx == nele_hess);
    } else {
        int parcel_idx = 0;
        for (int j(0); j < nvars;) {
            double alpha = alpha_vec[parcel_idx];
            double phi = phi_vec[parcel_idx][pollutant_idx];
            int s = s_h_u_vec[parcel_idx][0];
            int u = s_h_u_vec[parcel_idx][2];
            int acc_idx = 0;
            auto bmps = bmp_grp_src_links_vec[parcel_idx];
            double prod = 1.0;
            int ngroups = (int) bmp_grp_src_counter_vec[parcel_idx].size();
            std::vector<double> sigma_less(bmps.size(), 0.0);
            std::vector<double> sigma_full(ngroups, 0.0);
            int j_saved = j;
            int grp_counter = 0;
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
                double sigma = 0.0;
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    int bmp_idx = bmps[acc_idx + bc];
                    double pct = x[j];
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp_idx, s, u);
                    sigma += eta_dict[s_tmp][pollutant_idx] * pct;
                    ++j;
                }

                sigma_full[grp_counter] = sigma;
                j -= bmp_counter;

                for (int bc = 0; bc < bmp_counter; ++bc) {
                    int bmp_idx = bmps[acc_idx + bc];
                    double pct = x[j];
                    ++j;
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp_idx, s, u);
                    double eta_value = 0.0;
                    if (eta_dict.find(s_tmp) != eta_dict.end()) {
                        eta_value = eta_dict[s_tmp][pollutant_idx];
                    }
                    sigma_less[acc_idx + bc] = sigma - eta_value * pct;
                    //values[nvars+cnstr_counter] = sigma_less[acc_idx + bc];
                }
                acc_idx += bmp_counter;
                prod *= 1.0 - sigma;
                ++grp_counter;
            }
            grp_counter = 0;
            j = j_saved;
            acc_idx = 0;
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[parcel_idx]) {
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    int bmp_idx = bmps[acc_idx + bc];
                    ++j;
                }
                acc_idx += bmp_counter;
                ++grp_counter;
            }

            ++parcel_idx;
        }
        int nparcels = bmp_grp_src_counter_vec.size();
        int var_counter = 0;
        for (int my_parcel_idx(0); my_parcel_idx < nparcels; ++my_parcel_idx) {
            for (auto &&bmp_counter: bmp_grp_src_counter_vec[my_parcel_idx]) {
                for (int bc = 0; bc < bmp_counter; ++bc) {
                    ++var_counter;
                }
            }
        }
    }

    return false;
}


void EPA_NLP::save_files(
        Index                      n,
        const Number *x
) {
    //std::cout << "Cost(x*) = " << obj_value << std::endl;
    auto filename = fmt::format("{}/output/nsga3/{}/config/ipopt_{}.csv", msu_cbpo_path, prefix,0);
    auto json_filename = fmt::format("{}/output/nsga3/{}/config/ipopt.json", msu_cbpo_path, prefix);
    auto directory = fmt::format("{}/output/nsga3/{}/config", msu_cbpo_path, prefix);
    if (fs::exists(directory) == false) {
        fs::create_directories(directory);
    }

    json json_ipopt = {};
    if (fs::exists(json_filename) == true) {
        std::ifstream ifs_json(json_filename);
        json_ipopt = json::parse(ifs_json);
    }

    json this_json_ipopt = {};
    json empty_json = {};

    std::ofstream ofile(filename);
    ofile.precision(10);
    var_t v;

    final_x.resize(n);
    for(int i(0); i<n; ++i) {
       final_x[i] = x[i];
    }

    int idx = 0;
    int parcel_idx = 0;
    // Iterate through each line and split the content using delimeter
    bool flag = false;
    for (auto &&bmp_per_parcel: bmp_grp_src_links_vec) {
        int s = s_h_u_vec[parcel_idx][0];
        int h = s_h_u_vec[parcel_idx][1];
        int u = s_h_u_vec[parcel_idx][2];
        int unit_id = s_h_u_vec[parcel_idx][4];
        for (auto &&bmp_idx: bmp_per_parcel) {
            if(idx > n) {std::cerr<<"There is a problem n greater than maximum allowed\n";}
            if (x[idx] * alpha_vec[parcel_idx] > 1.0) {
                double amount = x[idx];
                ofile<<fmt::format("{},{},{},{},{},{}\n", s, h, u, bmp_idx, unit_id, amount);
                v.name = fmt::format("{}_{}_{}_{}", s, h, u, bmp_idx);
                v.location = {s,h,u};
                v.amount = (amount>1)?1.0:((amount<0)?0.0:amount);
                v.bmp = bmp_idx;
                this_json_ipopt.emplace_back(v);
                flag = true;
            }
            ++idx;
        }
        ++parcel_idx;
    }
    has_content = flag;
    if (flag) {
        json_ipopt.emplace_back(this_json_ipopt);
    }
    else {
        std::cout<<"Failed! No solution!\n";
    }
    std::ofstream ofs(json_filename);
    ofs<<json_ipopt<<std::endl;
    ofile.close();
}

void EPA_NLP::write_files(
        Index n,
        const Number *x,
        Index m,
        Number obj_value
) {

    std::vector<double> g_constr(3, 0.0);
    my_eval_g(nvars, x, g_constr);

    std::cout.precision(10);
    std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
    // std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
    std::cout << std::endl << std::endl << "Objective value" << std::endl;

    std::cout << "Cost(x*) = " << obj_value << std::endl;
    std::cout << std::endl << "Final value of the constraints:" << std::endl;


    fs::create_directories(fmt::format("{}/output/nsga3/{}/config/", msu_cbpo_path, prefix));
    std::string fileName = fmt::format("{}/output/nsga3/{}/config/ipopt_output.csv", msu_cbpo_path, prefix);
    std::ofstream file(fileName);
    file.precision(15);

    fileName = fmt::format("{}/output/nsga3/{}/config/ipopt_output.txt", msu_cbpo_path, prefix);
    std::ofstream file2(fileName);
    file2.precision(15);
    std::cout << "g_constr[NLoadEos] = " << g_constr[0] + sum_load_invalid[0] << "\n";
    std::cout << "g_constr[PLoadEos] = " << g_constr[1] + sum_load_invalid[1] << "\n";
    std::cout << "g_constr[SLoadEos] = " << g_constr[2] + sum_load_invalid[2] << "\n";
    std::cout << "Original NLoadEos = " << sum_load_valid[0] + sum_load_invalid[0] << "\n";
    std::cout << "Original PLoadEos = " << sum_load_valid[1] + sum_load_invalid[1] << "\n";
    std::cout << "Original SLoadEos = " << sum_load_valid[2] + sum_load_invalid[2] << "\n";

    file2 << obj_value << "," << g_constr[0] + sum_load_invalid[0] << "," << g_constr[1] + sum_load_invalid[1] << ","
          << g_constr[2] + sum_load_invalid[2] << "," << sum_load_valid[0] + sum_load_invalid[0] << ","
          << sum_load_valid[1] + sum_load_invalid[1] << "," << sum_load_valid[2] + sum_load_invalid[2] << ",";
    file2.close();
    // Iterate through each line and split the content using delimeter
    int idx = 0;
    int parcel_idx = 0;

    file<<"BmpSubmittedId,AgencyId,StateUniqueIdentifier,StateId,BmpId,GeographyId,LoadSourceId,UnitId,Amount,isValid,ErrorMessage,RowIndex,Cost,LrsegId,LoadSourceIdOriginal,TotalAmount,Capital\n";
    int counter = 0;
    int my_counter = 1;
    std::unordered_map<int, double> bmp_sum;

    for (auto &&bmp_per_parcel: bmp_grp_src_links_vec) {
        for (auto &&bmp_idx: bmp_per_parcel) {
            if (x[idx] > 0.001) {
                std::string tmp = "";
                int s = s_h_u_vec[parcel_idx][0];
                int h = s_h_u_vec[parcel_idx][1];
                int u = s_h_u_vec[parcel_idx][2];
                int unit_id = s_h_u_vec[parcel_idx][4];
                double amount = x[idx] * 100.0;
                double cost = 0.0;
                double capital = 0.0;

                unit_id = bmp_unit_dict[bmp_idx];
                if (unit_id == 1) { //acres
                    amount = x[idx] * alpha_vec[parcel_idx];
                    cost = amount * tau_dict[bmp_idx];
                    capital = amount * capital_tau_dict[bmp_idx];
                } else {
                    cost = x[idx] * alpha_vec[parcel_idx] * tau_dict[bmp_idx];
                    capital = x[idx] * alpha_vec[parcel_idx] * capital_tau_dict[bmp_idx];
                }

                file << counter + 1 << "," << h << ",SU" << counter << "," << s_state_dict[s] << "," << bmp_idx << ","
                     << s_geography_dict[s] << "," << u_u_group_dict[u] << "," << unit_id << "," << amount << ",True,,"
                     << counter + 1 << "," << cost << "," << s << "," << u << "," << alpha_vec[parcel_idx] << ","
                     << capital << "\n";
                ///std::cout << h << std::endl;
                //std::cout<<agencies_dict[h]<<std::endl;
                my_counter++;
                counter++;
                if (bmp_sum.find(bmp_idx) != bmp_sum.end()) {
                    bmp_sum[bmp_idx] += amount;
                } else {
                    bmp_sum[bmp_idx] = amount;
                }
                //file<<alpha_vec[idx_parcel]<< "," << bmp_vec[bmp][0]<<","<<x[idx] << ","<<alpha_vec[idx_parcel] *x[idx] << "\n";

            }
            ++idx;
        }
        ++parcel_idx;
    }

    fileName = "msu_cbpo_path + /output/" + prefix + "_bmp_summary.txt";

    std::ofstream file3(fileName);

    file3.precision(10);
    file3 << "Bmp_id,Acres,Cost,Capital" << std::endl;
    for (auto const&[key, sum]: bmp_sum) {
        file3 << key << "," << sum << "," << tau_dict[key] * sum << "," << capital_tau_dict[key] * sum << '\n';
    }
    file3.close();
    std::cout << "# Bmps: " << counter << std::endl;
    file.close();

    fs::create_directories(fmt::format("{}/output/nsga3/{}/config/", msu_cbpo_path, prefix));
    std::string filename_funcs= fmt::format("{}/output/nsga3/{}/config/ipopt_funcs.txt", msu_cbpo_path, prefix);
    std::string filename_vars = fmt::format("{}/output/nsga3/{}/config/ipopt_vars.txt", msu_cbpo_path, prefix);
    std::ofstream file_vars(filename_vars, std::ios_base::app);
    std::ofstream file_funcs(filename_funcs, std::ios_base::app);

    file_vars.precision(10);
    file_funcs.precision(10);

    for (int i = 0; i < n; i++){
	    double amount = (double) x[i];
	    if (amount < 0) amount = 0.0;
	    if (i+1 == n){
		file_vars<<amount<<"\n";
	    }
	    else{
	    	file_vars<< amount << ",";
	    }
    }

    file_vars.close();
    file_funcs<< obj_value << "," << g_constr[0] << "\n";// << g_constr[1] <<","<<g_constr[2] << "\n";
    file_funcs.close();
}

void EPA_NLP::write_formar_cast_file() {

    char filename[2000];
    std::unordered_map<int, std::string> usa_states;
    std::unordered_map<int, std::string> agencies;
    std::unordered_map<int, std::string> geography_name;
    std::unordered_map<int, std::string> load_source_group;

    std::string fileName = msu_cbpo_path + "/output/442_land_former_core_cast_output.txt";
    std::ofstream land_file(fileName);

    std::string path(msu_cbpo_path + "/output");
    std::string csvs_path(msu_cbpo_path + "/csvs");
    std::string county_name("442");

    sprintf(filename, "%s/csvs/%s", msu_cbpo_path, "usa_states.csv");
    CSVReader reader_usa_states(filename);

    sprintf(filename, "%s/csvs/%s", msu_cbpo_path, "agencies.csv");
    CSVReader reader_agencies(filename);

    sprintf(filename, "%s/csvs/%s", msu_cbpo_path, "geography_name.csv");
    CSVReader reader_geography_name(filename);

    sprintf(filename, "%s/csvs/%s", msu_cbpo_path, "load_source_group.csv");
    CSVReader reader_load_source_group(filename);
    sprintf(filename, "%s/output/%s", msu_cbpo_path, "442_output.csv");
    CSVReader reader_output(filename);

    std::vector<std::vector<std::string> > dataList = reader_usa_states.getData();

    for (std::vector<std::string> vec: dataList) {
        std::string tmp_str(vec[1]);
        //usa_states.insert(std::pair<int,std::string>(stoi(vec[0]), vec[1]));
        usa_states.insert({stoi(vec[0]), vec[1]});
        //usa_states[stoi(vec[0])] = tmp_str;
    }

    std::vector<std::vector<std::string> > dataList2 = reader_agencies.getData();

    for (std::vector<std::string> vec: dataList2) {
        std::string tmp_str(vec[1]);
        //agencies[stoi(vec[0])] = tmp_str;
        //agencies.insert(std::pair<int,std::string>(stoi(vec[0]), vec[1]));
        agencies.insert({stoi(vec[0]), vec[1]});
    }


    std::vector<std::vector<std::string> > dataList3 = reader_geography_name.getData();

    for (std::vector<std::string> vec: dataList3) {
        std::string tmp_str(vec[1]);
        //geography_name[stoi(vec[0])] = tmp_str;
        //geography_name.insert(std::pair<int,std::string>(stoi(vec[0]), vec[1]));
        geography_name.insert({stoi(vec[0]), vec[1]});
    }


    std::vector<std::vector<std::string> > dataList4 = reader_load_source_group.getData();

    for (std::vector<std::string> vec: dataList4) {
        std::string tmp_str(vec[1]);
        //load_source_group[stoi(vec[0])] = tmp_str;
        //load_source_group.insert(std::pair<int,std::string>(stoi(vec[0]), vec[1]));
        load_source_group.insert({stoi(vec[0]), vec[1]});
    }

    land_file<< "StateUniqueIdentifier\tAgencyCode\tStateAbbreviation\tBmpShortname\tGeographyName\tLoadSourceGroup\tAmount\tUnit\r\n";
    std::vector<std::vector<std::string> > ifile = reader_output.getData();
    FILE *fp;
    fp = fopen("/home/gtoscano/output/gtp.txt", "w");
    FILE *o_file = fopen("/home/gtoscano/output/gtp2.txt", "w");
    for (std::vector<std::string> row: ifile) {
        int counter = stoi(row[0]);
        int h = stoi(row[1]);
//land_file<< agencies[h].c_str() << " ";
        int StateId = stoi(row[3]);
        land_file << usa_states[StateId].c_str() << std::ends;
        fwrite((usa_states[StateId]).c_str(), 1, (usa_states[StateId]).size(), o_file);
        int BmpId = stoi(row[4]);
        //land_file << bmp_short_name_dict[BmpId].c_str()<< std::ends;
        //fwrite((bmp_short_name_dict[BmpId]).c_str(), 1, (bmp_short_name_dict[BmpId]).size(), o_file);
        int GeographyId = stoi(row[5]);
        land_file << geography_name[GeographyId].c_str() << std::ends;
        fwrite((geography_name[GeographyId]).c_str(), 1, (geography_name[GeographyId]).size(), o_file);
        int LoadSourceId = stoi(row[6]);
        land_file << load_source_group[LoadSourceId].c_str() << std::endl;
        int UnitId = stoi(row[7]);
        double Amount = stof(row[8]);
        //std::string RowIndex = row[11];
        char unit_type[100];

        if (UnitId == 0) {
            strcpy(unit_type, "percent");
        } else {
            strcpy(unit_type, "acres");
        }
        std::cout << agencies[h] << std::endl;
        //std::cout<<counter<<" "<<(agencies[h]).c_str()<<" "<<(usa_states[StateId]).c_str()<<" "<<(bmp_short_name_dict[BmpId]).c_str()<<" "<<(geography_name[GeographyId]).c_str()<<" "<<(load_source_group[LoadSourceId]).c_str()<<" "<<Amount<<" "<<unit_type<<std::endl;
        std::cout << counter << " " << (agencies[h]).c_str() << " " << (usa_states[StateId]).c_str() << " "
                  << (geography_name[GeographyId]).c_str() << " " << (load_source_group[LoadSourceId]).c_str() << " "
                  << Amount << " " << unit_type << std::endl;
        fprintf(fp, "%i\t%s\t%s\t%s\t%s\t%s\t%f\t%s\r\n", counter, (agencies[h]).c_str(), (usa_states[StateId]).c_str(),
                (bmp_short_name_dict[BmpId]).c_str(), (geography_name[GeographyId]).c_str(),
                (load_source_group[LoadSourceId]).c_str(), Amount, unit_type);
        printf("%i\t%s\t%s\t%s\t%s\t%s\t%f\t%s\r\n", counter, (agencies[h]).c_str(), (usa_states[StateId]).c_str(),
               (bmp_short_name_dict[BmpId]).c_str(), (geography_name[GeographyId]).c_str(),
               (load_source_group[LoadSourceId]).c_str(), Amount, unit_type);

    }
    fclose(fp);
    fclose(o_file);
    land_file.close();

}

// [TNLP_finalize_solution]
void EPA_NLP::finalize_solution(
        SolverReturn status,
        Index n,
        const Number *x,
        const Number *z_L,
        const Number *z_U,
        Index m,
        const Number *g,
        const Number *lambda,
        Number obj_value,
        const IpoptData *ip_data,
        IpoptCalculatedQuantities *ip_cq
) {
    final_x.resize(n);
    status_result= (int) status;
    std::cout<<"JAJAJAJA: "<<status<<std::endl;

    for (size_t i = 0; i < n; i++)
    {
        final_x[i] = x[i];
    }
    
    write_files(n, x, m, obj_value);
    save_files(n, x);
    // write_formar_cast_file();
}

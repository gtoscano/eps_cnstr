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
#include <random>
#include <nlohmann/json.hpp>

#include <boost/algorithm/string.hpp>
#include <misc_utilities.h>

using json = nlohmann::json;
namespace fs = std::filesystem;

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

namespace {
    std::string REDIS_HOST = misc_utilities::get_env_var("REDIS_HOST", "127.0.0.1");
    std::string REDIS_PORT = misc_utilities::get_env_var("REDIS_PORT", "6379");
    std::string REDIS_DB_OPT = misc_utilities::get_env_var("REDIS_DB_OPT", "1");
    std::string REDIS_URL = fmt::format("tcp://{}:{}/{}", REDIS_HOST, REDIS_PORT, REDIS_DB_OPT);

}
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
EPA_NLP::EPA_NLP(std::string filename_in, std::string filename_out, double max_constr, int pollutant_id) {
    this->filename_in = filename_in;
    load(filename_in);
    this->filename_out = filename_out;
    this-> pollutant_idx = pollutant_id;
    this->total_cost = 1.0;
    this->total_acres = 1.0;
    //read_global_files(prefix_file, pollutant_id);
    update_reduction(max_constr);
    fmt::print("Max constr: {}, {}\n", this->max_constr, max_constr);
    has_content = false;
}

void EPA_NLP::update_reduction(double max_constr) {
    this->max_constr = (1.0 - max_constr) * sum_load_valid_[pollutant_idx];
}

EPA_NLP::~EPA_NLP() {}

bool EPA_NLP::get_nlp_info(
        Index &n,
        Index &m,
        Index &nnz_jac_g,
        Index &nnz_h_lag,
        IndexStyleEnum &index_style
) {

    
    n = nvars_;
    m = ncons_;

    int tmp_limit_bmp_counter = 0;
    /*
    for (auto const&[key, val]: limit_bmps_dict)
        tmp_limit_bmp_counter += limit_vars_dict[key].size();
        */

    nnz_jac_g = nvars_ * 2 + tmp_limit_bmp_counter;
    nnz_h_lag = (int) (nvars_ * (nvars_ + 1)) / 2.0;
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

void EPA_NLP::filter_ef_keys() {
    auto redis = sw::redis::Redis(REDIS_URL);

    std::vector<std::string> keys_to_remove;
    size_t bmps_removed = 0;
    size_t bmps_sum = 0;
    size_t bmps_sum2 = 0;
    for (const auto &key: ef_keys_) {
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto load_src = out[2];
        auto state_id = lrseg_.at(lrseg)[1];
        auto& bmp_groups =  efficiency_[key];
        if (amount_.find(key) == amount_.end() ||
            phi_dict_.find(key) == phi_dict_.end() ||
            bmp_groups.size() <= 0 ) {
            keys_to_remove.push_back(key);
            continue;
        }
        auto alpha = amount_[key];
        auto phi = phi_dict_[key];
        
        if (alpha <= 1.0){
            keys_to_remove.push_back(key);
            continue;
        }
        for (auto &bmp_group : bmp_groups) {
            std::vector<int> bmps_to_remove;
            for (const auto &bmp : bmp_group) {
                std::string s_tmp = fmt::format("{}_{}_{}", bmp, lrseg, load_src);
                auto bmp_cost_key = fmt::format("{}_{}", state_id, bmp);
                if (!redis.hexists("ETA", s_tmp) ||
                    bmp_cost_.find(bmp_cost_key) == bmp_cost_.end() ||
                    bmp_cost_[bmp_cost_key] <= 0.0 ) {
                    //remove bmp from bmp_group
                    bmps_to_remove.push_back(bmp);
                }
            }
            bmps_removed += bmps_to_remove.size();
            bmps_sum += bmp_group.size();

            for (const auto &bmp : bmps_to_remove) {
                //fmt::print("Removing bmp {} from bmp_group\n", bmp);
                bmp_group.erase(std::remove(bmp_group.begin(), bmp_group.end(), bmp), bmp_group.end());
            }
            bmps_sum2 += bmp_group.size();
        }
    }
    // Sort keys_to_remove for efficient searching
    std::sort(keys_to_remove.begin(), keys_to_remove.end());
    
    // Use erase-remove idiom with a lambda function
    ef_keys_.erase(std::remove_if(ef_keys_.begin(), ef_keys_.end(), 
        [&keys_to_remove](const std::string& key) {
            return std::binary_search(keys_to_remove.begin(), keys_to_remove.end(), key);
        }), 
        ef_keys_.end());
    fmt::print("bmps_remove {} out from {} total {}\n", bmps_removed, bmps_sum, bmps_sum2);
    fmt::print("keys_to_remove size: {}\n", keys_to_remove.size());
}


void EPA_NLP::compute_ef_size() {
    nvars_ = 0;
    ncons_ = 1;

    int counter = 0;

    for (const auto &key: ef_keys_) {
        auto& bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                ++counter;
                ++nvars_;
            }
            ++ncons_;
        }
    }

    ef_size_ = counter;
}


void EPA_NLP::compute_eta() {
    auto redis = sw::redis::Redis(REDIS_URL);

    for (const auto &key: ef_keys_) {
        auto& bmp_groups =  efficiency_[key];
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto s = out[0];
        auto u = out[2];
        for (const auto &bmp_group: bmp_groups) {
            for (const auto &bmp: bmp_group) {
                std::string s_tmp = fmt::format("{}_{}_{}", bmp, s, u);

	        	std::vector<std::string> eta_tmp;
                std::string eta_str = "0.0_0.0_0.0";
                if (redis.hexists("ETA", s_tmp)) {
                    eta_str = *redis.hget("ETA", s_tmp);
                }
                else {
                    std::cout << "No ETA for " << s_tmp << std::endl;
                }

                boost::split(eta_tmp, eta_str, boost::is_any_of("_"));

                if(!eta_tmp.empty()) {
                   std::vector<double> content_eta({stof(eta_tmp[0]), stof(eta_tmp[1]), stof(eta_tmp[2])});
                   eta_dict_[s_tmp] = content_eta;
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
    std::vector<std::string> keys_to_check = {"amount", "phi", "efficiency", "lrseg", "bmp_cost", "u_u_group", "sum_load_valid", "sum_load_invalid" };
    for (const auto& key : keys_to_check) {
        if (!json_obj.contains(key)) {
            std::cout << "The JSON object does not contain the key '" << key << "'\n";
            exit(-1);
        }
    }

    // Access the JSON data
    amount_ = json_obj["amount"].get<std::unordered_map<std::string, double>>();
    phi_dict_ = json_obj["phi"].get<std::unordered_map<std::string, std::vector<double>>>();
    efficiency_ = json_obj["efficiency"].get<std::unordered_map<std::string, std::vector<std::vector<int>>>>();
    /*
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
            efficiency_[key] = bmp_group;
        }

    }
    */
    
    bmp_cost_ = json_obj["bmp_cost"].get<std::unordered_map<std::string, double>>();
    lrseg_ = json_obj["lrseg"].get<std::unordered_map<std::string, std::vector<int>>>();


    u_u_group_dict = json_obj["u_u_group"].get<std::unordered_map<std::string, int>>();


    sum_load_valid_ = json_obj["sum_load_valid"].get<std::vector<double>>();
    sum_load_invalid_ = json_obj["sum_load_invalid"].get<std::vector<double>>();
    compute_ef_keys();
    filter_ef_keys();
    compute_ef_size();
    compute_eta();
}

// EFICCIENCY END 

bool EPA_NLP::get_n_and_m(
        Index &n,
        Index &m
) {
    n = nvars_;
    m = ncons_;

    std::cout << "Number of parcels " << ef_keys_.size() << std::endl;
    std::cout << "Number of variables " << nvars_ << std::endl;
    std::cout << "Number of constraints " << ncons_ << std::endl;
    return true;
}

bool EPA_NLP::random_x(
        Index &n,
        Index &m,
        Number *x
) {
    size_t nvariables = 0;

    std::vector<double> my_initial_x;
    for (int i = 0; i < n; ++i) {
        x[i] = 0.0;
    }

    for (const auto &key: ef_keys_) {
        auto& bmp_groups =  efficiency_[key];
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
    for (Index i = 0; i < n; i++) {
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
    for (int i = 0; i < nvars_; ++i) {
        x[i] = initial_x[i];
    }

    return true;
}

bool EPA_NLP::set_starting_point(
        const std::vector<double>& my_x
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
    for (int i = 0; i < nvars_; ++i) {
        //x[i] = initial_x[i];
        x[i] = 0.0; 
    }

    //start with greedy solution
    size_t idx = 0;
    for (const auto &key: ef_keys_) {
        auto& bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            bool flag = true;
            for (const auto &bmp : bmp_group) {
                if(flag == true) {
                    x[idx] = 1.0;
                    flag = false;
                }
                ++idx;
            }
        }
    }

    return true;
}

bool EPA_NLP::eval_f(
        Index n,
        const Number *x,
        bool new_x,
        Number &obj_value
) {
    assert(n == nvars_);
    double fitness = 0.0;
    int idx = 0;

    for (const auto &key: ef_keys_) {
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto load_src = out[2];
        auto state_id = lrseg_.at(lrseg)[1];
        auto alpha = amount_[key];
        auto& bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                double pct = x[idx];
                ++idx;
                if (pct < 0.0) {
                    pct = 0.0;
                }

                auto bmp_cost_key = fmt::format("{}_{}", state_id, bmp);
                double cost = pct * alpha * bmp_cost_[bmp_cost_key];
                if(cost<0.0){
                    fmt::print("cost: {}\n", cost);
                    fmt::print("pct: {}\n", pct);
                    fmt::print("alpha: {}\n", alpha);
                    fmt::print("bmp_cost_key: {} = {} \n", bmp_cost_key, bmp_cost_[bmp_cost_key]);
                    exit(0);
                }

                fitness += cost;
            }
        }
    }


    obj_value = fitness;
    assert(obj_value>0);
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
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto state_id = lrseg_.at(lrseg)[1];
        auto alpha = amount_[key];
        auto& bmp_groups =  efficiency_[key];
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                auto bmp_cost_key = fmt::format("{}_{}", state_id, bmp);
                grad_f[idx] =  alpha* bmp_cost_[bmp_cost_key];
                ++idx;
            }
        }
    }


    return true;
}

bool EPA_NLP::eval_g(
        Index n,
        const Number *x,
        bool new_x,
        Index m,
        Number *g
) {
    return eval_g_proxy( n, x, new_x, m, g, false );
}

bool EPA_NLP::eval_g_proxy(
        Index n,
        const Number *x,
        bool new_x,
        Index m,
        Number *g,
        bool is_final 
) {
    assert(n == nvars_);

    double total_cost = 0;
    std::vector<double> fx(2, 0.0);
    int nconstraints = 1;

    std::vector<double> pt_load(3, 0.0);
    size_t idx = 0;
    for (const auto &key: ef_keys_) {
        auto& bmp_groups =  efficiency_[key];
        std::vector<double> prod(3, 1);
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto load_src = out[2];
        for (const auto &bmp_group: bmp_groups) {
            std::vector<double> sigma(3, 0.0);
            auto sigma_cnstr = 0.0;
            for (const auto &bmp: bmp_group) {
                auto pct = x[idx];
                idx++;
                std::string s_tmp = fmt::format("{}_{}_{}", bmp, lrseg, load_src);
                std::vector<double> eta(3,0.0);
                if (eta_dict_.find(s_tmp) != eta_dict_.end()) {
                    eta[0] = eta_dict_[s_tmp][0];
                    eta[1] = eta_dict_[s_tmp][1];
                    eta[2] = eta_dict_[s_tmp][2];
                }
                sigma[0] += eta[0] * pct;
                sigma[1] += eta[1] * pct;
                sigma[2] += eta[2] * pct;
                sigma_cnstr += pct;
            }


            prod[0] *= (1.0 - sigma[0]);
            prod[1] *= (1.0 - sigma[1]);
            prod[2] *= (1.0 - sigma[2]);

            g[nconstraints] = sigma_cnstr;
            ++nconstraints;
        }
        if(!phi_dict_.contains(key)) {
            std::cout<<"nooooooo.... no key here\n"<<key<<"\n";
        }
        auto phi = phi_dict_[key];
        auto alpha = amount_[key];
        pt_load[0] += phi[0] * alpha * prod[0];
        pt_load[1] += phi[1] * alpha * prod[1];
        pt_load[2] += phi[2] * alpha * prod[2];

    }
    g[0] = pt_load[pollutant_idx];

    if (is_final == true) {
        g[0] = pt_load[0];
        g[1] = pt_load[1];
        g[2] = pt_load[2];
    }
    /*

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
    */

    /*
    //g[0] =  (sum_load_invalid_[0] + pt_load[0]) - 0.8*(sum_load_invalid_[0] + sum_load_valid_[0]) ;
    double load_ub = sum_load_invalid_[load_to_opt_] + sum_load_valid_[load_to_opt_];
    double load_obtained = sum_load_invalid_[load_to_opt_] + pt_load[load_to_opt_];

    double load_lb = 0.8 * load_ub; //20% reduction
    g[0] = (load_obtained < load_lb)?(load_obtained-load_lb)/load_ub:0;

    double cost_ub = 1000000.0;
    double cost_steps = 1000000.0;
    g[1] = (total_cost > cost_ub)? (cost_ub-total_cost)/cost_steps:0;
    fx[1] = sum_load_invalid_[load_to_opt_] + pt_load[load_to_opt_];
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
    assert(n == nvars_);
    if (values == NULL) {
        // return the structure of the Jacobian
        // this particular Jacobian is not dense
        
        for (int i = 0; i < nvars_; ++i) {
            iRow[i] = 0;
            jCol[i] = i;
        }

        int jac_index = nvars_;
        int jac_row = 1;
        for (const auto &key: ef_keys_) {
            auto& bmp_groups =  efficiency_[key];
            for (const auto &bmp_group: bmp_groups) {
                for (const auto &bmp: bmp_group) {
                        iRow[jac_index] = jac_row;
                        jCol[jac_index] = jac_index - nvars_;
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
    } else {
        // return the values of the Jacobian of the constraints
        size_t parcel_idx = 0;
        size_t idx = 0;
        size_t j = 0;

        for (const auto &key: ef_keys_) {
            std::vector<double> prod(3, 1);
            std::vector <std::string> out;
            misc_utilities::split_str(key, '_', out);
            std::vector<std::vector<double>> sigma_less;
            std::vector<std::vector<double>> sigma_full;
            auto lrseg = out[0];
            auto load_src = out[2];
            auto alpha = amount_[key];
            auto phi = phi_dict_[key];
            auto& bmp_groups =  efficiency_[key];

            size_t grp_counter = 0;
            for (const auto &bmp_group: bmp_groups) {
                std::vector<double> sigma(3, 0.0);
                auto saved_idx = idx;
                for (const auto &bmp: bmp_group) {
                    double pct = x[idx];
                    ++idx;
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp, lrseg, load_src);

                    if (eta_dict_.find(s_tmp) != eta_dict_.end()) {
                        auto eta = eta_dict_[s_tmp];
                        sigma[0] += eta[0] * pct;
                        sigma[1] += eta[1] * pct;
                        sigma[2] += eta[2] * pct;
                    }
                }
                idx = saved_idx;
                sigma_full.push_back({sigma[0], sigma[1], sigma[2]});
                for (const auto &bmp: bmp_group) {
                    double pct = x[idx];
                    ++idx;
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp, lrseg, load_src);
                    std::vector<double> eta(3,0.0);
                    if (eta_dict_.find(s_tmp) != eta_dict_.end()) {
                        eta[0] = eta_dict_[s_tmp][0];
                        eta[1] = eta_dict_[s_tmp][1];
                        eta[2] = eta_dict_[s_tmp][2];
                    }
                    sigma_less.push_back({sigma[0] - eta[0] * pct, sigma[1] - eta[1] * pct, sigma[2] - eta[2] * pct});
                }

                prod[0] *= 1.0 - sigma[0];
                prod[1] *= 1.0 - sigma[1];
                prod[2] *= 1.0 - sigma[2];
                ++grp_counter;
            }
            grp_counter = 0;
            size_t k = 0;
            for (const auto &bmp_group: bmp_groups) {
                for (const auto &bmp: bmp_group) {
                    std::string s_tmp = fmt::format("{}_{}_{}", bmp, lrseg, load_src);
                    std::vector<double> eta (3, 0.0);
                    if (eta_dict_.find(s_tmp) != eta_dict_.end()) {
                        eta[0] = eta_dict_[s_tmp][0];
                        eta[1] = eta_dict_[s_tmp][1];
                        eta[2] = eta_dict_[s_tmp][2];
                    }
                    values[j] = prod[pollutant_idx] / (1.0 - sigma_full[grp_counter][pollutant_idx]);
                    values[j] *= (1.0 - sigma_less[k][pollutant_idx]);
                    values[j] *= alpha * phi[pollutant_idx] * (-eta[pollutant_idx]);
                    
                    ++j;
                    ++k;
                }
                ++grp_counter;
            }
        }

        int var_counter = 0;

        for (const auto &key: ef_keys_) {
            auto& bmp_groups =  efficiency_[key];
            for (const auto &bmp_group: bmp_groups) {
                for (const auto &bmp: bmp_group) {
                    values[nvars_+ var_counter] = 1;
                    ++var_counter;
                }
            }
        }
        int limit_cnstr_counter = nvars_ * 2;

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
    assert(n == nvars_);

    if (values == NULL) {
        Index idx = 0;
        for (Index row = 0; row < nvars_; row++) {
            for (Index col = 0; col <= row; col++) {
                iRow[idx] = row;
                jCol[idx] = col;
                idx++;
            }
        }

        assert(idx == nele_hess);
    } else {
        int parcel_idx = 0;
        /*
        for (int j(0); j < nvars_;) {
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
                    sigma += eta_dict_[s_tmp][pollutant_idx] * pct;
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
                    if (eta_dict_.find(s_tmp) != eta_dict_.end()) {
                        eta_value = eta_dict_[s_tmp][pollutant_idx];
                    }
                    sigma_less[acc_idx + bc] = sigma - eta_value * pct;
                    //values[nvars_+cnstr_counter] = sigma_less[acc_idx + bc];
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
        */
    }

    return false;
}


void EPA_NLP::save_files(
        Index                      n,
        const Number *x
) {
    auto filename = fmt::format("{}_reduced.csv", filename_out);
    auto json_filename = fmt::format("{}.json", filename_out);

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

    for (const auto &key: ef_keys_) {
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto agency = out[1];
        auto load_src = out[2];
        auto state_id = lrseg_.at(lrseg)[1];
        auto alpha = amount_[key];
        auto& bmp_groups =  efficiency_[key];
        auto unit_id = 1;
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                if (x[idx] * alpha > 1.0) {
                    double amount = x[idx];
                    ofile<<fmt::format("{},{},{},{},{},{}\n", lrseg, agency, load_src, bmp, unit_id, amount);
                    v.name = fmt::format("{}_{}_{}_{}", lrseg, agency, load_src, bmp);
                    v.location = {std::stoi(lrseg), std::stoi(agency), std::stoi(load_src)};
                    v.amount = (amount>1)?1.0:((amount<0)?0.0:amount);
                    v.bmp = bmp;
                    this_json_ipopt.emplace_back(v);
                    flag = true;
                }
                ++idx;
            }
        }
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

    bool new_x = true;
    Number *g_constr = new Number[ncons_];
    eval_g_proxy(n, x, new_x, m, g_constr, true);

    std::cout.precision(10);
    std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
    std::cout << std::endl << std::endl << "Objective value" << std::endl;

    std::cout << "Cost(x*) = " << obj_value << std::endl;
    std::cout << std::endl << "Final value of the constraints:" << std::endl;


    std::string file_name = fmt::format("{}.csv", filename_out);
    std::ofstream file(file_name);
    file.precision(15);

    file_name = fmt::format("{}_results.txt", filename_out);
    std::ofstream file2(file_name);
    file2.precision(15);
    std::cout << "g_constr[NLoadEos] = " << g_constr[0] + sum_load_invalid_[0] << "\n";
    std::cout << "g_constr[PLoadEos] = " << g_constr[1] + sum_load_invalid_[1] << "\n";
    std::cout << "g_constr[SLoadEos] = " << g_constr[2] + sum_load_invalid_[2] << "\n";
    std::cout << "Original NLoadEos = " << sum_load_valid_[0] + sum_load_invalid_[0] << "\n";
    std::cout << "Original PLoadEos = " << sum_load_valid_[1] + sum_load_invalid_[1] << "\n";
    std::cout << "Original SLoadEos = " << sum_load_valid_[2] + sum_load_invalid_[2] << "\n";

    file2 << obj_value << "," << g_constr[0] + sum_load_invalid_[0] << "," << g_constr[1] + sum_load_invalid_[1] << ","
          << g_constr[2] + sum_load_invalid_[2] << "," << sum_load_valid_[0] + sum_load_invalid_[0] << ","
          << sum_load_valid_[1] + sum_load_invalid_[1] << "," << sum_load_valid_[2] + sum_load_invalid_[2] << "\n";
    file2.close();
    // Iterate through each line and split the content using delimeter
    int idx = 0;

    file<<"BmpSubmittedId,AgencyId,StateUniqueIdentifier,StateId,BmpId,GeographyId,LoadSourceGroupId,UnitId,Amount,isValid,ErrorMessage,RowIndex,Cost,LrsegId,LoadSourceIdOriginal,Acreage\n";
    int counter = 0;
    std::unordered_map<std::string, double> bmp_sum;
    total_cost = 0.0;

    for (const auto &key: ef_keys_) {
        std::vector <std::string> out;
        misc_utilities::split_str(key, '_', out);
        auto lrseg = out[0];
        auto agency = out[1];
        auto load_src = out[2];
        auto alpha = amount_[key];
        auto& bmp_groups =  efficiency_[key];
        auto state_id = lrseg_.at(lrseg)[1];
        auto geography = lrseg_.at(lrseg)[2];
        auto unit_id = 1;
        for (const auto &bmp_group : bmp_groups) {
            for (const auto &bmp : bmp_group) {
                if (x[idx] * alpha > 1.0) {
                    double amount = x[idx] * alpha;
                    auto bmp_cost_key = fmt::format("{}_{}", state_id, bmp);
                    double cost = amount * bmp_cost_[bmp_cost_key];

                    file << counter + 1 << "," << agency << ",SU" << counter << "," << state_id << "," << bmp << ","
                         << geography << "," << u_u_group_dict[load_src] << "," << unit_id << "," << amount << ",True,,"
                         << counter + 1 << "," << cost << "," << lrseg << "," << load_src << "," << alpha <<"\n";
                    counter++;
                    if (bmp_sum.find(bmp_cost_key) != bmp_sum.end()) {
                        bmp_sum[bmp_cost_key] += amount;
                    } else {
                        bmp_sum[bmp_cost_key] = amount;
                    }

                }
                ++idx;
            }
        }
    }

    file_name = fmt::format("{}_summary.txt", filename_out);

    std::ofstream file3(file_name);

    file3.precision(10);
    file3 << "Bmp_id,Acres,Cost" << std::endl;
    for (auto const&[key, sum]: bmp_sum) {
        file3 << key << "," << sum << "," <<  bmp_cost_[key] * sum <<  '\n';
    }
    file3.close();
    std::cout << "# Bmps: " << counter << std::endl;
    file.close();

    std::string filename_funcs= fmt::format("{}_funcs.txt", filename_out);
    std::ofstream file_funcs(filename_funcs, std::ios_base::app);
    file_funcs.precision(10);
    file_funcs<< obj_value << "," << g_constr[0] << " "<<g_constr[1]<<" "<< g_constr[2]<< "\n";
    file_funcs.close();
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
    std::cout<<"Status: "<<status<<std::endl;

    for (size_t i = 0; i < n; i++)
    {
        final_x[i] = x[i];
    }
    
    write_files(n, x, m, obj_value);
    save_files(n, x);
    // write_formar_cast_file();
}

// Author: Gregorio Toscano Pulido 
//epa exec_id 1.0-target_reduction target_pollutant algorithm_var 1
//epa 45 0.90 0 0 1

//[MAIN]
#include "IpIpoptApplication.hpp"
#include "nlp.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <unistd.h>
#include <cstdlib>
#include <string>
#include <fmt/core.h>
#include <coin-or/IpSmartPtr.hpp>
#include <unordered_map>
#include <memory>
#include <misc_utilities.h>

using namespace Ipopt;

extern int nvars;
extern int ncons;

// https://dev.to/zigabrencic/parameters-object-in-c-1an9
class Parameters {
public:
    std::string filename_in;
    std::string filename_out;
    double max_constr;
    int pollutant_idx;
    std::string log_filename;
    Parameters() {
    }
};

class EpsConstraint {
    // https://pdimov.github.io/blog/2020/09/07/named-parameters-in-c20/
    struct Params {
        std::string filename_in;
        std::string filename_out;
        double max_constr;
        int pollutant_idx;
        std::string log_filename;
    };
    //SmartPtr <TNLP> mynlp;
    SmartPtr <EPA_NLP> mynlp;
    SmartPtr <IpoptApplication> app;
    std::string log_filename;
    //EPA_NLP *mynlp;
public:
    EpsConstraint(const std::shared_ptr<Parameters> &var);
    bool constr_eval(double, int);
    bool evaluate(double);
    void get_final_x(std::vector<double> &final_x);
    void set_init_x(std::vector<double> &init_x);
};

void EpsConstraint::get_final_x(std::vector<double> &final_x) {
    final_x.resize(mynlp->final_x.size());
    std::copy(mynlp->final_x.begin(), mynlp->final_x.end(), final_x.begin());
}

void EpsConstraint::set_init_x(std::vector<double> &init_x) {
    mynlp->init_x.resize(init_x.size());
    std::copy(init_x.begin(), init_x.end(), mynlp->init_x.begin());
}

EpsConstraint::EpsConstraint(const std::shared_ptr<Parameters> &var) {
    std::string filename_in = var->filename_in;
    std::string filename_out = var->filename_out;
    double max_constr = var->max_constr;
    int pollutant_idx = var->pollutant_idx;
    log_filename = var->log_filename;


    mynlp = new EPA_NLP(filename_in, filename_out, max_constr, pollutant_idx);
    app = IpoptApplicationFactory();
}

bool EpsConstraint::evaluate(double reduction) {
    //app->Options()->SetNumericValue("tol", 1e-8);

    mynlp->update_reduction(reduction);
    app->Options()->SetIntegerValue("max_iter", 1000);
    app->Options()->SetStringValue("linear_solver", "ma57");

    app->Options()->SetStringValue("output_file", log_filename.c_str());
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

    status = app->OptimizeTNLP(mynlp);
    //mynlp->save_files(n, x);
    return status;
}

bool EpsConstraint::constr_eval(double reduction, int nsteps){
    double step_size = reduction/nsteps;
    for (int i(1); i<= nsteps; ++i) {
        double max_constr = step_size*i;
        auto result = evaluate(max_constr);
    }
    return true;
}

int main(
        int argc,
        char **argv
) {
    //        1             2           3  4  5   6
    //run Jefferson.json   output.json  0 0.9 1.0 5 

    std::string filename_in;
    std::string filename_out;
    filename_in = argv[1];
    filename_out = argv[2];
    int pollutant_idx = atoi(argv[3]); // 0
    double max_constr = atof(argv[4]); // 0.9
    int nsteps = atoi(argv[5]); //5

    auto var = std::make_shared<Parameters>();
    std::string log_filename = fmt::format("{}_ipopt.out", filename_out);
    var->filename_in = filename_in;
    var->filename_out = filename_out;
    var->log_filename = log_filename;
    var->max_constr = max_constr;
    var->pollutant_idx = pollutant_idx;

    int option = 0;
    if (option == 0) { //ipopt Opt3
        EpsConstraint eps_constr(var);
        auto r = eps_constr.evaluate(var->max_constr);
        /*
        std::cout<<r<<std::endl;

        //std::vector<double> init_x;
        //eps_constr.get_final_x(init_x);
        //eps_constr.set_init_x(init_x);

        r = eps_constr.evaluate(0.7);
        std::cout<<r<<std::endl;
        */
    }
    if (option == 1) { //ipopt Opt3
        EpsConstraint eps_constr(var);
        eps_constr.constr_eval(var->max_constr, nsteps);
    }
    if (option == 6) { //ipopt Opt3
        SmartPtr <TNLP> mynlp = new EPA_NLP(filename_in, filename_out, max_constr, pollutant_idx);
        SmartPtr <IpoptApplication> app = IpoptApplicationFactory();

        //app->Options()->SetNumericValue("tol", 1e-8);
        app->Options()->SetIntegerValue("max_iter", 10000);
        app->Options()->SetStringValue("linear_solver", "ma57");

        //app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue("output_file", log_filename.c_str());
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");

        ApplicationReturnStatus status;
        status = app->Initialize();
        if (status != Solve_Succeeded) {
            std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
            return (int) status;
        }

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);
        std::ofstream file(log_filename, std::ios_base::app);

        if (status == Solve_Succeeded) {
            std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
            file << "Solved\n";
        } else {
            std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
            file << "Failed\n";
        }
        file.close();

    } else if (option == 2) { //random Opt1
        EPA_NLP *mynlp = new EPA_NLP(filename_in, filename_out, max_constr, pollutant_idx);
        Index n;
        Index m;
        mynlp->get_n_and_m(n, m);

        Number *x = new Number[n];
        Number obj_value;

        mynlp->random_x(n, m, x);
        //get_info(n, m);
        mynlp->my_get_starting_point(x);
        for(int i(0); i<n; i++){
            if(x[i]>0.1)
                std::cout<<i<<":"<<x[i]<<" ";
        }
        std::cout<<std::endl;
        mynlp->eval_f(n, x, true, obj_value);
        mynlp->write_files(n, x, m, obj_value);

        std::ofstream file(log_filename, std::ios_base::app);

        file << "Random Test\n";
        file.close();
        delete x;

    }
    return 0;
}


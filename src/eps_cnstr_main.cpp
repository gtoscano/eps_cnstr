// Author: Gregorio Toscano Pulido 
//epa exec_id 1.0-target_reduction target_pollutant algorithm_var 1
//epa 45 0.90 0 0 1

//[MAIN]

#include "IpIpoptApplication.hpp"
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
#include "nlp.hpp"
#include <eps_cnstr.h>


using namespace Ipopt;

int main(
        int argc,
        char **argv
) {
    //        1             2           3  4  5   6
    //run Jefferson.json   output.json  0 0.9 5 
    //run input_base_scenario.json input_scenario.json output.json pollutant_id  target_reduction nsteps 
    std::string filename_in;
    std::string filename_out;
    std::string filename_scenario;
    filename_in = argv[1];
    filename_scenario = argv[2];
    fmt::print("filename_scenario: {}\n", filename_scenario);
    filename_out = argv[3];
    int pollutant_idx = atoi(argv[4]); // 0
    double reduction= atof(argv[5]); // 0.9
    int nsteps = atoi(argv[6]); //5

    auto var = std::make_shared<Parameters>();
    std::string log_filename = fmt::format("{}_ipopt.out", filename_out);
    var->filename_in = filename_in;
    var->filename_out = filename_out;
    var->filename_scenario = filename_scenario;
    var->log_filename = log_filename;
    var->reduction = reduction;
    var->pollutant_idx = pollutant_idx;

    int option = 1;
    /*
    if (option == 0) { //ipopt Opt3
        EpsConstraint eps_constr(var);
        auto r = eps_constr.evaluate(var->reduction);
        //std::cout<<r<<std::endl;

        //std::vector<double> init_x;
        //eps_constr.get_final_x(init_x);
        //eps_constr.set_init_x(init_x);

        //r = eps_constr.evaluate(0.7);
        //std::cout<<r<<std::endl;
    }
    */
    if (option == 1) { //ipopt Opt3
        EpsConstraint eps_constr(var);
        eps_constr.constr_eval(var->reduction, nsteps);
    }
    /*
    if (option == 6) { //ipopt Opt3
        SmartPtr <TNLP> mynlp = new EPA_NLP(filename_in, filename_scenario, filename_out, reduction, pollutant_idx);
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
        EPA_NLP *mynlp = new EPA_NLP(filename_in, filename_scenario, filename_out, reduction, pollutant_idx);
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

    }*/
    return 0;
}


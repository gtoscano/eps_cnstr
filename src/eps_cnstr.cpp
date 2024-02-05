// Created by: Gregorio Toscano Pulido
#include <iostream>
#include <string>
#include <vector>
#include "eps_cnstr.h"

EpsConstraint::EpsConstraint(const std::shared_ptr<Parameters> &var) {
    std::string filename_in = var->filename_in;
    std::string filename_scenario = var->filename_scenario;
    std::string filename_out = var->filename_out;
    double reduction= var->reduction;
    int pollutant_idx = var->pollutant_idx;
    log_filename = var->log_filename;


    mynlp = new EPA_NLP(filename_in, filename_scenario, filename_out, reduction, pollutant_idx);
    app = IpoptApplicationFactory();
}

bool EpsConstraint::evaluate(double reduction, int current_iteration=1) {
    //app->Options()->SetNumericValue("tol", 1e-8);

    mynlp->update_reduction(reduction, current_iteration);
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
        double lower_bound = step_size*i;
        auto result = evaluate(lower_bound, i);
    }
    return true;
}

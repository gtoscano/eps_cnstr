// Created by gtoscano on 4/1/23.
//

#ifndef CBO_EVALUATION_MISC_UTILITIES_H
#define CBO_EVALUATION_MISC_UTILITIES_H

#include <vector>

namespace misc_utilities {
    /**
    * reads environment variables
    *
    * This function reads an environment variable.
    * If the variable is not set it returns a default value
    *
    * @param key the environment variable's name
    * @param  default_value The default value
    * @return either the environment variable or the default value.
    */
    std::string get_env_var(std::string const &key, std::string const &default_value);

    std::string find_file(std::string path,
                          std::string prefix);

    std::vector<std::string> find_files(std::string path, std::string prefix);
    void split_str(std::string const &str,
                   const char delim,
                   std::vector<std::string> &out);

    bool copy_full_directory(const std::string& source, const std::string& destination);
    bool copy_file(const std::string& source,
                   const std::string& destination);
    std::string current_time();

    double rand_double(double lower_bound, double upper_bound);
    void mkdir(std::string dir_path);
}

#endif //CBO_EVALUATION_MISC_UTILITIES_H

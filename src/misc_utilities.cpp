//
// Created by Gregorio Toscano on 4/1/23.
//

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include <random>
#include <fstream>


namespace fs = std::filesystem;

namespace misc_utilities {
    bool is_parquet_file(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
    
        if (!file) {
            std::cerr << "Failed to open the file: "<<filename<<".\n";
            return false;
        }
    
        char header[4];
        char footer[4];
    
        file.read(header, sizeof(header));
    
        if (!file.seekg(-4, std::ios::end)) {
            std::cerr << "File too small.\n";
            return false;
        }
    
        file.read(footer, sizeof(footer));
    
        return std::string_view(header, 4) == "PAR1" && std::string_view(footer, 4) == "PAR1";
    }

    std::string get_env_var(std::string const &key, std::string const &default_value) {
        const char *val = std::getenv(key.c_str());
        return val == nullptr ? std::string(default_value) : std::string(val);
    }

    std::string find_file(std::string path, std::string prefix) {
        std::string found_filename;
        for (const auto& entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file()
                && entry.path().extension() == ".csv"
                && entry.path().stem().string().starts_with(prefix)) {
                found_filename = entry.path();
                break;
            }
        }
        return found_filename;
    }
    std::vector<std::string> find_files(std::string path, std::string prefix) {
        std::vector<std::string> found_filenames;
        for (const auto& entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file()
                //&& entry.path().extension() == ".csv"
                && entry.path().stem().string().starts_with(prefix)) {
                found_filenames.emplace_back(entry.path().filename());
            }
        }
        return found_filenames;
    }


    void split_str(std::string const &str, const char delim,
                   std::vector<std::string> &out) {
        out.clear();

        std::stringstream s(str);

        std::string s2;

        while (std::getline(s, s2, delim)) {
            out.push_back(s2); // store the string in s2
        }
    }

    bool copy_full_directory(const std::string& source, const std::string& destination) {
        try {
            // Check if the source directory exists
            if (!fs::exists(source)) {
                std::cerr << "Error: Source directory does not exist" << std::endl;
                return false;
            }
    
            // Check if the source is actually a directory
            if (!fs::is_directory(source)) {
                std::cerr << "Error: Source is not a directory" << std::endl;
                return false;
            }
    
            // Create the destination directory if it doesn't exist
            if (!fs::exists(destination)) {
                fs::create_directories(destination);
            }
    
            // Copy the entire directory
            fs::copy(source, destination, fs::copy_options::recursive | fs::copy_options::overwrite_existing);
            
            return true;
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << std::endl;
            return false;
        }
    }

    bool copy_file(const std::string& source, const std::string& destination) {
        try {
            // Check if the source file exists
            std::cout << "source: " << source << std::endl;
            if (!fs::exists(source)) {
                std::cerr << "Error: Source file does not exist" << std::endl;
                return false;
            }

            // Copy the file to the destination
            fs::copy_file(source, destination, fs::copy_options::overwrite_existing);
            return true;
        } catch (const std::exception& ex) {
            std::cerr << "Error: " << ex.what() << std::endl;
            return false;
        }
    }
    std::string current_time() {
        // Get the current time
        auto now = std::chrono::system_clock::now();

        // Convert the time to a tm struct
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm now_tm = *std::localtime(&now_c);

        // Format the time as a string
        std::ostringstream oss;
        oss << std::put_time(&now_tm, "%H:%M:%S");
        return oss.str();
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    double rand_double(double lower_bound, double upper_bound) {
        return lower_bound + dis(gen) * (upper_bound - lower_bound);
    }

    void mkdir(std::string dir_path) {
    
        if (!fs::exists(dir_path)) {
            if (fs::create_directories(dir_path)) {
                std::cout << "Directory created successfully." << std::endl;
            } else {
                std::cerr << "Failed to open the directory: "<<dir_path<<".\n";
            }
        } else {
            std::cout << "Directory already exists. Doing nothing." << std::endl;
        }
    }

    /*
    double compute_cost(const std::string& loads_filename, int cost_profile_id) {
        std::shared_ptr<arrow::io::ReadableFile> infile;
        PARQUET_ASSIGN_OR_THROW(infile, arrow::io::ReadableFile::Open(loads_filename, arrow::default_memory_pool()));
        std::unique_ptr<parquet::arrow::FileReader> reader;
        PARQUET_THROW_NOT_OK(parquet::arrow::OpenFile(infile, arrow::default_memory_pool(), &reader));
    
        std::shared_ptr<arrow::ChunkedArray> amountArray;
        std::shared_ptr<arrow::ChunkedArray> bmpIdArray;
    
        // Get the schema from the Parquet file
        std::shared_ptr<arrow::Schema> schema;
        PARQUET_THROW_NOT_OK(reader->GetSchema(&schema));
        
        // Find the index of the "Amount" column
        int amount_index = schema->GetFieldIndex("Amount");
        if (amount_index == -1) {
            throw std::runtime_error("Amount column not found in Parquet file.");
        }
        
        // Find the index of the "BmpId" column
        int bmpId_index = schema->GetFieldIndex("BmpId");
        if (bmpId_index == -1) {
            throw std::runtime_error("BmpId column not found in Parquet file.");
        }
        
        // Read the "Amount" column by index
        PARQUET_THROW_NOT_OK(reader->ReadColumn(amount_index, &amountArray));
        
        // Read the "BmpId" column by index
        PARQUET_THROW_NOT_OK(reader->ReadColumn(bmpId_index, &bmpIdArray));
        //PARQUET_THROW_NOT_OK(reader->GetColumnByName("Amount", &amountArray));
        //PARQUET_THROW_NOT_OK(reader->GetColumnByName("BmpId", &bmpIdArray));
    
        if (amountArray->length() != bmpIdArray->length()) {
            throw std::runtime_error("Length mismatch between Amount and BmpId columns.");
        }
    
        double total_cost = 0.0;
    
        // Assuming the columns are split into chunks of equal size
        for (int i = 0; i < amountArray->num_chunks(); i++) {
            auto amountValues = std::static_pointer_cast<arrow::DoubleArray>(amountArray->chunk(i));
            auto bmpIdValues = std::static_pointer_cast<arrow::Int32Array>(bmpIdArray->chunk(i));
    
            for (int64_t j = 0; j < amountValues->length(); j++) {
                if (!amountValues->IsNull(j) && !bmpIdValues->IsNull(j)) {
                    int bmpId = bmpIdValues->Value(j);
                    double amount = amountValues->Value(j);
    
                    auto cost_per_unit = data_.get_bmp_cost(cost_profile_id, bmpId);
                    total_cost += cost_per_unit * amount;
                }
            }
        }
    
        return total_cost;
    }
    */

}

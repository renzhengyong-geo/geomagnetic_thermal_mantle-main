#include <iostream>
#include <string>
#include <vector>
#include <dirent.h>
#include <cstdlib>  // for exit
#include <cstdio>   // for remove

// Function to remove the ".dat" suffix if present
std::string remove_dat_suffix(const std::string& base_name) {
    const std::string suffix = ".dat";
    if (base_name.size() >= suffix.size() && base_name.substr(base_name.size() - suffix.size()) == suffix) {
        // Remove the ".dat" suffix
        return base_name.substr(0, base_name.size() - suffix.size());
    }
    return base_name; // Return unchanged if no ".dat" suffix
}


// Function to delete additional file types
void delete_additional_file_types() {
    // Define the list of extensions to delete
    std::vector<std::string> extensions = {".tab", ".plt", ".arf", ".blk", ".prn", ".tim", ".tof", "auto_refine.txt", "seismic_data.txt"};

    DIR *dir;
    struct dirent *ent;

    // Open the current directory
    if ((dir = opendir(".")) != NULL) {
        // Iterate over all files in the directory
        while ((ent = readdir(dir)) != NULL) {
            std::string file_name = ent->d_name;

            // Check if the file has any of the specified extensions
            for (const auto& ext : extensions) {
                if (file_name.size() >= ext.size() && file_name.substr(file_name.size() - ext.size()) == ext) {
                    // Delete the file
                    std::cout << "Deleting file: " << file_name << std::endl;
                    if (remove(file_name.c_str()) != 0) {
                        std::cerr << "Error deleting file: " << file_name << std::endl;
                    }
                }
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Error opening directory!" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return; 
}

// Function to find the latest .tab file and delete older ones
std::string get_latest_tab_file_and_delete_older(const std::string& base_name) {
    int max_index = 0;
    std::string latest_file_name;

    // Remove the ".dat" suffix from the base_name
    std::string stripped_base_name = remove_dat_suffix(base_name);

    DIR *dir;
    struct dirent *ent;

    // Open the current directory
    if ((dir = opendir(".")) != NULL) {
        // Iterate over all files in the directory
        while ((ent = readdir(dir)) != NULL) {
            std::string file_name = ent->d_name;

            // Check if the file name starts with the stripped base_name and has the format "a_N.tab"
            if (file_name.find(stripped_base_name) == 0 && file_name.substr(stripped_base_name.size(), 1) == "_") {
                // Find the index part of the filename
                size_t start_pos = stripped_base_name.size() + 1; // Skip the base name and underscore
                size_t end_pos = file_name.find(".tab");
                if (end_pos != std::string::npos) {
                    try {
                        int file_index = std::stoi(file_name.substr(start_pos, end_pos - start_pos));
                        if (file_index > max_index) {
                            max_index = file_index;
                            latest_file_name = file_name;  // Update latest file
                        }
                    } catch (const std::invalid_argument&) {
                        // Handle cases where conversion to int fails
                        continue;
                    }
                }
            }
        }
        closedir(dir);
    } else {
        // Could not open directory
        std::cerr << "Error opening directory!" << std::endl;
        exit(EXIT_FAILURE);
    }

    // If no latest file found, return empty string
    if (latest_file_name.empty()) {
        return latest_file_name;
    }

/*
    // Now delete the older files (files with index < max_index)
    if ((dir = opendir(".")) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            std::string file_name = ent->d_name;

            // Check if the file name starts with stripped_base_name and has the format "a_N.tab"
            if (file_name.find(stripped_base_name) == 0 && file_name.substr(stripped_base_name.size(), 1) == "_") {
                size_t start_pos = stripped_base_name.size() + 1;
                size_t end_pos = file_name.find(".tab");
                if (end_pos != std::string::npos) {
                    try {
                        int file_index = std::stoi(file_name.substr(start_pos, end_pos - start_pos));
                        if (file_index < max_index) {
                            // Delete the old file
                            std::cout << "Deleting old file: " << file_name << std::endl;
                            if (remove(file_name.c_str()) != 0) {
                                std::cerr << "Error deleting file: " << file_name << std::endl;
                            }
                        }
                    } catch (const std::invalid_argument&) {
                        // Ignore files that do not match the pattern
                        continue;
                    }
                }
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Error opening directory!" << std::endl;
        exit(EXIT_FAILURE);
    }
*/
    // Delete files with following list of extensions
    // ".tab", ".plt", ".arf", ".blk", ".prn", ".tim", ".tof"
    // delete_additional_file_types();
    return latest_file_name;
}



/*
int main() {
    // Base name of the files (e.g., "a.dat")
    std::string base_name = "lower_data_file_v1.dat";

    // Get the latest tab file and delete older files
    std::string latest_file = get_latest_tab_file_and_delete_older(base_name);

    if (!latest_file.empty()) {
        std::cout << "Latest .tab file: " << latest_file << std::endl;
    } else {
        std::cout << "No .tab files found!" << std::endl;
    }

    return 0;
}
*/

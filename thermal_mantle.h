// THERMAL_MANTLE.h
#ifndef thermal_mantle_h
#define thermal_mantle_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <cassert>
#include <iomanip>  // For std::setw


struct Layer {
    double upper_bound;
    double lower_bound;
    // Add other member variables as needed
};

double findPressureForDepth(double depth); 
void modify_perplex_datafile(const std::string &inputFileName, const std::string &outputFileName, double basalt_fraction);
void run_vertex_werami_with_parameters_file(const std::string &programName, const std::string &inputFileName);
void read_matrix_from_perplex_tabfile(const std::string &file_name, std::vector<std::vector<double>> &matrix, std::vector<std::string> &column_names);
std::map<std::string, std::vector<double>> combine_data_to_map(const std::vector<std::vector<double>> &matrix, const std::vector<std::string>&column_names,std::vector<std::string> &keywords);
void write_map_to_file(const std::string &output_file_name, const std::map<std::string, std::vector<double>> &data_map, const std::vector<std::string> &keywords);
std::vector<Layer> generateEarthLayers(int n_um, int n_tz, int n_lm);
void printEarthLayers(const std::vector<Layer>& layers, int n_um, int n_tz, int n_lm);
std::string get_latest_tab_file_and_delete_older(const std::string& base_name);
void combineMantleLayers(const std::vector<std::vector<double>>& matrix_upper, const std::vector<std::string>& column_names_upper,  const std::vector<std::vector<double>>& matrix_middle,   const std::vector<std::string>& column_names_middle,   const std::vector<std::vector<double>>& matrix_lower,   const std::vector<std::string>& column_names_lower,   std::vector<std::vector<double>>& combined_matrix, std::vector<std::string>& combined_column_names);
void delete_additional_file_types(); 

#endif

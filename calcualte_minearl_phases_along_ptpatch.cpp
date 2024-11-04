#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <cassert>
#include <iomanip>  // For std::setw

// Declare this function and data type defined in other files
#include "thermal_mantle.h"

// Function to generate PT patch and write to mantle geothermal files
void generate_PT_patch(int n_um, int n_tz, int n_lm,
                       std::vector<Layer>& layer,
                       std::vector<double>& temperature_profile,
                       std::string upper_PT,
                       std::string middle_PT,
                       std::string lower_PT) {
    // Check if the layer and temperature profiles have correct sizes
    assert(layer.size() == temperature_profile.size());
    assert(layer.size() == n_um + n_tz + n_lm);

    // Open output files for writing
    std::ofstream upper_PT_outputfile(upper_PT);
    std::ofstream middle_PT_outputfile(middle_PT);
    std::ofstream lower_PT_outputfile(lower_PT);

    if (!upper_PT_outputfile.is_open() || !middle_PT_outputfile.is_open() || !lower_PT_outputfile.is_open()) {
        std::cerr << "Error: One or more output files could not be opened." << std::endl;
        return;
    }

    // Write the pressure-temperature values for each layer to the corresponding files
    for (size_t i = 0; i < layer.size(); ++i) {
        double depth = (layer[i].upper_bound + layer[i].lower_bound) * 0.5; // Calculate midpoint depth
        double pressure = findPressureForDepth(depth)*10000;                // Calculate pressure based on depth
                                                                            // Gpa to Bar, 1Gpa=10000Bar
        double temperature = temperature_profile[i];                        // Get temperature from profile

        // Ensure temperature is within valid range (Kelvin)
        assert(temperature > 273.15);           // Ensure temperature is above 0°C (273.15 K)
        assert(temperature < 3500 + 273.15);    // Ensure temperature is below 3500°C (3773.15 K)

        // Write pressure-temperature values to appropriate file based on layer index
        if (i < n_um) {
            upper_PT_outputfile << pressure << "\t" << temperature << "\n";
        } else if (i < n_um + n_tz) {
            middle_PT_outputfile << pressure << "\t" << temperature << "\n";
        } else {
            lower_PT_outputfile << pressure << "\t" << temperature << "\n";
        }
    }

    // Close the output files
    upper_PT_outputfile.close();
    middle_PT_outputfile.close();
    lower_PT_outputfile.close();
}

// Function to compute mineral phases along the PT patch
void compute_mineralphases_along_PT_patch(double Xum, double Xtz, double Xlm, 
                                          std::string upper_data_file_v0, 
                                          std::string middle_data_file_v0,
                                          std::string lower_data_file_v0,
                                          std::string upper_data_file_v1,
                                          std::string middle_data_file_v1,
                                          std::string lower_data_file_v1,
                                          std::string upper_mantle_vertex_parameters_file,
                                          std::string upper_mantle_werami_parameters_file,
                                          std::string middle_mantle_vertex_parameters_file,
                                          std::string middle_mantle_werami_parameters_file,
                                          std::string lower_mantle_vertex_parameters_file,
                                          std::string lower_mantle_werami_parameters_file) {

    // Modify the Perplex data files for each mantle region
    modify_perplex_datafile(upper_data_file_v0, upper_data_file_v1, Xum);
    modify_perplex_datafile(middle_data_file_v0, middle_data_file_v1, Xtz);    
    modify_perplex_datafile(lower_data_file_v0, lower_data_file_v1, Xlm);

    std::cout<<"\nmodify_perplex_datafile good!\n";
    
    // Run thermodynamic phase diagram calculations using Vertex and Werami programs
    std::string vertexProgram = "vertex";
    std::string weramiProgram = "werami";

    run_vertex_werami_with_parameters_file(vertexProgram, upper_mantle_vertex_parameters_file);
    run_vertex_werami_with_parameters_file(weramiProgram, upper_mantle_werami_parameters_file);

    run_vertex_werami_with_parameters_file(vertexProgram, middle_mantle_vertex_parameters_file);
    run_vertex_werami_with_parameters_file(weramiProgram, middle_mantle_werami_parameters_file);

    run_vertex_werami_with_parameters_file(vertexProgram, lower_mantle_vertex_parameters_file);
    run_vertex_werami_with_parameters_file(weramiProgram, lower_mantle_werami_parameters_file);

    // Define keywords and process data for each mantle region
    std::vector<std::string> keywords = {"node#", "P(bar)", "T(K)", "O", "Cpx", "Opx", "Gt", "C2/c", "Wad", "st", "Ring", 
                                         "ca-pv", "Aki", "Pv", "per", "CF", "wus","sp", "Ppv"};

    // Process upper mantle data
    std::vector<std::vector<double>> matrix_upper;
    std::vector<std::string> column_names_upper;
    read_matrix_from_perplex_tabfile(upper_data_file_v1, matrix_upper, column_names_upper);  
    std::map<std::string, std::vector<double>> data_map_upper = combine_data_to_map(matrix_upper, column_names_upper, keywords);
    // write_map_to_file("mineral_phase_upper.dat", data_map_upper, keywords);

    // Process middle mantle data
    std::vector<std::vector<double>> matrix_middle;
    std::vector<std::string> column_names_middle;
    read_matrix_from_perplex_tabfile(middle_data_file_v1, matrix_middle, column_names_middle);
    std::map<std::string, std::vector<double>> data_map_middle = combine_data_to_map(matrix_middle, column_names_middle, keywords);
    // write_map_to_file("mineral_phase_middle.dat", data_map_middle, keywords);
    
    // Process lower mantle data
    std::vector<std::vector<double>> matrix_lower;
    std::vector<std::string> column_names_lower;
    read_matrix_from_perplex_tabfile(lower_data_file_v1, matrix_lower, column_names_lower);
    std::map<std::string, std::vector<double>> data_map_lower = combine_data_to_map(matrix_lower, column_names_lower, keywords);
    // write_map_to_file("mineral_phase_lower.dat", data_map_lower, keywords);

    // Final combined matrix and column names
    std::vector<std::vector<double>> combined_matrix;
    std::vector<std::string> combined_column_names;
    combineMantleLayers(matrix_upper, column_names_upper, matrix_middle, column_names_middle, matrix_lower, column_names_lower, combined_matrix, combined_column_names);  
    std::string combined_tab_file   = "combined_tab_file.tab";
    std::ofstream output_file(combined_tab_file.c_str());
   // Set a fixed width for the output columns
   int column_width = 12;
   output_file << std::left;  // Align text to the left
    // Output column headers with proper spacing
   for (size_t i = 0; i < combined_column_names.size(); ++i) {
     output_file << std::setw(column_width) << combined_column_names[i];  // Set width for each column
   }
   output_file << std::endl;
   // Output table rows with the same formatting
   for (size_t i = 0; i < combined_matrix.size(); ++i) {
      for (size_t j = 0; j < combined_matrix[i].size(); ++j) {
          output_file << std::setw(column_width) << combined_matrix[i][j];  // Set width for each value
      }
    output_file << std::endl;
   }
    
   return;
}



int main() {
    // Define the number of layers for each region
    int n_um = 8;  // Number of layers in the upper mantle (0 km to 410 km)
    int n_tz = 5;  // Number of layers in the transition zone (410 km to 660 km)
    int n_lm = 12; // Number of layers in the lower mantle (660 km to 2900 km)

    // Generate Earth's layer structure
    std::vector<Layer> earth_layers = generateEarthLayers(n_um, n_tz, n_lm);
    // Print the Earth's layer structure in a formatted way
    printEarthLayers(earth_layers, n_um, n_tz, n_lm);
    
    std::vector<double> temperature_profile; 
    // Generate the temparuture prifle with a constant value
    for(int i=0; i<n_um+n_tz+n_lm; i++) temperature_profile.push_back(1500); // 1500 K
    
    // A section of defing the file names
    // P-T patches
    std::string upper_PT   = "upper_PT.dat";
    std::string middle_PT  = "middle_PT.dat";
    std::string lower_PT   = "lower_PT.dat";
    // input dat file (Generated by Build program)
    std::string upper_data_file_v0   = "upper_data_file_v0.dat";
    std::string middle_data_file_v0  = "middle_data_file_v0.dat";
    std::string lower_data_file_v0   = "lower_data_file_v0.dat";
    // output .dat file (will be used by Vertex and Werami programs)
    std::string upper_data_file_v1   = "upper_data_file_v1.dat";
    std::string middle_data_file_v1  = "middle_data_file_v1.dat";
    std::string lower_data_file_v1   = "lower_data_file_v1.dat";
    // input parameters for Vertex and Werami programs
    std::string upper_mantle_vertex_parameters_file = "upper_mantle_vertex_parameters_file";
    std::string upper_mantle_werami_parameters_file = "upper_mantle_werami_parameters_file";
    std::string middle_mantle_vertex_parameters_file= "middle_mantle_vertex_parameters_file";
    std::string middle_mantle_werami_parameters_file= "middle_mantle_werami_parameters_file";
    std::string lower_mantle_vertex_parameters_file = "lower_mantle_vertex_parameters_file";
    std::string lower_mantle_werami_parameters_file = "lower_mantle_werami_parameters_file";
                                                     
    // Generate the three PT patch (geothermal files)                                         
    generate_PT_patch(n_um, n_tz, n_lm, earth_layers, temperature_profile,
                      upper_PT, middle_PT, lower_PT);
    double Xum=0.2;   // fraction of baslate in upper mantle
    double Xtz=0.2;   // fraction of basalt in the middle mantle
    double Xlm=0.2;   // fraction of basalt in the lower mantle
    // Call Vertex to do the thermaldynatic calcuation along P-T patches
    // Call Werami to extrac the phase diagram or model abundance of minerals along P-T pathces
    
    // Clear all previously generated internal files generated by runnning vertex and werami
    // Delete files with following list of extensions
    // ".tab", ".plt", ".arf", ".blk", ".prn", ".tim", ".tof"
    delete_additional_file_types();
    
    // do the calcuation with different percentages of baslat in the mantle
    compute_mineralphases_along_PT_patch(Xum, Xtz, Xlm, 
    	  				upper_data_file_v0, middle_data_file_v0,lower_data_file_v0,
                                        upper_data_file_v1, middle_data_file_v1,lower_data_file_v1,
                                        upper_mantle_vertex_parameters_file, upper_mantle_werami_parameters_file, 
                                        middle_mantle_vertex_parameters_file,middle_mantle_werami_parameters_file,
                                        lower_mantle_vertex_parameters_file, lower_mantle_werami_parameters_file);    
    return 0;
}


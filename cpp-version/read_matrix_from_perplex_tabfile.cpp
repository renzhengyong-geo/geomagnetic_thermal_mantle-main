#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>    // For atof()
#include <iomanip>    // For formatting output


std::string get_latest_tab_file_and_delete_older(const std::string& base_name);


// Function to split a line of text into a vector of strings
std::vector<std::string> split_line(const std::string &line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string word;
    while (stream >> word) {
        result.push_back(word);
    }
    return result;
}

// Function to read matrix data and column names from a Perplex tab file
void read_matrix_from_perplex_tabfile(const std::string &file_name, std::vector<std::vector<double>> &matrix, std::vector<std::string> &column_names) {
    
    // read Perple_X table file
    // Get the latest tab file and delete older files
    std::string latest_tab_file = get_latest_tab_file_and_delete_older(file_name);
    if (!latest_tab_file.empty()) {
        std::cout << "Latest .tab file: " << latest_tab_file << std::endl;
    } else {
        std::cout << "No .tab files found!" << std::endl;
    }
    std::ifstream file(latest_tab_file.c_str());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << latest_tab_file << std::endl;
        return;
    }

    std::string line;

    // Skip the first 6 lines (metadata and other information)
    for (int i = 0; i < 6; ++i) std::getline(file, line);

    // Read the number of rows and columns
    int row_count = 0, col_count = 0;
    std::getline(file, line); row_count = std::atoi(line.c_str());
    std::getline(file, line); col_count = std::atoi(line.c_str());

    // Read column names
    std::getline(file, line);
    column_names = split_line(line);

    // Initialize matrix with zero values
    matrix.assign(row_count, std::vector<double>(col_count, 0.0));

    // Read and parse each row of the matrix
    for (int i = 0; i < row_count; ++i) {
        std::getline(file, line);
        std::vector<std::string> columns = split_line(line);
        for (int j = 0; j < col_count; ++j) {
            if (columns[j] != "NaN")
                matrix[i][j] = std::atof(columns[j].c_str());
        }
    }
    file.close();
}

// Combine column names and matrix data into a map based on predefined keywords
std::map<std::string, std::vector<double>> combine_data_to_map(
    const std::vector<std::vector<double>> &matrix,
    const std::vector<std::string> &column_names,
    std::vector<std::string> &keywords) {  // Pass keywords by reference to modify it

    std::map<std::string, std::vector<double>> data_map;

    // Iterate through column names to check if they are in keywords
    for (size_t col = 0; col < column_names.size(); ++col) {
        bool is_keyword = false;

        // Check if the current column name is in the keywords list
        for (size_t k = 0; k < keywords.size(); ++k) {
            if (column_names[col].find(keywords[k]) == 0) {  // Check if column name starts with a keyword
                is_keyword = true;
                break;
            }
        }

        // If the column is not in the keywords, add it to keywords
        if (!is_keyword) {
            std::cout << "Warning: New mineral found: " << column_names[col] << ". Adding it to keywords list." << std::endl;
            keywords.push_back(column_names[col]);  // Add the new column name to keywords
        }

        // Add column data to the map
        data_map[column_names[col]] = std::vector<double>(matrix.size());
        for (size_t row = 0; row < matrix.size(); ++row) {
            data_map[column_names[col]][row] = matrix[row][col];
        }
    }

    return data_map;
}

// Write the combined data map to an output file
void write_map_to_file(const std::string &output_file_name, const std::map<std::string, std::vector<double>> &data_map, const std::vector<std::string> &keywords) {
    std::ofstream output_file(output_file_name.c_str());
    if (!output_file.is_open()) {
        std::cerr << "Error: Could not open file " << output_file_name << std::endl;
        return;
    }

    // Write column headers in the order of the updated keywords
    for (size_t i = 0; i < keywords.size(); ++i) {
        if (data_map.find(keywords[i]) != data_map.end()) {
            output_file << keywords[i] << "\t";
        }
    }
    output_file << std::endl;

    // Write each row of data in the order of the keywords
    size_t row_count = data_map.begin()->second.size();  // All columns have the same number of rows
    for (size_t i = 0; i < row_count; ++i) {
        for (size_t j = 0; j < keywords.size(); ++j) {
            if (data_map.find(keywords[j]) != data_map.end()) {
                output_file << data_map.find(keywords[j])->second[i] << "\t";
            }
        }
        output_file << std::endl;
    }

    output_file.close();
    std::cout << "Data map successfully written to " << output_file_name << std::endl;
}

/*
int main() {
    std::string input_file_name, output_file_name;
    std::vector<std::vector<double>> matrix;   // Store data as a matrix
    std::vector<std::string> column_names;     // Store column names

    // Predefined keywords with the new additions (order matters)
    std::vector<std::string> keywords = {
        "node#", "P(bar)", "T(K)", "O", "Cpx", "Opx", "Gt", "C2/c", "Wad", "st", "Ring", 
        "ca-pv", "Aki", "Pv", "per", "CF", "wus"
    };

    // Get input file name
    std::cout << "Enter the name of the input file: ";
    std::getline(std::cin, input_file_name);

    // Read matrix data and column names
    read_matrix_from_perplex_tabfile(input_file_name, matrix, column_names);

    // Combine data and column names into a map, and update keywords if new minerals are found
    std::map<std::string, std::vector<double>> data_map = combine_data_to_map(matrix, column_names, keywords);

    // Get output file name
    std::cout << "Enter the name of the output file: ";
    std::getline(std::cin, output_file_name);

    // Write data map to output file
    write_map_to_file(output_file_name, data_map, keywords);

    return 0;
}
*/

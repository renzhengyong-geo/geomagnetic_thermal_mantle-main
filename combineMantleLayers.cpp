#include <iostream>
#include <vector>
#include <string>
#include <map>

// Function to combine three tables into one
void combineMantleLayers(
    const std::vector<std::vector<double>>& matrix_upper,
    const std::vector<std::string>& column_names_upper,
    const std::vector<std::vector<double>>& matrix_middle,
    const std::vector<std::string>& column_names_middle,
    const std::vector<std::vector<double>>& matrix_lower,
    const std::vector<std::string>& column_names_lower,
    std::vector<std::vector<double>>& combined_matrix,
    std::vector<std::string>& combined_column_names
) {
    // Step 1: Create a map to store all unique column names
    std::map<std::string, int> column_index_map;
    int col_index = 0;
    
    // Add column names from upper mantle
    for (size_t i = 0; i < column_names_upper.size(); ++i) {
        if (column_index_map.find(column_names_upper[i]) == column_index_map.end()) {
            column_index_map[column_names_upper[i]] = col_index++;
            combined_column_names.push_back(column_names_upper[i]);
        }
    }
    
    // Add column names from middle mantle
    for (size_t i = 0; i < column_names_middle.size(); ++i) {
        if (column_index_map.find(column_names_middle[i]) == column_index_map.end()) {
            column_index_map[column_names_middle[i]] = col_index++;
            combined_column_names.push_back(column_names_middle[i]);
        }
    }
    
    // Add column names from lower mantle
    for (size_t i = 0; i < column_names_lower.size(); ++i) {
        if (column_index_map.find(column_names_lower[i]) == column_index_map.end()) {
            column_index_map[column_names_lower[i]] = col_index++;
            combined_column_names.push_back(column_names_lower[i]);
        }
    }
    
    // Step 2: Create empty rows for the combined matrix based on the combined column size
    int total_columns = combined_column_names.size();
    
    // Helper function to add rows with the correct columns from an input table
    for (size_t row_idx = 0; row_idx < matrix_upper.size(); ++row_idx) {
        std::vector<double> new_row(total_columns, 0.0);  // Initialize all columns to 0
        for (size_t i = 0; i < column_names_upper.size(); ++i) {
            const std::string& col_name = column_names_upper[i];
            int col_idx = column_index_map[col_name];
            new_row[col_idx] = matrix_upper[row_idx][i];  // Assign the correct value to the correct column
        }
        combined_matrix.push_back(new_row);
    }

    // Fix for the middle mantle
    for (size_t row_idx = 0; row_idx < matrix_middle.size(); ++row_idx) {
        std::vector<double> new_row(total_columns, 0.0);  // Initialize all columns to 0
        for (size_t i = 0; i < column_names_middle.size(); ++i) {
            const std::string& col_name = column_names_middle[i];
            int col_idx = column_index_map[col_name];  // Map the middle mantle columns to the correct index in combined columns
            new_row[col_idx] = matrix_middle[row_idx][i];  // Assign the correct value to the correct column
        }
        combined_matrix.push_back(new_row);
    }

    // Fix for the lower mantle
    for (size_t row_idx = 0; row_idx < matrix_lower.size(); ++row_idx) {
        std::vector<double> new_row(total_columns, 0.0);  // Initialize all columns to 0
        for (size_t i = 0; i < column_names_lower.size(); ++i) {
            const std::string& col_name = column_names_lower[i];
            int col_idx = column_index_map[col_name];  // Map the lower mantle columns to the correct index in combined columns
            new_row[col_idx] = matrix_lower[row_idx][i];  // Assign the correct value to the correct column
        }
        combined_matrix.push_back(new_row);
    }
    
    return; 
}

/*
int main() {
    // Upper mantle data
    std::vector<std::vector<double>> matrix_upper = {
        {1, 7487.5, 1500, 60.2741, 20.1101, 19.2278, 0, 0, 0.387986},
        {2, 23468.8, 1500, 61.1773, 20.377, 9.79626, 8.64938, 0, 0},
        {3, 40562.5, 1500, 61.3091, 21.341, 5.88518, 11.4647, 0, 0},
        {4, 57593.8, 1500, 61.3275, 20.8162, 5.78332, 12.073, 0, 0},
        {5, 74712.5, 1500, 61.3347, 19.9194, 6.22192, 12.524, 0, 0},
        {6, 92337.5, 1500, 61.4053, 18.5893, 0, 13.2276, 6.77773, 0},
        {7, 110094, 1500, 61.441, 16.1966, 0, 14.4061, 7.95636, 0},
        {8, 128031, 1500, 61.524, 14.3515, 0, 17.5332, 6.59135, 0}
    };
    std::vector<std::string> column_names_upper = {"node#", "P(bar)", "T(K)", "O", "Cpx", "Opx", "Gt", "C2/c", "sp"};
    
    // Middle (transition) mantle data
    std::vector<std::vector<double>> matrix_middle = {
        {1, 146650, 1500, 12.0296, 27.3889, 60.5814, 0, 0, 0, 0},
        {2, 165600, 1500, 5.21022, 34.0364, 52.7618, 0, 7.99152, 0, 0},
        {3, 184850, 1500, 0, 32.3383, 0, 1.8475, 65.8142, 0, 0},
        {4, 204500, 1500, 0, 20.4313, 0, 4.10655, 72.5043, 2.95783, 0},
        {5, 224400, 1500, 0, 16.6955, 0, 0, 61.5704, 4.15406, 17.58}
    };
    std::vector<std::string> column_names_middle = {"node#", "P(bar)", "T(K)", "Cpx", "Gt", "Wad", "st", "Ring", "ca-pv", "Aki"};
    
    // Lower mantle data
    std::vector<std::vector<double>> matrix_lower = {
        {1, 275100, 1500, 6.19297, 72.2375, 16.754, 4.13968, 0.675885, 0},
        {2, 358800, 1500, 6.19201, 72.305, 16.6547, 4.23166, 0.616643, 0},
        {3, 444367, 1500, 6.19038, 72.35, 16.5932, 4.33429, 0.532224, 0},
        {4, 531700, 1500, 6.18835, 72.3807, 16.5771, 4.44423, 0.409659, 0},
        {5, 620900, 1500, 6.18609, 72.4039, 16.6139, 4.55861, 0.237458, 0},
        {6, 711867, 1500, 6.18374, 72.4213, 16.7093, 4.6801, 0.00556044, 0},
        {7, 804867, 1500, 6.18121, 72.3924, 16.605, 4.82149, 0, 0},
        {8, 900000, 1500, 6.17861, 72.3554, 16.5012, 4.96483, 0, 0},
        {9, 997533, 1500, 6.17597, 72.3108, 16.4014, 5.11191, 0, 0},
        {10, 1.09773e+06, 1500, 6.18241, 63.0279, 16.3109, 5.33468, 0, 9.14414},
        {11, 1.201e+06, 1500, 6.24847, 0, 15.2234, 5.88387, 1.2483, 71.3959},
        {12, 1.30773e+06, 1500, 6.25475, 0, 16.3371, 5.7002, 0, 71.708}
    };
    std::vector<std::string> column_names_lower = {"node#", "P(bar)", "T(K)", "ca-pv", "Pv", "per", "CF", "wus", "Ppv"};
    
    // Final combined matrix and column names
    std::vector<std::vector<double>> combined_matrix;
    std::vector<std::string> combined_column_names;

    // Combine the mantle layers
    combineMantleLayers(matrix_upper, column_names_upper, matrix_middle, column_names_middle, matrix_lower, column_names_lower, combined_matrix, combined_column_names);

    // Output the combined table (column names and values)
    for (size_t i = 0; i < combined_column_names.size(); ++i) {
        std::cout << combined_column_names[i] << "\t";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < combined_matrix.size(); ++i) {
        for (size_t j = 0; j < combined_matrix[i].size(); ++j) {
            std::cout << combined_matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    
    return 0;
}
*/



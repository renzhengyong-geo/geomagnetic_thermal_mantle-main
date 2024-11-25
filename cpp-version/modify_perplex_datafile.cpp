#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>  // For controlling precision and width

// Function to modify the values of the third column based on the given formula
void modify_perplex_datafile(const std::string &inputFileName, const std::string &outputFileName, double basalt_fraction) {
    // Ensure that basalt_fraction is within the range [0, 1]
    if (basalt_fraction < 0 || basalt_fraction > 1) {
        std::cerr << "Error: The value of basalt_fraction must be in the range [0, 1]. Provided value = " << basalt_fraction << std::endl;
        return;
    }

    // Define compositions for Basalt and Harzburgite
    std::map<std::string, double> basaltValues = {{"CAO", 13.05}, {"FEO", 7.68}, {"MGO", 10.49}, {"AL2O3", 16.08}, {"SIO2", 50.39}, {"NA2O", 1.87}};
    std::map<std::string, double> harzburgiteValues = {{"CAO", 0.5}, {"FEO", 7.83}, {"MGO", 46.36}, {"AL2O3", 0.65}, {"SIO2", 43.64}, {"NA2O", 0.01}};

    // Open input and output files
    std::ifstream inputFile(inputFileName);
    std::ofstream outputFile(outputFileName);

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error: Could not open input or output file!" << std::endl;
        std::cout << inputFileName <<" or \t"<<outputFileName<<"\n";
        return;
    }

    std::string line;
    bool inComponentList = false;  // Track if we are within the "thermodynamic component list" section

    // Read and process each line in the file
    while (getline(inputFile, line)) {
        // Detect the start and end of the component list section
        if (line.find("begin thermodynamic component list") != std::string::npos) {
            inComponentList = true;  // Set flag to true after this line
            outputFile << line << std::endl;  // Write this line unchanged
            continue;  // Skip processing this line further
        }

        if (line.find("end thermodynamic component list") != std::string::npos) {
            inComponentList = false;  // Set flag to false after this line
            outputFile << line << std::endl;  // Write this line unchanged
            continue;  // Skip processing this line further
        }

        // Only modify lines inside the component list section
        if (inComponentList) {
            std::stringstream ss(line);
            std::string component;
            std::vector<std::string> columns;

            // Split the line into columns
            while (ss >> component) columns.push_back(component);

            // If the component exists in both maps, modify the third column value
            if (columns.size() > 2 && basaltValues.count(columns[0]) && harzburgiteValues.count(columns[0])) {
                // Calculate new value: (Composition of Basalt × basalt_fraction) + ((1 − basalt_fraction) × Composition of Harzburgite)
                double newValue = basaltValues[columns[0]] * basalt_fraction + harzburgiteValues[columns[0]] * (1 - basalt_fraction);

                // Convert the new double value to string with fixed precision
                std::ostringstream out;
                out << std::fixed << std::setw(12) << std::setprecision(6) << newValue;
                columns[2] = out.str();  // Update the third column with new value
            }

            // Reconstruct the line with modified columns and correct spacing
            outputFile <<                  columns[0] << "  " << std::setw(2) << columns[1] << "  "
                       << std::setw(12) << columns[2] << "  " << std::setw(12) << columns[3] << "  "
                       << std::setw(12) << columns[4] << "  " << columns[5] << "  amount" << std::endl;
        } else {
            // Write lines outside the component list section unchanged
            outputFile << line << std::endl;
        }
    }

    std::cout << "File processing completed. Check the output file: " << outputFileName << std::endl;
}

/*
int main() {
    std::string inputFileName = "khan2016pyrolite.dat";  // Input file name
    std::string outputFileName = "modified_khan2016pyrolite.dat";  // Output file name
    double basalt_fraction;

    // Prompt the user to enter the value of basalt_fraction
    std::cout << "Enter the value of basalt_fraction (0 <= basalt_fraction <= 1): ";
    std::cin >> basalt_fraction;

    // Call the function with the specified parameters
    modify_perplex_datafile(inputFileName, outputFileName, basalt_fraction);

    return 0;
}
*/

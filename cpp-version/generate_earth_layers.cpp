#include <iostream>
#include <vector>
#include <iomanip> // For std::setw


#include "thermal_mantle.h"

// Function to generate Earth's layer structure based on the number of layers in each region
std::vector<Layer> generateEarthLayers(int n_um, int n_tz, int n_lm) {
    std::vector<Layer> layers;

    // 1. Surface to 410 km
    double um_layer_height = 410.0 / n_um;
    for (int i = 0; i < n_um; ++i) {
        double upper = i * um_layer_height;
        double lower = (i + 1) * um_layer_height;
        layers.push_back({upper, lower});
    }

    // 2. 410 km to 660 km
    double tz_layer_height = (660.0 - 410.0) / n_tz;
    for (int i = 0; i < n_tz; ++i) {
        double upper = 410.0 + i * tz_layer_height;
        double lower = 410.0 + (i + 1) * tz_layer_height;
        layers.push_back({upper, lower});
    }

    // 3. 660 km to 2900 km
    double lm_layer_height = (2900.0 - 660.0) / n_lm;
    for (int i = 0; i < n_lm; ++i) {
        double upper = 660.0 + i * lm_layer_height;
        double lower = 660.0 + (i + 1) * lm_layer_height;
        layers.push_back({upper, lower});
    }

    return layers;
}

// Function to print the Earth's layers in a formatted way
void printEarthLayers(const std::vector<Layer>& layers, int n_um, int n_tz, int n_lm) {
    std::cout << std::fixed << std::setprecision(2); // Set precision to 2 decimal places

    std::cout << "------------------------------------------" << std::endl;
    std::cout << "|  Layer No.  |  Section Name    |  Bounds (km)           |" << std::endl;
    std::cout << "------------------------------------------" << std::endl;

    // Upper Mantle (0 - 410 km)
    std::cout << "|  Upper Mantle (0 km to 410 km)                             |" << std::endl;
    for (int i = 0; i < n_um; ++i) {
        std::cout << "|     " << std::setw(3) << i + 1 << "      |  Upper Mantle  |  "
                  << std::setw(6) << layers[i].upper_bound << " - " << std::setw(6) << layers[i].lower_bound << "    |" << std::endl;
    }

    // Transition Zone (410 - 660 km)
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "|  Transition Zone (410 km to 660 km)                         |" << std::endl;
    for (int i = n_um; i < n_um + n_tz; ++i) {
        std::cout << "|     " << std::setw(3) << i + 1 << "      |  Transition Z. |  "
                  << std::setw(6) << layers[i].upper_bound << " - " << std::setw(6) << layers[i].lower_bound << "    |" << std::endl;
    }

    // Lower Mantle (660 - 2900 km)
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "|  Lower Mantle (660 km to 2900 km)                           |" << std::endl;
    for (int i = n_um + n_tz; i < n_um + n_tz + n_lm; ++i) {
        std::cout << "|     " << std::setw(3) << i + 1 << "      |  Lower Mantle  |  "
                  << std::setw(6) << layers[i].upper_bound << " - " << std::setw(6) << layers[i].lower_bound << "    |" << std::endl;
    }

    std::cout << "------------------------------------------" << std::endl;
}

/*
int main() {
    // Define the number of layers for each region
    int n_um = 8;  // Number of layers in the upper mantle (0 km to 410 km)
    int n_tz = 5;  // Number of layers in the transition zone (410 km to 660 km)
    int n_lm = 12; // Number of layers in the lower mantle (660 km to 2900 km)

    // Generate Earth's layer structure
    std::vector<Layer> earth_layers = generateEarthLayers(n_um, n_tz, n_lm);

    // Print the Earth's layer structure in a formatted way
    printEarthLayers(earth_layers, n_um, n_tz, n_lm);

    return 0;
}
*/

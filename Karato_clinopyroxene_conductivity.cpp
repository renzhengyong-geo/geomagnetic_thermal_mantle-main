#include <cmath>
#include <stdexcept>
#include <iostream>

// Function to calculate conductivity
double Karato_orthopyroxene_conductivity(double temperature, double water_content, double pressure) {
    // Constants
    const double R = 8.314;       // Gas constant, J/(K⋅mol)
    const double A1 = std::pow(10, 2.7);   // Pre-exponential factor for polaron conduction, S/m
    const double A2 = std::pow(10, 2.6);   // Pre-exponential factor for proton conduction, S/m
    const double r1 = 0.0;        // Empirical constant for polaron conduction
    const double r2 = 0.62;       // Empirical constant for proton conduction
    const double E1 = 147000;     // Activation energy for polaron conduction, J/mol
    const double E2 = 82000;      // Activation energy for proton conduction, J/mol
    const double V1 = 0.0;        // Activation volume for polaron conduction, cc/mol
    const double V2 = 0.0;        // Activation volume for proton conduction, cc/mol

    // Check if temperature is greater than zero
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be greater than zero Kelvin.");
    }

    // Check if pressure is non-negative
    if (pressure < 0) {
        throw std::invalid_argument("Pressure must be non-negative.");
    }

    // Convert pressure from Pa to J/cc (1 Pa = 1 J/m³ = 1e-6 J/cc)
    pressure *= 1e3;

    // Check if water_content is within the range of 0 to 1
    if (water_content < 0 || water_content > 1) {
        throw std::invalid_argument("Water content must be in the range of 0 to 1.");
    }

    // Calculate conductivities
    double sigma_polaron = A1 * std::pow(water_content, r1) * std::exp(-(E1 + pressure * V1) / (R * temperature));
    double sigma_proton =  A2 * std::pow(water_content, r2) * std::exp(-(E2 + pressure * V2) / (R * temperature));

    // Sum the conductivities
    double sigma = sigma_polaron + sigma_proton;

    return sigma;
}

// Example usage
int main() {
    double temperature = 1000.0;  // Example temperature in Kelvin
    double water_content = 0.1;   // Example water content
    double pressure = 1e5;        // Example pressure in Pascal

    try {
        double conductivity = Karato_orthopyroxene_conductivity(temperature, water_content, pressure);
        std::cout << "Conductivity: " << conductivity << " S/m" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

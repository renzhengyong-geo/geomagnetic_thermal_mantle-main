#include <cmath>
#include <stdexcept>
#include <iostream>

double Yoshino_ca_pv_conductivity(double temperature, double water_content, double pressure) {
    // Constants
    const double k_boltzmann = 1.380649e-23;  // Boltzmann constant in J/K
    const double J2eV = 1.0 / 1.60217733e-19; // 1 eV = 1.60217733e-19 J
    const double k = k_boltzmann * J2eV;      // Boltzmann constant in eV/K
    const double sigma_i = 631;               // Intrinsic conductivity in S/m
    const double Hi = 1.03;                   // Activation enthalpy for intrinsic conduction in eV
    const double sigma_h = 0;                 // Hydrogen-related conductivity in S/m
    const double Hh = 0;                      // Activation enthalpy for hydrogen in eV
    const double sigma_p = 0;                 // Proton-related conductivity in S/m
    const double Hp = 0;                      // Activation enthalpy for proton conduction in eV
    const double alpha = 0;                   // Constant for proton conduction

    // Input validation
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be greater than zero Kelvin.");
    }

    if (pressure < 0) {
        throw std::invalid_argument("Pressure must be non-negative.");
    }

    if (water_content < 0 || water_content > 1) {
        throw std::invalid_argument("Water content must be in the range of 0 to 1.");
    }

    // Convert pressure from Pa to J/cc (1 Pa = 1 J/m³ = 1e-6 J/cc)
    pressure *= 1e3;

    // Sum the conductivities
    double sigma = sigma_i * std::exp(-Hi / (k * temperature)) +
                   sigma_h * std::exp(-Hh / (k * temperature)) +
                   sigma_p * water_content * std::exp(-(Hp - alpha * std::pow(water_content, 1.0 / 3.0)) / (k * temperature));

    return sigma;
}

// Example usage
int main() {
    double temperature = 1000.0;  // Example temperature in Kelvin
    double water_content = 0.1;   // Example water content
    double pressure = 1e5;        // Example pressure in Pascal

    try {
        double conductivity = Yoshino_ca_pv_conductivity(temperature, water_content, pressure);
        std::cout << "Conductivity: " << conductivity << " S/m" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

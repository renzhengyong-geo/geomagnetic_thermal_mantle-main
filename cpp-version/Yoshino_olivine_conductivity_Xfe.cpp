#include <cmath>
#include <stdexcept>
#include <iostream>

double Yoshino_olivine_conductivity_Xfe(double temperature, double water_content, double pressure, double Xfe) {
    // Constants
    const double k_boltzmann = 1.380649e-23;  // Boltzmann constant in J/K
    const double J2eV = 1.0 / 1.60217733e-19; // 1 eV = 1.60217733e-19 J
    const double k = k_boltzmann * J2eV;      // Boltzmann constant in eV/K
    const double sigma_i = std::pow(10, 4.73); // Intrinsic conductivity in S/m
    const double Hi = 2.31;                    // Activation enthalpy for intrinsic conduction in eV
    const double sigma_h = 525;                // Hydrogen-related conductivity in S/m
    const double Hh = 1.96;                    // Activation enthalpy for hydrogen in eV
    const double alpha_h = 1.49;               // Constant for hydrogen conduction
    const double V0 = -0.01;                   // Activation volume
    const double belta_h = 0.22;               // Beta constant for hydrogen
    const double sigma_p = std::pow(10, 1.9);  // Proton-related conductivity in S/m
    const double Hp = 0.92;                    // Activation enthalpy for proton conduction in eV
    const double alpha = 0.16;                 // Constant for proton conduction

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

    if (Xfe < 0 || Xfe > 1) {
        throw std::invalid_argument("Iron content (Xfe) must be in the range of 0 to 1.");
    }

    // Sum the conductivities
    double sigma = sigma_i * std::exp(-Hi / (k * temperature)) +
                   sigma_h * Xfe * std::exp(-((Hh - alpha_h * std::pow(Xfe, 1.0 / 3.0) + 0.0104 * pressure * (V0 - belta_h * Xfe))) / (k * temperature)) +
                   sigma_p * water_content * std::exp(-(Hp - alpha * std::pow(water_content, 1.0 / 3.0)) / (k * temperature));

    return sigma;
}

// Example usage
int main() {
    double temperature = 1000.0;  // Example temperature in Kelvin
    double water_content = 0.1;   // Example water content
    double pressure = 1e5;        // Example pressure in Pascal
    double Xfe = 0.1;             // Example iron content fraction

    try {
        double conductivity = Yoshino_olivine_conductivity_Xfe(temperature, water_content, pressure, Xfe);
        std::cout << "Conductivity: " << conductivity << " S/m" << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

#include <iostream>
#include <cmath>
#include <cstdio>  // For printf

// Constants
const double k = 1.380649e-23;  // Boltzmann constant, J/K

// Function to compute the conductivity of Wadsleyite using Yoshino's model
double Yoshino_wadsleyite_conductivity(double T, double Cw, double Xfe) {
    // Constants for the equation
    const double sigma_0i = 0.0;              // S/m (Pre-exponential factor for ionic conduction)
    const double delta_Hi = 0.0;              // J (Converted from eV to J)
    const double sigma_0h = pow(10, 2.46);    // S/m (Pre-exponential factor for hopping conduction)
    const double delta_Hh = 1.45 * 1.60218e-19;  // J (Converted from eV to J)
    const double sigma_0p = pow(10, 1.40);    // S/m (Pre-exponential factor for proton conduction)
    const double delta_H0_p = 0.83 * 1.60218e-19; // J (Converted from eV to J)
    const double alpha_p = 0.20 * 1.60218e-19;   // J (Converted from eV to J)

    // Check if temperature is greater than zero
    if (T <= 0) {
        printf("Error: Temperature must be greater than zero Kelvin.\n");
        std::exit(1);
    }

    // Check if water content is within the range of 0 to 1
    if (Cw < 0 || Cw > 1) {
        printf("Error: Water content must be in the range of 0 to 1.\n");
        std::exit(1);
    }

    // Check if iron content is within the range of 0 to 1
    if (Xfe < 0 || Xfe > 1) {
        printf("Error: Iron content must be in the range of 0 to 1.\n");
        std::exit(1);
    }

    // Sum the conductivities
    double sigma = sigma_0i * exp(-delta_Hi / (k * T)) +
                   sigma_0h * Xfe * exp(-delta_Hh / (k * T)) +
                   sigma_0p * Cw * exp((-delta_H0_p + alpha_p * pow(Cw, 1.0 / 3.0)) / (k * T));

    return sigma;
}

int main() {
    // Example inputs
    double T = 1200;      // Temperature in Kelvin
    double Cw = 0.01;     // Water content (dimensionless)
    double Xfe = 0.1;     // Iron content (dimensionless)

    // Compute conductivity
    double sigma = Yoshino_wadsleyite_conductivity(T, Cw, Xfe);

    // Output the result using printf
    printf("Electrical conductivity: %e S/m\n", sigma);

    return 0;
}

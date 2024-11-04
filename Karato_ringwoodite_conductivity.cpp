#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstdio>  // For printf

// Constants
const double R = 8.314;  // Gas constant, J/(K⋅mol)

// Function to compute the conductivity of Ringwoodite using Karato's model
double Karato_ringwoodite_conductivity(double T, double Cw, double P) {
    // Constants for the equation
    const double A1 = 0.0;              // Pre-exponential factor for polaron conduction, S/m
    const double A2 = pow(10, 3.6);     // Pre-exponential factor for proton conduction, S/m
    const double r1 = 0;                // Empirical constant for polaron conduction
    const double r2 = 0.69;             // Empirical constant for proton conduction
    const double E1 = 0;                // Activation energy for polaron conduction, J/mol
    const double E2 = 104000;           // Activation energy for proton conduction, J/mol
    const double V1 = 0.0;              // Activation volume for polaron conduction, cc/mol
    const double V2 = 0.0;              // Activation volume for proton conduction, cc/mol

    // Check if temperature is greater than zero
    if (T <= 0) {
        printf("Error: Temperature must be greater than zero Kelvin.\n");
        std::exit(1);
    }

    // Check if pressure is non-negative
    if (P < 0) {
        printf("Error: Pressure must be non-negative.\n");
        std::exit(1);
    }

    // Convert pressure from Pa to J/cc
    P = 1e-6 * P;  // 1 Pa = 1 J/m³ = 1e-6 J/cc

    // Check if water content is within the range of 0 to 1
    if (Cw < 0 || Cw > 1) {
        printf("Error: Water content must be in the range of 0 to 1.\n");
        std::exit(1);
    }

    // Calculate polaron conductivity
    double sigma_polaron = A1 * pow(Cw, r1) * exp(-(E1 + P * V1) / (R * T));

    // Calculate proton conductivity
    double sigma_proton = A2 * pow(Cw, r2) * exp(-(E2 + P * V2) / (R * T));

    // Sum the conductivities
    return sigma_polaron + sigma_proton;
}

int main() {
    // Example inputs
    double T = 1200;      // Temperature in Kelvin
    double Cw = 0.01;     // Water content (dimensionless)
    double P = 1e9;       // Pressure in Pascal

    // Compute conductivity
    double sigma = Karato_ringwoodite_conductivity(T, Cw, P);

    // Output the result using printf
    printf("Electrical conductivity: %e S/m\n", sigma);

    return 0;
}

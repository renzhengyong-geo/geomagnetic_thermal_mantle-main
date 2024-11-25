#include <iostream>
#include <cmath>
#include <stdexcept>
#include <cstdio>  // For printf

/* ------------------------------------------------------------------------
% FUNCTION: Karato_Olivine_conductivity
% ------------------------------------------------------------------------
% DESCRIPTION:
% Compute the conductivity of olivine using Karato's group high temperature
% and pressure experiments. There are two conductive mechanisms which are
% the hopping of electrons between ferric and ferrous iron (Fe3+->Fe2+,
% also referred to polaron conduction) or the migration of protons (H+,OH-)
%
% FORMULA:
% The conductivity is computed using the formula:
%   sigma = A1*Cw^{r1}*exp(-(E1+PV1)/(RT)) + A2*Cw^{r2}*exp(-(E2+PV2)/(RT))
% where:
%   - Subscript 1 stands for polaron conduction.
%   - Subscript 2 stands for proton conduction (e.g., H).
%   - A1, A2 are the pre-exponential terms which already contain the influence
%     of the Mg# and oxygen fugacity (A/s)
%   - Cw is the water content.
%   - r1, r2 are constants that depend on the mechanism of electrical conduction.
%   - E1, E2 are the activation energies.
%   - V1, V2 are the activation volumes.
%   - T is absolute temperature (in Kelvin)
%   - P is the pressure (in Pascal)
%
% REFERENCES:
% Values of A, r, E, V are from:
%   Karato, S. I. (2011). Water distribution across the mantle transition zone
%   and its implications for global material circulation. Earth and Planetary
%   Science Letters, 301(3–4), 413–423. https://doi.org/10.1016/j.epsl.2010.11.038
%
%
% INPUTS:
%   T   - Temperature in Kelvin.
%   Cw  - Water content (dimensionless).
%   P   - Pressure in Pascal.
%
% OUTPUT:
%   sigma - Electrical conductivity (S/m).
% ------------------------------------------------------------------------
*/

// Function to compute the conductivity of olivine using Karato's model
double Karato_olivine_conductivity(double T, double Cw, double P) {
    // temperature: Kelvin (K); 0 K is equivalent to -273.15°C.
    // water_content: dimensionless, typically in the range 0-1.
    // pressure: Pascal (Pa); 1 Pa = 1 J/m³ = 1e-6 J/cc.
    const double R = 8.314;           // Gas constant, J/(K⋅mol)
    const double A1 = pow(10, 2.4);   // Pre-exponential factor for polaron conduction, S/m
    const double A2 = pow(10, 3.1);   // Pre-exponential factor for proton conduction, S/m
    const double r1 = 0;              // Empirical constant for polaron conduction
    const double r2 = 0.62;           // Empirical constant for proton conduction
    const double E1 = 154000;         // Activation energy for polaron conduction, J/mol
    const double E2 = 87000;          // Activation energy for proton conduction, J/mol
    const double V1 = 2.4;            // Activation volume for polaron conduction, cc/mol
    const double V2 = 0.0;            // Activation volume for proton conduction, cc/mol

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
    double sigma_proton =  A2 * pow(Cw, r2) * exp(-(E2 + P * V2) / (R * T));

    // Sum the conductivities
    return sigma_polaron + sigma_proton;
}

int main() {
    // Example inputs
    double T = 1200;      // Temperature in Kelvin
    double Cw = 0.01;     // Water content (dimensionless)
    double P = 1e9;       // Pressure in Pascal

    // Compute conductivity
    double sigma = Karato_olivine_conductivity(T, Cw, P);

    // Output the result using printf
    printf("Electrical conductivity: %e S/m\n", sigma);

    return 0;
}

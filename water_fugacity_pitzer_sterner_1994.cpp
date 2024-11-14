#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

// Coefficients for the Pitzer and Sterner EOS
std::vector<std::vector<double>> coeff = {
    {0, 0, 0.24657688e6, 0.51359951e2, 0, 0},
    {0, 0, 0.58638965e0, -0.28646939e-2, 0.31375577e-4, 0},
    {0, 0, -0.62783840e1, 0.14791599e-1, 0.35779579e-3, 0.15432925e-7},
    {0, 0, 0, -0.42719875e0, -0.16325155e-4, 0},
    {0, 0, 0.56654978e4, -0.16580167e2, 0.76560762e-1, 0},
    {0, 0, 0, 0.10917883e0, 0, 0},
    {0.38878656e13, -0.13494878e9, 0.30916564e6, 0.75591105e1, 0, 0},
    {0, 0, -0.65537898e5, 0.18810675e3, 0, 0},
    {-0.14182435e14, 0.18165390e9, -0.19769068e6, -0.23530318e2, 0, 0},
    {0, 0, 0.92093375e5, 0.12246777e3, 0, 0}
};

struct PSeos_params {
    double temperature_C;
    double targetP_GPa;
};

// Function to calculate pressure using EOS
double PSeos(double volume, void *params) {
    PSeos_params *p = (PSeos_params *)params;
    double temperature_C = p->temperature_C;
    const double targetP_GPa = p->targetP_GPa;

    // Convert target pressure from GPa to bars
    const double targetP = targetP_GPa * 1e4; // GPa to bars (1 GPa = 10,000 bars)
    // Convert temperature from Celsius to Kelvin
    const double temperature = temperature_C + 273.15;
    const double R = 8314510.0; // Pa*cc/(K*mol)

    // Calculate density
    const double den = 1.0 / volume; // mol/cc

    // Calculate temperature-dependent coefficients
    std::vector<double> c(10);
    for (int i = 0; i < 10; ++i) {
        c[i] = coeff[i][0] * pow(temperature, -4) + coeff[i][1] * pow(temperature, -2) +
               coeff[i][2] * pow(temperature, -1) + coeff[i][3] +
               coeff[i][4] * temperature + coeff[i][5] * pow(temperature, 2);
    }

    // Calculate pressure based on Pitzer and Sterner EOS
    const double numerator = (c[2] + 2 * c[3] * den + 3 * c[4] * pow(den, 2) + 4 * c[5] * pow(den, 3));
    const double denominator = pow((c[1] + c[2] * den + c[3] * pow(den, 2) + c[4] * pow(den, 3) + c[5] * pow(den, 4)), 2);
    const double pressure = (den + c[0] * pow(den, 2) - pow(den, 2) * (numerator / denominator) +
                       c[6] * pow(den, 2) * exp(-c[7] * den) + c[8] * pow(den, 2) * exp(-c[9] * den)) *
                      (R * temperature) / 1e5; // Convert to bars

    return pressure - targetP; // Difference used for root finding
}

// Function to calculate molar volume for a given pressure and temperature
double PSvolume(double pressure_GPa, double temperature_C) {
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    const double x_lo = 1.0, x_hi = 100.0;
    gsl_function F;
    PSeos_params params = {temperature_C, pressure_GPa};

    F.function = &PSeos;
    F.params = &params;

    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    int status;
    int iter = 0, max_iter = 100;
    double r = 0;
    double x_lo_current = x_lo;
    double x_hi_current = x_hi;

    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        x_lo_current = gsl_root_fsolver_x_lower(s);
        x_hi_current = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo_current, x_hi_current, 0, 1e-8);
    } while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free(s);

    if (status != GSL_SUCCESS) {
        throw std::runtime_error("Accurate molar volume not found: pressure difference exceeds tolerance.");
    }

    return r;
}

// Function to calculate fugacity for a given pressure and temperature
double PSfugacity(double pressure_GPa, double temperature_C) {
    // Calculate the molar volume
    double volume = PSvolume(pressure_GPa, temperature_C);

    // Convert temperature from Celsius to Kelvin
    const double temperature = temperature_C + 273.15;
    const double R = 8314510.0; // Pa*cc/(K*mol)

    // Convert pressure from GPa to bars
    const double pressure_bar = pressure_GPa * 1e4; // 1 GPa = 10,000 bars

    // Calculate density
    const double den = 1.0 / volume; // mol/cc

    // Calculate temperature-dependent coefficients
    std::vector<double> c(10);
    for (int i = 0; i < 10; ++i) {
        c[i] = coeff[i][0] * pow(temperature, -4) + coeff[i][1] * pow(temperature, -2) +
               coeff[i][2] * pow(temperature, -1) + coeff[i][3] +
               coeff[i][4] * temperature + coeff[i][5] * pow(temperature, 2);
    }

    // Fugacity calculation
    // Adjust for correct unit handling: pressure in bars, R in Pa*cc/(K*mol)
    const double fug = exp(
        log(den) + c[0] * den +
        (1 / (c[1] + c[2] * den + c[3] * pow(den, 2) + c[4] * pow(den, 3) + c[5] * pow(den, 4)) - 1 / c[1]) -
        c[6] / c[7] * (exp(-c[7] * den) - 1) -
        c[8] / c[9] * (exp(-c[9] * den) - 1) +
        (pressure_bar * 1e5) / (den * R * temperature) + // Convert pressure from bars to Pascals for consistency with R
        log(R * temperature) - 1
    ) / 1e5; // Convert from Pa to bar

    // Convert fugacity from bars to GPa (since output should be in GPa)
    return fug * 1e-4; // 1 bar = 1e-4 GPa
}


int main() {
    try {
        double pressure_GPa = 8; // Pressure in GPa
        double temperature_C = 3000; // Temperature in Celsius

        double volume = PSvolume(pressure_GPa, temperature_C);
        double fugacity = PSfugacity(pressure_GPa, temperature_C);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Molar Volume: " << volume << " cc/mol" << std::endl;
        std::cout << "Fugacity: " << fugacity << " GPa" << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

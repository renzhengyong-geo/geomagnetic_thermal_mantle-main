#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <limits>

using namespace std;

// Function to compute the Hashin-Shtrikman bounds of electrical conductivity
void compute_bounds(const vector<double>& x_i, const vector<double>& sigma_i, double& sigma_hs_min, double& sigma_hs_max) {
    assert(x_i.size() == sigma_i.size());  // Ensure both vectors are of the same size
    int N = x_i.size();  // Number of mineral phases
    
    // Find the minimum and maximum conductivities
    double sigma_min = sigma_i[0];
    double sigma_max = sigma_i[0];
    
    for (int i = 1; i < N; ++i) {
        if (sigma_i[i] < sigma_min) {
            sigma_min = sigma_i[i];
        }
        if (sigma_i[i] > sigma_max) {
            sigma_max = sigma_i[i];
        }
    }

    // Compute the lower bound using sigma_- (minimum conductivity)
    double temp_min = 0.0;
    for (int i = 0; i < N; ++i) {
        temp_min += x_i[i] / (sigma_i[i] + 2 * sigma_min);
    }
    sigma_hs_min = 1.0 / temp_min - 2 * sigma_min;

    // Compute the upper bound using sigma_+ (maximum conductivity)
    double temp_max = 0.0;
    for (int i = 0; i < N; ++i) {
        temp_max += x_i[i] / (sigma_i[i] + 2 * sigma_max);
    }
    sigma_hs_max = 1.0 / temp_max - 2 * sigma_max;
}

// Function to compute the equation f(sigma_sc) based on the given formula
double f(double sigma_sc, const vector<double>& x_i, const vector<double>& sigma_i) {
    int N = x_i.size();
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += x_i[i] * (sigma_i[i] - sigma_sc) / (sigma_i[i] + 2 * sigma_sc);
    }
    return sum;
}

// Function to compute the self-consistent conductivity sigma_sc using the Bisection method
double compute_sigma_sc(const vector<double>& x_i, const vector<double>& sigma_i) {
    // Compute the Hashin-Shtrikman bounds
    double sigma_hs_min, sigma_hs_max;
    compute_bounds(x_i, sigma_i, sigma_hs_min, sigma_hs_max);

    // Bisection method for root-finding
    double lower_bound = sigma_hs_min;  // Lower bound for sigma_sc
    double upper_bound = sigma_hs_max;  // Upper bound for sigma_sc
    double tolerance = 1e-8;           // Desired accuracy for sigma_sc
    double sigma_sc = lower_bound;

    // Ensure f(lower_bound) and f(upper_bound) have opposite signs
    double f_lower = f(lower_bound, x_i, sigma_i);
    double f_upper = f(upper_bound, x_i, sigma_i);
    
    // Check if the function at the bounds have opposite signs
    if (f_lower * f_upper > 0) {
        cout << "Error: Function values at the bounds have the same sign. Cannot apply Bisection method." << endl;
        return -9999;  // Use a specific value to indicate failure
    }

    // Bisection method to find the root
    while (upper_bound - lower_bound > tolerance) {
        sigma_sc = (lower_bound + upper_bound) / 2.0;  // Midpoint
        double f_val = f(sigma_sc, x_i, sigma_i);  // Calculate f(sigma_sc)

        // Update bounds based on the sign of f(sigma_sc)
        if (f_val > 0) {
            lower_bound = sigma_sc;  // Root is between sigma_sc and upper_bound
        } else {
            upper_bound = sigma_sc;  // Root is between lower_bound and sigma_sc
        }
    }

    return sigma_sc;  // Return the computed sigma_sc
}

int main() {
    // Example parameters (replace with your actual values)
    vector<double> x_i = {0.3, 0.5, 0.2};      // Volume fractions of each mineral
    vector<double> sigma_i = {1.5, 2.0, 3.5};  // Conductivities of each mineral
    
    // Compute self-consistent conductivity
    double sigma_sc = compute_sigma_sc(x_i, sigma_i);
    
    // Display the result
    if (sigma_sc != -9999) {  // Check for failure
        cout << "The self-consistent conductivity sigma_sc is: " << sigma_sc << endl;
    } else {
        cout << "Error: Unable to compute self-consistent conductivity." << endl;
    }

    return 0;
}

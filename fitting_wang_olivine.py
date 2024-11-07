import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt

# Data
Cw = np.array([0.027, 0.027, 0.080, 0.060, 0.019, 0.019, 0.019, 0.019, 0.013, 0.013, 0.013, 0.013, 0.010, 0.010, 0.010, 0.010, 0.010])
T = np.array([960, 1273, 881, 873, 1273, 1073, 973, 873, 1286, 1171, 1068, 963, 1273, 1173, 1073, 973, 873])
s = np.array([0.0025, 0.0386, 0.0010, 0.0009, 0.0273, 0.0032, 0.0014, 0.0004, 0.0133, 0.0055, 0.0026, 0.0007, 0.0108, 0.0060, 0.0036, 0.0017, 0.0004])

# Compute ln(s) and ln(Cw)
ln_s = np.log(s)
ln_Cw = np.log(Cw)

# Compute 1/T
inverse_T = 1 / T

# MCMC Model using PyMC3
with pm.Model() as model:
    # Priors for A, r, H
    A = pm.Normal('A', mu=1e-3, sigma=1, shape=1)  # Normal prior for A (S/m)
    r = pm.Normal('r', mu=1, sigma=0.5, shape=1)    # Normal prior for r
    H = pm.Normal('H', mu=30, sigma=10, shape=1)    # Normal prior for H (kJ/mol)
    
    # Expected value based on the model
    expected_ln_s = pm.math.log(A) + r * ln_Cw - H * inverse_T
    
    # Likelihood (assuming Gaussian noise)
    sigma_obs = pm.HalfNormal('sigma_obs', sigma=0.1)
    obs = pm.Normal('obs', mu=expected_ln_s, sigma=sigma_obs, observed=ln_s)

    # Sampling from the posterior using MCMC
    trace = pm.sample(2000, return_inferencedata=False, tune=1000)

# Results
pm.traceplot(trace)
plt.show()

# Get the posterior samples
A_posterior = trace['A']
r_posterior = trace['r']
H_posterior = trace['H']

# Calculate the mean of the posterior samples
A_mean = np.mean(A_posterior)
r_mean = np.mean(r_posterior)
H_mean = np.mean(H_posterior)

print(f"Estimated A: {A_mean:.4e} S/m")
print(f"Estimated r: {r_mean:.4f}")
print(f"Estimated H: {H_mean:.4f} kJ/mol")


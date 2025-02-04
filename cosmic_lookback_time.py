import numpy as np
from scipy.integrate import quad

# Define constants
f_Mpc_to_km = 3.086e19  # Conversion factor from Mpc to km
H0 = 70  # Hubble constant in km/s/Mpc
Omega_m_0 = 0.3  # Matter density parameter
Omega_r_0 = 2e-4  # Radiation density parameter
Omega_Lambda_0 = 0.7  # Dark energy density parameter
w_Lambda = -1  # Equation of state parameter for dark energy
z_infinity = 1.6e5  # Upper limit for redshift
f_yr_to_sec = 365 * 24 * 60 * 60  # Seconds per year

# Hubble time in seconds (inverse of Hubble constant)
Hubble_time_sec = f_Mpc_to_km / H0

# Define the integrand for the lookback time
def integrand(z_prime):
    return 1 / ((1 + z_prime) * np.sqrt(
        Omega_m_0 * (1 + z_prime)**3 +
        Omega_r_0 * (1 + z_prime)**4 +
        Omega_Lambda_0 * (1 + z_prime)**(3 * (1 + w_Lambda))
    ))

# Perform the integration
t_clb_s, _ = quad(integrand, 0, z_infinity)

print('# Lookback time in seconds')
t_clb_s *= Hubble_time_sec

print('# Lookback time in years')
t_clb_yr = t_clb_s / f_yr_to_sec

print(t_clb_yr)

print('# Lookback time in billion years')
t_clb_Gyr = t_clb_s / f_yr_to_sec / 1e9

print(t_clb_Gyr)
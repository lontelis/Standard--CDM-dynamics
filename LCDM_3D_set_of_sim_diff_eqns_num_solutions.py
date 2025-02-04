"""
Author: Pierros Ntelis, 14 June 2024
"""

########## Solving and plotting the solution functions of time of the system of differential equations ###########

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the system of differential equations allowing for complex numbers
def system(N, y):
    x, y_comp, z = y  # Note: y_comp is used to distinguish from the input y
    dx_dN = -3 * x + x * (3 + y_comp - 3 * z)
    dy_dN = -4 * y_comp + y_comp * (3 + y_comp - 3 * z)
    dz_dN = z * (3 + y_comp - 3 * z)
    return [dx_dN, dy_dN, dz_dN]

# Initial conditions (allowing complex numbers if necessary)
N0 = -13
x0 = 0.001 + 2.0j  # Initial condition for x(N)
y0 = 0.998 + 2.0j  # Initial condition for y(N)
z0 = 0.001 + 2.0j  # Initial condition for z(N)

# Set the initial values and range for N
initial_conditions = [x0, y0, z0]
Nfinal = 2  # Define the final value of N
N_span = [N0, Nfinal]  # The range of N

# Solve the differential equations numerically allowing for complex values
solution = solve_ivp(system, N_span, initial_conditions, method='RK45', t_eval=np.linspace(N0, Nfinal, 400))

# Extract solutions for plotting
N_values = solution.t
x_values = solution.y[0]
y_values = solution.y[1]
z_values = solution.y[2]

# Plotting
plt.ion()

# Plotting real parts
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(N_values, x_values.real, label='Re[x(N)]', color='b')
plt.plot(N_values, y_values.real, label='Re[y(N)]', color='r')
plt.plot(N_values, z_values.real, label='Re[z(N)]', color='g')
plt.xlabel('N')
plt.ylabel('Real Parts')
plt.title('Real Parts of x(N), y(N), and z(N)')
plt.legend()
plt.grid()

# Plotting imaginary parts
plt.subplot(1, 2, 2)
plt.plot(N_values, x_values.imag, label='Im[x(N)]', color='b', linestyle='--')
plt.plot(N_values, y_values.imag, label='Im[y(N)]', color='r', linestyle='--')
plt.plot(N_values, z_values.imag, label='Im[z(N)]', color='g', linestyle='--')
plt.xlabel('N')
plt.ylabel('Imaginary Parts')
plt.title('Imaginary Parts of x(N), y(N), and z(N)')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()

# Print numerical solutions at specific N values
print("Numerical output at specific N values:")
N_specific_values = [N0, Nfinal/2, Nfinal]
for N in N_specific_values:
    x_at_N = np.interp(N, N_values, x_values)
    y_at_N = np.interp(N, N_values, y_values)
    z_at_N = np.interp(N, N_values, z_values)
    print(f"N = {N}")
    print(f"  x(N) = {x_at_N}")
    print(f"  y(N) = {y_at_N}")
    print(f"  z(N) = {z_at_N}")
    print()

"""
Author: Pierros Ntelis, 14 June 2024
"""

########## Solving and plotting the solution functions of time of the system of differential equations ###########

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the system of differential equations
def system(N, variables):
    x, y, z = variables
    dx_dN = -3 * x + x * (3 + y - 3 * z)
    dy_dN = -4 * y + y * (3 + y - 3 * z)
    dz_dN = z * (4 - x - 4 * z)
    return [dx_dN, dy_dN, dz_dN]

# Initial conditions
N0 = 0  # Starting point of the independent variable N
x0 = 0.01
y0 = 0.98
z0 = 0.01
initial_conditions = [x0, y0, z0]

# Define the range for the independent variable N
N_final = 10  # End point of the independent variable N
N_span = [N0, N_final]  # The range of N to solve over
N_eval = np.linspace(N0, N_final, 400)  # Points at which to store the computed solutions

# Solve the differential equations numerically
solution = solve_ivp(system, N_span, initial_conditions, method='RK45', t_eval=N_eval)

# Extract solutions for plotting
N_values = solution.t
x_values = solution.y[0]
y_values = solution.y[1]
z_values = solution.y[2]

# Plot the solutions
plt.ion()
plt.figure(figsize=(12, 6))

plt.plot(N_values, x_values, label='x(N)', color='b')
plt.plot(N_values, y_values, label='y(N)', color='r')
plt.plot(N_values, z_values, label='z(N)', color='g')

plt.xlabel('N')
plt.ylabel('Functions')
plt.title('Numerical Solutions of the Differential Equations')
plt.legend()
plt.grid()

plt.show()

plt.savefig('./savefigs/LCDM_2D_set_of_sim_diff_eqns_num_solutions_correctly_1_tmin_'+str(N0)+'_tmax'+str(N_final)+'_zero_is_'+str(x0)+'.pdf')

# Print numerical solutions at specific N values
print("Numerical output at specific N values:")
N_specific_values = [N0, N_final/2, N_final]
for N in N_specific_values:
    x_at_N = np.interp(N, N_values, x_values)
    y_at_N = np.interp(N, N_values, y_values)
    z_at_N = np.interp(N, N_values, z_values)
    print(f"N = {N:.2f}")
    print(f"  x(N) = {x_at_N:.4f}")
    print(f"  y(N) = {y_at_N:.4f}")
    print(f"  z(N) = {z_at_N:.4f}")
    print()

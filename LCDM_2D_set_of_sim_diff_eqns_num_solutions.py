"""
Author: Pierros Ntelis, 14 June 2024
"""

########## Solving and plotting the solution functions of time of the system of differential equations ###########

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Initialization of variables
t_min, t_max = -13.0, 1.0
N_step_size = 10000
var_zero_is = 1e-2  # 1e-3

# Define the system of differential equations
def system(t, Z):
    x, y = Z
    dxdt = x * (-3. + 3. * x + 4. * y)
    dydt = y * (-4. + 3. * x + 4. * y)
    return [dxdt, dydt]

# Define the initial conditions and time span
initial_conditions = [[var_zero_is, 1.0 - var_zero_is]]
t_span = (t_min, t_max)
t = np.linspace(t_min, t_max, N_step_size)

# Solve the system for each initial condition
plt.ion()
fig, ax1 = plt.subplots(figsize=(14, 7))

initial = initial_conditions[0]
sol = solve_ivp(system, t_span, initial, t_eval=t, dense_output=True)
ax1.plot(sol.t, sol.y[0], '-', label='$x(N)=\\Omega_{\\rm m}(N, z)$')
ax1.plot(sol.t, sol.y[1], '-', label='$y(N)=\\Omega_{\\rm r}(N, z)$')
ax1.plot(sol.t, 1. - sol.y[0] - sol.y[1], '-', label='$z(N)=\\Omega_{\\rm \\Lambda}(N, z)$')
ax1.plot(sol.t, - 1. + sol.y[0] + 4./3.* sol.y[1], '-', label='$w_{eff}(N)$')
ax1.set_title(f'Initial Condition: {initial}')
ax1.set_xlabel('Lapse function $N = \log[a(t)]$')
ax1.set_ylabel('Solution')
ax1.set_ylim(-1.5, 1.5)
ax1.legend()
ax1.grid()

# Add a secondary x-axis for the redshift z
def N_to_redshift(N):
    return (np.exp(-N) - 1)

def redshift_to_N(z):
    return -np.log(z + 1)

# Create a secondary x-axis
secax = ax1.secondary_xaxis('top', functions=(N_to_redshift, redshift_to_N))
secax.set_xlabel('Redshift $z$')

# Set the ticks on the secondary x-axis to match the primary axis at 3 significant figures
z_ticks = [0, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000]
secax.set_xticks(z_ticks)
redshift_from_N_extra_ticks = np.round(N_to_redshift(np.array([-10, -8, -6, -4, -2])),0)
extra_ticks = list(redshift_from_N_extra_ticks)
secax.set_xticks( sorted(list(secax.get_xticks()) + extra_ticks)[::-1] )

# Create labels in scientific notation
#labels = ['$0$', '$1 \\times 10^0$', '$1 \\times 10^1$', '$1 \\times 10^2$', '$1 \\times 10^3$', '$1 \\times 10^4$', '$1 \\times 10^5$', '$1 \\times 10^6$', '$1 \\times 10^7$', '$1 \\times 10^8$']
#labels = ['$8 \\times 10^5$', '$7 \\times 10^5$', '$6 \\times 10^5$', '$5 \\times 10^5$', '$4 \\times 10^5$', '$3 \\times 10^5$', '$2 \\times 10^5$', '$1 \\times 10^5$', '$0$']
secax.set_xticklabels(sorted(list(secax.get_xticks()) + extra_ticks)[::-1], rotation=45, ha='left', rotation_mode='default')


plt.tight_layout()
plt.show()
plt.savefig('./savefigs/LCDM_2D_set_of_sim_diff_eqns_num_solutions_tmin_'+str(t_min)+'_tmax'+str(t_max)+'_zero_is_'+str(var_zero_is)+'.pdf')

########### solution as a numeric output ########### 
results = {}
initial = initial_conditions[0]
sol = solve_ivp(system, t_span, initial, t_eval=t)
results[tuple(initial)] = sol

# Print the results
for initial, sol in results.items():
    print(f"Initial Condition {initial}:")
    print("t: ", sol.t[:10], "...")  # Print first 10 time points for brevity
    print("x(t): ", sol.y[0][:10], "...")
    print("y(t): ", sol.y[1][:10], "...")
    print("z(t): ", 1 - sol.y[0][:10] - sol.y[1][:10], "...")
    print("\n")



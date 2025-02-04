"""
Author: Pierros Ntelis, 9 October 2024
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

_lambda_ = 1  # example value for lambda#0 #10

# Define the system of differential equations
def system(t, variables):
    m, x, y = variables
    f_mxy = (4 - m + 2 * x**2 - 4 * y**2)
    dm_dt = -3 * m + m * f_mxy
    dx_dt = -3 * x + (np.sqrt(3)/np.sqrt(2)) * _lambda_  * y**2 + (1/2) * x * f_mxy
    dy_dt = - _lambda_ * (np.sqrt(3)/np.sqrt(2)) * y * x + (1/2) * y * f_mxy
    return [dm_dt, dx_dt, dy_dt]

# Set the initial conditions and time span
var_zero_is = 2.5e-9 #2.5e-9
initial_conditions = [0.009, 0.0, var_zero_is]  # initial values for m, x, y # 0.009

N_min = -12   #-13
N_max = 1     #2
N_span = (N_min, N_max)  # lapse time function range
N_eval = np.linspace(*N_span, 10000)  # lapse time function values to evaluate

# Solve the system of equations
solution = solve_ivp(system, N_span, initial_conditions, t_eval=N_eval, method='RK45', dense_output=True)

# Extract the results
N = solution.t
m, x, y = solution.y

r = 1-m-x**2-y**2

# Plot the results
plt.ion()
plt.figure(1, figsize=(14, 7))
plt.clf()
plt.plot(N, m, label='$\Omega_m(N)$', color='blue')
plt.plot(N, x, label='$x(N)$', color='red')
plt.plot(N, y, label='$y(N)$', color='green')
plt.plot(N, r, label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
plt.title('Numerical solution for the system of differential equations \n Energy density ratios, $\Omega_Z[N(t)]$ vs lapse function, $N(t) = \ln[a(t)]$, $\lambda=$%1.1e'%(_lambda_))
plt.xlabel('lapse function, $N(t) = \ln[a(t)]$')
plt.ylabel('$\Omega_Z(N)$')
plt.ylim([-1.1,1.1])
plt.legend()
plt.grid(True)
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mxyr_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')

# Plot the results with the effective equation of state

weff = 1./3.*r + x**2. - y**2.

plt.ion()
plt.figure(2, figsize=(14, 7))
plt.clf()
plt.plot(N, m, label='$\Omega_m(N)$', color='blue')
plt.plot(N, x, label='$x(N)$', color='red')
plt.plot(N, y, label='$y(N)$', color='green')
plt.plot(N, r, label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
plt.plot(N, weff, '--', label='$w_{\\rm eff}(N)=(\\frac{1}{3}\Omega_r+x^2-y^2)(N)$', color='purple')
plt.title('Numerical solution for the system of differential equations \n Energy density ratios, $\Omega_Z[N(t)]$ vs lapse function, $N(t) = \ln[a(t)]$, $\lambda=$%1.1e'%(_lambda_))
plt.xlabel('lapse function, $N(t) = \ln[a(t)]$')
plt.ylabel('$\Omega_Z(N)$')
plt.ylim([-1.1,1.1])
plt.legend()
plt.grid(True)
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mxyr_weos_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')


# Plot the results with the effective equation of state

weff = 1./3.*r + x**2. - y**2.

plt.ion()
plt.figure(3, figsize=(14, 7))
plt.clf()
plt.plot(N, m, label='$\Omega_m(N)$', color='blue')
plt.plot(N, x**2 + y**2, label='$\Omega_{\phi}(N)$', color='green')
plt.plot(N, r, label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
plt.plot(N, weff, '--', label='$w_{\\rm eff}(N)=(\\frac{1}{3}\Omega_r+x^2-y^2)(N)$', color='purple')
plt.title('Numerical solution for the system of differential equations \n Energy density ratios, $\Omega_Z[N(t)]$ vs lapse function, $N(t) = \ln[a(t)]$, $\lambda=$%1.1e'%(_lambda_))
plt.xlabel('lapse function, $N(t) = \ln[a(t)]$')
plt.ylabel('$\Omega_Z(N)$')
plt.ylim([-1.1,1.1])
plt.legend()
plt.grid(True)
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mphir_weos_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')


# Add second upper axis the redshift for mxyr
plt.ion()
fig, ax1 = plt.subplots(figsize=(14, 7))

ax1.plot(N, m, label='$\Omega_m(N)$', color='blue')
ax1.plot(N, x, label='$x(N)$', color='red')
ax1.plot(N, y, label='$y(N)$', color='green')
ax1.plot(N, r, label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
ax1.plot(N, weff, '--', label='$w_{\\rm eff}(N)=(\\frac{1}{3}\Omega_r+x^2-y^2)(N)$', color='purple')
ax1.set_title(f'Initial Condition: {initial_conditions}, $\lambda=$%1.1e'%(_lambda_))
ax1.set_xlabel('Lapse function $N = \ln[a(t)]$')
ax1.set_ylabel('Solution, $\Omega_Z(N)$')
plt.ylim([-1.1,1.1])
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
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mxyr_weos_with_redshift_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')


# Add second upper axis the redshift for mxyr
plt.ion()
fig, ax1 = plt.subplots(figsize=(14, 7))

ax1.plot(N, m, label='$\Omega_m(N)$', color='blue')
ax1.plot(N, x**2 + y**2, label='$\Omega_{\phi}(N)$', color='green')
ax1.plot(N, r, label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
ax1.plot(N, weff, '--', label='$w_{\\rm eff}(N)=(\\frac{1}{3}\Omega_r+x^2-y^2)(N)$', color='purple')
ax1.set_title(f'Initial Condition: {initial_conditions}, $\lambda=$%1.1e'%(_lambda_))
ax1.set_xlabel('Lapse function $N = \ln[a(t)]$')
ax1.set_ylabel('Solution, $\Omega_Z(N)$')
ax1.set_ylim([-1.1,1.1])
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
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mphir_weos_with_redshift_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')


print('# It finds the intersection points, and plots a new plot with points printed')

# Function to find intersection points
def find_intersections(x, y1, y2):
    # Find where the sign of the difference changes
    diff = y1 - y2
    sign_change_indices = np.where(np.diff(np.sign(diff)))[0]

    intersections = []
    
    # Linearly interpolate to find more accurate intersection points
    for i in sign_change_indices:
        x1, x2 = x[i], x[i + 1]
        y1_diff, y2_diff = diff[i], diff[i + 1]
        
        # Linear interpolation to find intersection point
        intersection_x = x1 - y1_diff * (x2 - x1) / (y2_diff - y1_diff)
        intersections.append(intersection_x)
    
    return intersections

phi = x**2 + y**2

print('# Find intersections between m(N) and r(N)')
intersection_N_m_r = find_intersections(N, m, r)
intersection_m_r   = np.interp(intersection_N_m_r, N, r)

print('# Find intersections between m(N) and phi(N)')
intersection_N_m_phi = find_intersections(N, m, phi)
intersection_m_phi   = np.interp(intersection_N_m_phi, N, phi)

print('# Find intersections between r(N) and phi(N)')
intersection_N_r_phi = find_intersections(N, r, phi)
intersection_r_phi   = np.interp(intersection_N_r_phi, N, phi)


# Add second upper axis the redshift for mxyr
plt.ion()
fig, ax1 = plt.subplots(figsize=(14, 7))

ax1.plot(N, m, label='$\Omega_m(N)$', color='blue')
ax1.plot(N, x**2 + y**2, label='$\Omega_{\phi}(N)$', color='green')
ax1.plot(N, r, label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
ax1.plot(N, weff, '--', label='$w_{\\rm eff}(N)=(\\frac{1}{3}\Omega_r+x^2-y^2)(N)$', color='purple')
ax1.set_title(f'Initial Condition: {initial_conditions}, $\lambda=$%1.1f'%(_lambda_))
ax1.set_xlabel('Lapse function $N = \ln[a(t)]$')
ax1.set_ylabel('Solution, $\Omega_Z(N)$')
ax1.set_ylim([-1.1,1.1])
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



# Plotting the curves

# Plot intersection points for m(N) and r(N)
ax1.scatter(intersection_N_m_r, intersection_m_r, color='blue', label='$m(N) ∩ r(N)$=(%1.1f,%1.1f), $z_{mr}$=%1.1f'%(intersection_N_m_r[0], intersection_m_r[0],N_to_redshift(intersection_N_m_r[0])))

ax1.scatter(intersection_N_r_phi, intersection_r_phi, color='magenta', label='$r(N) ∩ \\tilde{\phi}(N)$=(%1.1f,%1.1f), $z_{r\phi}$=%1.1f'%(intersection_N_r_phi[0], intersection_r_phi[0],N_to_redshift(intersection_N_r_phi[0])))

ax1.scatter(intersection_N_m_phi, intersection_m_phi, color='red', label='$m(N) ∩ \\tilde{\phi}(N)$=(%1.1f,%1.1f), $z_{m\phi}$=%1.1f'%(intersection_N_m_phi[0], intersection_m_phi[0],N_to_redshift(intersection_N_m_phi[0])))


#plt.xlabel('$N = \ln[a(t)]$')
#plt.ylabel('Z(N)')
ax1.legend()
ax1.grid(True)

plt.tight_layout()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mphir_weos_with_redshift_intersections_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')

print("Intersection points (N, m(N) ∩ r(N)): ", list(zip(intersection_N_m_r, intersection_m_r)))


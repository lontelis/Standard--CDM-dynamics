"""
Author: Pierros Ntelis, 22 October 2024
"""

print("""
########## ########## ########## ########## ##########  
########## Solving and ploting 
########## the solutions of 
########## 3D system of simultaneous differential equations for LCDM
########## we add also the constratin 1 = x + y + z, on all 3 equations
########## The solution x(N), y(N), z(N), 
########## functions of lapse time $N$
########## one initial conditions is chosen
########## with intersection points, to describe each equality epoch.
########## with the effective equation of state weff = 1/3\Omega_r - \Omega_Lambda
########## adding the intersection points.
########## ########## ########## ########## ########## 
""")


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

N_min,N_max = -12.0,1.0
N_step_size = 10000
var_zero_is = 0 #1e-18#3e-18 # 1e-3

# Define the system of differential equations
def system(N, Z):
    x, y, z = Z
    dx_dN = x*( -3 + 3*x + 4*y)
    dy_dN = y*( -4 + 3*x + 4*y) 
    dz_dN = (1-x-y)*(3*x + 4*y)
    return [dx_dN, dy_dN, dz_dN]

# Define the initial conditions and eta span
#initial_conditions = [0.01, 1.0-0.01, 0.0]
initial_conditions = [0.009, 1-0.009-var_zero_is, var_zero_is]

N_span = (N_min,N_max)
N = np.linspace(N_min,N_max, N_step_size)

# Solve the system for each initial condition
plt.ion()
fig, axes = plt.subplots(1, 1, figsize=(14, 7))
#axes = axes.flatten()


sol = solve_ivp(system, N_span, initial_conditions, t_eval=N, dense_output=True)
ax1 = axes

N_values     = sol.t
Omega_m      = sol.y[0]
Omega_r      = sol.y[1]
Omega_Lambda = 1-sol.y[0]-sol.y[1]
weff_LCDM = 1/3*Omega_r-Omega_Lambda

ax1.plot(N_values, Omega_m, 'b--', label='$x(N)=\Omega_m(N)$')
ax1.plot(N_values, Omega_r, 'y--',label='$y(N)=\Omega_r(N)$')
ax1.plot(N_values, Omega_Lambda, 'g--',label='$z(N)=\Omega_\Lambda(N)$')

#ax1.plot(N_values, 1-sol.y[0]-sol.y[1], '--', label='$(1-x-y)(N)=\Omega_\Lambda(N)$')  
ax1.plot(N_values, weff_LCDM, 'm--', label='$w_{\\rm{eff}}(N)=\\frac{1}{3}\Omega_r(N)-\Omega_\Lambda(N)$')    

ax1.set_title(f'3D Solutions and Initial Condition: {initial_conditions}',fontsize=(15))
ax1.set_xlabel('$N = \\ln [a(t)]$',size=20)
ax1.set_ylabel('Solution',size=20)
ax1.set_yscale('linear')
ax1.legend(fontsize=(15))
ax1.grid()


# Add a secondary x-axis for the redshift z
def N_to_redshift(N):
    return (np.exp(-N) - 1)

def redshift_to_N(z):
    return -np.log(z + 1)

# Create a secondary x-axis
secax = ax1.secondary_xaxis('top', functions=(N_to_redshift, redshift_to_N))
secax.set_xlabel('Redshift $z$',size=15)

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


print('# Find intersections between m(N) and r(N)')
intersection_N_m_r = find_intersections(N_values, Omega_m, Omega_r)
intersection_m_r   = np.interp(intersection_N_m_r, N_values, Omega_r)

print('# Find intersections between m(N) and Lambda(N)')
intersection_N_m_L = find_intersections(N_values, Omega_m, Omega_Lambda)
intersection_m_L   = np.interp(intersection_N_m_L, N_values, Omega_Lambda)

print('# Find intersections between r(N) and L(N)')
intersection_N_r_L = find_intersections(N_values, Omega_r, Omega_Lambda)
intersection_r_L   = np.interp(intersection_N_r_L, N_values, Omega_Lambda)

print('# Plot intersection points')
plt.scatter(intersection_N_m_r, intersection_m_r, color='blue', \
    label='$\Omega_m(N) ∩ \Omega_r(N)$=(%1.1f,%1.1f), z=%1.1f'%(intersection_N_m_r[0], intersection_m_r[0],N_to_redshift(intersection_N_m_r[0])))
plt.scatter(intersection_N_r_L, intersection_r_L, color='green', \
    label='$\Omega_r(N) ∩ \Omega_{\Lambda}(N)$=(%1.1f,%1.1f), z=%1.1f'%(intersection_N_r_L[0], intersection_r_L[0],N_to_redshift(intersection_N_r_L[0])))
plt.scatter(intersection_N_m_L, intersection_m_L, color='red', \
    label='$\Omega_m(N) ∩ \Omega_{\Lambda}(N)$=(%1.1f,%1.1f), z=%1.1f'%(intersection_N_m_L[0], intersection_m_L[0],N_to_redshift(intersection_N_m_L[0])))

plt.legend(fontsize=(12))

plt.tight_layout()
plt.show()
plt.savefig('./savefigs/phase_portraits_for_cosmology_from_ChatGPT_solutions_3D_4eqns_1initcond_with_zaxis_weff_intersections_delta_eta_tmin'+str(N_min)+'_tmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'.pdf')

print("""
########## ########## ########## ########## ##########  
########## Solving and ploting 
########## the solutions of 
########## 3D system of simultaneous differential equations for phiLCDM
########## we add also the constratin 1 = x + y + z, on all 3 equations
########## The solution m(N), x(N), y(N)
########## functions of lapse time function $N$
########## ########## ########## ########## ########## 
""")

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


# Plot the results with the effective equation of state

weff = 1./3.*r + x**2. - y**2.



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

ax1.plot(N, m,   '--',label='$\Omega_m(N)$', color='blue')
ax1.plot(N, x**2 + y**2, '--',label='$\Omega_{\phi}(N)$', color='green')
ax1.plot(N, r, '--',label='$\Omega_r(N)=(1-\Omega_m-x^2-y^2)(N)$', color='orange')
ax1.plot(N, weff, '--', label='$w_{\\rm eff}(N)=(\\frac{1}{3}\Omega_r+x^2-y^2)(N)$', color='purple')

#ax1.set_title(f'Initial Condition: {initial_conditions}, $\lambda=$%1.1f'%(_lambda_),size=15)
ax1.set_title('Initial Condition: [$0.009, 0, 2.5 \\times 10^{-9}$], $\lambda=$%1.0f'%(_lambda_), size=15)
ax1.set_xlabel('Lapse function $N = \ln[a(t)]$',size=20)
ax1.set_ylabel('Solution, $\Omega_Z(N)$',size=20)
ax1.set_ylim([-1.1,1.1])
ax1.legend(fontsize=(10))
ax1.grid()

# Add a secondary x-axis for the redshift z
def N_to_redshift(N):
    return (np.exp(-N) - 1)

def redshift_to_N(z):
    return -np.log(z + 1)

# Create a secondary x-axis
secax = ax1.secondary_xaxis('top', functions=(N_to_redshift, redshift_to_N))
secax.set_xlabel('Redshift $z$',size=15)

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
ax1.scatter(intersection_N_m_r, intersection_m_r, color='blue',  \
    label='$m(N) ∩ r(N)$=(%1.1f,%1.1f), $z_{mr}$=%1.1f'%(intersection_N_m_r[0], intersection_m_r[0],N_to_redshift(intersection_N_m_r[0])))

ax1.scatter(intersection_N_r_phi, intersection_r_phi, color='magenta', \
    label='$r(N) ∩ \\tilde{\phi}(N)$=(%1.1f,%1.1f), $z_{r\phi}$=%1.1f'%(intersection_N_r_phi[0], intersection_r_phi[0],N_to_redshift(intersection_N_r_phi[0])))

ax1.scatter(intersection_N_m_phi, intersection_m_phi, color='red', \
    label='$m(N) ∩ \\tilde{\phi}(N)$=(%1.1f,%1.1f), $z_{m\phi}$=%1.1f'%(intersection_N_m_phi[0], intersection_m_phi[0],N_to_redshift(intersection_N_m_phi[0])))


#plt.xlabel('$N = \ln[a(t)]$')
#plt.ylabel('Z(N)')
ax1.legend(fontsize=(12))
ax1.grid(True)

plt.tight_layout()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_num_solutions_mphir_weos_with_redshift_intersections_delta_N_Nmin'+str(N_min)+'_Nmax'+str(N_max)+'_zero_is_'+str(var_zero_is)+'_lambda_'+str(_lambda_)+'.pdf')

print("Intersection points (N, m(N) ∩ r(N)): ", list(zip(intersection_N_m_r, intersection_m_r)))


print('# Compare the outputs of the two models')

plt.ion()
fig, ax3 = plt.subplots(1, 1, figsize=(14, 7))

ax3.plot(N, m/Omega_m, label='$m$', color='blue')
ax3.plot(N, (x**2 + y**2)/Omega_Lambda, label='$\\tilde{\phi}/\Lambda$', color='green')
#ax3.plot(N, (y**2)/Omega_Lambda, '--', label='$y^2/\Lambda$', color='red')
ax3.plot(N, r/Omega_r, label='$r$', color='orange')
ax3.plot(N, weff/weff_LCDM, '--', label='$w_{\\rm eff}$', color='purple')
#ax3.set_title(f'Initial Condition: {initial_conditions}, $\lambda=$%1.1f'%(_lambda_), size=15)
ax3.set_title('Initial Condition: [$0.009, 0, 2.5 \\times 10^{-9}$], $\lambda=$%1.0f'%(_lambda_), size=15)
ax3.set_xlabel('Lapse function $N = \ln[a(t)]$',size=20)
ax3.set_ylabel('Solution, $\Omega_s^{\phi{\\rm CDM}}(N)/\Omega_s^{\Lambda{\\rm CDM}}(N)$',size=20)
ax3.set_ylim([0.0,1.8])
ax3.legend(fontsize=(15))
ax3.grid()

plt.savefig('./savefigs/ChatGPT_6D_to_3D_phiLCDM_comparison_with_3D_4eqns_LCDM_set_of_sim_diff_eqns_num_solutions.pdf')
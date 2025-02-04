"""
Author: Pierros Ntelis, 11 August 2024
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the system of differential equations
def system(N, y):
    m, r, Lambda, x, y_, lambda_, Gamma = y
    
    # Differential equations
    dm_dN = -3*m + m * (3 + r - 3*Lambda + 3*x**2 - 3*y_**2)
    dr_dN = -4*r + r * (3 + r - 3*Lambda + 3*x**2 - 3*y_**2)
    dLambda_dN = Lambda * (3 + r - 3*Lambda + 3*x**2 - 3*y_**2)
    dx_dN = (-3/2)*x + (np.sqrt(3)/np.sqrt(2))*lambda_*y_**2 + (1/2)*x*(r - 3*Lambda + 3*x**2 - 3*y_**2)
    dy_dN = -(lambda_*np.sqrt(3)/np.sqrt(2))*y_*x + (1/2)*y_*(3 + r - 3*Lambda + 3*x**2 - 3*y_**2)
    dlambda_dN = -np.sqrt(6)*x*lambda_**2*(Gamma - 1)
    dGamma_dN = 0  # Since Gamma' = 0

    return [dm_dN, dr_dN, dLambda_dN, dx_dN, dy_dN, dlambda_dN, dGamma_dN]

# Initial conditions m, r, Λ, x, y, λ, Γ
N0 = -13
#y0 = [0.1, 1-(1)*0.01-(1+1+1)*0.1, 0.0, 0.1, 0.01, 0.1, 1.0]
#y0 = [1-0.11, 0.0, 0.0, 0.01, 0.01, 0.1, 1.0] # m,x,y model
#y0 = [0.01, 1-0.02, 0.01, 0.0, 0.0, 0.1, 1.0] # m,r,Λ model
y0 = [0.01, 1-0.01, 0.0 , 0.0, 1e-9, 0.1, 1.0] # m,r model

# Define the range of N over which to solve
N_final = 2  # Choose a suitable final N
N_span  = [N0, N_final]
N_eval = np.linspace(*N_span, 300) 

# Solve the system using solve_ivp
sol = solve_ivp(system, N_span, y0, t_eval=N_eval, method='RK45', dense_output=True)

# Extract the results
N_values = sol.t
m_values, r_values, Lambda_values, x_values, y_values, lambda_values, Gamma_values = sol.y

phi_values = x_values**2. + y_values**2.

# Plot the results
plt.ion()
plt.figure(1,figsize=(14, 10))
plt.clf()

plt.subplot(3, 3, 1)
plt.plot(N_values, m_values, label='m(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('m(N)')
#plt.title('m(N) vs N')
plt.title(' \n ')
plt.grid(True)

plt.subplot(3, 3, 2)
plt.plot(N_values, r_values, label='r(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('r(N)')
#plt.title('r(N) vs N')
plt.title(' \n ')
plt.grid(True)

plt.subplot(3, 3, 3)
plt.plot(N_values, Lambda_values, label='Λ(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('Λ(N)')
#plt.title('Λ(N) vs N')
plt.title(' \n ')
plt.grid(True)

plt.subplot(3, 3, 4)
plt.plot(N_values, x_values, label='x(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('x(N)')
plt.title('x(N) vs N')
plt.grid(True)

plt.subplot(3, 3, 5)
plt.plot(N_values, y_values, label='y(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('y(N)')
plt.title('y(N) vs N')
plt.grid(True)

plt.subplot(3, 3, 6)
plt.plot(N_values, lambda_values, label='λ(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('λ(N)')
plt.title('λ(N) vs N')
plt.grid(True)


plt.subplot(3, 3, 7)
plt.plot(N_values, Gamma_values, label='Γ(N)')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('Γ(N)')
plt.title('Γ(N) vs N')
plt.grid(True)

plt.subplot(3, 3, 8)
plt.plot(N_values, x_values**2. + y_values**2., label='$\\tilde{\phi}(N)$')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('$\\tilde{\phi}(N)$')
plt.title('$\\tilde{\phi}(N)$ vs N')
plt.grid(True)

plt.subplot(3, 3, 9)
plt.plot(N_values, m_values, label='m(N)')
plt.plot(N_values, r_values, label='r(N)')
plt.plot(N_values, Lambda_values, label='Λ(N)')
plt.plot(N_values, x_values**2. + y_values**2., label='$\\tilde{\phi}(N)$')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('Z(N)')
plt.legend()
#plt.title('Z(N) vs N')
plt.title(' \n ')
plt.grid(True)

plt.suptitle('Initial Conditions: [m,r,Λ,x,y,λ,Γ]=%s'%(y0))

plt.tight_layout()
plt.show()

plt.savefig('./ChatGPT_scalarLCDM_7D_plots/ChatGPT_scalarLCDM_7D_plots_N_%s_-_%s_initial_conditions_Gamma_%s.pdf'%(N0,N_final,y0))

plt.figure(2, figsize=(14, 7))
plt.clf()
plt.plot(N_values, m_values, label='m(N)')
plt.plot(N_values, r_values, label='r(N)')
plt.plot(N_values, Lambda_values, label='Λ(N)')
plt.plot(N_values, phi_values, label='$\\tilde{\phi}(N)$')
plt.plot(N_values, 1./3.*r_values - Lambda_values + x_values**2. - y_values**2., label='$w_{\\rm eff}(N)$')
plt.xlabel('$N = ln[a(t)]$')
plt.ylabel('Z(N)')
plt.legend()
plt.grid(True)
plt.title('Initial Conditions: [m,r,Λ,x,y,λ,Γ]=%s'%(y0))

plt.tight_layout()
plt.show()

plt.savefig('./ChatGPT_scalarLCDM_7D_plots/ChatGPT_scalarLCDM_7D_plots_N_%s_-_%s_initial_conditions_Gamma_%s_main.pdf'%(N0,N_final,y0))

# Print the numerical output at specific points
print("Numerical output at specific N values:")
N_specific_values = [N0, N_final, 10]
for N in N_specific_values:
    y_at_N = sol.sol(N)
    print(f"N = {N}")
    print(f"  m(N) = {y_at_N[0]}")
    print(f"  r(N) = {y_at_N[1]}")
    print(f"  Λ(N) = {y_at_N[2]}")
    print(f"  x(N) = {y_at_N[3]}")
    print(f"  y(N) = {y_at_N[4]}")
    print(f"  λ(N) = {y_at_N[5]}")
    print(f"  Γ(N) = {y_at_N[6]}")
    print()


print('# It finds the intersection points, and plots a new plot with points printed')

import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
#from scipy.optimize import fsolve

# Your existing data

"""
# Interpolation of the curves to continuous functions
m_interp = interp1d(N_values, m_values, kind='linear', fill_value="extrapolate")
r_interp = interp1d(N_values, r_values, kind='linear', fill_value="extrapolate")
Lambda_interp = interp1d(N_values, Lambda_values, kind='linear', fill_value="extrapolate")
phi_interp = interp1d(N_values, x_values**2 + y_values**2, kind='linear', fill_value="extrapolate")
w_eff_interp = interp1d(N_values, 1./3.*r_values - Lambda_values + x_values**2 - y_values**2, kind='linear', fill_value="extrapolate")
"""

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
intersection_N_m_r = find_intersections(N_values, m_values, r_values)
intersection_m_r   = np.interp(intersection_N_m_r, N_values, r_values)
"""
print('# Find intersections between m(N) and Λ(N')
intersection_N_m_L = find_intersections(N_values, m_values, Lambda_values)
intersection_m_L   = np.interp(intersection_N_m_L, N_values, Lambda_values)

print('# Find intersections between r(N) and Λ(N')
intersection_N_r_L = find_intersections(N_values, r_values, Lambda_values)
intersection_r_L   = np.interp(intersection_N_r_L, N_values, r_values)
"""
print('# Find intersections between m(N) and phi(N)')
intersection_N_m_phi = find_intersections(N_values, m_values, phi_values)
intersection_m_phi   = np.interp(intersection_N_m_phi, N_values, phi_values)

print('# Find intersections between r(N) and phi(N)')
intersection_N_r_phi = find_intersections(N_values, r_values, phi_values)
intersection_r_phi   = np.interp(intersection_N_r_phi, N_values, phi_values)
"""
print('# Find intersections between r(N) and phi(N)')
intersection_N_L_phi = find_intersections(N_values, Lambda_values, phi_values)
intersection_L_phi   = np.interp(intersection_N_L_phi, N_values, phi_values)
"""


# Plotting the curves
plt.figure(3, figsize=(14, 7))
plt.clf()
plt.plot(N_values, m_values, label='$\Omega_m(N)$')
plt.plot(N_values, r_values, label='$\Omega_r(N)$')
"""
plt.plot(N_values, Lambda_values, label='Λ(N)')
"""
plt.plot(N_values, x_values**2. + y_values**2., label='$\\tilde{\phi}(N)$')
plt.plot(N_values, 1./3.*r_values - Lambda_values + x_values**2. - y_values**2., label='$w_{\\rm eff}(N)$')

# Plot intersection points for m(N) and r(N)
plt.scatter(intersection_N_m_r, intersection_m_r, color='blue', label='$m(N) ∩ r(N)$=(%1.1e,%1.1e)'%(intersection_N_m_r[0], intersection_m_r[0]))
"""
plt.scatter(intersection_N_m_L, intersection_m_L, color='green', label='$m(N) ∩ \Lambda(N)$=(%1.1e,%1.1e)'%(intersection_N_m_L[0], intersection_m_L[0]))
plt.scatter(intersection_N_r_L, intersection_r_L, color='orange', label='$r(N) ∩ \Lambda(N)$=(%1.1e,%1.1e)'%(intersection_N_r_L[0], intersection_r_L[0]))
"""
plt.scatter(intersection_N_m_phi, intersection_m_phi, color='red', label='$m(N) ∩ \\tilde{\phi}(N)$=(%1.1e,%1.1e)'%(intersection_N_m_phi[0], intersection_m_phi[0]))

plt.scatter(intersection_N_r_phi, intersection_r_phi, color='magenta', label='$r(N) ∩ \\tilde{\phi}(N)$=(%1.1e,%1.1e)'%(intersection_N_r_phi[0], intersection_r_phi[0]))
"""
plt.scatter(intersection_N_L_phi, intersection_L_phi, color='brown', label='$\Lambda(N) ∩ \\tilde{\phi}(N)$')#=(%1.1e,%1.1e)'%(intersection_N_L_phi[0], intersection_L_phi[0]))
"""

plt.xlabel('$N = \ln[a(t)]$')
plt.ylabel('Z(N)')
plt.legend()
plt.grid(True)
plt.title('Initial Conditions: [m,r,Λ,x,y,λ,Γ]=%s' % (y0))

plt.tight_layout()
plt.show()
plt.savefig('./ChatGPT_scalarLCDM_7D_plots/ChatGPT_scalarLCDM_7D_plots_N_%s_-_%s_initial_conditions_Gamma_%s_main_with_intersections.pdf'%(N0,N_final,y0))

print("Intersection points (N, m(N) ∩ r(N)): ", list(zip(intersection_N_m_r, intersection_m_r)))


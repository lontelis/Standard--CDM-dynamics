"""
Author: Pierros Ntelis, 10 October 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq

# Define the symbolic variables
mu, l_, m = symbols('mu l_ m')

# Coefficients from the cubic equation
B = (3/2) - m - 4*l_**2 + 6*l_**-2
C = -53/2 + (123/4)*l_**-2 - 6*l_**2 - (117/4)*l_**-4 - 4*l_**2*m - (5/4)*m - (7/4)*m**2 - 13*m*l_**-2
D = -4 + 12*l_**2 + (99/4)*l_**-2 - (513/4)*l_**-4 + 117*l_**-6 + (19/2)*m - (5/4)*m**2 - (1/2)*m**3 - 84*m*l_**-2 + 117*m*l_**-4 - m**2*l_**2 - (21/2)*m**2*l_**-2

# Define the cubic equation: A*mu^3 + B*mu^2 + C*mu + D = 0
# cubic_eq = Eq(mu**3 + B*mu**2 + C*mu + D, 0)

# Solve the cubic equation symbolically for specific values of m and l_
def solve_eigenvalues(m_value, l_value):
    B_val = B.subs({m: m_value, l_: l_value})
    C_val = C.subs({m: m_value, l_: l_value})
    D_val = D.subs({m: m_value, l_: l_value})

    # Coefficients in numerical form
    coeffs = [1, B_val, C_val, D_val]

    # Solve the cubic equation
    roots = np.roots(coeffs)
    return roots

# Calculate for specific values of m and lambda
m_val = 0.5
lambda_values = [ -100, -10, -1, -0.1, -0.01, 1e-15, 0.01, 0.1, 1., 10., 100.]
eigenvalues = {l_val: solve_eigenvalues(m_val, l_val) for l_val in lambda_values}

print('# For m=%1.1f'%(m_val))
print('# and for \lambda={'+str(lambda_values)+'}')
print('# the eigenvalues are')
print(eigenvalues)

# Calculate for specific values of m and lambda
m_val = 0.0
lambda_values = [ -100, -10, -1, -0.1, -0.01, 1e-15, 0.01, 0.1, 1., 10., 100.]
eigenvalues = {l_val: solve_eigenvalues(m_val, l_val) for l_val in lambda_values}

print('# For m=%1.1f'%(m_val))
print('# and for \lambda={'+str(lambda_values)+'}')
print('# the eigenvalues are')
print(eigenvalues)

# Calculate for specific values of m and lambda
m_val = 1.0
lambda_values = [ -100, -10, -1, -0.1, -0.01, 1e-15, 0.01, 0.1, 1., 10., 100.]
eigenvalues = {l_val: solve_eigenvalues(m_val, l_val) for l_val in lambda_values}

print('# For m=%1.1f'%(m_val))
print('# and for \lambda={'+str(lambda_values)+'}')
print('# the eigenvalues are')
print(eigenvalues)

# Calculate for specific values of m and lambda
m_val = 0.3
lambda_values = [ -100, -10, -1, -0.1, -0.01, 1e-15, 0.01, 0.1, 1., 10., 100.]
eigenvalues_m_leq_05 = {l_val: solve_eigenvalues(m_val, l_val) for l_val in lambda_values}

print('# For m=%1.1f'%(m_val))
print('# and for \lambda={'+str(lambda_values)+'}')
print('# the eigenvalues are')
print(eigenvalues_m_leq_05)

# Calculate for specific values of m and lambda
m_val = 0.7
lambda_values = [ -100, -10, -1, -0.1, -0.01, 1e-15, 0.01, 0.1, 1., 10., 100.]
eigenvalues_m_geq_05 = {l_val: solve_eigenvalues(m_val, l_val) for l_val in lambda_values}

print('# For m=%1.1f'%(m_val))
print('# and for \lambda={'+str(lambda_values)+'}')
print('# the eigenvalues are')
print(eigenvalues_m_geq_05)


abs_lambda_vals = np.array([1e-15, 0.01, 0.1, 1., 10., 100.])
eigenvalues_m_leq_05_real_mxy = np.zeros((len(abs_lambda_vals),3))
eigenvalues_m_geq_05_real_mxy = np.zeros((len(abs_lambda_vals),3))
eigenvalues_m_leq_05_imag_mxy = np.zeros((len(abs_lambda_vals),3))
eigenvalues_m_geq_05_imag_mxy = np.zeros((len(abs_lambda_vals),3))

for i in range(len(abs_lambda_vals)): 
	print(i)
	print(abs_lambda_vals[i])
	eigenvalues_m_leq_05_real_mxy[i] = np.round(np.real(eigenvalues_m_leq_05[abs_lambda_vals[i]]),0)
	eigenvalues_m_geq_05_real_mxy[i] = np.round(np.real(eigenvalues_m_geq_05[abs_lambda_vals[i]]),0)
	eigenvalues_m_leq_05_imag_mxy[i] = np.round(np.imag(eigenvalues_m_leq_05[abs_lambda_vals[i]]),0)
	eigenvalues_m_geq_05_imag_mxy[i] = np.round(np.imag(eigenvalues_m_geq_05[abs_lambda_vals[i]]),0)

plt.ion()
plt.figure(1, figsize=(14,7))
plt.clf()

plt.plot(abs_lambda_vals, eigenvalues_m_leq_05_real_mxy[:,0] , '.-', label='$m\leq 1/2$')
plt.plot(abs_lambda_vals, eigenvalues_m_geq_05_real_mxy[:,0] , '.:', label='$m\geq 1/2$')

plt.xlabel('$\lambda$') 
plt.ylabel('$Re[\\vec{ev}_m](\lambda)$') 
plt.legend()
plt.title('Real part of eigenvalues for $\Omega_m(N)$')
plt.grid()
plt.xscale('log')
plt.yscale('symlog')
plt.draw()
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_Jacobian_eigenvalues_ml_real.pdf')

plt.ion()
plt.figure(2, figsize=(14,7))
plt.clf()

plt.plot(abs_lambda_vals, eigenvalues_m_leq_05_real_mxy[:,1] , '.-', label='$m\leq 1/2$')
plt.plot(abs_lambda_vals, eigenvalues_m_geq_05_real_mxy[:,1] , '.:', label='$m\geq 1/2$')

plt.xlabel('$\lambda$') 
plt.ylabel('$Re[\\vec{ev}_x](\lambda)$') 
plt.legend()
plt.title('Real part of eigenvalues for $\Omega_{x^2}(N)$')
plt.grid()
plt.xscale('log')
plt.yscale('symlog')
plt.draw()
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_Jacobian_eigenvalues_xl_real.pdf')

plt.ion()
plt.figure(3, figsize=(14,7))
plt.clf()

plt.plot(abs_lambda_vals, eigenvalues_m_leq_05_real_mxy[:,2] , '.-', label='$m\leq 1/2$')
plt.plot(abs_lambda_vals, eigenvalues_m_geq_05_real_mxy[:,2] , '.:', label='$m\geq 1/2$')

plt.xlabel('$\lambda$') 
plt.ylabel('$Re[\\vec{ev}_y](\lambda)$') 
plt.legend()
plt.title('Real part of eigenvalues for $\Omega_{y^2}(N)$')
plt.grid()
plt.xscale('log')
plt.yscale('symlog')
plt.draw()
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_Jacobian_eigenvalues_yl_real.pdf')


plt.ion()
plt.figure(4, figsize=(14,7))
plt.clf()

plt.plot(abs_lambda_vals, eigenvalues_m_leq_05_imag_mxy[:,0] , '.-', label='$m\leq 1/2$')
plt.plot(abs_lambda_vals, eigenvalues_m_geq_05_imag_mxy[:,0] , '.:', label='$m\geq 1/2$')

plt.xlabel('$\lambda$') 
plt.ylabel('$Im[\\vec{ev}_m](\lambda)$') 
plt.legend()
plt.title('Imaginary part of eigenvalues for $\Omega_m(N)$')
plt.grid()
plt.xscale('log')
plt.yscale('symlog')
plt.draw()
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_Jacobian_eigenvalues_ml_imag.pdf')

plt.ion()
plt.figure(5, figsize=(14,7))
plt.clf()

plt.plot(abs_lambda_vals, eigenvalues_m_leq_05_imag_mxy[:,1] , '.-', label='$m\leq 1/2$')
plt.plot(abs_lambda_vals, eigenvalues_m_geq_05_imag_mxy[:,1] , '.:', label='$m\geq 1/2$')

plt.xlabel('$\lambda$') 
plt.ylabel('$Im[\\vec{ev}_x](\lambda)$') 
plt.legend()
plt.title('Imaginary part of eigenvalues for $\Omega_{x^2}(N)$')
plt.grid()
plt.xscale('log')
plt.yscale('symlog')
plt.draw()
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_Jacobian_eigenvalues_xl_imag.pdf')

plt.ion()
plt.figure(6, figsize=(14,7))
plt.clf()

plt.plot(abs_lambda_vals, eigenvalues_m_leq_05_imag_mxy[:,2] , '.-', label='$m\leq 1/2$')
plt.plot(abs_lambda_vals, eigenvalues_m_geq_05_imag_mxy[:,2] , '.:', label='$m\geq 1/2$')

plt.xlabel('$\lambda$') 
plt.ylabel('$Im[\\vec{ev}_y](\lambda)$') 
plt.legend()
plt.title('Imaginary part of eigenvalues for $\Omega_{y^2}(N)$')
plt.grid()
plt.xscale('log')
plt.yscale('symlog')
plt.draw()
plt.show()

plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_Jacobian_eigenvalues_yl_imag.pdf')
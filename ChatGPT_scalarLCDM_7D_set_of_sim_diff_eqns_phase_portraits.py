"""
Author: Pierros Ntelis, 11 August 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the system of differential equations
def system(t, vars):
    m, r, Lambda, x, y, lambda_, Gamma = vars
    
    dm_dt = -3 * m + m * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dr_dt = -4 * r + r * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dLambda_dt = Lambda * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dx_dt = -3/2 * x + np.sqrt(3/2) * lambda_ * y**2 + 1/2 * x * (r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dy_dt = -np.sqrt(3/2) * lambda_ * y * x + 1/2 * y * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dlambda_dt = -np.sqrt(6) * x * lambda_**2 * (Gamma - 1)
    dGamma_dt = 0  # Gamma' = 0
    
    return [dm_dt, dr_dt, dLambda_dt, dx_dt, dy_dt, dlambda_dt, dGamma_dt]


# Parameter values
m       = 0.01  # Example parameter
r       = 1.01  # Example parameter
Lambda  = 0.01  # Example parameter
x       = 0.01  # Example parameter
y       = 0.01  # Example parameter
lambda_ = 0.01  # Example parameter
Gamma   = 1.00  # Example parameter

# Define the vector field for the phase portrait of x and y
def vector_fields(m, r, Lambda, x, y, lambda_, Gamma):
    # These are the relevant derivatives from the system for the (x, y) plane
    dm_dt = -3 * m + m * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dr_dt = -4 * r + r * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dLambda_dt = Lambda * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dx_dt = -3/2 * x + np.sqrt(3/2) * lambda_ * y**2 + 1/2 * x * (r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dy_dt = -np.sqrt(3/2) * lambda_ * y * x + 1/2 * y * (3 + r - 3 * Lambda + 3 * x**2 - 3 * y**2)
    dlambda_dt = -np.sqrt(6) * x * lambda_**2 * (Gamma - 1)
    dGamma_dt = 0  # Gamma' = 0
    return dm_dt, dr_dt, dLambda_dt, dx_dt, dy_dt, dlambda_dt, dGamma_dt

print('# Create a grid of points in the (m, r, Lambda, x, y, lambda_, Gamma) planes')

NNN=10
m_vals      = np.linspace(-2, 2, NNN)
r_vals      = np.linspace(-2, 2, NNN)
Lambda_vals = np.linspace(-2, 2, NNN)
x_vals      = np.linspace(-2, 2, NNN)
y_vals      = np.linspace(-2, 2, NNN)
lambda_vals_= np.linspace(-2, 2, NNN)
Gamma_vals  = np.linspace(-2, 2, NNN)

print('# grid')

"""
m, r = np.meshgrid(m_vals, r_vals)

x, y = np.meshgrid(x_vals, y_vals)

Lambda,  lambda_ = np.meshgrid(Lambda_vals, lambda_vals_)

Gamma,  lambda_  = np.meshgrid(Gamma_vals, lambda_vals_)
"""

print('# Compute the vector field at each point')
dm_dt, dr_dt, dLambda_dt, dx_dt, dy_dt, dlambda_dt_, dGamma_dt = vector_fields(m, r, Lambda, x, y, lambda_, Gamma)


# Setting up plotting
plt.ion()

# Plot the vector field
plt.figure(1,figsize=(8, 6))
plt.clf()
plt.quiver(x, y, dx_dt, dy_dt, color='blue', pivot='mid')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('Phase Portrait in ($x, y$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_xy.pdf')


# Plot the vector field
plt.figure(2,figsize=(8, 6))
plt.clf()
plt.quiver(m, r, dm_dt, dr_dt, color='blue', pivot='mid')
plt.xlabel('$m$')
plt.ylabel('$r$')
plt.title('Phase Portrait in ($m, r$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_mr.pdf')

# Plot the vector field
plt.figure(3,figsize=(8, 6))
plt.clf()
plt.quiver(m, Lambda, dm_dt, dLambda_dt, color='blue', pivot='mid')
plt.xlabel('$m$')
plt.ylabel('$\\Lambda$')
plt.title('Phase Portrait in ($m, \\Lambda$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_mL.pdf')

# Plot the vector field
plt.figure(4,figsize=(8, 6))
plt.clf()
plt.quiver(r, Lambda, dr_dt, dLambda_dt, color='blue', pivot='mid')
plt.xlabel('$r$')
plt.ylabel('$\\Lambda$')
plt.title('Phase Portrait in ($r, \\Lambda$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_rL.pdf')


# Plot the vector field
plt.figure(5,figsize=(8, 6))
plt.clf()
plt.quiver(Gamma, lambda_, dGamma_dt, dlambda_dt_, color='blue', pivot='mid')
plt.xlabel('$\\lambda$')
plt.ylabel('$\\Gamma$')
plt.title('Phase Portrait in ($\\Gamma, \\lambda$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_lG.pdf')


# Plot the vector field
plt.figure(6,figsize=(8, 6))
plt.clf()
plt.quiver(lambda_, Lambda, dlambda_dt_, dLambda_dt, color='blue', pivot='mid')
plt.xlabel('$\\lambda$')
plt.ylabel('$\\Lambda$')
plt.title('Phase Portrait in ($\\lambda, \\Lambda$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_lL.pdf')

# Plot the vector field
plt.figure(7,figsize=(8, 6))
plt.clf()
plt.quiver(x, Lambda, dx_dt, dLambda_dt, color='blue', pivot='mid')
plt.xlabel('$x$')
plt.ylabel('$\\Lambda$')
plt.title('Phase Portrait in ($x, \\Lambda$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_xL.pdf')

# Plot the vector field
plt.figure(8,figsize=(8, 6))
plt.clf()
plt.quiver(y, Lambda, dy_dt, dLambda_dt, color='blue', pivot='mid')
plt.xlabel('$y$')
plt.ylabel('$\\Lambda$')
plt.title('Phase Portrait in ($y, \\Lambda$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_yL.pdf')

# Plot the vector field
plt.figure(9,figsize=(8, 6))
plt.clf()
plt.quiver(m, x, dm_dt, dx_dt, color='blue', pivot='mid')
plt.xlabel('$m$')
plt.ylabel('$x$')
plt.title('Phase Portrait in ($m, x$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_mx.pdf')

# Plot the vector field
plt.figure(10,figsize=(8, 6))
plt.clf()
plt.quiver(m, y, dm_dt, dy_dt, color='blue', pivot='mid')
plt.xlabel('$m$')
plt.ylabel('$y$')
plt.title('Phase Portrait in ($m, y$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_my.pdf')

# Plot the vector field
plt.figure(11,figsize=(8, 6))
plt.clf()
plt.quiver(r, x, dr_dt, dx_dt, color='blue', pivot='mid')
plt.xlabel('$r$')
plt.ylabel('$x$')
plt.title('Phase Portrait in ($r, x$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_rx.pdf')


# Plot the vector field
plt.figure(12,figsize=(8, 6))
plt.clf()
plt.quiver(r, y, dr_dt, dy_dt, color='blue', pivot='mid')
plt.xlabel('$r$')
plt.ylabel('$y$')
plt.title('Phase Portrait in ($r, y$) Plane')
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_ry.pdf')


"""


plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_set_of_sim_diff_eqns_phase_portraits_xl.pdf')
"""


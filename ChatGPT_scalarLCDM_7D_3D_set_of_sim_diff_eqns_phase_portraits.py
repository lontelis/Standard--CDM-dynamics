"""
Author: Pierros Ntelis, 11 August 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sys

_lambda_ = int(sys.argv[1])
#_lambda_ = 10.0

# Define the system of differential equations
def system(t, vars):
    m, x, y = vars

    f_mxy = (4 - m + 2 * x**2 - 4 * y**2)
    dm_dt = -3 * m + m * f_mxy
    dx_dt = -3 * x + (np.sqrt(3)/np.sqrt(2)) * _lambda_  * y**2 + (1/2) * x * f_mxy
    dy_dt = - _lambda_ * (np.sqrt(3)/np.sqrt(2)) * y * x + (1/2) * y * f_mxy
    
    return [dm_dt, dx_dt, dy_dt]


# Parameter values
m       = 0.01  # Example parameter
x       = 0.01  # Example parameter
y       = 0.01  # Example parameter

# Define the vector field for the phase portrait of x and y
def vector_fields(m, x, y):
    # These are the relevant derivatives from the system for the (x, y) plane
    f_mxy = (4 - m + 2 * x**2 - 4 * y**2)
    dm_dt = -3 * m + m * f_mxy
    dx_dt = -3 * x + (np.sqrt(3)/np.sqrt(2)) * _lambda_  * y**2 + (1/2) * x * f_mxy
    dy_dt = - _lambda_ * (np.sqrt(3)/np.sqrt(2)) * y * x + (1/2) * y * f_mxy

    return dm_dt, dx_dt, dy_dt

# Define the vector field for the phase portrait of x and y
def vector_fields_xy(x, y):
    # These are the relevant derivatives from the system for the (x, y) plane
    f_mxy = (4 - m + 2 * x**2 - 4 * y**2)
    dx_dt = -3 * x + (np.sqrt(3)/np.sqrt(2)) * _lambda_  * y**2 + (1/2) * x * f_mxy
    dy_dt = - _lambda_ * (np.sqrt(3)/np.sqrt(2)) * y * x + (1/2) * y * f_mxy

    return dx_dt, dy_dt

# Define the vector field for the phase portrait of x and y
def vector_fields_mx(m, x):
    # These are the relevant derivatives from the system for the (x, y) plane
    f_mxy = (4 - m + 2 * x**2 - 4 * y**2)
    dm_dt = -3 * m + m * f_mxy
    dx_dt = -3 * x + (np.sqrt(3)/np.sqrt(2)) * _lambda_  * y**2 + (1/2) * x * f_mxy

    return dm_dt, dx_dt

# Define the vector field for the phase portrait of x and y
def vector_fields_my(m, y):
    # These are the relevant derivatives from the system for the (x, y) plane
    f_mxy = (4 - m + 2 * x**2 - 4 * y**2)
    dm_dt = -3 * m + m * f_mxy
    dy_dt = - _lambda_ * (np.sqrt(3)/np.sqrt(2)) * y * x + (1/2) * y * f_mxy

    return dm_dt, dy_dt

print('# Create a grid of points in the (m, x, y) planes')

NNN=10
extremes = 1
m_vals      = np.linspace(0, extremes, NNN)
x_vals      = np.linspace(-extremes, extremes, NNN)
y_vals      = np.linspace(-extremes, extremes, NNN)


print('# grid')

print('# Compute the vector field at each point')
#dm_dt, dx_dt, dy_dt = vector_fields(m, x, y)

"""
# Setting up plotting
plt.ion()

print('# Compute the vector field at each point, for x, y')
x, y = np.meshgrid(x_vals, y_vals)
dx_dt, dy_dt = vector_fields_xy(x, y)
print('# Plot the vector field')
plt.figure(1,figsize=(8, 6))
plt.clf()
plt.quiver(x, y, dx_dt, dy_dt, color='blue', pivot='mid')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('Phase Portrait in ($x, y$) Plane, $\lambda=$%1.2f'%_lambda_)
plt.grid()
plt.show()
#plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy_lambda_'+str(_lambda_)+'.pdf')

print('# Compute the vector field at each point, for m, x')
m, x = np.meshgrid(m_vals, x_vals)
dm_dt, dx_dt = vector_fields_mx(m, x)
print('# Plot the vector field')
plt.figure(2,figsize=(8, 6))
plt.clf()
plt.quiver(m, x, dm_dt, dx_dt, color='blue', pivot='mid')
plt.xlabel('$m$')
plt.ylabel('$x$')
plt.title('Phase Portrait in ($m, x$) Plane, $\lambda=$%1.2f'%_lambda_)
plt.grid()
plt.show()
#plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx_lambda_'+str(_lambda_)+'.pdf')

print('# Compute the vector field at each point, for m, y')
m, y = np.meshgrid(m_vals, y_vals)
dm_dt, dy_dt = vector_fields_my(m, y)
print('# Plot the vector field')
# Plot the vector field
plt.figure(3,figsize=(8, 6))
plt.clf()
plt.quiver(m, y, dm_dt, dy_dt, color='blue', pivot='mid')
plt.xlabel('$m$')
plt.ylabel('$y$')
plt.title('Phase Portrait in ($m, y$) Plane, $\lambda=$%1.2f'%_lambda_)
plt.grid()
plt.show()
#plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my_lambda_'+str(_lambda_)+'.pdf')
"""

# Setting up plotting
plt.ion()
figsize_var = 7
print('# Compute the vector field at each point, for x, y')
x, y = np.meshgrid(x_vals, y_vals)
#m = m_vals
m = 0
dm_dt, dx_dt, dy_dt = vector_fields(m, x, y)
print('# Plot the vector field')
plt.figure(4,figsize=(figsize_var, figsize_var))
plt.clf()
plt.quiver(x, y, dx_dt, dy_dt, color='blue', pivot='mid')
if _lambda_==0.0 : 
    plt.plot(0,0, 'bo', label='$O \\to S^{xy}_0(0,0)$ saddle' )
    plt.plot(1,0, 'go', label='$O_8 \\to R^{xy}_+(+1,0)$ unstable' )    
    plt.plot(-1,0, 'go', label='$O_8 \\to R^{xy}_-(-1,0)$ unstable' )
    plt.plot(0,+1, 'ro', label='$O_6 \\to A^{xy}_+(0,+1)$ stable' )    
    plt.plot(0,-1, 'ro', label='$O_6 \\to A^{xy}_-(0,-1)$ stable' )    
elif _lambda_==1.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{xy}_0(0,0)$ saddle' )
    plt.plot(1,0, 'go', label='$O_8 \\to R^{xy}_+(+1,0)$ unstable' )    
    plt.plot(-1,0, 'go', label='$O_8 \\to R^{xy}_-(-1,0)$ unstable' )
    plt.plot(+np.sqrt(1/6),+np.sqrt(5/6), 'ro', label='$O_{10} \\to A^{xy}_+(+\\sqrt{\\frac{1}{6}},+\\sqrt{\\frac{5}{6}})$ stable' )    #O_5
    plt.plot(+np.sqrt(1/6),-np.sqrt(5/6), 'ro', label='$O_{10} \\to A^{xy}_-(+\\sqrt{\\frac{1}{6}},-\\sqrt{\\frac{5}{6}})$ stable' )    #O_5
elif _lambda_==-1.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{xy}_0(0,0)$ saddle' )
    plt.plot(1,0, 'go', label='$O_8 \\to R^{xy}_+(+1,0)$ unstable' )    
    plt.plot(-1,0, 'go', label='$O_8 \\to R^{xy}_-(-1,0)$ unstable' )
    plt.plot(-np.sqrt(1/6),+np.sqrt(5/6), 'ro', label='$O_{11} \\to A^{xy}_+(-\\sqrt{\\frac{1}{6}},+\\sqrt{\\frac{5}{6}})$ stable' )   #O_5 
    plt.plot(-np.sqrt(1/6),-np.sqrt(5/6), 'ro', label='$O_{11} \\to A^{xy}_-(-\\sqrt{\\frac{1}{6}},-\\sqrt{\\frac{5}{6}})$ stable' )   #O_5     
elif _lambda_==10.0 :
    plt.plot(0,0, 'ko', label='$O_5 \\to S^{xy}_{p+}(0,0)$ saddle' ) 
    plt.plot(-1,0, 'go', label='$O_8 \\to R^{xy}_-(-1,0)$ unstable' )
    plt.plot(+1,0, 'bo', label='$O_9 \\to S^{xy}_+(+1,0)$ saddle' )  
    plt.plot(+np.sqrt(3/2)*_lambda_**-1,+ np.sqrt(3/2)*_lambda_**-1, 'ko', label='$O_4 \\to S^{xy}_+(+\\sqrt{\\frac{3}{2}}(%1.1f)^{-1},+\\sqrt{\\frac{3}{2}}(%1.1f)^{-1})$ saddle'%(_lambda_,_lambda_))       
    plt.plot(+np.sqrt(3/2)*_lambda_**-1,- np.sqrt(3/2)*_lambda_**-1, 'ko', label='$O_4 \\to S^{xy}_-(+\\sqrt{\\frac{3}{2}}(%1.1f)^{-1},-\\sqrt{\\frac{3}{2}}(%1.1f)^{-1})$ saddle'%(_lambda_,_lambda_))           
elif _lambda_==-10.0 :
    plt.plot(0,0, 'ko', label='$O_5 \\to S^{xy}_{p+}(0,0)$ saddle' ) 
    plt.plot(-1,0, 'bo', label='$O_9 \\to A^{xy}_-(-1,0)$ saddle' )
    plt.plot(+1,0, 'go', label='$O_8 \\to R^{xy}_+(+1,0)$ unstable' )  
    plt.plot(np.sqrt(3/2)*_lambda_**-1,+ np.sqrt(3/2)*_lambda_**-1, 'ko', label='$O_4 \\to S^{xy}_+(-\\sqrt{\\frac{3}{2}}(%1.1f)^{-1},+\\sqrt{\\frac{3}{2}}(%1.1f)^{-1})$ saddle'%(_lambda_,_lambda_))       
    plt.plot(np.sqrt(3/2)*_lambda_**-1,- np.sqrt(3/2)*_lambda_**-1, 'ko', label='$O_4 \\to S^{xy}_-(-\\sqrt{\\frac{3}{2}}(%1.1f)^{-1},-\\sqrt{\\frac{3}{2}}(%1.1f)^{-1})$ saddle'%(_lambda_,_lambda_))               
elif _lambda_==2 :
    plt.plot(0,0, 'go', label='$R^{xy}_0(0,0)$ unstable' )  
    plt.plot(np.sqrt(2/3),np.sqrt(1/3), 'ro', label='$A^{xy}_+(\\sqrt{\\frac{2}{3}},\\sqrt{\\frac{1}{3}})$ stable spiral' )  

elif _lambda_==3 :
    plt.plot(0,0, 'go', label='$R^{xy}_0(0,0)$ unstable' )  
    plt.plot(np.sqrt(1/6),np.sqrt(1/6), 'ro', label='$A^{xy}_+(\\sqrt{\\frac{1}{6}},\\sqrt{\\frac{1}{6}})$ stable spiral' )  

if _lambda_==10.0 or _lambda_==-10.0: plt.legend(loc='lower left',fontsize='small')
elif _lambda_==-1.0: plt.legend(loc='lower right',fontsize='large') 
else: plt.legend(loc='lower left',fontsize='large')
plt.xlabel('$x$', size=15)
plt.ylabel('$y$', size=15)
plt.title('Phase Portrait in ($x, y$) Plane, $\lambda=$%1.2f'%_lambda_)
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')
#                       ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy1_lambda_0

print('# Compute the vector field at each point, for m, x')
m, x = np.meshgrid(m_vals, x_vals)
y = 0
dm_dt, dx_dt, dy_dt = vector_fields(m, x, y)
print('# Plot the vector field')
plt.figure(5,figsize=(figsize_var, figsize_var))
plt.clf()
plt.quiver(m, x, dm_dt, dx_dt, color='blue', pivot='mid')
if _lambda_==0.0 : 
    plt.plot(0,0, 'bo', label='$O \\to S^{mx}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'go', label='$O_8 \\to R^{mx}_+(0,+1)$ unstable' )  #O_5  
    plt.plot(0,-1, 'go', label='$O_8 \\to R^{mx}_-(0,-1)$ unstable' )  #O_5
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{mx}_+(+1,0)$ saddle' )      
elif _lambda_==1.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{mx}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'go', label='$O_8 \\to R^{mx}_+(0,+1)$ unstable' )   #O_5 
    plt.plot(0,-1, 'go', label='$O_8 \\to R^{mx}_-(0,-1)$ unstable' )   #O_5
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{mx}_+(+1,0)$ saddle' )       
elif _lambda_==-1.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{mx}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'go', label='$O_8 \\to R^{mx}_+(0,+1)$ unstable' )    #O_5
    plt.plot(0,-1, 'go', label='$O_8 \\to R^{mx}_-(0,-1)$ unstable' )    #O_5
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{mx}_+(+1,0)$ saddle' )   
elif _lambda_==10.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{mx}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'go', label='$O_8 \\to R^{mx}_+(0,+1)$ unstable' )    
    plt.plot(0,-1, 'go', label='$O_8 \\to R^{mx}_-(0,-1)$ unstable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{mx}_+(+1,0)$ saddle' )  
    plt.plot(0, np.sqrt(3/2)*_lambda_**-1, 'bo', label='$O_5 \\to S^{mx}_+(0,+\\sqrt{\\frac{3}{2}}(%1.1f)^{-1})$ saddle'%(_lambda_))       
    plt.plot(0,-np.sqrt(3/2)*_lambda_**-1, 'bo', label='$O_5 \\to S^{mx}_-(0,-\\sqrt{\\frac{3}{2}}(%1.1f)^{-1})$ saddle'%(_lambda_))       
elif _lambda_==-10.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{mx}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'go', label='$O_8 \\to R^{mx}_+(0,+1)$ unstable' )    
    plt.plot(0,-1, 'go', label='$O_8 \\to R^{mx}_-(0,-1)$ unstable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{mx}_+(+1,0)$ saddle' )  
    #plt.plot( np.sqrt(3/2)*_lambda_**-1,0, 'bo', label='$S^{mx}_+(+\\sqrt{\\frac{3}{2}}(%1.1f)^{-1},0)$ saddle'%(_lambda_))       
    plt.plot(-np.sqrt(3/2)*_lambda_**-1,0, 'bo', label='$S^{mx}_-(-\\sqrt{\\frac{3}{2}}(%1.1f)^{-1},0)$ saddle'%(_lambda_))       

if _lambda_==0.0 or _lambda_==1.0 or _lambda_==-1.0 or _lambda_==10.0 or _lambda_==-10.0: 
    plt.legend(loc='lower right',fontsize='large') 
else: plt.legend(loc='lower left',fontsize='large') 
plt.xlabel('$\Omega_m$', size=15)
plt.ylabel('$x$', size=15)
plt.title('Phase Portrait in ($\Omega_m, x$) Plane, $\lambda=$%1.2f'%_lambda_)
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

print('# Compute the vector field at each point, for m, y')
m, y = np.meshgrid(m_vals, y_vals)
x = 0
dm_dt, dx_dt, dy_dt = vector_fields(m, x, y)
print('# Plot the vector field')
# Plot the vector field
plt.figure(6,figsize=(figsize_var, figsize_var))
plt.clf()
plt.quiver(m, y, dm_dt, dy_dt, color='blue', pivot='mid')
if _lambda_==0.0 : 
    plt.plot(0,0, 'bo', label='$O \\to S^{my}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'ro', label='$O_6 \\to A^{my}_+(0,+1)$ stable' )    
    plt.plot(0,-1, 'ro', label='$O_6 \\to A^{my}_-(0,-1)$ stable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{my}_+(+1,0)$ saddle' )     
elif _lambda_==1.0 :
    plt.plot(0,0, 'go', label='$O_8 \\to R^{my}_0(0,0)$ unstable' )
    #plt.plot(0,+1, 'ro', label='$A^{my}_+(0,+1)$ stable' )    
    #plt.plot(0,-1, 'ro', label='$A^{my}_-(0,-1)$ stable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{my}_+(+1,0)$ saddle' ) 
    plt.plot(0, np.sqrt(5/6)*_lambda_**-1, 'ro', label='$O_7 \\to A^{my}_+(0,+\\sqrt{\\frac{5}{6}}(%1.0f)^{-1})$ stable'%(_lambda_))       
    plt.plot(0,-np.sqrt(5/6)*_lambda_**-1, 'ro', label='$O_7 \\to A^{my}_-(0,-\\sqrt{\\frac{5}{6}}(%1.0f)^{-1})$ stable'%(_lambda_))                  
elif _lambda_==-1.0 :
    plt.plot(0,0, 'go', label='$O_8 \\to R^{my}_0(0,0)$ unstable' )
    #plt.plot(0,+1, 'ro', label='$A^{my}_+(0,+1)$ stable' )    
    #plt.plot(0,-1, 'ro', label='$A^{my}_-(0,-1)$ stable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{my}_+(+1,0)$ saddle' )  
    plt.plot(0, np.sqrt(5/6)*_lambda_**-1, 'ro', label='$O_{11} \\to A^{my}_+(0,+\\sqrt{\\frac{5}{6}}(%1.0f)^{-1})$ stable'%(_lambda_))  #O_5     
    plt.plot(0,-np.sqrt(5/6)*_lambda_**-1, 'ro', label='$O_{11} \\to A^{my}_-(0,-\\sqrt{\\frac{5}{6}}(%1.0f)^{-1})$ stable'%(_lambda_))  #O_5                     
elif _lambda_==10.0 :
    plt.plot(0,0, 'bo', label='$O \\to S^{my}_0(0,0)$ saddle' )
    plt.plot(0,+1, 'ro', label='$O_7 \\to A^{my}_+(0,+1)$ stable' )    
    plt.plot(0,-1, 'ro', label='$O_7 \\to A^{my}_-(0,-1)$ stable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{my}_+(+1,0)$ saddle' )     
elif _lambda_==-10.0 :
    plt.plot(0,0, 'go', label='$O_8 \\to R^{my}_0(0,0)$ unstable' )
    plt.plot(0,+1, 'ro', label='$O_7 \\to A^{my}_+(0,+1)$ stable' )    
    plt.plot(0,-1, 'ro', label='$O_7 \\to A^{my}_-(0,-1)$ stable' )
    plt.plot(+1,0, 'bo', label='$O_1 \\to S^{my}_+(+1,0)$ saddle' ) 
    #plt.plot( np.sqrt(3/2)*_lambda_**-1,0, 'go', label='$R^{my}_+(+\\sqrt{\\frac{3}{2}}(%1.0f)^{-1},0)$ saddle'%(_lambda_))       
    plt.plot(-np.sqrt(3/2)*_lambda_**-1,0, 'go', label='$O_4 \\to R^{my}_-(-\\sqrt{\\frac{3}{2}}(%1.0f)^{-1},0)$ saddle'%(_lambda_))                   


if _lambda_==0.0 or _lambda_==1.0 or _lambda_==-1.0 or _lambda_==10.0 or  _lambda_==-10.0: plt.legend(loc='lower right',fontsize='large')
else: plt.legend(loc='lower left',fontsize='large')
plt.xlabel('$\Omega_m$', size=15)
plt.ylabel('$y$', size=15)
plt.title('Phase Portrait in ($m, y$) Plane, $\lambda=$%1.2f'%_lambda_)
plt.grid()
plt.show()
plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')


print('# Creating phase portrait trajectories from one critical point to another')



# Time span for the integration
t_span = (-12,1)  # Adjust the end time if needed
t_eval = np.linspace(t_span[0], t_span[1], 1000)  # Points at which to evaluate the solution


if _lambda_==0:

    print('# plot the trajectories in $(x,y)$ ')
    plt.ion()
    plt.figure(4,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    initial_point = (0, -0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    initial_point = (0, -0.999, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    initial_point = (0,  0.999, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    initial_point = (0,  0.0001, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    initial_point = (0,  -0.0001, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,x)$ ')
    plt.ion()
    plt.figure(5,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0,  -0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0.0001,  0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0.0001,  -0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,y)$ ')
    plt.ion()
    plt.figure(6,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0.0001, 0, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])

    # Solve the system of ODEs
    initial_point = (-0.0001, 0, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])    

    # Solve the system of ODEs
    initial_point = (0.9999, 0, +0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])  

    # Solve the system of ODEs
    initial_point = (0.9999, 0,  -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])            
    
    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

if _lambda_==1 or _lambda_==-1:

    print('# plot the trajectories in $(x,y)$ ')
    plt.ion()
    plt.figure(4,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])  

    # Solve the system of ODEs
    initial_point = (0,  0.999, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])  

    # Solve the system of ODEs
    initial_point = (0,  -0.999, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2]) 

    # Solve the system of ODEs
    initial_point = (0,  0.999, -0.01) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2]) 

    # Solve the system of ODEs
    initial_point = (0,  +0.999, +0.01) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2]) 

    # Solve the system of ODEs
    initial_point = (0,  -0.0001, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])            

    # Solve the system of ODEs
    initial_point = (0,  +0.0001, +0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2]) 

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,x)$ ')
    plt.ion()
    plt.figure(5,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0,  -0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0.0001,  +0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1]) 

    # Solve the system of ODEs
    initial_point = (0.0001,  -0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])        

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')


    print('# plot the trajectories in $(m,y)$ ')
    plt.ion()
    plt.figure(6,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0.001,  0, 0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0.001,  0, -0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])    

    # Solve the system of ODEs
    initial_point = (+0.999,  0, +0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])  

    # Solve the system of ODEs
    initial_point = (+0.999,  0, -0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2]) 

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

if _lambda_==-1:

    print('# plot the trajectories in $(x,y)$ ')
    plt.ion()
    plt.figure(4,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -0.999, -0.01) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2]) 

    # Solve the system of ODEs
    initial_point = (0,  -0.999, +0.01) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,x)$ ')
    plt.ion()
    plt.figure(5,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,y)$ ')
    plt.ion()
    plt.figure(6,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0.001,  0, 0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[2])

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

if _lambda_==10:

    print('# plot the trajectories in $(x,y)$ ')
    plt.ion()
    plt.figure(4,figsize=(figsize_var, figsize_var))
    
    # Solve the system of ODEs
    initial_point = (0,  0, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -0.999, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  np.sqrt(3/2)/10, np.sqrt(3/2)/10) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'k-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  np.sqrt(3/2)/10, -np.sqrt(3/2)/10) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'k-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])    

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,x)$ ')
    plt.ion()
    plt.figure(5,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0,  -0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])    

    # Solve the system of ODEs
    initial_point = (0.001,  0.9999, 0.0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0.001,  -0.9999, 0.0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])   

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')


    print('# plot the trajectories in $(m,y)$ ')
    plt.ion()
    plt.figure(6,figsize=(figsize_var, figsize_var))

    def system_2(t, vars): 
        dm_dt, dx_dt, dy_dt = vector_fields(vars[0], 0.0, vars[1])
        return(dm_dt, dy_dt) 
    initial_point = (np.sqrt(3/2)/10,  -0.01) 

    initial_point = (0,  +0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])  

    initial_point = (0,  -0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])  

    initial_point = (1,  +0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])     

    initial_point = (1,  -0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])              

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

if _lambda_==-10:

    print('# plot the trajectories in $(x,y)$ ')
    plt.ion()
    plt.figure(4,figsize=(figsize_var, figsize_var))
    
    # Solve the system of ODEs
    initial_point = (0,  0, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -0.999, 0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'k-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -0.999, -0.0001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'k-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -np.sqrt(3/2)/10, np.sqrt(3/2)/10) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'k-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])

    # Solve the system of ODEs
    initial_point = (0,  -np.sqrt(3/2)/10, -np.sqrt(3/2)/10) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'k-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2]) 

    # Solve the system of ODEs
    initial_point = (0,  1.0, 0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])     

    # Solve the system of ODEs
    initial_point = (0,  1.0, -0.001) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[1], initial_point[2])            

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_xy'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,x)$ ')
    plt.ion()
    plt.figure(5,figsize=(figsize_var, figsize_var))

    # Solve the system of ODEs
    initial_point = (0,  0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0,  -0.9999, 0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])    

    # Solve the system of ODEs
    initial_point = (0.001,  0.9999, 0.0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])

    # Solve the system of ODEs
    initial_point = (0.001,  -0.9999, 0.0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])   

    # Solve the system of ODEs
    initial_point = (-np.sqrt(3/2)*_lambda_**-1,  +0.0001, 0.0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])  

    # Solve the system of ODEs
    initial_point = (-np.sqrt(3/2)*_lambda_**-1,  -0.0001, 0.0) 
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%(initial_point[0], initial_point[1], initial_point[2]))
    plt.plot(initial_point[0], initial_point[1])      

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_mx'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')

    print('# plot the trajectories in $(m,y)$ ')
    plt.ion()
    plt.figure(6,figsize=(figsize_var, figsize_var))
      

    def system_2(t, vars): 
        dm_dt, dx_dt, dy_dt = vector_fields(vars[0], 0.0, vars[1])
        return(dm_dt, dy_dt) 
    initial_point = (np.sqrt(3/2)/10,  -0.01) 

    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1]) 

    initial_point = (np.sqrt(3/2)/10,  +0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])  

    initial_point = (0,  +0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])  

    initial_point = (0,  -0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'r-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])  

    initial_point = (1,  +0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])     

    initial_point = (1,  -0.01) 
    sol = solve_ivp(system_2, t_span, initial_point, t_eval=t_eval, method='RK45')
    # Plot in the x-y plane
    plt.plot(sol.y[0], sol.y[1], 'b-', label='Trajectory from Point $(m,x,y)=$(%1.1f,0,%1.1f)'%(initial_point[0], initial_point[1] ))
    plt.plot(initial_point[0], initial_point[1])                     

    plt.savefig('./savefigs/ChatGPT_scalarLCDM_7D_3D_set_of_sim_diff_eqns_phase_portraits_my'+str(np.round(extremes,0))+'_lambda_'+str(_lambda_)+'.pdf')



"""
print('# General Creating phase portrait trajectories from one critical point to another')

# Initialize an empty array to store solutions
solutions = []



# Define the lists of critical points
list_1 = [
    (0, 0, 0),
    (1, 0, 0),
    (0,  0.999, 0.0001),
    (0, -0.999, 0.0001),    
    (0, -0.999, -0.0001), 
    (0,  0.999, -0.0001), 
    (0, 0, 1),
    (np.sqrt(1/6)+0.0001, np.sqrt(5/6)+0.0001, 0),
    (np.sqrt(1/6)-0.0001, np.sqrt(5/6)+0.0001, 0),   
    (np.sqrt(1/6)-0.0001, np.sqrt(5/6)-0.0001, 0),   
    (np.sqrt(1/6)+0.0001, np.sqrt(5/6)-0.0001, 0),
    (np.sqrt(2/3), np.sqrt(1/3), 0),
    (0.3,0,0),
    (0, 1,  1),    
    (0, 1, -1),        
]

plt.ion()
plt.figure(7,figsize=(figsize_var, figsize_var))
plt.clf()

# Solve the system for each pair of critical points in list_1

for i, initial_point in enumerate(list_1):

    # Solve the system of ODEs
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')
    solutions.append(sol)
    # Plot in the x-y plane
    plt.plot(sol.y[1], sol.y[2]) #, label='Trajectory from Point %1.0f $(m,x,y)=$(%1.1f,%1.1f,%1.1f)'%( i+1, initial_point[0], initial_point[1], initial_point[2] ))
    print(initial_point)

# Customize and show the plot
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectories in the $x-y$ Plane')
plt.legend()
plt.grid()
plt.show()

"""

print('# Create a 3D plot, so that we can see the tracking solutions easier')

from mpl_toolkits.mplot3d import Axes3D

# Time span and initial conditions
t_span = (-12, 1)  # time range
t_eval = np.linspace(t_span[0], t_span[1], 1000)  # evaluation points

# Plotting the 3D trajectory
plt.ion()
fig = plt.figure(8,figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

initial_point_list_2 = [
(0., 0., 0.),
(0.0001, 0, 0.9999),
(0.0001, 0.0001, 0.9999),
(0.0001, 0, 0.9999),
(0.0001, 0.0001, 0.9999),
(0.0001, 0, -0.9999),
(0.0001, 0.0001, -0.9999),
(0.0001, -0.0001, -0.9999),
(0.0001, -0.0001,  0.9999)
]

for i, initial_point in enumerate(initial_point_list_2): # initial conditions for (m, x, y)

    # Solve the system of ODEs
    sol = solve_ivp(system, t_span, initial_point, t_eval=t_eval, method='RK45')

    # Extract solution components
    m_values = sol.y[0]
    x_values = sol.y[1]
    y_values = sol.y[2]


    ax.plot(m_values, x_values, y_values, color='b', lw=2)

# Labels and titles
ax.set_xlabel('m')
ax.set_ylabel('x')
ax.set_zlabel('y')
ax.set_title('3D Plot of Solution Trajectory for System of Differential Equations')

# Show plot
plt.show()

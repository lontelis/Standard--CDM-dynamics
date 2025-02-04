"""
Author: Pierros Ntelis, 24 September 2024
"""

print('############ Jacobian 1, determinant ############')

import sympy as sp

# Define symbolic variables
m, k, l = sp.symbols('m k λ')
sqrt3_over_2 = sp.sqrt(3) / sp.sqrt(2)

# Define the matrix J using sp.Rational to enforce fractions
J = sp.Matrix([
    [0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0],
    [0,  0, 3, 0, 0],
    [0,  0, 0, sp.Rational(-3, 2), 0],
    [0,  0, 0, 0, sp.Rational(3, 2)],
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fraction form

# Print the solution
print(simplified_solution)


print('############ Jacobian 2 ############')

import sympy as sp

# Define the variables
m, lam, k = sp.symbols('m lambda k', real=True)
sqrt3_over_2 = sp.sqrt(3) / sp.sqrt(2)

# Define the Jacobian matrix, ensuring fractions are used
J = sp.Matrix([
    [0, m, -3*m, 0, 0],
    [0, -1, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, 0, 0, -sp.Rational(3, 2), 0],
    [0, 0, 0, 0, sp.Rational(3, 2)]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fraction form

# Print the solution
print(simplified_solution)


print('############ Jacobian 3 +, determinant ############')

import sympy as sp

# Define symbolic variables
m, k, l = sp.symbols('m k λ')
sqrt3_over_2 = sp.sqrt(3) / sp.sqrt(2)

# Define the Jacobian matrix J, ensuring fractions are used
J = sp.Matrix([
    [0, m, -3*m, 6*m*sqrt3_over_2, -6*m*sqrt3_over_2],
    [0, -1, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, sp.Rational(1, 2)*sqrt3_over_2, -sp.Rational(3, 2)*sqrt3_over_2, 3, 3*l - sp.Rational(9, 2)*sqrt3_over_2],
    [0, sp.Rational(1, 2)*sqrt3_over_2, -sp.Rational(3, 2)*sqrt3_over_2, -sp.Rational(3, 2)*l + sp.Rational(9, 2), -3 - sp.Rational(3, 2)*l]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Ensure fraction format

# Print the solution
print(simplified_solution)


print('############ Jacobian 3 -, determinant ############')

import sympy as sp

# Define symbolic variables
m, k, l = sp.symbols('m k λ')
sqrt3_over_2 = sp.sqrt(3) / sp.sqrt(2)

# Define the Jacobian matrix J, using sp.Rational to ensure fractions
J = sp.Matrix([
    [0, m, -3*m, 6*m*sqrt3_over_2, 6*m*sqrt3_over_2],
    [0, -1, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, sp.Rational(1, 2)*sqrt3_over_2, -sp.Rational(3, 2)*sqrt3_over_2, 3, -3*l + sp.Rational(9, 2)*sqrt3_over_2],
    [0, -sp.Rational(1, 2)*sqrt3_over_2, sp.Rational(3, 2)*sqrt3_over_2, sp.Rational(3, 2)*l - sp.Rational(9, 2), -3 - sp.Rational(3, 2)*l]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)


print('############ Jacobian 4 +, determinant ############')

import sympy as sp

# Define symbolic variables
m, k, l = sp.symbols('m k λ')
sqrt3_over_2 = sp.sqrt(3) / sp.sqrt(2)

# Define the Jacobian matrix J, using sp.Rational for fractional values
J = sp.Matrix([
    [0, m, -3*m, 6*m*sqrt3_over_2*l**-1, -6*m*sqrt3_over_2*l**-1],
    [0, -1, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, sp.Rational(1, 2)*sqrt3_over_2*l**-1, -sp.Rational(3, 2)*sqrt3_over_2*l**-1, 9*l**-2 - sp.Rational(3, 2), 3 - sp.Rational(9, 2)*sqrt3_over_2*l**-3],
    [0, sp.Rational(1, 2)*sqrt3_over_2*l**-1, -sp.Rational(3, 2)*sqrt3_over_2*l**-1, -sp.Rational(3, 2) + sp.Rational(9, 2)*l**-2, sp.Rational(9, 2)*l**-2]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)


print('############ Jacobian 4 -, determinant ############')

import sympy as sp

# Define symbolic variables
m, k, l = sp.symbols('m k λ')
sqrt3_over_2 = sp.sqrt(3) / sp.sqrt(2)

# Define the Jacobian matrix J, using sp.Rational for fractional values
J = sp.Matrix([
    [0, m, -3*m, 6*m*sqrt3_over_2*l**-1, 6*m*sqrt3_over_2*l**-1],
    [0, -1, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, sp.Rational(1, 2)*sqrt3_over_2*l**-1, -sp.Rational(3, 2)*sqrt3_over_2*l**-1, 9*l**-2 - sp.Rational(3, 2), - 3 + sp.Rational(9, 2)*sqrt3_over_2*l**-3],
    [0, -sp.Rational(1, 2)*sqrt3_over_2*l**-1, sp.Rational(3, 2)*sqrt3_over_2*l**-1, sp.Rational(3, 2) - sp.Rational(9, 2)*l**-2, sp.Rational(9, 2)*l**-2]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)


print('############ Jacobian 5, determinant ############')

import sympy as sp

# Define symbolic variables
r, k = sp.symbols('r k')

# Define the Jacobian matrix J, using sp.Rational for fractional entries
J = sp.Matrix([
    [r, 0, 0, 0, 0],
    [0, -1 + 2*r, -3*r, 0, 0],
    [0, 0, 3 + r, 0, 0],
    [0, 0, 0, (r - 3)*sp.Rational(1, 2), 0],
    [0, 0, 0, 0, (3 + r)*sp.Rational(1, 2)]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)


print('############ Jacobian 6 +, determinant ############')

import sympy as sp

# Define symbolic variables
r, k, l = sp.symbols('r k λ')
sqrt2_over_3 = sp.sqrt(2) / sp.sqrt(3)
sqrt3 = sp.sqrt(3)

# Define the Jacobian matrix J, using sp.Rational for fractional entries
J = sp.Matrix([
    [r - 4, 0, 0, 0, 0],
    [0, 2*r + 3, -3*r, 12*sqrt2_over_3*r, -4*sqrt3*r],
    [0, 0, r + 7, 0, 0],
    [0, sqrt2_over_3, -3*sqrt2_over_3, (r + 17)*sp.Rational(1, 2), (7 + r)*sp.Rational(1, 2) - 2*l],
    [0, sqrt3*sp.Rational(1, 3), -sqrt3, -sp.sqrt(2)*l + 4*sp.sqrt(6), (r - 3)*sp.Rational(1, 2) - 2*l]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)


print('############ Jacobian 6 -, determinant ############')

import sympy as sp

# Define symbolic variables
r, k, l = sp.symbols('r k λ')
sqrt2_over_3 = sp.sqrt(2) / sp.sqrt(3)
sqrt3 = sp.sqrt(3)

# Define the Jacobian matrix J, using sp.Rational for fractional entries
J = sp.Matrix([
    [r - 4, 0, 0, 0, 0],
    [0, 2*r + 3, -3*r, 12*sqrt2_over_3*r, 4*sqrt3*r],
    [0, 0, r + 7, 0, 0],
    [0, sqrt2_over_3, -3*sqrt2_over_3, (r + 17)*sp.Rational(1, 2), (7 + r)*sp.Rational(1, 2) - 2*l],
    [0, sqrt3*sp.Rational(-1, 3), sqrt3, sp.sqrt(2)*l + 4*sp.sqrt(6), (r - 3)*sp.Rational(1, 2) - 2*l]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)

print('############ Jacobian 7 +, determinant ############')

import sympy as sp

# Define symbolic variables
r, k, l = sp.symbols('r k λ')

# Use rational numbers and square roots without decimals
sqrt2 = sp.sqrt(2)
sqrt3 = sp.sqrt(3)
sqrt2_over_3 = sqrt2 / sp.Integer(3)

# Define the Jacobian matrix J
J = sp.Matrix([
    [r - 8*l**6 + 4*l**2, 0, 0, 0, 0],
    [0, -1 + 2*r + 8*l**6 - 4*l**2, -3*r, 12*sqrt2_over_3*r*l**3, -4*sqrt3*r*l],
    [0, 0, 3 + r + 8*l**6 - 4*l**2, 0, 0],
    [0, sqrt2_over_3*l**3, -3*sqrt2_over_3*l**3, r/2 + 12*l**6 - 2*l**2 - sp.Rational(3, 2), 
     2*sqrt2*l**2 - 16*sqrt3*l**7 / sp.Integer(3)],
    [0, sqrt3/3*l, -sqrt3*l, -sqrt2*l**2 + 4*sqrt2*l**4, 
     (3 + r)/2 + 4*l**6 - 2*l**4 - 9*l**2]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)



print('############ Jacobian 7 -, determinant ############')

import sympy as sp

# Define symbolic variables
r, k, l = sp.symbols('r k λ')

# Use rational numbers and square roots without decimals
sqrt2 = sp.sqrt(2)
sqrt3 = sp.sqrt(3)
sqrt2_over_3 = sqrt2 / sp.Integer(3)

# Define the Jacobian matrix J
J = sp.Matrix([
    [r - 8*l**6 + 4*l**2, 0, 0, 0, 0],
    [0, -1 + 2*r + 8*l**6 - 4*l**2, -3*r, 12*sqrt2_over_3*r*l**3, -4*sqrt3*r*l],
    [0, 0, 3 + r + 8*l**6 - 4*l**2, 0, 0],
    [0, sqrt2_over_3*l**3, -3*sqrt2_over_3*l**3, r/2 + 12*l**6 - 2*l**2 - sp.Rational(3, 2), 
     -2*sqrt2*l**2 + 16*sqrt3*l**7 / sp.Integer(3)],
    [0, -sqrt3/3*l, sqrt3*l, sqrt2*l**2 - 4*sqrt2*l**4, 
     (3 + r)/2 + 4*l**6 - 2*l**4 - 9*l**2]
])

# Define the identity matrix scaled by k
I = sp.eye(5) * k

# Compute the determinant of (J - k*I)
determinant = (J - I).det()

# Solve for k and ensure the result is in fractional form
solution = sp.solve(determinant, k)
simplified_solution = [sp.nsimplify(sol) for sol in solution]  # Simplify result to fractions

# Print the solution
print(simplified_solution)

print('# Polynomial Root finder on eigenvalues')

print('############ 3rd solution eigenvalues conditions explorations ############')

import sympy as sp
import numpy as np

# Define symbolic variables
r, k, l = sp.symbols('r k λ')

# Define the inequality expression
inequality_expr = 7*l**2 - 6*sp.sqrt(6)*l - 32*l + 18*sp.sqrt(6) - 16

# Solve the inequality
solution_intervals = sp.solveset(inequality_expr <= 0, l, domain=sp.S.Reals)

print(solution_intervals)

sqrt6 = sp.sqrt(6)
sqrt2 = sp.sqrt(2)

# Expression inside the square root
expr = sp.sqrt(211 - 15 * sqrt6)

# Calculate the interval bounds
lower_bound = (-sqrt2 * expr / 7 + 3 * sqrt6 / 7 + 16 / 7).evalf()
upper_bound = (sqrt2 * expr / 7 + 3 * sqrt6 / 7 + 16 / 7).evalf()

print('# Output the approximate decimal values of the interval bounds')
(lower_bound, upper_bound)
print('# The interval of $\lambda$ is:')
(round(lower_bound,3), round(upper_bound,3))

print('# solve the second inequality')

import sympy as sp
import numpy as np

# Define symbolic variables
r, k, l = sp.symbols('r k λ')

# Define the inequality expression
inequality_expr = 8*l**2 - 6*sp.sqrt(6)*l - 32*l + 18*sp.sqrt(6) - 16

# Solve the inequality
solution_intervals = sp.solveset(inequality_expr >= 0, l, domain=sp.S.Reals)

print(solution_intervals)

def sqrt(x): return(sp.sqrt(x))

sqrt6 = sp.sqrt(6)
sqrt2 = sp.sqrt(2)

# Expression inside the square root
expr = sp.sqrt(211 - 15 * sqrt6)

# Calculate the interval bounds
lower_bound = (-sqrt(6)*sqrt(73 - 8*sqrt(6))/8 + 3*sqrt(6)/8 + 2).evalf()
upper_bound = (3*sqrt(6)/8 + 2 + sqrt(6)*sqrt(73 - 8*sqrt(6))/8).evalf()

print('# Output the approximate decimal values of the interval bounds')
(lower_bound, upper_bound)
print('# The interval of $\lambda$ is:')
(round(lower_bound,3), round(upper_bound,3))



print('############ 4th solution eigenvalues conditions explorations ############')

import sympy as sp
import numpy as np

# Define symbolic variables
r, k, l = sp.symbols('r k l')

# Define the inequality expression
inequality_expr = (-7*l**10 
              + 18*l**8 
              + 6*sp.sqrt(6)*l**7 
              + 9*l**6 
              - 18*sp.sqrt(6)*l**5)

# Solve the inequality
solution_intervals = sp.solveset(inequality_expr >= 0, l, domain=sp.S.Reals)

print(solution_intervals)

def sqrt(x): return(sp.sqrt(x))

sqrt6 = sp.sqrt(6)
sqrt2 = sp.sqrt(2)

# Expression inside the square root
expr = sp.sqrt(211 - 15 * sqrt6)

# Calculate the interval bounds
lower_bound = (-1/(7*(3*sqrt(6)/7 + sqrt(2653)/49)**(1/3)) + (3*sqrt(6)/7 + sqrt(2653)/49)**(1/3)).evalf()
upper_bound = (sqrt(3)).evalf()

print('# Output the approximate decimal values of the interval bounds')
print( (lower_bound, upper_bound) )
print('# The interval of $\lambda$ is:')
print( sp.Union( sp.Interval(round(-sqrt(3)), 0) , sp.Interval(round(lower_bound,3), round(upper_bound,3)) ) )


print('# Solving ineqequality by factorising out lamda^5')

import sympy as sp

# Define the symbolic variable
lmbda = sp.symbols('lambda')

# Define the factored polynomial
polynomial = -7*lmbda**5 + 18*lmbda**3 + 6*sp.sqrt(6)*lmbda**2 + 9*lmbda - 18*sp.sqrt(6)

# Solve the polynomial equation numerically
solutions = sp.solve(polynomial, lmbda)

# Display the solutions
solutions

print('## plot ev_x and ev_y at [-\sqrt{3}, 0] \cup (1.169, \sqrt{3})')
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Define symbolic variable lambda
λ = sp.symbols('λ')

# Define f(λ) function
f_lambda = -7*λ**10 + 18*λ**8 + 6*sp.sqrt(6)*λ**7 + 9*λ**6 - 18*sp.sqrt(6)*λ**5

# Define the first eigenvalue equation
eigenvalue_eq1 = (-3*λ**5 + 27*λ**3 - 3*sp.sqrt(f_lambda)) / (4*λ**5)

# Define the second eigenvalue equation
eigenvalue_eq2 = (-3*λ**5 + 27*λ**3 + 3*sp.sqrt(f_lambda)) / (4*λ**5)

# Convert symbolic expressions into numerical functions for plotting
eigenvalue_func1 = sp.lambdify(λ, eigenvalue_eq1, 'numpy')
eigenvalue_func2 = sp.lambdify(λ, eigenvalue_eq2, 'numpy')

# Define the range for lambda (excluding zero to avoid division by zero)
λ_values = np.linspace(-5, -0.1, 200).tolist() + np.linspace(0.1, 5, 200).tolist()

# Calculate the corresponding eigenvalues for both functions
eigenvalues1 = eigenvalue_func1(np.array(λ_values))
eigenvalues2 = eigenvalue_func2(np.array(λ_values))

# Plotting the eigenvalue dependence on lambda for both equations
plt.ion()
plt.figure(figsize=(10, 6))

# Plot the first eigenvalue equation
plt.plot(λ_values, eigenvalues1, label=r'$\vec{ev}_x = \frac{-3\lambda^5 + 27\lambda^3 - 3\sqrt{f(\lambda)}}{4\lambda^5}$', color='b')

# Plot the second eigenvalue equation
plt.plot(λ_values, eigenvalues2, label=r'$\vec{ev}_y = \frac{-3\lambda^5 + 27\lambda^3 + 3\sqrt{f(\lambda)}}{4\lambda^5}$', color='r')

# Add title and labels
plt.title('Eigenvalue Dependence on λ')
plt.xlabel('λ')
plt.ylabel('Eigenvalue')
plt.grid(True)
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)

# Show the legend
plt.legend()

# Show the plot
plt.show()
plt.savefig('./ChatGPT_scalarLCDM_7D_plots/solution_4th_eigenvalue_characterisation.pdf')


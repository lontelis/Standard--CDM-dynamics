"""
Author: Pierros Ntelis, 11 August 2024
"""

import sympy as sp

print('# Define variables')
m, r, Lambda, x, y, lambda_ = sp.symbols('m r Lambda x y lambda')

print('# Define the system of differential equations')
m_prime = -3*m + m*(3 + r - 3*Lambda + 3*x**2 - 3*y**2)
r_prime = -4*r + r*(3 + r - 3*Lambda + 3*x**2 - 3*y**2)
Lambda_prime = Lambda*(3 + r - 3*Lambda + 3*x**2 - 3*y**2)
x_prime = -3/2*x + sp.sqrt(3)/sp.sqrt(2)*lambda_*y**2 + (1/2)*x*(r - 3*Lambda + 3*x**2 - 3*y**2)
y_prime = -lambda_*(sp.sqrt(3)/sp.sqrt(2))*y*x + (1/2)*y*(3 + r - 3*Lambda + 3*x**2 - 3*y**2)
lambda_prime = 0 #-sp.sqrt(6)*x*lambda_**2*(1 - 1)  # Since Gamma = 1, this term is zero

print('# Define the system')
equations = [m_prime, r_prime, Lambda_prime, x_prime, y_prime, lambda_prime]

print('# Compute the critical points by setting each derivative to zero')
critical_points = sp.solve([eq for eq in equations], (m, r, Lambda, x, y, lambda_))

print("Critical Points:")
for point in critical_points:
    print(point)

print('# Calculate the Jacobian matrix of the system')
variables = [m, r, Lambda, x, y, lambda_]
jacobian_matrix = sp.Matrix([[sp.diff(eq, var) for var in variables] for eq in equations])

print("\nJacobian Matrix:")
sp.pprint(jacobian_matrix)

print('# Evaluate the Jacobian at each critical point and compute eigenvalues and eigenvectors')
for idx, point in enumerate(critical_points):
    jacobian_at_point = jacobian_matrix.subs({m: point[0], r: point[1], Lambda: point[2], x: point[3], y: point[4], lambda_: point[5]})
    eigenvalues = jacobian_at_point.eigenvals()
    eigenvectors = jacobian_at_point.eigenvects()
    
    print(f"\nCritical Point {idx+1}: {point}")
    print("Jacobian Matrix at Critical Point:")
    sp.pprint(jacobian_at_point)
    print("Eigenvalues:")
    for eigenvalue in eigenvalues:
        print(eigenvalue)
    print("Eigenvectors:")
    for eigenvector in eigenvectors:
        print(eigenvector)

    print('# Determine stability')
    stable = all(sp.re(ev) < 0 for ev in eigenvalues)
    unstable = any(sp.re(ev) > 0 for ev in eigenvalues)
    saddle = any(sp.re(ev) > 0 for ev in eigenvalues) and any(sp.re(ev) < 0 for ev in eigenvalues)

    if stable:
        print("Classification: Stable (All eigenvalues have negative real parts)")
    elif unstable:
        print("Classification: Unstable (At least one eigenvalue has a positive real part)")
    elif saddle:
        print("Classification: Saddle Point (Mixed signs in eigenvalues)")
    else:
        print("Classification: Cannot be determined (Complex or zero eigenvalues)")

print("\nAnalysis complete.")

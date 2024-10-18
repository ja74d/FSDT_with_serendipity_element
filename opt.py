import time
import numpy as np
import sympy as sp
from sympy import Matrix
from input import *
from mesh import Nex, Ney, coordinations
from code_table import code
from gaussian_quad import RIP_Gauss



# Symbolic variables
k, e = sp.symbols('k e')

# Lagrangian Interpolation Functions
N = np.array([
    (1 / 4) * (1 - k) * (1 - e),
    (1 / 4) * (1 + k) * (1 - e),
    (1 / 4) * (1 + k) * (1 + e),
    (1 / 4) * (1 - k) * (1 + e)
])

# Material constants
db = (E * h**3) / (12 * (1 - nu**2))
ds = (E * h * ka) / (2 * (1 + nu))

# Material property matrices
Db = db * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
Ds = ds * np.array([[1, 0], [0, 1]])

# Pre-allocate Ke and Fe lists
Ke = []
Fe = []

# Pre-compute shape function derivatives
dNk = np.array([sp.diff(Ni, k) for Ni in N])
dNe = np.array([sp.diff(Ni, e) for Ni in N])

count = 0
# Loop through each element and compute stiffness and force matrices
for elemc in range(len(coordinations)):
    # Start the timer
    start_time = time.perf_counter()
    element_coordinates = np.array(coordinations[elemc])

    count += 1
    print(count)

    # Symbolic Jacobian matrix
    J = np.array([
        [sum(dNk[i] * element_coordinates[i][0] for i in range(4)),
         sum(dNk[i] * element_coordinates[i][1] for i in range(4))],
        [sum(dNe[i] * element_coordinates[i][0] for i in range(4)),
         sum(dNe[i] * element_coordinates[i][1] for i in range(4))]
    ], dtype=object)

    # Determinant and inverse of Jacobian
    det_J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    J_inv = np.array([[J[1, 1], -J[0, 1]], [-J[1, 0], J[0, 0]]]) / det_J

    # Precompute shape function derivatives in global coordinates
    DNx = J_inv[0, 0] * dNk + J_inv[0, 1] * dNe
    DNy = J_inv[1, 0] * dNk + J_inv[1, 1] * dNe

    # Bending strain-displacement matrix (Bb)
    Bb = np.zeros((3, 12), dtype=object)
    for i in range(4):
        Bb[:, i * 3 + 2] = [DNx[i], -DNy[i], DNy[i] - DNx[i]]

    # Shear strain-displacement matrix (Bs)
    Bs = np.zeros((2, 12), dtype=object)
    for i in range(4):
        Bs[:, i * 3] = [DNx[i], DNy[i]]
        Bs[0, i * 3 + 2] = N[i]
        Bs[1, i * 3 + 1] = -N[i]

    # Compute element stiffness matrix (Ke)
    Gb = Bb.T @ Db @ Bb
    Gs = Bs.T @ Ds @ Bs
    K_eb = np.zeros((12, 12))
    K_es = np.zeros((12, 12))

    for o in range(12):
        for p in range(12):
            K_eb[o, p] = RIP_Gauss(Gb[o, p] * det_J)
            K_es[o, p] = RIP_Gauss(Gs[o, p] * det_J)

    Ke.append(K_eb + K_es)

    # Compute element force vector (Fe)
    Nw = np.array([N[0], 0, 0, N[1], 0, 0, N[2], 0, 0, N[3], 0, 0])
    F_e = np.zeros((12, 1))
    for i in range(12):
        F_e[i, 0] = RIP_Gauss(p0 * Nw[i] * det_J)
    Fe.append(F_e)
    # End the timer
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Time taken: {elapsed_time} seconds")


# Assembling global stiffness matrix and force vector
num_dofs = np.max(code)
K = np.zeros((num_dofs, num_dofs))
F = np.zeros(num_dofs)

for elem in range(code.shape[0]):
    for i in range(12):
        if code[elem, i] != 0:
            for j in range(12):
                if code[elem, j] != 0:
                    K[code[elem, i] - 1, code[elem, j] - 1] += Ke[elem][i, j]
            F[code[elem, i] - 1] += Fe[elem][i, 0]

# Solving for displacement
Delta = np.linalg.solve(K, F)

# Finding displacement at midpoint
Wmid = np.max(Delta)
wmidND = Wmid / (p0 * Lx**4 / db)

# Output the result
print(f"number of elements: {code.shape[0]}")
print(f"displacement at midpoint: {Wmid}")
print(f"Non-dimensional displacement at midpoint: {wmidND}")


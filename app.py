import time
import numpy as np
import sympy as sp
from tqdm import tqdm
from sympy import Matrix
from input import *
from mesh import coordinations
from code_table import code
from gaussian_quad import RIP_Gauss
#from scipy.linalg import lu_factor, lu_solve
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

# Start the timer
start = time.perf_counter()

#Symboles
k, e = sp.symbols('k e')
#k for kesi
#e for eta

#Lagrangian Interpolation Functions
N1 = (1/4)*(1-k)*(1+e)*(-k+e-1)
N2 = (1/4)*(1+k)*(1+e)*(k+e-1)
N3 = (1/4)*(1+k)*(1-e)*(k-e-1)
N4 = (1/4)*(1-k)*(1-e)*(-k-e-1)

N5 = (1/2)*(1-k**2)*(1+e)
N7 = (1/2)*(1-k**2)*(1-e)

N6 = (1/2)*(1+k)*(1-e**2)
N8 = (1/2)*(1-k)*(1-e**2)

N = np.array([N1, N2, N3, N4, N5, N6, N7, N8])
Nw = np.array([N[0], 0, 0, N[1], 0, 0, N[2], 0, 0, N[3], 0, 0, N[4], 0, 0, N[5], 0, 0, N[6], 0, 0, N[7], 0, 0])

def DNx(N):
    return (J_inv[0,0])*( sp.diff(N, k) ) + (J_inv[0,1])*( sp.diff(N, e) )

def DNy(N):
    return (J_inv[1,0])*( sp.diff(N, k) ) + (J_inv[1,1])*( sp.diff(N, e) )

def matrix_in_list(matrix_to_check, matrix_list):
    matrix_to_check = (matrix_to_check) 

#D
#Db matrix(Bending)
db = (E*h**3)/(12*(1-(nu**2)))
Db = db*np.array([
    [1, nu, 0],
    [nu, 1, 0],
    [0, 0, (1-nu)/2]
])

#Ds matrix(shearing)
ds = (E*h*(ka))/(2*(1+nu))
Ds = ds*np.array([
    [1, 0],
    [0, 1]
])

Jacob = []
Ke = []
Fe = []
count = 0


for elemc in tqdm(range(len(coordinations)),desc="Calculating elements"):
    # Start the timer
    #start_time = time.perf_counter()
    element_coordinates = coordinations[elemc]
    count += 1
    n_elem = len(code)
    #if (count + 1) % (n_elem // 20) == 0:
    #    percentage_complete = (count + 1) / n_elem * 100
    #    print(f"{percentage_complete:.1f}% complete")


    #Saving Private Jacobian
    #symbolic in jacobian
    J = Matrix([
        [0, 0],
        [0, 0]
    ])

    j00 = 0
    for i in range(8):
        j00 += sp.diff(N[i], k) * element_coordinates[i][0]
    J[0, 0] = j00

    j01 = 0
    for i in range(8):
        j01 += sp.diff(N[i], k) * element_coordinates[i][1]
    J[0, 1] = j01

    j10 = 0
    for i in range(8):
        j10 += sp.diff(N[i], e) * element_coordinates[i][0]
    J[1, 0] = j10

    j11 = 0
    for i in range(8):
        j11 += sp.diff(N[i], e) * element_coordinates[i][1]
    J[1, 1] = j11

    #Det J
    det_J = (J[0,0])*(J[1,1]) - (J[0,1])*(J[1,0])
    for ij in Jacob:
        ij_check = np.zeros((2, 2)); J_check = np.zeros((2,2))
        J_check[0,0] = RIP_Gauss(J[0,0],3); J_check[0,1] = RIP_Gauss(J[0,1],3)
        J_check[1,0] = RIP_Gauss(J[1,0],3); J_check[1,1] = RIP_Gauss(J[1,1],3)
        ij_check[0,0] = RIP_Gauss(ij[0,0],3); ij_check[0,1] = RIP_Gauss(ij[0,1],3)
        ij_check[1,0] = RIP_Gauss(ij[1,0],3); ij_check[1,1] = RIP_Gauss(ij[1,1],3)
        
        if (np.linalg.norm(ij_check) - np.linalg.norm(J_check)) < tol:
            K_e = Ke[Jacob.index(ij)]
            Ke.append(K_e)
            break
    else:
        Jacob.append(J)

        #inverse of jacoian
        #J_inv = np.linalg.inv(J)
        J_inv = (1/det_J)*np.array([
            [J[1,1], -J[0,1]],
            [-J[0,1], J[0,0]]
        ])


        #B
        #Bb matrix(Bending)
        Bb = []
        for i in range(8):

            BB = np.array([

            [0, 0, DNx(N[i])],
            [0, -DNy(N[i]), 0],
            [0, -DNx(N[i]), DNy(N[i])]

            ])
            Bb.append(BB)
        Bb = np.array(Bb)


        #Bs matrix(Shearing)
        Bs = []
        for i in range(8):
            BS = np.array([

            [DNx(N[i]), 0, N[i]],
            [DNy(N[i]), -N[i], 0],

            ])
            Bs.append(BS)
        Bs = np.array(Bs)



        #Gaussian Integration Method *
        #Ke

        def compute_kb(i, j):
            return Bb[i].T @ Db @ Bb[j]


        Gb = np.block([[compute_kb(i, j) for j in range(8)] for i in range(8)])


        K_eb = np.zeros((12, 12))
        for o in range(0, 12):
            for p in range(0, 12):
               K_eb[o, p] = RIP_Gauss(Gb[o, p]*det_J)


        def calculate_ks(i, j):
            return Bs[i].T @ Ds @ Bs[j]


        Gs = np.block([[calculate_ks(i, j) for j in range(8)] for i in range(8)])


        K_eb = np.zeros((24, 24))
        for o in range(0, 24):
            for p in range(0, 24):
               #K_eb[o, p] = sp.integrate(sp.integrate((Gb[o, p]*det_J), (k, -1, 1)), (e, -1, 1))
               K_eb[o, p] = RIP_Gauss(Gb[o, p]*det_J)


        Gs = np.array(Gs)
        K_es = np.zeros((24, 24))
        for o in range(0, 24):
            for p in range(0, 24):
                #K_es[o, p] = sp.integrate(sp.integrate((Gs[o, p]*det_J), (k, -1, 1)), (e, -1, 1))
                K_es[o, p] = RIP_Gauss(Gs[o, p]*det_J)

        K_e = K_es + K_eb
        Ke.append(K_e)

        #Fe
    F_e = np.zeros((24, 1))
    for i in range(0,24):
        #F_e[i,0] = sp.integrate(sp.integrate( (p0*Nw[i]*det_J), (k,-1,1) ),(e,-1,1))
        F_e[i, 0] = RIP_Gauss(p0*Nw[i]*det_J)
    Fe.append(F_e)
    
    # End the timer
    #end_time = time.perf_counter()
    # Calculate elapsed time
    #elapsed_time = end_time - start_time
    #print(f"Time taken: {elapsed_time} seconds")


#Assembling
num_dofs = np.max(code)

K = np.zeros((num_dofs, num_dofs))
F = np.zeros(num_dofs)

num_elements = code.shape[0]

for elem in range(num_elements):
    for i in range(24):
        if code[elem, i] != 0:
            for j in range(24):
                if code[elem, j] != 0:
                    K[code[elem, i] - 1, code[elem, j] - 1] += Ke[elem][i, j]
            F[code[elem, i] - 1] += Fe[elem][i, 0]

K_sparse = csc_matrix(K)
Delta = spsolve(K_sparse, F)
#lu, piv = lu_factor(K)
#Delta = lu_solve((lu, piv), F)
#Delta = np.linalg.inv(K) @ F
#Delta = K**-1 @ F

Wmid = max(Delta)
wmidND = (Wmid / (p0 * (Lx)**4 / db))

# Output the result
print(f"number of elements: {num_elements}")
print(f"displacement at midpoint: {Wmid}")
print(f"Non-dimensional displacement at midpoint: {wmidND}")  

# End the timer
end = time.perf_counter()

# Calculate elapsed time
elapsed_time = end - start
print(f"Time taken: {elapsed_time} seconds")

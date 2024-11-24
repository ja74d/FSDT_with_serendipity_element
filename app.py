import time
import sys
import h5py
import numpy as np
import sympy as sp
from tqdm import tqdm
from interpolation_functions import *
from sympy import Matrix
from input import *
from mesh import coordinations
from code_table import code
from gaussian_quad import RIP_Gauss, apply_rip_gauss
#from scipy.linalg import lu_factor, lu_solve
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import eig
from isotropic_A_B_D_S import A, D, S, db, Db, Ds
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigsh
#from Single_FGM import A, D, S, db
#sys.path.append('/home/javad/FSDT_with_serendipity_element/Sandwich_FSDT')
#from Sandwich_FSDT.app import A, B, D, S, db

# Start the timer
start = time.perf_counter()

#S1 = ((np.pi**2)*db)/Lx**2
S1 = 1
Sigma = np.array([
    [S1, 0],
    [0, 0]
])

def DNx(i):
    return (J_inv[0,0])*( dNdk[i] ) + (J_inv[0,1])*( dNde[i] )

def DNy(i):
    return (J_inv[1,0])*( dNdk[i] ) + (J_inv[1,1])*( dNde[i] )


def Bm_function(i):
    return np.array([
        [0, 0, 0, DNx(i), 0],
        [0, 0, 0, 0, DNy(i)],
        [0, 0, 0, DNy(i), DNx(i)]])

def Bb_function(i):
    return np.array([
        [0, DNx(i), 0, 0, 0],
        [0, 0, DNy(i), 0, 0],
        [0,  DNy(i), DNx(i), 0, 0]])

def Bs_function(i):
    return np.array([
        [DNx(i), N[i], 0, 0, 0],
        [DNy(i), 0, N[i], 0, 0],])

def BGb_function(i):
    return  np.array([
            [DNx(i), 0, 0, 0, 0],
            [DNy(i), 0, 0, 0, 0]])

def Bgs1_function(i):
    return np.array([
            [0, 0, DNx(i), 0, 0],
            [0, 0, DNy(i), 0, 0]])

def Bgs2_function(i):
    return np.array([
            [0, DNx(i), 0, 0, 0],
            [0, DNy(i), 0, 0, 0]])

Jacob = []
Ke = []
Kge = []
Fe = []
count = 0


for elemc in tqdm(range(len(coordinations)),desc="Calculating elements"):
    element_coordinates = coordinations[elemc]
    count += 1
    n_elem = len(code)

    #Saving Private Jacobian
    #symbolic in jacobian
    J = Matrix([
        [0, 0],
        [0, 0]
    ])

    j00 = 0
    j01 = 0
    j10 = 0
    j11 = 0
    for i in range(8):
        j00 += dNdk[i] * float(element_coordinates[i][0])
        j01 += dNdk[i] * float(element_coordinates[i][1])
        j10 += dNde[i] * float(element_coordinates[i][0])
        j11 += dNde[i] * float(element_coordinates[i][1])
    J[0, 0] = j00
    J[0, 1] = j01
    J[1, 0] = j10
    J[1, 1] = j11

    #Det J
    det_J = (J[0,0])*(J[1,1]) - (J[0,1])*(J[1,0])
    if Jacob_cache == 1:
        pass
    #for ij in Jacob:
    #    J_check = apply_rip_gauss(J, 3)
    #    ij_check = apply_rip_gauss(ij, 3)
    #    
    #    if abs((np.linalg.norm(ij_check) - np.linalg.norm(J_check))) < tol:
    #        K_e = Ke[Jacob.index(ij)]
    #        K_eg = Kge[Jacob.index(ij)]
    #        Ke.append(K_e)
    #        Kge.append(K_eg)
    #        break
    
    else:
        Jacob.append(J)

        #inverse of jacoian
        #J_inv = np.linalg.inv(J)
        J_inv = (1/det_J)*np.array([
            [J[1,1], -J[0,1]],
            [-J[1,0], J[0,0]]
        ])

        #B
        #Bm matrix(membrane)
        Bm = [Bm_function(0), Bm_function(1), Bm_function(2), Bm_function(3), Bm_function(4), Bm_function(5), Bm_function(6), Bm_function(7)]

        #Bb matrix(Bending)
        Bb = [Bb_function(0), Bb_function(1), Bb_function(2), Bb_function(3), Bb_function(4), Bb_function(5), Bb_function(6), Bb_function(7)]

        #Bs matrix(Shearing)
        Bs = [Bs_function(0), Bs_function(1), Bs_function(2), Bs_function(3), Bs_function(4), Bs_function(5), Bs_function(6), Bs_function(7)]
        
        #Bg matrix(Geometric)
        Bgb = [BGb_function(0), BGb_function(1), BGb_function(2), BGb_function(3), BGb_function(4), BGb_function(5), BGb_function(6), BGb_function(7)]

        Bgs1 = [Bgs1_function(0), Bgs1_function(1), Bgs1_function(2), Bgs1_function(3), Bgs1_function(4), Bgs1_function(5), Bgs1_function(6), Bgs1_function(7)]

        Bgs2 = [Bgs2_function(0), Bgs2_function(1), Bgs2_function(2), Bgs2_function(3), Bgs2_function(4), Bgs2_function(5), Bgs2_function(6), Bgs2_function(7)]


        def calculate_kg(i, j):
            return h*( Bgb[i].T @ Sigma @ Bgb[j] ) + ((h**3)/12)*(Bgs1[i].T @ Sigma @ Bgs1[j])  + (((h**3)/12)*(Bgs2[i].T @ Sigma @ Bgs2[j]))
        Pg = np.block([[calculate_kg(i, j) for j in range(8)] for i in range(8)])

        k_eg = np.zeros((40, 40))
        for o in range(0, 40):
            for p in range(0, 40):
                k_eg[o, p] = RIP_Gauss(Pg[o, p]*det_J)
        Kge.append(k_eg)

        #Gaussian Integration Method *
        #Ke

        def compute_km(i, j):
            return Bm[i].T @ A @ Bm[j]
        Gm = np.block([[compute_km(i, j) for j in range(8)] for i in range(8)])

        K_em = np.zeros((40, 40))
        for o in range(0, 40):
            for p in range(0, 40):
                K_em[o, p] = RIP_Gauss(Gm[o, p]*det_J)

        #K_eb
        def compute_kb(i, j):
            return Bb[i].T @ D @ Bb[j]
        Gb = np.block([[compute_kb(i, j) for j in range(8)] for i in range(8)])

        K_eb = np.zeros((40, 40))
        for o in range(0, 40):
            for p in range(0, 40):
               K_eb[o, p] = RIP_Gauss(Gb[o, p]*det_J)

        def calculate_ks(i, j):
            return Bs[i].T @ S @ Bs[j]
        Gs = np.block([[calculate_ks(i, j) for j in range(8)] for i in range(8)])

        K_es = np.zeros((40, 40))
        for o in range(0, 40):
            for p in range(0, 40):
                K_es[o, p] = RIP_Gauss(Gs[o, p]*det_J)


        def compute_kb(i,j):
            return Bb[i].T @ D @ Bb[j]

        Gb = np.block([[compute_kb(i, j) for j in range(8)] for i in range(8)])

        K_e = K_es + K_eb + K_em
        Ke.append(K_e)

        #Fe
    F_e = np.zeros((40, 1))
    for i in range(0,40):
        #F_e[i,0] = sp.integrate(sp.integrate( (p0*Nw[i]*det_J), (k,-1,1) ),(e,-1,1))
        F_e[i, 0] = RIP_Gauss(p0*Nw[i]*det_J)
    Fe.append(F_e)

#Assembling
num_dofs = np.max(code)

K = np.zeros((num_dofs, num_dofs))
Kg = np.zeros((num_dofs, num_dofs))
F = np.zeros(num_dofs)

num_elements = code.shape[0]

for elem in range(num_elements):
    for i in range(40):
        if code[elem, i] != 0:
            for j in range(40):
                if code[elem, j] != 0:
                    K[code[elem, i] - 1, code[elem, j] - 1] += Ke[elem][i, j]
                    Kg[code[elem, i] - 1, code[elem, j] - 1] += Kge[elem][i, j]
            F[code[elem, i] - 1] += Fe[elem][i, 0]

#db = (Em*h**3)/(12*(1-(nu**2)))

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
print(f"Non-dimensional displacement at midpoint: {100*wmidND}")

# End the timer
end = time.perf_counter()

# Calculate elapsed time
elapsed_time = end - start
print(f"Time taken: {elapsed_time} seconds")

#Saving the Results
with h5py.File(f'{file_name}.h5', 'w') as hdf:

    # Create a group for matrices
    matrices_group = hdf.create_group('Matrices')
    matrices_group.create_dataset('Stiffness_Matrix', data=K)
    matrices_group.create_dataset('Force_Matrix', data=F)
    matrices_group.create_dataset('Geometrical_Matrix', data=Kg)

    # Create a group for nodal data
    nodal_group = hdf.create_group('Nodal_Data')
    nodal_group.create_dataset('Nodal_Displacement', data=Delta)

    # You can also add attributes for metadata
    hdf.attrs['Description'] = 'Finite Element Analysis Results'
    hdf.attrs['Version'] = '1.0'

print(f"Data has been saved in '{file_name}.h5'")
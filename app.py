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
from scipy.linalg import eig
#from isotropic_A_B_D_S import A, D, S, db, Db, Ds
from Single_FGM import A, D, S, db

# Start the timer
start = time.perf_counter()

#Symboles
k, e = sp.symbols('k e')
#k for kesi
#e for eta

S1 = ((np.pi**2)*db)/Lx**2

Sigma = np.array([
    [S1, 0],
    [0, 0]
])

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
Nw = np.array([N[0], 0, 0, 0, 0, N[1], 0, 0, 0, 0, N[2], 0, 0, 0, 0, N[3], 0, 0,0 , 0, N[4], 0, 0, 0 ,0 , N[5], 0, 0, 0, 0, N[6], 0, 0,0 , 0, N[7], 0, 0, 0, 0])

def DNx(N):
    return (J_inv[0,0])*( sp.diff(N, k) ) + (J_inv[0,1])*( sp.diff(N, e) )

def DNy(N):
    return (J_inv[1,0])*( sp.diff(N, k) ) + (J_inv[1,1])*( sp.diff(N, e) )

Jacob = []
Ke = []
Kge = []
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
            #K_eg = Kge[Jacob.index(ij)]
            Ke.append(K_e)
            #Kge.append(K_eg)
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
        Bm = []
        for i in range(8):
            BM = np.array([
                [0, 0, 0, DNx(N[i]), 0],
                [0, 0, 0, 0, DNy(N[i])],
                [0, 0, 0, DNy(N[i]), DNx(N[i])]
            ])
            Bm.append(BM)

        #Bb matrix(Bending)
        Bb = []
        for i in range(8):
            BB = np.array([
            [0, DNx(N[i]), 0, 0, 0],
            [0, 0, DNy(N[i]), 0, 0],
            [0,  DNy(N[i]), DNx(N[i]), 0, 0]
            ])
            Bb.append(BB)

        #Bs matrix(Shearing)
        Bs = []
        for i in range(8):
            BS = np.array([
            [DNx(N[i]), N[i], 0, 0, 0],
            [DNy(N[i]), 0, N[i], 0 , 0],
            ])
            Bs.append(BS)

        Bgb = []
        for i in range(8):
            BGB = np.array([
                [DNx(N[i]), 0, 0, 0, 0],
                [DNy(N[i]), 0, 0, 0, 0]
            ])
            Bgb.append(BGB)

        Bgs1 = []
        for i in range(8):
            BGS1 = np.array([
                [0, 0, DNx(N[i]), 0, 0],
                [0, 0, DNy(N[i]), 0, 0]
            ])
            Bgs1.append(BGS1)

        Bgs2 = []
        for i in range(8):
            BGS2 = np.array([
                [0, DNx(N[i]), 0, 0, 0],
                [0, DNy(N[i]), 0, 0, 0]
            ])
            Bgs2.append(BGS2)

        def calculate_kg(i, j):
            return h*( Bgb[i].T @ Sigma @ Bgb[j] ) + ((h**3)/12)*(Bgs1[i].T @ Sigma @ Bgs1[j])  + (((h**3)/12)*(Bgs2[i].T @ Sigma @ Bgs2[j]))
        #Pg = np.block([[calculate_kg(i, j) for j in range(8)] for i in range(8)])

        #k_eg = np.zeros((40, 40))
        #for o in range(0, 40):
        #    for p in range(0, 40):
        #        k_eg[o, p] = RIP_Gauss(Pg[o, p]*det_J)
        #Kge.append(k_eg)

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
    
    # End the timer
    #end_time = time.perf_counter()
    # Calculate elapsed time
    #elapsed_time = end_time - start_time
    #print(f"Time taken: {elapsed_time} seconds")


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
                    #Kg[code[elem, i] - 1, code[elem, j] - 1] += Kge[elem][i, j]
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
#eigenvalues, eigenvectors = eig(K, Kg)
#print(f"Buckling load factor: {min(eigenvalues)}")

# End the timer
end = time.perf_counter()

# Calculate elapsed time
elapsed_time = end - start
print(f"Time taken: {elapsed_time} seconds")

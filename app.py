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

        kb11 = (Bb[0].T@Db@Bb[0])
        kb12 = (Bb[0].T@Db@Bb[1])
        kb13 = (Bb[0].T@Db@Bb[2])
        kb14 = (Bb[0].T@Db@Bb[3])
        kb15 = (Bb[0].T@Db@Bb[4])
        kb16 = (Bb[0].T@Db@Bb[5])
        kb17 = (Bb[0].T@Db@Bb[6])
        kb18 = (Bb[0].T@Db@Bb[7])

        kb21 = (Bb[1].T@Db@Bb[0])
        kb22 = (Bb[1].T@Db@Bb[1])
        kb23 = (Bb[1].T@Db@Bb[2])
        kb24 = (Bb[1].T@Db@Bb[3])
        kb25 = (Bb[1].T@Db@Bb[4])
        kb26 = (Bb[1].T@Db@Bb[5])
        kb27 = (Bb[1].T@Db@Bb[6])
        kb28 = (Bb[1].T@Db@Bb[7])

        kb31 = (Bb[2].T@Db@Bb[0])
        kb32 = (Bb[2].T@Db@Bb[1])
        kb33 = (Bb[2].T@Db@Bb[2])
        kb34 = (Bb[2].T@Db@Bb[3])
        kb35 = (Bb[2].T@Db@Bb[4])
        kb36 = (Bb[2].T@Db@Bb[5])
        kb37 = (Bb[2].T@Db@Bb[6])
        kb38 = (Bb[2].T@Db@Bb[7])

        kb41 = (Bb[3].T@Db@Bb[0])
        kb42 = (Bb[3].T@Db@Bb[1])
        kb43 = (Bb[3].T@Db@Bb[2])
        kb44 = (Bb[3].T@Db@Bb[3])
        kb45 = (Bb[3].T@Db@Bb[4])
        kb46 = (Bb[3].T@Db@Bb[5])
        kb47 = (Bb[3].T@Db@Bb[6])
        kb48 = (Bb[3].T@Db@Bb[7])

        kb51 = (Bb[4].T@Db@Bb[0])
        kb52 = (Bb[4].T@Db@Bb[1])
        kb53 = (Bb[4].T@Db@Bb[2])
        kb54 = (Bb[4].T@Db@Bb[3])
        kb55 = (Bb[4].T@Db@Bb[4])
        kb56 = (Bb[4].T@Db@Bb[5])
        kb57 = (Bb[4].T@Db@Bb[6])
        kb58 = (Bb[4].T@Db@Bb[7])

        kb61 = (Bb[5].T@Db@Bb[0])
        kb62 = (Bb[5].T@Db@Bb[1])
        kb63 = (Bb[5].T@Db@Bb[2])
        kb64 = (Bb[5].T@Db@Bb[3])
        kb65 = (Bb[5].T@Db@Bb[4])
        kb66 = (Bb[5].T@Db@Bb[5])
        kb67 = (Bb[5].T@Db@Bb[6])
        kb68 = (Bb[5].T@Db@Bb[7])
        
        kb71 = (Bb[6].T@Db@Bb[0])
        kb72 = (Bb[6].T@Db@Bb[1])
        kb73 = (Bb[6].T@Db@Bb[2])
        kb74 = (Bb[6].T@Db@Bb[3])
        kb75 = (Bb[6].T@Db@Bb[4])
        kb76 = (Bb[6].T@Db@Bb[5])
        kb77 = (Bb[6].T@Db@Bb[6])
        kb78 = (Bb[6].T@Db@Bb[7])

        kb81 = (Bb[7].T@Db@Bb[0])
        kb82 = (Bb[7].T@Db@Bb[1])
        kb83 = (Bb[7].T@Db@Bb[2])
        kb84 = (Bb[7].T@Db@Bb[3])
        kb85 = (Bb[7].T@Db@Bb[4])
        kb86 = (Bb[7].T@Db@Bb[5])
        kb87 = (Bb[7].T@Db@Bb[6])
        kb88 = (Bb[7].T@Db@Bb[7])


        Gb = np.block([
            [kb11, kb12, kb13, kb14, kb15, kb16, kb17, kb18],
            [kb21, kb22, kb23, kb24, kb25, kb26, kb27, kb28],
            [kb31, kb32, kb33, kb34, kb35, kb36, kb37, kb38],
            [kb41, kb42, kb43, kb44, kb45, kb46, kb47, kb48],
            [kb51, kb52, kb53, kb54, kb55, kb56, kb57, kb58],
            [kb61, kb62, kb63, kb64, kb65, kb66, kb67, kb68],
            [kb71, kb72, kb73, kb74, kb75, kb76, kb77, kb78],
            [kb81, kb82, kb83, kb84, kb85, kb86, kb87, kb88],
        ])


        K_eb = np.zeros((24, 24))
        for o in range(0, 24):
            for p in range(0, 24):
               #K_eb[o, p] = sp.integrate(sp.integrate((Gb[o, p]*det_J), (k, -1, 1)), (e, -1, 1))
               K_eb[o, p] = RIP_Gauss(Gb[o, p]*det_J)

        ks11 = Bs[0].T@Ds@Bs[0]
        ks12 = Bs[0].T@Ds@Bs[1]
        ks13 = Bs[0].T@Ds@Bs[2]
        ks14 = Bs[0].T@Ds@Bs[3]
        ks15 = Bs[0].T@Ds@Bs[4]
        ks16 = Bs[0].T@Ds@Bs[5]
        ks17 = Bs[0].T@Ds@Bs[6]
        ks18 = Bs[0].T@Ds@Bs[7]

        ks21 = Bs[1].T@Ds@Bs[0]
        ks22 = Bs[1].T@Ds@Bs[1]
        ks23 = Bs[1].T@Ds@Bs[2]
        ks24 = Bs[1].T@Ds@Bs[3]
        ks25 = Bs[1].T@Ds@Bs[4]
        ks26 = Bs[1].T@Ds@Bs[5]
        ks27 = Bs[1].T@Ds@Bs[6]
        ks28 = Bs[1].T@Ds@Bs[7]

        ks31 = Bs[2].T@Ds@Bs[0]
        ks32 = Bs[2].T@Ds@Bs[1]
        ks33 = Bs[2].T@Ds@Bs[2]
        ks34 = Bs[2].T@Ds@Bs[3]
        ks35 = Bs[2].T@Ds@Bs[4]
        ks36 = Bs[2].T@Ds@Bs[5]
        ks37 = Bs[2].T@Ds@Bs[6]
        ks38 = Bs[2].T@Ds@Bs[7]

        ks41 = Bs[3].T@Ds@Bs[0]
        ks42 = Bs[3].T@Ds@Bs[1]
        ks43 = Bs[3].T@Ds@Bs[2]
        ks44 = Bs[3].T@Ds@Bs[3]
        ks45 = Bs[3].T@Ds@Bs[4]
        ks46 = Bs[3].T@Ds@Bs[5]
        ks47 = Bs[3].T@Ds@Bs[6]
        ks48 = Bs[3].T@Ds@Bs[7]

        ks51 = Bs[4].T@Ds@Bs[0]
        ks52 = Bs[4].T@Ds@Bs[1]
        ks53 = Bs[4].T@Ds@Bs[2]
        ks54 = Bs[4].T@Ds@Bs[3]
        ks55 = Bs[4].T@Ds@Bs[4]
        ks56 = Bs[4].T@Ds@Bs[5]
        ks57 = Bs[4].T@Ds@Bs[6]
        ks58 = Bs[4].T@Ds@Bs[7]

        ks61 = Bs[5].T@Ds@Bs[0]
        ks62 = Bs[5].T@Ds@Bs[1]
        ks63 = Bs[5].T@Ds@Bs[2]
        ks64 = Bs[5].T@Ds@Bs[3]
        ks65 = Bs[5].T@Ds@Bs[4]
        ks66 = Bs[5].T@Ds@Bs[5]
        ks67 = Bs[5].T@Ds@Bs[6]
        ks68 = Bs[5].T@Ds@Bs[7]

        ks71 = Bs[6].T@Ds@Bs[0]
        ks72 = Bs[6].T@Ds@Bs[1]
        ks73 = Bs[6].T@Ds@Bs[2]
        ks74 = Bs[6].T@Ds@Bs[3]
        ks75 = Bs[6].T@Ds@Bs[4]
        ks76 = Bs[6].T@Ds@Bs[5]
        ks77 = Bs[6].T@Ds@Bs[6]
        ks78 = Bs[6].T@Ds@Bs[7]

        ks81 = Bs[7].T@Ds@Bs[0]
        ks82 = Bs[7].T@Ds@Bs[1]
        ks83 = Bs[7].T@Ds@Bs[2]
        ks84 = Bs[7].T@Ds@Bs[3]
        ks85 = Bs[7].T@Ds@Bs[4]
        ks86 = Bs[7].T@Ds@Bs[5]
        ks87 = Bs[7].T@Ds@Bs[6]
        ks88 = Bs[7].T@Ds@Bs[7]

        Gs = np.block([
            [ks11, ks12, ks13, ks14, ks15, ks16, ks17, ks18],
            [ks21, ks22, ks23, ks24, ks25, ks26, ks27, ks28],
            [ks31, ks32, ks33, ks34, ks35, ks36, ks37, ks38],
            [ks41, ks42, ks43, ks44, ks45, ks46, ks47, ks48],
            [ks51, ks52, ks53, ks54, ks55, ks56, ks57, ks58],
            [ks61, ks62, ks63, ks64, ks65, ks66, ks67, ks68],
            [ks71, ks72, ks73, ks74, ks75, ks76, ks77, ks78],
            [ks81, ks82, ks83, ks84, ks85, ks86, ks87, ks88]
        ])

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

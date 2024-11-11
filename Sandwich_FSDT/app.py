import numpy as np
from sympy import Matrix, integrate, symbols
from Honeycomb import Ceramic_only
from FGM import E_z_top, E_z_bottom, h0, h1, h2, h3, nu, Em , h

z = symbols('z')

db = (Em*h**3)/(12*(1-(nu**2)))

#Locations of top and bottom of each layer of the laminate
#h = np.array([h0, h1, h2, h3])
#h=1/3
ka = 5/6
ka_FGM = 5/6

Shorhan = Ceramic_only
Ec =Shorhan['E1']

#isotropic test
#Sandwich = np.array([ [E1, E2, nu, G12, G13, G23],           #FGM_Facesheet_bottom
#                      [E1, E2, nu, G12, G13, G23],           #Honeycomb
#                      [E1, E2, nu, G12, G13, G23] ])         #FGM_Facesheet_top


#Reduced stiffness matrix
def Q_(E_z, nu):
    Q = Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ])

    Q[0,0] = Q[1,1] = E_z/(1-nu**2)
    Q[0,1] = Q[1,0] = (nu*E_z)/(1-nu**2)
    Q[2,2] = E_z/(2*(1+nu))

    return Q

def Qs_(E_z, nu):
    Qs = np.zeros((2,2))
    Qs = Matrix([
        [Q_(E_z, nu)[2,2], 0],
        [0, Q_(E_z, nu)[2,2]]
    ])
    return Qs

#print(Q(Sandwich[0]))

#A
A = np.zeros((3, 3))
A2 = np.zeros((3, 3))

for i in range(3):
    for j in range(3):
        A[i,j] = -1*( integrate( Q_(E_z_bottom, nu)[i, j], (z, h0, h1) ) + integrate( Q_(Ec, nu)[i, j], (z, h1, h2) ) + integrate( Q_(E_z_top, nu)[i, j], (z, h2, h3) ) )
        #A2[i,j] = integrate( Q_(Ec, nu)[i, j], (z, h0, h1) ) + integrate( Q_(Ec, nu)[i, j], (z, h1, h2) ) + integrate( Q_(Ec, nu)[i, j], (z, h2, h3) )

#A_prime = np.linalg.inv(A)

#B
B = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        B[i,j] = integrate( Q_(E_z_bottom, nu)[i, j]*z, (z, h0, h1) ) + integrate( Q_(Ec, nu)[i, j]*z, (z, h1, h2) ) + integrate( Q_(E_z_top, nu)[i, j]*z, (z, h2, h3) )

#D
D = np.zeros((3, 3))
D2 = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        D[i,j] = -1*( integrate( Q_(E_z_bottom, nu)[i, j]*z**2, (z, h0, h1) ) + integrate( Q_(Ec, nu)[i, j]*z**2, (z, h1, h2) ) + integrate( Q_(E_z_top, nu)[i, j]*z**2, (z, h2, h3) ) )
        #D2[i,j] = integrate( Q_(Ec, nu)[i, j]*z**2, (z, h0, h1) ) + integrate( Q_(Ec, nu)[i, j]*z**2, (z, h1, h2) ) + integrate( Q_(Ec, nu)[i, j]*z**2, (z, h2, h3) )

#D_prime = np.linalg.inv(D)
#print(D)

#S
S = np.zeros((2, 2))
S2 = np.zeros((2, 2))
for i in range(2):
    for j in range(2):
        S[i,j] = -1 * (ka_FGM * integrate( Qs_(E_z_bottom, nu)[i,j], (z, h0, h1) ) + ka * integrate( Qs_(Ec, nu)[i,j], (z, h1, h2) ) + ka_FGM * integrate( Qs_(E_z_top, nu)[i,j], (z, h2, h3) ) )
        #S2[i,j] = ka * ( integrate( Qs_(Ec, nu)[i,j], (z, h0, h1) ) + integrate( Qs_(Ec, nu)[i,j], (z, h1, h2) ) + integrate( Qs_(Ec, nu)[i,j], (z, h2, h3) ) )



#Tests

#print('A:\b',A)
#print()
#print('A2\b', A2)
#print('B:\b', B)
#print()
#print('D:\b',D)
#print()
#print('D2\b', D2)
#print('S:\b',S)
#print()
#print('S2\b', S2)

import math
from sympy import Matrix
import numpy as np
import scipy as sp
from input import h
from sympy import symbols, integrate

#volume fraction index "n"
n = 0

#Ceramic-Aluminum
Ec = 200e+09
Em = 70e+09
nu = 0.3

Gc = Ec/(2*(1+nu))
Gm = Em/(2*(1+nu))

z = symbols('z')

def powerlaw(n, pc, pm):
    return (pc-pm)*(((2*z+h)/(2*h))**n) + pm

E_z = powerlaw(n, Ec, Em)
#e = integrate(E_z*z, (z, -h/2, h/2))/integrate(E_z, (z, -h/2, h/2))
e = 0

Q = Matrix([
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
])

Q[0,0] = Q[1,1] = E_z/(1-nu**2)
Q[0,1] = Q[1,0] = (nu*E_z)/(1-nu**2)
Q[2,2] = E_z/(2*(1+nu))

Qs = Matrix([
    [Q[2,2], 0],
    [0, Q[2,2]]
])

A = np.zeros((3, 3))
D = np.zeros((3, 3))
S = np.zeros((2, 2))

for i in range(3):
    for j in range(3):
        A[i, j] = integrate(Q[i,j], (z, -h/2, h/2))
#A[0, 0] = integrate(Q[0, 0], (z, -h/2, h/2))

for ii in range(3):
    for jj in range(3):
        D[ii,jj] = integrate((Q[ii,jj]*((z-e)**2)), (z, -h/2, h/2))

for iii in range(2):
    for jjj in range(2):
        S[iii,jjj] = 5/6 * integrate(Qs[iii,jjj], (z, -h/2, h/2))

db = (Em*h**3)/(12*(1-(nu**2)))
#E = integrate(E_z, (z, -h/2, h/2))/h

#nu_z = powerlaw(1, nuc, num)
#nu = integrate(nu_z, (z, -h/2, h/2))/h

#G_z = powerlaw(1, Gc, Gm)
#G = integrate(G_z, (z, -h/2, h/2))/h

import sympy as sp
import numpy as np

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

dNdk = [ (-(0.25 - 0.25*k)*(e + 1) - 0.25*(e + 1)*(e - k - 1)), ((e + 1)*(0.25*k + 0.25) + 0.25*(e + 1)*(e + k - 1)), ((1 - e)*(0.25*k + 0.25) + 0.25*(1 - e)*(-e + k - 1)), (-(0.25 - 0.25*k)*(1 - e) - 0.25*(1 - e)*(-e - k - 1)), (-1.0*k*(e + 1)), (0.5 - 0.5*e**2), (-1.0*k*(1 - e)), (0.5*e**2 - 0.5) ]
dNde = [ ((0.25 - 0.25*k)*(e + 1) + (0.25 - 0.25*k)*(e - k - 1)), ((e + 1)*(0.25*k + 0.25) + (0.25*k + 0.25)*(e + k - 1)), (-(1 - e)*(0.25*k + 0.25) - (0.25*k + 0.25)*(-e + k - 1)), (-(0.25 - 0.25*k)*(1 - e) - (0.25 - 0.25*k)*(-e - k - 1)), (0.5 - 0.5*k**2), (-2*e*(0.5*k + 0.5)), (0.5*k**2 - 0.5), (-2*e*(0.5 - 0.5*k))]

N = np.array([N1, N2, N3, N4, N5, N6, N7, N8])
Nw = np.array([N[0], 0, 0, 0, 0, N[1], 0, 0, 0, 0, N[2], 0, 0, 0, 0, N[3], 0, 0,0 , 0, N[4], 0, 0, 0 ,0 , N[5], 0, 0, 0, 0, N[6], 0, 0,0 , 0, N[7], 0, 0, 0, 0])
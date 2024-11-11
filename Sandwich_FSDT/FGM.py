import numpy as np
from sympy import symbols, integrate

#volume fraction index
n = 0

t1, t2, t3 = 1, 1, 1

h_ = t1 + t2 + t3

t1 = t1/h_
t2 = t2/h_
t3 = t3/h_

h0 = (t1+t2+t3)/2
h1 = h0 - t1
h2 = h1 - t2
h3 = h2 - t3

h = h0 + h1 + h2 + h3

#==========bottom--1.5-metal
#                                                  |-z
#++++++++++core             -----------------------|
#                                                  |+z
#========== top-1.5-metal

#Ceramic-Aluminum
Ec = 151e+09
Em = 70e+09
nu = 0.3

#Linear Temperature Change Through the Thickness
def T(T0, deltaT):
    #return (Tbottom + ( (Ttop - Tbottom)/h )*z)
    return T0+deltaT
#Temperature at top
T = T(300, 600)

#Tempreture Dependent Material Properties
def TDP(P0, P_1, P1, P2, P3, T):
    P = P0*( P_1*(T**-1) +1 + P1*T + P2*(T**2) +P3*(T**3) )
    return P

z = symbols('z')

#Simple Power Law
#if h0<h<h1 ===> Vn = ((z-h0)/(h1-h0))^n
#if h2<h<h3 ===> Vn = ((z-h3)/(h2-h3))^n

def powerlaw_bottom(n, pc, pm):
    Vbottom = ((z-h0)/(h1-h0))**n
    return (pc-pm)*Vbottom + pm

E_z_bottom = powerlaw_bottom(n, Ec, Em)

#FGM_bottom = {'E':E_z_bottom, 'nu':nu}

def powerlaw_top(n, pc, pm):
    Vtop = ((z-h3)/(h2-h3))**n
    return (pc-pm)*Vtop + pm

E_z_top = powerlaw_top(n, Ec, Em)

FGM_top = {'E':E_z_top, 'nu':nu}


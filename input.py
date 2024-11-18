#File Name
file_name = 'test file'

#B.C.
BCleft, BCright, BCtop, BCbottom = 'S', 'Sym', 'S', 'Sym'

#volume fraction index "n"
n = 1

#Mechanincal Properties
#Ceramic-Aluminum
Ec = 380e+09
Em = 70e+09

E = 1
h = 1
nu = 0.3
K_dictionary = {
    0: 5/6,
    1: 0.831,
    2: 0.7949,
    4: 0.76,
    8: 0.7619,
    10: 0.7694
}

ka = K_dictionary[n]

#Distributed Load
p0 = 1

#Geometry
Lx = Ly = 10

#d
d = (E*h**3)/(12*(1-nu**2))

#tolerance
tol = 1e-4

#Jacobian cache status
Jacob_cache = 0
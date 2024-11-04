#B.C.
BCleft, BCright, BCtop, BCbottom = 'S', 'S', 'S', 'S'

#volume fraction index "n"
n = 0

#Mechanincal Properties
#Ceramic-Aluminum
Ec = 151e+09
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
Lx = 10

#tolerance
tol = 1e-10
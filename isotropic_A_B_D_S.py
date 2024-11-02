from input import *
import numpy as np

Q = np.zeros((3, 3))

Q[0,0] = Q[1,1] = E/(1-nu**2)
Q[0,1] = Q[1,0] = (nu*E)/(1-nu**2)
Q[2,2] = E/(2*(1+nu))

A = h*Q

B = np.zeros((3, 3))

D = (h**3/12) * Q

S = h*(ka)*np.array([
    [Q[2,2], 0],
    [0, Q[2,2]]
])

#D
#Db matrix(Bending)
db = (E*h**3)/(12*(1-(nu**2)))
Db = db*np.array([
    [1, nu, 0],
    [nu, 1, 0],
    [0, 0, (1-nu)/2]
])

#Ds matrix(shearing)
ds = (E*h*ka)/(2*(1+nu))
Ds = ds*np.array([
    [1, 0],
    [0, 1]
])

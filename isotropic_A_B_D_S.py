from input import *
import numpy as np

Q = np.zeros((3, 3))

Q[0,0] = Q[1,1] = E/(1-nu**2)
Q[0,1] = Q[1,0] = (nu*E)/(1-nu**2)
Q[2,2] = E/(2*(1+nu))

<<<<<<< HEAD
A = h*Q
=======
A = h*np.array([
    [Q[0,0], Q[0,1], Q[0,2]],
    [Q[1,0], Q[1,1], Q[1,2]],
    [Q[2,0], Q[2,1], Q[2,2]]
])
>>>>>>> 869aa3cc66a713b9baae50c1acfa04331fe3c93f

B = np.zeros((3, 3))

D = (h**3/12) * Q

<<<<<<< HEAD
S = h*(ka)*np.array([
=======
S = h*(5/6)*np.array([
>>>>>>> 869aa3cc66a713b9baae50c1acfa04331fe3c93f
    [Q[2,2], 0],
    [0, Q[2,2]]
])

#D
#Db matrix(Bending)
<<<<<<< HEAD
=======
#E = 70e+09
>>>>>>> 869aa3cc66a713b9baae50c1acfa04331fe3c93f
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
<<<<<<< HEAD
])
=======
])
>>>>>>> 869aa3cc66a713b9baae50c1acfa04331fe3c93f

import h5py
import numpy as np
from input import file_name, tol, d, Lx, p0

with h5py.File('test.h5', 'r') as hdf:
    # Access datasets
    K = hdf['Matrices/Stiffness_Matrix'][:]
    Kg = hdf['Matrices/Geometrical_Matrix'][:]
    F = hdf['Matrices/Force_Matrix'][:]
    Delta = hdf['Nodal_Data/Nodal_Displacement'][:]

#print(max(np.linalg.inv(K)@F) / (p0 * (Lx**4) / d))
#print(K)

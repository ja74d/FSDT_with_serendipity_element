import numpy as np
from input import *
from mesh import Nex, Ney, elements, left_dofs, right_dofs, top_dofs, bottom_dofs, num_nodes
from collections import Counter

#num_nodes = (Nex+1)*(Ney+1)
node_dofs = {}

#This loop gives every node, its DOFs without considering its Boundary Conditions
for node in range(1, num_nodes+1):
    start_dof = (node-1)*5+1
    node_dofs[node] = list(range(start_dof, start_dof+5))

#code = np.zeros((Nex*Ney, 12), dtype=int)

code = np.zeros((len(elements), 40), dtype=int)

# Fill the code table
for i, element in enumerate(elements):
    element_dofs = []
    for node in element:
        element_dofs.extend(node_dofs[node])
    code[i, :] = element_dofs
#print(code)

#EFFECT OF BOUNDARY CONDITIONS ON DOFs
resL = []
if BCleft == 'S':
    for i in left_dofs:
        resL.append(i[0])
        resL.append(i[2])
        resL.append(i[3])
        resL.append(i[4])
elif BCleft == 'C':
    for i in left_dofs:
        resL.append(i[0])
        resL.append(i[1])
        resL.append(i[2])
        resL.append(i[3])
        resL.append(i[4])
elif BCleft == 'F':
    pass

resR = []
if BCright == 'S':
    for i in right_dofs:
        resR.append(i[0])
        resR.append(i[2])
        resR.append(i[3])
elif BCright == 'C':
    for i in right_dofs:
        resR.append(i[0])
        resR.append(i[1])
        resR.append(i[2])
        resR.append(i[3])
        resR.append(i[4])
elif BCright == 'Sym':
    for i in right_dofs:
        #resR.append(i[0])
        resR.append(i[1])
        #resR.append(i[2])
        resR.append(i[3])
        #resR.append(i[4])
elif BCright == 'F':
    pass

resT = []
if BCtop == 'S':
    for i in top_dofs:
        resT.append(i[0])
        resT.append(i[1])
        resT.append(i[4])
elif BCtop == 'C':
    for i in top_dofs:
        resT.append(i[0])
        resT.append(i[1])
        resT.append(i[2])
        resT.append(i[3])
        resT.append(i[4])
elif BCtop == 'F':
    pass

resB = []
if BCbottom == 'S':
    for i in bottom_dofs:
        resB.append(i[0])
        resB.append(i[1])
        resB.append(i[4])
elif BCbottom == 'C':
    for i in bottom_dofs:
        resB.append(i[0])
        resB.append(i[1])
        resB.append(i[2])
        resB.append(i[3])
        resB.append(i[4])
elif BCbottom == 'Sym':
    for i in bottom_dofs:
        #resB.append(i[0])
        #resB.append(i[1])
        resB.append(i[2])
        #resB.append(i[3])
        resB.append(i[4])
elif BCbottom == 'F':
    pass

#converting lists to numpy arrays
resT = np.array(resT)
resR = np.array(resR)
resB = np.array(resB)
resL = np.array(resL)

res = np.sort(np.concatenate([resL, resT, resR, resB]))
res = np.unique(res)

sizeres = res.size

for k in range(sizeres - 1, -1, -1):
    for j in range(len(elements)):
        for i in range(40):
            if code[j, i] == res[k]:
                code[j, i] = 0
            elif code[j, i] > res[k]:
                code[j, i] -= 1
#print()

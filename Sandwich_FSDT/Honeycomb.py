import numpy as np
import math

#Section Properties

b = 1
l2 = 2
l = 1
t = 0.1
s = 0.1
teta = -30
teta = math.radians(teta)

Lx = 2*l*(math.cos(teta))
Ly = l2 + 2*l*(math.sin(teta))

alpha = l2/l
beta = s/l
gamma = b/l

#material properties
E = 70e+09
nu = 0.3
G = E/(2*(1+nu))

#correction Coefficients
k1 = 1/( 1 + beta**2*( 2.4+1.5*nu*(math.tan(teta))**-2 ) )
k2 = 1/( 1 + beta**2*( 2.4 + 1.5*nu + math.tan(teta)**2 + ((2*alpha)/(math.cos(teta)**2)) ) )
k12 = 1/ (  1 + 2*alpha + beta**2 *( ( ((2.4+1.5*nu)/(alpha))*(2+alpha+math.sin(teta)) ) + ( ((alpha+math.sin(teta))/(alpha**2))*(((2+alpha)*(math.tan(teta)**2))+(math.sin(teta))) ) ) )
c12 = k1*( 1 + beta**2*( 1.4 + 1.5*nu ) )
c21 = k2*( 1 + beta**2*( 1.4 + 1.5*nu ) )
G23u = G*(beta)*( (alpha+2*(math.sin(teta)**2))/(2*(alpha+math.sin(teta))*(math.cos(teta))) )
G23l = G*(beta)*( (alpha+math.sin(teta))/((1+2*alpha)*(math.cos(teta))) )
if teta>=0:
    alpha_=0.787
elif teta<0:
    alpha_=1.342

#in-plane mechanical properties

###gibson
#in-plan
E1 = E*((t/l)**3)*( (math.cos(teta))/((l2/l+math.sin(teta))*(math.sin(teta)**2)) )
E2 = ((t/l)**3)*E*( ( l2/l + math.sin(teta) )/(math.cos(teta)**3) )
G12 = ((t/l)**3)*E*( (l2/l+math.sin(teta))/( ((l2/l)**2)*(1+2*(l2/l))*(math.cos(teta)) ) )
nu12 = ( (math.cos(teta)**2)/( (l2/l+math.sin(teta))*(math.sin(teta)) ) )
nu21 = ( ((l2/l + math.sin(teta))*(math.sin(teta)))/(math.cos(teta)**2) )

#out of plan
G13 = (t/l)*G*( (math.cos(teta))/(l2/l + math.sin(teta)) )
G23 = (t/(2*l))*G*( ( l2/l + 2*(math.sin(teta))**2 )/((l2/l + math.sin(teta))*math.cos(teta)) )


Gibson = { 
    'E1': E1,
    'E2': E2,
    'G12': G12,
    'G13': G13,
    'G23': G23,
    'nu12': nu12,
    'nu21': nu21
         }

###shorhan
#in-plan
E1 = k1*E*(beta**3)*( (math.cos(teta))/((alpha+math.sin(teta))*(math.sin(teta)**2)) )
E2 = k2*E*(beta**3)*( (alpha+math.sin(teta))/(math.cos(teta)**3) )
G12 = k12*E*(beta**3)*( (alpha+math.sin(teta))/( (alpha**2)*(1+2*alpha)*(math.cos(teta)) ) )
nu12 = c12*( (math.cos(teta)**2)/( (alpha+math.sin(teta))*(math.sin(teta)) ) )
nu21 = c21*( ((alpha+math.sin(teta))*(math.sin(teta)))/(math.cos(teta)**2) )

#out of plan
G13 = G*(beta)*( (math.cos(teta))/(alpha+math.sin(teta)) )
G23 = G23l + (alpha_/gamma)*(G23u-G23l)


Shorhan = { 
    'E1': E1,
    'E2': E2,
    'G12': G12,
    'G13': G13,
    'G23': G23,
    'nu12': nu12,
    'nu21': nu21
         }


Ec = 151e+09
nuc = 0.3
Gc = Ec/(2*(1+nuc))

#print(Gibson)
#print(Shorhan)
Ceramic_only = {
    'E1': Ec,
    'E2': Ec,
    'G12': Gc,
    'G13': Gc,
    'G23': Gc,
    'nu12': nuc,
    'nu21': nuc
}
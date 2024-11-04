import numpy as np
import matplotlib.pyplot as plt
from FGM import h0, h1, h2, h3

#h0 = 0
#h1 = 1
#
#h2 = 2
#h3 = 3

#n_values = np.linspace(0.1, 7, 10)
n_values = [0.1, 1, 6, 3, 0.4]

plt.figure(figsize=(10, 6))

z1 = np.linspace(h0, h1, 100)
z2 = np.linspace(h1, h2, 100)
z3 = np.linspace(h2, h3, 100)


for n in n_values:
    vc1 = ((z1 - h0) / (h1 - h0)) ** n
    vc2 = np.ones_like(z2)
    vc3 = ( (z3-h3) / (h2-h3) )**n
    plt.plot(vc1, z1, label=f'n = {n}')
    plt.plot(vc2, z2, label=f'n = {n}')
    plt.plot(vc3, z3, label=f'n = {n}')

plt.title('Volume Fraction V(z) for Different n Values')
plt.xlabel('Thickness (z)')
plt.ylabel('Volume Fraction')
plt.axhline(0, color='black', linewidth=0.5, ls='--')  # X-axis
plt.axvline(0, color='black', linewidth=0.5, ls='--')  # Y-axis
plt.grid()
plt.xlim(-0.1, 1.1)  # Set y limits to show volume fractions
plt.legend()
plt.show()

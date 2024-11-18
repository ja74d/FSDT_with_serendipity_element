import sympy as sp
import numpy as np

k, e = sp.symbols('k e')

def gaussian_quadrature(n):
    x = sp.Symbol('x')
    # Get the nth Legendre polynomial
    Pn = sp.legendre(n, x)
    
    # Compute the roots of the Legendre polynomial (abscissas)
    abscissas = sp.solvers.solve(Pn, x)
    
    # Compute the weights
    weights = []
    for xi in abscissas:
        wi = 2 / ((1 - xi**2) * (sp.diff(Pn, x).subs(x, xi))**2)
        weights.append(wi)
    
    # Evaluate the symbolic expressions numerically
    abscissas = [sp.N(xi) for xi in abscissas]
    weights = [sp.N(wi) for wi in weights]
    
    return abscissas, weights


def RIP_Gauss(fn, n=3):
    if fn == 0:
        integ = 0
    else:
        #try:
        #    deg_k = sp.Poly(fn, k).degree()
        #    deg_e = sp.Poly(fn, e).degree()
        #    n = math.ceil(max(deg_k, deg_e) +1 / 2)
        #except:
            #n = 5
        if n==5:
            abscissas, weights =[0, -0.538469310105683, 0.538469310105683, -0.906179845938664, 0.906179845938664],[0.568888888888889, 0.478628670499366, 0.478628670499366, 0.236926885056189, 0.236926885056189]
        elif n==3:
            abscissas, weights = [0, -0.774596669241483, 0.774596669241483], [0.888888888888889, 0.555555555555556, 0.555555555555556]
        else:
            abscissas, weights = gaussian_quadrature(n)
        integ = 0
        for i in range(len(abscissas)):
           for j in range(len(abscissas)):
                f_val = fn.subs( {k:abscissas[i], e:abscissas[j]} )
                integ += weights[i]*weights[j]*f_val
    return integ

def apply_rip_gauss(matrix, n):
    # this function applies the gauss numerical integration over all entities of a matrix
    result = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            result[i, j] = RIP_Gauss(matrix[i, j], n)
    return result
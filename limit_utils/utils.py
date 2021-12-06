from math import e
import numpy as np

def branching_ratio(x):
    '''Computes branching ratio from xsec dataframe.'''
    gbgm_ratio = x.gb/x.gmu
    return 2/3*(1 + gbgm_ratio**2 * (1 + 2*x.dbs**2))**-1

def total_xsec(x):
    '''Computes xsec from xsec dataframe with branching ratio.'''
    return x.xsec/branching_ratio(x)

# formula to fit for c1 and k per mass point 
def mm_xsec_func(x, c1, k):
    '''Formula for xsec from gb, dbs, gmu'''
    gb = x.gb
    return c1 * gb**2 * (1 + k * x.dbs**2) * x.gmu**2
# formula to fit for c1 and k per mass point 
def total_xsec_func(x, c1, k):
    '''Formula for xsec from gb, dbs'''
    gb = x.gb
    return c1 * gb**2 * (1 + k * x.dbs**2) 


def power_law(x, c, *params):
    x = x.astype('float')
    y = x*0.
    for i, par in enumerate(params):
        y += np.power(x,i)*par
    return c*np.exp(y)

def linear(x, c1, c2):
    return c1*x + c2

def exponet_500(x,c1,c2,c3,c4,c5,c6,c7):
    '''Custom fit function to fit the c1 as a function of mass curve'''
    x = x/500
    return np.exp(c2*x+c3*x**-1+c4*x**-2+c5*x**-3)*(c1)

def compute_c1(x, popt_c1): return exponet_500(x, *popt_c1)
def compute_k(x, popt_k): return linear(x, *popt_k)

def compute_expected_xsec(x, popt_c1, popt_k):
    '''compute xsec using c1 and k from params'''
    mass = x.mass
    c1, k  = compute_c1(mass, popt_c1), compute_k(mass, popt_k)
    return total_xsec_func(x, c1, k)


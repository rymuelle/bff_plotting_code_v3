from glob import glob
import os
import uproot
import numpy as np


def log_norm_np(x, norm, sigma, theta, mean):
    return norm/((x-theta)*sigma*2*3.14159)*np.exp(-(np.log((x-theta)/mean))**2/(2*sigma**2))

def log_norm_unp(x, norm, sigma, theta, mean):
    import uncertainties.unumpy as unp
    return norm/((x-theta)*sigma*2*3.14159)*unp.exp(-(unp.log((x-theta)/mean))**2/(2*sigma**2))

def parabola(x,offset,m1,m2,m3):
    x = x - offset
    return m1+m2*x+m3*x**2

def power_func(x, c, p):
    return c*x**p

def power_law(x, c, *params):
    x = x.astype('float')
    y = x*0.
    for i, par in enumerate(params):
        y += np.power(x,i)*par
    return c*np.exp(y)

def heaviside(x, cutoff, scale):
    return np.heaviside(x-cutoff, 0)*scale

def significance(sig,bck):
    return sig/(sig+bck+1e-12)**.5

def chiSquared(unc2,unc1,dof=0):
    deltasquared = (unc1-unc2)**2
    return np.sum(deltasquared/(unc2) )/(len(unc2)-dof)   

def constant(x, b):
    return linear(x, b, 0)

def linear(x, m, b):
    return m * x + b

linear_old=linear

#def linear_old(x,m1,m2):
#    return m1*x+m2

def quad(x, b, m, m2):
    return m2 * x ** 2 + m * x + b

def double_crystalball(x, norm, mean, sigma, b1, m1, b2, m2):
    from scipy.stats import crystalball
    xprime = (x-mean)/sigma
    
    left_x = xprime[xprime<0]
    left_y = crystalball.pdf(left_x, b1, m1)
    #even out normalization
    left_cdf = crystalball.cdf(0, b1, m1)
    scale_from_norm = (1-left_cdf)/.5
    left_y=left_y/scale_from_norm
    new_left_norm = left_cdf/scale_from_norm
    
    right_x = -xprime[xprime>=0]
    right_y = crystalball.pdf(right_x, b2, m2)
    #even out normalization
    right_cdf = crystalball.cdf(0, b2, m2)
    scale_from_norm = (1-right_cdf)/.5
    right_y=right_y/scale_from_norm    
    new_right_norm = right_cdf/scale_from_norm
    
    #renormalize 
    new_norm = new_left_norm+new_right_norm
    return np.concatenate([left_y,right_y])/new_norm*2


def make_bpoly(x, *constants, x_range=[105,900]):
    from scipy.interpolate import BPoly
    constants = [[c] for c in constants]
    bp = BPoly(constants, x_range )
    return bp(x)


def make_bpoly_exp(x, *args, **kwargs):
    ''' for fitting a bernstein poly on log scale'''
    return np.exp(make_bpoly(x, *args, **kwargs))

def lognorm(x, norm, sigma, mean, theta, module=np):
    import math
    pi = math.pi
    ln = 1./((x-theta)*sigma*2*pi)*module.exp(-(module.log((x-theta)/mean))**2/(2*sigma**2))
    ln = ln/ln.sum()
    return norm*ln


### limit plotting

def branching_ratio(x):
    '''Computes branching ratio from xsec dataframe.'''
    gbgm_ratio = x.gb/x.gmu
    return 2/3*(1 + gbgm_ratio**2 * (1 + 2*x.dbs**2))**-1
    
def total_xsec(x):
    '''Computes xsec from xsec dataframe with branching ratio.'''
    return x.xsec/branching_ratio(x)

# formula to fit for c1 and k per mass point 
def total_xsec_func(x, c1, k):
    '''Formula for xsec from gb, dbs'''
    gb = x.gb
    return c1 * gb**2 * (1 + k * x.dbs**2) 


import pickle
with open('fits/limit_setting/popt_c1.pkl', 'rb') as f:
    popt_c1 = pickle.load(f)
with open('fits/limit_setting/popt_k.pkl', 'rb') as f:
    popt_k = pickle.load(f)
    



import pickle
with open('fits/limit_setting/popt_c1.pkl', 'rb') as f:
    popt_c1 = pickle.load(f)
with open('fits/limit_setting/popt_k.pkl', 'rb') as f:
    popt_k = pickle.load(f)
    



#y = 2/3/(1+x/m^2 * (1+2*s^2)) * c * x * (1+k*s^2) solve for x
# y = xsec, m = gmu, s = dbs, gb**2 = x
def gb_from_xsec_gmu_dbs_br(xsec, df):
    c1, k = compute_c1(df.mass, popt_c1), compute_k(df.mass, popt_k)
    return ((3 * df.gmu ** 2 * xsec)/(2 * c1 * df.gmu ** 2 * (1 + k * df.dbs**2) - 3 * (1 + 2 * df.dbs**2)*xsec) )**.5

# formula to fit for c1 and k per mass point to mumu br xsec
def mumu_xsec_func(x, c1, k):
    '''Formula for xsec from gb, dbs'''
    return c1 * x.gb**2 * (1 + k * x.dbs**2) * 2 * x.gmu**2

def gb_from_mm_xsec(tx, dbs, gmu, c1, k):
    '''Formula for gb from mumu br xsec and dbs'''
    return (tx/(c1 * (1 + k * dbs**2) *2 * gmu ** 2 ))**.5

# formula to fit for c1 and k per mass point 
def gb_from_total_xsec(tx, dbs, c1, k):
    '''Formula for gb from total xsec and dbs'''
    return (tx/(c1 * (1 + k * dbs**2)))**.5


def exponet_500(x,c1,c2,c3,c4,c5,c6,c7):
    '''Custom fit function to fit the c1 as a function of mass curve'''
    x = x/500
    return np.exp(c2*x+c3*x**-1+c4*x**-2+c5*x**-3)*(c1)


# not needed really, but making a linear fit to mass and gmu values from paper values to show validity of min gmu value for mass
def make_mass_to_gmu():
    from scipy.optimize import curve_fit
    gmu = [0.08, 0.14, 0.20]
    mass = np.array([200,350,500])
    popt_mass_gmu ,_ = curve_fit(linear, mass, gmu)
    return lambda x: linear(x, *popt_mass_gmu)

mass_to_gmu = make_mass_to_gmu()

def compute_c1(x, popt_c1): return exponet_500(x, *popt_c1)
def compute_k(x, popt_k): return linear(x, *popt_k)
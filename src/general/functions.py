from glob import glob
import os
import uproot
import numpy as np


def log_norm(x, norm, sigma, theta, mean):
    return norm/((x-theta)*sigma*2*3.14159)*np.exp(-(np.log((x-theta)/mean))**2/(2*sigma**2))

def log_norm_unp(x, norm, sigma, theta, mean):
    import uncertainties.unumpy as unp
    return norm/((x-theta)*sigma*2*3.14159)*unp.exp(-(unp.log((x-theta)/mean))**2/(2*sigma**2))

def parabola(x,offset,m1,m2,m3):
    x = x - offset
    return m1+m2*x+m3*x**2

def power_func(x, c, p):
    return c*x**p

def heaviside(x, cutoff, scale):
    return np.heaviside(x-cutoff, 0)*scale

def significance(sig,bck):
    return sig/(sig+bck+1e-12)**.5

def chiSquared(unc2,unc1,dof=0):
    deltasquared = (unc1-unc2)**2
    return np.sum(deltasquared/(unc2) )/(len(unc2)-dof)   

def constant(x, b):
    return linear(x, b, 0)

def linear(x, b, m):
    return m * x + b

def linear_old(x,m1,m2):
    return m1*x+m2

def quad(x, b, m, m2):
    return m2 * x ** 2 + m * x + b
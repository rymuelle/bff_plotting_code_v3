import numpy as np

def poly(x, *params):
    x = np.asarray(x)
    y = np.full(len(x),0)
    for i,par in enumerate(params):
        y += par*x**i
    return y
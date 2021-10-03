import numpy as np

def array_center(arr):
	return np.array([(arr[i]+arr[i+1])/2 for i in range(len(arr[:-1]))])

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).sum(-1).sum(1)

def chi_squared(x, x1, error, nDOF=0):
     return np.sum((x-x1)**2/(error)**2)/(len(x)-nDOF)

def rebin_simul(x,y, min_val=10):
    rebinned_x = []
    rebinned_y = []
    tempx=0
    tempy=0
    for valuex, valuey in zip(x,y):
        tempx += valuex
        tempy +=valuey
        if tempx > min_val: 
            rebinned_x.append(tempx)
            rebinned_y.append(tempy)
            tempx = 0
            tempy = 0
    if tempx > 0: 
        rebinned_x.append(tempx)
        rebinned_y.append(tempy)
    return np.array(rebinned_x), np.array(rebinned_y)


def hist_chi2(hist1, hist2, nDOF=0, systematics=False):
    hist1_rebin, hist2_rebin = rebin_simul(hist1.nominal, hist2.nominal, min_val=20)
    _, hist1_var_rebin = rebin_simul(hist1.nominal, hist1.std**2, min_val=20)
    if systematics:
        _, hist1_sys_var_rebin = rebin_simul(hist1.nominal, hist1.up**2, min_val=20)
        hist1_var_rebin+=hist1_sys_var_rebin
    return chi_squared(hist2_rebin, hist1_rebin, hist1_var_rebin**.5, nDOF=nDOF)
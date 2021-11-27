from glob import glob
import os
import uproot
import numpy as np

def toVector(vtype, list):
    from ROOT import vector
    v = vector(vtype)()
    for o in list:
        v.push_back(o)
    return v

def get_nEvents(fileglob, key = 'nEventsGenWeighted'):
    '''
    Finds total events pre cut.
    Input:
        fileglob: list of files
        key: optional, key for one bin hist containing counts.
        
    Output:
        scalar of total events added.
    '''
    file_path = fileglob
    total = 0
    for file in fileglob:
        up_f = uproot.open(file)
        # Try two different methods to get value
        try:
            total += up_f[key].numpy()[0][0]
        except:
            total += up_f[key].values()[0]
    return total

def make_view(mass_cut = [-np.inf,np.inf], HTLT = np.inf, RelMET = np.inf, SBM = 0, MET_filter = 1, region=0):
    view =  {
        'lt':{
            'DiLepMass': mass_cut[1],
            'HTLT_nominal': HTLT,
            'RelMET_nominal': RelMET
        },
        'gt':{
            'DiLepMass': mass_cut[0],
            'SBM_nominal': SBM
        },
        'eq':{
            'Flag_METFilters': 1,
        }
    }
    if region: view['eq'][region] = 1
    return view

def ratio_plot_template(**kwargs):
    import matplotlib.pyplot as plt
    fig = plt.figure(constrained_layout=True, **kwargs)
    gs = fig.add_gridspec(10,10)
    top = fig.add_subplot(gs[0:8,:])
    bottom = fig.add_subplot(gs[8:,:])
    return fig, top, bottom

def nratio_plot_template(nPlots=[2,2],rps = 2,**kwargs):
    '''
This produces a grid of ratio plots to work with of arbitrary size.
Caveat: gridspec uses the x corrdinate to refer to the height, while figsize uses x for width.
Here, x_bins then refers to height in the grid spec.
The grid kwarg is flipped so the first corrdiate is width.
'''
    import matplotlib.pyplot as plt
    import itertools
    fig = plt.figure(constrained_layout=True, **kwargs)
    
    size = np.flip(nPlots)
    one_plot_size = np.array([10,10])
    gs_size = np.array(one_plot_size)*size
    gs = fig.add_gridspec(*gs_size)
    #calculate boundaries for each cell
    x_bins = np.linspace(0,gs_size[0],size[0]+1)
    x_bins = [[x_bins[i],x_bins[i+1]] for i in range(size[0])]
    y_bins = np.linspace(0,gs_size[1],size[1]+1)
    y_bins = [[y_bins[i],y_bins[i+1]] for i in range(size[1])]
    #conver to int
    x_bins = np.array(x_bins,dtype=int)
    y_bins = np.array(y_bins,dtype=int)
    #ratio plot size
    assert rps < one_plot_size[0], "rps size must be smaller than plot size"
    axes = []
    print(axes)
    for j,y in enumerate(y_bins):
        axes.append([])
        for i,x in enumerate(x_bins):
            tt,tb, bb = x[0], x[1]-rps, x[1]
            l,r = y
            top = fig.add_subplot(gs[tt:tb,l:r])
            bottom = fig.add_subplot(gs[tb:bb,l:r])
            axes[j].append([top,bottom])
    return fig, axes

def color_map(n_things):
    import matplotlib.colors
    import matplotlib.pyplot as plt
    cmap = plt.cm.rainbow
    colors = np.linspace(0,1,n_things)
    return [cmap(c) for c in colors]

def log_norm(x, norm, sigma, theta, mean):
    return norm/((x-theta)*sigma*2*3.14159)*np.exp(-(np.log((x-theta)/mean))**2/(2*sigma**2))

def log_norm_unp(x, norm, sigma, theta, mean):
    import uncertainties.unumpy as unp
    return norm/((x-theta)*sigma*2*3.14159)*unp.exp(-(unp.log((x-theta)/mean))**2/(2*sigma**2))
import uncertainties

def unzip_unp(unp_arr):
    import uncertainties
    return np.array(list(zip(*[(x.nominal_value, x.std_dev) for x in unp_arr])))

def time_func(func):
    from functools import wraps
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = perf_counter()
        res = func(*args, **kwargs)
        print("{:.2e}".format(perf_counter()-start))
        return res
    return wrapper


def hist2unc(hist):
    import uncertainties
    from uncertainties import ufloat
    val = hist.values()
    var = hist.variances()**.5
    return np.array([ufloat(vl,vr) for vl, vr in zip(val,var)])

def hist2array(hist):
    return hist.values()

def sum_in_quad(np_arr, **kwargs):
    return np.sum(np_arr**2, **kwargs)**.5

vhist2unc = np.vectorize(hist2unc)
vhist2array = np.vectorize(hist2array)
vunc2nom = np.vectorize(lambda x: x.nominal_value)
vunc2std = np.vectorize(lambda x: x.std_dev)

def parabola(x,offset,m1,m2,m3):
    x = x - offset
    return m1+m2*x+m3*x**2



def power_func(x, c, p):
    return c*x**p

def heaviside(x, cutoff, scale):
    return np.heaviside(x-cutoff, 0)*scale

def significance(sig,bck):
    return sig/(sig+bck+1e-12)**.5

def apply_multiple_filters(df, filter_list):
    filter_prod = np.array(filter_list).prod(axis=0)
    return df[filter_prod==1]


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
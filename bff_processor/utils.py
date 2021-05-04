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

def get_nEvents(fileglob, key = 'sumPosMinusNegGenWeight'):
    file_path = fileglob.replace('*.root','')
    files = ["{}/{}".format(file_path,f) for f in os.listdir(file_path) if '.root' in f]
    total = 0
    for file in files:
        up_f = uproot.open(file)
        total += up_f[key].numpy()[0][0]
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
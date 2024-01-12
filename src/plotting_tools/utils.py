import matplotlib.pyplot as plt
import numpy as np

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


def pair_wise_list_operation(arr, operation = lambda x,y: (x+y)/2): 
    return [operation(arr[i], arr[i+1]) for i in range(len(arr)-1)]

def calc_bin_widths(arr):
    return pair_wise_list_operation(arr, operation=lambda x, y: -x+y)

def calc_bin_centers(arr):
    return pair_wise_list_operation(arr)

def rebin_np(oldbins, newbins, data, drawplot=False):
    values = np.histogram(oldbins, newbins, weights = data)
    return values[0]

# this is much better than the np method
def rebin_integrate(x,y, xnew_binedges, keep_norm=False, usequad=False, straight_interpolate = False, sample_density=100):
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    f = interp1d(x,y, fill_value=0, bounds_error=False)
     
    centers = np.array(Bins(bin_edges).calc_bin_centers())
    widths = np.array(Bins(bin_edges).calc_bin_widths())
    if usequad:
        ynew_int = []
        for i, center in enumerate(centers):
            integral = quad(f, xnew_binedges[i], xnew_binedges[i+1])[0]
            ynew_int.append(integral)
        ynew_int = np.array(ynew_int)
    elif straight_interpolate:
        ynew_int = f(centers)*widths
    else:
        ynew_int = []
        for i, (center, widht) in enumerate(zip(centers, widths)):
            sample_points = np.linspace(xnew_binedges[i], xnew_binedges[i+1],int(sample_density*widht))
            integral = f(sample_points).sum()/sample_density
            ynew_int.append(integral)
        ynew_int = np.array(ynew_int)
        
    
    
    if keep_norm:  
        y_total = y[(x>=np.min(bin_edges)) & (x<=np.max(bin_edges))].sum()
        return centers, ynew_int/ynew_int.sum()*y_total
    return centers, ynew_int
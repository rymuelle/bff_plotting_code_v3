import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import boost_histogram as bh
from uncertainties import ufloat

def combine_sys(hist_dict, sys_cat_key, nominal):
        sys_list = []
        for sys_key in hist_dict[sys_cat_key]:
            up, down = hist_dict[sys_cat_key][sys_key]
            up,down = up[0]-nominal, down[0]-nominal
            #sort the bins
            up_down = np.array([up,down]).T
            up_down.sort(axis=1)
            up_down = up_down.T
            sys_list.append(up_down)
        return np.sum(sys_list,axis=0)
    
def double_list(l):
    doubled = np.array([l,l]).T
    return doubled.reshape(-1)

def filter_bins(x_bins,min_val):
    b = (x_bins>=min_val)
    bin_filter = np.array([float(b[i])+float(b[2]) for i in range(len(b)-1)])
    bin_filter = bin_filter==2.0
    return bin_filter, b

def make_plot(hist_dict):
    x_nom,x_stat,x_bins,x_center = hist_dict['nom']
    x_center = np.array(x_center)
    x_nom = np.array(x_nom) 
    sys = np.full([2,len(x_nom)],0.)
    if 'jv' in hist_dict: 
        sys += combine_sys(hist_dict, 'jv', x_nom)
    if 'weights' in hist_dict: 
        weight_sys = combine_sys(hist_dict, 'weights', x_nom)
        sys += weight_sys
    return x_nom,x_stat,sys[0],sys[1], x_bins,x_center

td_dict = {
    'data':{'label': 'data', 'type': 'data', 'kwargs':{'color': 'black'}},
    'BFFZprimeToMuMu_M_175':{'label': '175 $\delta_{bs}=0.04$', 'type': 'signal', 'kwargs':{'color': '#3380ff'}},
    'BFFZprimeToMuMu_M_200':{'label': '200 $\delta_{bs}=0.04$', 'type': 'signal', 'kwargs':{'color': '#6633ff'}},
    'BFFZprimeToMuMu_M_200_dbs0p5':{'label': '200 $\delta_{bs}=0.5$', 'type': 'signal', 'kwargs':{'color': '#3900e6','linestyle':'--'}},
    'BFFZprimeToMuMu_M_200_dbs1p0':{'label': '200 $\delta_{bs}=1.0$', 'type': 'signal', 'kwargs':{'color': '#260099','linestyle':':'}},
    'BFFZprimeToMuMu_M_350':{'label': '350 $\delta_{bs}=0.04$', 'type': 'signal', 'kwargs':{'color': '#ff33ff'}},
    'BFFZprimeToMuMu_M_350_dbs0p5':{'label': '350 $\delta_{bs}=0.5$', 'type': 'signal', 'kwargs':{'color': '#e600e6','linestyle':'--'}},
    'BFFZprimeToMuMu_M_350_dbs1p0':{'label': '350 $\delta_{bs}=1.0$', 'type': 'signal', 'kwargs':{'color': '#990099','linestyle':':'}},
    'BFFZprimeToMuMu_M_500':{'label': '500 $\delta_{bs}=0.04$', 'type': 'signal', 'kwargs':{'color': '#ff3366'}},
    'BFFZprimeToMuMu_M_500_dbs0p5':{'label': '500 $\delta_{bs}=0.5$', 'type': 'signal', 'kwargs':{'color': '#e60039','linestyle':'--'}},
    'BFFZprimeToMuMu_M_500_dbs1p0':{'label': '500 $\delta_{bs}=1.0$', 'type': 'signal', 'kwargs':{'color': '#990026','linestyle':':'}},

    'WW/WZ/ZZ':{'label': 'WW/WZ/ZZ', 'type': 'background', 'kwargs':{'color': '#ccffcc'}},
    'ST':{'label': 'ST', 'type': 'background', 'kwargs':{'color': '#ffccff'}},
    'TT':{'label': 'TT', 'type': 'background', 'kwargs':{'color': '#ffff99'}},
    'DY':{'label': 'DY', 'type': 'background', 'kwargs':{'color': '#99ffff'}},
    }
td_dict = [{'hname':key, **td_dict[key]} for key in td_dict]
for row in td_dict:
    row['color'] = row['kwargs']['color']
    row.pop('kwargs')
meta_df = pd.DataFrame(td_dict) 

def center(l_list): return [(l_list[i]+l_list[i+1])/2. for i in range(len(l_list)-1)]

def make_stack(df, reg_dict_temp, region, sum_hist=True,**kwargs):
    hnames = df.hname
    colors = df.color
    labels = df.label
    hists = [make_plot(reg_dict_temp[region][hname], **kwargs) for hname in hnames if hname in reg_dict_temp[region]]
    if len(hists)==0: return [], [], [], []
    #x_nom,x_stat,sys, sys_p_stats , x_bins,x_center, doubled_x, doubled_sys_p_stats
    hists = np.array(hists, dtype=object)
    nom = hists[:,0]
    stats = hists[:,1]
    sysup = hists[:,2]
    sysdown = hists[:,3]
    x_bin = hists[0,4]
    if sum_hist:
        nom = np.sum(nom,axis=0)
        stats = np.sum(stats**2,axis=0)**.5
        sysup = np.sum(sysup,axis=0)
        sysdown = np.sum(sysdown,axis=0)
        return nom, stats, np.array([sysup,sysdown]), x_bin
    else:
        return nom, stats, np.array([sysup,sysdown]), center(x_bin), colors, labels

def plot_w_error(ax, nom, stat, sys, bins, color='black', alpha=.25, **kwargs):
    comb_sys_stat = [sys[0]-stat+nom, sys[1]+stat+nom]
    double_edge = double_list(bins)[1:-1]
    double_comb_sys_stat = [double_list(comb_sys_stat[0]), double_list(comb_sys_stat[1])]
    dnom = double_list(nom)
    bin_center = [(bins[i]+bins[i+1])/2. for i in range(len(nom))]
    ax.plot(double_edge, dnom, color=color, **kwargs)
    ax.fill_between(double_edge, *double_comb_sys_stat, color=color, alpha=alpha)

    
def make_plot_boost(hist_dict):
    x_nom,x_stat,x_bins,x_center = hist_dict['nom']
    x_center = np.array(x_center)
    x_nom = np.array(x_nom) 
    sys = np.full([2,len(x_nom)],0.)
    if 'jv' in hist_dict: 
        sys += combine_sys(hist_dict, 'jv', x_nom)
    if 'weights' in hist_dict: 
        weight_sys = combine_sys(hist_dict, 'weights', x_nom)
        sys += weight_sys
    return x_nom,x_stat,sys[0],sys[1], x_bins,x_center

    
def make_stack_boost(df, reg_dict_temp, region, sum_hist=True,**kwargs):
    hnames = df.hname
    colors = df.color
    labels = df.label
    hists = [make_plot(reg_dict_temp[region][hname], **kwargs) for hname in hnames if hname in reg_dict_temp[region]]
    if len(hists)==0: return [], [], [], []
    #x_nom,x_stat,sys, sys_p_stats , x_bins,x_center, doubled_x, doubled_sys_p_stats
    hists = np.array(hists, dtype=object)
    nom = hists[:,0]
    stats = hists[:,1]
    sysup = hists[:,2]
    sysdown = hists[:,3]
    x_bin = hists[0,4]
    if sum_hist:
        nom = np.sum(nom,axis=0)
        stats = np.sum(stats**2,axis=0)**.5
        sysup = np.sum(sysup,axis=0)
        sysdown = np.sum(sysdown,axis=0)
        return nom, stats, np.array([sysup,sysdown]), x_bin
    else:
        return nom, stats, np.array([sysup,sysdown]), center(x_bin), colors, labels
    
def make_plot_boost(hist_dict):
    nom = hist_dict['nom']
    sysup = []
    sysdown = []
    if 'jv' in hist_dict:
        for key in  hist_dict['jv']:
            hists = hist_dict['jv'][key]
            hists = sorted(hists,key= lambda x: -x.sum().value)
            sysup.append(-1*nom+hists[0])
            sysdown.append(-1*nom+hists[1])
    if 'weights' in hist_dict:
        for key in  hist_dict['weights']:
            hists = hist_dict['weights'][key]
            hists = sorted(hists,key= lambda x: -x.sum().value)
            sysup.append(-1*nom+hists[0])
            sysdown.append(-1*nom+hists[1])
    return nom, sum(sysup), sum(sysdown)
    
def make_stack_boost(df, reg_dict_temp, region, sum_hist=True):
    hnames = df.hname
    colors = df.color
    labels = df.label
    hists = [make_plot_boost(reg_dict_temp[region][hname]) for hname in hnames if hname in reg_dict_temp[region]]
    nom_l=[]
    sysup_l = []
    sysdown_l =[]
    for nom,sysup, sysdown in hists:
        nom_l.append(nom)
        sysup_l.append(sysup)
        sysdown_l.append(sysdown)
    if sum_hist:
        nom = sum(nom_l)
        sysup = sum(sysup_l)
        sysdown = sum(sysdown_l)
        return nom, [sysup,sysdown]
    else:
        return nom_l, [sysup_l,sysdown_l], colors, labels
    
def plot_w_error_boost(ax, nom, sys, color='black', alpha=.25, **kwargs):
    nom_values = nom.values()
    nom_var = nom.variances()
    if sys[0]: 
        sysup = np.array([sys[0].values() ** 2 + nom_var]) ** .5 + nom_values
        sysdown = -np.array([sys[1].values() ** 2 + nom_var]) ** .5 + nom_values 
    else:
        sysup = np.array([nom_var]) ** .5 + nom_values
        sysdown = -np.array([nom_var]) ** .5 + nom_values 
    double_edge = double_list(nom.axes[0].edges)[1:-1]
    double_comb_sys_stat = [double_list(sysdown), double_list(sysup)]
    bin_center = nom.axes[0].centers
    ax.plot(double_edge, double_list(nom_values), color=color, **kwargs)
    ax.fill_between(double_edge, *double_comb_sys_stat, color=color, alpha=alpha)
    return double_list(nom_values), double_comb_sys_stat


def boost_from_unc(unc_array, edges, flip=0):
    step = 1
    if flip: step = -1
    b_nom = np.array(([x.nominal_value for x in unc_array]))[::step]
    b_std = np.array(([x.std_dev for x in unc_array]))[::step]
    b_hist = bh.Histogram(bh.axis.Variable(edges[::step]))
    return b_nom, b_std
    #print(b_hist)
    b_hist.values()[:] = b_nom
    #print(b_hist)
    b_hist.variances()[:] = b_std**2 
    # print(b_hist)
    return b_hist

def boost2unc(bhist):
    from uncertainties import ufloat
    return np.array([ufloat(x,y) for x, y in np.array(bhist)])

def boost_ratio_hist(b1,b2, **kwargs):
    edges = b1.axes[0].edges
    b1_unc = boost2unc(b1)
    b2_unc = boost2unc(b2)
    print(b2_unc[b2_unc==0])
    b1_rat = b1_unc/b2_unc
    return boost_from_unc(b1_rat, edges, **kwargs)

def produce_bff_hists(df, name, columns, weight='Weight'):
    '''This function helps produce plots for bff cut selection'''
    from itertools import combinations
    hist_1d_dict = {}
    for column, bin_meta in columns:
        hist_1d_dict[column] = bh.Histogram(
            bh.axis.Regular(*bin_meta, metadata="{} {}".format(column,name)),
            storage=bh.storage.Weight()
            )
        hist_1d_dict[column].fill(df[column], weight=df[weight])
    
    hist_2d_dict= {}
    for (c1,bm1), (c2,bm2) in combinations(columns, 2):
        hist_2d_dict["{}_{}".format(c1,c2)] = bh.Histogram(
            bh.axis.Regular(*bm1, metadata="{} {}".format(c1,name)),
            bh.axis.Regular(*bm2, metadata="{} {}".format(c2,name)),
            storage=bh.storage.Weight()
        )
        hist_2d_dict["{}_{}".format(c1,c2)].fill(df[c1],df[c2], weight=df[weight])
    
    return hist_1d_dict, hist_2d_dict

def boost_plot(ax, bh, **kwargs):
    val, var = bh.values(), bh.variances()
    center = bh.axes[0].centers
    ax.errorbar(center, val, yerr=var**.5, **kwargs)
    
def boost_plot2d(ax, h, lock_aspect=0, log=0, min_val=.1, **kwargs):
    w, x, y = h.to_numpy()
    # Draw the count matrix
    if not log:
        ax.pcolormesh(x, y, w.T)
    else:
        import matplotlib.colors as colors
        vmin = max(w.T.min(),min_val)
        vmax = w.T.max()
        ax.pcolor(x, y, w.T,
                   norm=colors.LogNorm(vmin=vmin, vmax=vmax),
                   cmap='PuBu_r', shading='auto')
    ax.set_xlabel(h.axes[0].metadata)
    ax.set_ylabel(h.axes[1].metadata)
    if lock_aspect: ax.set_aspect("equal")

def unc_plot(ax, centers, unc, fill_between=False, **kwargs):
    from bff_processor.utils import vunc2nom, vunc2std
    val = vunc2nom(unc)
    std = vunc2std(unc)
    if not fill_between:
        ax.errorbar(centers, val, std, **kwargs)
    if fill_between:
        ax.plot(centers, val, **kwargs)
        fb_kwargs = {**kwargs, 'label':None}
        ax.fill_between(centers, val+std, val-std, alpha=.25, **fb_kwargs)


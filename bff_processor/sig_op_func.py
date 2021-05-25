from bff_processor.utils import significance, hist2unc
from bff_processor.plotting_utils import produce_bff_hists, boost_plot, unc_plot
from bff_processor.bff_meta import preselection, band_cut
import numpy as np
from scipy.stats import norm

from scipy.optimize import minimize
from bff_processor.utils import vunc2nom, vunc2std
import uncertainties
import matplotlib.pyplot as plt
import re 

columns = [
    ['DiLepMass', [139,105,800]],
    ['TMB_nom',       [80,0, 800]],
    ['HTLT_nom',      [100,-500,500]],
    ['RelMET_nom',    [100, 0,1]]
]

#optimize RelMET cut:
def column_1d_sig(sig_df, bck_df, mass_band=1, filter_func=lambda x: x):
    #select for region
    _sig_df = filter_func(sig_df)
    _bck_df =  filter_func(bck_df)
    mean, std = sig_df['DiLepMass'].mean(), sig_df['DiLepMass'].std()
    widht = std*mass_band
    #mass band cut
    _sig_df = band_cut('DiLepMass',mean-widht, mean+widht)(_sig_df)
    _bck_df = band_cut('DiLepMass',mean-widht, mean+widht)(_bck_df)
    bck_1d_hist, bck_2d_hist = produce_bff_hists(_bck_df, "", columns, weight='Weight')
    sig_1d_hist, sig_2d_hist = produce_bff_hists(_sig_df, "", columns, weight='Weight')
    return bck_1d_hist, sig_1d_hist, bck_2d_hist, sig_2d_hist

def calc_sig_cut(s,b, direction=1):
    from itertools import accumulate
    s,b = np.array(list(accumulate(s[::direction])))[::direction], np.array(list(accumulate(b[::direction])))[::direction]
    return significance(s,b)

def map_signf(value, signf, centers):
    diff = np.abs(centers - value)
    delta  = (centers[1] - centers[0])
    weight = norm.pdf(diff, loc=0, scale=delta*1)
    wval = np.dot(signf,weight)/np.sum(weight+1e-12)
    return wval.nominal_value

def fit_sig(masses, sigs, centers, fit_function, *popt):
    signfs = []
    for m, sig in zip(masses,sigs):
        cut_val = fit_function(m,*popt)
        signf = map_signf(cut_val, sig, centers)
        signfs.append(signf+1e-12)
    return np.dot(-np.array(signfs), np.power(masses, 0))

def plot_opt_sig(column, bff_dict, background_df, bff_samples, func2fit, minimize_func=True,  filter_func=lambda x: x, postfix="",direction=1, p0=[0,.5], title_size=20):
    fig,ax = plt.subplots(2,len(bff_samples), figsize=[38,15])
    fit_x = []
    fit_y = []
    fit_y_unc = []
    sigs = []
    masses = []
    
    for j, m in enumerate(bff_samples):
        mass = int(re.findall('([0-9]+)',m)[0])
        if m not in bff_dict:
            masses.append(mass)
            sigs.append([])
            continue
        bck_1d_hist, sig_1d_hist, bck_2d_hist, sig_2d_hist = column_1d_sig(bff_dict[m], background_df, mass_band=2,filter_func=filter_func)
        
        top_ax = ax[0,j]
        top_ax.set_xlabel(column.replace('_',' '))
        boost_plot(top_ax, bck_1d_hist[column] ,label='bck')
        boost_plot(top_ax, sig_1d_hist[column] ,label='sig')
        title = r'{} GeV $\delta_{{bs}}$ {}'.format(*m.split(' '))
        top_ax.set_title(title, fontsize=title_size)
        top_ax.legend()
        
        s,b = hist2unc(sig_1d_hist[column]), hist2unc(bck_1d_hist[column])
        signf_val = calc_sig_cut(s,b, direction=direction)
        
        centers = sig_1d_hist[column].axes[0].centers
        unc_plot(ax[1][j], signf_val,centers, zorder=1)
        peak_index = np.argmax(signf_val)
        peak_center = centers[peak_index]
        peak_height = signf_val[peak_index]
        
        sigs.append(signf_val)
        masses.append(mass)
            
        #uncertainty of center
        upper_limit = vunc2nom(signf_val)+vunc2std(signf_val)
        in_range = centers[upper_limit > peak_height]
        sigma_range = np.asarray((in_range[0], in_range[-1]))
        unc = np.max(np.abs(sigma_range-peak_center))
        fit_x.append(mass)
        fit_y.append(peak_center)
        fit_y_unc.append(unc)
        
        ax[1][j].set_title("max: {:.2f}+/-{:.2f}".format(peak_center, unc), fontsize=title_size)
    
    ax[0][0].set_ylabel('Counts')
    ax[1][0].set_ylabel(r'Significance $\frac{s}{\sqrt{s+b}}$')
    
    if minimize_func:
        fit_func = lambda popt: fit_sig(masses, sigs, centers, func2fit, *popt)
        min_val = minimize(fit_func, p0)
        print(min_val)
        popt = min_val.x
    else:
        popt = p0
    masses = np.array(masses)
    #cut_points = linear(masses, *popt)
    cut_points = func2fit(masses, *popt)
    
    for j, (m, c, s) in enumerate(zip(masses, cut_points, sigs)):

        if len(s)<10: continue
        sig_at_cut = map_signf(c, s, centers)
        print(sig_at_cut, c)
        #ax[1][j].plot([c,c], [0,sig_at_cut], zorder=1)
        ax[1][j].scatter([c], [sig_at_cut], zorder=1,s=150, c='red', marker=(5, 1), label='cut loc.')
        ax[1][j].legend(loc=4)
    return popt, fig
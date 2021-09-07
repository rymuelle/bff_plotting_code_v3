from bff_signal_model.bff_signal_model import reset_params, bff_signal_model, mu
import zfit
from zfit import z
import numpy as np
from bff_processor.utils import nratio_plot_template, hist2unc, vunc2nom, chiSquared, color_map, quad, linear, constant
import mplhep as hep

def get_unique_masses(df, x_range=[110,800]):
    bff_data = df[df.name.str.contains("BFF")]
    df  = bff_data[(bff_data.DiLepMass > x_range[0]) & (bff_data.DiLepMass < x_range[1]) & (bff_data.dbs==0.04)]
    masses = df.mass.unique()
    masses = sorted(masses)
    return masses

def sigma_from_mass(mass):
    return quad(mass, *[-8.46335525e-01,  1.66590567e-02,  1.51511530e-05])

def chi2(xm, xp, error, ndf=0):
    y = (xm - xp) ** 2 / error **2
    return np.sum(y)/(len(y) - ndf)

def make_plot_dict(df, obs, masses, compute_hesse=True, regions = ['SR1', 'SR2']):
    tail_sys = 1
    tail_sys_start = 2
    constant_sys = .05
    x_range = obs.limit1d
    param_list = []
    plot_dict = {}
    plot_dict_centered = {}
    df = df[df.dbs==0.04]
    for reg in regions:
        plot_dict[reg] = {}
        plot_dict_centered[reg] = {}
        for mass in masses:
            print(mass)
            if reg == "Both":
                tdf = df[(df.mass==mass)]
            else:
                tdf = df[(df.mass==mass) & (df[reg+'_nom']==1)]
            data, weights = tdf.DiLepMass.to_numpy(),tdf.Weight.to_numpy()
            mean, std = reset_params(data)
    
            doublecb = bff_signal_model(obs=obs, mu=mu)
    
            #set up fit
            data_zfit = zfit.Data.from_numpy(obs=obs, array=data, weights=weights)
            nll = zfit.loss.UnbinnedNLL(model=doublecb, data=data_zfit)
            minimizer = zfit.minimize.Minuit()
            result = minimizer.minimize(nll)
            if compute_hesse:
                x = result.hesse()
            else:
                x = {}
            param_dict = {'reg': reg, 'mass': mass, **{p.name:p.value().numpy() for p in result.params}, **{p.name+"_error":x[p]['error'] for p in x}}
            param_list.append(param_dict)
    
            #make the plot
            bins = np.linspace(*x_range, int((x_range[1]-x_range[0])/5 + 1))
            print(len(bins))
            y, y_unc = doublecb.tail_sys(bins, np.sum(weights), width=tail_sys_start, supersample=100,
                                           tail_percent=tail_sys, constant_percent=constant_sys, stat_unc=0,)
            hist, _ = np.histogram(data, weights=weights, bins=bins)
            hist_var, _ = np.histogram(data, weights=weights **2 , bins=bins)
            hist_std = hist_var ** .5
    
            plot_dict[reg][mass] = {"fit": y, "fit_unc": y_unc, "hist": hist, "hist_std": hist_std, "bins":bins}
    
            #make tail sys plots
            width = 10
            bins_centered_peak = np.linspace(mean-width*std, mean+width*std, 2*5*width+1)
            y, y_unc = doublecb.tail_sys(bins_centered_peak, np.sum(weights), width=tail_sys_start, supersample=100, area=width*std*2,
                                        tail_percent=tail_sys, constant_percent=constant_sys, stat_unc=0,)
            y2 = doublecb.fill_bins(bins_centered_peak,  np.sum(weights), supersample=100, area=width*std*2)
            hist, _ = np.histogram(data, weights=weights, bins=bins_centered_peak)
            hist_var, _ = np.histogram(data, weights=weights **2 , bins=bins_centered_peak)
            hist_std = hist_var ** .5
            bins_centered_peak = (bins_centered_peak-mean)/std
    
            plot_dict_centered[reg][mass] = {"fit": y, "fit_unc": y_unc, "hist": hist, "hist_std": hist_std, "bins":bins_centered_peak}
    return plot_dict, plot_dict_centered, param_list

def compute_bin_centers(bins):
    return [(bins[i]+bins[i+1])/2 for i in range(len(bins) - 1)]

def prepare_plots(mass_dict, mass, systematics=0, width = 100):
    fit_plot = mass_dict['fit']
    fit_unc_plot = mass_dict['fit_unc']
    hist = mass_dict['hist']
    hist_std = mass_dict['hist_std']
    bins = mass_dict['bins']
    bin_centers = compute_bin_centers(bins)
    # remove low content bins:
    total = np.sum(fit_plot)
    max_bin, min_bin = mass+sigma_from_mass(mass)*width,mass-sigma_from_mass(mass)*width
    filter_array =np.logical_and(bin_centers < max_bin, bin_centers > min_bin)
    hist = hist[filter_array]
    hist_std = ((hist_std[filter_array] ** 2) + (hist*systematics ** 2)) ** .5
    bin_centers_temp =np.array(bin_centers)[filter_array]
    fit_plot = fit_plot[filter_array]
    fit_unc_plot = fit_unc_plot[filter_array]
    return bins, bin_centers, bin_centers_temp, hist, hist_std, fit_plot, fit_unc_plot
    
def make_stack_plot(plot_dict, masses, pdf, lumi, era, compute_hesse, systematics=0, width=100, legend=True,
bottom_limit=[-3,3],postfix="",
                   yscale='log'):
    residual_dict = {}
    colors = color_map(len(masses))
    for reg, reg_dict in plot_dict.items():
        residual_dict[reg] = {}
        fig, ax = nratio_plot_template(nPlots=[1,1])
        (top, bottom) = ax[0][0]
        reg_dict = {k:i for k,i in reg_dict.items() if k in masses}
        for color, (mass, mass_dict) in zip(colors,reg_dict.items()):
            #make a different color for histogram
            hist_color = np.power(color, 1)*.75
            #set alpha
            hist_color[-1] = 1
            bins, bin_centers, bin_centers_temp, hist, hist_std, fit_plot, fit_unc_plot = prepare_plots(mass_dict, mass, systematics=systematics, width=width)

            top.errorbar(bin_centers_temp, hist, yerr=hist_std, label='{} GeV'.format( mass), color=hist_color, linestyle="None", marker='o')
            top.errorbar(bin_centers_temp, fit_plot, yerr=fit_unc_plot, color=color)
            bottom.errorbar(bin_centers_temp, (fit_plot-hist)/fit_unc_plot, yerr=fit_unc_plot/fit_unc_plot, color=color)
            residual_dict[reg][mass] = {"residual": (fit_plot-hist)/fit_plot, "std":hist_std/fit_plot, "bin_centers": bin_centers}
            
            total_unc = (fit_unc_plot**2+hist_std**2)**.5
            pdf.loc[(pdf.mass==mass) & (pdf.reg==reg), 'chi2'] = chi2(fit_plot, hist, total_unc , ndf=1)
    
        if legend: top.legend(title="{}, {} sys, {} sigma".format(reg, systematics, width))
        top.set_yscale(yscale)
        top.set_ylim(bottom=.1e-2,top=1e4)
        bottom.set_ylim(*bottom_limit)
        bottom.plot(top.get_xlim(), np.full(len(top.get_xlim()), -1), color='black', linestyle=':')
        bottom.plot(top.get_xlim(), np.full(len(top.get_xlim()), 0), color='black', linestyle='--')
        bottom.plot(top.get_xlim(), np.full(len(top.get_xlim()), 1), color='black', linestyle=':')
        bottom.set_xlabel('DiLepMass [GeV]')
        top.set_ylabel('Count per 5 GeV')
        hep.cms.label(loc=0,ax=top,lumi=lumi,year=era, data=False)
        if compute_hesse:
            fig.savefig('fits/bff/bff_mass_only_model_{}_{}{}.png'.format(reg, era, postfix))
        else:
            fig.savefig('fits/bff/bff_mass_only_model_{}_{}_nohess{}.png'.format(era, reg, postfix))
    return residual_dict

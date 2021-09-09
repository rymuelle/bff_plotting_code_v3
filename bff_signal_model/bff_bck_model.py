import zfit
from zfit import z
from zfit.models.physics import DoubleCB
from scipy.optimize import curve_fit
from zfit.models.physics import double_crystalball_func
from bff_signal_model.functions import *
import numpy as np
from bff_plotting_tools.make_hists import SysHist
from math import pi 

try:
    sigma    = zfit.Parameter("sigma"   , .8, 0, 1.5)
    theta = zfit.Parameter("theta", 80,  0, 100)
    mean = zfit.Parameter("mean", 91,  0, 800)
except:
    print("already defined")

class bff_bck_model(zfit.pdf.ZPDF):
    """1-dimensional PDF implementing the exp(alpha * x) shape."""
    _PARAMS = ['sigma', 'theta', 'mean']  # specify which parameters to take

    def lognorm(self, data, module, sigma=.5, theta=100, mean=90):
        return 1./((data-theta)*sigma*2*pi)*module.exp(-(module.log((data-theta)/mean))**2/(2*sigma**2))
        
    def _unnormalized_pdf(self, x):  # implement function
        data = z.unstack_x(x)
        sigma = self.params['sigma']
        theta = self.params['theta']
        mean = self.params['mean']
        return self.lognorm(data, z.numpy, sigma=sigma, theta=theta, mean=mean)

    def scaled_pdf(self, x, sumW, area=0):
        if area == 0: area = self.norm_range.area()
        y = self.pdf(x)
        n_bins = len(x)
        plot_scaling = sumW / n_bins * area
        return (y * plot_scaling).numpy()
    
    def par2ufloat(self, name): 
        assert self.result, "Need to add result object"
        from uncertainties import ufloat
        return  ufloat(self.result.params[name]['value'], self.result.params[name]['minuit_hesse']['error'])

    def super_sample(self, bins, supersample):
        bin_edges, nBins = bins.bin_edges, bins.calc_nBins()
        return [np.linspace(bin_edges[i], bin_edges[i+1], supersample) for i in range(nBins)]

    def fill_bins(self, bins, *args, **kwargs):
        y_unc = self.super_sample_unc(bins, *args, **kwargs)
        y  = np.array(list(map(lambda x: x.nominal_value,y_unc)))
        unc  = np.array(list(map(lambda x: x.std_dev,y_unc)))
        return SysHist(y, -y*0, y*0, unc, bins.bin_edges)
    
    def super_sample_unc(self, bins, norm, supersample=100, area=0):
        bin_widths= bins.calc_bin_widths()
        supersampled_bins = self.super_sample(bins, supersample)
        percent_weight = bin_widths/np.sum(bin_widths)
        supersampled_bin_values = list(map(lambda x: self.scaled_pdf_uncertainty(x, 1, area=area), supersampled_bins))
        y_unc = np.sum(supersampled_bin_values, axis=1)*percent_weight
        y_unc = y_unc * norm / np.sum(y_unc)
        return y_unc

    def scaled_pdf_uncertainty(self, data, norm, area=0):
        data = np.array(data)
        import uncertainties.unumpy as unp
        name = []
        values = []
        for x in self.result.params:
            name.append(x.name)
            values.append(x.value().numpy())
        import uncertainties
        popt = uncertainties.correlated_values(values, self.result.covariance())
        y = self.lognorm(data, unp, **{n:p for n, p in zip(name,popt)})
        n_bins = len(x) 
        if area == 0: area = self.norm_range.area()
        #plot_scaling = float(sumW2 / n_bins * area)
        return y / norm


def logistics_func(x, offset, widht):
    return 1/(1+10.**(-(x-offset)/widht))
    
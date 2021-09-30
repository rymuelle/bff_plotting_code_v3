import zfit
from zfit import z
from zfit.models.physics import DoubleCB
from scipy.optimize import curve_fit
from zfit.models.physics import double_crystalball_func
from bff_signal_model.functions import *
import numpy as np
from bff_plotting_tools.make_hists import SysHist

try:
    mu = zfit.Parameter("mu", 300,  100, 600)
    sigma = zfit.Parameter("sigma", 20,  0, 100)
    alphal = zfit.Parameter("alphal", 2,  0, 100)
    nl = zfit.Parameter("nl", 3,  -10, 10)
    alphar = zfit.Parameter("alphar", 1.5,  0, 100)
    nr = zfit.Parameter("nr", 6,  -100, 100)
    Nsig = zfit.Parameter("Nsig", 1., -20., 1e8)
except:
    print("already defined")

class bff_signal_model(zfit.pdf.ZPDF):
    """1-dimensional PDF implementing the double cb signal model"""
    _PARAMS = ['mu']  # specify which parameters to take
    def mu(self):
        return linear(self.params['mu'], *[0.57839198, 0.99561578])
    def nr(self):
        return linear(self.params['mu'], *[1.99831802, 0.01178409])
    def sigma(self):
        return quad(self.params['mu'], *[-8.46335525e-01,  1.66590567e-02,  1.51511530e-05])
    def nl(self):
        return linear(self.params['mu'], *[2.42, 0])
    def alphal(self):
        return linear(self.params['mu'], *[1.38, 0])
    def alphar(self):
        return linear(self.params['mu'], *[1.8, 0])

    def _unnormalized_pdf(self, x):  # implement function
        data = z.unstack_x(x)
        mu = self.params['mu']
        nr = self.nr()
        sigma = self.sigma()
        nl = self.nl()
        alphal = self.alphal()
        alphar = self.alphar()
        return double_crystalball_func(data,mu=mu, sigma=sigma, 
                                       alphal=alphal, alphar=alphar, 
                                       nl=nl, nr=nr)
    def scaled_pdf(self, x, sumW, area=0):
        if area == 0: area = self.norm_range.area()
        y = self.pdf(x)
        n_bins = len(x)
        plot_scaling = sumW / n_bins * area
        return (y * plot_scaling).numpy()

    def super_sample(self, bins, supersample):
        bin_edges, nBins = bins.bin_edges, bins.calc_nBins()
        return [np.linspace(bin_edges[i], bin_edges[i+1], supersample) for i in range(nBins)]
    
    def fill_bins(self, bins, sumW2, supersample=100, area=0):
        bin_widths= bins.calc_bin_widths()
        supersampled_bins = self.super_sample(bins, supersample)
        percent_weight = bin_widths/np.sum(bin_widths)
        supersampled_bin_values = list(map(lambda x: self.scaled_pdf(x, sumW2, area=area), supersampled_bins))
        return np.sum(supersampled_bin_values, axis=1)*percent_weight

    def tail_sys(self, bins, sumW2, 
                 supersample=100, width=2.5, tail_percent=.1, constant_percent=.05, stat_unc=0,
                 **kwargs):
        y = self.fill_bins(bins, sumW2, supersample=supersample,
                           **kwargs)
        supersampled_bins = self.super_sample(bins, supersample)
        mu = self.mu()
        sigma = self.sigma() 
        sigma_from_mu = list(map(lambda x: abs(x-mu)/sigma, supersampled_bins))
        tail_array = list(map(lambda x: logistics_func(x, width, .5), sigma_from_mu))
        subsampled_systematics = np.mean(tail_array, axis=1)*tail_percent*y
        constant_percent =  y*constant_percent
        uncertainty = (subsampled_systematics ** 2 + constant_percent ** 2 + y*stat_unc) ** .5
        return y, uncertainty
    def make_hist(self, bins, sumW2, **kwargs):
        y, unc = self.tail_sys(bins, sumW2, **kwargs)
        return SysHist(y, -unc, +unc, y*0, bins.bin_edges)

def reset_params(data):
    mean = np.mean(data)
    std = np.std(data)
    mu.set_value(mean)
    sigma.set_value(std)
    alphal.set_value(2)
    nl.set_value(3)
    alphar.set_value(1.5)
    nr.set_value(6)
    return mean, std

def logistics_func(x, offset, widht):
    return 1/(1+10.**(-(x-offset)/widht))
    
def sys_func(x, p0, p1, p2):
  fit = ((p0/x.mass)**2+(p1)**2)**.5 + p2*x.dbs
  #fit = (p0-p1*x.dbs)*np.log(x.mass)*np.log(p4)+p2+p3*x.dbs
  return fit 

def sys_func_offset(x, offset, *args):
  return sys_func(x, *args) + offset


popt_sys_dict = {
    "2016":{
    'SR1': [-24.54276342,   0.13869324,  -0.04319781],
    'SR2': [-16.95213727,  0.14670199, -0.03403927]
    },
    "2017":{
    'SR1': [-25.66070814,   0.11570483,  -0.04050246],
    'SR2': [-18.77099018,   0.1226328,   -0.0306517]
    },
    "2018":{
    'SR1': [-4.11823379e+01,  1.38593459e-01, -3.94183979e-02],
    'SR2': [-2.68983655e+01,  1.19450655e-01, -1.55283528e-02]
    },

} 
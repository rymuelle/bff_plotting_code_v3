import zfit
from zfit import z
from zfit.models.physics import DoubleCB
from scipy.optimize import curve_fit
from zfit.models.physics import double_crystalball_func
from bff_signal_model.functions import *
import numpy as np

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

    def super_sample(self, bin_edges, supersample):
        nbins = len(bin_edges) - 1
        return [np.linspace(bin_edges[i], bin_edges[i+1], supersample) for i in range(nbins)]
    
    def fill_bins(self, bin_edges, sumW2, supersample=100, area=0):
        nbins = len(bin_edges) - 1
        supersampled_bins = self.super_sample(bin_edges, supersample)
        supersampled_bin_values = list(map(lambda x: self.scaled_pdf(x, sumW2/nbins, area=area), supersampled_bins))
        return np.sum(supersampled_bin_values, axis=1)

    def tail_sys(self, bin_edges, sumW2, 
                 supersample=100, width=2.5, tail_percent=.1, constant_percent=.05, stat_unc=0,
                 **kwargs):
        y = self.fill_bins(bin_edges, sumW2, supersample=supersample,
                           **kwargs)
        supersampled_bins = self.super_sample(bin_edges, supersample)
        mu = self.mu()
        sigma = self.sigma() 
        sigma_from_mu = list(map(lambda x: abs(x-mu)/sigma, supersampled_bins))
        tail_array = list(map(lambda x: logistics_func(x, width, .5), sigma_from_mu))
        subsampled_systematics = np.mean(tail_array, axis=1)*tail_percent*y
        constant_percent =  y*constant_percent
        uncertainty = (subsampled_systematics ** 2 + constant_percent ** 2 + y*stat_unc) ** .5
        return y, uncertainty
        
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
    
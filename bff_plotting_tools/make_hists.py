import numpy as np
from plotting_meta.plotting_meta import Bins, bins

class SysHist(Bins):
    def __init__(self, nominal, down, up, std, bin_edges):
        Bins.__init__(self, bin_edges)
        self.nominal = nominal
        self.down = down
        self.up = up
        self.std = std
        self.bins = bins
    def calc_ratio(self, divisor):
        return SysHist(self.nominal/divisor, 
            self.down/divisor, 
            self.up/divisor, 
            self.std/divisor, 
            self.bin_edges)
    def make_density_hist(self, scale=1):
        width = self.calc_bin_widths()*scale
        return self.calc_ratio(width)
    def draw(self, ax, color='blue', **kwargs):
        ax.errorbar(self.calc_bin_centers(), self.nominal, yerr=self.std, drawstyle='steps-mid',color=color, **kwargs)
        ax.fill_between(self.calc_bin_centers(), self.up+self.nominal, self.down+self.nominal, step='mid', alpha=.5,color=color)
    def calc_sum(self):
        return np.sum(self.nominal)

def make_hist(df, column, region, bin_edges=bins.bin_edges, weights='Weight', region_sys='nom', std=0):
    '''Produce a binned count from a histogram'''
    region = "{}_{}".format(region, region_sys)
    df_temp = df[df[region]==1]
    hist = np.histogram(df_temp[column],
                 bins=bin_edges,
                 weights=df_temp[weights])[0]
    if std:
        std_hist = np.histogram(df_temp[column],
                 bins=bin_edges,
                 weights=df_temp[weights]**2)[0]**.5
        return hist, std_hist
    return hist

def make_sys(df, column, regionname, bin_edges=bins.bin_edges): 
    '''Produces up down systematics, linearly added.'''
    nominal, std = make_hist(df, column, regionname, bin_edges=bin_edges, std=1)
    sys_array = []
    for jetcorr in ['jer', 'jesTotal']:
        cor = []
        sig_cor = []
        for direction in ['Up','Down']:
            reg_sys = '{}{}'.format(jetcorr, direction)
            sys_hist = make_hist(df, column, regionname, bin_edges=bin_edges, region_sys=reg_sys)
            cor.append(sys_hist - nominal) 
        cor = sorted(cor, key=lambda x: np.sum(x))
        sys_array.append(cor)   
    weightsys = ["Weight_Pu","Weight_BTag","Weight_PUID","Weight_PDF_ISRFSR_","Weight_MuonSF","Weight_ElectronSF"]
    for sys in weightsys:
        cor = []
        sig_cor = []
        for direction in ['Up','Down']:
            weightnamesys = sys+direction
            sys_hist = make_hist(df, column, regionname, bin_edges=bin_edges, weights=weightnamesys)
            cor.append(sys_hist - nominal) 
        cor = sorted(cor, key=lambda x: np.sum(x))
        sys_array.append(cor)

    #add in quadrature, but preserve signs:
    signs = np.sign(sys_array)
    square_keep_sign = np.power(sys_array,2)*signs
    sum_keep_sign = np.sum(square_keep_sign,axis=0)
    sum_signs = np.sign(sum_keep_sign)
    down, up = abs(np.sum(square_keep_sign,axis=0))**.5*sum_signs
    # one line to sort up/down bins
    # list(zip(*map(sorted,zip(up,down))))
    return SysHist(nominal, down, up, std, bin_edges)

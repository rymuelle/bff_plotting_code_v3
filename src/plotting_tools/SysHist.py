import numpy as np
from src.plotting_tools.Bins import Bins, bins
from bff_plotting_tools.utils import rebin
import uncertainties

class SysHist(Bins):
    def __init__(self, nominal, down, up, std, bin_edges, sys={}):
        Bins.__init__(self, bin_edges)
        self.nominal = nominal
        self.down = down
        self.up = up
        self.std = std
        self.bins = Bins(bin_edges)
        self.sys = sys
    def add_sys(self, key, down, up):
        self.sys[key] = [down, up]
    def sys_summary(self):
        nom = self.nominal.sum()
        sys_string = "nominal: {:.1f}:".format(nom)
        for sys, (down, up) in self.sys.items():
            sys_string+="\n\t{}: {:.2f} to {:.2f}".format(sys, np.sum(down)/nom, np.sum(up)/nom)
        return sys_string
    
    
    def sys_from_sys_dict(self):
        ''' updates up down values from the sys dict object'''
        up_list, down_list = [], []
        for sys, (down, up) in self.sys.items():
            up_list.append(up)
            down_list.append(down)
        def sum_in_quad(sys_list): return (np.stack(sys_list)**2).sum(axis=0)**.5
        up, down, = sum_in_quad(up_list), -sum_in_quad(down_list)
        self.up = up
        self.down = down        
    def sys_per_from_add_sys(self):
        sys_tot = 0
        nom = np.sum(self.nominal)
        for sys, (down, up) in self.sys.items():
            sys_tot += (np.sum(up)**2+np.sum(down)**2)/2
        return sys_tot**.5/nom
    def sys_pers(self):
        nom = np.sum(self.nominal)
        systematics = {}
        for sys, (down, up) in self.sys.items():
            sys_tot = (np.sum(up)**2+np.sum(down)**2)**.5/(2*nom)   
            systematics[sys.replace('Up','Comb')] = sys_tot
        systematics['tot'] = self.sys_per_from_add_sys()
        return systematics
    def sys_string(self):
        string = ""
        nom = np.sum(self.nominal)
        for sys, (down, up) in self.sys.items():
            string += "{} {:.2f} to {:.2f}\n".format(sys, (np.sum(down)/nom), (np.sum(up)/nom))
        string += "nom: {:.2f} sys: {:.2f}".format(nom, self.sys_per_from_add_sys())
        return string
    @classmethod
    def from_ufloats(cls, bins, ufloats):
        y  = np.array(list(map(lambda x: x.nominal_value,ufloats)))
        unc  = np.array(list(map(lambda x: x.std_dev,ufloats)))
        return cls(y, -y*0, y*0, unc, bins.bin_edges)
    def rebin(self, nGroup):
        def rebin1D(x): return rebin(x.shape(1,-1), (1,nGroup))
        nominal = rebin1D(self.nominal)
        down = rebin1D(self.down)
        up = rebin1D(up)
        var = std**2
        var = rebin1D(var)
        std = var**.5
        
    def reduce_range(self, top=np.inf, bottom=-np.inf):
        bin_edges = self.bin_edges[(self.bin_edges < top) & (self.bin_edges>bottom)]
        new_bins = Bins(bin_edges)
        centers = new_bins.calc_bin_centers()
        x = np.array(self.calc_bin_centers())
        select = np.in1d(x, centers)
        nom = self.nominal[select]
        up = self.up[select]
        down = self.down[select]
        std = self.std[select]
        bin_edges = self.bin_edges[(self.bin_edges < top) & (self.bin_edges>bottom)]
        return SysHist(nom, 
            down, 
            up, 
            std, 
            bin_edges)        
        
    def calc_ratio(self, divisor):
        sys = self.sys
        new_sys = {}
        for key, (down, up) in sys.items():
            new_sys[key] = [down/divisor, up/divisor]
        return SysHist(self.nominal/divisor, 
            self.down/divisor, 
            self.up/divisor, 
            self.std/divisor, 
            self.bin_edges, 
            sys=new_sys)
    def inverse_make_density_hist(self, scale=1):
        width = np.array(self.calc_bin_widths()*scale)
        return self.calc_ratio(1./width)
    def make_density_hist(self, scale=1):
        width = np.array(self.calc_bin_widths()*scale)
        return self.calc_ratio(width)
    def draw(self, ax, color='blue', error_scale=1, draw_sys=True, **kwargs):
        ax.errorbar(self.calc_bin_centers(), self.nominal, yerr=self.std*error_scale, drawstyle='steps-mid',color=color, **kwargs)
        if draw_sys: ax.fill_between(self.calc_bin_centers(), self.up+self.nominal, self.down+self.nominal, step='mid', alpha=.5,color=color)
    def calc_sum(self):
        return np.sum(self.nominal)

    def uncertainty_std_dev(self):
        return np.array([uncertainties.ufloat(nom, std) for nom, std in zip(self.nominal, self.std)])
    
    def normalize(self):
        return self.calc_ratio(self.calc_sum())
    def to_dict(self):
        return {'nom': self.nominal, 'up':self.up, 'down':self.down,
                'std':self.std, 'bins': self.bin_edges}
    
    @classmethod
    def from_dict(cls, hist_dict):
        down = hist_dict['down']
        up = hist_dict['up']
        down[np.isnan(down)]=0
        up[np.isnan(up)]=0
        return cls(hist_dict['nom'],hist_dict['down'],hist_dict['up'],
                   hist_dict['std'],hist_dict['bins'])
    
    def calc_total_unc(self):
        return (((self.up-self.down)/2)**2 + self.std**2)**.5
    def chisquared(hist1, hist2, count_var=True, ndof=0):
        diff = hist1.nominal-hist2.nominal
        #uc1 = hist1.calc_total_unc()
        #uc2 = hist2.calc_total_unc()
        #print(uc1.sum(), uc2.sum(), diff.sum())
        #chi2 = diff**2/(uc1**2+uc2**2)
        if count_var:
            uc1 = hist1.nominal
        else: 
            uc1 = hist1.calc_total_unc()**2
        return (((diff)**2/uc1)/(len(hist1.nominal)-ndof)).sum()
    
    
    def __add__(self,other):
        ''' down is an approx'''
        sys = self.sys
        sys2 = other.sys
        new_sys = {}
        for key, (down, up) in sys.items():
            new_sys[key] = [(sys[key]['down']**2 + sys2[key]['down']**2)**.5,
                            (sys[key]['up']**2 + sys2[key]['up']**2)**.5]
            
        assert (self.bin_edges == other.bin_edges).all(), "must be same binning" 
        nom = self.nominal + other.nominal
        std = (self.std**2 + other.std**2)**.5
        up = (self.up**2 + other.up**2)**.5
        down= -(self.down**2 + other.down**2)**.5
        return SysHist(nom, down, up, std, self.bin_edges, sys=new_sys)
    
    def __mul__(self, mul):
        sys = self.sys
        new_sys = {}
        for key, (down, up) in sys.items():
            new_sys[key] = [down*mul, up*mul]
        return SysHist(self.nominal * mul, 
            self.down * mul, 
            self.up * mul, 
            self.std * mul, 
            self.bin_edges,
            sys = new_sys)        

def make_hist(x, bin_edges=bins.bin_edges, weights=1, std=0):
    hist = np.histogram(x,
                 bins=bin_edges,
                 weights=weights)[0]
    if std:
        std_hist = np.histogram(x,
                 bins=bin_edges,
                 weights=weights**2)[0]**.5
        return hist, std_hist
    return hist

def isin(df,region, select_level=1): return df[df[region]>=select_level]


def make_sys_hist(mdf, column, reg, replace_dict = {}, bin_edges=bins.bin_edges,
                 ind_sys_hist=False, select_level=1):
    
    
    nom_mdf = isin(mdf,'{}_jet_nom_muon_corrected_pt_ele_pt'.format(reg), select_level=select_level)
    nom_hist, nom_std =  make_hist(
        nom_mdf[column+'_jet_nom_muon_corrected_pt_ele_pt'],
        weights=nom_mdf.Weight,
        std=1,
        bin_edges=bin_edges
    )
    sys_list = []
    sys_weight_columns = mdf.filter(regex='^Weight.*Up').columns
    sys_weight_columns_down = mdf.filter(regex='^Weight.*Down').columns
    
    sys_weights = {}
    for down, up in list(zip(sys_weight_columns_down, sys_weight_columns)):
        down_hist =  make_hist(
            nom_mdf[column+'_jet_nom_muon_corrected_pt_ele_pt'],
            weights=nom_mdf[down],
        bin_edges=bin_edges)-nom_hist
        up_hist =  make_hist(
            nom_mdf[column+'_jet_nom_muon_corrected_pt_ele_pt'],
            weights=nom_mdf[up],
        bin_edges=bin_edges)-nom_hist
        sorted_sys_hist = list(sorted([down_hist,up_hist], key=lambda x: np.sum(x)))
        sys_weights[up]=sorted_sys_hist
    
    region_sys_columns = mdf.filter(regex='{}.*Up'.format(reg)).columns
    region_sys_columns_down = mdf.filter(regex='{}.*Down'.format(reg)).columns
    
    sys_region = {}
    for down, up in list(zip(region_sys_columns_down, region_sys_columns)):
        column_down, column_up = down.replace(reg, column), up.replace(reg, column)
        sys_mdf = isin(mdf,down, select_level=select_level)
        down_hist =  make_hist(
            sys_mdf[column_down],
            weights=sys_mdf.Weight,
        bin_edges=bin_edges)-nom_hist
        sys_mdf = isin(mdf,up, select_level=select_level)
        up_hist =  make_hist(
            sys_mdf[column_up],
            weights=sys_mdf.Weight,
        bin_edges=bin_edges)-nom_hist
        sorted_sys_hist = list(sorted([down_hist,up_hist], key=lambda x: np.sum(x)))
        sys_region[up]=sorted_sys_hist
        
    all_sys = {**sys_region, **sys_weights}
    #def replace bad sys for some years
    for sys in all_sys.keys():
        for w,v in replace_dict.items():
            if w in sys: 
                all_sys[sys] = [
                    nom_hist*(-v),
                    nom_hist*(v)
                ]
                
    all_sys_arr = np.array([x for _, x in all_sys.items()])
    sign = (all_sys_arr+1e-4)/abs(all_sys_arr+1e-4)
    all_sys_arr = (all_sys_arr**2).sum(axis=0)**.5
    all_sys_arr[0] = all_sys_arr[0]*-1
    if ind_sys_hist:
        syshist = SysHist(nom_hist, all_sys_arr[0], all_sys_arr[1], nom_std, bin_edges), all_sys
    else: syshist = SysHist(nom_hist, all_sys_arr[0], all_sys_arr[1], nom_std, bin_edges, sys={})
    for key, (down, up) in all_sys.items():
        syshist.add_sys(key, down, up)
    return syshist


def make_sys_hist_v2(*args, **kwargs):
    return make_sys_hist(*args, **kwargs)

#def era_corrected_sys_hist(era, *args, **kwargs):
#    if era=='2016': kwargs['replace_dict'] = {'Weight_MuonSF':.16}
#    return make_sys_hist(*args, **kwargs)
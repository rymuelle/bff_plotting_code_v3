import numpy as np

class Bins():
    def __init__(self, bin_edges):
        self.bin_edges = np.array(bin_edges)
    @classmethod
    def from_bin_stack(cls, bin_stack):
        bin_edges = np.unique(bin_stack.reshape(-1))
        return cls(bin_edges)
    def calc_bin_stack(self):
        return np.array([(self.bin_edges[i],self.bin_edges[i+1]) for i in range(self.calc_nBins())])
    def rebin(self):
        grouped_bins = self.calc_bin_stack().reshape(-1,2,2)
        rebinned_grouped_bins = np.array(list(map(lambda x: (np.min(x), np.max(x)), grouped_bins)))
        return Bins.from_bin_stack(rebinned_grouped_bins)
    def calc_nBins(self):
        return len(self.bin_edges)-1
    def calc_bin_centers(self):
        return [(self.bin_edges[i]+self.bin_edges[i+1])/2 for i in range(self.calc_nBins())]
    def calc_bin_widths(self):
        return [(self.bin_edges[i+1]-self.bin_edges[i]) for i in range(self.calc_nBins())] 
    def calc_bin_range(self):
        return (self.bin_edges[0], self.bin_edges[-1])
    def __repr__(self):
        return "{}".format(self.bin_edges)


def make_bins(binning_type="split", minVal=110, maxVal=800):
    if binning_type=="constant":
        bin_edges = np.linspace(minVal,maxVal, int((maxVal-minVal)/5 + 1))
    elif binning_type=="split":
        # 110-200 1 gev, 200-800, 5 gev. Skip the first bin edge of the latter
        bin_edges = np.append(np.linspace(minVal,200, int((200-minVal)/1+1)),
                         np.linspace(200,maxVal, int((maxVal-200)/5+1))[1:])
    elif binning_type=="linear":
        min_width = 1
        max_width = 5
        total_range = maxVal-minVal
        nBins = int(total_range/(min_width/2+max_width/2))
        bin_widths = np.linspace(min_width, max_width, nBins)
        bin_edges = np.array([np.sum(bin_widths[:i]) for i in range(len(bin_widths)+1)]) + minVal
    elif binning_type=="200":
        bin_edges = np.linspace(150,250, 101)
    else: return 0
    return Bins(bin_edges)
    
binning_type = "split"
bins = make_bins(binning_type=binning_type, minVal=110)

color_cycle = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"]
lumi_dict = {"2016": 36.33 , "2017": 41.48 , "2018":59.83, "16-18":  137.65}
    
def cms_format_fig(era, ax, **kwargs):
    import mplhep as hep
    import matplotlib.pyplot as plt
    hep.set_style(hep.style.CMS)
    plt.rcParams.update({
        "text.usetex": True,
    })
    lumi = lumi_dict[era]
    lumi = hep.cms.label(ax=ax, lumi=lumi,year=era, **kwargs);
    
    
    
import math

vtb = 1.021
vts = 0.04
wma = .8768 # weak mixing angle cos(theta_w)

def trident(mass): return 540/mass * 1.306 * 10 ** -9 * 200 ** 2
def BS_BS(mass,vtb,vts,wma): 
    greater_than = 0.16
    greater_than = greater_than * (1 - .029 * math.log(mass / 1000)) ** -1
    greater_than = greater_than * (mass / 1000) ** 2
    greater_than = greater_than ** .5
    greater_than = greater_than * (vtb*vts) / (2 * wma) 
    return greater_than

def xs_ys(x,c):
    return c*(1/x)


def trident_y(mass, x):
    return xs_ys(x,trident(mass))

def BS_BS_y(mass, x):
    return xs_ys(x,BS_BS(mass, vtb,vts,wma))

def get_plot_values(mass):
    trident_dict = {350:trident(350),
        200:trident(200),
        500:trident(500)}
    bs_dict  = {350:BS_BS(350,vtb,vts,wma),
        200:BS_BS(200,vtb,vts,wma),
        500:BS_BS(500,vtb,vts,wma)}
    inc_xsection = {350:9,
        200:53,
        500:2.45}
    return_obj = [bs_dict,trident_dict,inc_xsection]
    return [x[mass] for x in return_obj]


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
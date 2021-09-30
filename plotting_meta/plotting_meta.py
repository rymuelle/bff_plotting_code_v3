import numpy as np

class Bins():
    def __init__(self, bin_edges):
        self.bin_edges = np.array(bin_edges)
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


def make_bins(binning_type="split"):
    if binning_type=="constant":
        bin_edges = np.linspace(110,800, int((800-110)/5 + 1))
    elif binning_type=="split":
        # 110-200 1 gev, 200-800, 5 gev. Skip the first bin edge of the latter
        bin_edges = np.append(np.linspace(120,200, int((200-120)/1+1)),
                         np.linspace(200,800, int((800-200)/5+1))[1:])
    elif binning_type=="linear":
        min_width = 1
        max_width = 5
        total_range = 800-110
        nBins = int(total_range/(min_width/2+max_width/2))
        bin_widths = np.linspace(min_width, max_width, nBins)
        bin_edges = np.array([np.sum(bin_widths[:i]) for i in range(len(bin_widths)+1)]) + 110
    elif binning_type=="200":
        bin_edges = np.linspace(150,250, 101)
    else: return 0
    return Bins(bin_edges)
    
binning_type = "split"
bins = make_bins(binning_type=binning_type)

color_cycle = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"]

def cms_format_fig(era, ax):
    import mplhep as hep
    import matplotlib.pyplot as plt
    hep.set_style(hep.style.CMS)
    plt.rcParams.update({
        "text.usetex": True,
    })
    lumi_dict = {"2016": 35.50, "2017": 41.85, "2018":58.88}
    lumi = lumi_dict[era]
    lumi = hep.cms.label(ax=ax, lumi=lumi,year=era);
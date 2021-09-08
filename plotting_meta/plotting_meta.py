import numpy as np

def calc_bin_center(bin_edges):
    return [(bin_edges[i]+bin_edges[i+1])/2 for i in range(len(bin_edges[:-1]))]

def calc_bin_widths(bin_edges):
    return [(bin_edges[i+1]-bin_edges[i]) for i in range(len(bin_edges[:-1]))]

def make_bins(binning_type="split"):
    if binning_type=="constant":
        bin_edges = np.linspace(110,800, int((800-110)/5 + 1))
    elif binning_type=="split":
        # 110-200 1 gev, 200-800, 5 gev. Skip the first bin edge of the latter
        bin_edges = np.append(np.linspace(110,200, int((200-110)/1+1)),
                         np.linspace(200,800, int((800-200)/5+1))[1:])
    elif binning_type=="linear":
        min_width = 1
        max_width = 5
        total_range = 800-110
        nBins = int(total_range/(min_width/2+max_width/2))
        bin_widths = np.linspace(min_width, max_width, nBins)
        bin_edges = np.array([np.sum(bin_widths[:i]) for i in range(len(bin_widths)+1)]) + 110
    else: return 0,0,0
    bin_centers = calc_bin_center(bin_edges)
    bin_widths = calc_bin_widths(bin_edges)
    x_range = (bin_edges[0],bin_edges[-1])
    return bin_edges, bin_centers, bin_widths, x_range
    



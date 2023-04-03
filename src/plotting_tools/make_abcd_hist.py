from src.plotting_tools.Bins import bins
from src.plotting_tools.SysHist import SysHist

def make_abcd_hist(_abcd_arrays, bins=bins):
    nom = _abcd_arrays['nom']
    std = _abcd_arrays['std']
    bin_centers = bins.calc_bin_centers()
    return SysHist(nom, -std, std, std*0, bins.bin_edges)
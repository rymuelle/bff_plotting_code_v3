import mplhep as hep
import matplotlib.pyplot as plt

def cms_style():
    hep.style.use(hep.style.CMS)
    # This sets up matplotlib to use latex
    plt.rcParams.update({
        "text.usetex": True,
        })
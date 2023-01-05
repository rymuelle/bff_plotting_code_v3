import mplhep as hep
import matplotlib.pyplot as plt

def cms_style():
    hep.style.use(hep.style.CMS)
    # This sets up matplotlib to use latex
    plt.rcParams.update({
        "text.usetex": True,
        })
    

lumi_dict = {"2016": 36.33 , "2017": 41.48 , "2018":59.83, "16-18":  137.65}
lumi_dict["201X"] =  lumi_dict["16-18"]  
lumi_dict["Run 2"] =  lumi_dict["16-18"]  

def cms_format_fig(era, ax, **kwargs):
    lumi = lumi_dict[era]
    lumi = hep.cms.label(ax=ax, lumi=lumi,year=era, **kwargs);
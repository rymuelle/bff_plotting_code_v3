from plotting_meta.plotting_meta import make_bins
import numpy as np
import hist 

def make_hist(
              df, 
              axismeta, axisname, 
              weightname,
              regionname,
              typename,
              variablename,
              uniquelabels=[]
             ):

    ldf = df[(df[regionname]==1) & (df.type==typename)]
    values = ldf[variablename]
    weights = ldf[weightname]
    labels = ldf['labels']
    if uniquelabels==[]:
        uniquelabels =  np.unique(labels)
    
    dr_hist = (hist.Hist.new
               .StrCat(uniquelabels, name = "labels")
               .Reg(*axismeta, name=axisname)
               .Weight())
    
    dr_hist.fill(labels, values, weight=weights)
    return dr_hist


def make_sys(df, column, regionname): 
	bins, bins_center, bins_widht, x_range  = make_bins()
    nominal, std = make_hist(df, column, regionname,bins=bins, std=1)
    sys_array = []
    for jetcorr in ['jer', 'jesTotal']:
        cor = []
        sig_cor = []
        for direction in ['Up','Down']:
            reg_sys = '{}{}'.format(jetcorr, direction)
            sys_hist = make_hist(df, column, regionname, region_sys=reg_sys,bins=bins)
            cor.append(sys_hist - nominal) 
        cor = sorted(cor, key=lambda x: np.sum(x))
        sys_array.append(cor)   
    weightsys = ["Weight_Pu","Weight_BTag","Weight_PUID","Weight_PDF_ISRFSR_","Weight_MuonSF","Weight_ElectronSF"]
    for sys in weightsys:
        cor = []
        sig_cor = []
        for direction in ['Up','Down']:
            weightnamesys = sys+direction
            sys_hist = make_hist(df, column, regionname, weights=weightnamesys,bins=bins)
            cor.append(sys_hist - nominal) 
        cor = sorted(cor, key=lambda x: np.sum(x))
        sys_array.append(cor)
    down, up = np.sum(sys_array,axis=0)
    # one line to sort up/down bins
    # list(zip(*map(sorted,zip(up,down))))
    return nominal, down, up, std
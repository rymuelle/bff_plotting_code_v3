import numpy as np


def make_hist(df, column, region, bins=np.linspace(110,800, 100), weights='Weight', region_sys='nom', std=0):
    '''Produce a binned count from a histogram'''
    region = "{}_{}".format(region, region_sys)
    df_temp = df[df[region]==1]
    hist = np.histogram(df_temp[column],
                 bins=bins,
                 weights=df_temp[weights])[0]
    if std:
        std_hist = np.histogram(df_temp[column],
                 bins=bins,
                 weights=df_temp[weights]**2)[0]**.5
        return hist, std_hist
    return hist

def make_sys(df, column, regionname, bins=np.linspace(110,800,100)): 
    '''Produces up down systematics, linearly added.'''
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

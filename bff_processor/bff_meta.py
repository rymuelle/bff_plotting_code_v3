import numpy as np

nJets = [1,2]
regions = ['SR{}','CR{}0','CR{}3','CR{}4']
weights = ['Weight_Pu{}',
 'Weight_BTag{}',
 'Weight_PUID{}',
 'Weight_PDF_ISRFSR_{}',
 'Weight_MuonSF{}',
 'Weight_ElectronSF{}']
jet_variations = ['jer{}', 'jesTotal{}']
var_types = ['Up','Down']
all_reg = [reg.format(nJet) for reg in regions for nJet in nJets]
weights = [[weight.format(v) for v in var_types] for weight in weights]
jet_variations = [[jv.format(v) for v in var_types] for jv in jet_variations]

BFF_cut_values ={1:{'HTLT':-18,'RelMET':.3125,'SBM':18.75}, 
                    2:{'HTLT':-10,'RelMET':.2625,'SBM':1.25}}

def bff_cuts_pd(df, reg, jv, bff_cv=BFF_cut_values):
    import re
    import pandas as pd
    nJets = int(re.findall('(?:SR|CR)([0-9])', reg)[0])
    _bff_cv = bff_cv[nJets]
    less_than, greater_than, equal_to = {},{},{}
    less_than['HTLT_'+jv] = _bff_cv['HTLT']
    less_than['RelMET_'+jv] = _bff_cv['RelMET']
    less_than['DiLepMass'] = 800
    greater_than['SBM_'+jv] = _bff_cv['SBM']
    greater_than['DiLepMass'] = 105
    equal_to[reg+'_'+jv] = 1
    equal_to['Flag_METFilters'] = 1
    return df.loc[(df[list(less_than)] < pd.Series(less_than)).all(axis=1) & 
                 (df[list(greater_than)] > pd.Series(greater_than)).all(axis=1) & 
                 (df[list(equal_to)] == pd.Series(equal_to)).all(axis=1)
                 ]
def band_cut(column, low=-np.inf,high=np.inf):
    import pandas as pd
    def cut_df(df):
        return df[(df[column] > low) & ( df[column]< high)]
    return cut_df

def band_cut2d(column1, column2, low=[0,-np.inf], high=[0,np.inf]):
    import pandas as pd 
    from bff_processor.utils import linear
    def cut_df(df):
        return df[(linear(df[column1],*low) < df[column2]) & (linear(df[column1],*high) > df[column2])]
    return cut_df    

def isin(string):
    import pandas as pd
    def cut_df(df):
        return df[df[string]==1]
    return cut_df

def preselection():
    import pandas as pd
    def cut_df(df):
        return isin('Flag_METFilters')(band_cut('DiLepMass',low=105,high=800)(df))
    return cut_df

#identity wrapper for df
def identity(*args, **kwargs):
    def wrapper(x):
        return x
    return wrapper


def make_bff_cuts():
    import yaml
    with open('../bff_cut.yml', 'r') as f:
        param_values = yaml.load(f, yaml.SafeLoader)
    popt_relmet_1 = param_values['popt_relmet_1']
    popt_relmet_2 = param_values['popt_relmet_2']
    RelMET_ff1 = lambda x: band_cut2d('DiLepMass','RelMET_nom', high=popt_relmet_1)(isin('SR1_nom')(x))
    RelMET_ff2 = lambda x: band_cut2d('DiLepMass','RelMET_nom', high=popt_relmet_2)(isin('SR2_nom')(x))
    
    popt_htlt_1 = param_values['popt_htlt_1']
    popt_htlt_2 = param_values['popt_htlt_2']
    HTLT_ff1 = lambda x: band_cut2d('DiLepMass','HTLT_nom', high=popt_htlt_1)(RelMET_ff1(x))
    HTLT_ff2 = lambda x: band_cut2d('DiLepMass','HTLT_nom', high=popt_htlt_2)(RelMET_ff2(x))
    
    popt_TMB_1 = param_values['popt_TMB_1']
    popt_TMB_2 = param_values['popt_TMB_2']
    TMB_ff1 = lambda x: band_cut2d('DiLepMass','TMB_nom', low=popt_TMB_1)(HTLT_ff1(x))
    TMB_ff2 = lambda x: band_cut2d('DiLepMass','TMB_nom', low=popt_TMB_2)(HTLT_ff2(x))
    
import pandas as pd
import numpy as np
from src.assets.lumi import lumi_dict

def get_data_from_csv(file):
    df = pd.read_csv(file)
    clean_data(df)
    return df


def clean_data(df):
    df.replace([-np.inf, np.inf],np.nan, inplace=True)
    weight_columns = df.filter(regex='Weight.*').columns
    df.loc[df.type=='data', weight_columns] = 1
    df.fillna(0, inplace=True)

def delta_r_cut(df, dr):
    deltar_columns = ['minGoodJetElDR_jet_nom_muon_corrected_pt_ele_pt', 'minGoodJetMuDR_jet_nom_muon_corrected_pt_ele_pt']
    return df[df[deltar_columns].min(axis=1) > dr]


def read_feather_parquet(path):
    import pyarrow.feather as feather
    import pyarrow.parquet as parquet
    import pandas as pd
    try:
         df = pd.read_parquet(path+'.parquet')
    except:
         df = feather.read_feather(path+'.feather')
    return df

import pyarrow.feather as feather
def get_data(era, path, df_filter=lambda x: x.DiLepMass_jet_nom_muon_corrected_pt_ele_pt>0, stitch_dy=1, verbose=0,
            blinded=True, filterdeltar=True):
    if verbose: print("loading")
    if era=='2016':
        lumi=lumi_dict['2016']
    if era=='2017':
        lumi=lumi_dict['2017']
    if era=='2018':
        lumi=lumi_dict['2018']
    if era=='16-18':
        lumi = lumi_dict['2016']+lumi_dict['2017']+lumi_dict['2018']
        df = read_feather_parquet('{}/data/combined_2016'.format(path,era))
        print(df.shape)
        df = df.append( read_feather_parquet('{}/data/combined_2017'.format(path,era)))
        print(df.shape)
        df = df.append( read_feather_parquet('{}/data/combined_2018'.format(path,era)))
        print(df.shape)
        print("loaded all df")
    else:
        df = read_feather_parquet('{}/data/combined_{}'.format(path,era))
    clean_data(df)
    df = df[df_filter(df)]
    if verbose: print("loaded")
    if stitch_dy:
        #4 dy samples, weight samples evenly, remove DYLL as it contains some ISR sys issues
        dy_sel = (df.category.str.contains('DY')) & (df.category != 'DYLL')
        dy_types = df[dy_sel].category.unique()
        ntypes = len(dy_types)
        weight_columns = df.filter(regex='^Weight').columns
        total_dy = np.sum(dy_sel)
        for dy_type in dy_types:
            current_dy_sel = (df.category == dy_type)
            nevents = np.sum(current_dy_sel)
            df.loc[current_dy_sel,weight_columns] *= 1.0/ntypes #nevents/total_dy
        # set dyll weights to 0:
        df.loc[df.category == 'DYLL', weight_columns] *= 0
    else:
        dy_sel = df.type.str.contains('bck_ext')
        df = df[~dy_sel]
        
    #correct some issues from skimming
    # tops seem to use a different ISRFSR format, until rerun, just don't use
    # does not affect results
    bad_pdf_isr_sys = ['mc_santitop','mc_stop',]
    df.loc[df.name.isin(bad_pdf_isr_sys) ,'Weight_ISRFSR_Up'] = df.loc[df.name.isin(bad_pdf_isr_sys)]['Weight']
    df.loc[df.name.isin(bad_pdf_isr_sys) ,'Weight_ISRFSR_Down'] = df.loc[df.name.isin(bad_pdf_isr_sys)]['Weight']
    # Muon Trigger SFs have the standard weight subtracted out accidentally
    df['Weight_MuonTriggerDown'] = df['Weight'] - df['Weight_MuonTriggerDown'] 
    df['Weight_MuonTriggerUp'] += df['Weight']
    
    # filter out delta r
    if filterdeltar: df = delta_r_cut(df, 0.4)
            
    # blind SR
    if blinded:
        in_sr = df.filter(regex='SR[1|2]').sum(axis=1) > 0
        is_data = df.type.str.contains('data')
        nob_blinded =  (is_data  & in_sr) != 1
        df = df[nob_blinded]
        

    return df, lumi
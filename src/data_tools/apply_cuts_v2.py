# contains functions to build dataframe for use, such as applying bff cuts
import numpy as np
import dill as pickle
import pandas as pd
import re
from src.general.functions import linear_old, heaviside, power_func
# new linear has opposite order
linear = linear_old


def apply_multiple_filters(df, filter_list):
    filter_prod = np.array(filter_list).prod(axis=0)
    return df[filter_prod==1]

with open('bff_cut.dill', 'rb') as f:
    param_dict = pickle.load(f)

def RelMET_filter(df, col, jv, *popt):
    return power_func(df['DiLepMass_{}'.format(jv)], *popt) > df[col]

def HTLT_filter(df, col, jv, *popt):
    return linear(df['DiLepMass_{}'.format(jv)], *popt) > df[col]

def TMB_filter(df, col, jv, *popt):
    return heaviside(df['DiLepMass_{}'.format(jv)], *popt) < df[col]

RelMET_filter1 = lambda df, col, jv: RelMET_filter(df, col, jv, *param_dict['popt_relmet_1'])
RelMET_filter2 = lambda df, col, jv: RelMET_filter(df, col, jv, *param_dict['popt_relmet_2'])

HTLT_filter1 = lambda df, col, jv: HTLT_filter(df, col, jv, *param_dict['popt_htlt_1'])
HTLT_filter2 = lambda df, col, jv: HTLT_filter(df, col, jv, *param_dict['popt_htlt_2'])

TMB_filter1 = lambda df, col, jv: TMB_filter(df, col, jv, *param_dict['popt_TMB_1'])
TMB_filter2 = lambda df, col, jv: TMB_filter(df, col, jv, *param_dict['popt_TMB_2'])




def bff1_value(df, jv): return np.array([RelMET_filter1(df, 'RelMET_{}'.format(jv), jv), HTLT_filter1(df, 'HTLT_{}'.format(jv), jv), TMB_filter1(df, 'TMB_{}'.format(jv))]).prod(axis=0)
def bff2_value(df, jv): return np.array([RelMET_filter2(df, 'RelMET_{}'.format(jv), jv), HTLT_filter2(df, 'HTLT_{}'.format(jv), jv), TMB_filter2(df, 'TMB_{}'.format(jv))]).prod(axis=0)

bff_value = {"1":bff1_value,"2":bff2_value}

def bff1_no_tmb_value(df, jv): return np.array([RelMET_filter1(df, 'RelMET_{}'.format(jv), jv)*2, HTLT_filter1(df, 'HTLT_{}'.format(jv), jv)*4]).sum(axis=0)
def bff2_no_tmb_value(df, jv): return np.array([RelMET_filter2(df, 'RelMET_{}'.format(jv), jv)*2, HTLT_filter2(df, 'HTLT_{}'.format(jv), jv)*4]).sum(axis=0)

bff_no_tmb_value = {"1":bff1_no_tmb_value,"2":bff2_no_tmb_value}



def is_muon_reg(reg):
    return ('SR' in reg) or ('CR10' in reg) or ('CR20' in reg)

def trigger_filter(df, reg, era):
    era = int(era)
    if era == 2016:
        if is_muon_reg(reg): triggers =  ['HLT_Mu50', 'HLT_TkMu50']
        else: triggers =  ['HLT_DoubleEle33_CaloIdL_MW'], #'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL']
    if era == 2017:
        if is_muon_reg(reg): triggers =  ['HLT_Mu50', 'HLT_OldMu100', 'HLT_TkMu100']
        else: triggers =  ['HLT_DoubleEle33_CaloIdL_MW', 'HLT_DoubleEle25_CaloIdL_MW']
    if era == 2018:
        if is_muon_reg(reg): triggers =  ['HLT_Mu50', 'HLT_OldMu100', 'HLT_TkMu100']
        else: triggers =  ['HLT_DoubleEle25_CaloIdL_MW']
    tdf = df[df[reg]==1]
    filter_list = np.full(df.shape[0], 0)
    for trig in triggers:
        filter_list = filter_list | (df[trig].to_numpy().reshape(-1)==1).astype(bool)
        #filter_list = filter_list | df[trig].to_numpy().reshape(-1)==1
    return filter_list & df[reg]==1

def reg_filter(df, reg): return df[df[reg]==1]


def process_sample(row, era, verbose=False, trigger_fix=False):
    '''Applies bff cuts optimized previously'''
    #get stuff for bff samples
    name = re.findall('data\/tw_[0-9]{4}_(.+).csv', row.file)[0]
    #open file and filter out events with bff selection
    df = pd.read_csv(row.file)
    #loop over regions
    for reg in df.filter(regex='(?:SR|CR)\d+_.+'):
        nJets, jv = re.findall('(?:SR|CR)(\d)\d*_(.+)', reg)[0]
        # selected events are ==1, events pre bff selection are .5
        # trigger=1, relmet==2, htlt ==4, mudr=8, eledr = 16
        bff = bff_no_tmb_value[nJets](df, jv)
        trigger = trigger_filter(df, reg, era)
        mudr = (df['minGoodJetMuDR_{}'.format(jv)] > 0.4)*8
        eldr = (df['minGoodJetElDR_{}'.format(jv)] > 0.4)*16
        df[reg] = trigger_filter(df, reg, era)
        df[reg] += bff
        df[reg] += mudr
        df[reg] += eldr
        df[reg] *= trigger_filter(df, reg, era)/31
        
    selected_events = df.filter(regex='(?:SR|CR)\d+_.+').sum(axis=1)>0
    in_at_least_one_region_post_bff = (df.filter(regex='(?:SR|CR)\d+_.+')==1).sum(axis=1)>0
    if verbose:
        print(name)
        print("\t{} remaining".format(  in_at_least_one_region_post_bff.mean()))
    df = df[selected_events]
    #remove unnamed column from index, probably a better way, but ok for now
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df.reset_index(inplace=True)
    df['name'] = row.sample_name
    df['category'] = row.category
    df['mass'] = row.mass
    df['type'] = row.type
    df['gmu'] = row.gmu
    df['gb'] = row.gb
    df['dbs'] = row.dbs
    df['xsec'] = row.xsec
    df = df.replace([np.inf, -np.inf], 0)
    # make hlt columns bool:
    for hlt in df.filter(regex='HLT+_.+'):
        df[hlt] = df[hlt].astype(bool)
    #trigger sfs are not relative to nominal, make them relative
    if trigger_fix:
        df['Weight_MuonTriggerDown'] = df['Weight'] - df['Weight_MuonTriggerDown'] 
        df['Weight_MuonTriggerUp'] += df['Weight']
    return df


def process_sample_from_file(file, verbose=False):
    '''Applies bff cuts optimized previously'''
    #get stuff for bff samples
    name = re.findall('data\/tw_[0-9]{4}_(.+).csv', file)[0]
    #open file and filter out events with bff selection
    df = pd.read_csv(file)
    #loop over regions
    for reg in df.filter(regex='(?:SR|CR)\d+_.+'):
        nJets, jv = re.findall('(?:SR|CR)(\d)\d*_(.+)', reg)[0]
        # selected events are ==1, events pre bff selection are .5
        df[reg] = df[reg]*(1+bff_no_tmb_value[nJets](df, jv))/2
    selected_events = df.filter(regex='(?:SR|CR)\d+_.+').sum(axis=1)>0
    in_at_least_one_region_post_bff = (df.filter(regex='(?:SR|CR)\d+_.+')==1).sum(axis=1)>0
    if verbose:
        print(name)
        print("\t{} remaining".format(  in_at_least_one_region_post_bff.mean()))
    df = df[selected_events]
    return df

if __name__=="__main__":
    print(param_dict)

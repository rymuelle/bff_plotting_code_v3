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

def RelMET_filter(df, col, *popt):
    return power_func(df.DiLepMass, *popt) > df[col]

def HTLT_filter(df, col, *popt):
    return linear(df.DiLepMass, *popt) > df[col]

def TMB_filter(df, col, *popt):
    return heaviside(df.DiLepMass, *popt) < df[col]

RelMET_filter1 = lambda df, col: RelMET_filter(df, col, *param_dict['popt_relmet_1'])
RelMET_filter2 = lambda df, col: RelMET_filter(df, col, *param_dict['popt_relmet_2'])

HTLT_filter1 = lambda df, col: HTLT_filter(df, col, *param_dict['popt_htlt_1'])
HTLT_filter2 = lambda df, col: HTLT_filter(df, col, *param_dict['popt_htlt_2'])

TMB_filter1 = lambda df, col: TMB_filter(df, col, *param_dict['popt_TMB_1'])
TMB_filter2 = lambda df, col: TMB_filter(df, col, *param_dict['popt_TMB_2'])

RelMET_SR1_nom = lambda df: apply_multiple_filters(df, [df.SR1_nom, RelMET_filter1(df, 'RelMET_nom')])
RelMET_SR2_nom = lambda df: apply_multiple_filters(df, [df.SR2_nom, RelMET_filter2(df, 'RelMET_nom')])

HTLT_SR1_nom = lambda df: apply_multiple_filters(df, [df.SR1_nom, RelMET_filter1(df, 'RelMET_nom'), HTLT_filter1(df, 'HTLT_nom')])
HTLT_SR2_nom = lambda df: apply_multiple_filters(df, [df.SR2_nom, RelMET_filter2(df, 'RelMET_nom'), HTLT_filter2(df, 'HTLT_nom')])

TMB_SR1_nom = lambda df: apply_multiple_filters(df, [df.SR1_nom, RelMET_filter1(df, 'RelMET_nom'), HTLT_filter1(df, 'HTLT_nom'), TMB_filter1(df, 'TMB_nom') ])
TMB_SR2_nom = lambda df: apply_multiple_filters(df, [df.SR2_nom, RelMET_filter2(df, 'RelMET_nom'), HTLT_filter2(df, 'HTLT_nom'), TMB_filter2(df, 'TMB_nom')])



def bff1_value(df, jv): return np.array([RelMET_filter1(df, 'RelMET_{}'.format(jv)), HTLT_filter1(df, 'HTLT_{}'.format(jv)), TMB_filter1(df, 'TMB_{}'.format(jv))]).prod(axis=0)
def bff2_value(df, jv): return np.array([RelMET_filter2(df, 'RelMET_{}'.format(jv)), HTLT_filter2(df, 'HTLT_{}'.format(jv)), TMB_filter2(df, 'TMB_{}'.format(jv))]).prod(axis=0)

bff_value = {"1":bff1_value,"2":bff2_value}

def bff1_no_tmb_value(df, jv): return np.array([RelMET_filter1(df, 'RelMET_{}'.format(jv))*5/8, HTLT_filter1(df, 'HTLT_{}'.format(jv))*3/8]).sum(axis=0)
def bff2_no_tmb_value(df, jv): return np.array([RelMET_filter2(df, 'RelMET_{}'.format(jv))*5/8, HTLT_filter2(df, 'HTLT_{}'.format(jv))*3/8]).sum(axis=0)

bff_no_tmb_value = {"1":bff1_no_tmb_value,"2":bff2_no_tmb_value}


bff_1 = lambda df, jv: apply_multiple_filters(df, [RelMET_filter1(df, 'RelMET_{}'.format(jv)), HTLT_filter1(df, 'HTLT_{}'.format(jv)), TMB_filter1(df, 'TMB_{}'.format(jv))])
bff_2 = lambda df, jv: apply_multiple_filters(df, [RelMET_filter2(df, 'RelMET_{}'.format(jv)), HTLT_filter2(df, 'HTLT_{}'.format(jv)), TMB_filter2(df, 'TMB_{}'.format(jv))])

bff_1_no_TMB = lambda df, jv: apply_multiple_filters(df, [RelMET_filter1(df, 'RelMET_{}'.format(jv)), HTLT_filter1(df, 'HTLT_{}'.format(jv))])
bff_2_no_TMB = lambda df, jv: apply_multiple_filters(df, [RelMET_filter2(df, 'RelMET_{}'.format(jv)), HTLT_filter2(df, 'HTLT_{}'.format(jv))])


def reg_filter(df, reg): return df[df[reg]==1]

def process_sample(row, verbose=False):
    '''Applies bff cuts optimized previously'''
    #get stuff for bff samples
    name = re.findall('data\/tw_[0-9]{4}_(.+).csv', row.file)[0]
    #open file and filter out events with bff selection
    df = pd.read_csv(row.file)
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
    return df

if __name__=="__main__":
    print(param_dict)

import dill as pickle
from bff_processor.utils import linear_old, heaviside, power_func, apply_multiple_filters
# new linear has opposite order
linear = linear_old
import numpy as np

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

if __name__=="__main__":
    print(param_dict)

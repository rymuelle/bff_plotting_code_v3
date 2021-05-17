#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from pandas import DataFrame
import numpy as np
from bff_processor.SampleManager import SampleManager, SampleManagerPlotting
from bff_processor.SampleStack import SampleStack
from bff_processor.utils import make_view
import matplotlib.pyplot as plt
import mplhep as hep
import pprint
import matplotlib.gridspec as gridspec
hep.set_style(hep.style.CMS)
import argparse
import pickle
import itertools

parser = argparse.ArgumentParser(description='Produced flattened hists')
parser.add_argument('-e','--era', help='Era of sample.', required=True)
parser.add_argument('-o','--outname', help='Out name for pickle', required=True)
parser.add_argument('-c','--column', help='Column to plot', default='DiLepMass')
parser.add_argument('--min', help='min bin val', default=105, type=float)
parser.add_argument('--max', help='max bin val', default=800, type=float)
parser.add_argument('--nBins', help='n bin val', default=(800-105)/5+1, type=int)
parser.add_argument('--no_bff', help="don't apply bff cuts", default=False, action='store_true')
parser.add_argument('-t','--threads', help='n threads boost', default=0, type=int)

args = parser.parse_args()

b_kwargs = {}
if args.threads:
    print("use {} threads".format(args.threads))
    b_kwargs['threads'] = args.threads

no_bff = args.no_bff
era = args.era
outname = args.outname
column = args.column
meta_bin = [int(args.nBins),args.min,args.max]

sample_list =[
['data/tw_{}_BFFZprimeToMuMu_M_200.csv'.format(era), 'BFFZprimeToMuMu_M_200', 'BFFZprimeToMuMu_M_200', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_200_dbs0p5.csv'.format(era), 'BFFZprimeToMuMu_M_200_dbs0p5', 'BFFZprimeToMuMu_M_200_dbs0p5', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_200_dbs1p0.csv'.format(era), 'BFFZprimeToMuMu_M_200_dbs1p0', 'BFFZprimeToMuMu_M_200_dbs1p0', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_350.csv'.format(era), 'BFFZprimeToMuMu_M_350', 'BFFZprimeToMuMu_M_350', 'sig', 'BFFZprimeToMuMu_M_350', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_350_dbs0p5.csv'.format(era), 'BFFZprimeToMuMu_M_350_dbs0p5', 'BFFZprimeToMuMu_M_350_dbs0p5', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_350_dbs1p0.csv'.format(era), 'BFFZprimeToMuMu_M_350_dbs1p0', 'BFFZprimeToMuMu_M_350_dbs1p0', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_500.csv'.format(era), 'BFFZprimeToMuMu_M_500', 'BFFZprimeToMuMu_M_500', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_500_dbs0p5.csv'.format(era), 'BFFZprimeToMuMu_M_500_dbs0p5', 'BFFZprimeToMuMu_M_500_dbs0p5', 'sig', 'blue'],
['data/tw_{}_BFFZprimeToMuMu_M_500_dbs1p0.csv'.format(era), 'BFFZprimeToMuMu_M_500_dbs1p0', 'BFFZprimeToMuMu_M_500_dbs1p0', 'sig', 'blue'],

['data/tw_{}_ZToEE_M_120_200.csv'.format(era), 'ZToEE_M_120_200', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToEE_M_200_400.csv'.format(era), 'ZToEE_M_200_400', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToEE_M_400_800.csv'.format(era), 'ZToEE_M_400_800', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToEE_M_50_120.csv'.format(era), 'ZToEE_M_50_120', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToEE_M_800_1400.csv'.format(era), 'ZToEE_M_800_1400', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToMuMu_M_120_200.csv'.format(era), 'ZToMuMu_M_120_200', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToMuMu_M_200_400.csv'.format(era), 'ZToMuMu_M_200_400', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToMuMu_M_400_800.csv'.format(era), 'ZToMuMu_M_400_800', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToMuMu_M_50_120.csv'.format(era), 'ZToMuMu_M_50_120', 'DY', 'bkc', 'red'],
['data/tw_{}_ZToMuMu_M_800_1400.csv'.format(era), 'ZToMuMu_M_800_1400', 'DY', 'bkc', 'red'],

['data/tw_{}_mc_santitop.csv'.format(era), 'mc_santitop', 'ST', 'bkc', 'green'],
['data/tw_{}_mc_stop.csv'.format(era), 'mc_stop', 'ST', 'bkc', 'green'],
['data/tw_{}_mc_ttbar.csv'.format(era), 'mc_ttbar', 'TT', 'bkc', 'orange'],
['data/tw_{}_mc_ww.csv'.format(era), 'mc_ww', 'WW/WZ/ZZ', 'bkc', 'purple'],
['data/tw_{}_mc_wz.csv'.format(era), 'mc_wz', 'WW/WZ/ZZ', 'bkc', 'purple'],
['data/tw_{}_mc_zz.csv'.format(era), 'mc_zz', 'WW/WZ/ZZ', 'bkc', 'purple']
]
if era=='2017':
    sample_list.append(['data/tw_{}_BFFZprimeToMuMu_M_175.csv'.format(era), 'BFFZprimeToMuMu_M_175', 'BFFZprimeToMuMu_M_175', 'sig', 'blue'],)

sample_list_data =[
['data/tw_{}_data_el.csv'.format(era), 'data_el', 'data', 'data', 'black'],
['data/tw_{}_data_mu.csv'.format(era), 'data_mu', 'data', 'data', 'black'],
]

labels = [x[2] for x in sample_list]
labels = np.unique(labels)

stack = SampleStack(sample_list=sample_list)
stack_data = SampleStack(sample_list=sample_list_data)


weights = [
['Weight_PuUp','Weight_PuDown'],
['Weight_BTagUp','Weight_BTagDown'],
['Weight_PUIDUp','Weight_PUIDDown'],
['Weight_PDF_ISRFSR_Up','Weight_PDF_ISRFSR_Down'],
['Weight_MuonSFUp','Weight_MuonSFDown'],
['Weight_ElectronSFUp','Weight_ElectronSFDown'],]
regions = ['SR{}_{{}}','CR{}0_{{}}', 'CR{}3_{{}}', 'CR{}4_{{}}']
nJets = [1,2]
jes_var = [['jesDown', 'jesUp'], ['jerDown', 'jerUp']]
BFF_cuts =  ['HTLT_{}', 'RelMET_{}', 'SBM_{}']
def def_BFF_cuts(jv, nJet):
    t1 =  [bff.format(jv) for bff in BFF_cuts]
    return(list(zip(t1,BFF_cut_values[nJet])))
if no_bff:
    print("WARNING: no bff cuts applied")
    BFF_cut_values ={1:[np.inf, np.inf,0], 2:[np.inf,np.inf,0]}
else:
    BFF_cut_values ={1:[-18,.3125,18.75], 2:[-10,.2625,1.25]}
def make_view(jv,nJets, region):
    bff = def_BFF_cuts(jv, nJet)
    view_dict = {'lt': {'DiLepMass': 800, bff[0][0]: bff[0][1], bff[1][0]: bff[1][1]},
         'gt': {'DiLepMass': 105, bff[2][0]: bff[2][1]},
         'eq': {'Flag_METFilters': 1, region.format(nJets).format(jv): 1}}
    return view_dict

def make_view_no_reg():
    bff = def_BFF_cuts('nominal', 1)
    view_dict = {'lt': {'DiLepMass': 800, bff[0][0]: bff[0][1]+20, bff[1][0]: bff[1][1]+.2},
         'gt': {'DiLepMass': 105, bff[2][0]: 0},
         'eq': {'Flag_METFilters': 1}}
    return view_dict


#reduce num events
stack.view = make_view_no_reg()
stack = stack.SMP(clone=True)
stack_data.view = make_view_no_reg()
stack_data = stack_data.SMP(clone=True)


def make_sys_plots(stack, column, reg, nJet, ismc=0, **kwargs):
        n_jet_region  = reg.format(nJet)
        nom_region = n_jet_region.format('nominal')
        hist_dict = {}
        weight_dict = {}
        jv_dict = {}
        
        stack.view = make_view_no_reg()

        stack.view = make_view('nominal',nJet, reg)
        nom_stack = stack.SMP(clone=True)
        nom_hist = nom_stack.sum_boost(column, *meta_bin, **kwargs, w_kwargs={'weight_names':['Weight']}, b_kwargs=b_kwargs)
        if ismc: 
            for up,down in weights:
                up_hist = nom_stack.sum_boost(column, *meta_bin, **kwargs, w_kwargs={'weight_names':[up]}, b_kwargs=b_kwargs)
                down_hist = nom_stack.sum_boost(column, *meta_bin, **kwargs, w_kwargs={'weight_names':[down]}, b_kwargs=b_kwargs)
                weight_dict[up] = [up_hist,down_hist]
                
            for up,down in jes_var:
                stack.view = make_view(up,nJet, reg)
                up_hist = stack.sum_boost(column, *meta_bin, **kwargs, w_kwargs={'weight_names':['Weight']}, b_kwargs=b_kwargs)
                stack.view = make_view(down,nJet, reg)
                down_hist = stack.sum_boost(column, *meta_bin, **kwargs, w_kwargs={'weight_names':['Weight']}, b_kwargs=b_kwargs)
                jv_dict[up] = [up_hist,down_hist]
            
        hist_dict['nom'] = nom_hist
        if ismc:
            hist_dict['jv'] = weight_dict
            hist_dict['weights'] = jv_dict
        return hist_dict

reg_list =  [[reg, nJet] for reg in regions for nJet in nJets]

from time import perf_counter
reg_dict = {}
for reg, nJet in reg_list:
        start = perf_counter()
        hist_dict = {}
        if 'SR' not in reg:
            hist_dict["data"] = make_sys_plots(stack_data, column, reg, nJet)
        for label in labels:
            print(label)
            hist_dict[label] = make_sys_plots(stack, column, reg, nJet, ismc=1, label=label)

        reg_dict[reg.format(nJet)] = hist_dict



with open(outname,'wb') as f:
    pickle.dump(reg_dict,f)


with open(outname,'rb') as f:
    reg_dict_2016_temp = pickle.load(f)

import pprint
pprint.pprint(reg_dict_2016_temp)


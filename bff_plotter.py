#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
from glob import glob
from ROOT import RDataFrame, TFile, TH1F, TH2F
import ROOT
import yaml
from bff_processor.cpp_function import def_cpp
from bff_processor.utils import toVector, get_nEvents
from bff_processor.df_definitions import *
import pandas as pd
from time import perf_counter
import argparse

parser = argparse.ArgumentParser(description='Produce Root/CSV skims of nanoskims for plotting and analysis')
parser.add_argument('-e','--era', help='Era of sample.', required=True)
parser.add_argument('-y','--yaml', help='Yaml file with samples', required=True)
parser.add_argument('-m','--multi', help='Enable multithreading', default=False, action='store_true')

args = parser.parse_args()

print(args)

#define gint c++
def_cpp()

#for multithreading
if args.multi:
    ROOT.ROOT.EnableImplicitMT()
RDFrame = RDataFrame


# In[ ]:


#fname = "samplesCR_2016_Apr2020.yml"
#fname = "sampleCR_2017_v6_lo.yml"
#fname = "samplesCR_2018_v6_lo.yml"

era = str(args.era)
if era == "2016":
    bDiscValue = 0.6321
elif era == "2017":
    bDiscValue = 0.4941
elif era == "2018":
    bDiscValue = 0.4184
else: 
    assert (era in ['2016','2017','2018']), '{} not in 2016,2017,2018'.format(era)

fname = args.yaml
outname = fname.replace('.yml','.root')
outname

class sample_processor():
    def __init__(self,file_name,outname,bDiscValue,is_inclusive=0):
        #load config
        self.file_name = file_name
        with open(file_name,'r') as f:
            self.sample_dict = yaml.load(f, Loader=yaml.FullLoader)
        #setup outfile
        self.outname = outname
        self.out = TFile(outname, 'recreate')
        self.outdirs_dict = {}
        for sample in self.samples():
            name = sample['name']
            self.outdirs_dict[name] = self.out.mkdir(name)
        self.lumi = self.sample_dict['lumi']
        #get and write lumi info
        hlumi = TH1F("lumi", "lumi", 1, 0, 1)
        hlumi.SetDirectory(self.out)
        hlumi.SetBinContent(1, self.lumi)
        hlumi.Write()
        self.bDiscValue = bDiscValue
        self.is_inclusive = is_inclusive
    def samples(self):
        return self.sample_dict['samples']
    def sample_names(self):
        return [s['name'] for s in self.samples()]
    def close(self):
        self.out.Close()
    def __repr__(self):
        text_dict = {"fn":self.file_name,
                     "on":self.outname, 
                     "lumi":self.lumi,
                    "samples": self.sample_names()}
        return '''from {fn} to {on}\nlumi: {lumi}\nSamples {samples}'''.format(**text_dict)

sp = sample_processor(fname, outname, bDiscValue)
print(sp)


def create_regions(df, ismc):
    # create regions
    if int(ismc):
        JERC_var = ['nominal','jerUp','jerDown','jesUp','jesDown']
    else:
        JERC_var = ['nominal']
    rs = ["CR10", "CR11", "CR12", "CR13", "CR14", "CR20", "CR21", "CR22", "CR23", "CR24", "SR1", "SR2"]
    rs = [(r,var) for r in rs for var in JERC_var]
    for reg,var in rs:
        r = '{}_{}'.format(reg,var)
        HTLT_string = 'HTLT_{}'.format(var)
        RelMET_string = 'RelMET_{}'.format(var)
        SBM_string = 'SBM_{}'.format(var)

        HTLT,RelMET,SBM = -120,0.22,0
        if r[2] == "2":
            HTLT,RelMET,SBM = -60,0.22,150
        
        region_string = "{}pre_bff".format(r)
        format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
        regions.append((df.Filter("{r} && DiLepMass>54".format(**format_dict), region_string), format_dict))

        region_string = "{}".format(r)
        format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
        regions.append((df.Filter("{r} && DiLepMass>54 && {HTLT_string}<{HTLT} && {RelMET_string}<{RelMET} && {SBM_string}>{SBM}".format(**format_dict), region_string), format_dict))

        region_string = "{}_200_GeV_htlt_sig".format(r)
        format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
        regions.append((df.Filter("{r} && DiLepMass>54 && {RelMET_string}<{RelMET} && {SBM_string}>{SBM}".format(**format_dict), region_string), format_dict))

        region_string = "{}_200_GeV_sig".format(r)
        format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
        regions.append((df.Filter("{r} && DiLepMass>54 && {RelMET_string}<{RelMET} && {SBM_string}>{SBM}".format(**format_dict), region_string), format_dict))

        
columns = ['Weight_PuUp','Weight_PuDown','Weight_BTagUp','Weight_BTagDown','Weight_PUIDUp','Weight_PUIDDown','Weight_PDF_ISRFSR_Up','Weight_PDF_ISRFSR_Down','Weight_MuonSFUp','Weight_MuonSFDown','Weight_ElectronSFUp','Weight_ElectronSFDown','Weight','sample_weight','DiLepMass','TriggerWeight','Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_eeBadScFilter','Flag_METFilters','HTLT_nominal','RelMET_nominal','SBM_nominal','SR2_nominal','SR1_nominal','CR10_nominal','CR11_nominal','CR12_nominal','CR13_nominal','CR14_nominal','SR2_nominal','CR20_nominal','CR21_nominal','CR22_nominal','CR23_nominal','CR24_nominal','HTLT_jerDown','RelMET_jerDown','SBM_jerDown','SR2_jerDown','SR1_jerDown','CR10_jerDown','CR11_jerDown','CR12_jerDown','CR13_jerDown','CR14_jerDown','SR2_jerDown','CR20_jerDown','CR21_jerDown','CR22_jerDown','CR23_jerDown','CR24_jerDown','HTLT_jerUp','RelMET_jerUp','SBM_jerUp','SR2_jerUp','SR1_jerUp','CR10_jerUp','CR11_jerUp','CR12_jerUp','CR13_jerUp','CR14_jerUp','SR2_jerUp','CR20_jerUp','CR21_jerUp','CR22_jerUp','CR23_jerUp','CR24_jerUp','HTLT_jesDown','RelMET_jesDown','SBM_jesDown','SR2_jesDown','SR1_jesDown','CR10_jesDown','CR11_jesDown','CR12_jesDown','CR13_jesDown','CR14_jesDown','SR2_jesDown','CR20_jesDown','CR21_jesDown','CR22_jesDown','CR23_jesDown','CR24_jesDown','HTLT_jesUp','RelMET_jesUp','SBM_jesUp','SR2_jesUp','SR1_jesUp','CR10_jesUp','CR11_jesUp','CR12_jesUp','CR13_jesUp','CR14_jesUp','SR2_jesUp','CR20_jesUp','CR21_jesUp','CR22_jesUp','CR23_jesUp','CR24_jesUp']

columns_data = ['Weight','sample_weight','DiLepMass','TriggerWeight','Flag_goodVertices','Flag_globalSuperTightHalo2016Filter','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadPFMuonFilter','Flag_eeBadScFilter','Flag_METFilters','HTLT_nominal','RelMET_nominal','SBM_nominal','SR2_nominal','SR1_nominal','CR10_nominal','CR11_nominal','CR12_nominal','CR13_nominal','CR14_nominal','SR2_nominal','CR20_nominal','CR21_nominal','CR22_nominal','CR23_nominal','CR24_nominal',]
        
def process_sample(sp,sample,era,verbose=1):
    #get metadata
    name,xsec,nevts = sample['name'],sample['xsec'],sample['nevts']
    ismc,fileglob,bTagEff = int(sample['ismc']),sample['fileglob'],sample['bTagEff']
    if not nevts:
        nevts = get_nEvents(fileglob)
    sample_weight = float(xsec)*sp.lumi/float(nevts)
    if verbose: print("name: {} , xsec: {}, nevents: {} ismc: {}".format(name,xsec,nevts,ismc))
    #make file glob
    files = toVector('string', glob(fileglob))
    #set up btagging and puid sf files
    bTagFile, PUIDSFfile = setup_btag_puid(ismc, era, bTagEff)
    #make rdf
    df = RDFrame('Events', files)
    df = df.Filter("DiLepMass>54", "mass_cut")
    df = df.Filter("SR1_nominal or CR10_nominal or CR11_nominal or CR12_nominal or CR13_nominal or CR14_nominal or SR2_nominal or CR20_nominal or CR21_nominal or CR22_nominal or CR23_nominal or CR24_nominal", "in_region")
    df = def_good_jet(df,ismc, bDiscValue)
    df = def_good_leptons(df, ismc)
    df = def_HLT(df, ismc, era)
    df = def_sf_and_weight(df,ismc, sp.is_inclusive, name, sample_weight)
    df = def_lep_selections(df)
    lcolumns = columns
    if not ismc:
        lcolumns = columns_data
    df.Snapshot("Events", "data/tw_{}_{}.root".format(era,name), lcolumns)
    df_np = df.AsNumpy(lcolumns)
    df_df = pd.DataFrame(df_np)
    df_df.to_csv('data/tw_{}_{}.csv'.format(era,name))
    return name,df, fileglob

#####
##### BTag SF not included
#####
for sample in sp.samples():
    start_time = perf_counter()
    name,df,fileglob = process_sample(sp,sample,era)
    #count = df.Count()
    end_time = perf_counter()
    try:
        print(name,end_time-start_time)
    except:
        print(name,end_time-start_time)
sp.close()

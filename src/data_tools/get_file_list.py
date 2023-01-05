import pandas as pd
import numpy as np
import os
import re
from src.data_tools.make_sample_df import make_sample_df


def get_value_from_df(row, sample_df, column):
    file = row.file
    matched_row = sample_df[sample_df.file==file]
    if matched_row.shape[0]==1: return matched_row[column].values[0]
    else: 
        return -1

def regex_select(era):
    return "tw_{}.+\.csv".format(era)

def one_in(list_of_strings, test_string):
    for string_in in list_of_strings:
        if string_in in test_string: return True
    return False

def get_file_df(regex_select='tw_(?:2016|2017|2018).+\.csv', path='data'):
    files = ["{}/{}".format(path,x) for x in os.listdir(path) if re.match(regex_select, x)]
    DY = [{'type': 'bck', 'category': 'DY', 'file': x} for x in files if re.match('.+ZTo(?:Mu|EE).+', x)]
    ST = [{'type': 'bck','category':'ST', 'file':x} for x in files if re.match('.+top.csv', x)]
    VB = [{'type': 'bck','category':'VB', 'file':x} for x in files if re.match('.+mc_(?:ww|wz|zz)', x)]
    TT = [{'type': 'bck','category':'TT', 'file':x} for x in files if re.match('.+ttbar', x)]
    BFF = [{'type': 'sig','category':'BFF', 'file':x} for x in files if re.match('.+BFFZp', x)]
    Y3 = [{'type': 'sig','category':'Y3', 'file':x} for x in files if re.match('.+y3.+', x)]
    data = [{'type': 'data','category':'data', 'file':x} for x in files if re.match('.+_data_', x)]    
    ST_extra  = [{'type': 'bck_ext', 'category': 'ST_extra', 'file': x} for x in files if re.match('.+ST_.+', x)]
    TB  = [{'type': 'bck_ext', 'category': 'TB', 'file': x} for x in files if re.match('.+[W,Z]{3}.+', x)]
    WJ  = [{'type': 'bck_ext', 'category': 'WJ', 'file': x} for x in files if re.match('.+_WJets.+', x)]
    TTV  = [{'type': 'bck_ext', 'category': 'TTV', 'file': x} for x in files if re.match('.+TT[W,Z].+', x)]
    higgs  = [{'type': 'bck_ext', 'category': 'higgs', 'file': x} for x in files if one_in(['ggH','VBF', 'Wp', 'Wm', 'ZH', 'ttH', 'GluGluHToMuMu'],x)]
    DYLL  = [{'type': 'bck_ext', 'category': 'DYLL', 'file': x} for x in files if re.match('.+DYJetsToLL.+[0-9]+to[0-9]+.+', x)]
    DYAMad  = [{'type': 'bck_ext', 'category': 'DYAMad', 'file': x} for x in files if re.match('.+DYJetsToLL.+M-50_.+mad.+', x)]
    DYAMC  = [{'type': 'bck_ext', 'category': 'DYAMC', 'file': x} for x in files if re.match('.+DYJetsToLL.+M-50_.+amc.+', x)]
    all_lists = Y3 + DY + ST + VB + TT + BFF + data + ST_extra + TB + WJ + TTV + higgs + DYLL + DYAMad + DYAMC
    all_list_list = [ x['file'] for x in all_lists]
    not_in_all_lists = list(set(files) - set(all_list_list))
    if len(not_in_all_lists)!=0: print("not in all lists: {}".format(not_in_all_lists))
    assert len(files) == len(all_lists), "duplicate or uncaught file nfiles: {} vs ncaptured: {}".format(len(files), len(all_lists))
    df = pd.DataFrame(all_lists)
    df['mass'] = df.file.apply(get_mass)
    df['dbs'] = df.file.apply(get_dbs)
    df['gmu'] = df.file.apply(get_gmu)
    df['gb'] = df.apply(get_gb, axis=1)
    df['era'] = df.file.apply(lambda x: int(re.findall('.*tw_([0-9]+)_.*', x)[0]))
    #get xsec
    sample_df = make_sample_df()
    df['xsec'] = df.apply(lambda x: get_value_from_df(x,sample_df, 'xsec'), axis=1)
    df['sample_name'] = df.apply(lambda x: get_value_from_df(x,sample_df, 'name'), axis=1)
    return df


def linearCut2d(x, y, x0, y0, x1, y1, lessThan=False, greaterThan=False):
    '''x, y, x0, y0, x1, yx
    make a line between x0,y0 and x1, y1
    return y above this line
    (y-y0) > m*(x-x0)'''
    m = (y1-y0)/(x1-x0)
    if greaterThan and not lessThan:
        return  y-y0 > m*(x-x0)
    if lessThan and not greaterThan:
        return  y-y0 < m*(x-x0)
    
    
def get_mass(string):
    if not 'BFF' in string: return None
    mass = re.findall('.*M_([0-9]+).*', string)
    if len(mass)==1: return float(mass[0])
    else: return None
    
    
def get_dbs(string):
    if not 'BFF' in string: return None
    value = re.findall('.*dbs([0-9,p]+).*', string)
    if len(value)==1: return float(value[0].replace('p','.'))
    else: return None
    
def get_gmu(string):
    if not 'BFF' in string: return None
    value = re.findall('.*gmu([0-9,p]+).*', string)
    if len(value)==1: return float(value[0].replace('p','.'))
    return 0.17
   #else: 
   #    mass = get_mass(string)
   #    #experimental linear fit to paper values of gmu (.08, .14, .2) to (200, 350, 500)
   #    return 4.00000000e-04*mass 
    
    
def get_gb(row):
    if row.type!="sig": return None
    value = re.findall('.*gb([0-9,p]+).*', row.file)
    if len(value)==1: return float(value[0].replace('p','.'))
    return 0.02
    #else: 
    #    constant = 1.3*1e-5
    #    return constant*(row.mass/100)**2/(row.gmu*row.dbs)  
           
  
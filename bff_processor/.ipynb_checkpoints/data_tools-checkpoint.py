import pandas as pd
import numpy as np
import os
import re

def regex_select(era):
    return "tw_{}.+\.csv".format(era)

def get_files(regex_select, path='data'):
    files = ["{}/{}".format(path,x) for x in os.listdir(path) if re.match(regex_select, x)]
    DY = [x for x in files if re.match('.+ZTo(?:Mu|EE).+', x)]
    ST = [x for x in files if re.match('.+top.csv', x)]
    VB = [x for x in files if re.match('.+mc_(?:ww|wz|zz)', x)]
    TT = [x for x in files if re.match('.+ttbar', x)]
    BFF = [x for x in files if re.match('.+BFFZp', x)]
    data = [x for x in files if re.match('.+_data_', x)]
    assert len(files) == len(DY+ST+TT+VB+data+BFF), "duplicate or uncaught file"
    return {'DY': DY,'ST': ST,'VB': VB,'TT': TT,'BFF': BFF,'data': data}

def make_df(df_list, verbose=0):
    t_df = pd.DataFrame()
    for fname in df_list:
        if verbose: print(fname)
        t_df = pd.concat([t_df, pd.read_csv(fname)])
    return t_df

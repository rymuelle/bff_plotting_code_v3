import pandas as pd
import numpy as np
import os
import re

def regex_select(era):
    return "tw_{}.+\.csv".format(era)

def get_file_dict(regex_select='tw_(?:2016|2017|2018).+\.csv', path='data'):
    files = ["{}/{}".format(path,x) for x in os.listdir(path) if re.match(regex_select, x)]
    DY = [{'category': 'bck', 'type': 'DY', 'file': x} for x in files if re.match('.+ZTo(?:Mu|EE).+', x)]
    ST = [{'category': 'bck','type':'ST', 'file':x} for x in files if re.match('.+top.csv', x)]
    VB = [{'category': 'bck','type':'VB', 'file':x} for x in files if re.match('.+mc_(?:ww|wz|zz)', x)]
    TT = [{'category': 'bck','type':'TT', 'file':x} for x in files if re.match('.+ttbar', x)]
    BFF = [{'category': 'sig','type':'BFF', 'file':x} for x in files if re.match('.+BFFZp', x)]
    data = [{'category': 'data','type':'data', 'file':x} for x in files if re.match('.+_data_', x)]
    assert len(files) == len(DY+ST+TT+VB+data+BFF), "duplicate or uncaught file"
    df = pd.DataFrame(DY+ST+TT+VB+data+BFF)
    df['mass'] 
    return 


def read_csv(row):
    return pd.read_csv(row.file)

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
    mass = re.findall('.*M_([0-9]+).*', string)
    if len(mass)==1: return int(mass[0])
    else: return None
from glob import glob
from ROOT import vector
import os
import uproot
import numpy as np

def toVector(vtype, list):
    v = vector(vtype)()
    for o in list:
        v.push_back(o)
    return v

def get_nEvents(fileglob, key = 'sumPosMinusNegGenWeight'):
    file_path = fileglob.replace('*.root','')
    files = ["{}/{}".format(file_path,f) for f in os.listdir(file_path) if '.root' in f]
    total = 0
    for file in files:
        up_f = uproot.open(file)
        total += up_f[key].numpy()[0][0]
    return total

def make_view(mass_cut = [-np.inf,np.inf], HTLT = np.inf, RelMET = np.inf, SBM = 0, MET_filter = 1, region=0):
    view =  {
        'lt':{
            'DiLepMass': mass_cut[1],
            'HTLT_nominal': HTLT,
            'RelMET_nominal': RelMET
        },
        'gt':{
            'DiLepMass': mass_cut[0],
            'SBM_nominal': SBM
        },
        'eq':{
            'Flag_METFilters': 1,
        }
    }
    if region: view['eq'][region] = 1
    return view


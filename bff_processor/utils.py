from glob import glob
from ROOT import vector
import os
import uproot

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

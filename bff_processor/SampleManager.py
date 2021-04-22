import pandas as pd
from pandas import DataFrame
import numpy as np

class SampleManager(DataFrame):
    Flag_METFilters = {"Flag_goodVertices":1,
        "Flag_globalSuperTightHalo2016Filter":1,
        "Flag_HBHENoiseFilter":1,
        "Flag_HBHENoiseIsoFilter":1,
        "Flag_EcalDeadCellTriggerPrimitiveFilter":1,
        "Flag_BadPFMuonFilter":1,
        "Flag_eeBadScFilter":1,
        "Flag_METFilters":1}
    @classmethod
    def merge(cls,df_list):
        name = ''
        for i, df in enumerate(df_list):
            if i == 0: 
                df_temp = df
                name += '{}'.format(df.name)
            else: 
                df_temp = pd.merge(df_temp, df, 'outer')
                name += ' + {}'.format(df.name)
        return cls(df_temp, name=name)
    def __init__(self, df, name="", *args, **kwargs):
        super().__init__(*args, **kwargs)
        for key in df:
            if "Unnamed" in key: continue
            self[key] = df[key]
        self.name = name
    @classmethod
    def from_file(cls, file_name, *args, **kwargs):
        df = pd.read_csv(file_name)
        return cls(df, name=file_name, *args, **kwargs)
    def view(self, lt={}, gt={}, eq={}):
        df = self
        for key in lt:
            df = df[df[key] < lt[key]] 
        for key in gt:
            df = df[df[key] > gt[key]]
        for key in eq:
            df = df[df[key] == eq[key]]
        return self.__class__(df)
    def mass_window(self,min_mass,max_mass):
        return self.view(gt={'DiLepMass':min_mass},lt={'DiLepMass':max_mass})
    def bff(self,HTLT=np.inf,RelMET=np.inf,SBM=0):
        return self.view(gt={'SBM_nominal':SBM},lt={'RelMET_nominal':RelMET,'HTLT_nominal':HTLT})
    def met_filter(self):
        return self.view(eq=Flag_METFilters)
    def flag(self,flag_name, istrue=1):
        return self.view(eq={flag_name:istrue})
    def total(self):
        return np.sum(self.sample_weight)
    def histogram(self,column_name, **kwargs):
        pd_series = self[column_name]
        weights = self.sample_weight
        n,bins  = np.histogram(pd_series,weights=weights,**kwargs)
        stat_error = np.sqrt(np.histogram(pd_series,weights=weights**2, bins=bins)[0])
        bin_centers = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
        return n, stat_error, bins, bin_centers
    def __repr__(self):
        mean,std = self.DiLepMass.mean(),self.DiLepMass.std()
        return '''Name: {} \n shape: {}\n total: {:.2e}\n mass mean: {:.2f} std: {:.2f}'''.format(self.name,self.shape,self.total(),mean,std)
    def print(self):
        print(repr(self))

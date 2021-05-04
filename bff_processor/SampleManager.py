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
    def view(self, lt={}, gt={}, eq={}, **kwargs):
        df = self
        for key in lt:
            df = df[df[key] < lt[key]] 
        for key in gt:
            df = df[df[key] > gt[key]]
        for key in eq:
            df = df[df[key] == eq[key]]
        return self.__class__(df, **kwargs)
    def mass_window(self,min_mass,max_mass):
        return self.view(gt={'DiLepMass':min_mass},lt={'DiLepMass':max_mass})
    def bff(self,HTLT=np.inf,RelMET=np.inf,SBM=0):
        return self.view(gt={'SBM_nominal':SBM},lt={'RelMET_nominal':RelMET,'HTLT_nominal':HTLT})
    def met_filter(self):
        return self.view(eq=Flag_METFilters)
    def flag(self,flag_name, istrue=1):
        return self.view(eq={flag_name:istrue})
    def weights(self,weight_names=['Weight']):
        #this is really weird, but if you return a df with only one row in it, this forms a 1d list, unless you conver to an iterable like numpy array
        weights = [self[weight].to_numpy() for weight in weight_names]
        return np.prod(weights,axis=0)
    def total(self, sumW2=0, **kwargs):
        weights = self.weights(**kwargs)
        if sumW2: return np.sum(weights), np.sum(weights**2)**.5
        return np.sum(weights)
    def histogram(self,column_name,w_kwargs={}, **kwargs):
        pd_series = self[column_name]
        weights = self.weights(**w_kwargs)
        n,bins  = np.histogram(pd_series,weights=weights,**kwargs)
        stat_error = np.sqrt(np.histogram(pd_series,weights=weights**2, bins=bins)[0])
        bin_centers = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
        return n, stat_error, bins, bin_centers
    def boost(self, column_name, *meta_bin, b_kwargs={}, w_kwargs={}, **kwargs):
        import boost_histogram as bh
        pd_series = self[column_name]
        weights = self.weights(**w_kwargs)
        h = bh.Histogram(
            bh.axis.Regular(*meta_bin, metadata=column_name),
            storage=bh.storage.Weight(),
        )
        h.fill(pd_series,weight=weights, **b_kwargs)
        return h
    def __repr__(self):
        mean,std = self.DiLepMass.mean(),self.DiLepMass.std()
        return '''Name: {} \n shape: {}\n total: {:.2e}\n mass mean: {:.2f} std: {:.2f}'''.format(self.name,self.shape,self.total(),mean,std)
    def print(self):
        print(repr(self))

class SampleManagerPlotting(SampleManager):
    #@classmethod
    #def merge(cls,df_list): maybe a good idea, slow
    #    temp_cls = super().merge(df_list)
    #    print(type(temp_cls))
    def __init__(self, df, *args, name="", label="", category="", color="", **kwargs):
        super().__init__(df, *args, name=name, **kwargs)
        self.label = label
        self.category = category
        self.color = color
    @classmethod
    def from_file(cls,file_name,name,label,category,color, *args, **kwargs):
        df = pd.read_csv(file_name)
        return cls(df,name=name,label=label,category=category,color=color)
    def view(self, lt={}, gt={}, eq={}, **kwargs):
        cls_tmp = super().view(lt=lt, gt=gt, eq=eq, **kwargs)
        cls_tmp.label = self.label
        cls_tmp.category = self.category
        cls_tmp.color = self.color
        cls_tmp.name = self.name
        return cls_tmp

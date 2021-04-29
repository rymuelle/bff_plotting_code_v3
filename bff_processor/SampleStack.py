import numpy as np
from bff_processor.SampleManager import SampleManager, SampleManagerPlotting
import pprint

class SampleStack():
    def __init__(self, view={}, sample_list=[]):
        self._SMP = []
        for smp in sample_list:
            self.add_sample(*smp)
        self.view = {}
    def add_sample(self,*args,**kwargs):
        print(args)
        self._SMP.append(SampleManagerPlotting.from_file(*args, **kwargs))
    def SMP(self,view={}):
        view = {**self.view, **view}
        return [smp.view(**view, name=smp.name) for smp in self._SMP]
    def select_smp(self, label='',name='',category=''):
        eval_smp = lambda smp, l, n, c: (l == smp.label or not l) and (n == smp.name or not n) and (c == smp.category or not c)
        return [smp for smp in self.SMP() if eval_smp(smp, label, name, category)]
    def sum_hists(self,column_name,bins, label='',name='',category='',**kwargs):
        hists =  [smp.histogram(column_name, bins=bins,**kwargs) for smp in self.select_smp(label=label,name=name,category=category)]
        bins_c = np.array(hists)[:,3][0]
        bins = np.array(hists)[:,2][0]
        nom = np.sum(np.array(hists)[:,0],axis=0)
        std = np.sum( np.array(hists)[:,1]**2 ,axis=0)**.5
        return nom, std, bins, bins_c
    def __repr__(self):
        return pprint.pformat(self.SMP())

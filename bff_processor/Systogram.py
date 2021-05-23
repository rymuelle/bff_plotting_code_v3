import boost_histogram as bh
import numpy as np
from bff_processor.utils import sum_in_quad

class Systogram(bh.Histogram):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.args = args
        self.kwargs = kwargs
        self.systematics = np.zeros(0,dtype=[('name', object), ('up', object), ('down', object)])
    def add_sys(self, sys_name, verbose=False):
        if sys_name not in self.systematics['name']:
            hist_up = bh.Histogram(*self.args, **self.kwargs)
            hist_down = bh.Histogram(*self.args, **self.kwargs)
            self.systematics = np.append(self.systematics, np.array([(sys_name,hist_up,hist_down)], dtype=self.systematics.dtype))
            if verbose: print('Created {}!'.format(sys_name))
            return 1
        else:
            if verbose: print('{} already exists!'.format(sys_name))
            return 0
    def get_sys(self,sys_name, create=True, **kwargs):
        if create: self.add_sys(sys_name, **kwargs)
        sys_list = self.systematics[self.systematics['name']==sys_name]
        n_sys_list = len(sys_list)
        assert n_sys_list == 1, "Systematic '{}' should have length 1, has length {}".format(sys_name,n_sys_list)
        return sys_list[0]
    def nom_to_uncertainty(self):
        return hist2unc(self)
    def sys_to_array(self, relative=True):
        up = np.array(list(map(lambda x: x.values(), self.systematics['up'])))
        down = np.array(list(map(lambda x: x.values(), self.systematics['down']))) 
        if relative:
            up = np.subtract(up,self.values())
            down = np.subtract(down,self.values())
        up_down_sorted = map(lambda x: sorted(x, key=lambda y: -np.sum(y)), zip(up,down))
        up,down = list(zip(*up_down_sorted))
        return {'up': np.asarray(up), 'down': np.asarray(down)}
    def sys_sum(self, linear=True):
        sys_array = self.sys_to_array()
        if linear: return {'up': sum(sys_array['up']), 'down':sum(sys_array['down'])}
        return {'up': sum_in_quad(sys_array['up'], axis=0), 'down': -sum_in_quad(sys_array['down'], axis=0)}
    def nom_std_arrays(self):
        rs = self.sys_sum()
        nom = self.values()
        return nom, (self.variances()+rs['up']**2)**.5+nom, -(self.variances() +rs['down']**2)**.5+nom
     

from src.plotting_tools.colors import color_fader
from src.assets.file_groups import bck_dict, bck_list, bck_colors
import numpy as np
import pandas as pd
from src.general.utils import is_equal
from src.plotting_tools.SysHist import SysHist
from src.plotting_tools.Bins import Bins, bins
from src.plotting_tools.latexAssets import signal_type_dict
split_bins = bins

class StackPlotter():
    def __init__(self, plot_df, era, rebin=0, x_range=(-np.inf, np.inf)):
        self.plot_df = plot_df
        self.era = era
        self.rebin = rebin
        self.x_range = x_range
        self.scale=1
        
    def make_hist(self, row, preserve_binning=False):
        sh = SysHist.from_dict(row)
        sh *= self.scale
        if (not is_equal(self.rebin, 0) and not preserve_binning): sh = sh.rebin(self.rebin)
        if not preserve_binning: sh = sh.reduce_range(bottom=self.x_range[0], top=self.x_range[1])
        # ST isr fsr issue
        if row.sample_name in ['mc_santitop','mc_stop']:
            if 'Weight_ISRFSR_Up' in sh.sys: 
                sh.sys['Weight_ISRFSR_Up'] = [sh.sys['Weight_ISRFSR_Up'][0]*0, sh.sys['Weight_ISRFSR_Up'][0]*0]
        if sh.sys != {}: sh.sys_from_sys_dict()
        return sh 
    
    def feature_reg_df(self, feature, reg):
        return self.plot_df[(self.plot_df.reg==reg) & (self.plot_df.feature==feature) & (self.plot_df.era==int(self.era))]
    
    def stich_dy(self, feature, reg):
        '''stitched dy samples and replaces old seperate ones'''
        
        dy_sel = (self.plot_df.category.str.contains('DY')) & (self.plot_df.category != 'DYLL') & (self.plot_df.sample_name != 'dy_stitched')
        feat_reg = (self.plot_df.feature==feature) & (self.plot_df.reg==reg)
        ncategorys = len(self.plot_df[dy_sel].category.unique())
        
        dy_df = self.plot_df[dy_sel & feat_reg]
        comb_hist  = self.combine_hists(dy_df, preserve_binning=True).calc_ratio(ncategorys)
        
        #now label prestitched so it is not used in stack
        self.plot_df.loc[dy_sel & feat_reg,"type"] = 'bck_prestitched'
        #add new stitched row
        hist_dict = {"type": "bck", "category": "DY",  "sample_name": "dy_stitched", 'file':'dy_stitched',
                     "era": int(self.era), "reg": reg, "feature": feature,
                     **comb_hist.to_dict()}
        self.plot_df = self.plot_df.append(hist_dict, ignore_index=True)
        #avoid duplicates if run multiple times
        self.plot_df.drop_duplicates(subset=['type', 'category', 'file', 'mass', 'dbs', 'gmu', 'gb', 
                                             'era', 'xsec', 'sample_name','reg', 'feature'], inplace=True)
        return comb_hist
    
    def bck_df(self, feature, reg):
        tdf = self.feature_reg_df(feature, reg)
        return tdf[tdf.type=="bck"]
    def combine_hists(self, df, **kwargs):
        nhists = 0
        for i, row in df.iterrows():
            if not nhists:
                chist = self.make_hist(row, **kwargs)
            else: chist += self.make_hist(row, **kwargs)
            nhists += 1
        return chist
    def combine_back(self, feature, reg):
        self.stich_dy(feature, reg)
        bdf = self.bck_df(feature, reg)
        return self.combine_hists(bdf)     
    def draw_background(self, ax, feature, reg, ratio = -1, draw_sys=1,  error_scale=1, make_density=1, scale = 1, sys_label = None, nom_color='red', **kwargs):
        self.stich_dy(feature, reg)
        bdf = self.bck_df(feature, reg)
        hist_dict = {}
        nhists = 0
        _nominal_values = []
        bck_dict['DY'] = ['dy_stitched']
        for cat in bck_dict:
            _name_list = bck_dict[cat]
            _name_tdf = bdf[bdf.sample_name.isin(_name_list)]
            _chist = self.combine_hists(_name_tdf)
            if make_density: _chist = _chist.make_density_hist()
            _nominal_values.append(_chist.nominal) 
            
        labels = [signal_type_dict[x] for x in bck_dict]

        ax.stackplot(_chist.calc_bin_centers(), _nominal_values,
                     labels=labels, alpha=1, step='mid', colors=bck_colors)
        bhist = self.combine_back(feature, reg)
        if make_density: bhist = bhist.make_density_hist()
        bhist.draw(ax, alpha=.5, draw_sys=draw_sys, error_scale=error_scale, color=nom_color, sys_label=sys_label, label="MC background", **kwargs)
        return bhist
    def select_hists(self,**kwargs):
        tdf = self.plot_df
        for key, value in kwargs.items():
            tdf = tdf[tdf[key] == value]   
        tdf['hist'] = tdf.apply(self.make_hist, axis=1)
        return tdf
    def get_signal_hist(self, feature, reg, **kwargs):
        tdf = self.feature_reg_df(feature, reg)
        tdf = tdf[tdf.type=='sig']
        for key, value in kwargs.items():
            tdf = tdf[tdf[key] == value]
        shists = []
        for i, row in tdf.iterrows():
            shists.append({"hist":self.make_hist(row), **row.to_dict()})
        return pd.DataFrame(shists)

    def draw_signals(self, ax, feature, reg,
                     dbs_values = [0.04], mass_values = [125., 150., 175., 200., 250, 300, 350.],
                    c1='#ff2f00', c2='#0486ff',
                    ratio=-1,  draw_sys=1, make_density=1, **kwargs):
        nmass = len(mass_values)
        colors = [color_fader(c1,c2,mix=(i+.0)/nmass) for i in range(nmass)]
        
        tdf = self.feature_reg_df(feature, reg)
        sdf =  tdf[tdf.type=='sig']
        for dbs in dbs_values:
            for color, mass in zip(colors, mass_values): 
                _sdf = sdf[(sdf.mass==mass) & (sdf.dbs==dbs)]
                _shist = self.make_hist(_sdf.iloc[0])
                if make_density: _shist = _shist.make_density_hist()
                _shist.draw(ax, color=color, label='{} GeV'.format(int(mass)), draw_sys=draw_sys)
                
    def draw_signals_compare_dbs(self, ax, feature, reg,
                     dbs_values = [0.04, 1.0], mass_values = [125., 350.], 
                                 lss =['solid', 'dashed'],
                                c1='#ff2f00', c2='#0486ff',
                                ratio=-1,  draw_sys=1, make_density=1, **kwargs):
        nmass = len(mass_values)
        colors = [color_fader(c1,c2,mix=(i+.0)/(nmass-1)) for i in range(nmass)]
        
        tdf = self.feature_reg_df(feature, reg)
        sdf =  tdf[tdf.type=='sig']
        for ls, dbs in zip(lss,dbs_values):
            for color, mass in zip(colors, mass_values): 
                _sdf = sdf[(sdf.mass==mass) & (sdf.dbs==dbs)]
                _shist = self.make_hist(_sdf.iloc[0])
                if make_density: _shist = _shist.make_density_hist()
                _shist.normalize().draw(ax, color=color, label='{} GeV $\delta_{{bs}}$ = {}'.format(int(mass), dbs), draw_sys=draw_sys, ls=ls)
        return _shist
                
    def make_data_hist(self, feature, reg, blinded=True):
        mu_regions = ['SR1', 'SR2', 'CR10', 'CR20', 'CRA', 'CRD', 'CRA2', 'CRD2', 
                      'CR20_median', 'CRA_median', 'CRD_median', 'CRA2_median', 'CRD2_median']
        
        tdf = self.feature_reg_df(feature, reg)
        data_string = '_el.csv'
        if reg in mu_regions: data_string = '_mu.csv'
        ddf = tdf[(tdf.type=='data') & (tdf.file.str.contains(data_string))]
        _dhist = self.combine_hists(ddf)  
        if blinded and 'SR' in reg:
            _dhist*=0
        return _dhist
    
    def draw_data(self, ax, feature, reg, return_hist=0, make_density=1, **kwargs):
        _dhist = self.make_data_hist(feature, reg)
        if make_density: _dhist = _dhist.make_density_hist()
        ax.errorbar(_dhist.calc_bin_centers(), _dhist.nominal, yerr=_dhist.std, color='black',
                    label='Observed', ls='', marker='o')
        if return_hist: return _dhist
        return _dhist
    
    
def get_stack_plotter(output_dir, era, bins = 'split', x_range=[100,900]):
    plot_df = pd.read_pickle('{}/data/combined_{}_flat_hist.pkl'.format(output_dir, era))
    if bins=='split':
        bin_edges = split_bins.bin_edges
    else:
        bin_edges = 0
    return StackPlotter(plot_df, era, rebin=bin_edges, x_range=x_range)

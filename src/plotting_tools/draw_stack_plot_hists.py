import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.assets.file_groups import bck_dict, bck_list, bck_colors

from src.plotting_tools.Bins import Bins, bins, binning_type
from src.plotting_tools.cms_format import cms_format_fig, cms_style
from src.plotting_tools.SysHist import make_hist, make_sys_hist, era_corrected_sys_hist, SysHist
from src.plotting_tools.colors import color_fader

import mplhep as hep
hep.style.use(hep.style.CMS)
plt.rcParams.update({
    "text.usetex": True,
})


def draw_bckground(ax, tdf, feature, reg, era, ratio = -1, draw_sys=1,  error_scale=1, make_density=1, scale = 1, **kwargs):
    ##
    ## draw mc backgrounds
    ##
    #combined backgrounds
    _bhist = era_corrected_sys_hist(era, tdf, feature, reg, **kwargs)*scale
    # histograms of each category:
    _nominal_values = []
    _cats = []
    for cat in bck_dict:
        _name_list = bck_dict[cat]
        _name_tdf = tdf[tdf.name.isin(_name_list)]
        _thist = era_corrected_sys_hist(era,_name_tdf, feature, reg, **kwargs)*scale
        if ratio > 0:  _thist = _thist.calc_ratio(1./ratio)
        if make_density: _thist = _thist.make_density_hist()
        _nominal_values.append(_thist.nominal)
        _cats.append(cat)
        
    if ratio > 0: _bhist = _bhist.calc_ratio(1./ratio)
    ax.stackplot(_bhist.calc_bin_centers(), _nominal_values, labels=_cats, alpha=1, step='mid', colors=bck_colors)
    if make_density:
        _bhist.make_density_hist().draw(ax, color='gray', alpha=.5, draw_sys=draw_sys, error_scale=error_scale)
    else:
        _bhist.draw(ax, color='gray', alpha=.5, draw_sys=draw_sys, error_scale=error_scale)
    return _bhist
    
    
    
def draw_signals(ax, tdf, feature, reg, era,
                dbs_values = [0.04], mass_values = [125., 150., 175., 200., 250, 300, 350., 400, 450, 500., 750.],
                c1='#ff2f00', c2='#0486ff',
                ratio=-1,  draw_sys=1, make_density=1, **kwargs):
    nmass = len(mass_values)
    colors = [color_fader(c1,c2,mix=(i+.0)/nmass) for i in range(nmass)]
    
    for dbs in dbs_values:
        for color, mass in zip(colors, mass_values):
            _signal_tdf = tdf[(tdf.dbs==dbs) & (tdf.mass==mass)]
            _shist = era_corrected_sys_hist(era, _signal_tdf, feature, reg, **kwargs)
            if ratio > 0: _shist = _shist.calc_ratio(1./ratio)
            if make_density: _shist = _shist.make_density_hist()
            _shist.draw(ax, color=color, label='{} GeV'.format(int(mass)), draw_sys=draw_sys)
            
            
            
def draw_data(ax, tdf, feature, reg, era, return_hist=0, make_density=1, **kwargs):
    _dhist = era_corrected_sys_hist(era, tdf, feature, reg, **kwargs)
    if make_density: _dhist = _dhist.make_density_hist()
    ax.errorbar(_dhist.calc_bin_centers(), _dhist.nominal, yerr=_dhist.std, color='black', label='data',
           ls='', marker='o')
    if return_hist: return _dhist
    return _dhist.nominal.sum()



def draw_stackplot(ax, df, feature, reg, era, lumi,  draw_sig_sys=1, draw_bck_sys=1, error_scale=1, top=1e4, xlabel='$m_{\ell\ell}$ [GeV]',
                   **kwargs):
    background_df = df[df.type.str.contains('bck')]
    #removes DYJLL samples and so on as they duplicate ZtoLL samples
    #background_df = background_df[background_df.name.isin(bck_list)]
    #draw data
    if 'CR' in reg:
        total = draw_data(ax, df[df.type=='data'], feature, reg, era, **kwargs)
        _bhist = era_corrected_sys_hist(era, background_df, feature, reg, **kwargs)
        background_total = _bhist.nominal.sum()       
    #draws background
    draw_bckground(ax,  background_df, feature, reg, era, draw_sys=draw_bck_sys, error_scale=error_scale, **kwargs)
    #draw signals
    draw_signals(ax, df[df.type=='sig'], feature, reg, era, draw_sys=draw_sig_sys, **kwargs)
    
    ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Events per GeV')
    ax.legend(ncol=2)
    ax.set_ylim(1e-2, top)
    lumilabel_ = hep.cms.label(ax=ax, lumi=lumi,year=era);
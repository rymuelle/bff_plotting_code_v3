from src.assets.lumi import lumi_dict
from src.physics.ModelParams import ModelParams

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from src.plotting_tools.cms_format import cms_format_fig
from src.general.limit_utils import draw_other_experiments, draw_gmu_curve

##
## Create dfs and grids for gbs and dbses plotting
##
def make_gb_dbs_grid(reg, era, mass,
                     gbs, dbses, 
                     acceptance_model
                    ):
    lumi = lumi_dict[str(era)]     
    def figures_of_merit(_gb, _dbs): 
        mp = ModelParams(mass, _dbs, _gb, 0)
        acceptance = acceptance_model(mp, reg)
        return {"acceptance": acceptance,
         "xsec": mp.mumu_xsec_from_br(),
         #"xsec_acc": mp.mumu_xsec_from_br()*acceptance,
         "mumu_br": mp.mumu_br(),
         "nevents": acceptance*lumi*mp.mumu_xsec_from_br(),
         "nevents_inc": lumi*mp.mumu_xsec_from_br(),
         "inclusive_xsec": mp.inclusive_xsec(),
         "gmu": mp.gmu,
        }
    df = pd.DataFrame([{'gb':gb, 'dbs': dbs, **figures_of_merit(gb, dbs)} for gb in gbs for dbs in dbses])
    return df.sort_values(['gb', 'dbs'], ascending=True)

def make_fb_grid(reg, era, mass, gbs, dbses, acceptance_model, show = False, weights = 'xsec'):
    gb_dbs_grid = make_gb_dbs_grid(reg, era, mass, gbs=gbs.calc_bin_centers(), dbses=dbses.calc_bin_centers(), acceptance_model=acceptance_model)
    
    fig, ax = plt.subplots()
    counts,ybins,xbins,image = ax.hist2d(gb_dbs_grid.gb, gb_dbs_grid.dbs, weights=gb_dbs_grid[weights], 
                                     bins = (gbs.bin_edges, dbses.bin_edges))
    # hack to elimnate splitting of contours
    #counts[0,:]=1
    if not show: plt.close(fig)
    return counts, ybins, xbins, image, gb_dbs_grid
##
## Functions to plot our observed limits on db dbs space
##

def get_counts(limit_df, reg, mass, era):
    limit  = limit_df[(limit_df.nJets==reg) & 
            (limit_df.mass==mass) & 
            (limit_df.era==era)]
    lim_counts = limit[['n_16.0*acc','n_50.0*acc','n_84.0*acc']].mean()
    return lim_counts.to_numpy()


def connect_segments(contours):
    connected_segs = []
    for i in range(len(contours.allsegs)):
        connected_segs.append(np.concatenate(contours.allsegs[i]))
    return connected_segs

def make_contors(counts, x, y, limits, show=False):
    fig, ax = plt.subplots()
    contours = ax.contour(counts.T,
                             extent=[x.min(),x.max(),y.min(),y.max()],
                             linewidths=3, 
                             levels=limits)

    if not show: plt.close(fig)
    return connect_segments(contours)


def draw_three_level_limit(ax, contours, color='red', **kwargs):
    x0, y0 = contours[0].T
    x2, y2 = contours[2].T
    ax.fill(
        [*x0,*x2[::-1]],
        [*y0,*y2[::-1]],
        color=color,
        alpha=.1
    )
    ax.plot(*contours[1].T, color=color, **kwargs)
    
    
def draw_gb_dbs_limits(ax, reg, mass, era, gbs, dbses, limit_df, acceptance_model, **kwargs):
    lim_counts = get_counts(limit_df, reg, mass, era)
    count_landscape, ybins, xbins, image, gb_dbs_df = make_fb_grid(reg, era, mass, gbs, dbses, acceptance_model, show=False, weights='nevents')

    contours = make_contors(count_landscape, gbs.bin_edges, dbses.bin_edges, lim_counts, show=False)
    draw_three_level_limit(ax, contours, **kwargs)
    return contours

##
## Functions to plot inclusive limits on db dbs space
##
def draw_inclusive_gb_dbs_lim(ax, reg, mass, era, gbs, dbses, acceptance_model, predict_inc):
    fb_landscape, ybins, xbins, image, gb_dbs_df = make_fb_grid(reg, era, mass, gbs, dbses, acceptance_model, show=False, weights='xsec')
    if mass >=200:
        contors = ax.contour(fb_landscape.T,
                             extent=[gbs.bin_edges.min(),gbs.bin_edges.max(),
                         dbses.bin_edges.min(),dbses.bin_edges.max()],
                             linewidths=3, 
                             levels=[predict_inc(mass)],
                             linestyles=[':'],
                             colors=['black'], 
                             label='JHEP 07 (2021) 208'
                        )
##
## Functions to plot width broadening
##        
def draw_gmu_width_curve(ax, mass):
    _ = draw_gmu_curve(ax, mass, percent=10, color='black',linestyle=':', 
                                    zorder=1, label='width $<$ {}\%'.format(10))
    gmu_cons_curve = draw_gmu_curve(ax, mass, percent=1, color='black',linestyle='-.', 
                                    zorder=1, label='width $<$ {}\%'.format(1))
    return gmu_cons_curve
##
## intersection plots
##
from src.general.intersect import intersection

def find_intersection(x,y, one_lim_contour, gb=0.02, dbs=1.0):
    xc,yc = one_lim_contour.T
    x,y = intersection(x, y, xc, yc)
    filter_for_area = (x<gb) & (y<dbs)
    x = x[filter_for_area]
    y = y[filter_for_area]
    if len(x)==1: return x[0], y[0]
    return 0,0

def find_intersections(x,y, curves):
    xd, yd = find_intersection(x,y, curves[0])
    xn, yn = find_intersection(x,y, curves[1])
    xu, yu = find_intersection(x,y, curves[2])
    return np.array([[xd,yd], [xn,yn], [xu,yu]])

def plot_intersections(ax, reg, mass, _era, other_exp, width_curve, lim_contour, sys_dir):
        intesection_ty = find_intersections(other_exp['gb'], other_exp['ty'], lim_contour)
        intesection_gmu = find_intersections(width_curve['gb'], width_curve['dbs'], lim_contour)
        gb_limits = np.min([intesection_ty[:,0], intesection_gmu[:,0]], axis=0)
        dbs_limits = np.max([intesection_ty[:,1], intesection_gmu[:,1]], axis=0)
        
        ax.scatter(gb_limits, dbs_limits, marker="*", color='black')
        gb_limits_dict_ty = {'gb_'+k: v for k,v in zip(sys_dir,np.sort(intesection_ty[:,0]))}
        dbs_limits_dict_ty = {'dbs_'+k: v for k,v in zip(sys_dir,np.sort(intesection_ty[:,1]))}
        gb_limits_dict_gmu = {'gb_gmu_'+k: v for k,v in zip(sys_dir,np.sort(intesection_gmu[:,0]))}
        dbs_limits_dict_gmu = {'dbs_gmu_'+k: v for k,v in zip(sys_dir,np.sort(intesection_gmu[:,1]))}
        ty = (intesection_ty[:,1] > intesection_gmu[:,1]).all()
        return {"ty": ty, "era": _era, "reg": reg, "mass": mass, 
                                       **gb_limits_dict_ty, **dbs_limits_dict_ty, 
                                       **gb_limits_dict_gmu, **dbs_limits_dict_gmu}
    
    
##
## Functions to format gb dbs plot
##
def format_gb_dbs_plot(ax, mass):
    legend_opts = {
                "facecolor": 'white',
                "framealpha": 1,
                "frameon": True,
                "title":"$m_{{Z'}} = {}$ GeV".format(mass),
                "loc": 'upper right'
            }
    cms_format_fig("Run 2", ax, "\emph{Preliminary}")
    ax.legend(**legend_opts)
    ax.set_xlabel('$g_b$', fontsize=45)
    ax.set_ylabel('$\\delta_{bs}$', fontsize=45)
    ax.set_ylim([0,1])
    ax.set_xlim([0,.02])
    
def imshow(ax, counts, ybins, xbins):
    '''placeholder'''
    plt.imshow(count_landscape.T[::-1])
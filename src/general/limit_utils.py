import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from math import e
import numpy as np
from src.general.constants import c9

def compute_acceptance(x_df, reg):
    '''From region, df, use previous fit params to compute acceptance'''
    popt_acc = acceptance_dict[reg]
    return dbs_acceptance_fit(x_df,*popt_acc)

#compute expected xsec from mass


def branching_ratio_df(x):
    '''Computes branching ratio from xsec dataframe.'''
    gbgm_ratio = x.gb/x.gmu
    return 2/3*(1 + gbgm_ratio**2 * (1 + 2*x.dbs**2))**-1

def total_xsec_df(x):
    '''Computes xsec from xsec dataframe with branching ratio.'''
    return x.xsec/branching_ratio_df(x)

# formula to fit for c1 and k per mass point 
def mm_xsec_func(x, c1, k):
    '''Formula for xsec from gb, dbs, gmu'''
    gb = x.gb
    return c1 * gb**2 * (1 + k * x.dbs**2) * x.gmu**2
# formula to fit for c1 and k per mass point 
def total_xsec_func_df(x, c1, k):
    '''Formula for xsec from gb, dbs'''
    gb = x.gb
    return c1 * gb**2 * (1 + k * x.dbs**2) 


def power_law(x, c, *params):
    x = x.astype('float')
    y = x*0.
    for i, par in enumerate(params):
        y += np.power(x,i)*par
    return c*np.exp(y)

def linear(x, c1, c2):
    return c1*x + c2



def exponet_500(x,c1,c2,c3,c4,c5,c6,c7):
    '''Custom fit function to fit the c1 as a function of mass curve'''
    x = x/500
    return np.exp(c2*x+c3*x**-1+c4*x**-2+c5*x**-3)*(c1)


import pickle
with open('fits/limit_setting/popt_c1.pkl', 'rb') as f:
    popt_c1 = pickle.load(f)
with open('fits/limit_setting/popt_k.pkl', 'rb') as f:
    popt_k = pickle.load(f)

def compute_c1(x, popt_c1=popt_c1): return exponet_500(x, *popt_c1)
def compute_k(x, popt_k=popt_k): return linear(x, *popt_k)

def compute_expected_xsec(x, popt_c1, popt_k):
    '''compute xsec using c1 and k from params'''
    mass = x.mass
    c1, k  = compute_c1(mass, popt_c1), compute_k(mass, popt_k)
    return total_xsec_func_df(x, c1, k)


def produce_xsec_model(popt_c1, popt_k):
    '''load xsec model using c1 and k'''
    return lambda x: compute_expected_xsec(x, popt_c1, popt_k)
xsec_model = produce_xsec_model(popt_c1, popt_k)

def make_xsec_df_producion_values(mass, nSamples = 100):
    '''creates a dataframe of a varienty of gmu, gb, dbs, mass values'''
    dbs_gb_df = pd.DataFrame()
    dbs_gb_df['dbs'] = np.linspace(0,1, nSamples)
    dbs_gb_df['gb'] = np.full(nSamples, 0.02)
    dbs_gb_df['gmu'] = np.full(nSamples, 0.17)
    dbs_gb_df['mass'] = np.full(nSamples, mass)
    return dbs_gb_df



def compute_xsec(df):
    from src.assets.lumi import lumi_dict
    df['comp_total_xsec'] = xsec_model(df)
    df['branching_ratio'] = branching_ratio(df)
    df['comp_xsec'] = df['branching_ratio']*df['comp_total_xsec']
    for reg in ['SR1', 'SR2', 'SRX']:
        df['comp_acc_'+reg] = compute_acceptance(df, reg)
        df['comp_xsec_'+reg] = df['comp_xsec']*df['comp_acc_'+reg]
        for era in [2016, 2017, 2018]:
            lumi = lumi_dict[str(era)]
            df['comp_nevt_{}_{}'.format(reg,era)] = df['comp_xsec_'+reg]*lumi
    return df

def make_mass_x_df(dbs, nSamples = 100):
    dbs_gb_df = pd.DataFrame()
    dbs_gb_df['mass'] = np.linspace(125,750, nSamples)
    dbs_gb_df['dbs'] = np.full(nSamples, dbs)
    return dbs_gb_df

def make_dbs_mass_df():
    mass = np.array([125,150,175,200,350,500])
    gmu = mass_to_gmu(mass)
    #gb = np.linspace(2e-4, 2e-2, int((2e-2-2e-4)/2e-4+1))
    dbs = np.linspace(1e-2, 1., int((1-1e-2)/(1e-2)+1))
    df =  pd.DataFrame([{"mass": m, "gmu": mass_to_gmu(m), "dbs": d} for m in mass for d in dbs])
    df['gb'] = df.apply(calc_gb, axis=1)
    return compute_xsec(df)
def make_gb_dbs_mass_df():
    mass = np.array([125,150,175,200,350,500])
    gmu = mass_to_gmu(mass)
    gb = np.linspace(2e-4, 2e-2, int((2e-2-2e-4)/2e-4+1))
    dbs = np.linspace(1e-2, 1., int((1-1e-2)/(1e-2)+1))
    df = pd.DataFrame([{"mass": m, "gmu": mass_to_gmu(m),"gb":b, "dbs": d} for m in mass for b in gb for d in dbs])
    return compute_xsec(df)

def make_mass_df(mass):
    gmu = mass_to_gmu(mass)
    #gb = np.linspace(2e-4, 2e-2, int((2e-2-2e-4)/2e-4+1))
    dbs = np.linspace(1e-4, 1., int((1-1e-2)/(1e-2)+1))
    df =  pd.DataFrame([{"mass": mass, "gmu": mass_to_gmu(mass), "dbs": d} for d in dbs])
    return df

def calc_gb(row):
    from src.general.constants import c9
    return c9*(row.mass/100)**2/(row.gmu*row.dbs)

def calc_gmu(row):
    from src.general.constants import c9
    return c9*(row.mass/100)**2/(row.gb*row.dbs)

def calc_gmu(mass, gb, dbs):
    from src.general.constants import c9
    return c9*(mass/100)**2/(gb*dbs)

def branching_ratio(gb, dbs, gmu):
    '''Computes branching ratio from xsec dataframe.'''
    gbgm_ratio = gb/gmu
    return 2/3*(1 + gbgm_ratio**2 * (1 + 2*dbs**2))**-1

def total_xsec_func(gb, dbs, c1, k):
    '''Formula for xsec from gb, dbs'''
    return c1 * gb**2 * (1 + k * dbs**2) 

# not needed really, but making a linear fit to mass and gmu values from paper values to show validity of min gmu value for mass
gmu = [0.08, 0.14, 0.20]
mass = np.array([200,350,500])
def linear(x,c1,c2): 
    return x*c1+c2
popt_mass_gmu ,_ = curve_fit(linear, mass, gmu)
def mass_to_gmu(x): return linear(x, *popt_mass_gmu)


#used for acceptance fit
def logisitic(x, c1, c2):
    return 1/(1+10**(-(x-c1)/c2))

def dbs_func(x, c1, c2):
    return c1*(1+c2*x**.5)

def dbs_acceptance_fit(x, c1, c2, c3, c4, c5, *params):
    #y = poln(x.mass, *params)
    y = logisitic(x.mass, *params)
    yp = c3 + c4*y + c5*y**2
    return yp*dbs_func(x.dbs, c1, c2)
##
#gb dbs 
##
def draw_other_experiments(ax, mass, gb = np.linspace(1e-8, 0.25,1000)):
    from src.physics.gb_dbs_constraings import  BS_BS_y, trident_y
    #this draws exclusion curves from other experiments
    #draw BS-BS, neutrino triden
    
    bsbsy = BS_BS_y(mass, gb)
    ty = trident_y(mass, gb)
    #ax.fill_between(gb,bsbsy,bsbsy+999, color='#c2c2c2', label=r'$B_s-\bar{B_s}$', zorder=0)
    ax.fill_between(gb,bsbsy,bsbsy+999, color='#F3E5AB', label=r'$B_s-\bar{B_s}$', zorder=0)
    ax.fill_between(gb,ty, color='#ffd7ff', label='$\\nu$ Trident', zorder=0)
    return {"gb": gb, "bsbsy":bsbsy, "ty":ty}
def draw_gmu_curve(ax, mass, percent=1, **kwargs):
    #draw gmu width curve
    gb, dbs = curve_of_const_gmu(mass, width_to_gmu(percent))
    ax.plot(gb, dbs, **kwargs)
    return {"gb": gb, "dbs":dbs}
def format_gb_gmu_plot(ax, mass=""):
    ax.set_ylim(0,1.01)
    ax.set_xlim(0,0.02)
    legend_opts = {
            "facecolor": 'white',
            "framealpha": 1,
            "frameon": True
        }
    if mass!="": legend_opts["title"] = "$m_{{Z'}} = {}$ GeV".format(mass)
    legend = ax.legend( **legend_opts)
    ax.set_xlabel('$g_b$')
    ax.set_ylabel('$\\delta_{bs}$')
    cms_format_fig(str(lera), ax,  data=True, label="Work in Progress")
#sigma[%] = exp[2ln(g_mu) + 1.381] + 0.00038

def width_to_gmu(width):
    return np.exp((np.log(width-0.00038)-1.381)/2)
def curve_of_const_gmu(zmass, gmu):
    dbs = np.linspace(1e-4, 1., int((1-1e-4)/(1e-2)+1)) 
    gb = c9*(zmass/100)**2/(gmu*dbs)
    return gb, dbs
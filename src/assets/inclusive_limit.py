import pandas as pd
import matplotlib.pyplot as plt
from src.general.functions import make_bpoly, power_func
from scipy.interpolate import interp1d
import numpy as np

# from figure 4 https://arxiv.org/pdf/2103.02708.pdf
# raw data: https://www.hepdata.net/record/ins1849964

hepdata_limit_file = 'limits/HEPData-ins1849964-v2-Observed_limit_in_narrow_width,_spin_1.csv'
inclusive_data = pd.read_csv(hepdata_limit_file)
inclusive_data['xsec'] = inclusive_data['xsec(pb)']*1000
inclusive_data['xsec_log'] = np.log(inclusive_data['xsec'])

def make_observed_limit_inter():
    plt.scatter(inclusive_data.mass, inclusive_data['xsec'])

    f = interp1d(inclusive_data.mass, inclusive_data['xsec'])

    def predict_inc(mass):
        return f(mass)

    x= np.linspace(200,5000, 1000)
    y = predict_inc(x)
    plt.plot(x,y, color='orange')
    plt.yscale('log')
    plt.ylabel('xsec [fb]')
    plt.xlabel('DiLepMass [GeV]')
    return predict_inc
    
predict_inc = make_observed_limit_inter()


def make_mass_df(mass):
    gmu = mass_to_gmu(mass)
    #gb = np.linspace(2e-4, 2e-2, int((2e-2-2e-4)/2e-4+1))
    dbs = np.linspace(1e-4, 1., int((1-1e-4)/(1e-4)+1))
    df =  pd.DataFrame([{"mass": mass, "gmu": mass_to_gmu(mass), "dbs": d} for d in dbs])
    return df

from src.general.functions import gb_from_xsec_gmu_dbs_br


def draw_inclusive_dbs_curve(ax, _df, **kwargs):
    excluded_fb = predict_inc(_df.mass)
    gb = gb_from_xsec_gmu_dbs_br(excluded_fb, _df)
    ax.plot(gb,_df.dbs, **kwargs)
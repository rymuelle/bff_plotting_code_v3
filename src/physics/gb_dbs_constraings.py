import math

vtb = 1.021
vts = 0.04
wma = .8768 # weak mixing angle cos(theta_w)
from src.general.constants import c9


#def trident(mass): return 540/mass * 1.306 * 10 ** -9 * 200 ** 2
#corrected version?
def trident(mass): return 540/mass * c9 * 10**-4 * 200 ** 2

def BS_BS(mass,vtb,vts,wma): 
    greater_than = 0.16
    greater_than = greater_than * (1 - .029 * math.log(mass / 1000)) ** -1
    greater_than = greater_than * (mass / 1000) ** 2
    greater_than = greater_than ** .5
    greater_than = greater_than * (vtb*vts) / (2 * wma) 
    return greater_than

def xs_ys(x,c):
    return c*(1/x)


def trident_y(mass, x):
    return xs_ys(x,trident(mass))

def BS_BS_y(mass, x):
    return xs_ys(x,BS_BS(mass, vtb,vts,wma))

def get_plot_values(mass):
    trident_dict = {350:trident(350),
        200:trident(200),
        500:trident(500)}
    bs_dict  = {350:BS_BS(350,vtb,vts,wma),
        200:BS_BS(200,vtb,vts,wma),
        500:BS_BS(500,vtb,vts,wma)}
    inc_xsection = {350:9,
        200:53,
        500:2.45}
    return_obj = [bs_dict,trident_dict,inc_xsection]
    return [x[mass] for x in return_obj]
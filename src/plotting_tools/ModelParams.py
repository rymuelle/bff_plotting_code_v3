from src.general.constants import c9
from src.general.limit_utils import branching_ratio, total_xsec_func, compute_c1, compute_k

class ModelParams:
    def __init__(self, mass, dbs, gb, gmu):
        self.mass = mass
        self.dbs = dbs
        self.gb = gb
        if gmu:
            self.gmu = gmu
        else:
            self.gmu = c9*(mass/100)**2/(gb*dbs)
        self.c1 = compute_c1(mass)
        self.k = compute_k(mass)
    def inclusive_xsec(self):
        return total_xsec_func(self.gb, self.dbs, self.c1, self.k)
    def mumu_xsec_from_br(self):
        return self.inclusive_xsec()*self.mumu_br()
    def is_valid(self):
        calc_c9 = self.gmu*self.dbs*self.gb*(100/self.mass)**2
        percent_off = (calc_c9 - c9)/c9
        return percent_off < 1e-5
    
    def mumu_br(self):
        return branching_ratio(self.gb, self.dbs, self.gmu)
    def _mumu_xsec(self):
        return 2*self.inclusive_xsec()*self.gmu**2
    def _nunu_xsec(self):
        return self.inclusive_xsec()*self.gmu**2
    def _bb_xsec(self):
        return 3*self.inclusive_xsec()*self.gb**2
    def _bs_xsec(self):
        return 6*self.inclusive_xsec()*self.gb**2*self.dbs**2
    
    def __repr__(self):
        text = '''mass: {} dbs: {:.2e} gb: {:.2e} gmu: {:.2e} c1: {:.2e} k : {:.2f}
is valid: {} 
inc. xsec: {:.2e} mumu br: {:.2f}
mumu xsec: {:.2e}'''.format(
        self.mass, self.dbs, self.gb, self.gmu, self.c1, self.k,
        self.is_valid(), self.inclusive_xsec(), self.mumu_br(),
        self.mumu_xsec_from_br())
        return text    
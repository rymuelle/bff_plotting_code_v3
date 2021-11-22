# Making skims:
bff_plotter.ipynb

# calculating sig cuts:
creates bff_dict of cut values

optimize_2d_cuts.ipynb

# produce combined df
This should be run for each era. It produces a single df of all eras for ease of next steps (and it applies cuts)

produce_combined_df.ipynb

# Bins 

plotting_meta includes the shared bin value function.

```
from plotting_meta.plotting_meta import bin_edges, bin_centers, bin_widhts, x_range
```

# Make histogram with systematics

```
from bff_plotting_tools.make_hists import make_sys
hist = make_sys(df[df.mass==125], 'DiLepMass', "SR1")

fig, ax = plt.subplots(figsize=(20,10))
ax.scatter(hist.calc_bin_center(), hist.nominal)
```


example drawing a ratio hist:

```
hist = make_sys(df[df.mass==125], 'DiLepMass', "SR1")
dhist = hist.make_density_hist()
rhist = hist.calc_ratio(hist.nominal)

fig, ax = plt.subplots(figsize=(15,10))
dhist.draw(ax, color='red', label='350')
ax.set_yscale('log')
ax.set_ylabel('Events per GeV')
ax.legend()
```

# Make ABCD pred

fit_and_abcd.ipynb

# Make a fit to signal

fit_signal_restarting.ipynb

```
import zfit
from zfit import z
from bff_signal_model.bff_signal_model import bff_signal_model, reset_params, mu
...
def make_plot(mass):
    #mc hist:
    hist = make_sys(df[(df.mass==mass) & (df.dbs==0.04)], 'DiLepMass', "SR1")
    
    #make fit hist
    obs = zfit.Space("x", limits=bins.calc_bin_range())
    sm = bff_signal_model(obs, mu=mass)
    fit_hist = sm.make_hist(bins, np.sum(hist.nominal), tail_percent=1)
    
    fig, ax = plt.subplots(figsize=(15,10))
    hist.make_density_hist().draw(ax, color='red', label='MC {}'.format(mass))
    fit_hist.make_density_hist().draw(ax, color='blue', label='parametric {}'.format(mass))
    ax.set_yscale('log')
    ax.set_ylabel('Events per GeV')
    ax.set_ylim(bottom=1e-3)
    ax.legend()
```
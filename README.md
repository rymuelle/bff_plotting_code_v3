# Setting up on lxplus:

I set up in CMSSW_12_1_0

```
cmsenv
git checkout <repo path>
# create virtual enviroment for installing pacakges
python3 -m venv env
# activate env
. env/bin/activate
#set it up so you can use the kernel in jupyter
pip install ipykernel
python -m ipykernel install --user --name=bff_12_1
pip install <insert packages here>
```

Now, you should be able to log in to lxplus with port forwarding:
```
 ssh -L 8080:localhost:8080 rymuelle@lxplus.cern.ch
 cd <path to git repo>
 jupyter-notebook --no-browser --port=8080
```

You should now be able to open up the jupyter notebook launch page from the links printed out on a local browser.

# Making skims:
Takes yml files and produces csv skims. Tested in 12_1_0.

0_bff_skimmer.ipynb

# calculating sig cuts:
Creates bff_dict of cut values. Needed to produced .dill file of cut values used in the next step. Tested in 12_1_0.

1_optimize_2d_cuts.ipynb

# produce combined df
This should be run for each era. It produces a single df of all eras for ease of next steps (and it applies cuts)
Tested in 12_1_0.

2_produce_combined_df.ipynb

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
Makes abcd fits using lognorm function. 
Tested in 12_1_0

3_fit_and_abcd.ipynb

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
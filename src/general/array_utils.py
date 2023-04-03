import numpy as np

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def moving_sum(a, n=3, module=np) :
    ret = module.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:]

def moving_sum(a, n=3, module=np):
    return module.array([a[i:i+n].sum() for i in range(len(a)-n+1)])
    


def super_sample(func, bin_edges, n_samples=3e4):
    x_super_sample = np.linspace(bin_edges[0], bin_edges[-1], int(n_samples))
    y_super_sample = func(x_super_sample)
    counts, _ = np.histogram(x_super_sample, weights=y_super_sample, bins=hist.bin_edges)
    y = counts/n_samples
    return y

def super_sample_function(function, bin_edges, n_samples=3e4):
    return lambda x, *popt: super_sample(lambda x: function(x,*popt), bin_edges)

def moving_avg_func(n, function):
    return lambda *args: moving_average(function(*args), n=n)

def unp_array_to_nom_std(arr):
    nom = np.array(list(map(lambda x: x.nominal_value, arr)))
    std = np.array(list(map(lambda x: x.std_dev, arr)))
    return nom, std

def unzip(arr):
    return list(zip(*arr))
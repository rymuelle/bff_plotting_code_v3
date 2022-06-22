import pyarrow.feather as feather
import pandas as pd

def get_data(era):
    if era=='2016':
        lumi=35.50
    if era=='2017':
        lumi=41.85
    if era=='2018':
        lumi=58.88
    if era=='16-18':
        lumi = 35.5+41.85+58.88
        df = feather.read_feather('data/combined_2016.feather'.format(era))
        print(df.shape)
        df = df.append( feather.read_feather('data/combined_2017.feather'.format(era)))
        print(df.shape)
        df = df.append( feather.read_feather('data/combined_2018.feather'.format(era)))
        print(df.shape)
        print("loaded all df")
    else:
        df = feather.read_feather('data/combined_{}.feather'.format(era))
    df.replace([np.inf, -np.inf], 0, inplace=True)
    return df, lumi

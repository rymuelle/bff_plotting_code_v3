import yaml
import pandas as pd

samplefname = {}
samplefname[2016] = "samples/samplesCR_2016_Apr2020.yml"
samplefname[2017] = "samples/samplesCR_2017_Apr2020.yml"
samplefname[2018] = "samples/samplesCR_2018_Apr2020.yml"

def make_sample_df():
    sample_list = []  
    for era, file_name in samplefname.items():
        with open(file_name,'r') as f:
            sample_dict = yaml.load(f, Loader=yaml.FullLoader)
            lumi = sample_dict['lumi']
            sample_list += [{**x, 'era':era, 'lumi': lumi} for x in sample_dict['samples']]
    sample_df =  pd.DataFrame(sample_list)
    sample_df['file'] = sample_df.apply(lambda x: "data/tw_{}_{}.csv".format(x.era, x['name']), axis=1)
    return sample_df
from __future__ import print_function
import argparse
import glob
import re
import subprocess
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('era', type=int)

args = parser.parse_args()

print("processing", args.era)

file_list = [f for f in glob.glob('combine_cfgs/*') if str(args.era) in f]

combine_command_template = "combine -M AsymptoticLimits {} --run blind"

limit_list = []
for i, f in enumerate(file_list):
    command = combine_command_template.format(f)
    print('file: {}, {} out of {}'.format(f, i, len(file_list)))
    p = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
    out, err = p.communicate()
    
    nJets = re.findall(r'SR(\d)', f)[0]
    mass = re.findall(r'M_(\d+)_', f)[0]
    dbs = re.findall(r'dbs(\d)p(\d+)', f)[0]
    dbs = '{}.{}'.format(*dbs)
    limits = re.findall(r'Expected +(\d+.\d+)%: r < (\d+.\d+)',out)
    lim_dict = {k:float(v) for k,v in limits}
    lim_dict['mass'] = int(mass)
    lim_dict['nJets'] = int(nJets)
    lim_dict['dbs'] = float(dbs)
    limit_list.append(lim_dict)

df = pd.DataFrame(limit_list)

df.to_csv('lim_{}.csv'.format(args.era))

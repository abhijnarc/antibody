#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Created on Tue Mar 5 15:46:51 2024

script to cript to get top-10 models from capri_ss file
top 10 models in terms of irmsd (<= 2) and haddock score ( < 0)

Then get the binding energy by prodigy

To be run from prodigy environment

@author: user """

import subprocess
#import glob
import sys
import gzip
import shutil
import pandas as pd
import os

# get prefix from argument
pref = sys.argv[1]
#pref = 'V2F9'

# capri file
capri_file = pref + '/12_caprieval/capri_ss.tsv'
# output file
out_file = 'energy/energy_' + pref +  '.csv'
# directory
clust_dir = pref +'/11_seletopclusts/'
# distance cutoff value
dcut = 4.5

# Initialize a list
pdres = []

# read the capri file
capri_df = pd.read_csv(capri_file, sep='\s+')
#capri_df = pd.read_csv(capri_file, delim_whitespace = True)

# replace the model name with split value
capri_df.model = capri_df['model'].str.split('/').str[2]

# filter with irmsd <= 2 and haddock score < 0

capri_filt = capri_df[(capri_df['irmsd'] <= 2) & (capri_df['score'] < 0)]
capri_filt2 = capri_filt[['model', 'score', 'irmsd', 'fnat', 'lrmsd','dockq']]

# sort the dataframe on irmsd
capri_sort = capri_filt2.sort_values(by = ['irmsd'])

# get the top 10
capri_top = capri_filt2.head(n = 10)

# run a loop over the models in dataframe
for ind in capri_top.index:
    print('starting ' + capri_top['model'][ind])
    # pdb file with full path
    pdb = capri_top['model'][ind] + '.gz'
    uzpdb = clust_dir + capri_top['model'][ind]
    gzpdb = clust_dir + pdb
    # get model name
    temp = capri_top['model'][ind].split('.')[0]
    model = temp.split('_')[3]
    print('analyzing model ' + model)

    # unzip the filename - only if it exists
    if os.path.isfile(gzpdb):
        with gzip.open(gzpdb, 'rb') as f_in:
            with open(uzpdb, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    # run prodigy as a system command
    args = ('/home/user/miniconda3/envs/prodigy/bin/prodigy ' + uzpdb +
            ' --selection A B --contact_list --distance-cutoff ' + str(dcut))
    pres = subprocess.run(args, capture_output = True, shell = True,
                         encoding = "utf-8").stdout
    perr = subprocess.run(args, capture_output = True, shell = True,
                          encoding = "utf-8").stderr
    print(perr)

    # error trapping
    if len(perr) > 0:
        nintf = 0
        energy = kd = 9999
    else:
        p_array = pres.splitlines()
        kd = float(p_array[-1].split(':')[1])
        energy = float(p_array[-2].split(':')[1])
        nintf = float(p_array[-11].split(':')[1])
 
    pdres.append({'model': temp,
                  'haddock': capri_top['score'][ind],
                  'irmsd': capri_top['irmsd'][ind],
                  'fnat': capri_top['fnat'][ind],
                  'lrmsd': capri_top['lrmsd'][ind],
                  'dockq': capri_top['dockq'][ind],
                  'n_intf': nintf, 'Kd': kd, 'energy': energy})
    
pdres_df = pd.DataFrame(pdres)

    # sort the dataframe
    #pdres_sort = pdres_df.sort_values(by = ['energy'])

# save to a csv file
pdres_df.to_csv(out_file, encoding = 'utf-8', index = False)

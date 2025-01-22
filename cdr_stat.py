# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 17:12:18 2024

Script to  get CDRs and their lengths

@author: sujan
"""

import pandas as pd
import glob
import re

# define chain type
chain = 'heavy'

# output filename
out_file = chain + '_CDRstat.csv'

# cluster filename
clus_file = './clus_pseudo/' + chain + 'Clus_cluster.tsv'

# read the cluster file
clus_df = pd.read_csv(clus_file, sep = '\t', header = None)

# assign names
clus_df.columns = ['clid', 'contig']

# get the unique clid's
clus_list = clus_df.clid.unique()

# initiatlize a list
cdr_stat = []

# run a loop 
for cl in clus_list:
    airr_file = 'airr/' + chain + '_' + cl + '.airr'
    print('Processing ... ' + cl + '\n')
    airr_data = pd.read_csv(airr_file, sep = '\t')

    # get each of the CDR - skip the record if its empty
    # CDR1
    if airr_data['cdr1_aa'].notna().sum() > 0 :
        cdr1 = airr_data['cdr1_aa'].to_list()[0]
        cdr1 = cdr1.replace('*', '')
    else:
        continue
    # CDR2
    if airr_data['cdr2_aa'].notna().sum() > 0 :
        cdr2 = airr_data['cdr2_aa'].to_list()[0]
        cdr2 = cdr2.replace('*', '')
    else:
        continue
    # CDR3
    if airr_data['cdr3_aa'].notna().sum() > 0 :
        cdr3 = airr_data['cdr3_aa'].to_list()[0]
        cdr3 = cdr3.replace('*', '')
    else:
        continue

 
    #append to the list
    cdr_stat.append({'contig': cl,
                    'chain': chain,
                    'CDR1': cdr1,
                    'CDR1_len': len(cdr1),
                    'CDR2': cdr2,
                    'CDR2_len': len(cdr2),
                    'CDR3': cdr3,
                    'CDR3_len': len(cdr3)}
                    )

# convert the list to dataframe
cdr_df = pd.DataFrame(cdr_stat)

# sort the dataframe
cdr_sort = cdr_df.sort_values(by = ['CDR3_len'], ascending = False)

# print to outfile
cdr_sort.to_csv(out_file, encoding = 'utf-8', index = False)


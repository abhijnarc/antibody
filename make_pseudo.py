# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 17:12:18 2024

Script to parse the aiir files and make pseudo sequences

@author: sujan
"""

import pandas as pd
import glob
import re

# define chain type
chain = 'light'

# output fasta filename
out_file = chain + '_pseudo.fa'
of = open(out_file, 'w')

# run through the airr files
for filename in glob.glob('/data/Abhijna/RS0566/LCN2/LCN2/LCN2_airr.tsv'):
    temp = filename.split('_')[1]
    contig = temp.split('.')[0]

    print('Processing ... ' + contig + '\n')
    airr_data = pd.read_csv(filename, sep = '\t')

    # get each of the CDR - skip the record if its empty
    # CDR1
    if airr_data['cdr1'].notna().sum() > 0 :
        cdr1 = airr_data['cdr1'].to_list()[0]
    else:
        continue
    # CDR2
    if airr_data['cdr2'].notna().sum() > 0 :
        cdr2 = airr_data['cdr2'].to_list()[0]
    else:
        continue
    # CDR3
    if airr_data['cdr3'].notna().sum() > 0 :
        cdr3 = airr_data['cdr3'].to_list()[0]
    else:
        continue

    pseudo = cdr1 + cdr2 + cdr3
    of.write('>' + contig + '\n')
    of.write(pseudo + '\n')

of.close()

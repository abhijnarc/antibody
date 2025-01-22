# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 17:12:18 2024

Script to parse the aiir files and make Fv  sequences

@author: sujan
"""

import pandas as pd
import glob
import re

# define chain type
chain = 'heavy'

# output fasta filename
out_file = chain + '_fv.fa'
of = open(out_file, 'w')

# run through the airr files
for filename in glob.glob('airr/' + chain + '_*.airr'):
    temp = filename.split('_')[1]
    contig = temp.split('.')[0]

    print('Processing ... ' + contig + '\n')
    airr_data = pd.read_csv(filename, sep = '\t')

    # get each of the CDR - skip the record if its empty
    #FWR1
    if airr_data['fwr1_aa'].notna().sum() > 0 :
        fwr1 = airr_data['fwr1_aa'].to_list()[0]
    else:
        continue
    # CDR1
    if airr_data['cdr1_aa'].notna().sum() > 0 :
        cdr1 = airr_data['cdr1_aa'].to_list()[0]
    else:
        continue
    #FWR2
    if airr_data['fwr2_aa'].notna().sum() > 0 :
        fwr2 = airr_data['fwr2_aa'].to_list()[0]
    else:
        continue
    # CDR2
    if airr_data['cdr2_aa'].notna().sum() > 0 :
        cdr2 = airr_data['cdr2_aa'].to_list()[0]
    else:
        continue
    #FWR3
    if airr_data['fwr3_aa'].notna().sum() > 0 :
        fwr3 = airr_data['fwr3_aa'].to_list()[0]
    else:
        continue
    # CDR3
    if airr_data['cdr3_aa'].notna().sum() > 0 :
        cdr3 = airr_data['cdr3_aa'].to_list()[0]
    else:
        continue
    #FWR4
    if airr_data['fwr4_aa'].notna().sum() > 0 :
        fwr4 = airr_data['fwr4_aa'].to_list()[0]
    else:
        continue

    fv = fwr1 + cdr1 + fwr2 + cdr2 + fwr3 + cdr3 + fwr4
    # replace the *
    fv2 = fv.replace('*', '')

    of.write('>' + contig + '\n')
    of.write(fv2 + '\n')

of.close()

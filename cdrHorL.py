#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# # code to clean interface data obtained from prodigy
# for a full VH-VL antibody from all ic files in a directory

#from functools import reduce
import pandas as pd
import glob
import os
#import glob
import sys

# One to three letter AA dict
aa_123 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
          'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
          'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
          'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
          'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}


# Invert the dict to make it three to one letter
aa_321 = {v: k for k, v in aa_123.items()}


# Prefix for the sequence
pref = sys.argv[1]


# Function to get CDR mapping
def cdr_map(row):
    cdr_file = '/data/Abhijna/scripts/cdr_coord.csv'
    cdr_df = pd.read_csv(cdr_file)
    cdr_df.set_index("seq", drop=True, inplace=True)
    cdr_dict = cdr_df.to_dict(orient="index")


    # Extract the CDR coordinates
    h1s, h1e = cdr_dict[pref]['hcdr1_start'], cdr_dict[pref]['hcdr1_end']
    h2s, h2e = cdr_dict[pref]['hcdr2_start'], cdr_dict[pref]['hcdr2_end']
    h3s, h3e = cdr_dict[pref]['hcdr3_start'], cdr_dict[pref]['hcdr3_end']
    hlen = cdr_dict[pref]['hc_len']
    l1s, l1e = cdr_dict[pref]['lcdr1_start'] + hlen, cdr_dict[pref]['lcdr1_end'] + hlen
    l2s, l2e = cdr_dict[pref]['lcdr2_start'] + hlen, cdr_dict[pref]['lcdr2_end'] + hlen
    l3s, l3e = cdr_dict[pref]['lcdr3_start'] + hlen, cdr_dict[pref]['lcdr3_end'] + hlen


    # Check which domain the position belongs to
    if h1s <= row['ab_pos'] <= h1e or h2s <= row['ab_pos'] <= h2e or h3s <= row['ab_pos'] <= h3e:
        return 'H'
    elif l1s <= row['ab_pos'] <= l1e or l2s <= row['ab_pos'] <= l2e or l3s <= row['ab_pos'] <= l3e:
        return 'L'
    else:
        return 'FWR'


# Function to classify ag_pos
def classify_ag_pos(df):
    grouped = df.groupby('ag_pos')['ab_dom'].apply(set).reset_index()
    grouped['binding'] = grouped['ab_dom'].apply(lambda x: 'both' if {'H', 'L'}.issubset(x) else ('H' if 'H' in x else 'L'))
    return grouped[['ag_pos', 'binding']]


# get the prefix and model number and set file
for fname in glob.glob(pref + '/11_seletopclusts/cluster_1_model_1.ic'):
    bdir, ic_file = os.path.split(fname)
    model = ic_file.split('.')[0]
    seqdir = bdir.split('/')[0]
    print("processing {}, model {} ....".format(bdir, model))
    out_file = os.path.join(bdir, model + '_intf.csv')


    # Read the IC file
    ic_df = pd.read_csv(fname, sep=r'\s+', header=None)
    ic_df.columns = ['ab_aa3', 'ab_pos', 'ab_chain', 'ag_aa3', 'ag_pos', 'ag_chain']


    # Replace three-letter code with one-letter code
    ic_df['ab_aa'] = ic_df['ab_aa3'].map(aa_321)
    ic_df['ag_aa'] = ic_df['ag_aa3'].map(aa_321)


    # Reduce dataframe to key columns
    red_df = ic_df[['ab_pos', 'ab_aa', 'ag_pos', 'ag_aa']].drop_duplicates().sort_values(by=['ag_pos'])


    # Map binding domain (H or L)
    red_df['ab_dom'] = red_df.apply(cdr_map, axis=1)


    # Classify ag_pos interactions
    classified_df = classify_ag_pos(red_df)


    # Merge the classification back with the original dataframe
    final_df = pd.merge(red_df, classified_df, on='ag_pos', how='left')


    # Save the output
    final_df.to_csv(out_file, index=False)
    print(f"Output saved to {out_file}")

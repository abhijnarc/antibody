# code to clean interface data obtained from prodigy
# for a full VH-VL antibody from all ic files in a directory

#from functools import reduce
import pandas as pd
import glob
import os

# one to three letter AA dict
aa_123 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
          'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
          'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N':'ASN',
          'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S':'SER',
          'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y':'TYR'}

# invert the dict to make it three to one letter
aa_321 = {v: k for k,v in aa_123.items()}

# prefix for the sequence
pref = 'RS09H0'
#%%

# function to get CDR mapping
def cdr_map(row):
    # read the CDR coord file
    cdr_file = 'RS09H0_coord.csv'

    # read the CDR coordinates and store the coordinates in a list
    cdr_df = pd.read_csv(cdr_file)

    # convert this dataframe to a dict with seq id as key
    cdr_df.set_index("seq", drop=True, inplace=True)
    cdr_dict = cdr_df.to_dict(orient = "index")

    # make easy names for the coodinates
    h1s = cdr_dict[pref]['hcdr1_start']
    h1e = cdr_dict[pref]['hcdr1_end']
    h2s = cdr_dict[pref]['hcdr2_start']
    h2e = cdr_dict[pref]['hcdr2_end']
    h3s = cdr_dict[pref]['hcdr3_start']
    h3e = cdr_dict[pref]['hcdr3_end']
    hlen = cdr_dict[pref]['hc_len']
    l1s = cdr_dict[pref]['lcdr1_start'] + hlen
    l1e = cdr_dict[pref]['lcdr1_end'] + hlen
    l2s = cdr_dict[pref]['lcdr2_start'] + hlen
    l2e = cdr_dict[pref]['lcdr2_end'] + hlen
    l3s = cdr_dict[pref]['lcdr3_start'] + hlen
    l3e = cdr_dict[pref]['lcdr3_end'] + hlen

    # convert the prefix to integer for quering the dict

    # check the abpos value and assign
    if row['ab_pos'] >= h1s and row['ab_pos'] <= h1e:
        val = 'H1'
    elif row['ab_pos'] >= h2s and row['ab_pos'] <= h2e:
        val = 'H2'
    elif row['ab_pos'] >= h3s and row['ab_pos'] <= h3e:
        val = 'H3'
    elif row['ab_pos'] >= l1s and row['ab_pos'] <= l1e:
        val = 'L1'
    elif row['ab_pos'] >= l2s and row['ab_pos'] <= l2e:
        val = 'L2'
    elif row['ab_pos'] >= l3s and row['ab_pos'] <= l3e:
        val = 'L3'
    else:
        val = 'FWR'

    return val

# %%
# get the prefix and model number and set file
for fname in glob.glob(pref + '/11_seletopclusts/cluster_1_model_1.ic'):
    bdir, ic_file = os.path.split(fname)
    model = ic_file.split('.')[0]
    seqdir = bdir.split('/')[0]
    print("processing {}, model {} ....".format(bdir, model))
    out_file = os.path.join(bdir, model + '_intf.csv')

    ic_df = pd.read_csv(fname, sep = '\s+', header = None)
    # assign column headers
    ic_df.columns = ['ab_aa3', 'ab_pos', 'ab_chain', 'ag_aa3', 'ag_pos', 'ag_chain']

    # replace three letter code with one letter
    ic_df['ab_aa'] = ic_df['ab_aa3'].map(aa_321)
    ic_df['ag_aa'] = ic_df['ag_aa3'].map(aa_321)

    # make a reduced df with lesser columns
    red_df = ic_df[['ab_pos', 'ab_aa', 'ag_pos', 'ag_aa']]


    # drop duplicates and sort
    red_df = red_df.drop_duplicates()

    # sort by pos
    red_df = red_df.sort_values(by = ['ag_pos'])

    # add the binding interface
    red_df['ab_dom'] = red_df.apply(cdr_map, axis = 1)

    # number of residues
    print("Number of interface pairs = {}".format(red_df.shape[0]))

    # write to file
    red_df.to_csv(out_file, index = False)

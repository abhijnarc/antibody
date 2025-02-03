import pandas as pd
import os 
import sys

pref = sys.argv[1]

# Load the input CSV files with error handling
try:
    input_df = pd.read_csv(pref + '/11_seletopclusts/cluster_1_model_1_intf.csv')
    lcn2_df = pd.read_csv('lcn2.csv')
except pd.errors.EmptyDataError:
    print("One of the input files is empty or not formatted correctly.")
    exit()

# Filter out rows with 'FWR' in the 'ab_dom' column
filtered_input_df = input_df[input_df['ab_dom'] != 'FWR']

# Create sets for M, B, and S categories from lcn2.csv
M_positions = set(lcn2_df['M'].dropna().astype(int))
B_positions = set(lcn2_df['B'].dropna().astype(int))
S_positions = set(lcn2_df['S'].dropna().astype(int))

# Initialize counters for each category
M_count = 0
B_count = 0
S_count = 0

# Get unique positions from the filtered input
unique_positions = set(filtered_input_df['ag_pos'])

# Categorize and count positions
for pos in unique_positions:
    if pos in M_positions:
        M_count += 1
    elif pos in B_positions:
        B_count += 1
    elif pos in S_positions:
        S_count += 1

# Print the counts
print(f"Only MMP9: {M_count}")
print(f"Both MMP9 and SLC17A22: {B_count}")
print(f"Only SLC17A22 : {S_count}")

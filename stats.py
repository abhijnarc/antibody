import subprocess
import pandas as pd
import os
import sys
pref = sys.argv[1]
# Define the scripts to run
scripts = ['/data/Abhijna/scripts/get_top10energy.py', '/data/Abhijna/scripts/intf_fullab.py', '/data/Abhijna/scripts/bindingsurface.py']

# Define the arguments for each script
args = [pref] * len(scripts)

# Run each script sequentially
for script, arg in zip(scripts, args):
    subprocess.run(['python3', script , arg])

# Read the results from the scripts
energy_file = pd.read_csv('energy_' + pref +  '.csv')
binding_surface_output = pd.read_csv(pref + '/binding_surface_counts.csv')

# Extract required values
model = energy_file.iloc[0, 0]  
haddock = energy_file.iloc[0, 1]  
only_mmp9 = binding_surface_output.iloc[0, 0]  
both = binding_surface_output.iloc[0, 1]  
only_slc22a17 = binding_surface_output.iloc[0, 2]   
energy = energy_file.iloc[0, -1] 

# Create a DataFrame with two rows
df = pd.DataFrame({
    "model": [model],  # Same value for both rows
    "haddock": [haddock],  # Same value for both rows
    "only MMP9": [only_mmp9],
    "both MMP9 and SLC22A17": [both],
    "only SLC22A17": [only_slc22a17],
    "Column6": [energy]
})

# Save to CSV
df.to_csv(pref + "_results.csv", index=False)

print(f"CSV file '{pref}_results.csv' created successfully!")
print(pd.read_csv(pref + "_results.csv"))

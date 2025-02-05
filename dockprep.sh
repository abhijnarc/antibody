#!/bin/bash

# Check if a prefix is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <prefix>"
    exit 1
fi

# Set variables
prefix=$1
fasta_file="${prefix}.fa"
igf_pdb="${prefix}_igf.pdb"
clean_pdb="${prefix}_igf_clean.pdb"
coord_file="${prefix}_coord.csv"
res_file="${prefix}_res.txt"
ambig_file="ambig${prefix}.tbl"
unambig_file="unambig${prefix}.tbl"

source "$(conda info --base)/etc/profile.d/conda.sh"
# Activate IgFold environment
conda activate igfold

# Step 1: Run igf_full.py to generate the antibody model
echo "Running igf_full.py..."
python3 /data/Abhijna/scripts/igf_full.py $prefix

#Activate haddock envirnoment
conda activate haddock3

# Step 2: Clean the antibody PDB file
echo "Running clean_ab.sh..."
bash /data/Abhijna/scripts/clean_ab.sh $prefix

# Step 3: Parse the CDR regions from the FASTA file
echo "Running cdr_parse.py..."
python3 /data/Abhijna/scripts/cdr_parse.py --fasta_file $fasta_file --prefix $prefix

# Step 4: Generate the restraints file
echo "Running ab_res.py..."
python3 /data/Abhijna/scripts/ab_res.py $coord_file $res_file

# Step 5: Generate the ambiguous restraints file
echo "Running ambig.sh..."
bash /data/Abhijna/scripts/ambig.sh $prefix

# Step 6: Generate the unambiguous restraints file
echo "Running unambig.sh..."
bash /data/Abhijna/scripts/unambig.sh $prefix

# Step 7: Generate the configuration file
echo "Running generate_config.sh..."
bash /data/Abhijna/scripts/generate_config.sh $prefix

echo "Pipeline completed successfully."

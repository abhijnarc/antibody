#!/bin/bash

# Check if a prefix is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <prefix>"
    exit 1
fi

# Set variables
prefix=$1
filename="ambig${prefix}.tbl"

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate haddock3

# Run haddock3-restraints to generate the file
haddock3-restraints active_passive_to_ambig "${prefix}_res.txt" ag_res.txt > "$filename"

# Check if the file was created successfully
if [ ! -f "$filename" ]; then
    echo "Error: Failed to generate $filename"
    exit 1
fi

# Remove the first line and overwrite the file
tail -n +2 "$filename" > temp_file && mv temp_file "$filename"

echo "Successfully generated and cleaned $filename."

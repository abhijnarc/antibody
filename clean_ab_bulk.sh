#!/bin/bash

# Directory containing input PDB files
input_dir="pdb"  # Update this if necessary
output_dir="igf_clean"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Loop through each PDB file in the input directory
for pdb_file in $input_dir/*.pdb; do
    # Extract the prefix (file name without extension)
    base_name=$(basename "$pdb_file" .pdb)
    
    # Set input and output filenames
    input=$pdb_file
    output=${output_dir}/${base_name}_clean.pdb
    hc=${output_dir}/${base_name}_H.pdb
    lc=${output_dir}/${base_name}_L.pdb

    echo "Processing $pdb_file..."
    
    # Clean the heavy chain
    cat $input | pdb_tidy -strict | pdb_selchain -H | pdb_delhetatm | \
        pdb_fixinsert | pdb_keepcoord | pdb_tidy -strict > $hc
    
    # Clean the light chain
    cat $input | pdb_tidy -strict | pdb_selchain -L | pdb_delhetatm | \
        pdb_fixinsert | pdb_keepcoord | pdb_tidy -strict > $lc
    
    # Merge the two chains and renumber
    pdb_merge $hc $lc | pdb_chain -A | pdb_chainxseg | pdb_reres | \
        pdb_tidy -strict > $output

    echo "Cleaned PDB saved to $output"
done

echo "All PDB files have been processed."

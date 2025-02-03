#/bin/bash
# script for cleanup of antibody pdb file
# get file prefix from command line
pref=$1
# set filenames
input=${pref}.pdb
output=${pref}_clean.pdb
hc=${pref}_H.pdb
lc=${pref}_L.pdb
# first the heavy chain
cat $input | pdb_tidy -strict | pdb_selchain -H | pdb_delhetatm | \
         pdb_fixinsert | pdb_keepcoord | pdb_tidy -strict > $hc
# next the light chain
cat $input | pdb_tidy -strict | pdb_selchain -L | pdb_delhetatm | \
         pdb_fixinsert | pdb_keepcoord | \
         pdb_tidy -strict > $lc
# merge the two and renumber
pdb_merge $hc $lc  | pdb_chain -A | \
        pdb_chainxseg | pdb_reres | pdb_tidy -strict > $output

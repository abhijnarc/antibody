import os
from pyfasta import Fasta
from igfold import IgFoldRunner
import sys

# Define the directory containing FASTA files
input_dir = "sequences"  # Update to your actual directory
output_dir = "pdb"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Create an igfold object
igfold = IgFoldRunner()

# Loop through each FASTA file in the input directory
for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith(".fa"):  # Ensure we're processing only FASTA files
        pref = os.path.splitext(fasta_file)[0]  # Get the file prefix (without .fa)
        input_path = os.path.join(input_dir, fasta_file)
        output_path = os.path.join(output_dir, f"{pref}_igf.pdb")
        
        print(f"Processing {fasta_file}...")
        
        # Read the FASTA file
        this_seq = Fasta(input_path)
        
        # Set the Fasta to an antibody dict
        abseq = {
            'H': str(this_seq['VH']),
            'L': str(this_seq['VL'])
        }

        # Run IgFold
        out = igfold.fold(
            output_path,       # Output PDB file
            sequences=abseq,   # Antibody sequences
            do_refine=True,    # Refine the antibody structure with PyRosetta
            use_openmm=True,   # Use OpenMM for refinement
            do_renum=False     # Do not renumber predicted antibody structure
        )
        
        print(f"Prediction for {fasta_file} completed. Output saved to {output_path}")

print("All predictions are completed.")

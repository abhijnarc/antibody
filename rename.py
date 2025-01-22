from Bio import SeqIO
import sys

# Input and output FASTQ file names
input_fasta = "S1C1_light.fa"  # Replace with your input FASTQ file
output_fasta = "S1C1_pseudo_light.fa"  # Output file for sequences with renamed IDs

# Define the old and new prefix
old_prefix = "assemble"  # Replace with the current prefix
new_prefix = "lcn2_"  # Replace with the new prefix

# Define a function to modify the ID (rename prefix while keeping suffix)
def modify_id(old_id):
    # Check if the ID starts with the old prefix
    if old_id.startswith(old_prefix):
        # Replace only the old prefix with the new prefix, leaving the rest of the ID (suffix) intact
        new_id = old_id.replace(old_prefix, new_prefix, 1)
        return new_id
    else:
        # Return the original ID if it doesn't match the old prefix
        return old_id

# Open the input FASTQ file for reading
with open(input_fasta, "r") as input_handle:
    # Parse the FASTQ file
    records = SeqIO.parse(input_handle, "fasta")
    
    # Open the output FASTQ file for writing
    with open(output_fasta, "w") as output_handle:
        # Loop through each sequence record
        for record in records:
            # Modify the ID using the function
            record.id = modify_id(record.id)
            SeqIO.write(record, output_handle, "fasta")

print(f"Renamed sequences saved to {output_fasta}")

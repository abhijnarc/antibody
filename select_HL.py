# -*- coding: utf-8 -*-
"""
Modified Script to save the H chain and L chain with the highest coverage in seq_summary.csv.
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys

pref = sys.argv[1]
# Input and output file setup
input_fasta = pref +  "_filtered.fa"  # Provided filtered FASTA file
output_summary = "antibody_seq_summary.csv"
output_filtered_fasta = "antibody_filtered.fa"
output_vh_vl_fasta = "antibody_vh_vl.fasta"
aa_vh_vl = "aa_vh_vl.fa"

# Parse the FASTA file
sequences = SeqIO.parse(input_fasta, "fasta")

# Initialize containers
seq_summary = []
filtered_sequences = []
filter_count = 0

# Process sequences
for record in sequences:
    header = record.description
    contig = header.split()[0]
    coverage = float(header.split()[2])  # Assume 3rd element is coverage
    
    # Get V loci information
    v_loci_raw = header.split()[3]
    v_loci = v_loci_raw.split('*')[0]
    
    # Determine chain type
    chain = "heavy" if "IGH" in v_loci else "light"

    # Extract CDR scores using regex
    cdr1_score = int(re.search(r"CDR1\(\d+-\d+\):(\d+)", header).group(1))
    cdr2_score = int(re.search(r"CDR2\(\d+-\d+\):(\d+)", header).group(1))
    cdr3_score = int(re.search(r"CDR3\(\d+-\d+\):(\d+)", header).group(1))
    
    # Calculate total score and check filtering criteria
    total_score = cdr1_score + cdr2_score + cdr3_score
    is_valid = (cdr1_score > 0) and (cdr2_score > 0) and (cdr3_score > 0) and (coverage >= 10)
    
    if is_valid:
        # Increment filter count and append to filtered list
        filter_count += 1
        seq_summary.append({
            "contig": contig,
            "chain": chain,
            "coverage": coverage,
            "cdr1_score": cdr1_score,
            "cdr2_score": cdr2_score,
            "cdr3_score": cdr3_score,
            "total_score": total_score,
            "sequence": str(record.seq)
        })
        filtered_sequences.append(record)

# Convert summary to DataFrame
summary_df = pd.DataFrame(seq_summary)

# Sort by chain type and coverage, descending
sorted_df = summary_df.sort_values(by=["chain", "coverage"], ascending=[True, False])

# Save the filtered sequences to CSV and FASTA
sorted_df.to_csv(output_summary, index=False)
SeqIO.write(filtered_sequences, output_filtered_fasta, "fasta")

# Select the H and L chains with the highest coverage
heavy_chain = sorted_df[sorted_df["chain"] == "heavy"].iloc[0]
light_chain = sorted_df[sorted_df["chain"] == "light"].iloc[0]

# Create SeqRecord objects for VH and VL
vh_record = SeqRecord(
    Seq(heavy_chain["sequence"]),
    id="VH",
    description=f"Chain: {heavy_chain['chain']}, Coverage: {heavy_chain['coverage']}, Total Score: {heavy_chain['total_score']}"
)

vl_record = SeqRecord(
    Seq(light_chain["sequence"]),
    id="VL",
    description=f"Chain: {light_chain['chain']}, Coverage: {light_chain['coverage']}, Total Score: {light_chain['total_score']}"
)

# Save VH and VL to a new FASTA file
SeqIO.write([vh_record, vl_record], output_vh_vl_fasta, "fasta")


# Print summary
print(f"Filtered sequences: {filter_count}")
print(f"Summary saved to: {output_summary}")
print(f"Filtered sequences saved to: {output_filtered_fasta}")
print(f"VH and VL saved to: {output_vh_vl_fasta}")


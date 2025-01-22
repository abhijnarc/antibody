# -*- coding: utf-8 -*-
"""
Script to parse, filter sequences for antibody construction,
save first heavy and light chains as H and L chains.
"""

import pandas as pd
from Bio import SeqIO, Seq
import re
import sys

pref = sys.argv[1]
# Input and output file setup
input_fasta = pref +"_filtered.fa"  # Provided filtered FASTA file
output_summary = "seq_summary.csv"
output_filtered_fasta = "filtered.fa"
output_hl_fasta = "H_L_chains.fasta"
aa_hl = "aa_H_L.fa"

# Parse the FASTA file
sequences = SeqIO.parse(input_fasta, "fasta")

# Initialize containers
seq_summary = []
filtered_sequences = []
filter_count = 0

# Variables to store the first heavy and light chain
first_heavy_chain = None
first_light_chain = None

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

        # Save the first heavy and light chain
        if chain == "heavy" and first_heavy_chain is None:
            first_heavy_chain = record
        if chain == "light" and first_light_chain is None:
            first_light_chain = record

# Convert summary to DataFrame
summary_df = pd.DataFrame(seq_summary)

# Save the filtered sequences to CSV and FASTA
summary_df.to_csv(output_summary, index=False)
SeqIO.write(filtered_sequences, output_filtered_fasta, "fasta")

# Save H and L chains to a new FASTA file
if first_heavy_chain and first_light_chain:
    h_record = SeqIO.SeqRecord(
        first_heavy_chain.seq,
        id="H",
        description=f"First Heavy Chain, Coverage: {coverage}, Total Score: {total_score}"
    )
    l_record = SeqIO.SeqRecord(
        first_light_chain.seq,
        id="L",
        description=f"First Light Chain, Coverage: {coverage}, Total Score: {total_score}"
    )
    SeqIO.write([h_record, l_record], output_hl_fasta, "fasta")

    
# Print summary
print(f"Filtered sequences: {filter_count}")
print(f"Summary saved to: {output_summary}")
print(f"Filtered sequences saved to: {output_filtered_fasta}")
print(f"H and L chains saved to: {output_hl_fasta}")


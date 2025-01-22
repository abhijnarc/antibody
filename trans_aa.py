from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

pref = sys.argv[1]
# Input and output FASTA file names
input_fasta = pref +".fa"  # Replace with your nucleotide FASTA file
output_fasta = pref +"_aa.fa"  # Output file for amino acid sequences

# Parse the input FASTA file
nucleotide_records = SeqIO.parse(input_fasta, "fasta")

# Create a list to store translated sequences
amino_acid_records = []

# Iterate over nucleotide sequences and translate them
for record in nucleotide_records:
    # Translate the nucleotide sequence to amino acids (stop at stop codon)
    translated_seq = record.seq.translate(to_stop=True)
    
    # Create a new SeqRecord for the amino acid sequence
    translated_record = SeqRecord(
        translated_seq,
        id=record.id,  # Keep the same ID as the original sequence
    )
    
    # Append the translated sequence to the list
    amino_acid_records.append(translated_record)

# Save the translated sequences to the output FASTA file
SeqIO.write(amino_acid_records, output_fasta, "fasta")

print(f"Translated sequences saved to {output_fasta}")

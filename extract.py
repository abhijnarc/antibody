import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_cdr_pseudo_sequences(input_file, output_file):
    """
    Extract CDR1, CDR2, and CDR3 sequences from the FASTA headers (description lines),
    concatenate them, and create pseudo-sequences in single-line FASTA format.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
    """
    # Compile regex patterns for extracting CDR regions from the description
    cdr_patterns = {
        'CDR1': re.compile(r'CDR1\((\d+)-(\d+)\):([A-Za-z]+)'),
        'CDR2': re.compile(r'CDR2\((\d+)-(\d+)\):([A-Za-z]+)'),
        'CDR3': re.compile(r'CDR3\((\d+)-(\d+)\):([A-Za-z]+)')
    }

    output_records = []

    # Process each sequence in the FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        header = record.description

        # Extract CDR regions from the header
        cdr_sequences = []
        for cdr_name, pattern in cdr_patterns.items():
            match = pattern.search(header)
            if match:
                # Extract the sequence and adjust for any potential errors
                cdr_sequence = match.group(3).strip()
                cdr_sequences.append(cdr_sequence)

        # Concatenate CDR sequences into a pseudo-sequence
        if cdr_sequences:
            pseudo_sequence = ''.join(cdr_sequences)
            output_record = SeqRecord(
                Seq(pseudo_sequence),
                id=record.id,
                description="Pseudo-sequence of concatenated CDRs"
            )
            output_records.append(output_record)

    # Write the output FASTA file
    with open(output_file, "w") as f:
        SeqIO.write(output_records, f, "fasta")

# Example usage
# extract_cdr_pseudo_sequences("input.fasta", "output.fasta")
if __name__ == "__main__":
    input_file = "LCN2_filtered.fa"
    output_file = "LCN2_CDRs.fa"
    extract_cdr_pseudo_sequences(input_file, output_file)

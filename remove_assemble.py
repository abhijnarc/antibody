from Bio import SeqIO

def edit_descriptions(input_file, output_file, prefix="assemble"):
    """
    Edits out instances of a specific prefix from the description lines in a FASTA file.
    
    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file.
        prefix (str): Prefix to remove from descriptions (default: "assemble").
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # Remove all occurrences of the prefix from the description
            record.description = " ".join(
                word for word in record.description.split() if not word.startswith(prefix)
            )
            SeqIO.write(record, outfile, "fasta")

# Example Usage
input_fasta = "renamed_heavy.fa"
output_fasta = "lcn2_heavy.fa"
edit_descriptions(input_fasta, output_fasta)

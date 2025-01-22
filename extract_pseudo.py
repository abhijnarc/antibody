import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_cdr_sequences(input_file, output_file):
    """
    Extract CDR1, CDR2, and CDR3 sequences from FASTA headers and concatenate them.
    Output is written in single-line FASTA format.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output FASTA file
    """
    # Compile regex patterns for finding CDR regions
    cdr_patterns = {
        'CDR1': re.compile(r'CDR1\((\d+)-(\d+)\):([A-Za-z]+)'),
        'CDR2': re.compile(r'CDR2\((\d+)-(\d+)\):([A-Za-z]+)'),
        'CDR3': re.compile(r'CDR3\((\d+)-(\d+)\):([A-Za-z]+)')
        }
    
    output_records = []
    
    # Process each sequence in the FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        header = record.description
        sequence = str(record.seq)
        
        # Extract CDR regions
        cdr_sequences = []
        for cdr_name, pattern in cdr_patterns.items():
            match = pattern.search(header)
            if match:
                start, end = int(match.group(1)), int(match.group(2))
                cdr_seq = sequence[start-1:end]  # Convert to 0-based indexing
                cdr_sequences.append(cdr_seq)
        
        if len(cdr_sequences) == 3:  # Only process if all CDRs are found
            # Concatenate CDR sequences
            combined_seq = ''.join(cdr_sequences)
            
            # Create new record with concatenated sequence
            new_record = SeqRecord(
                Seq(combined_seq),
                id=record.id,
                description=f"concatenated_CDRs_{len(combined_seq)}bp"
            )
            output_records.append(new_record)
    
    # Write sequences in single-line FASTA format
    with open(output_file, 'w') as f:
        for record in output_records:
            header = f">{record.id} {record.description}\n"
            sequence = f"{str(record.seq)}\n"
            f.write(header)
            f.write(sequence)
            
    print(f"Processed sequences written to {output_file}")

# Example usage
if __name__ == "__main__":
    input_file = "LCN2_filtered.fa"
    output_file = "LCN2_CDRs.fa"
    extract_cdr_sequences(input_file, output_file)

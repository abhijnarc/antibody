import pandas as pd
from Bio import SeqIO
import argparse
import os

def parse_fasta_and_process(fasta_file, prefix):
    heavy_sequence = None
    light_sequence = None

    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"Processing record: {record.id}")  # Debugging
        if "VH" in record.id.upper():
            heavy_sequence = str(record.seq)
        elif "VL" in record.id.upper():
            light_sequence = str(record.seq)

    if heavy_sequence:
        append_heavy_chain(heavy_sequence, prefix)
    else:
        print("No heavy chain sequence found in the FASTA file.")

    if light_sequence:
        append_light_chain(light_sequence, prefix)
    else:
        print("No light chain sequence found in the FASTA file.")

def append_heavy_chain(heavy_sequence, prefix):
    print(f"Processing heavy chain sequence: {heavy_sequence[:10]}...")  # Debugging
    url = f"http://www.abysis.org/abysis/sequence_input/key_annotation/key_annotation.cgi?aa_sequence={heavy_sequence}&nuc_sequence=&translation=sixft&humanorganism=on"
    try:
        tables = pd.read_html(url, match='Sequence Fragment', index_col=None)
        if not tables:
            print("No matching tables found.")
            return
        extracted_table = tables[5]
        cdr_table = extracted_table[extracted_table['Region'].isin(['CDR-H1', 'CDR-H2', 'CDR-H3', 'NaN'])]
        cdr_coords = {}
        for region in ['CDR-H1', 'CDR-H2', 'CDR-H3']:
            row = cdr_table[cdr_table['Region'] == region]
            if not row.empty:
                start, end = map(int, row['Residues'].iloc[0].split(' - '))
                cdr_coords[f'{region.lower()}_start'] = start
                cdr_coords[f'{region.lower()}_end'] = end
            else:
                cdr_coords[f'{region.lower()}_start'] = None
                cdr_coords[f'{region.lower()}_end'] = None

        # Get the total length
        total_length = extracted_table['Length'].dropna().iloc[7]

        hc_data = {
            'seq': [prefix],
            'hcdr1_start': [cdr_coords.get('cdr-h1_start')],
            'hcdr1_end': [cdr_coords.get('cdr-h1_end')],
            'hcdr2_start': [cdr_coords.get('cdr-h2_start')],
            'hcdr2_end': [cdr_coords.get('cdr-h2_end')],
            'hcdr3_start': [cdr_coords.get('cdr-h3_start')],
            'hcdr3_end': [cdr_coords.get('cdr-h3_end')],
            'hc_len': [total_length]
        }
        hc_df = pd.DataFrame(hc_data)
        output_file = f'{prefix}_coord.csv'
        hc_df.to_csv(output_file, index=False)
        print("\nHeavy Chain Data:")
        print(hc_df)
        print(f"\nData saved to {output_file}")
    except ValueError as ve:
        print(f"ValueError: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred while processing heavy chain: {e}")

def append_light_chain(light_sequence, prefix):
    print(f"Processing light chain sequence: {light_sequence[:10]}...")  # Debugging
    url = f"http://www.abysis.org/abysis/sequence_input/key_annotation/key_annotation.cgi?aa_sequence={light_sequence}&nuc_sequence=&translation=sixft&humanorganism=on"
    try:
        tables = pd.read_html(url, match='Sequence Fragment', index_col=None)
        if not tables:
            print("No matching tables found.")
            return
        extracted_table = tables[5]
        cdr_table = extracted_table[extracted_table['Region'].isin(['CDR-L1', 'CDR-L2', 'CDR-L3'])]
        cdr_coords = {}
        for region in ['CDR-L1', 'CDR-L2', 'CDR-L3']:
            row = cdr_table[cdr_table['Region'] == region]
            if not row.empty:
                start, end = map(int, row['Residues'].iloc[0].split(' - '))
                cdr_coords[f'{region.lower()}_start'] = start
                cdr_coords[f'{region.lower()}_end'] = end
            else:
                cdr_coords[f'{region.lower()}_start'] = None
                cdr_coords[f'{region.lower()}_end'] = None
        lc_data = {
            'lcdr1_start': [cdr_coords.get('cdr-l1_start')],
            'lcdr1_end': [cdr_coords.get('cdr-l1_end')],
            'lcdr2_start': [cdr_coords.get('cdr-l2_start')],
            'lcdr2_end': [cdr_coords.get('cdr-l2_end')],
            'lcdr3_start': [cdr_coords.get('cdr-l3_start')],
            'lcdr3_end': [cdr_coords.get('cdr-l3_end')]
        }
        lc_df = pd.DataFrame(lc_data)
        output_file = f'{prefix}_coord.csv'
        try:
            existing_df = pd.read_csv(output_file)
        except FileNotFoundError:
            print(f"{output_file} not found. Please ensure the heavy chain script runs first.")
            return
        combined_df = pd.concat([existing_df, lc_df], axis=1)
        combined_df.to_csv(output_file, index=False)
        print("\nAppended Data:")
        print(combined_df)
        print(f"\nData appended to {output_file}")
    except ValueError as ve:
        print(f"ValueError: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred while processing light chain: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse FASTA file and process sequences.")
    parser.add_argument("--fasta_file", type=str, required=True, help="The path to the FASTA file containing sequences.")
    parser.add_argument("--prefix", type=str, required=True, help="The prefix for the output coord CSV file.")
    args = parser.parse_args()
    parse_fasta_and_process(args.fasta_file, args.prefix)

import pandas as pd

def append_light_chain(sequence):
    """
    Scrapes a table from the Abysis page using the provided light chain sequence,
    extracts start and end residues for LCDR1, LCDR2, and LCDR3, 
    and appends this information to an existing CSV file in the specified format.

    Parameters:
        sequence (str): The amino acid sequence for the light chain.
        seq_name (str): A name or identifier for the sequence.

    Returns:
        None: Appends the data to the CSV and prints the appended data.
    """
    # Construct the URL dynamically with the user's input
    url = f"http://www.abysis.org/abysis/sequence_input/key_annotation/key_annotation.cgi?aa_sequence={sequence}&nuc_sequence=&translation=sixft&humanorganism=on"

    try:
        # Use pandas to read the table containing 'Sequence Fragment'
        tables = pd.read_html(url, match='Sequence Fragment', index_col=None)

        if not tables:
            print("No matching tables found.")
            return

        # Extract the desired table
        extracted_table = tables[5]

        # Filter rows for LCDR1, LCDR2, and LCDR3
        cdr_table = extracted_table[extracted_table['Region'].isin(['CDR-L1', 'CDR-L2', 'CDR-L3'])]

        # Parse the start and end residues
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

        # Prepare the light chain data
        lc_data = {
            'lcdr1_start': [cdr_coords.get('cdr-l1_start')],
            'lcdr1_end': [cdr_coords.get('cdr-l1_end')],
            'lcdr2_start': [cdr_coords.get('cdr-l2_start')],
            'lcdr2_end': [cdr_coords.get('cdr-l2_end')],
            'lcdr3_start': [cdr_coords.get('cdr-l3_start')],
            'lcdr3_end': [cdr_coords.get('cdr-l3_end')]
        }

        lc_df = pd.DataFrame(lc_data)

        # Load the existing CSV file
        output_file = 'RS09H0_coord.csv'
        try:
            existing_df = pd.read_csv(output_file)
        except FileNotFoundError:
            print(f"{output_file} not found. Please ensure the heavy chain script runs first.")
            return

        # Append the light chain data to the existing data
        combined_df = pd.concat([existing_df, lc_df], axis=1)

        # Save the updated file
        combined_df.to_csv(output_file, index=False)

        # Print the combined data
        print("\nAppended Data:")
        print(combined_df)
        print(f"\nData appended to {output_file}")

    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

import argparse

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Scrape CDR table from Abysis.")
    parser.add_argument("--sequence", type=str, required=True, help="The amino acid sequence for the chain.")

    # Parse arguments
    args = parser.parse_args()
    sequence = args.sequence

    # Call the scrape function with the provided arguments
    append_light_chain(sequence)

import pandas as pd

def scrape_and_format_table(sequence, seq_name):
    """
    Scrapes a table from the Abysis page using the provided amino acid sequence,
    formats the data to include start and end residues for CDR-H1, CDR-H2, and CDR-H3,
    and saves the data in the desired CSV format.

    Parameters:
        sequence (str): The amino acid sequence to query.
        seq_name (str): A name or identifier for the sequence.

    Returns:
        None: Prints the processed data and saves it to a CSV file.
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

        # Filter rows for CDR-H1, CDR-H2, and CDR-H3
        cdr_table = extracted_table[extracted_table['Region'].isin(['CDR-H1', 'CDR-H2', 'CDR-H3'])]

        # Parse the start and end residues
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
        total_length = extracted_table['Length'].dropna().astype(int).sum()

        # Prepare the final dataframe
        formatted_data = {
            'seq': [seq_name],
            'hcdr1_start': [cdr_coords.get('cdr-h1_start')],
            'hcdr1_end': [cdr_coords.get('cdr-h1_end')],
            'hcdr2_start': [cdr_coords.get('cdr-h2_start')],
            'hcdr2_end': [cdr_coords.get('cdr-h2_end')],
            'hcdr3_start': [cdr_coords.get('cdr-h3_start')],
            'hcdr3_end': [cdr_coords.get('cdr-h3_end')],
            'hc_len': [total_length]
        }

        df = pd.DataFrame(formatted_data)

        # Save the processed data to a CSV file
        output_file = 'cdr_coord.csv'
        df.to_csv(output_file, index=False)

        # Print the formatted data
        print("\nFormatted Data:")
        print(df)
        print(f"\nData saved to {output_file}")

    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# Main function to take user input and call the scraper
import argparse

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Scrape CDR table from Abysis.")
    parser.add_argument("--sequence", type=str, required=True, help="The amino acid sequence for the chain.")
    parser.add_argument("--name", type=str, required=True, help="A name or identifier for the sequence.")

    # Parse arguments
    args = parser.parse_args()
    sequence = args.sequence
    seq_name = args.name

    # Call the scrape function with the provided arguments
    scrape_and_format_table(sequence, seq_name)

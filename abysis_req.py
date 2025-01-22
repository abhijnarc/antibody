import pandas as pd

def scrape_table(sequence):
    """
    Scrapes a table containing 'Sequence Fragment' from the Abysis page using the provided amino acid sequence.

    Parameters:
        sequence (str): The amino acid sequence to query.

    Returns:
        None: Prints the extracted table or an error message.
    """
    # Construct the URL dynamically with the user's input
    url = f"http://www.abysis.org/abysis/sequence_input/key_annotation/key_annotation.cgi?aa_sequence={sequence}&nuc_sequence=&translation=sixft&humanorganism=on"

    try:
        # Use pandas to read the table containing 'Sequence Fragment'
        tables = pd.read_html(url, match='Sequence Fragment', index_col=None)

        # Display the number of tables found
        print(f"Number of tables found: {len(tables)}")

        # Extract and print the desired table (assuming the first match)
        if tables:
            desired_table = tables[4]  # Retrieve the first matching table
            print("\nExtracted Table:")
            print(desired_table)
        else:
            print("No matching tables found.")

    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# Main function to take user input and call the scraper
if __name__ == "__main__":
    # Prompt the user for the amino acid sequence
    user_sequence = input("Enter the amino acid sequence: ").strip()

    # Call the scrape function with the user's input
    scrape_table(user_sequence)

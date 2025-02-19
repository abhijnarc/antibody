import csv
import sys

def generate_restraints(input_csv, output_txt):
    try:
        with open(input_csv, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            with open(output_txt, 'w') as txtfile:
                for row in reader:
                    # Extract CDR coordinates and generate ranges
                    cdr_ranges = [
                        range(int(row['hcdr1_start']), int(row['hcdr1_end']) + 1),
                        range(int(row['hcdr2_start']), int(row['hcdr2_end']) + 1),
                        range(int(row['hcdr3_start']), int(row['hcdr3_end']) + 1),
                        range(int(row['lcdr1_start']), int(row['lcdr1_end']) + 1),
                        range(int(row['lcdr2_start']), int(row['lcdr2_end']) + 1),
                        range(int(row['lcdr3_start']), int(row['lcdr3_end']) + 1)
                    ]
                    # Flatten the list of ranges and write to the output file
                    all_numbers = sorted([num for cdr_range in cdr_ranges for num in cdr_range])
                    unique_numbers = sorted(set(all_numbers))  # Get only unique numbers
                    txtfile.write(' '.join(map(str, unique_numbers)) + '\n\n')  # Add an extra empty line
        print(f"Restraint file generated: {output_txt}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ab_res.py <input_csv> <output_txt>")
    else:
        generate_restraints(sys.argv[1], sys.argv[2])

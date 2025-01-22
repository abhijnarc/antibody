import csv
import sys

def generate_restraints(input_csv, output_txt):
    try:
        with open(input_csv, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            with open(output_txt, 'w') as txtfile:
                for row in reader:
                    # Extract CDR coordinates
                    cdr_coords = [
                        row['hcdr1_start'], row['hcdr1_end'],
                        row['hcdr2_start'], row['hcdr2_end'],
                        row['hcdr3_start'], row['hcdr3_end'],
                        row['lcdr1_start'], row['lcdr1_end'],
                        row['lcdr2_start'], row['lcdr2_end'],
                        row['lcdr3_start'], row['lcdr3_end']
                    ]
                    # Write coordinates to the output file
                    txtfile.write(' '.join(cdr_coords) + '\n')
                    txtfile.write('\n')  # Add an empty line
        print(f"Restraint file generated: {output_txt}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_restraints.py <input_csv> <output_txt>")
    else:
        input_csv = sys.argv[1]
        output_txt = sys.argv[2]
        generate_restraints(input_csv, output_txt)

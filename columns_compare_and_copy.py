# -*- coding: cp1250 -*-

# This script copies rows from the first CSV file to a new file.
# It only copies rows where the first column contains values 
# that also appear in the second CSV file (reference file).
# Values in the second file must be in a single column, with each value 
# in a separate row.

import csv

def compare_and_remove(input_file, reference_file, output_file):
    # Create a set to store the values from the reference file
    molecules = set()

    # Read the MOLECULE column from the second file (reference file) and store in the set
    with open(reference_file, 'r') as ref_csv:
        ref_reader = csv.reader(ref_csv)
        next(ref_reader)  # Skip the header
        for row in ref_reader:
            molecule = row[0]  # Assuming the relevant values are in the first column
            molecules.add(molecule)

    # Read the input file and compare the MOLECULE column values to the reference set.
    # Write only the rows where the MOLECULE matches a value in the reference set.
    with open(input_file, 'r') as input_csv, open(output_file, 'w', newline='') as output_csv:
        input_reader = csv.reader(input_csv)
        output_writer = csv.writer(output_csv)

        # Copy the header from the input file to the output file
        header = next(input_reader)
        output_writer.writerow(header)

        # Process each row and check if the first column matches a value from the reference file
        for row in input_reader:
            molecule = row[0]
            if molecule in molecules:
                output_writer.writerow(row)  # Write the row to the output file if it matches

    print("Finished! Relevant rows have been saved to:", output_file)

# Prompt the user to input file names for processing
input_file = input("Enter the name of the input file: ")
reference_file = input("Enter the name of the reference file: ")
output_file = input("Enter the name of the output file: ")

# Run the function to compare and filter the rows
compare_and_remove(input_file, reference_file, output_file)

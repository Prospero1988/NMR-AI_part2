# -*- coding: utf-8 -*-

import os
import csv
import sys

def add_two_zero_columns_to_csv(file_path):
    temp_file_path = file_path + '.tmp'

    with open(file_path, 'r', newline='') as infile, open(temp_file_path, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        # Add headers "zero1" and "zero2" if the CSV has headers
        headers = next(reader, None)
        if headers is not None:
            writer.writerow(headers + ['zero1','zero2','zero3', 'zero4', 'zero5', 'zero6', 'zero7', 'zero8', 'zero9', 'zero10', 'zero11', 'zero12'])
        else:
            # If there are no headers, just add "zero1" and "zero2" headers
            writer.writerow(['zero1','zero2','zero3', 'zero4', 'zero5', 'zero6', 'zero7', 'zero8', 'zero9', 'zero10', 'zero11', 'zero12'])

        # Add the two columns with "0" to each row
        for row in reader:
            writer.writerow(row + ['0', '0' ,'0', '0', '0', '0', '0', '0', '0', '0', '0', '0'])
    
    os.replace(temp_file_path, file_path)

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1)

    current_directory = sys.argv[1]

    if not os.path.isdir(current_directory):
        print(f"The path {current_directory} is not a valid directory.")
        sys.exit(1)

    csv_files = [f for f in os.listdir(current_directory) if f.endswith('.csv')]

    for csv_file in csv_files:
        full_path = os.path.join(current_directory, csv_file)
        print(f"Processing file: {csv_file}")
        add_two_zero_columns_to_csv(full_path)

    print("All CSV files have been processed. Twelve columns with headers 'zero','zero2', etc. have been added to the end of each file.")

if __name__ == "__main__":
    main()

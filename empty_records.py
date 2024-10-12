# -*- coding: utf-8 -*-

import csv

def check_empty_records(file_name):
    # Initialize a counter for empty records and a list to store rows with empty records
    total_empty_records = 0
    rows_with_empty_records = []

    # Open the CSV file for reading
    with open(file_name, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)

        # Loop through each row and check for empty cells
        for row_idx, row in enumerate(csv_reader, start=1):
            # Count the number of empty records (cells) in the current row
            empty_records_in_row = sum(not record.strip() for record in row)
            total_empty_records += empty_records_in_row

            # If there are any empty records, add the row index and count to the list
            if empty_records_in_row > 0:
                rows_with_empty_records.append((row_idx, empty_records_in_row))

    # Return the total count of empty records and the list of rows with empty records
    return total_empty_records, rows_with_empty_records

# Get the file name from the user
file_name = input("Enter the name of the CSV file: ")
print()

# Call the function to check for empty records
total, rows = check_empty_records(file_name)

# Display the total number of empty records
print(f"Total number of empty records: {total}")
print()

# If there are no empty records, exit the script
if total == 0:
    exit()
else:
    # If there are empty records, display the row numbers and the count of empty records in each row
    print("Rows with empty records:")
    for row_num, empty_count in rows:
        print(f"Row: {row_num}, Number of empty records: {empty_count}")
print()

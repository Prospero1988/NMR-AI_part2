# -*- coding: utf-8 -*-

import csv

def remove_rows_with_empty_cells(csv_file):
    # Opening the CSV file for reading
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        rows = list(reader)  # Reading all rows into a list

    # Checking each row and removing those that contain empty cells
    cleaned_rows = [row for row in rows if all(cell.strip() for cell in row)]

    # Opening the CSV file for writing and saving the cleaned rows
    with open(csv_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(cleaned_rows)  # Writing the cleaned rows back to the file

    print("Rows with empty cells have been removed from the CSV file.")

# Example usage
csv_file = input("Enter the name of the CSV file: ")
remove_rows_with_empty_cells(csv_file)

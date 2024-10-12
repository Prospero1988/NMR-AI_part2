# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 18:30:17 2024

@author: aleniak
"""

import pandas as pd
import random
import os
import sys

# Check if a CSV file is provided as an argument
if len(sys.argv) < 2:
    raise ValueError("Please provide a CSV file as an argument!")

# Load the CSV file with pairs of file paths to be mixed
csv_file_with_paths = sys.argv[1]
file_pairs = pd.read_csv(csv_file_with_paths)

# Iterate through each pair of files in the CSV file
for index, row in file_pairs.iterrows():
    file1_path = row[0]
    file2_path = row[1]

    # Check if both files exist
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        raise ValueError(f"Files do not exist: {file1_path}, {file2_path}")
    
    # Load the CSV files
    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)

    # Check if the number of rows and columns match between the two files
    num_rows, num_cols = df1.shape
    if df2.shape != (num_rows, num_cols):
        raise ValueError(f"Files {file1_path} and {file2_path} have different numbers of rows or columns.")

    # Copy the column headers from the first file to the second
    df2.columns = df1.columns

    # Randomly select rows and combine the DataFrames
    percentages = [0.2, 0.4, 0.6, 0.8]
    output_files = []

    for pct in percentages:
        num_rows_sample = int(num_rows * pct)

        for draw in range(1, 11):  # Repeat 10 times for each percentage
            # Randomly select the appropriate percentage of rows from the first file
            indices1 = random.sample(range(num_rows), num_rows_sample)

            # Remove the selected rows from the second file
            df2_sample = df2.drop(indices1)

            # Combine the randomly selected rows from the first file with the remaining rows from the second file
            df1_sample = df1.loc[indices1]
            df_combined = pd.concat([df1_sample, df2_sample])

            # Ensure the resulting DataFrame has the same number of rows and columns as the original files
            assert df_combined.shape == (num_rows, num_cols), "The number of rows or columns in the resulting file does not match the original."

            # Save the combined data to a CSV file
            folder_name = f"mixed_files_{index+1}"
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)

            output_filename = f'{folder_name}/mixed_{int(pct*100)}_percent_draw{draw:02}.csv'
            df_combined.to_csv(output_filename, index=False)
            output_files.append(output_filename)

            # Information about the mixing process
            print(f"Creating file: {output_filename}")
            print(f"Number of rows from file '{file1_path}': {num_rows_sample}")
            print(f"Number of rows from file '{file2_path}': {num_rows - num_rows_sample}")

    # Display information about the processed files
    print(f"Processed files: {file1_path}, {file2_path}")
    print(f"Number of rows: {num_rows}, Number of columns: {num_cols}")
    print(f"Created files: {', '.join(output_files)}")

print("Script completed successfully.")

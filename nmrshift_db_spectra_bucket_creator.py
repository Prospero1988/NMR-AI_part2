# -*- coding: utf-8 -*-

import os
import csv
import argparse
import numpy as np

def create_buckets():
    # Create 500 buckets from -1 to 13 inclusive
    bucket_range = np.linspace(-1, 13, 501)
    buckets = {i: 0 for i in range(500)}
    return bucket_range, buckets

def find_bucket_index(value, bucket_range):
    # Find the bucket index for the given value
    if value < -1 or value > 13:
        return None
    for i in range(len(bucket_range) - 1):
        if bucket_range[i] <= value < bucket_range[i + 1]:
            return i
    return len(bucket_range) - 2  # Last bucket if value is exactly 13

def process_file(file_path, output_dir, bucket_range):
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        buckets = {i: 0 for i in range(500)}
        error_values = []
        
        for row in reader:
            try:
                value = float(row[0])
                bucket_index = find_bucket_index(value, bucket_range)
                if bucket_index is not None:
                    buckets[bucket_index] += 1
                else:
                    error_values.append(value)
            except ValueError:
                error_values.append(row[0])
        
        output_file_path = os.path.join(output_dir, os.path.basename(file_path))
        with open(output_file_path, 'w') as out_f:
            for i in range(500):
                out_f.write(f"{buckets[i]}\n")
        
        return error_values

def main(directory):
    bucket_range, _ = create_buckets()
    processed_dir = os.path.join(directory, 'processed')
    os.makedirs(processed_dir, exist_ok=True)
    
    success_count = 0
    error_files = {}
    
    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            file_path = os.path.join(directory, filename)
            error_values = process_file(file_path, processed_dir, bucket_range)
            
            if error_values:
                error_files[filename] = error_values
            else:
                success_count += 1
    
    print(f"Successfully processed {success_count} files.")
    if error_files:
        print(f"Files with errors: {len(error_files)}")
        for fname, errors in error_files.items():
            print(f"{fname}: {errors}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CSV file parser and processor.')
    parser.add_argument('directory', type=str, help='Path to the directory containing CSV files.')
    args = parser.parse_args()
    main(args.directory)

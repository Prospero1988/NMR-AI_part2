import os
import sys
import glob
import pandas as pd
import hashlib

def hash_to_200_bits(input_str: str):
    """
    Converts a string into a 200-bit binary vector.
    Uses sha256 (256 bits) and truncates it to 200 bits.
    Returns a list of 200 elements (0/1).
    """
    # Compute the hash (256 bits = 32 bytes)
    hash_bytes = hashlib.sha256(input_str.encode('utf-8')).digest()
    # Convert to a large integer
    big_int = int.from_bytes(hash_bytes, 'big')  
    # Truncate to 200 bits: shift left by (256 - 200) = 56 bits
    truncated = big_int >> 56  
    # Create a binary string (without the "0b" prefix), padded to 200 bits
    bit_string = bin(truncated)[2:].zfill(200)
    # Convert to a list of [0,1]
    return [int(bit) for bit in bit_string]

def process_csv_file(csv_path: str):
    """
    Processes a single CSV file according to the requirements:
    1. Reads the CSV file.
    2. Moves the third column to the first position in the new file (keeping the column name).
    3. Generates 200-bit vectors from the MOLECULE_NAME column.
    4. Saves the new file with the suffix '_hashed_IDs'.
    """
    # Read the file
    df = pd.read_csv(csv_path)
    
    # Report column names and rows (optional)
    print(f"  -> Successfully read file: {csv_path}")
    print(f"     Columns: {list(df.columns)}")

    # Assuming MOLECULE_NAME is the 1st column (index 0),
    # and the third column (index 2) is the label to be moved to the start.
    third_col_name = df.columns[2]   # e.g., 'CHILOGD_7.4'
    molecule_col_name = df.columns[0]  # e.g., 'MOLECULE_NAME'

    # Prepare the new DataFrame
    new_df = pd.DataFrame()

    # 1) Move the third column to the first position
    new_df[third_col_name] = df.iloc[:, 2]

    # 2) Generate 200-bit vectors from the MOLECULE_NAME column
    bit_vectors = []
    for val in df[molecule_col_name].astype(str):
        bit_vec = hash_to_200_bits(val)
        bit_vectors.append(bit_vec)

    # Create a DataFrame with columns "bit_1" ... "bit_200"
    bit_col_names = [f"hash_bit_{i+1:03d}" for i in range(200)]
    bit_df = pd.DataFrame(bit_vectors, columns=bit_col_names)

    # 3) Add the 200-bit binary vector to new_df
    new_df = pd.concat([new_df, bit_df], axis=1)

    # Generate the output path
    base_name, ext = os.path.splitext(csv_path)
    output_path = f"{base_name}_hashed_IDs.csv"

    # Save the new file
    new_df.to_csv(output_path, index=False)

    print(f"  -> Saved new file: {output_path}")

def process_folder(folder_path: str):
    """
    Processes all CSV files in the given folder.
    """
    csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
    if not csv_files:
        print(f"No CSV files found in the folder: {folder_path}")
        return

    print(f"Found {len(csv_files)} CSV file(s) in the folder: {folder_path}")
    for csv_file in csv_files:
        try:
            process_csv_file(csv_file)
        except Exception as e:
            print(f"  [ERROR] Problem with file {csv_file}: {e}")
        else:
            print(f"  [OK] Successfully processed: {csv_file}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <path_to_folder>")
        sys.exit(1)
    
    folder_path = sys.argv[1]
    print(f"Starting processing of folder: {folder_path}")
    try:
        process_folder(folder_path)
    except Exception as e:
        print(f"  [ERROR] Failed to process the folder: {e}")
    else:
        print("Finished processing all files.")

if __name__ == "__main__":
    main()

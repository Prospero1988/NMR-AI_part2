import os
import csv
import shutil
import sys

def create_directory_if_not_exists(directory_path):
    """Creates a directory if it does not exist."""
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def search_and_copy_files(source_dir, target_dir, molecule_names):
    """Search for files in source_dir containing phrases from molecule_names and copy them to target_dir."""
    files_copied = 0
    found_molecule_names = set()

    for file_name in os.listdir(source_dir):
        for molecule_name in molecule_names:
            if molecule_name in file_name:
                source_file_path = os.path.join(source_dir, file_name)
                target_file_path = os.path.join(target_dir, file_name)
                shutil.copy2(source_file_path, target_file_path)
                files_copied += 1
                found_molecule_names.add(molecule_name)

    return files_copied, found_molecule_names

def process_csv_file(csv_file_path, source_dir):
    """Processes a CSV file, searches for files and copies them to appropriate folders."""
    folder_name = os.path.splitext(os.path.basename(csv_file_path))[0]
    target_dir = os.path.join(source_dir, folder_name)

    create_directory_if_not_exists(target_dir)

    molecule_names = []
    with open(csv_file_path, mode='r', encoding='utf-8') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            molecule_names.append(row['MOLECULE_NAME'])

    print(f"Processing CSV file: {csv_file_path}")
    print(f"Target folder: {target_dir}")
    print(f"Number of records in the 'MOLECULE_NAME' column: {len(molecule_names)}")

    files_copied, found_molecule_names = search_and_copy_files(source_dir, target_dir, molecule_names)

    print(f"Number of files found and copied: {files_copied}")

    # Determine which molecule names were not found
    not_found_molecule_names = set(molecule_names) - found_molecule_names
    if not_found_molecule_names:
        print(f"Molecule names not found during search: {', '.join(not_found_molecule_names)}")
    else:
        print("All molecule names were found in the search.")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <source_directory> <csv_directory>")
        sys.exit(1)

    source_dir = sys.argv[1]
    csv_dir = sys.argv[2]

    if not os.path.exists(source_dir) or not os.path.isdir(source_dir):
        print(f"Directory {source_dir} does not exist or is not a directory.")
        sys.exit(1)

    if not os.path.exists(csv_dir) or not os.path.isdir(csv_dir):
        print(f"Directory {csv_dir} does not exist or is not a directory.")
        sys.exit(1)

    csv_files = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]

    if not csv_files:
        print("No CSV files found in the directory.")
        sys.exit(1)

    for csv_file in csv_files:
        csv_file_path = os.path.join(csv_dir, csv_file)
        process_csv_file(csv_file_path, source_dir)

if __name__ == "__main__":
    main()

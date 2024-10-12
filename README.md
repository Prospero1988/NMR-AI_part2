
# NMR Data Processing and Machine Learning Scripts

This repository contains several Python scripts developed for processing Nuclear Magnetic Resonance (NMR) data, managing and manipulating CSV files, and performing machine learning (ML) regression tasks. Below is a detailed overview of each script, its purpose, usage instructions, and additional relevant information.

## Input data availability

This repository also contains input data files that are essential for the proper functioning of several scripts. 
The archive **Input_data.rar** contains all the CSV files with input data for ML computations, 
including reduced spectra combined with corresponding logD parameters in the rows of the feature space matrix and labels.
These input files are crucial for running the calculations and must be placed in the correct directories for successful execution.

## Table of Contents
- [add_12_zero_columns_at_the_end.py](#add_12_zero_columns_at_the_endpy)
- [bucket_integration_dir.py](#bucket_integration_dirpy)
- [columns_compare_and_copy.py](#columns_compare_and_copypy)
- [create_file_list.py](#create_file_listpy)
- [delete_empty_rows.py](#delete_empty_rowspy)
- [extract_first_column_from_csv.py](#extract_first_column_from_csvpy)
- [header_and_column_insert_dir.py](#header_and_column_insert_dirpy)
- [input_percentage_mixer_verbose_list.py](#input_percentage_mixer_verbose_listpy)
- [MOL_REP_generation.py](#mol_rep_generationpy)
- [nmrshift_db_spectra_bucket_creator.py](#nmrshift_db_spectra_bucket_creatorpy)
- [preparation_total.py](#preparation_totalpy)
- [Regressor_10CV_query_model.py](#regressor_10cv_query_modelpy)
- [Regressor_10CV_with_model_saving.py](#regressor_10cv_with_model_savingpy)
- [search_files_by_list_csv.py](#search_files_by_list_csvpy)
- [smiles_to_mol.py](#smiles_to_molpy)
- [total_reduction_from_list.py](#total_reduction_from_listpy)

---

## Script Descriptions

### add_12_zero_columns_at_the_end.py

**Purpose**: 
This script adds 12 additional columns at the end of a CSV file and fills them with zeros. It is useful when padding datasets with extra information or preparing a dataset for certain algorithms that require a fixed number of features.

**How it works**:
- The script iterates over all CSV files in a directory and appends 12 zero-filled columns to each file.
- The new columns are named C1 to C12.
  
**How to use**:
```bash
python add_12_zero_columns_at_the_end.py --dir DIRECTORY_PATH
```

Example:
```bash
python add_12_zero_columns_at_the_end.py --dir ./data
```

---

### bucket_integration_dir.py

**Purpose**: 
This script performs bucket integration on NMR datasets. It splits the data into a specified number of ranges (buckets) and calculates the sum of values in each range, reducing the dimensionality of the dataset.

**Key Features**:
- Accepts CSV files with NMR data stored as intensity values per spectral point.
- Allows the user to define how many "buckets" to create for integration.
- Outputs a new CSV file with integrated data.

**How to use**:
```bash
python bucket_integration_dir.py --dir DIRECTORY_PATH --ranges N
```
- **`--dir`**: Specify the directory containing the CSV files.
- **`--ranges`**: Number of ranges (buckets) to divide the data into.

---

### columns_compare_and_copy.py

**Purpose**: 
This script copies rows from one CSV file to another if the values in the first column of both files match. This is particularly useful for filtering datasets based on a reference list.

**Key Features**:
- Compares the first column of an input CSV file with a reference CSV file.
- Copies matching rows into a new CSV file.

**How to use**:
```bash
python columns_compare_and_copy.py
```
The script prompts the user to provide input files (the source CSV file, the reference CSV file, and the output file).

---

### create_file_list.py

**Purpose**: 
Scans a directory and its subdirectories, then creates a CSV file listing all file paths.

**Key Features**:
- Recursively scans a directory and its subdirectories.
- Outputs a CSV file containing the full paths of all files.

**How to use**:
```bash
python create_file_list.py
```
The script prompts the user to input the directory path and output file name.

---

### delete_empty_rows.py

**Purpose**: 
This script removes rows from a CSV file if any of the cells in the row are empty.

**Key Features**:
- Cleans a CSV file by removing rows that contain empty cells.
- Overwrites the original file with the cleaned data.

**How to use**:
```bash
python delete_empty_rows.py
```
The script prompts for the name of the CSV file to clean.

---

### extract_first_column_from_csv.py

**Purpose**: 
Extracts the first column from a CSV file and saves it to a new CSV file.

**Key Features**:
- Saves the first column of a given CSV file into a new file.
- Can be used for creating datasets with only labels or identifiers.

**How to use**:
```bash
python extract_first_column_from_csv.py
```
The script prompts the user for the input CSV file and output file paths.

---

### header_and_column_insert_dir.py

**Purpose**: 
Manipulates CSV files by copying columns from one file to another and inserting headers.

**Key Features**:
- Inserts a column from a source CSV file into another file as the first column.
- Allows the user to add custom headers to the columns.

**How to use**:
```bash
python header_and_column_insert_dir.py --dir DIRECTORY_PATH --headers CUSTOM_HEADERS --second-file SECOND_FILE_PATH --column-number COLUMN_NUMBER
```
The script processes the directory containing CSV files, adding headers and columns based on user input.

---

### input_percentage_mixer_verbose_list.py

**Purpose**: 
Creates mixed datasets by randomly combining rows from two CSV files. This is useful for creating training datasets with a mix of experimental and theoretical data.

**Key Features**:
- Mixes rows between two CSV files based on user-defined percentages (e.g., 20%, 40%, 60%, 80%).
- Outputs multiple files with different row combinations.

**How to use**:
```bash
python input_percentage_mixer_verbose_list.py <csv_file_with_paths>
```
- `<csv_file_with_paths>`: A CSV file containing pairs of file paths to mix.

---

### MOL_REP_generation.py

**Purpose**: 
Generates molecular fingerprints and descriptors from SMILES strings stored in a CSV file. These fingerprints can be used for cheminformatics tasks or machine learning models predicting molecular properties.

**Key Features**:
- Generates multiple types of molecular fingerprints (e.g., RDKit, ECFP4, MACCS).
- Outputs separate CSV files for each fingerprint type.

**How to use**:
```bash
python MOL_REP_generation.py input_file.csv
```
The script processes the input CSV file and generates fingerprint files for each molecular structure.

---

### nmrshift_db_spectra_bucket_creator.py

**Purpose**: 
Generates pseudo-NMR spectra from chemical shift data and integrates the results into buckets. It is particularly useful for datasets from DFT calculations or NMR predictor databases.

**Key Features**:
- Reads CSV files containing chemical shift data and applies bucket integration.
- Outputs pseudo-NMR spectra in bucketed form.

**Usage Details**: 
Refer to the `nmrshift_db_spectra_bucket_creator_info.docx` file for more detailed instructions.

---

### preparation_total.py

**Purpose**: 
Processes NMR datasets stored in CSV files, allowing for normalization and selective removal of columns. The script also supports visualization of the data before and after processing.

**Special Feature**: 
The script integrates with `List_for_point_cutting_solvents.txt` to define ranges for removing columns from the dataset during preprocessing.

**How to use**:
```bash
python preparation_total.py
```
The script guides the user through the process with prompts for input, parameters, and optional visualization.

---

### Regressor_10CV_query_model.py

**Purpose**: 
Applies a pre-trained machine learning model to new datasets and evaluates the predictions. This script is particularly useful for validating regression models on new data.

**Key Features**:
- Loads a pre-trained model and applies it to new input data.
- Calculates metrics such as RMSE, MAE, and R² for each prediction.

**How to use**:
```bash
python Regressor_10CV_query_model.py <csv_file> <model_file>
```

---

### Regressor_10CV_with_model_saving.py

**Purpose**: 
Performs regression analysis using Support Vector Regression (SVR), AdaBoost, and Gradient Boosting with 10-fold cross-validation. This script saves the trained model for future use.

**Key Features**:
- Trains and validates regression models using cross-validation.
- Saves the trained models for future inference.

**How to use**:
```bash
python Regressor_10CV_with_model_saving.py <directory_path_to_CSV_files>
```
The script processes all CSV files in the specified directory and logs performance metrics.

---

### search_files_by_list_csv.py

**Purpose**: 
Searches for files in a directory based on lists from CSV files and copies the matching files into corresponding directories.

**Key Features**:
- Uses CSV files as a reference for searching and copying files.
- Automatically organizes found files into directories based on the CSV file names.

**How to use**:
```bash
python search_files_by_list_csv.py <katalog_1> <katalog_2>
```

---

### smiles_to_mol.py

**Purpose**: 
Converts SMILES strings from a CSV file into `.mol` files containing 3D molecular structures.

**How to use**:
```bash
python smiles_to_mol.py molecules.csv
```

---

### total_reduction_from_list.py

**Purpose**: 
Performs dimensionality reduction on high-dimensional datasets using various techniques (PCA, ICA, NMF, UMAP). The script allows you to specify target dimensionality and the number of decimal places to round values to.

**How to use**:
```bash
python total_reduction_from_list.py
```

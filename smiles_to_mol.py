# -*- coding: utf-8 -*-

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Function to generate .mol files from a CSV containing molecule names and SMILES
def generate_mol_files(csv_path):
    try:
        # Load the CSV file
        data = pd.read_csv(csv_path)
        
        # Check if the necessary columns exist in the CSV file
        if 'MOLECULE_NAME' not in data.columns or 'SMILES' not in data.columns:
            raise ValueError("CSV file must contain 'MOLECULE_NAME' and 'SMILES' columns.")
        
        file_count = 0  # Counter to track the number of files generated
        
        # Iterate through each row in the CSV file
        for index, row in data.iterrows():
            try:
                name = row['MOLECULE_NAME']  # Get the molecule name
                smiles = row['SMILES']  # Get the SMILES string
                
                # Generate the molecule from the SMILES string
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError(f"Invalid SMILES string: {smiles}")
                
                # Add hydrogens to the molecule
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates for the molecule
                AllChem.EmbedMolecule(mol)
                
                # Optimize the molecule's 2D coordinates
                AllChem.Compute2DCoords(mol)
                
                # Save the molecule structure to a .mol file
                mol_file = f"{name}.mol"  # Use the molecule name for the .mol file
                with open(mol_file, 'w') as f:
                    f.write(Chem.MolToMolBlock(mol))  # Write the molecule block to the file
                
                file_count += 1  # Increment the file counter after successful save
                
            except Exception as e:
                print(f"Error processing row {index}: {e}")
        
        print(f"Generated {file_count} .mol files.")  # Output the number of files generated
    
    except FileNotFoundError:
        print(f"File not found: {csv_path}")  # Handle case when the CSV file is missing
    except pd.errors.EmptyDataError:
        print(f"File is empty: {csv_path}")  # Handle case when the CSV file is empty
    except pd.errors.ParserError:
        print(f"Error parsing CSV file: {csv_path}")  # Handle case when there is a CSV parsing issue
    except Exception as e:
        print(f"An unexpected error occurred: {e}")  # Catch any other unexpected errors

# Main script execution
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_csv>")  # Inform user how to run the script
    else:
        csv_path = sys.argv[1]  # Get the CSV file path from command-line arguments
        generate_mol_files(csv_path)  # Call the function to generate .mol files

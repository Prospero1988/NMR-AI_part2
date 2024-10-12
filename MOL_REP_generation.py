# -*- coding: utf-8 -*-

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import rdMolDescriptors
import numpy as np

# Function to generate specific types of molecular fingerprints
def generate_fingerprint(smiles, fp_type):
    mol = Chem.MolFromSmiles(smiles)
    if fp_type == "RDKit":
        # Generate RDKit fingerprint (Morgan with radius 0)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 0, nBits=2048)
    elif fp_type == "ECFP4":
        # Generate ECFP4 fingerprint (Morgan with radius 2)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    elif fp_type == "MACCS":
        # Generate MACCS fingerprint
        fp = MACCSkeys.GenMACCSKeys(mol)
    elif fp_type == "Klekota-Roth":
        # Generate Klekota-Roth (Atom Pair) fingerprint
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=2048)
    return list(fp)

# Function to generate and save fingerprints
def save_fingerprints(df, fp_type):
    # Apply the fingerprint generation to the 'SMILES' column
    fp_data = df['SMILES'].apply(lambda x: generate_fingerprint(x, fp_type))
    fp_df = pd.DataFrame(fp_data.tolist(), index=df.index)
    # Concatenate the original ID with the fingerprints
    output_df = pd.concat([df.iloc[:, 0], fp_df], axis=1)
    output_file = f"{fp_type}_fingerprints.csv"
    output_df.to_csv(output_file, index=False)
    print(f"Saved {fp_type} fingerprints to file {output_file}")

# Function to generate molecular descriptors
def generate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    # Calculate all available descriptors
    descriptors = [desc(mol) for _, desc in Descriptors._descList]
    return descriptors

# Function to save molecular descriptors
def save_descriptors(df):
    # Apply descriptor generation to the 'SMILES' column
    descriptors_data = df['SMILES'].apply(generate_descriptors)
    descriptors_df = pd.DataFrame(descriptors_data.tolist(), index=df.index)
    # Normalize descriptor values
    descriptors_df = (descriptors_df - descriptors_df.min()) / (descriptors_df.max() - descriptors_df.min())
    # Concatenate the original ID with the descriptors
    output_df = pd.concat([df.iloc[:, 0], descriptors_df], axis=1)
    output_file = "Molecular_Descriptors.csv"
    output_df.to_csv(output_file, index=False)
    print(f"Saved molecular descriptors to file {output_file}")

# Main script logic
if __name__ == "__main__":
    # Input CSV file from command-line argument
    input_file = sys.argv[1]
    df = pd.read_csv(input_file)

    # Generate and save fingerprints for each type
    for fp_type in ["RDKit", "ECFP4", "MACCS", "Klekota-Roth"]:
        save_fingerprints(df, fp_type)
    
    # Generate and save molecular descriptors
    save_descriptors(df)

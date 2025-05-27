import os
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# Define input and output directories
input_dir = "reference_adata"
output_dir = "reference_adata"


# Function to process each file
def process_adata(file_path):
    print(f"Processing {file_path.name}...")
    
    # Load the AnnData object
    adata = sc.read_h5ad(file_path)
    
    # Print basic information
    print(f"Original shape: {adata.shape}")
    
    # Filter genes with min count of 2
    sc.pp.filter_genes(adata, min_counts=2)
    
    # Filter cells with min count of 2
    sc.pp.filter_cells(adata, min_counts=2)
    
    # Remove mitochondrial genes
    adata = adata[:, [g for g in adata.var_names if not g.startswith(('mt-', 'MT-', 'Mt-'))]]
    
    print(f"Shape after filtering: {adata.shape}")
    
    # Show available observations for renaming
    print("\nAvailable observations (columns in .obs):")
    for i, col in enumerate(adata.obs.columns):
        print(f"{i+1}. {col}")
    
    # Ask which observation to rename to cellType
    rename_choice = input("\nEnter the number of the observation to rename as 'cellType' (or 'skip' to proceed without renaming): ")
    
    if rename_choice.lower() != 'skip':
        try:
            col_idx = int(rename_choice) - 1
            if 0 <= col_idx < len(adata.obs.columns):
                col_name = adata.obs.columns[col_idx]
                print(f"Renaming '{col_name}' to 'cellType'")
                
                # Create a copy of the selected column as 'cellType'
                adata.obs['cellType'] = adata.obs[col_name].copy()
                
                print(f"Unique cell types: {adata.obs['cellType'].unique()}")
            else:
                print("Invalid selection. No renaming performed.")
        except ValueError:
            print("Invalid input. No renaming performed.")
    
    # Save the processed file
    output_path = os.path.join(output_dir, file_path.name)
    adata.write(output_path)
    print(f"Saved to {output_path}\n")
    
    return adata

# List all h5ad files in the input directory
h5ad_files = list(Path(input_dir).glob("*.h5ad"))

# Process files one by one
for file_path in h5ad_files:
    process_or_skip = input(f"Process {file_path.name}? (y/n): ")
    if process_or_skip.lower() == 'y':
        process_adata(file_path)
    else:
        print(f"Skipping {file_path.name}\n")

print("Processing complete!")
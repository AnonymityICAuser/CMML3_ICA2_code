import scipy.sparse as sparse
import scipy.io as sio
import numpy as np
import scanpy as sc
import os
import pandas as pd

def save_data_to_csv_and_mtx(h5ad_file_path, output_dir):
    # Check if the output directory exists, if not create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load the .h5ad file
    all_data = sc.read_h5ad(h5ad_file_path)

    # filter cell and gene
    sc.pp.filter_cells(all_data, min_genes=1)
    sc.pp.filter_genes(all_data, min_cells=1)
    
    # Extract cell info, gene info, and count matrix
    cellinfo = all_data.obs
    geneinfo = all_data.var
    mtx = all_data.X.T  # Transpose to get genes as rows, cells as columns

    # set cellID as index for cellinfo (currently not unified)
    cellinfo.index = all_data.obs_names
    # reset name as cellID
    cellinfo.index.name = 'cell_ID'

    geneinfo.index = all_data.var_names
    geneinfo.index.name = 'gene_identifier'


    # Define the file paths for the output files
    cellinfo_csv_path = os.path.join(output_dir, "cellinfo.csv")
    geneinfo_csv_path = os.path.join(output_dir, "geneinfo.csv")
    mtx_mm_path = os.path.join(output_dir, "sparse_matrix.mtx")

    # Handle different types of data (spatial or single-cell)
    if 'spatial' in all_data.obsm:
        # Save spatial information if available
        spatial = all_data.obsm['spatial']
        spatial_df = pd.DataFrame(spatial, columns=['x', 'y'])
        spatial_df.to_csv(os.path.join(output_dir, "spatial.csv"))
        print(f"Spatial information saved for {os.path.basename(h5ad_file_path)}")
    else:
        print(f"No spatial information found in {os.path.basename(h5ad_file_path)} - treating as scRNA-seq data")
        

    # Export to CSV
    cellinfo.to_csv(cellinfo_csv_path)
    geneinfo.to_csv(geneinfo_csv_path)
    
    # Make sure the matrix is in csr format
    if not sparse.issparse(mtx):
        mtx = sparse.csr_matrix(mtx)
    else:
        mtx = mtx.tocsr()
        
    # Export count matrix in Matrix Market format
    sio.mmwrite(mtx_mm_path, mtx)

    return all_data

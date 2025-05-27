import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import cell2location
from cell2location.models import RegressionModel
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
import scipy.sparse as sparse
import torch # Keep torch for checking availability
import pandas as pd
import gc
# Determine accelerator and devices for scvi-tools/PyTorch Lightning
if torch.cuda.is_available():
    accelerator = "cuda"
    devices = 1 # Use 1 GPU (typically the first one, ID 0)
    # If you want to specify a particular GPU or multiple, e.g., [0] or [0, 1]
    # devices = [0] # To use GPU 0
else:
    accelerator = "cpu"
    devices = "auto" # or 1 for one CPU core, "auto" is usually fine

def run_cell2loc(adata_ref,adata_decon):
    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
    labels_key='cellType',
    )

    mod = RegressionModel(adata_ref)

    # Pass accelerator and devices arguments
    mod.train(accelerator=accelerator) # Default max_epochs is 250

    adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
    )

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    if sparse.issparse(adata_decon.X):
        # Ensure adata_ref.X is also sparse if adata_decon.X is sparse, and convert to CSR
        if not sparse.issparse(adata_ref.X) or not sparse.isspmatrix_csr(adata_ref.X):
             adata_ref.X = sparse.csr_matrix(adata_ref.X) # Convert to sparse CSR
        if not sparse.isspmatrix_csr(adata_decon.X):
            adata_decon.X = adata_decon.X.tocsr()
    else: # If adata_decon.X is dense, ensure adata_ref.X is also treated as dense then converted
        adata_ref.X = sparse.csr_matrix(adata_ref.X) # Convert to sparse CSR, assumes it can be dense or sparse
        adata_decon.X = sparse.csr_matrix(adata_decon.X) # Convert to sparse CSR

    common_genes = adata_decon.var_names.intersection(inf_aver.index)
    if len(common_genes) == 0:
        raise ValueError("No common genes found between reference signatures (inf_aver) "
                         "and spatial data (adata_decon). Ensure gene names are consistent.")
    print(f"Found {len(common_genes)} common genes. "
          f"Original genes in reference: {len(inf_aver.index)}, "
          f"Original genes in spatial data: {len(adata_decon.var_names)}.")
    adata_decon = adata_decon[:, common_genes].copy()
    inf_aver = inf_aver.loc[common_genes, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata=adata_decon)
    mod = cell2location.models.Cell2location(
    adata_decon, cell_state_df=inf_aver,
    N_cells_per_location=10,
    detection_alpha=20
    )

    # Pass accelerator and devices arguments
    mod.train(max_epochs=20000, accelerator=accelerator)

    adata_decon = mod.export_posterior(
    adata_decon, sample_kwargs={'num_samples': 1000, 'batch_size': adata_decon.shape[0]}
    )
    if 'q05_cell_abundance_w_sf' not in adata_decon.obsm:
        raise KeyError("Key 'q05_cell_abundance_w_sf' not found in adata_decon.obsm. ")
    if 'mod' not in adata_decon.uns or 'factor_names' not in adata_decon.uns['mod']:
        raise KeyError("Key 'factor_names' not found in adata_decon.uns['mod']. ")
    adata_decon.obs[adata_decon.uns['mod']['factor_names']] = adata_decon.obsm['q05_cell_abundance_w_sf']
    return get_cell_abundance_from_obs(adata_decon)

def get_cell_abundance_from_obs(adata_decon):
    potential_metadata_cols = ['n_counts', '_indices', '_scvi_batch', '_scvi_labels', 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'total_counts']

    if 'mod' in adata_decon.uns and 'factor_names' in adata_decon.uns['mod']:
        cell_type_cols = list(adata_decon.uns['mod']['factor_names'])
        cell_type_cols = [col for col in cell_type_cols if col in adata_decon.obs.columns]
    else:
        cell_type_cols = [col for col in adata_decon.obs.columns if col not in potential_metadata_cols and not col.startswith('_')]

    if not cell_type_cols:
        print("Warning: No cell type abundance columns identified in `adata_decon.obs`. ")
        return pd.DataFrame(index=adata_decon.obs.index)
    abundance_df = adata_decon.obs[cell_type_cols].copy()

    # clean up cuda memory
    torch.cuda.empty_cache()
    gc.collect()
    return abundance_df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run cell2location deconvolution')
    parser.add_argument('--adata_ref', type=str, required=True, help='Path to reference h5ad file')
    parser.add_argument('--adata_decon', type=str, required=True, help='Path to deconvolution h5ad file')
    parser.add_argument('--output_csv', type=str, required=True, help='Path to output file')
    args = parser.parse_args()
    adata_ref = sc.read_h5ad(args.adata_ref)
    adata_decon = sc.read_h5ad(args.adata_decon)
    if not adata_ref.var_names.is_unique:
        print("Warning: Reference AnnData var_names are not unique. Making them unique.")
        adata_ref.var_names_make_unique()
    if not adata_decon.var_names.is_unique:
        print("Warning: Deconvolution AnnData var_names are not unique. Making them unique.")
        adata_decon.var_names_make_unique()
    df = run_cell2loc(adata_ref, adata_decon)
    df.to_csv(args.output_csv)
    print(f"Cell abundance table saved to {args.output_csv}")
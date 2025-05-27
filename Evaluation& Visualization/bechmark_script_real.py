import numpy as np
import pandas as pd
import os
import scipy.stats as stats
from scipy.spatial.distance import jensenshannon
import scanpy as sc
import scipy.sparse
import anndata # Explicitly import for type hinting
from typing import Union, Dict, List, Tuple, Optional # Optional is also good for `X | None`
# --- I. Helper Functions ---
def normalize_proportions_df(proportions_df: pd.DataFrame) -> pd.DataFrame:
    """Row-normalizes a DataFrame of proportions to sum to 1 per row (spot)."""
    if proportions_df.empty: return proportions_df.copy() # Return a copy if empty
    df_norm = proportions_df.copy()
    row_sums = df_norm.sum(axis=1)
    non_zero_rows = row_sums > 1e-9 # Avoid division by zero for very small sums
    df_norm.loc[non_zero_rows] = df_norm.loc[non_zero_rows].div(row_sums[non_zero_rows], axis=0)
    # For rows that summed to zero (or near zero), their proportions will remain zero or become NaN if all were zero.
    # If NaNs are an issue later, consider df_norm.fillna(0, inplace=True) after division.
    return df_norm

# --- II. Indicator Calculation Functions ---

def calculate_gt_accuracy_metrics(predicted_norm_df: pd.DataFrame, ground_truth_df: pd.DataFrame) -> dict:
    """Calculates RMSE, JSD, and per-cell-type PCC against ground truth."""
    overall_metrics = {
        'Overall_RMSE_vs_GT': np.nan, 'median_JSD_vs_GT': np.nan,
        'median_PCC_per_cell_type_vs_GT': np.nan, 'mean_PCC_per_cell_type_vs_GT': np.nan
    }
    by_cell_type_metrics = pd.DataFrame()
    by_spot_metrics = pd.DataFrame() # Retained for potential future use, not currently returned to main pipeline
    
    if predicted_norm_df.empty or ground_truth_df.empty:
        return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics, 'by_spot': by_spot_metrics}

    common_spots = predicted_norm_df.index.intersection(ground_truth_df.index)
    if not common_spots.any(): return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics, 'by_spot': by_spot_metrics}
    
    pred_aligned = predicted_norm_df.loc[common_spots]
    gt_aligned = ground_truth_df.loc[common_spots]
    
    common_cts = pred_aligned.columns.intersection(gt_aligned.columns)
    if not common_cts.any(): return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics, 'by_spot': by_spot_metrics}

    pred_aligned = pred_aligned[common_cts]
    gt_aligned = gt_aligned[common_cts] # gt_aligned proportions are used for RMSE and PCC directly
    
    pred_matrix = pred_aligned.values
    gt_matrix = gt_aligned.values

    overall_metrics['Overall_RMSE_vs_GT'] = np.sqrt(np.mean((pred_matrix - gt_matrix)**2))

    # For JSD, ensure both predicted and ground truth are valid probability distributions (row-normalized)
    # pred_matrix is already from predicted_norm_df, so it's row-normalized.
    # gt_aligned might or might not be normalized from source, ensure it is for JSD.
    gt_aligned_norm_for_jsd = normalize_proportions_df(gt_aligned.copy()) 
    gt_matrix_norm_for_jsd = gt_aligned_norm_for_jsd.values
    
    jsd_values = np.zeros(len(common_spots))
    for i in range(pred_matrix.shape[0]):
        p_spot = pred_matrix[i, :]
        q_spot = gt_matrix_norm_for_jsd[i, :]

        sum_p = np.sum(p_spot)
        sum_q = np.sum(q_spot)

        if sum_p < 1e-9 and sum_q < 1e-9: 
            jsd_values[i] = 0.0
        elif sum_p < 1e-9 or sum_q < 1e-9: 
            jsd_values[i] = 1.0 # Max JSD (base 2 log(2)=1) if one dist is zero and other is not
        else:
            jsd_values[i] = jensenshannon(p_spot, q_spot, base=2)
            
    by_spot_metrics['JSD_vs_GT'] = jsd_values # For potential detailed by-spot output
    overall_metrics['median_JSD_vs_GT'] = np.nanmedian(jsd_values)

    pcc_ct = pd.Series(index=common_cts, dtype=float, name="PCC_vs_GT")
    rmse_ct = pd.Series(index=common_cts, dtype=float, name="RMSE_vs_GT")
    for i, ct in enumerate(common_cts): 
        p_vec, t_vec = pred_matrix[:, i], gt_matrix[:, i] # Use original (potentially unnormalized for GT) proportions for PCC/RMSE
        rmse_ct[ct] = np.sqrt(np.mean((p_vec - t_vec)**2))
        if len(np.unique(p_vec)) > 1 and len(np.unique(t_vec)) > 1:
            pcc_ct[ct], _ = stats.pearsonr(p_vec, t_vec)
        else: pcc_ct[ct] = np.nan
    by_cell_type_metrics['PCC_vs_GT'] = pcc_ct
    by_cell_type_metrics['RMSE_vs_GT'] = rmse_ct
    overall_metrics['median_PCC_per_cell_type_vs_GT'] = pcc_ct.median()
    overall_metrics['mean_PCC_per_cell_type_vs_GT'] = pcc_ct.mean()
    
    return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics, 'by_spot': by_spot_metrics}

def calculate_intrinsic_metrics(predicted_norm_df: pd.DataFrame) -> dict:
    overall_metrics = {'median_Entropy_Predicted': np.nan, 'mean_Num_Detected_CellTypes_Predicted': np.nan}
    by_spot_metrics = pd.DataFrame(index=predicted_norm_df.index if not predicted_norm_df.empty else None)
    if predicted_norm_df.empty: return {'overall': overall_metrics, 'by_spot': by_spot_metrics}
    pred_matrix = predicted_norm_df.values

    entropy_values = np.zeros(pred_matrix.shape[0])
    for i in range(pred_matrix.shape[0]):
        p_spot = pred_matrix[i, :]
        p_spot_nz = p_spot[p_spot > 1e-9] # Ensure positive for entropy
        if len(p_spot_nz) > 0: entropy_values[i] = stats.entropy(p_spot_nz, base=2)
        else: entropy_values[i] = 0.0 # Entropy of a zero vector or single point is 0
    by_spot_metrics['Entropy_Predicted'] = entropy_values
    overall_metrics['median_Entropy_Predicted'] = np.nanmedian(entropy_values)

    detection_threshold = 0.01 # Cell type considered "detected" if proportion > this value
    num_detected = np.sum(pred_matrix > detection_threshold, axis=1)
    by_spot_metrics['Num_Detected_CellTypes_Predicted'] = num_detected
    overall_metrics['mean_Num_Detected_CellTypes_Predicted'] = np.mean(num_detected)
    return {'overall': overall_metrics, 'by_spot': by_spot_metrics}

def calculate_marker_correlation_metric(
    predicted_norm_df: pd.DataFrame, adata_spatial_expression: anndata.AnnData, 
    marker_gene_map: dict, expression_layer: str = None, min_spots_expressed: int = 5
    ) -> dict:
    overall_metrics = {'Mean_Marker_PCC_all_cts': np.nan, 'Median_Marker_PCC_all_cts': np.nan}
    by_cell_type_metrics = pd.DataFrame(columns=['Marker_PCC']) # Initialize to ensure it exists
    if adata_spatial_expression is None or not marker_gene_map or predicted_norm_df.empty: 
        return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics}

    common_spots = predicted_norm_df.index.intersection(adata_spatial_expression.obs_names)
    if not common_spots.any(): return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics}
        
    proportions_aligned = predicted_norm_df.loc[common_spots]
    adata_aligned = adata_spatial_expression[common_spots, :].copy()

    marker_pccs_ct = pd.Series(dtype=float, name="Marker_PCC") 
    for cell_type, markers in marker_gene_map.items():
        if cell_type not in proportions_aligned.columns: continue
        prop_vector = proportions_aligned[cell_type].values
        ct_marker_corrs = []
        for marker in markers:
            if marker not in adata_aligned.var_names: continue
            
            if expression_layer and expression_layer in adata_aligned.layers:
                expr_data = adata_aligned[:, marker].layers[expression_layer]
            else:
                expr_data = adata_aligned[:, marker].X

            marker_expr_vector = expr_data.toarray().flatten() if scipy.sparse.issparse(expr_data) else np.array(expr_data).flatten()
            
            if np.sum(marker_expr_vector > 1e-6) < min_spots_expressed: continue # Ensure marker is expressed in min spots
            if len(np.unique(prop_vector)) > 1 and len(np.unique(marker_expr_vector)) > 1: # Check for variance
                corr, _ = stats.pearsonr(prop_vector, marker_expr_vector)
                if not np.isnan(corr): ct_marker_corrs.append(corr)
        marker_pccs_ct[cell_type] = np.nanmean(ct_marker_corrs) if ct_marker_corrs else np.nan # Use nanmean
    
    by_cell_type_metrics['Marker_PCC'] = marker_pccs_ct
    overall_metrics['Mean_Marker_PCC_all_cts'] = marker_pccs_ct.mean() # Mean of per-celltype average marker PCCs
    overall_metrics['Median_Marker_PCC_all_cts'] = marker_pccs_ct.median()
    return {'overall': overall_metrics, 'by_cell_type': by_cell_type_metrics}

def generate_deg_marker_map(
    adata_ref: anndata.AnnData, cell_type_obs_key: str, 
    n_top_genes: int = 20, positive_logfc_threshold: float = 0.25
    ) -> dict:
    print(f"  Generating top {n_top_genes} DEGs using '{cell_type_obs_key}' from AnnData with shape {adata_ref.shape}...")
    adata_temp = adata_ref.copy()
    if cell_type_obs_key not in adata_temp.obs:
        print(f"  Error: Cell type key '{cell_type_obs_key}' not found in reference AnnData .obs. Cannot generate DEGs.")
        return {}
        
    if not pd.api.types.is_categorical_dtype(adata_temp.obs[cell_type_obs_key]):
        adata_temp.obs[cell_type_obs_key] = adata_temp.obs[cell_type_obs_key].astype('category')
    
    min_cells_per_group = 3 
    type_counts = adata_temp.obs[cell_type_obs_key].value_counts()
    valid_groups_for_deg = type_counts[type_counts >= min_cells_per_group].index.tolist()
    
    if len(valid_groups_for_deg) < 2 : 
        print(f"  Warning: Not enough groups ({len(valid_groups_for_deg)}) with sufficient cells (>= {min_cells_per_group}) for DEG analysis. Returning empty marker map.")
        return {}
        
    adata_temp = adata_temp[adata_temp.obs[cell_type_obs_key].isin(valid_groups_for_deg)].copy()
    adata_temp.obs[cell_type_obs_key] = adata_temp.obs[cell_type_obs_key].astype('category') # Re-categorize after filtering

    # Request more genes initially to have enough after filtering
    sc.tl.rank_genes_groups(adata_temp, groupby=cell_type_obs_key, method='wilcoxon', 
                            use_raw=False, n_genes=max(n_top_genes * 5, 100), corr_method='benjamini-hochberg') 
    marker_map = {}
    if 'rank_genes_groups' not in adata_temp.uns: 
        print("  Warning: 'rank_genes_groups' not found in .uns after sc.tl.rank_genes_groups. DEG analysis failed.")
        return marker_map
        
    groups = adata_temp.obs[cell_type_obs_key].cat.categories.tolist() 
    for group in groups:
        try:
            deg_results = sc.get.rank_genes_groups_df(adata_temp, group=group)
            positive_degs = deg_results[
                (deg_results['logfoldchanges'] > positive_logfc_threshold) & 
                (deg_results['pvals_adj'] < 0.05) # Example significance filter
            ]
            top_degs_for_group = positive_degs['names'].head(n_top_genes).tolist()
            
            if top_degs_for_group: 
                marker_map[group] = top_degs_for_group
            else: # Fallback if no significant positive DEGs
                fallback_degs = deg_results[deg_results['logfoldchanges'] > 0] # Just positive logFC
                marker_map[group] = fallback_degs['names'].head(n_top_genes).tolist() if not fallback_degs.empty else []
                if not marker_map[group]:
                    print(f"    No markers found for '{group}' even with fallback (logFC > {positive_logfc_threshold}, p_adj < 0.05 or logFC > 0).")
        except KeyError: 
            print(f"    Warning: Cell type group '{group}' not found in DEG results (should not happen if using .cat.categories). Skipping.")
            marker_map[group] = []
    return marker_map

def calculate_aggregated_score(summary_df: pd.DataFrame, metrics_for_as: list, higher_is_better_map: dict) -> pd.DataFrame:
    as_df = summary_df.copy()
    rank_score_cols_for_mean = []
    
    # Determine N (number of methods being ranked) from the input summary_df
    # This assumes summary_df is for a single dataset, with methods either in index or 'Method' column
    if isinstance(summary_df.index, pd.MultiIndex) and 'Method' in summary_df.index.names:
         # If methods are part of a MultiIndex (e.g., before reset_index)
        num_methods_total = len(summary_df.index.get_level_values('Method').unique())
    elif 'Method' in summary_df.columns:
        num_methods_total = summary_df['Method'].nunique()
    else: # Assume index contains method names
        num_methods_total = summary_df.index.nunique()

    if num_methods_total == 0 and not as_df.empty : # Should not happen if df not empty
        num_methods_total = len(as_df)


    for metric_col in metrics_for_as:
        rank_col_name = f'Rank_Score_{metric_col}'
        if metric_col not in as_df.columns or as_df[metric_col].isnull().all():
            as_df[rank_col_name] = np.nan 
            continue
        
        is_higher_better = higher_is_better_map.get(metric_col, True)
        
        # Step 1: Get ranks where 1 is "best" (e.g., lowest RMSE is rank 1, highest PCC is rank 1)
        if is_higher_better:
            ranks_1_is_best = as_df[metric_col].rank(method='min', ascending=False, na_option='bottom')
        else:
            ranks_1_is_best = as_df[metric_col].rank(method='min', ascending=True, na_option='bottom')

        # Step 2: Convert these ranks to scores where higher score is better.
        # Score = N_valid - rank_1_is_best + 1. Best gets N_valid points, worst gets 1 point.
        # N_valid is the number of non-NaN values for this specific metric column.
        num_valid_for_metric = ranks_1_is_best.notna().sum()
        
        if num_valid_for_metric > 0:
            final_rank_scores = num_valid_for_metric - ranks_1_is_best + 1
            # Ensure that original NaNs in ranks_1_is_best result in NaN scores
            final_rank_scores[ranks_1_is_best.isna()] = np.nan 
        else: # All values for this metric were NaN
            final_rank_scores = pd.Series(np.nan, index=as_df.index) # Assign NaNs of correct length

        as_df[rank_col_name] = final_rank_scores
        rank_score_cols_for_mean.append(rank_col_name)

    if not rank_score_cols_for_mean: 
        as_df['Aggregated_Score'] = np.nan
    else:
        as_df['Aggregated_Score'] = as_df[rank_score_cols_for_mean].mean(axis=1, skipna=True) 
    return as_df

# --- III. Modified Main Orchestration for a Single Dataset ---
def run_benchmark_for_single_dataset(
    input_mode: str,
    dataset_id_tag: str,
    spatial_adata_path: str,
    reference_adata_path: Optional[str],  # MODIFIED HERE
    method_prediction_files: Dict[str, str], # MODIFIED HERE (if you want to be explicit)
    cell_type_obs_key_in_ref: Optional[str], # MODIFIED HERE
    n_top_deg_markers: int,
    metrics_for_as: List, # MODIFIED HERE (or List[str])
    higher_is_better_map: Dict # MODIFIED HERE (or Dict[str, bool])
) :
    print(f"\n===== Processing dataset: {dataset_id_tag} (Mode: {input_mode}) =====")

    marker_gene_map_generated = {}
    if (input_mode == "real_paired" or input_mode == "simulated") and \
       cell_type_obs_key_in_ref and reference_adata_path:
        if os.path.exists(reference_adata_path):
            try:
                print(f"  Loading reference AnnData for DEG markers for {dataset_id_tag} from: {reference_adata_path}")
                adata_ref_for_degs = sc.read_h5ad(reference_adata_path)
                if 'counts' in adata_ref_for_degs.layers:
                    adata_ref_for_degs.X = adata_ref_for_degs.layers['counts'].copy()
                # Basic check if .X might be raw counts (if 'counts' layer not present)
                elif scipy.sparse.issparse(adata_ref_for_degs.X) and adata_ref_for_degs.X.max() > 100 and np.issubdtype(adata_ref_for_degs.X.dtype, np.integer):
                     pass # Assume if sparse, large integer values, it's counts-like.
                sc.pp.normalize_total(adata_ref_for_degs, target_sum=1e4)
                sc.pp.log1p(adata_ref_for_degs)
                marker_gene_map_generated = generate_deg_marker_map(
                    adata_ref_for_degs, cell_type_obs_key_in_ref, n_top_genes=n_top_deg_markers
                )
            except Exception as e:
                print(f"    ERROR generating DEG marker map for {dataset_id_tag}: {e}")
        else:
            print(f"    INFO: Reference AnnData for DEG markers for {dataset_id_tag} not found at {reference_adata_path}. Marker map will be empty.")

    ground_truth_df = None
    adata_spatial_expression_norm = None

    if not os.path.exists(spatial_adata_path):
        print(f"  ERROR: Spatial AnnData file not found: {spatial_adata_path}. Cannot process dataset {dataset_id_tag}.")
        return None, None, None
    
    try:
        adata_spatial_raw = sc.read_h5ad(spatial_adata_path)
        adata_spatial_expression_norm = adata_spatial_raw.copy()
        if 'counts' in adata_spatial_expression_norm.layers:
            adata_spatial_expression_norm.X = adata_spatial_expression_norm.layers['counts'].copy()
        elif scipy.sparse.issparse(adata_spatial_expression_norm.X) and adata_spatial_expression_norm.X.max() > 100 and np.issubdtype(adata_spatial_expression_norm.X.dtype, np.integer):
            pass # Heuristic for raw counts
        sc.pp.normalize_total(adata_spatial_expression_norm, target_sum=1e4)
        sc.pp.log1p(adata_spatial_expression_norm)

        if input_mode == "simulated":
            if 'cell_type_proportions' in adata_spatial_raw.obsm and \
               'cell_type_names' in adata_spatial_raw.uns and \
               isinstance(adata_spatial_raw.uns['cell_type_names'], (list, np.ndarray)) and \
               len(adata_spatial_raw.uns['cell_type_names']) > 0 :
                gt_array = adata_spatial_raw.obsm['cell_type_proportions']
                gt_ct_names = list(adata_spatial_raw.uns['cell_type_names'])
                ground_truth_df = pd.DataFrame(gt_array, index=adata_spatial_raw.obs_names.astype(str), columns=gt_ct_names).fillna(0)
                # Ensure ground truth proportions sum to 1 (or handle minor float precision issues)
                ground_truth_df = normalize_proportions_df(ground_truth_df)
            else:
                print(f"  Warning: Ground truth data incomplete in {spatial_adata_path} for {dataset_id_tag}. GT metrics will be NaN.")
    except Exception as e:
        print(f"  ERROR loading spatial/GT AnnData {spatial_adata_path} for {dataset_id_tag}: {e}. Skipping dataset.")
        return None, None, None
    
    loaded_pred_dfs = {}
    for method, fpath in method_prediction_files.items():
        if os.path.exists(fpath):
            try:
                df = pd.read_csv(fpath, header=0, index_col=0)
                if not df.empty:
                    df.index = df.index.astype(str)
                    loaded_pred_dfs[method] = df
                else:
                    print(f"    INFO: File for {method} is empty: {fpath}")
            except Exception as e:
                print(f"    ERROR loading {method} from {fpath}: {e}")
        else:
            print(f"    INFO: File for {method} not found: {fpath}")
    
    if not loaded_pred_dfs:
        print(f"  No valid predictions loaded for {dataset_id_tag}. Skipping.")
        return None, None, None

    ref_spots = adata_spatial_raw.obs_names.astype(str)
    all_cts_set = set(ground_truth_df.columns) if ground_truth_df is not None else set()
    for df_ in loaded_pred_dfs.values():
        all_cts_set.update(df_.columns)
    ref_cts_list = sorted(list(all_cts_set))
    
    if not ref_cts_list:
        print(f"  ERROR: No cell types identified for alignment for {dataset_id_tag}. Skipping.")
        return None, None, None
    
    aligned_pred_dfs = {
        name: df.reindex(index=ref_spots, columns=ref_cts_list).fillna(0.0).replace([np.inf, -np.inf], 0.0)
        for name, df in loaded_pred_dfs.items()
    }

    current_dataset_main_summaries_list = []
    current_dataset_detailed_gt_dfs_list = []
    current_dataset_detailed_marker_dfs_list = []

    for method_name, pred_df_aligned in aligned_pred_dfs.items():
        print(f"    Evaluating {method_name} for {dataset_id_tag}...")
        if pred_df_aligned.shape[0] == 0 or pred_df_aligned.shape[1] == 0:
            print(f"      Skipping {method_name} due to empty aligned DataFrame.")
            continue

        pred_df_norm = normalize_proportions_df(pred_df_aligned)
        
        overall_summary_stats_method = {}
        
        if ground_truth_df is not None:
            gt_metrics_results = calculate_gt_accuracy_metrics(pred_df_norm, ground_truth_df)
            overall_summary_stats_method.update(gt_metrics_results['overall'])
            if not gt_metrics_results['by_cell_type'].empty:
                df_detail_gt = gt_metrics_results['by_cell_type'].copy()
                df_detail_gt['Method'] = method_name
                df_detail_gt['DatasetID'] = dataset_id_tag
                current_dataset_detailed_gt_dfs_list.append(df_detail_gt.reset_index().rename(columns={'index':'CellType'}))

        intrinsic_metrics_results = calculate_intrinsic_metrics(pred_df_norm)
        overall_summary_stats_method.update(intrinsic_metrics_results['overall'])

        if marker_gene_map_generated and adata_spatial_expression_norm is not None:
            marker_metrics_results = calculate_marker_correlation_metric(
                pred_df_norm, adata_spatial_expression_norm, marker_gene_map_generated
            )
            overall_summary_stats_method.update(marker_metrics_results['overall'])
            if not marker_metrics_results['by_cell_type'].empty:
                df_detail_marker = marker_metrics_results['by_cell_type'].copy()
                df_detail_marker['Method'] = method_name
                df_detail_marker['DatasetID'] = dataset_id_tag
                current_dataset_detailed_marker_dfs_list.append(df_detail_marker.reset_index().rename(columns={'index':'CellType'}))
        
        main_summary_row_df = pd.DataFrame([overall_summary_stats_method], index=[method_name])
        # 'DatasetID' will be added after concatenating all methods for this dataset
        current_dataset_main_summaries_list.append(main_summary_row_df)

    # --- Consolidate results for the current dataset ---
    final_summary_df_this_dataset = None
    if current_dataset_main_summaries_list:
        # Concatenate summaries from all methods for this dataset
        summary_df_this_dataset_all_methods = pd.concat(current_dataset_main_summaries_list)
        summary_df_this_dataset_all_methods = summary_df_this_dataset_all_methods.rename_axis('Method').reset_index()
        summary_df_this_dataset_all_methods['DatasetID'] = dataset_id_tag
        
        # Calculate aggregated score if applicable (simulated mode with GT)
        if input_mode == "simulated" and ground_truth_df is not None:
            # Set index for calculate_aggregated_score if it expects 'Method' in index
            temp_summary_for_as = summary_df_this_dataset_all_methods.set_index(['Method', 'DatasetID'])
            temp_summary_for_as = calculate_aggregated_score(temp_summary_for_as, metrics_for_as, higher_is_better_map)
            final_summary_df_this_dataset = temp_summary_for_as.reset_index()
        else:
            final_summary_df_this_dataset = summary_df_this_dataset_all_methods

        print(f"\n--- Summary for Dataset {dataset_id_tag} (Mode: {input_mode}) ---")
        # Ensure DatasetID is one of the first columns for readability
        cols = ['DatasetID', 'Method'] + [c for c in final_summary_df_this_dataset.columns if c not in ['DatasetID', 'Method']]
        print(final_summary_df_this_dataset[cols].to_string())
    else:
        print(f"  No method summaries generated for {dataset_id_tag}.")

    final_detailed_gt_df_this_dataset = pd.concat(current_dataset_detailed_gt_dfs_list, ignore_index=True) if current_dataset_detailed_gt_dfs_list else None
    final_detailed_marker_df_this_dataset = pd.concat(current_dataset_detailed_marker_dfs_list, ignore_index=True) if current_dataset_detailed_marker_dfs_list else None
            
    return final_summary_df_this_dataset, final_detailed_gt_df_this_dataset, final_detailed_marker_df_this_dataset


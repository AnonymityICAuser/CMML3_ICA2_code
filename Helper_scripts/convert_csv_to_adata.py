import scipy.sparse as sparse
import numpy as np
import scanpy as sc
import os
import pandas as pd

# --- Helper for unique index ---
def make_index_unique_if_needed(index, name_for_log=""):
    if not isinstance(index, pd.Index): # Ensure it's a pandas Index
        index = pd.Index(index)
    if index.duplicated().any():
        print(f"Warning: Duplicate values found in {name_for_log} index. Making them unique using sc.utils.make_index_unique().")
        return sc.utils.make_index_unique(pd.Index(index.astype(str)))
    return index

# --- Robust CSV/TSV Reader ---
def read_csv_or_tsv(file_path, index_col=0, header=0, is_expr_matrix_no_header=False, **kwargs):
    """
    Attempts to read a file as CSV, then as TSV.
    is_expr_matrix_no_header: Special flag for expression matrices where genes are index,
                              cells are columns, and there's NO header for cell columns.
    """
    filename = os.path.basename(file_path)
    # print(f"Attempting to read: {filename}, index_col={index_col}, header={header}, is_expr_matrix_no_header={is_expr_matrix_no_header}")
    try:
        if is_expr_matrix_no_header and (header is None or header==0):
            try:
                # For "Gene val val" with no cell header, header=None is critical.
                # index_col=0 takes the first column as index (genes).
                df = pd.read_csv(file_path, index_col=index_col, header=None, sep=r'\s+', engine='python', **kwargs)
                print(f"Successfully read {filename} as space/tab-separated (expr_matrix_no_header mode). Shape: {df.shape}")
                # If first row was mistaken for header, it might be in df.columns. Fix if needed.
                # This case usually applies if index_col was NOT 0.
                # If index_col=0, the first column is index, rest are data.
                # If after this, df.shape[1] is 0 but df.shape[0] > 0, it means only index was read.
                if df.shape[1] == 0 and df.shape[0] > 0 and index_col is not None:
                    print(f"Warning: {filename} (expr_matrix_no_header mode) read with 0 data columns. Check file content/format.")
                return df
            except Exception as e_tsv_first:
                print(f"Attempt to read {filename} as space/tab (expr_matrix_no_header mode) failed: {e_tsv_first}. Trying CSV mode for it.")
                # Fallback to trying as CSV even in this special mode, though less likely to be correct
                df = pd.read_csv(file_path, index_col=index_col, header=None, sep=',', engine='python', **kwargs)
                print(f"Fallback CSV read for {filename} (expr_matrix_no_header mode). Shape: {df.shape}")
                return df

        # Standard CSV attempt
        df = pd.read_csv(file_path, index_col=index_col, header=header, sep=',', engine='python', **kwargs)
        # Heuristic for misparsed TSV (single data column after index)
        is_single_data_column = (df.shape[1] == 1 and index_col is None) or \
                                (df.shape[1] == 0 and index_col is not None and df.index.name is not None and df.index.name in df.columns) # index took the only data col name
        
        if is_single_data_column and os.path.getsize(file_path) > 30 :
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                first_line = f.readline().strip()
                # Check if the first line (potentially header) or second line (data) contains tabs/multiple spaces
                second_line = f.readline().strip() if os.path.getsize(file_path) > len(first_line) + 5 else ""
                if ('\t' in first_line or '  ' in first_line) or \
                   ('\t' in second_line or '  ' in second_line) :
                    print(f"CSV read of {filename} resulted in a single data column. Attempting TSV/space-separated read.")
                    raise ValueError("Likely TSV/space-separated")
        print(f"Successfully read {filename} as CSV. Shape: {df.shape}")
        return df
        
    except (pd.errors.ParserError, ValueError) as e_csv:
        # print(f"CSV parsing for {filename} failed or triggered heuristic: {e_csv}. Trying TSV/space-separated.")
        try:
            effective_header = None if is_expr_matrix_no_header else header
            # For files like "Jam2 1 0 1...", index_col=0 and header=None (if no cell headers) is crucial for sep=r'\s+'
            df = pd.read_csv(file_path, index_col=index_col, header=effective_header, sep=r'\s+', engine='python', **kwargs)
            print(f"Successfully read {filename} as space/tab-separated. Shape: {df.shape}")
            if df.shape[1] == 0 and df.shape[0] > 0 and index_col is not None:
                 print(f"Warning: {filename} (TSV mode) read with 0 data columns after index. Check file content/format.")
            return df
        except Exception as e_tsv:
            print(f"Error: Could not parse {filename} as CSV or TSV/space-separated: {e_tsv}")
            raise

# --- Function 1: Load Spatial Data ---
def csv_to_spatial_adata(input_dir,
                         count_file="st_count.csv", 
                         location_file="st_location.csv", 
                         count_index_col=0, count_header=0,
                         location_index_col=0, location_header=0,
                         spot_id_col_in_loc_df=None, 
                         coord_col_names=None):
    print(f"\n--- Loading Spatial Data from: {input_dir} ---")
    
    # 1. Load Spatial Locations
    loc_path = os.path.join(input_dir, location_file)
    if not os.path.exists(loc_path):
        print(f"Error: Spatial location file not found: {loc_path}"); return None
    try:
        loc_df = read_csv_or_tsv(loc_path, index_col=location_index_col, header=location_header)
        if spot_id_col_in_loc_df:
            if spot_id_col_in_loc_df not in loc_df.columns:
                print(f"Error: spot_id_col_in_loc_df '{spot_id_col_in_loc_df}' not in {location_file}."); return None
            loc_df = loc_df.set_index(spot_id_col_in_loc_df)
        
        loc_df.index = make_index_unique_if_needed(loc_df.index.astype(str), f"{location_file} spot IDs")
        print(f"Loaded {location_file}: {loc_df.shape[0]} spots, {loc_df.shape[1]} attributes.")
    except Exception as e:
        print(f"Error loading {location_file}: {e}"); return None

    final_obs_names = loc_df.index # These are the unique, string-based spot IDs
    print(f"Debug (Spatial): Using final_obs_names (length {len(final_obs_names)}), first 5: {final_obs_names.tolist()[:5]}")

    # Identify coordinate columns
    actual_coord_cols = []
    if coord_col_names and len(coord_col_names) == 2 and all(c in loc_df.columns for c in coord_col_names):
        actual_coord_cols = list(coord_col_names)
    else:
        coord_options = [('array_row', 'array_col'), ('x', 'y'), ('imagecol', 'imagerow'),
                         ('pxl_col_in_fullres', 'pxl_row_in_fullres'), ('spatial1', 'spatial2'),
                         ('coord_x', 'coord_y'), ('X', 'Y')]
        for c1, c2 in coord_options:
            if c1 in loc_df.columns and c2 in loc_df.columns:
                actual_coord_cols = [c1, c2]; break
        if not actual_coord_cols:
            numeric_cols = loc_df.select_dtypes(include=np.number).columns
            potential_cols = [col for col in numeric_cols if col != loc_df.index.name] # Exclude index if it was named and numeric
            if len(potential_cols) >=2:
                actual_coord_cols = potential_cols[:2]
                print(f"Auto-detected numeric coordinate columns: {actual_coord_cols}")
            elif len(loc_df.columns) >=2 : # Last resort, take first two columns if all else fails
                actual_coord_cols = loc_df.columns[:2].tolist()
                print(f"Warning: Using first two columns as coordinates: {actual_coord_cols}")


    if not actual_coord_cols or len(actual_coord_cols) < 2:
        print(f"Error: Could not reliably determine 2 coordinate columns in {location_file}. Found: {actual_coord_cols}"); return None
    print(f"Using coordinate columns: {actual_coord_cols}")
    
    try:
        # loc_df's index is already final_obs_names
        spatial_coords_df = loc_df[actual_coord_cols].copy().astype(float)
    except Exception as e:
        print(f"Error extracting/converting coordinate columns {actual_coord_cols} from loc_df: {e}"); return None

    # 2. Load Spatial Counts
    count_path = os.path.join(input_dir, count_file)
    if not os.path.exists(count_path):
        print(f"Error: Spatial count file not found: {count_path}"); return None
    try:
        # count_header=None if genes are index and spots have no header.
        counts_df_raw = read_csv_or_tsv(count_path, index_col=count_index_col, header=count_header, 
                                        is_expr_matrix_no_header=(count_header is None))
        print(f"Loaded {count_file}: shape {counts_df_raw.shape}")
    except Exception as e:
        print(f"Error loading {count_file}: {e}"); return None

    # Determine orientation
    idx_overlap = np.isin(counts_df_raw.index.astype(str), final_obs_names).mean()
    col_overlap = np.isin(counts_df_raw.columns.astype(str), final_obs_names).mean()
    X_st_oriented = None; st_var_names_pd_idx = None; overlap_threshold = 0.75

    if counts_df_raw.empty:
        print(f"Error: {count_file} loaded as an empty DataFrame."); return None

    if idx_overlap >= overlap_threshold and (counts_df_raw.shape[1] > 0 and col_overlap < overlap_threshold / 2):
        print(f"Interpreting {count_file} as Spots x Genes."); X_st_oriented = counts_df_raw; st_var_names_pd_idx = X_st_oriented.columns
    elif col_overlap >= overlap_threshold and (counts_df_raw.shape[0] > 0 and idx_overlap < overlap_threshold / 2):
        print(f"Interpreting {count_file} as Genes x Spots. Transposing."); X_st_oriented = counts_df_raw.T; st_var_names_pd_idx = X_st_oriented.columns
    elif idx_overlap >= overlap_threshold and col_overlap >= overlap_threshold : # Both high, check shape
        print(f"Warning: Both rows/cols of {count_file} highly overlap spot IDs. Assuming Spots x Genes if rows match spot count better.")
        if abs(counts_df_raw.shape[0] - len(final_obs_names)) <= abs(counts_df_raw.shape[1] - len(final_obs_names)):
             X_st_oriented = counts_df_raw; st_var_names_pd_idx = X_st_oriented.columns
        else: X_st_oriented = counts_df_raw.T; st_var_names_pd_idx = X_st_oriented.columns; print("  Transposed.")
    else: # Low overlap, try to guess
        if counts_df_raw.shape[0] > 0 and counts_df_raw.shape[1] > 0:
            if col_overlap > idx_overlap + 0.1 : X_st_oriented = counts_df_raw.T; st_var_names_pd_idx = X_st_oriented.columns; print(f"Low overlap ({idx_overlap:.2f} vs {col_overlap:.2f}), assuming Genes x Spots due to col match. Transposed.")
            elif idx_overlap > col_overlap + 0.1 : X_st_oriented = counts_df_raw; st_var_names_pd_idx = X_st_oriented.columns; print(f"Low overlap ({idx_overlap:.2f} vs {col_overlap:.2f}), assuming Spots x Genes due to row match.")
            else: print(f"Error: Cannot confidently determine orientation of {count_file}. idx_overlap={idx_overlap:.2f}, col_overlap={col_overlap:.2f}"); return None
        else: print(f"Error: Count matrix {count_file} is empty or has zero dimension after orientation check."); return None

    X_st_oriented.index = make_index_unique_if_needed(X_st_oriented.index.astype(str), f"{count_file} oriented spot IDs")
    if not isinstance(st_var_names_pd_idx, pd.Index): st_var_names_pd_idx = pd.Index(st_var_names_pd_idx)
    st_var_names_final = make_index_unique_if_needed(st_var_names_pd_idx.astype(str), f"{count_file} gene names")
    X_st_oriented.columns = st_var_names_final

    # 3. Align and Create AnnData
    X_aligned = X_st_oriented.reindex(index=final_obs_names, columns=st_var_names_final)
    if X_aligned.isnull().values.any():
        print(f"Warning: NaNs found in aligned expression matrix X_aligned. Filling with 0. Spot IDs in counts might not perfectly match location IDs.")
        X_aligned = X_aligned.fillna(0)
    
    adata_obs = pd.DataFrame(index=final_obs_names)
    other_loc_cols = loc_df.columns.difference(pd.Index(actual_coord_cols)) 
    if not other_loc_cols.empty:
        adata_obs = adata_obs.join(loc_df[other_loc_cols]) # loc_df index is already final_obs_names
        
    adata_var = pd.DataFrame(index=st_var_names_final)
    
    # spatial_coords_df's index is final_obs_names.
    # We need to ensure the order matches final_obs_names for adata.obsm['spatial']
    # This reindex also ensures that if spatial_coords_df somehow had extra spots, they are dropped.
    spatial_coords_for_adata = spatial_coords_df.reindex(final_obs_names)

    if spatial_coords_for_adata.isnull().values.any():
        print(f"CRITICAL ERROR: NaNs found in spatial_coords_for_adata after reindexing to final_obs_names. This indicates a mismatch.")
        print(f"  spatial_coords_df index (first 5): {spatial_coords_df.index.tolist()[:5]}")
        print(f"  final_obs_names (first 5): {final_obs_names.tolist()[:5]}")
        # Check for specific mismatches
        missing_in_spatial_df = final_obs_names[~final_obs_names.isin(spatial_coords_df.index)]
        if len(missing_in_spatial_df) > 0:
            print(f"  Spot IDs in final_obs_names but NOT in spatial_coords_df.index (first 5): {missing_in_spatial_df.tolist()[:5]}")
        return None

    try:
        adata_st = sc.AnnData(X=sparse.csr_matrix(X_aligned.values.astype(np.float32)), obs=adata_obs, var=adata_var)
        adata_st.obsm['spatial'] = spatial_coords_for_adata[actual_coord_cols].values.astype(np.float32)
        print(f"Successfully created spatial AnnData: {adata_st.n_obs} spots, {adata_st.n_vars} genes.")
        if np.isnan(adata_st.obsm['spatial']).any():
            print("WARNING: Final adata_st.obsm['spatial'] still contains NaNs!")
        return adata_st
    except Exception as e:
        print(f"Error creating spatial AnnData object: {e}"); return None

def csv_to_single_cell_adata(input_dir,
                             expression_file="sc_mousebrain.csv", 
                             celltype_file=None, 
                             expr_index_col=0, expr_header=0, 
                             expr_cells_are_rows=None, 
                             celltype_index_col=0, celltype_header=0, 
                             celltype_id_col_in_df=None, 
                             celltype_annot_col_name=None 
                             ):
    print(f"\n--- Loading Single-Cell Data from: {input_dir} ---")

    # 1. Load Expression Data
    expr_path = os.path.join(input_dir, expression_file)
    if not os.path.exists(expr_path):
        print(f"Error: SC expression file not found: {expr_path}"); return None
    try:
        is_expr_g_x_c_no_header = (expr_header is None and (expr_cells_are_rows is False or expr_cells_are_rows is None))
        # Read all expression data as string first to handle potential mixed types before numeric conversion
        expr_df_raw = read_csv_or_tsv(expr_path, index_col=expr_index_col, header=expr_header, 
                                      is_expr_matrix_no_header=is_expr_g_x_c_no_header, dtype=str) 
        print(f"Loaded {expression_file}: shape {expr_df_raw.shape} (read as string)")
    except Exception as e:
        print(f"Error loading {expression_file}: {e}"); return None

    # 2. Load Cell Type Data (as before)
    sc_cell_ids_from_types_pd_idx = None
    cell_annotations_df = None
    if celltype_file:
        ctype_path = os.path.join(input_dir, celltype_file)
        if not os.path.exists(ctype_path): print(f"Warning: SC celltype file '{ctype_path}' not found.")
        else:
            try:
                ctype_df_loaded = None
                if celltype_index_col is None and celltype_header is None:
                    try: ctype_df_loaded = pd.read_csv(ctype_path, header=None, sep=r'\s+', engine='python')
                    except pd.errors.EmptyDataError: ctype_df_loaded = pd.read_csv(ctype_path, header=None, sep=',', engine='python')
                    except Exception: 
                        try: ctype_df_loaded = pd.read_csv(ctype_path, header=None, sep=',', engine='python')
                        except Exception as e_c: print(f"Could not read {celltype_file} (no-header/no-index): {e_c}"); raise
                    if ctype_df_loaded.shape[1] > 1: cell_annotations_df = ctype_df_loaded.iloc[:, [0]].copy(); print(f"Warning: {celltype_file} read with >1 col, using first.")
                    elif ctype_df_loaded.shape[1] == 1: cell_annotations_df = ctype_df_loaded.copy()
                    else: print(f"Error: {celltype_file} loaded with 0 columns."); return None
                    annot_col_name_to_use = celltype_annot_col_name if celltype_annot_col_name else 'cell_type_auto'
                    cell_annotations_df.columns = [annot_col_name_to_use]
                    sc_cell_ids_from_types_pd_idx = cell_annotations_df.index # RangeIndex
                else: 
                    ctype_df_loaded = read_csv_or_tsv(ctype_path, index_col=celltype_index_col, header=celltype_header)
                    if celltype_id_col_in_df: 
                        if celltype_id_col_in_df not in ctype_df_loaded.columns: print(f"Error: {celltype_id_col_in_df} not in {celltype_file}."); return None
                        ctype_df_loaded = ctype_df_loaded.set_index(celltype_id_col_in_df)
                    ctype_df_loaded.index = make_index_unique_if_needed(ctype_df_loaded.index, f"{celltype_file} cell IDs")
                    sc_cell_ids_from_types_pd_idx = ctype_df_loaded.index
                    if celltype_annot_col_name and celltype_annot_col_name in ctype_df_loaded.columns:
                        cell_annotations_df = ctype_df_loaded[[celltype_annot_col_name]].copy()
                    elif ctype_df_loaded.shape[1] == 1: cell_annotations_df = ctype_df_loaded.copy()
                    else: 
                        non_numeric_cols = ctype_df_loaded.select_dtypes(exclude=np.number).columns
                        if len(non_numeric_cols) > 0: cell_annotations_df = ctype_df_loaded[[non_numeric_cols[0]]].copy(); print(f"Auto-selected '{non_numeric_cols[0]}' from {celltype_file} as annotation.")
                        elif ctype_df_loaded.shape[1] > 0: cell_annotations_df = ctype_df_loaded[[ctype_df_loaded.columns[0]]].copy(); print(f"Warning: Using first col '{ctype_df_loaded.columns[0]}' from {celltype_file}.")
                if cell_annotations_df is not None: print(f"Loaded {celltype_file}: {cell_annotations_df.shape[0]} entries, annotation col '{cell_annotations_df.columns[0]}'. Index type: {type(sc_cell_ids_from_types_pd_idx)}")
                else: print(f"Loaded {celltype_file} but no annotations extracted.")
            except Exception as e: print(f"Error processing {celltype_file}: {e}"); sc_cell_ids_from_types_pd_idx = None; cell_annotations_df = None
    
    # 3. Orient and Align Expression Matrix
    if expr_df_raw.empty or (expr_df_raw.shape[1] == 0 and expr_index_col is not None) :
        print(f"Error: Expression DataFrame {expression_file} is empty or has 0 data columns after loading. Cannot proceed.")
        return None

    expr_df_oriented = None; sc_var_names_pd_idx = None; sc_obs_names_pd_idx = None
    if expr_cells_are_rows is None: # Auto-detect
        if sc_cell_ids_from_types_pd_idx is not None and not isinstance(sc_cell_ids_from_types_pd_idx, pd.RangeIndex): 
            idx_overlap_mean = np.isin(expr_df_raw.index.map(str), sc_cell_ids_from_types_pd_idx.map(str)).mean()
            col_is_str_like = expr_df_raw.columns.inferred_type in ['string', 'unicode', 'mixed'] or all(isinstance(c, str) for c in expr_df_raw.columns)
            if col_is_str_like:
                col_overlap_mean = np.isin(expr_df_raw.columns.map(str), sc_cell_ids_from_types_pd_idx.map(str)).mean()
            else: col_overlap_mean = 0.0
            if idx_overlap_mean > col_overlap_mean + 0.1: expr_cells_are_rows = True
            elif col_overlap_mean > idx_overlap_mean + 0.1: expr_cells_are_rows = False
            else: expr_cells_are_rows = expr_df_raw.shape[0] >= expr_df_raw.shape[1]; print(f"Ambiguous orientation for {expression_file} via celltypes. Guessed cells_are_rows={expr_cells_are_rows} by shape.")
        else: 
            col_is_numeric = all(isinstance(c, int) or isinstance(c, float) for c in expr_df_raw.columns)
            idx_is_str_like = expr_df_raw.index.inferred_type in ['string', 'unicode', 'mixed'] or all(isinstance(i, str) for i in expr_df_raw.index)
            if col_is_numeric and idx_is_str_like: expr_cells_are_rows = False; print(f"{expression_file} has string-like index & num cols. Assuming Genes x Cells.")
            else: expr_cells_are_rows = expr_df_raw.shape[0] >= expr_df_raw.shape[1]; print(f"No usable celltype IDs for {expression_file} orientation. Assuming cells_are_rows={expr_cells_are_rows} by shape.")
    
    if expr_cells_are_rows:
        print(f"Interpreting {expression_file} as Cells x Genes."); expr_df_oriented = expr_df_raw
        sc_obs_names_pd_idx = expr_df_oriented.index; sc_var_names_pd_idx = expr_df_oriented.columns
    else: 
        print(f"Interpreting {expression_file} as Genes x Cells. Transposing."); expr_df_oriented = expr_df_raw.T
        sc_obs_names_pd_idx = expr_df_oriented.index ; sc_var_names_pd_idx = expr_df_oriented.columns
    
    # Define final_sc_obs_names
    if sc_cell_ids_from_types_pd_idx is not None and not isinstance(sc_cell_ids_from_types_pd_idx, pd.RangeIndex):
        final_sc_obs_names_pd_idx = sc_cell_ids_from_types_pd_idx
    else: # Includes case where celltype file gave RangeIndex or no celltype file
        final_sc_obs_names_pd_idx = sc_obs_names_pd_idx 
        # If sc_obs_names_pd_idx is also RangeIndex (e.g. from transposed GxP matrix with no P header)
        # AND celltype file also gave RangeIndex, make sure they have same length for positional alignment
        if isinstance(sc_obs_names_pd_idx, pd.RangeIndex) and \
           isinstance(sc_cell_ids_from_types_pd_idx, pd.RangeIndex) and \
           cell_annotations_df is not None and \
           len(sc_obs_names_pd_idx) != len(sc_cell_ids_from_types_pd_idx):
            print(f"Warning: Length mismatch between expression cells ({len(sc_obs_names_pd_idx)}) and celltype labels ({len(sc_cell_ids_from_types_pd_idx)}) when both are order-based. Annotations might be misaligned.")
            # Default to expression cell count for AnnData shape
    
    final_sc_obs_names = make_index_unique_if_needed(final_sc_obs_names_pd_idx, "final SC cell IDs")
    final_sc_var_names = make_index_unique_if_needed(pd.Index(sc_var_names_pd_idx), "final SC gene names")

    expr_df_oriented.index = expr_df_oriented.index.map(str)
    expr_df_oriented.columns = expr_df_oriented.columns.map(str) # Ensure columns are also string before reindex
    
    X_aligned_str = expr_df_oriented.reindex(index=final_sc_obs_names, columns=final_sc_var_names)
    
    # Convert to numeric, coercing errors to NaN, then fill NaNs with 0
    # This is where the '1772071015_C02' would become NaN then 0
    X_numeric = X_aligned_str.apply(pd.to_numeric, errors='coerce').fillna(0)
    num_genes_final = X_numeric.shape[1]

    if num_genes_final != len(final_sc_var_names):
         print(f"Warning: Number of columns in numeric matrix ({num_genes_final}) does not match final_sc_var_names ({len(final_sc_var_names)}). This could be due to all-NaN columns being dropped by apply(pd.to_numeric). Re-creating with all expected genes.")
         X_numeric = X_numeric.reindex(columns=final_sc_var_names, fill_value=0)


    adata_obs_sc = pd.DataFrame(index=final_sc_obs_names)
    if cell_annotations_df is not None:
        if isinstance(sc_cell_ids_from_types_pd_idx, pd.RangeIndex) and \
           isinstance(final_sc_obs_names_pd_idx, pd.RangeIndex) and \
           len(sc_cell_ids_from_types_pd_idx) == len(final_sc_obs_names_pd_idx): # Both order-based and same length
            cell_annotations_df.index = final_sc_obs_names 
        else: # Align by actual ID values if available
             cell_annotations_df = cell_annotations_df.reindex(final_sc_obs_names)
        for col_annot in cell_annotations_df.columns:
             adata_obs_sc[col_annot] = cell_annotations_df[col_annot]
    
    adata_var_sc = pd.DataFrame(index=final_sc_var_names)

    try:
        adata_sc = sc.AnnData(X=sparse.csr_matrix(X_numeric.values.astype(np.float32)), obs=adata_obs_sc, var=adata_var_sc)
        print(f"Successfully created single-cell AnnData: {adata_sc.n_obs} cells, {adata_sc.n_vars} genes.")
        if adata_sc.n_vars == 0 and len(final_sc_var_names) > 0 :
            print(f"CRITICAL WARNING: Single-cell AnnData created with 0 genes, but expected {len(final_sc_var_names)}. Check expression file parsing, orientation, and numeric conversion.")
        return adata_sc
    except Exception as e:
        print(f"Error creating single-cell AnnData object: {e}"); return None

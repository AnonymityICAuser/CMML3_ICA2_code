# Load required libraries
library(Giotto)
library(Matrix)
library(data.table)
library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("-s", "--st_dir"), type="character", default=NULL, help="Directory containing spatial transcriptomics data"),
  make_option(c("-c", "--sc_dir"), type="character", default=NULL, help="Directory containing single-cell data"),
  make_option(c("-t", "--cell_type_path"), type="character", default=NULL, help="Path to cell type ground truth CSV file for single-cell data"),
  make_option(c("-o", "--output_path"), type="character", default="cell_type_dwls_output.csv", help="Output file path [default=%default]"),
  make_option(c("-r", "--resolution"), type="numeric", default=0.4,  help="Resolution for Leiden clustering on ST data [default=%default]"),
  make_option(c("-f", "--scale_factor"), type="numeric", default=6000, help="Scale factor for normalization [default=%default]"),
  make_option(c("-m", "--markers"), type="numeric", default=10, help="Number of top markers per cell type for signature matrix [default=%default]"),
  # Removed hvf_method and skip_hvf options
  make_option(c("--pca_dims"), type="integer", default=10, help="Number of principal components for clustering [default=%default]"),
  make_option(c("--nn_k"), type="integer", default=15, help="Number of nearest neighbors for graph construction [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list, description="Perform spatial deconvolution using DWLS method with actual ST Leiden clustering. PCA uses all common genes.")
opt <- parse_args(opt_parser)

if (is.null(opt$st_dir) || is.null(opt$sc_dir) || is.null(opt$cell_type_path)) {
  print_help(opt_parser)
  stop("Missing required arguments: st_dir, sc_dir, and cell_type_path must be provided", call.=FALSE)
}

run_spatial_dwls <- function(st_dir, sc_dir, cell_type_path, output_path,
                             resolution, scale_factor, top_n_markers,
                             pca_dims_to_use, nn_k_val) { # Removed hvf_method and skip_hvf_calc
  print("Starting run_spatial_dwls function...")
  instrs <- createGiottoInstructions(save_plot = FALSE, show_plot = FALSE, return_plot = FALSE, python_path = NULL)

  # --- Read and Prepare Single-Cell Data ---
  print("Reading and preparing single-cell data...")
  sc_cell_types <- fread(cell_type_path, header=TRUE)
  if("cell_label" %in% colnames(sc_cell_types) && "class" %in% colnames(sc_cell_types)) {
    setnames(sc_cell_types, c("cell_label", "class"), c("cell_ID", "cellType"))
  } else if (!("cell_ID" %in% colnames(sc_cell_types) && "cellType" %in% colnames(sc_cell_types))) {
    warning("Cell type CSV does not have expected 'cell_ID'/'cellType' or 'cell_label'/'class' columns. Using first two columns.")
    if(ncol(sc_cell_types) >= 2) setnames(sc_cell_types, colnames(sc_cell_types)[1:2], c("cell_ID", "cellType"))
    else stop("Cell type CSV must have at least two columns.")
  }
  sc_cell_info <- fread(file.path(sc_dir, "cellinfo.csv"), skip = 1, header = FALSE); sc_gene_info <- fread(file.path(sc_dir, "geneinfo.csv"), skip = 1, header = FALSE)
  sc_cell_ids <- sc_cell_info[[1]]; sc_gene_ids <- sc_gene_info[[1]]

  cell_type_dt <- sc_cell_types[, .(cell_ID = as.character(cell_ID), cellType = as.character(cellType))]
  sc_meta_dt <- data.table(cell_ID = as.character(sc_cell_ids))
  sc_meta_dt <- merge(sc_meta_dt, cell_type_dt, by = "cell_ID", all.x = TRUE)
  sc_meta_dt <- sc_meta_dt[!is.na(cellType)]
  sc_meta_dt[, sampleInfo := "sample1"]
  setkey(sc_meta_dt, cell_ID)

  sc_meta <- as.data.frame(sc_meta_dt)
  rownames(sc_meta) <- sc_meta$cell_ID

  if(nrow(sc_meta) == 0) stop("No cells remaining in single-cell metadata after merging and filtering for cell types.")
  sc_cell_ids_filtered <- rownames(sc_meta)

  sc_sparse_matrix <- readMM(file.path(sc_dir, "sparse_matrix.mtx"))
  dim_sc_genes <- min(nrow(sc_sparse_matrix), length(sc_gene_ids)); dim_sc_cells <- min(ncol(sc_sparse_matrix), length(sc_cell_ids))
  sc_count_pre <- sc_sparse_matrix[1:dim_sc_genes, 1:dim_sc_cells, drop=FALSE]
  rownames(sc_count_pre) <- sc_gene_ids[1:dim_sc_genes]; colnames(sc_count_pre) <- sc_cell_ids[1:dim_sc_cells]

  cells_to_keep_sc <- intersect(colnames(sc_count_pre), sc_cell_ids_filtered)
  if(length(cells_to_keep_sc) == 0) stop("No common cells between SC count matrix and filtered metadata.")
  sc_count_all_genes_orig <- sc_count_pre[, cells_to_keep_sc, drop=FALSE]

  sc_meta <- sc_meta[cells_to_keep_sc, , drop=FALSE]

  if(ncol(sc_count_all_genes_orig) == 0 || nrow(sc_count_all_genes_orig) == 0) stop("SC count matrix is empty after filtering.")

  print("Reading and preparing spatial data (raw)...")
  st_cell_info <- fread(file.path(st_dir, "cellinfo.csv"), skip = 1, header = FALSE); st_gene_info <- fread(file.path(st_dir, "geneinfo.csv"), skip = 1, header = FALSE)
  st_cell_ids <- st_cell_info[[1]]; st_gene_ids <- st_gene_info[[1]]
  st_sparse_matrix <- readMM(file.path(st_dir, "sparse_matrix.mtx"))
  dim_st_genes <- min(nrow(st_sparse_matrix), length(st_gene_ids)); dim_st_cells <- min(ncol(st_sparse_matrix), length(st_cell_ids))
  spatial_count_all_genes_orig_pre <- st_sparse_matrix[1:dim_st_genes, 1:dim_st_cells, drop=FALSE]
  rownames(spatial_count_all_genes_orig_pre) <- st_gene_ids[1:dim_st_genes]; colnames(spatial_count_all_genes_orig_pre) <- st_cell_ids[1:dim_st_cells]
  spatial_count_all_genes_orig <- spatial_count_all_genes_orig_pre

  if (ncol(st_cell_info) >= 3 && is.numeric(st_cell_info[[2]]) && is.numeric(st_cell_info[[3]])) {
    spatial_location_df <- data.frame(x = st_cell_info[[2]], y = st_cell_info[[3]])
    rownames(spatial_location_df) <- st_cell_info[[1]]
    cells_to_keep_st_locs <- intersect(rownames(spatial_location_df), colnames(spatial_count_all_genes_orig))
    if(length(cells_to_keep_st_locs) == 0) stop("No common cells between ST count matrix and spatial locations.")
    spatial_location_for_giotto <- spatial_location_df[cells_to_keep_st_locs, , drop = FALSE]
    spatial_count_all_genes_orig <- spatial_count_all_genes_orig[, cells_to_keep_st_locs, drop = FALSE]
    print("Using x and y coordinates from st/cellinfo.csv")
  } else {
    print("Generating placeholder grid coordinates for spatial data.")
    n_spots_loc <- ncol(spatial_count_all_genes_orig)
    n_cols_grid_loc <- ceiling(sqrt(n_spots_loc))
    spatial_location_for_giotto <- data.frame(
      x = (seq_len(n_spots_loc) - 1) %% n_cols_grid_loc,
      y = (seq_len(n_spots_loc) - 1) %/% n_cols_grid_loc,
      row.names = colnames(spatial_count_all_genes_orig),
      stringsAsFactors = FALSE )
  }
  if(nrow(spatial_location_for_giotto) != ncol(spatial_count_all_genes_orig)) stop("Mismatch between spatial locations and spatial_count spots after reconciling.")

  common_genes <- intersect(rownames(sc_count_all_genes_orig), rownames(spatial_count_all_genes_orig))
  if(length(common_genes) == 0) {
    print("No common genes found. Skipping deconvolution.")
    empty_results <- data.table(cell_ID = character(0), Note = "No common genes found between SC and ST")
    write.csv(empty_results, output_path, row.names = FALSE)
    return(NULL)
  }
  print(paste("Number of common genes:", length(common_genes)))

  sc_count_common <- sc_count_all_genes_orig[common_genes, , drop=FALSE]
  spatial_count_common <- spatial_count_all_genes_orig[common_genes, , drop=FALSE]
  spatial_location_for_giotto <- spatial_location_for_giotto[colnames(spatial_count_common), , drop = FALSE]
  current_sc_meta <- sc_meta[colnames(sc_count_common), , drop=FALSE]

  print(paste("Dimensions of sc_count_common (for signature):", paste(dim(sc_count_common), collapse="x")))
  print(paste("Dimensions of current_sc_meta (for signature):", paste(dim(current_sc_meta), collapse="x")))
  print(paste("Dimensions of spatial_count_common (for ST object):", paste(dim(spatial_count_common), collapse="x")))
  print(paste("Dimensions of spatial_location_for_giotto (for ST object):", paste(dim(spatial_location_for_giotto), collapse="x")))

  print("Processing single-cell data for signature matrix...")
  giotto_SC <- createGiottoObject(expression = sc_count_common, cell_metadata = current_sc_meta, instructions = instrs)
  giotto_SC <- normalizeGiotto(giotto_SC, scalefactor = scale_factor, verbose = FALSE)

  print("Finding markers on common genes for single-cell data...")
  markers_scran <- findMarkers_one_vs_all(gobject=giotto_SC, method="scran", expression_values="normalized", cluster_column = "cellType", min_feats=3)
  if(is.null(markers_scran) || nrow(markers_scran) == 0) stop("No SC marker genes found.")
  top_markers_dt <- markers_scran[, head(.SD, top_n_markers), by="cluster"]
  unique_signature_genes <- unique(top_markers_dt$feats)
  if(length(unique_signature_genes) == 0) stop("No unique SC signature genes selected.")
  print(paste("Number of unique signature genes for DWLS (derived from common genes):", length(unique_signature_genes)))

  sc_normalized_matrix_for_sig <- getExpression(giotto_SC, values = "normalized", output = "matrix")
  unique_signature_genes <- intersect(unique_signature_genes, rownames(sc_normalized_matrix_for_sig))
  if(length(unique_signature_genes) == 0) stop("Signature genes not found in SC normalized matrix after marker selection intersection.")

  sc_norm_expr_subset_for_sig <- sc_normalized_matrix_for_sig[unique_signature_genes, , drop = FALSE]

  if(any(!is.finite(as.matrix(sc_norm_expr_subset_for_sig)))){
    warning("Non-finite values (NA, NaN, Inf) found in sc_norm_expr_subset_for_sig. Replacing with 0.")
    sc_norm_expr_subset_for_sig[which(!is.finite(sc_norm_expr_subset_for_sig), arr.ind = TRUE)] <- 0
  }

  DWLS_matrix <- makeSignMatrixDWLSfromMatrix(matrix = sc_norm_expr_subset_for_sig,
                                            cell_type = pDataDT(giotto_SC)$cellType,
                                            sign_gene = unique_signature_genes)
  if(is.null(DWLS_matrix) || nrow(DWLS_matrix) == 0 || ncol(DWLS_matrix) == 0) stop("Failed to create DWLS signature matrix.")

  if(any(!is.finite(as.matrix(DWLS_matrix)))){
    warning("Non-finite values (NA, NaN, Inf) found in DWLS_matrix after creation. Replacing with 0.")
    DWLS_matrix[which(!is.finite(DWLS_matrix), arr.ind = TRUE)] <- 0
  }

  DWLS_matrix <- DWLS_matrix[unique_signature_genes, , drop = FALSE]

  all_zero_cols <- apply(DWLS_matrix, 2, function(col) all(col == 0))
  if(any(all_zero_cols)){
    warning(paste("Signature matrix has columns (cell types) that are all zero:",
                  paste(colnames(DWLS_matrix)[all_zero_cols], collapse=", ")))
    DWLS_matrix <- DWLS_matrix[, !all_zero_cols, drop = FALSE]
    if(ncol(DWLS_matrix) == 0) stop("Signature matrix empty after removing all-zero columns.")
    print(paste("Signature matrix reduced to", ncol(DWLS_matrix), "cell types after removing all-zero columns."))
  }

  if(is.null(DWLS_matrix) || nrow(DWLS_matrix) == 0 || ncol(DWLS_matrix) == 0) stop("DWLS signature matrix became empty after cleaning.")
  print("DWLS signature matrix created and cleaned.")

  print("Processing ST data: Normalization and Clustering (Giotto filtering on ST data is SKIPPED).")
  spot_sums_initial <- Matrix::colSums(spatial_count_common)
  zero_spots_indices <- which(spot_sums_initial == 0)
  if(length(zero_spots_indices) > 0) {
    print(paste("Removing", length(zero_spots_indices), "spots with zero total counts from ST data (common gene set)."))
    spatial_count_common <- spatial_count_common[, -zero_spots_indices, drop = FALSE]
    spatial_location_for_giotto <- spatial_location_for_giotto[colnames(spatial_count_common), , drop = FALSE]
  }
  if(ncol(spatial_count_common) == 0) {
    print("All ST spots had zero counts on common genes or were removed. Skipping.")
    empty_results <- data.table(cell_ID = character(0), Note = "All ST spots removed (zero counts on common genes)")
    write.csv(empty_results, output_path, row.names = FALSE)
    return(NULL)
  }

  giotto_ST <- createGiottoObject(expression = spatial_count_common,
                                  spatial_locs = spatial_location_for_giotto,
                                  instructions = instrs)
  print("Giotto 'filterGiotto' step on ST data is SKIPPED.")
  print("Normalizing ST data...")
  giotto_ST <- normalizeGiotto(giotto_ST, scalefactor = scale_factor, verbose = FALSE)
  print("Adding statistics to ST data...")
  giotto_ST <- addStatistics(gobject = giotto_ST, expression_values = "normalized", verbose = FALSE)

  # --- HVF CALCULATION REMOVED ---
  print("HVF calculation is SKIPPED. Using all common genes for PCA.")
  # Ensure 'common_genes' are those present in the ST Giotto object's expression matrix
  # The `spatial_count_common` already subsetted to common_genes, so `rownames(giotto_ST@expression$rna$normalized)` will be these.
  features_for_pca <- rownames(getExpression(giotto_ST, values="normalized", output="matrix"))
  if(is.null(features_for_pca) || length(features_for_pca) == 0) {
      stop("No features available for PCA. This means the ST normalized matrix is empty.")
  }
  print(paste("Using", length(features_for_pca), "genes (all common genes in ST normalized matrix) for PCA."))


  num_feats_for_pca_calc <- length(features_for_pca)
  num_cells_for_pca_calc <- ncol(getExpression(giotto_ST, values="normalized", output="matrix"))

  if(num_cells_for_pca_calc <=1 || num_feats_for_pca_calc <=1) {
    print("Not enough cells or features remaining for PCA. Skipping deconvolution.")
    empty_results <- data.table(cell_ID = character(0), Note = "Not enough data for PCA")
    write.csv(empty_results, output_path, row.names = FALSE)
    return(NULL)
  }
  max_allowable_pcs <- min(num_cells_for_pca_calc -1, num_feats_for_pca_calc); if (max_allowable_pcs < 1) max_allowable_pcs <- 1
  capped_pca_dims <- min(pca_dims_to_use, 50, max_allowable_pcs)
  if(capped_pca_dims < 1) {
      if (max_allowable_pcs >=1) {
          warning("PCA dimensions capped at 1. Clustering may not be meaningful.")
          actual_pca_dims_to_use <- 1
      } else {
        print("PCA cannot be run (0 dimensions). Skipping deconvolution.")
        empty_results <- data.table(cell_ID = character(0), Note = "PCA cannot be run, 0 dimensions")
        write.csv(empty_results, output_path, row.names = FALSE)
        return(NULL)
      }
  } else {
      actual_pca_dims_to_use <- capped_pca_dims
  }

  print(paste("Running PCA on ST data using up to", actual_pca_dims_to_use, "dimensions, on", length(features_for_pca), "features."))
  giotto_ST <- runPCA(gobject = giotto_ST, expression_values = "normalized", feats_to_use = features_for_pca, ncp = actual_pca_dims_to_use, verbose = FALSE)
  pca_results <- getDimReduction(giotto_ST, reduction = 'cells', reduction_method = 'pca', name = 'pca', output = 'matrix')
  if(is.null(pca_results) || ncol(pca_results) < 1) {
      print(paste("PCA did not produce any dimensions. Skipping deconvolution."))
      empty_results <- data.table(cell_ID = character(0), Note = "PCA failed or produced 0 dimensions")
      write.csv(empty_results, output_path, row.names = FALSE)
      return(NULL)
  }
  actual_pca_dims_obtained <- ncol(pca_results)
  if(actual_pca_dims_obtained < 1) {
      print("PCA resulted in 0 usable dimensions. Skipping deconvolution.")
      empty_results <- data.table(cell_ID = character(0), Note = "PCA resulted in 0 usable dimensions")
      write.csv(empty_results, output_path, row.names = FALSE)
      return(NULL)
  }


  print(paste("Creating nearest neighbor network on ST data using k =", nn_k_val, "and top", actual_pca_dims_obtained, "PCs."))
  actual_nn_k_val <- min(nn_k_val, num_cells_for_pca_calc - 1);
  if(actual_nn_k_val < 1 && num_cells_for_pca_calc > 1) {
      actual_nn_k_val <- 1
  } else if (num_cells_for_pca_calc <= 1 || actual_nn_k_val < 1) {
      print(paste("Cannot create NN, k < 1 or not enough cells (", num_cells_for_pca_calc, "). Skipping deconvolution."))
      empty_results <- data.table(cell_ID = character(0), Note = "Cannot create NN, k < 1 or too few cells")
      write.csv(empty_results, output_path, row.names = FALSE)
      return(NULL)
  }
  if(actual_nn_k_val != nn_k_val) print(paste("Adjusted k for NN from", nn_k_val, "to", actual_nn_k_val))
  giotto_ST <- createNearestNetwork(gobject = giotto_ST, type = 'sNN', dimensions_to_use = 1:actual_pca_dims_obtained, k = actual_nn_k_val, name = "sNN.network", verbose = FALSE)

  st_actual_cluster_colname <- "st_leiden_cluster"
  print(paste("Running Leiden clustering on ST data with resolution =", resolution))
  giotto_ST <- doLeidenCluster(gobject = giotto_ST, name = st_actual_cluster_colname, resolution = resolution, n_iterations = 1000, nn_network_to_use = 'sNN', network_name = 'sNN.network', verbose = FALSE)
  if(!st_actual_cluster_colname %in% colnames(pDataDT(giotto_ST))) { stop(paste("Actual ST cluster column '",st_actual_cluster_colname,"' not added.",sep="")) }
  num_st_clusters_found <- length(unique(pDataDT(giotto_ST)[[st_actual_cluster_colname]]))
  print(paste("Found", num_st_clusters_found, "ST Leiden clusters."))
  if(num_st_clusters_found == 0 ) stop("Leiden clustering resulted in 0 clusters.")
  if(num_st_clusters_found <= 1 && num_cells_for_pca_calc > 10) { warning("Leiden clustering resulted in <=1 cluster.") }

  sign_matrix_to_use <- DWLS_matrix
  print(paste("Running DWLS Deconvolution using ST actual cluster column:", st_actual_cluster_colname))
  print(paste("DWLS Signature matrix for deconv has", nrow(sign_matrix_to_use), "genes and", ncol(sign_matrix_to_use), "cell types."))

  st_norm_expr_genes_for_deconv <- rownames(getExpression(giotto_ST, values='normalized', output='matrix'))
  print(paste("ST object (normalized) for deconv has", length(st_norm_expr_genes_for_deconv), "genes."))

  common_genes_final_deconv <- intersect(rownames(sign_matrix_to_use), st_norm_expr_genes_for_deconv)
  if(length(common_genes_final_deconv) == 0){
    stop("No common genes between final signature matrix and ST normalized expression for deconvolution.")
  }
  if(length(common_genes_final_deconv) < nrow(sign_matrix_to_use)){
    warning(paste("Signature matrix has", nrow(sign_matrix_to_use), "genes, but only", length(common_genes_final_deconv),
                  "are shared with ST normalized expression. Subsetting signature matrix to these common genes."))
    sign_matrix_to_use <- sign_matrix_to_use[common_genes_final_deconv, , drop = FALSE]
    if(nrow(sign_matrix_to_use) == 0 || ncol(sign_matrix_to_use) == 0) stop("Signature matrix empty after final intersection with ST genes for deconvolution.")
  }
  print(paste("Final signature matrix for DWLS has", nrow(sign_matrix_to_use), "genes and", ncol(sign_matrix_to_use), "cell types."))

  if(any(!is.finite(as.matrix(sign_matrix_to_use)))) {
      stop("Non-finite values found in final signature matrix before DWLS. This should not happen.")
  }
  if(nrow(sign_matrix_to_use) == 0 || ncol(sign_matrix_to_use) == 0){
      stop("Final signature matrix is empty before DWLS.")
  }


  giotto_ST <- runDWLSDeconv(
        gobject = giotto_ST,
        sign_matrix = sign_matrix_to_use,
        expression_values = "normalized",
        cluster_column = st_actual_cluster_colname,
        return_gobject = TRUE
    )

  print("DWLS Deconvolution completed successfully.")
  deconv_results <- tryCatch({
      getSpatialEnrichment(gobject = giotto_ST, name = "DWLS", output = "data.table")
  }, error = function(e_get) {
      warning(paste("Error in getSpatialEnrichment for DWLS:", e_get$message))
      return(NULL)
  })

  if(is.null(deconv_results) || nrow(deconv_results) == 0) {
    warning("DWLS ran but getSpatialEnrichment returned no results (empty table or error).")
    deconv_results <- data.table(cell_ID = character(0))
    expected_ct_cols <- if(exists("sign_matrix_to_use") && !is.null(colnames(sign_matrix_to_use))) colnames(sign_matrix_to_use) else character(0)
    for(ct_name in expected_ct_cols) {
        deconv_results[, (ct_name) := numeric(0)]
    }
    deconv_results[, Note := "DWLS getSpatialEnrichment returned 0 rows or failed"]
  }

  write.csv(deconv_results, output_path, row.names = FALSE)
  print(paste("Deconvolution results saved to:", output_path))
  print("run_spatial_dwls function completed.")
  return(deconv_results)
}

print("Attempting to run deconvolution script...")
tryCatch({
    run_spatial_dwls(
      st_dir = opt$st_dir, sc_dir = opt$sc_dir, cell_type_path = opt$cell_type_path,
      output_path = opt$output_path, resolution = opt$resolution,
      scale_factor = opt$scale_factor, top_n_markers = opt$markers,
      # Removed hvf_method and skip_hvf from call
      pca_dims_to_use = opt$pca_dims, nn_k_val = opt$nn_k
    )
    print("Deconvolution script finished successfully.")
}, error = function(e) {
    print(paste("Error during script execution: ", e$message))
    if(exists("sys.calls")) { print("Traceback (simplified):"); print(head(sys.calls(), n=10L)) }
    print("Session Info:")
    print(sessionInfo())
    if(exists("giotto_ST")) save(giotto_ST, file = "giotto_ST_at_error.RData")
    if(exists("sign_matrix_to_use")) save(sign_matrix_to_use, file = "sign_matrix_to_use_at_error.RData")
    if(exists("DWLS_matrix")) save(DWLS_matrix, file = "DWLS_matrix_initial_at_error.RData")
    stop(e)
})

print("Script execution finished.")
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
  make_option(c("-r", "--resolution"), type="numeric", default=0.5,  help="Resolution for Leiden clustering (NOT USED in this version for ST) [default=%default]"),
  make_option(c("-f", "--scale_factor"), type="numeric", default=10000, help="Scale factor for normalization [default=%default]"),
  make_option(c("-m", "--markers"), type="numeric", default=10, help="Number of top markers per cell type for signature matrix [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list, description="Perform spatial deconvolution using DWLS method (with a dummy single cluster for ST data)")
opt <- parse_args(opt_parser)

if (is.null(opt$st_dir) || is.null(opt$sc_dir) || is.null(opt$cell_type_path)) {
  print_help(opt_parser)
  stop("Missing required arguments: st_dir, sc_dir, and cell_type_path must be provided", call.=FALSE)
}

run_spatial_dwls <- function(st_dir, sc_dir, cell_type_path, output_path, resolution, scale_factor, top_n_markers) {
  print("Starting run_spatial_dwls function...")

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
  cell_type_map <- sc_cell_types$cellType; names(cell_type_map) <- sc_cell_types$cell_ID
  sc_meta <- data.frame(cell_ID = sc_cell_ids, cellType = cell_type_map[sc_cell_ids], sampleInfo = "sample1", stringsAsFactors = FALSE)
  rownames(sc_meta) <- sc_cell_ids
  sc_meta <- sc_meta[!is.na(sc_meta$cellType), , drop=FALSE]
  if(nrow(sc_meta) == 0) stop("No cells remaining in single-cell metadata after filtering for cell types.")
  sc_cell_ids_filtered <- rownames(sc_meta)
  sc_sparse_matrix <- readMM(file.path(sc_dir, "sparse_matrix.mtx"))
  dim_sc_genes <- min(nrow(sc_sparse_matrix), length(sc_gene_ids)); dim_sc_cells <- min(ncol(sc_sparse_matrix), length(sc_cell_ids))
  sc_count_pre <- sc_sparse_matrix[1:dim_sc_genes, 1:dim_sc_cells, drop=FALSE]
  rownames(sc_count_pre) <- sc_gene_ids[1:dim_sc_genes]; colnames(sc_count_pre) <- sc_cell_ids[1:dim_sc_cells]
  sc_count <- sc_count_pre[, intersect(colnames(sc_count_pre), sc_cell_ids_filtered), drop=FALSE]
  if(ncol(sc_count) == 0 || nrow(sc_count) == 0) stop("SC count matrix is empty after filtering.")

  print("Reading and preparing spatial data...")
  st_cell_info <- fread(file.path(st_dir, "cellinfo.csv"), skip = 1, header = FALSE); st_gene_info <- fread(file.path(st_dir, "geneinfo.csv"), skip = 1, header = FALSE)
  st_cell_ids <- st_cell_info[[1]]; st_gene_ids <- st_gene_info[[1]]
  st_sparse_matrix <- readMM(file.path(st_dir, "sparse_matrix.mtx"))
  dim_st_genes <- min(nrow(st_sparse_matrix), length(st_gene_ids)); dim_st_cells <- min(ncol(st_sparse_matrix), length(st_cell_ids))
  spatial_count_pre <- st_sparse_matrix[1:dim_st_genes, 1:dim_st_cells, drop=FALSE]
  rownames(spatial_count_pre) <- st_gene_ids[1:dim_st_genes]; colnames(spatial_count_pre) <- st_cell_ids[1:dim_st_cells]
  spatial_count <- spatial_count_pre
  if (ncol(st_cell_info) >= 3 && is.numeric(st_cell_info[[2]]) && is.numeric(st_cell_info[[3]])) {
    spatial_location <- data.frame(x = st_cell_info[[2]], y = st_cell_info[[3]], row.names = st_cell_ids, stringsAsFactors = FALSE)
    spatial_location <- spatial_location[intersect(rownames(spatial_location), colnames(spatial_count)), , drop = FALSE]
    print("Using x and y coordinates from st/cellinfo.csv")
  } else {
    print("Generating placeholder grid coordinates for spatial data.")
    n_spots_loc <- ncol(spatial_count); n_cols_grid_loc <- ceiling(sqrt(n_spots_loc))
    spatial_location <- data.frame(x = (seq_len(n_spots_loc) - 1) %% n_cols_grid_loc, y = (seq_len(n_spots_loc) - 1) %/% n_cols_grid_loc, row.names = colnames(spatial_count), stringsAsFactors = FALSE)
  }
  if(nrow(spatial_location) != ncol(spatial_count)) stop("Mismatch between spatial locations and spatial_count spots.")

  common_genes <- intersect(rownames(sc_count), rownames(spatial_count))
  if(length(common_genes) == 0) stop("No common genes found.")
  print(paste("Number of common genes:", length(common_genes)))
  sc_count <- sc_count[common_genes, , drop=FALSE]; spatial_count <- spatial_count[common_genes, , drop=FALSE]
  print(paste("Dimensions of sc_count (common genes, filtered cells):", paste(dim(sc_count), collapse="x")))
  print(paste("Dimensions of spatial_count (common genes):", paste(dim(spatial_count), collapse="x")))

  print("Processing single-cell data for signature matrix...")
  current_sc_meta <- sc_meta[colnames(sc_count), , drop=FALSE]
  giotto_SC <- createGiottoObject(expression = sc_count, cell_metadata = current_sc_meta)
  giotto_SC <- normalizeGiotto(giotto_SC, scalefactor = scale_factor)
  print("Finding markers on all (common) genes for single-cell data...")
  markers_scran <- findMarkers_one_vs_all(gobject=giotto_SC, method="scran", expression_values="normalized", cluster_column = "cellType", min_feats=3, feats_to_use = NULL)
  if(is.null(markers_scran) || nrow(markers_scran) == 0) stop("No SC marker genes found.")
  top_markers_dt <- markers_scran[, head(.SD, top_n_markers), by="cluster"]
  unique_signature_genes <- unique(top_markers_dt$feats)
  if(length(unique_signature_genes) == 0) stop("No unique SC signature genes selected.")
  print(paste("Number of unique signature genes for DWLS (from all SC genes):", length(unique_signature_genes)))
  DWLS_matrix <- makeSignMatrixDWLSfromMatrix(matrix = getExpression(giotto_SC, values = "normalized", output = "matrix"), cell_type = pDataDT(giotto_SC)$cellType, sign_gene = unique_signature_genes)
  if(is.null(DWLS_matrix) || nrow(DWLS_matrix) == 0 || ncol(DWLS_matrix) == 0) stop("Failed to create DWLS signature matrix.")
  print("DWLS signature matrix created.")

  print("Processing ST data: Filtering, Normalization, and adding dummy cluster...")
  spot_sums <- Matrix::colSums(spatial_count)
  zero_spots_indices <- which(spot_sums == 0)
  if(length(zero_spots_indices) > 0) {
    print(paste("Removing", length(zero_spots_indices), "spots with zero total counts from ST data."))
    spatial_count <- spatial_count[, -zero_spots_indices, drop = FALSE]
    spatial_location <- spatial_location[colnames(spatial_count), , drop = FALSE]
  }
  if(ncol(spatial_count) == 0) stop("All ST spots removed. Deconvolution cannot proceed.")
  st_cell_ids_filtered_st <- colnames(spatial_count) # These are the final spot IDs to use

  giotto_ST <- createGiottoObject(expression = spatial_count, spatial_locs = spatial_location)
  # Add initial metadata (e.g., sample_ID)
  initial_st_meta_df <- data.frame(sample_ID = rep("sample1", length(st_cell_ids_filtered_st)), 
                                   row.names = st_cell_ids_filtered_st, 
                                   stringsAsFactors = FALSE)
  giotto_ST <- addCellMetadata(giotto_ST, new_metadata = initial_st_meta_df, by_column = FALSE)
  
  giotto_ST <- normalizeGiotto(giotto_ST, scalefactor = scale_factor)

  # Create and add a dummy cluster column where all spots belong to one cluster
  st_dummy_cluster_colname <- "all_spots_in_one_cluster"
  print(paste("Creating dummy cluster column '", st_dummy_cluster_colname, "' for ST data.", sep=""))
  
  # Create a new data.table for the dummy cluster metadata. 
  # It must have a column with cell_IDs that match those in giotto_ST.
  dummy_cluster_metadata <- data.table(
    cell_ID = st_cell_ids_filtered_st, # Use the filtered spot IDs
    temp_cluster_col = factor(rep("global_cluster", length(st_cell_ids_filtered_st)))
  )
  # Rename the temporary column to the desired final column name
  setnames(dummy_cluster_metadata, "temp_cluster_col", st_dummy_cluster_colname)
  
  giotto_ST <- addCellMetadata(
    gobject = giotto_ST,
    new_metadata = dummy_cluster_metadata,
    by_column = TRUE, # TRUE because we are matching based on the 'cell_ID' column in dummy_cluster_metadata
    column_cell_ID = "cell_ID" # The name of the column in dummy_cluster_metadata that holds the cell/spot IDs
  )
  
  print(paste("Dummy cluster column '", st_dummy_cluster_colname, "' added to ST metadata.", sep=""))
  if(!st_dummy_cluster_colname %in% colnames(pDataDT(giotto_ST))) {
      stop(paste("Dummy cluster column '", st_dummy_cluster_colname, "' was not successfully added.", sep=""))
  }

  # --- Run DWLS Deconvolution using the dummy ST cluster column ---
  print(paste("Running DWLS Deconvolution using ST dummy cluster column:", st_dummy_cluster_colname))
  giotto_ST <- runDWLSDeconv(
        gobject = giotto_ST, 
        sign_matrix = DWLS_matrix, 
        expression_values = "normalized",
        cluster_column = st_dummy_cluster_colname, 
        return_gobject = TRUE
    )
  
  print("DWLS Deconvolution completed successfully.")
  deconv_results <- getSpatialEnrichment(gobject = giotto_ST, name = "DWLS", output = "data.table")
  if(is.null(deconv_results) || nrow(deconv_results) == 0) {
    stop("DWLS ran but getSpatialEnrichment returned no results.")
  }

  # --- Save Results ---
  write.csv(deconv_results, output_path, row.names = FALSE)
  print(paste("Deconvolution results saved to:", output_path))
  print("run_spatial_dwls function completed.")
  return(deconv_results)
}

print("Attempting to run deconvolution script...")
run_spatial_dwls(
  st_dir = opt$st_dir, sc_dir = opt$sc_dir, cell_type_path = opt$cell_type_path,
  output_path = opt$output_path, resolution = opt$resolution,
  scale_factor = opt$scale_factor, top_n_markers = opt$markers
)
print("Deconvolution script finished successfully (if this message is reached).")
print("Script execution finished.")
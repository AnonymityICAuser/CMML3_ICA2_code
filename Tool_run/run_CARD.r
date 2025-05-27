library(Matrix)
library(data.table)
library(CARD)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-s", "--st_dir"), type="character", default=NULL,
              help="Directory containing spatial transcriptomics data"),
  make_option(c("-c", "--sc_dir"), type="character", default=NULL,
              help="Directory containing single-cell data"),
  make_option(c("-t", "--cell_type_file"), type="character", default=NULL,
              help="Path to cell type annotation file"),
  make_option(c("-o", "--output"), type="character", default="cell_type_proportions.csv",
              help="Output file path for cell type proportions [default= %default]"),
  make_option(c("-m", "--min_count_gene"), type="integer", default=1,
              help="Minimum count per gene [default= %default]"),
  make_option(c("-p", "--min_count_spot"), type="integer", default=1,
              help="Minimum count per spot [default= %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$st_dir) || is.null(opt$sc_dir) || is.null(opt$cell_type_file)) {
  print_help(opt_parser)
  stop("Missing required arguments. Please provide st_dir, sc_dir, and cell_type_file.")
}

# Read cell type information
sc_cell_types <- fread(opt$cell_type_file, header=TRUE)

# Rename columns if needed
if("cell_label" %in% colnames(sc_cell_types) && "class" %in% colnames(sc_cell_types)) {
  setnames(sc_cell_types, c("cell_label", "class"), c("cell_ID", "cellType"))
}

# Read cell and gene information
st_cell_info <- fread(file.path(opt$st_dir, "cellinfo.csv"), skip = 1, header = FALSE)
sc_cell_info <- fread(file.path(opt$sc_dir, "cellinfo.csv"), skip = 1, header = FALSE)
st_gene_info <- fread(file.path(opt$st_dir, "geneinfo.csv"), skip = 1, header = FALSE)
sc_gene_info <- fread(file.path(opt$sc_dir, "geneinfo.csv"), skip = 1, header = FALSE)

# Extract IDs
st_cell_ids <- st_cell_info[[1]]
sc_cell_ids <- sc_cell_info[[1]]
st_gene_ids <- st_gene_info[[1]]
sc_gene_ids <- sc_gene_info[[1]]

# Create sc_meta with proper alignment
sc_meta <- data.frame(
  cellID = sc_cell_ids,
  cellType = sc_cell_types$cellType[match(sc_cell_ids, sc_cell_types$cell_ID)],
  sampleInfo = "sample1",
  stringsAsFactors = FALSE
)
rownames(sc_meta) <- sc_cell_ids

# Read the sparse matrices
st_sparse_matrix <- readMM(file.path(opt$st_dir, "sparse_matrix.mtx"))
sc_sparse_matrix <- readMM(file.path(opt$sc_dir, "sparse_matrix.mtx"))

# For spatial coordinates
n_cols <- ceiling(sqrt(length(st_cell_ids)))
st_x_coords <- (seq_along(st_cell_ids) - 1) %% n_cols
st_y_coords <- (seq_along(st_cell_ids) - 1) %/% n_cols

# Create spatial_location
spatial_location <- data.frame(
  x = st_x_coords,
  y = st_y_coords,
  row.names = st_cell_ids,
  stringsAsFactors = FALSE
)

# Handle matrix size differences
if (nrow(st_sparse_matrix) == length(st_gene_ids) + 1) {
  st_sparse_matrix <- st_sparse_matrix[-1, ]
}
if (ncol(st_sparse_matrix) == length(st_cell_ids) + 1) {
  st_sparse_matrix <- st_sparse_matrix[, -1]
}

if (nrow(sc_sparse_matrix) == length(sc_gene_ids) + 1) {
  sc_sparse_matrix <- sc_sparse_matrix[-1, ]
}
if (ncol(sc_sparse_matrix) == length(sc_cell_ids) + 1) {
  sc_sparse_matrix <- sc_sparse_matrix[, -1]
}

# Assign matrices and set row/column names
sc_count <- sc_sparse_matrix
rownames(sc_count) <- sc_gene_ids[1:nrow(sc_count)]
colnames(sc_count) <- sc_cell_ids[1:ncol(sc_count)]

spatial_count <- st_sparse_matrix
rownames(spatial_count) <- st_gene_ids[1:nrow(spatial_count)]
colnames(spatial_count) <- st_cell_ids[1:ncol(spatial_count)]

# Find common genes
common_genes <- intersect(rownames(sc_count), rownames(spatial_count))

# Filter matrices to only include common genes
sc_count <- sc_count[common_genes, ]
spatial_count <- spatial_count[common_genes, ]

print(dim(sc_count))
print(dim(spatial_count))

cat(sprintf("Dimensions of sc_count: %d genes × %d cells\n", nrow(sc_count), ncol(sc_count)))
cat(sprintf("Dimensions of spatial_count: %d genes × %d spots\n", nrow(spatial_count), ncol(spatial_count)))

# Create CARD object
CARD_obj <- createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 0,
  minCountSpot = 0
)

CARD_obj <- CARD_deconvolution(CARD_obj)

# Extract proportions and save
proportions <- CARD_obj@Proportion_CARD
write.csv(proportions, opt$output)
cat(sprintf("Results saved to: %s\n", opt$output))
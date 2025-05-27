
# Base directory for your processed data pairs
DATA_ROOT_DIR="Simulation_merfish"

# Base output directory
OUTPUT_DIR_BASE="Simulation_merfish/results_data_pairs"

# Base script directory (assuming your scripts are in this structure relative to where you submit from)
SCRIPT_BASE_DIR="scripts" # Adjust if your scripts are elsewhere

# Create base output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR_BASE}"

# Array of descriptive names for each pair (used for output file naming)
# These are derived from the MERFISH filenames.
declare -a PAIR_NAMES=(
   # Add your pair names here
)

# Corresponding ST H5AD files (Spatial Transcriptomics data for cell2loc or similar tools)
# These files are located in DATA_ROOT_DIR and end with .h5ad.h5ad
declare -a ST_H5AD_FILES=(
  # Add your ST H5AD files here
)

# Corresponding SC H5AD files for cell2loc (Single-Cell data from reference_adata directory)
declare -a SC_H5AD_FILES=(
  # Add your SC H5AD files here
)

# Corresponding SC directories for R scripts (Single-Cell data from for_R directory)
declare -a SC_FOR_R_DIRS=(
    # Add your SC R-dir files here
)

# Corresponding ST files for R scripts (Spatial Transcriptomics data formatted for R)
# These files are located in DATA_ROOT_DIR/for_R and end with .h5ad
# Note: Original variable name was ST_FOR_R_DIRS; these are files.
# Ensure your R scripts can handle these paths correctly.
declare -a ST_FOR_R_DIRS=( # Renamed from ST_FOR_R_DIRS to ST_FOR_R_FILES for clarity
  # Add your ST R-files here
)

NUM_PAIRS=${#PAIR_NAMES[@]}

# Note: If you plan to run each pair as a separate SLURM array task,
# you would modify this loop to use ${SLURM_ARRAY_TASK_ID}.
# For a single job processing all pairs sequentially, this loop is fine.

for (( i=0; i<${NUM_PAIRS}; i++ )); do
    CURRENT_PAIR_NAME="${PAIR_NAMES[$i]}"

    # Specific output directory for this pair
    CURRENT_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${CURRENT_PAIR_NAME}"
    mkdir -p "${CURRENT_OUTPUT_DIR}"

    # Define paths for the current pair
    SC_H5AD_PATH="${SC_H5AD_FILES[$i]}"
    ST_H5AD_PATH="${ST_H5AD_FILES[$i]}"
    SC_R_DIR_PATH="${SC_FOR_R_DIRS[$i]}"
    ST_R_DIR_PATH="${ST_FOR_R_DIRS[$i]}"

    # Cell type file for R scripts (assumed to be cellinfo.csv in the SC_R_DIR_PATH)
    CELL_TYPE_FILE_R="${SC_R_DIR_PATH}/cellinfo.csv"

    # Output paths for tools
    CARD_OUTPUT="${CURRENT_OUTPUT_DIR}/CARD_proportions.csv"
    DWLS_OUTPUT="${CURRENT_OUTPUT_DIR}/DWLS_proportions.csv"
    CELL2LOC_OUTPUT="${CURRENT_OUTPUT_DIR}/cell2loc_proportions.csv"

    # --- Sanity checks for input files/dirs ---
    echo "Checking input files for ${CURRENT_PAIR_NAME}:"
    error_found=0
    if [ ! -f "${SC_H5AD_PATH}" ]; then echo "  ERROR: SC H5AD not found: ${SC_H5AD_PATH}"; error_found=1; fi
    if [ ! -f "${ST_H5AD_PATH}" ]; then echo "  ERROR: ST H5AD not found: ${ST_H5AD_PATH}"; error_found=1; fi
    if [ ! -d "${SC_R_DIR_PATH}" ]; then echo "  ERROR: SC R-dir not found: ${SC_R_DIR_PATH}"; error_found=1; fi
    if [ ! -f "${CELL_TYPE_FILE_R}" ]; then echo "  ERROR: Cell type file for R not found: ${CELL_TYPE_FILE_R}"; error_found=1; fi
    if [ ! -d "${ST_R_DIR_PATH}" ]; then echo "  ERROR: ST R-dir not found: ${ST_R_DIR_PATH}"; error_found=1; fi

    if [ ${error_found} -eq 1 ]; then
        echo "  Skipping ${CURRENT_PAIR_NAME} due to missing input files."
        continue # Skip to the next iteration of the loop
    fi
    echo "  All input files/dirs seem to exist for ${CURRENT_PAIR_NAME}."


    echo "----------------------------------------------------"
    echo "Running CARD for ${CURRENT_PAIR_NAME}..."
    echo "Start time: $(date)"
    conda run -n r_env Rscript "${SCRIPT_BASE_DIR}/deconv/run_CARD.r" \
        --st_dir "${ST_R_DIR_PATH}" \
        --sc_dir "${SC_R_DIR_PATH}" \
        --cell_type_file "${CELL_TYPE_FILE_R}" \
        --output "${CARD_OUTPUT}"
    if [ $? -ne 0 ]; then echo "ERROR: CARD failed for ${CURRENT_PAIR_NAME}"; else echo "CARD completed."; fi
    echo "End time: $(date)"
    echo "----------------------------------------------------"

    # --- Run spatialDWLS using r_env ---
    echo "----------------------------------------------------"
    echo "Running spatialDWLS for ${CURRENT_PAIR_NAME}..."
    echo "Start time: $(date)"
    conda run -n r_env Rscript "${SCRIPT_BASE_DIR}/deconv/run_spatalDWLS.r" \
        --st_dir "${ST_R_DIR_PATH}" \
        --sc_dir "${SC_R_DIR_PATH}" \
        --cell_type_path "${CELL_TYPE_FILE_R}" \
        --output_path "${DWLS_OUTPUT}"
    if [ $? -ne 0 ]; then echo "ERROR: spatialDWLS failed for ${CURRENT_PAIR_NAME}"; else echo "spatialDWLS completed."; fi
    echo "End time: $(date)"
    echo "----------------------------------------------------"

    # --- Run cell2loc using cell2loc_env ---
    echo "----------------------------------------------------"
    echo "Running cell2loc for ${CURRENT_PAIR_NAME}..."
    echo "Start time: $(date)"
    conda run -n cell2loc_env python "${SCRIPT_BASE_DIR}/deconv/run_cell2loc.py" \
        --adata_ref "${SC_H5AD_PATH}" \
        --adata_decon "${ST_H5AD_PATH}" \
        --output_csv "${CELL2LOC_OUTPUT}"
    if [ $? -ne 0 ]; then echo "ERROR: cell2loc failed for ${CURRENT_PAIR_NAME}"; else echo "cell2loc completed."; fi
    echo "End time: $(date)"
    echo "----------------------------------------------------"

    echo "Completed processing for data pair: ${CURRENT_PAIR_NAME}"
    echo "" # Add a blank line for better readability
done

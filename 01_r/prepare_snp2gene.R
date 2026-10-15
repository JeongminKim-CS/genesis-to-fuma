# ============================================================
# Step 1. Configuration
# ============================================================

# Input and output directories
input_dir <- "03_data"
output_dir <- "02_output/snp2gene"

# P-value threshold (configurable)
p_threshold <- 1e-3

# Generate P-value threshold label
p_label <- format(p_threshold, scientific = TRUE)
p_label <- sub("e-0", "e-", p_label)
p_label <- sub("e\\+0?", "e", p_label)

# Genome build
genome_build <- "GRCh38"

# Required GENESIS columns
required_cols <- c(
  "CHR",
  "POS",
  "EFFECT_ALLELE",
  "OTHER_ALLELE",
  "BETA",
  "PVAL"
)

# Validate input directory
if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}

# Create output directory
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Step 2. Discover GENESIS CSV Files
# ============================================================

csv_files <- list.files(
  path = input_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Stop if no CSV files are found
if (length(csv_files) == 0) {
  stop("No CSV files found in: ", input_dir)
}

# Display discovered files
cat("CSV files found:", length(csv_files), "\n")
print(basename(csv_files))

# ============================================================
# Step 3. Process GENESIS CSV Files
# ============================================================

for (csv_file in csv_files) {
  
  # Extract dataset name
  csv_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Keep only gene ID and phenotype
  dataset_name <- sub(
    "^([^_]+_[^_]+).*",
    "\\1",
    csv_name
  )
  
  # Load CSV
  genesis_data <- read.csv(
    csv_file,
    check.names = FALSE
  )
  
  # Check required columns
  missing_cols <- setdiff(required_cols, names(genesis_data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ", basename(csv_file),
      ": ", paste(missing_cols, collapse = ", ")
    )
  }
  
  # Filter SNPs
  genesis_filtered <- genesis_data[
    !is.na(genesis_data$PVAL) &
      genesis_data$PVAL <= p_threshold,
  ]
  
  # Select FUMA columns
  fuma_data <- genesis_filtered[, required_cols]
  
  # Rename columns
  colnames(fuma_data) <- c(
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "p_value"
  )
  
  # Check filtered data
  if (nrow(fuma_data) == 0) {
    warning("No SNPs passed threshold: ", dataset_name)
    next
  }
  
  stopifnot(
    all(!is.na(fuma_data$p_value)),
    all(fuma_data$p_value <= p_threshold)
  )
  
  output_name <- paste0(
    dataset_name,
    "_FUMA_snp2gene_p",
    p_label,
    "_",
    genome_build,
    ".txt"
  )
  
  output_file <- file.path(output_dir, output_name)
  
  # Export TXT
  write.table(
    fuma_data,
    file = output_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  
  cat(
    "Created:", output_name,
    "| SNPs:", nrow(fuma_data),
    "\n"
  )
  
} # End of for loop

# ============================================================
# Step 4. Validate Exported FUMA Files
# ============================================================

validation_results <- data.frame(
  Dataset = character(),
  SNPs = integer(),
  Status = character()
)

expected_cols <- c(
  "chromosome",
  "base_pair_location",
  "effect_allele",
  "other_allele",
  "beta",
  "p_value"
)

for (csv_file in csv_files) {
  
  # Reconstruct output filename
  csv_name <- tools::file_path_sans_ext(basename(csv_file))
  
  # Extract gene ID and phenotype
  dataset_name <- sub(
    "^([^_]+_[^_]+).*",
    "\\1",
    csv_name
  )
  
  output_name <- paste0(
    dataset_name,
    "_FUMA_snp2gene_p",
    p_label,
    "_",
    genome_build,
    ".txt"
  )
  
  output_file <- file.path(output_dir, output_name)
  
  # Check whether output exists
  if (!file.exists(output_file)) {
    status <- "MISSING"
    snp_count <- 0L
    
  } else {
    
    # Read exported TXT
    data_check <- read.delim(
      output_file,
      check.names = FALSE
    )
    
    snp_count <- nrow(data_check)
    
    # Validate data structure and values
    valid <- identical(names(data_check), expected_cols) &&
      snp_count > 0 &&
      all(!is.na(data_check)) &&
      all(is.finite(data_check$p_value)) &&
      all(data_check$p_value >= 0 &
            data_check$p_value <= p_threshold)
    
    status <- if (valid) "PASS" else "FAIL"
  }
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      Dataset = dataset_name,
      SNPs = snp_count,
      Status = status
    )
  )
}

print(validation_results)

# Stop execution if any validation fails
if (!all(validation_results$Status == "PASS")) {
  stop("One or more FUMA output files failed validation.")
}

cat("\nAll FUMA output files passed validation.\n")
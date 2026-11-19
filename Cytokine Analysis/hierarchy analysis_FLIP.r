# ---- Packages ----
if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

suppressPackageStartupMessages({
  library(factoextra)
})

# ---- Load data (CSV or XLSX) ----
file_path <- "/Users/allisonpickle/Desktop/Composite Timed LPS.xlsx"  # <-- CHANGE ME
xlsx_sheet <- "Inflam_Combo"  # e.g., 1 OR "Sheet1"

# ==============================================================================
# NEW OPTION: Preserve sample order from input file?
# Set to TRUE to keep samples in left-to-right order as they appear in your file
# Set to FALSE to let hierarchical clustering determine sample order
# ==============================================================================
preserve_sample_order <- TRUE  # <-- CHANGE ME (TRUE or FALSE)

if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
  data_raw <- readxl::read_excel(path = file_path, sheet = xlsx_sheet)
  data_raw <- as.data.frame(data_raw, check.names = FALSE)
} else {
  data_raw <- read.csv(file_path, check.names = FALSE)
}

# Quick checks
head(data_raw)
tail(data_raw)
str(data_raw)

# ---- Expect FLIPPED layout ----
# Column 1 = feature names (cytokines), remaining columns = samples (headers in row 1)
feature_col <- names(data_raw)[1]
if (ncol(data_raw) < 2) stop("Expected at least 2 columns: feature column + >=1 sample column.")

data_raw[[feature_col]] <- as.character(data_raw[[feature_col]])

# EDIT: make feature labels plot-safe (avoid α/β conversion errors)
data_raw[[feature_col]] <- gsub("α", "alpha", data_raw[[feature_col]], fixed = TRUE)
data_raw[[feature_col]] <- gsub("β", "beta",  data_raw[[feature_col]], fixed = TRUE)

if (anyNA(data_raw[[feature_col]]) || any(data_raw[[feature_col]] == "")) {
  stop("First column (feature names) contains missing/blank values.")
}

# Use first column as rownames (feature names), then drop it
feature_names <- make.unique(data_raw[[feature_col]])
rownames(data_raw) <- feature_names
mat <- data_raw[, setdiff(names(data_raw), feature_col), drop = FALSE]

# STORE ORIGINAL SAMPLE ORDER from input file
original_sample_order <- colnames(mat)

# Coerce all sample columns to numeric safely (Excel often imports as character)
mat_num <- as.data.frame(lapply(mat, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE)

# EDIT: carry over Column A labels as rownames for cytokines/features
rownames(mat_num) <- feature_names

# Drop features that are all NA after coercion
keep_feature <- rowSums(!is.na(mat_num)) > 0
mat_num <- mat_num[keep_feature, , drop = FALSE]

# ---- Build a sample-by-feature data frame for sample clustering ----
# df: rows = samples, cols = features
df <- as.data.frame(t(as.matrix(mat_num)), check.names = FALSE)
rownames(df) <- make.unique(as.character(rownames(df)))

# Keep only numeric columns (extra safety)
numeric_cols <- vapply(df, is.numeric, logical(1))
df <- df[, numeric_cols, drop = FALSE]

# Remove rows with any NA (same behavior as your original script)
df <- na.omit(df)

# ---- Output folder for exported plots ----
sheet_tag <- if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
  gsub("[^A-Za-z0-9]+", "_", as.character(xlsx_sheet))
} else {
  "csv"
}

out_dir <- file.path(dirname(file_path), paste0("hc_outputs_", sheet_tag))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

base_tag <- paste0(tools::file_path_sans_ext(basename(file_path)), "__", sheet_tag, "__FLIP")


# =============================================================================
# Heatmap vs a reference sample: cytokines on Y, samples on X
# =============================================================================

# NEW OPTION: Toggle between log2 fold-change and raw values
use_log2fc <- FALSE  # Set to FALSE to use raw values instead of log2 fold-change

if (use_log2fc) {
  ref_sample <- "Control (6hr)"  # <-- CHANGE ME
  if (!(ref_sample %in% colnames(mat_num))) {
    stop(paste0(
      "Reference sample '", ref_sample, "' not found among sample columns. ",
      "Available columns include: ", paste(head(colnames(mat_num), 10), collapse = ", "), " ..."
    ))
  }

  pseudocount <- 1e-9
  ref_vals_feat <- as.numeric(mat_num[[ref_sample]])

  # Calculate log2 fold-change
  fc_heat <- log2((as.matrix(mat_num) + pseudocount) / (ref_vals_feat + pseudocount))
  fc_heat[, ref_sample] <- 0

  # Determine sample order based on user preference
  if (preserve_sample_order) {
    # Keep samples in their original left-to-right order from input file
    sample_order <- intersect(original_sample_order, colnames(fc_heat))
    fc_heat <- fc_heat[, sample_order, drop = FALSE]
    cluster_cols_setting <- FALSE
    order_tag <- "_ORIGINAL_ORDER"
    message("Using ORIGINAL sample order from input file (left to right)")
  } else {
    # Use hierarchical clustering to determine sample order
    cluster_cols_setting <- res.hc.sample
    order_tag <- "_CLUSTERED"
    message("Using HIERARCHICAL CLUSTERING to determine sample order")
  }

  heatmap_data <- fc_heat
  heatmap_title <- paste0("log2 fold-change vs ", ref_sample, " (samples on X axis)")
} else {
  # Use raw values directly
  heatmap_data <- as.matrix(mat_num)

  # Determine sample order based on user preference
  if (preserve_sample_order) {
    # Keep samples in their original left-to-right order from input file
    sample_order <- intersect(original_sample_order, colnames(heatmap_data))
    heatmap_data <- heatmap_data[, sample_order, drop = FALSE]
    cluster_cols_setting <- FALSE
    order_tag <- "_ORIGINAL_ORDER"
    message("Using ORIGINAL sample order from input file (left to right)")
  } else {
    # Use hierarchical clustering to determine sample order
    cluster_cols_setting <- res.hc.sample
    order_tag <- "_CLUSTERED"
    message("Using HIERARCHICAL CLUSTERING to determine sample order")
  }

  heatmap_title <- "Raw values (samples on X axis)"
}

heatmap_png <- file.path(out_dir, paste0(base_tag, "__heatmap_", if (use_log2fc) "log2FC" else "raw", "_samples_on_x", order_tag, ".png"))
heatmap_pdf <- file.path(out_dir, paste0(base_tag, "__heatmap_", if (use_log2fc) "log2FC" else "raw", "_samples_on_x", order_tag, ".pdf"))

lim <- max(abs(heatmap_data), na.rm = TRUE)

# PNG via ragg (Unicode-safe)
ragg::agg_png(heatmap_png, width = 2600, height = 1700, res = 200)
pheatmap::pheatmap(
  mat = heatmap_data,
  color = grDevices::colorRampPalette(c("blue", "white", "red"))(101),
  breaks = seq(-lim, lim, length.out = 102),
  cluster_rows = res.hc.feature,
  cluster_cols = cluster_cols_setting,
  gaps_row = NULL,
  gaps_col = NULL,
  fontsize_row = 7,
  fontsize_col = 9,
  angle_col = 45,
  main = heatmap_title
)
dev.off()

# PDF output
pdf(heatmap_pdf, width = 12, height = 8.5)
pheatmap::pheatmap(
  mat = heatmap_data,
  color = grDevices::colorRampPalette(c("blue", "white", "red"))(101),
  breaks = seq(-lim, lim, length.out = 102),
  cluster_rows = res.hc.feature,
  cluster_cols = cluster_cols_setting,
  gaps_row = NULL,
  gaps_col = NULL,
  fontsize_row = 7,
  fontsize_col = 9,
  angle_col = 45,
  main = heatmap_title
)
dev.off()

message("Reference heatmap exported to:\n  ", heatmap_png, "\n  ", heatmap_pdf)
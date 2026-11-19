# ============================================================
# Simplified Treatment Analysis with Significance Testing
# + Bar plots with significance stars
# + No random effects modeling
# + Simple t-tests or ANOVA
# ============================================================

# ----------------- PACKAGES -----------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

options(warn = -1)

# ============================================================
# SETTINGS
# ============================================================

input_dir <- '/Users/allisonpickle/Desktop/260603_3_26Diss_Excitotox/MM Stats'
file_pattern <- "\\.csv$"

rename_cols_by_position <- NULL

baseline_treatment <- "Tri_Ctrl "

# Match Treatments in CSV
treatment_levels <- c(
  "Tri_Ctrl",
  "Tri_w3",
  "Tri_Glu",
  "Tri_w3Glu"
)
# Match Treatments to strip plot Labels. Have enough \n- per number of treatment groups/levels. Align +/- with row labels below.
treatment_labels <- c(
  "Tri_ctrl" = "-\n-",
  "Tri_w3" = "+\n-",
  "Tri_Glu" = "-\n+",
  "Tri_w3Glu" = "+\n+"
)

# Override strip plot labels with these in order of treatment_levels
row_labels <- c(
  "20 µM w3",
  "5 µM Glutamate"
)

plot_title_base <- "Treatment"
x_axis_label <- ""
y_axis_label_raw <- "Results (mean ± SE)"
y_axis_label_log <- "log(Results + 0.5) (mean ± SE)"

# ============================================================
# HELPERS
# ============================================================

add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

make_aligned_tick_labels <- function(levels_vec, labels_named_vec) {
  out <- rep("", length(levels_vec))
  names(out) <- levels_vec
  
  common <- intersect(levels_vec, names(labels_named_vec))
  out[common] <- as.character(labels_named_vec[common])
  
  out
}

# ============================================================
# STATISTICAL TESTING
# ============================================================

compute_pairwise_pvalues <- function(data, response_col, baseline) {
  
  treatments <- levels(data$Treatment)
  
  p_df <- data.frame(
    Treatment = treatments,
    p.value = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Baseline gets NA (no comparison with itself)
  p_df$p.value[p_df$Treatment == baseline] <- NA_real_
  
  # Compare each treatment to baseline
  for (trt in setdiff(treatments, baseline)) {
    
    baseline_vals <- data[[response_col]][data$Treatment == baseline]
    treatment_vals <- data[[response_col]][data$Treatment == trt]
    
    # Remove NAs
    baseline_vals <- baseline_vals[!is.na(baseline_vals)]
    treatment_vals <- treatment_vals[!is.na(treatment_vals)]
    
    if (length(baseline_vals) >= 2 && length(treatment_vals) >= 2) {
      
      # Perform t-test
      test_result <- tryCatch(
        t.test(treatment_vals, baseline_vals, var.equal = FALSE),
        error = function(e) NULL
      )
      
      if (!is.null(test_result)) {
        p_df$p.value[p_df$Treatment == trt] <- test_result$p.value
      }
    }
  }
  
  p_df
}

# ============================================================
# SUMMARY TABLES
# ============================================================

build_summary_with_sig <- function(data, response_col, baseline) {
  
  summary_df <- data %>%
    dplyr::group_by(Treatment) %>%
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se = sd(.data[[response_col]], na.rm = TRUE) /
        sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )
  
  summary_df$Treatment <- factor(
    summary_df$Treatment,
    levels = levels(data$Treatment)
  )
  
  # Get p-values
  p_df <- compute_pairwise_pvalues(data, response_col, baseline)
  
  p_df$Treatment <- factor(
    p_df$Treatment,
    levels = levels(data$Treatment)
  )
  
  out <- dplyr::left_join(
    summary_df,
    p_df,
    by = "Treatment"
  )
  
  out$stars <- vapply(
    out$p.value,
    add_stars,
    character(1)
  )
  
  out
}

# ============================================================
# PLOTTING
# ============================================================

plot_bar_p_values <- function(summ_df,
                              raw_data,
                              response_col,
                              y_label,
                              title_text) {
  
  # Check if Biological_Rep column exists
  has_bio_rep <- "Biological_Rep" %in% names(raw_data)
  show_legend <- has_bio_rep && dplyr::n_distinct(raw_data$Biological_Rep, na.rm = TRUE) > 1
  
  aligned_treatment_labels <-
    make_aligned_tick_labels(
      treatment_levels,
      treatment_labels
    )
  
  aligned_treatment_labels <-
    stats::setNames(
      paste0("\n", aligned_treatment_labels),
      names(aligned_treatment_labels)
    )
  
  combined_row_label <- paste0(
    "\n",
    paste(row_labels, collapse = "\n")
  )
  
  p <- ggplot2::ggplot(
    summ_df,
    ggplot2::aes(
      x = Treatment,
      y = mean
    )
  ) +
    
    ggplot2::geom_col(
      fill = "#6c6c6c",
      alpha = 0.7
    ) +
    
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean - se,
        ymax = mean + se
      ),
      width = 0.2
    ) +
    
    ggplot2::geom_text(
      ggplot2::aes(
        y = mean + se,
        label = stars
      ),
      vjust = -0.5,
      size = 7,
      fontface = "bold"
    ) +
    
    ggplot2::scale_x_discrete(
      labels = aligned_treatment_labels,
      drop = FALSE
    ) +
    
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    
    ggplot2::labs(
      title = title_text,
      x = x_axis_label,
      y = y_label
    ) +
    
    ggplot2::theme_bw() +
    
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        vjust = -0.2,
        lineheight = 0.8,
        size = 9
      ),
      plot.margin = ggplot2::margin(
        5, 5, 45, 55
      ),
      legend.position = "none"
    ) +
    
    ggplot2::coord_cartesian(clip = "off")
  
  # Add individual points if Biological_Rep exists
  if (has_bio_rep && show_legend) {
    p <- p + ggplot2::geom_jitter(
      data = raw_data,
      ggplot2::aes(
        x = Treatment,
        y = .data[[response_col]],
        shape = Biological_Rep
      ),
      width = 0.15,
      alpha = 0.75,
      size = 2,
      color = "#1f5fbf"
    )
  } else {
    p <- p + ggplot2::geom_jitter(
      data = raw_data,
      ggplot2::aes(
        x = Treatment,
        y = .data[[response_col]]
      ),
      width = 0.15,
      alpha = 0.75,
      size = 2,
      color = "#1f5fbf"
    )
  }
  
  p <- p + ggplot2::annotate(
    "text",
    x = 0.4,
    y = -Inf,
    label = combined_row_label,
    hjust = 1,
    vjust = 1.3,
    size = 3,
    lineheight = 0.8
  )
  
  p
}

# ============================================================
# MAIN PER-FILE RUNNER
# ============================================================

run_one_csv <- function(input_file) {
  
  cat(
    "\nProcessing:",
    basename(input_file),
    "\n"
  )
  
  process_data <- read.csv(
    input_file,
    stringsAsFactors = FALSE
  )
  
  # --------------------------------------------------------
  # DEBUG: Print column names
  # --------------------------------------------------------
  cat("Column names found:", paste(names(process_data), collapse = ", "), "\n")
  
  # --------------------------------------------------------
  # Optional rename
  # --------------------------------------------------------
  if (!is.null(rename_cols_by_position) && length(rename_cols_by_position) > 0) {
    for (new_name in names(rename_cols_by_position)) {
      col_pos <- rename_cols_by_position[[new_name]]
      if (col_pos <= ncol(process_data)) {
        colnames(process_data)[col_pos] <- new_name
        cat("Renamed column", col_pos, "to", new_name, "\n")
      }
    }
  }
  
  # --------------------------------------------------------
  # Verify required columns exist
  # --------------------------------------------------------
  required_cols <- c("Treatment", "Results")
  missing_cols <- setdiff(required_cols, names(process_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(names(process_data), collapse = ", "))
  }
  
  # --------------------------------------------------------
  # Factor setup
  # --------------------------------------------------------
  process_data$Results <- suppressWarnings(
    as.numeric(process_data$Results)
  )
  
  process_data$Treatment <- factor(
    process_data$Treatment,
    levels = treatment_levels
  )
  
  # Optional: convert Biological_Rep to factor if it exists
  if ("Biological_Rep" %in% names(process_data)) {
    process_data$Biological_Rep <- factor(process_data$Biological_Rep)
  }
  
  # Remove rows with NA in required columns
  analysis_data <- process_data[
    !is.na(process_data$Results) & !is.na(process_data$Treatment),
  ]
  
  cat("Rows after removing NAs:", nrow(analysis_data), "\n")
  
  analysis_data$logResults <- log(analysis_data$Results + 0.5)
  
  # --------------------------------------------------------
  # Output folder
  # --------------------------------------------------------
  csv_base <- tools::file_path_sans_ext(
    basename(input_file)
  )
  
  out_dir <- file.path(
    dirname(input_file),
    paste0("Simple_Outputs_", csv_base)
  )
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  # --------------------------------------------------------
  # PLOTS - RAW
  # --------------------------------------------------------
  
  summ_raw <- build_summary_with_sig(
    analysis_data,
    "Results",
    baseline_treatment
  )
  
  # Save summary table
  write.csv(
    summ_raw,
    file.path(out_dir, "Summary_Raw.csv"),
    row.names = FALSE
  )
  
  p_raw <- plot_bar_p_values(
    summ_raw,
    analysis_data,
    "Results",
    y_axis_label_raw,
    paste0(plot_title_base, " (raw)")
  )
  
  ggplot2::ggsave(
    file.path(
      out_dir,
      "Barplot_Raw.png"
    ),
    p_raw,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  # --------------------------------------------------------
  # PLOTS - LOG
  # --------------------------------------------------------
  
  summ_log <- build_summary_with_sig(
    analysis_data,
    "logResults",
    baseline_treatment
  )
  
  # Save summary table
  write.csv(
    summ_log,
    file.path(out_dir, "Summary_Log.csv"),
    row.names = FALSE
  )
  
  p_log <- plot_bar_p_values(
    summ_log,
    analysis_data,
    "logResults",
    y_axis_label_log,
    paste0(plot_title_base, " (log)")
  )
  
  ggplot2::ggsave(
    file.path(
      out_dir,
      "Barplot_Log.png"
    ),
    p_log,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  cat(
    "Done:",
    out_dir,
    "\n"
  )
}

# ============================================================
# RUN ALL FILES
# ============================================================

csv_files <- list.files(
  input_dir,
  pattern = file_pattern,
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop(
    "No CSV files found in: ",
    input_dir
  )
}

cat(
  "Found ",
  length(csv_files),
  " file(s)\n",
  sep = ""
)

results <- lapply(
  csv_files,
  function(f) {
    
    tryCatch(
      
      run_one_csv(f),
      
      error = function(e) {
        
        message(
          "FAILED: ",
          basename(f),
          " -> ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
  }
)

cat("\n=== All files processed ===\n")
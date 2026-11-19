# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","emmeans","readxl")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Disable pager to prevent "Selection:" prompts
options(pager = "console")

# Working directory
setwd('/Users/allisonpickle/Desktop/Sadi_Data/Crack Width')

# Data
input_file  <- "Book6.xlsx"  # Can be .csv or .xlsx

# Read CSV or XLSX and combine data from columns A-D and F-I
read_and_combine_columns <- function(file_path) {
  # Determine file type
  file_ext <- tolower(tools::file_ext(file_path))
  
  cat("\n=== PROCESSING FILE ===\n")
  cat("File:", basename(file_path), "\n")
  cat("Type:", file_ext, "\n")
  
  # Read the file based on extension
  if (file_ext == "csv") {
    raw_data <- read.csv(file_path, stringsAsFactors = FALSE, 
                        fill = TRUE, check.names = FALSE)
  } else if (file_ext %in% c("xlsx", "xls")) {
    raw_data <- as.data.frame(readxl::read_excel(file_path))
  } else {
    stop("ERROR: Unsupported file type. Please use .csv, .xlsx, or .xls files.")
  }
  
  # Get column names
  col_names <- names(raw_data)
  
  cat("Detected", ncol(raw_data), "columns\n")
  cat("Column names:", paste(col_names, collapse = ", "), "\n")
  
  # Check if we have at least 9 columns (A-D and F-I = columns 1-4 and 6-9)
  if (ncol(raw_data) < 9) {
    # For CSV, try reading without header
    if (file_ext == "csv") {
      raw_data_no_header <- read.csv(file_path, header = FALSE, 
                                     stringsAsFactors = FALSE, 
                                     fill = TRUE, check.names = FALSE)
      cat("\nActual number of columns in file:", ncol(raw_data_no_header), "\n")
      
      if (ncol(raw_data_no_header) < 9) {
        stop("ERROR: CSV file needs at least 9 columns. Found ", 
             ncol(raw_data_no_header), " columns.")
      }
      
      # Use the no-header version and create column names
      raw_data <- raw_data_no_header
      names(raw_data) <- paste0("Col", 1:ncol(raw_data))
      col_names <- names(raw_data)
      
      cat("Using first row as data (no header detected)\n")
      cat("Assigned column names:", paste(col_names, collapse = ", "), "\n")
    } else {
      stop("ERROR: File needs at least 9 columns. Found ", ncol(raw_data), " columns.")
    }
  }
  
  # Extract first set of data (columns 1-4: A-D)
  part1 <- raw_data[, 1:4]
  names(part1) <- col_names[1:4]
  
  # Extract second set of data (columns 6-9: F-I, skipping E which is column 5)
  part2 <- raw_data[, 6:9]
  names(part2) <- col_names[1:4]  # Use same names as part1
  
  # Remove rows with all NA values from both parts
  part1 <- part1[rowSums(!is.na(part1)) > 0, ]
  part2 <- part2[rowSums(!is.na(part2)) > 0, ]
  
  # Combine the two parts vertically
  combined_data <- rbind(part1, part2)
  
  cat("Part 1 (columns 1-4):", nrow(part1), "rows\n")
  cat("Part 2 (columns 6-9):", nrow(part2), "rows\n")
  cat("Combined total:", nrow(combined_data), "rows\n\n")
  
  return(combined_data)
}

# Read and combine the data
process_data <- read_and_combine_columns(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "as-dealloyed"

# Treatment order and labels - just use the treatment names as they are
treatment_levels <- c(
  "as-dealloyed",
  "5760s",
  "9600s"
)

process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)
process_data$Image <- factor(process_data$Image)

# CHECK IF SAMPLE COLUMN EXISTS
if ("Sample" %in% names(process_data)) {
  process_data$Sample <- factor(process_data$Sample)
  message("Sample column found and converted to factor")
} else {
  stop("ERROR: No 'Sample' column found in CSV. Please check column names above and either:\n",
       "  1. Add a 'Sample' column to your CSV, or\n",
       "  2. Uncomment one of the options above to create/rename it")
}

# ---- Apply Log Transformation ----
cat("\n=== APPLYING LOG TRANSFORMATION ===\n")
process_data$logResults <- log(process_data$Results + 0.5)
cat("Log transformation: log(Results + 0.5)\n")

# ---- Fit Models ----

# Model with NESTED random effects (Image nested within Sample)
RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Sample/Image), 
                        data = process_data)

# Log-transformed model with same random effects structure
LogFit <- lmerTest::lmer(logResults ~ Treatment + (1|Sample/Image), 
                        data = process_data)

# Check for singular fit
if (isSingular(RawFit)) {
  cat("\n⚠️  WARNING: Singular fit detected in RawFit! Sample variance is zero.\n")
  cat("Refitting with simplified model: Results ~ Treatment + (1|Image)\n\n")
  RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Image), 
                          data = process_data)
}

if (isSingular(LogFit)) {
  cat("\n⚠️  WARNING: Singular fit detected in LogFit! Sample variance is zero.\n")
  cat("Refitting with simplified model: logResults ~ Treatment + (1|Image)\n\n")
  LogFit <- lmerTest::lmer(logResults ~ Treatment + (1|Image), 
                          data = process_data)
}

cat("\n=== RAW FIT MODEL ===\n")
anova(RawFit); summary(RawFit)

cat("\n=== LOG FIT MODEL ===\n")
anova(LogFit); summary(LogFit)

# ---- Sequential Pairwise Comparisons Function ----
perform_sequential_comparisons <- function(fit, model_name) {
  cat("\n=== SEQUENTIAL PAIRWISE COMPARISONS (", model_name, ") ===\n", sep = "")
  cat("Comparing consecutive treatment groups (excluding baseline)\n\n")
  
  # Get estimated marginal means
  emm <- emmeans(fit, ~ Treatment)
  
  # Dynamically generate sequential comparisons based on treatment levels
  n_treatments <- length(treatment_levels)
  
  # Create contrast vectors for sequential comparisons (excluding baseline)
  sequential_contrasts <- list()
  contrast_names <- c()
  
  for (i in 2:(n_treatments - 1)) {
    contrast_vec <- rep(0, n_treatments)
    contrast_vec[i] <- 1
    contrast_vec[i + 1] <- -1
    
    contrast_name <- paste(treatment_levels[i], "vs", treatment_levels[i + 1])
    contrast_names <- c(contrast_names, contrast_name)
    sequential_contrasts[[contrast_name]] <- contrast_vec
  }
  
  # Print what comparisons will be made
  cat("Generating the following sequential comparisons:\n")
  for (name in contrast_names) {
    cat("  -", name, "\n")
  }
  cat("\n")
  
  # Perform contrasts
  sequential_results <- contrast(emm, method = sequential_contrasts, adjust = "none")
  cat("Sequential Treatment Comparisons:\n")
  print(summary(sequential_results))
  
  # Calculate effect sizes (differences in means)
  cat("\n\nEffect Sizes (Mean Differences):\n")
  seq_summary <- summary(sequential_results)
  for (i in 1:nrow(seq_summary)) {
    sig_level <- if (seq_summary$p.value[i] < 0.001) "***" else 
                 if (seq_summary$p.value[i] < 0.01) "**" else
                 if (seq_summary$p.value[i] < 0.05) "*" else "ns"
    
    cat(sprintf("%s: %.4f (SE = %.4f, t = %.2f, p = %.4f %s)\n",
                seq_summary$contrast[i],
                seq_summary$estimate[i],
                seq_summary$SE[i],
                seq_summary$t.ratio[i],
                seq_summary$p.value[i],
                sig_level))
  }
  
  return(list(emm = emm, results = sequential_results))
}

# Run sequential comparisons for both models
raw_seq <- perform_sequential_comparisons(RawFit, "RawFit")
log_seq <- perform_sequential_comparisons(LogFit, "LogFit")

# Helper plotting function: ONLY std residuals vs fitted + QQ
plot_diag <- function(fit, tag) {
  par(mfrow = c(1, 2))
  fitted_vals <- fitted(fit)
  std_resid   <- resid(fit, type = "pearson")
  plot(fitted_vals, std_resid,
       xlab = "Fitted Values",
       ylab = "Standardized (Pearson) Residuals",
       main = paste(tag, "Std Residuals vs Fitted"))
  abline(h = 0, col = "red", lwd = 2)
  
  # Custom QQ plot without labels
  qq_data <- qqnorm(resid(fit), plot.it = FALSE)
  plot(qq_data$x, qq_data$y, 
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles",
       main = paste(tag, "QQ Plot"),
       pch = 16, col = "black")
  qqline(resid(fit), col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
}

# ---- Save ANOVA and summary outputs ----
save_model_summaries <- function() {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)

  out_file <- file.path(out_dir, "Mixed_Model_Summary.txt")

  sink(out_file)
  cat("================================================================\n")
  cat("MIXED EFFECTS MODEL SUMMARIES\n")
  cat("================================================================\n\n")
  
  # Extract formulas
  raw_formula <- formula(RawFit)
  log_formula <- formula(LogFit)
  
  # RawFit
  cat("=== RawFit Model (Untransformed) ===\n")
  cat("Formula: "); print(raw_formula)
  cat("\nANOVA:\n"); print(anova(RawFit))
  cat("\nSummary:\n"); print(summary(RawFit))
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  print(as.data.frame(VarCorr(RawFit)))
  
  # Add sequential pairwise comparisons for RawFit
  cat("\n\n================================================================\n")
  cat("SEQUENTIAL PAIRWISE COMPARISONS (RawFit)\n")
  cat("================================================================\n\n")
  
  emm_raw <- emmeans(RawFit, ~ Treatment)
  n_treatments <- length(treatment_levels)
  sequential_contrasts <- list()
  contrast_names <- c()
  
  for (i in 2:(n_treatments - 1)) {
    contrast_vec <- rep(0, n_treatments)
    contrast_vec[i] <- 1
    contrast_vec[i + 1] <- -1
    
    contrast_name <- paste(treatment_levels[i], "vs", treatment_levels[i + 1])
    contrast_names <- c(contrast_names, contrast_name)
    sequential_contrasts[[contrast_name]] <- contrast_vec
  }
  
  sequential_results_raw <- contrast(emm_raw, method = sequential_contrasts, adjust = "none")
  
  cat("Comparing consecutive treatment groups (excluding baseline):\n")
  cat("Treatment levels:", paste(treatment_levels, collapse = " → "), "\n\n")
  print(summary(sequential_results_raw))
  
  cat("\n\nInterpretation:\n")
  cat("----------------------------------------------------------------\n")
  seq_summary_raw <- summary(sequential_results_raw)
  for (i in 1:nrow(seq_summary_raw)) {
    sig_level <- if (seq_summary_raw$p.value[i] < 0.001) "***" else 
                 if (seq_summary_raw$p.value[i] < 0.01) "**" else
                 if (seq_summary_raw$p.value[i] < 0.05) "*" else "ns"
    
    cat(sprintf("%s:\n", seq_summary_raw$contrast[i]))
    cat(sprintf("  Mean difference: %.2f (SE = %.2f)\n", 
                seq_summary_raw$estimate[i], seq_summary_raw$SE[i]))
    cat(sprintf("  t(%.1f) = %.2f, p = %.4f %s\n\n",
                seq_summary_raw$df[i], seq_summary_raw$t.ratio[i], 
                seq_summary_raw$p.value[i], sig_level))
  }
  
  # LogFit
  cat("\n\n================================================================\n")
  cat("=== LogFit Model (Log Transformed) ===\n")
  cat("Formula: "); print(log_formula)
  cat("Transformation: log(Results + 0.5)\n")
  cat("\nANOVA:\n"); print(anova(LogFit))
  cat("\nSummary:\n"); print(summary(LogFit))
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  print(as.data.frame(VarCorr(LogFit)))
  
  # Add sequential pairwise comparisons for LogFit
  cat("\n\n================================================================\n")
  cat("SEQUENTIAL PAIRWISE COMPARISONS (LogFit)\n")
  cat("================================================================\n\n")
  
  emm_log <- emmeans(LogFit, ~ Treatment)
  sequential_results_log <- contrast(emm_log, method = sequential_contrasts, adjust = "none")
  
  cat("Comparing consecutive treatment groups (excluding baseline):\n")
  cat("Treatment levels:", paste(treatment_levels, collapse = " → "), "\n\n")
  print(summary(sequential_results_log))
  
  cat("\n\nInterpretation:\n")
  cat("----------------------------------------------------------------\n")
  seq_summary_log <- summary(sequential_results_log)
  for (i in 1:nrow(seq_summary_log)) {
    sig_level <- if (seq_summary_log$p.value[i] < 0.001) "***" else 
                 if (seq_summary_log$p.value[i] < 0.01) "**" else
                 if (seq_summary_log$p.value[i] < 0.05) "*" else "ns"
    
    cat(sprintf("%s:\n", seq_summary_log$contrast[i]))
    cat(sprintf("  Mean difference: %.4f (SE = %.4f)\n", 
                seq_summary_log$estimate[i], seq_summary_log$SE[i]))
    cat(sprintf("  t(%.1f) = %.2f, p = %.4f %s\n\n",
                seq_summary_log$df[i], seq_summary_log$t.ratio[i], 
                seq_summary_log$p.value[i], sig_level))
  }
  
  sink()

  message("Saved ANOVA and summary outputs to: ", out_file)
}

# ---- Save diagnostics ----
save_diagnostics <- function(save_data = TRUE) {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  models <- list(
    RawFit = RawFit,
    LogFit = LogFit
  )
  
  all_resid <- list()
  
  for (nm in names(models)) {
    fit <- models[[nm]]
    fitted_vals <- fitted(fit)
    std_resid   <- resid(fit, type = "pearson")
    
    ## 1) STANDARDIZED (Pearson) residuals vs fitted
    png(file.path(out_dir, paste0(nm, "_residuals.png")),
        width = 1200, height = 1000, res = 150)
    plot(fitted_vals, std_resid,
         xlab = "Fitted Values",
         ylab = "Standardized (Pearson) Residuals",
         main = paste(nm, "Std Residuals vs Fitted"))
    abline(h = 0, col = "red", lwd = 2)
    dev.off()
    
    ## 2) QQ plot - WITHOUT LABELS
    png(file.path(out_dir, paste0(nm, "_qq.png")),
        width = 1200, height = 1000, res = 150)
    qq_data <- qqnorm(resid(fit), plot.it = FALSE)
    plot(qq_data$x, qq_data$y, 
         xlab = "Theoretical Quantiles", 
         ylab = "Sample Quantiles",
         main = paste(nm, "QQ Plot"),
         pch = 16, col = "#0068f9")
    qqline(resid(fit), col = "red", lwd = 2)
    dev.off()
    
    if (save_data) {
      all_resid[[nm]] <- data.frame(
        Model       = nm,
        Fitted      = fitted_vals,
        StdResidual = std_resid,
        RawResidual = resid(fit)
      )
    }
  }
  
  if (save_data) {
    write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
  
  message("Saved diagnostic plots for RawFit and LogFit to: ", out_dir)
}

# ---------- Helpers: extract p-values and make barplots with stars ----------

extract_treatment_pvals <- function(fit, baseline_name = "as-dealloyed") {
  fe <- coef(summary(fit))

  rn <- rownames(fe)
  is_trt <- grepl("^Treatment", rn)
  if (!any(is_trt)) {
    stop("No Treatment* rows found in fixed effects.")
  }

  trt_rows  <- fe[is_trt, , drop = FALSE]
  trt_names <- sub("^Treatment", "", rownames(trt_rows))

  p_col <- grep("Pr\\(", colnames(trt_rows), value = TRUE)[1]
  if (is.na(p_col)) stop("Could not find p-value column in fixed effects.")

  pvals <- trt_rows[, p_col]
  names(pvals) <- trt_names

  data.frame(
    Treatment = c(baseline_name, trt_names),
    p.value   = c(NA, as.numeric(pvals)),
    stringsAsFactors = FALSE
  )
}

add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

build_bar_summary_with_sig <- function(data, fit, response_col, baseline_name = baseline_treatment) {
  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) /
             sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  # Re-apply original factor order to keep plotting order
  summary_df$Treatment <- factor(
    summary_df$Treatment,
    levels = levels(data$Treatment)
  )

  p_df <- extract_treatment_pvals(fit, baseline_name)
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars <- vapply(out$p.value, add_stars, character(1))
  out
}

# Simplified plotting function - all samples as circles with different colors
plot_bar_p_values <- function(summ_df, raw_data, response_col, y_label, title_text) {
  
  ggplot(summ_df, aes(x = Treatment, y = mean)) +
    geom_col(fill = "#6c6c6c", alpha = 0.7) +
    geom_jitter(data = raw_data, 
                aes(x = Treatment, y = .data[[response_col]], color = Sample), 
                width = 0.2, 
                alpha = 0.6,
                size = 1.0,
                shape = 16) +  # All circles (shape 16)
    scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00", "#FF1493", "#8B008B", "#00CED1", "#FF4500", "#FFD700")) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(
      aes(y = mean + se, label = stars),
      vjust = -0.5,
      size = 7,
      fontface = "bold"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      x = "Treatment",
      y = y_label, 
      title = title_text,
      color = "Sample"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

# ---------- Build and save barplots ----------

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base)
dir.create(out_dir, showWarnings = FALSE)

# Raw scale barplot
raw_bar_summary <- build_bar_summary_with_sig(
  data         = process_data,
  fit          = RawFit,
  response_col = "Results",
  baseline_name = baseline_treatment
)

raw_bar_plot <- plot_bar_p_values(
  raw_bar_summary,
  raw_data     = process_data,
  response_col = "Results",
  y_label      = "Results (mean ± SE)",
  title_text   = "Treatment Effects (RawFit)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_RawFit.png"),
  plot     = raw_bar_plot,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# Log scale barplot
log_bar_summary <- build_bar_summary_with_sig(
  data         = process_data,
  fit          = LogFit,
  response_col = "logResults",
  baseline_name = baseline_treatment
)

log_bar_plot <- plot_bar_p_values(
  log_bar_summary,
  raw_data     = process_data,
  response_col = "logResults",
  y_label      = "log(Results) (mean ± SE)",
  title_text   = "Treatment Effects (LogFit)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_LogFit.png"),
  plot     = log_bar_plot,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# ---- Run saves at the end ----
save_diagnostics()
save_model_summaries()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Comparing RawFit vs LogFit transformations\n")
cat("All outputs saved to:", out_dir, "\n")
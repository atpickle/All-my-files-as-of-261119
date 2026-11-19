# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","emmeans")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Working directory
setwd("/Users/allisonpickle/Desktop/Sadi_Data")

# Data
input_file  <- "Lig data_raw_reform_Time_Temp_Fixed Effects.csv"
process_data <- read.csv(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "as-dealloyed"

# Treatment order and labels - just use the treatment names as they are
treatment_levels <- c(
  "as-dealloyed",
  "1920s",
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

# CHECK IF TEMPERATURE COLUMN EXISTS
if ("Temperature" %in% names(process_data)) {
  process_data$Temperature <- factor(process_data$Temperature)
  message("Temperature column found and converted to factor")
  
  # Print temperature levels for verification
  cat("\nTemperature levels found:\n")
  print(levels(process_data$Temperature))
  cat("\n")
} else {
  warning("WARNING: No 'Temperature' column found in CSV. Model will only include Treatment.\n")
  # Set a flag to determine which model to use
  has_temperature <- FALSE
}

# Set flag if temperature exists
has_temperature <- "Temperature" %in% names(process_data)

# Model with NESTED random effects (Image nested within Sample)
# Include Temperature if available
if (has_temperature) {
  cat("Fitting model with Treatment and Temperature as fixed effects\n\n")
  RawFit <- lmerTest::lmer(Results ~ Treatment + Temperature + (1|Sample/Image), 
                          data = process_data)
} else {
  cat("Fitting model with Treatment only as fixed effect\n\n")
  RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Sample/Image), 
                          data = process_data)
}
                        
anova(RawFit); summary(RawFit)

# ---- Sequential Pairwise Comparisons ----
cat("\n=== SEQUENTIAL PAIRWISE COMPARISONS ===\n")
cat("Comparing consecutive treatment groups (excluding baseline)\n\n")

# Get estimated marginal means
emm <- emmeans(RawFit, ~ Treatment)

# Define sequential comparisons (1920s vs 5760s, then 5760s vs 9600s)
sequential_contrasts <- list(
  "1920s vs 5760s" = c(0, 1, -1, 0),
  "5760s vs 9600s" = c(0, 0, 1, -1)
)

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
  
  cat(sprintf("%s: %.2f nm (SE = %.2f, t = %.2f, p = %.4f %s)\n",
              seq_summary$contrast[i],
              seq_summary$estimate[i],
              seq_summary$SE[i],
              seq_summary$t.ratio[i],
              seq_summary$p.value[i],
              sig_level))
}

# If Temperature exists, also compare temperatures
if (has_temperature && length(levels(process_data$Temperature)) > 1) {
  cat("\n\n=== TEMPERATURE COMPARISONS ===\n")
  cat("Comparing different temperature levels\n\n")
  
  emm_temp <- emmeans(RawFit, ~ Temperature)
  temp_contrasts <- contrast(emm_temp, method = "pairwise", adjust = "none")
  
  cat("Pairwise Temperature Comparisons:\n")
  print(summary(temp_contrasts))
}

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
  cat("MIXED EFFECTS MODEL SUMMARY\n")
  cat("================================================================\n\n")
  
  # Extract formula
  raw_formula <- formula(RawFit)
  
  # RawFit
  cat("=== RawFit Model (Untransformed) ===\n")
  cat("Formula: "); print(raw_formula)
  cat("\nANOVA:\n"); print(anova(RawFit))
  cat("\nSummary:\n"); print(summary(RawFit))
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  print(as.data.frame(VarCorr(RawFit)))
  
  # Add sequential pairwise comparisons
  cat("\n\n================================================================\n")
  cat("SEQUENTIAL PAIRWISE COMPARISONS\n")
  cat("================================================================\n\n")
  
  emm <- emmeans(RawFit, ~ Treatment)
  
  sequential_contrasts <- list(
    "1920s vs 5760s" = c(0, 1, -1, 0),
    "5760s vs 9600s" = c(0, 0, 1, -1)
  )
  
  sequential_results <- contrast(emm, method = sequential_contrasts, adjust = "none")
  
  cat("Comparing consecutive treatment groups (excluding baseline):\n\n")
  print(summary(sequential_results))
  
  cat("\n\nInterpretation:\n")
  cat("----------------------------------------------------------------\n")
  seq_summary <- summary(sequential_results)
  for (i in 1:nrow(seq_summary)) {
    sig_level <- if (seq_summary$p.value[i] < 0.001) "***" else 
                 if (seq_summary$p.value[i] < 0.01) "**" else
                 if (seq_summary$p.value[i] < 0.05) "*" else "ns"
    
    cat(sprintf("%s:\n", seq_summary$contrast[i]))
    cat(sprintf("  Mean difference: %.2f nm (SE = %.2f)\n", 
                seq_summary$estimate[i], seq_summary$SE[i]))
    cat(sprintf("  t(%.1f) = %.2f, p = %.4f %s\n\n",
                seq_summary$df[i], seq_summary$t.ratio[i], 
                seq_summary$p.value[i], sig_level))
  }
  
  # Add temperature comparisons if available
  if (has_temperature && length(levels(process_data$Temperature)) > 1) {
    cat("\n\n================================================================\n")
    cat("TEMPERATURE COMPARISONS\n")
    cat("================================================================\n\n")
    
    emm_temp <- emmeans(RawFit, ~ Temperature)
    temp_contrasts <- contrast(emm_temp, method = "pairwise", adjust = "none")
    
    cat("Pairwise Temperature Comparisons:\n\n")
    print(summary(temp_contrasts))
    
    cat("\n\nInterpretation:\n")
    cat("----------------------------------------------------------------\n")
    temp_summary <- summary(temp_contrasts)
    for (i in 1:nrow(temp_summary)) {
      sig_level <- if (temp_summary$p.value[i] < 0.001) "***" else 
                   if (temp_summary$p.value[i] < 0.01) "**" else
                   if (temp_summary$p.value[i] < 0.05) "*" else "ns"
      
      cat(sprintf("%s:\n", temp_summary$contrast[i]))
      cat(sprintf("  Mean difference: %.2f nm (SE = %.2f)\n", 
                  temp_summary$estimate[i], temp_summary$SE[i]))
      cat(sprintf("  t(%.1f) = %.2f, p = %.4f %s\n\n",
                  temp_summary$df[i], temp_summary$t.ratio[i], 
                  temp_summary$p.value[i], sig_level))
    }
  }
  
  sink()

  message("Saved ANOVA and summary outputs to: ", out_file)
}

# ---- Save diagnostics ----
save_diagnostics <- function(save_data = TRUE) {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  fit <- RawFit
  fitted_vals <- fitted(fit)
  std_resid   <- resid(fit, type = "pearson")
  
  ## 1) STANDARDIZED (Pearson) residuals vs fitted
  png(file.path(out_dir, "RawFit_residuals.png"),
      width = 1200, height = 1000, res = 150)
  plot(fitted_vals, std_resid,
       xlab = "Fitted Values",
       ylab = "Standardized (Pearson) Residuals",
       main = "RawFit Std Residuals vs Fitted")
  abline(h = 0, col = "red", lwd = 2)
  dev.off()
  
  ## 2) QQ plot - WITHOUT LABELS
  png(file.path(out_dir, "RawFit_qq.png"),
      width = 1200, height = 1000, res = 150)
  qq_data <- qqnorm(resid(fit), plot.it = FALSE)
  plot(qq_data$x, qq_data$y, 
       xlab = "Theoretical Quantiles", 
       ylab = "Sample Quantiles",
       main = "RawFit QQ Plot",
       pch = 16, col = "#0068f9")
  qqline(resid(fit), col = "red", lwd = 2)
  dev.off()
  
  if (save_data) {
    all_resid <- data.frame(
      Model       = "RawFit",
      Fitted      = fitted_vals,
      StdResidual = std_resid,
      RawResidual = resid(fit)
    )
    
    write.csv(
      all_resid,
      file.path(out_dir, "Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
  
  message("Saved standardized (Pearson) residual vs fitted and QQ plots to: ", out_dir)
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
      x = "Treatment Duration",  # Updated label
      y = y_label, 
      title = title_text,
      color = "Sample",
      caption = if(has_temperature) paste0("Adjusted for Temperature (", 
                                          paste(levels(process_data$Temperature), collapse = ", "), 
                                          "°C)") else NULL  # Add caption showing temperature adjustment
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, color = "black"),
      plot.caption = element_text(hjust = 0.5, size = 8, color = "gray30"),  # Style caption
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}

# ---------- Build and save barplot ----------

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

# Updated title to reflect model structure
plot_title <- if(has_temperature) {
  "Treatment Effects (Controlling for Temperature)"
} else {
  "Treatment Effects"
}

raw_bar_plot <- plot_bar_p_values(
  raw_bar_summary,
  raw_data     = process_data,
  response_col = "Results",
  y_label      = "Ligament Thickness (nm, mean ± SE)",  # More descriptive
  title_text   = plot_title  # Updated title
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_RawFit.png"),
  plot     = raw_bar_plot,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# ---- Run saves at the end ----
save_diagnostics()
save_model_summaries()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All outputs saved to:", out_dir, "\n")
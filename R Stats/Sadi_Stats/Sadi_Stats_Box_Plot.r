# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","pbkrtest")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Disable pager to prevent "Selection:" prompts
options(pager = "console")

# Working directory
setwd("/Users/allisonpickle/Desktop/Sadi_Data")

# Data
input_file  <- "Lig data_raw_reform_375C.csv"
process_data <- read.csv(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "as-dealloyed"

# Treatment order and labels - just use the treatment names as they are
treatment_levels <- c(
  "as-dealloyed",
  "15s",
  "45s",
  "75s"
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

# Model with NESTED random effects (Image nested within Sample)
RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Sample/Image), 
                        data = process_data)          
anova(RawFit); summary(RawFit)

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
  vc <- as.data.frame(VarCorr(RawFit))
  print(vc)
  
  # Check for singular fit
  if (isSingular(RawFit)) {
    cat("\n⚠️  WARNING: Singular fit detected!\n")
    cat("One or more random effect variances are zero or near-zero.\n")
    cat("Consider simplifying random effects structure.\n")
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

# ---------- Helpers: extract p-values and make boxplots with stars ----------

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

# Box plot function
plot_bar_p_values <- function(summ_df, raw_data, response_col, y_label, title_text) {
  
  # Calculate max value for positioning stars
  max_val <- max(raw_data[[response_col]], na.rm = TRUE)
  
  ggplot(raw_data, aes(x = Treatment, y = .data[[response_col]])) +
    geom_boxplot(fill = "#6c6c6c", alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    geom_text(data = summ_df,
              aes(x = Treatment, y = max_val, label = stars),
              vjust = -0.5,
              size = 7,
              fontface = "bold",
              inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      x = "Treatment Duration",
      y = y_label, 
      title = title_text
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, color = "black")
    )
}

# ---------- Build and save boxplot ----------

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base)
dir.create(out_dir, showWarnings = FALSE)

# Raw scale boxplot
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
  y_label      = "Ligament Thickness (nm)",
  title_text   = "Treatment Effects"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Boxplot_RawFit.png"),
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
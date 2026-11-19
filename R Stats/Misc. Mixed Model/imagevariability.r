# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Working directory
setwd("/Users/allisonpickle/Desktop")

# Data
input_file  <- "PD_Microglia Area_RMixedFormated _BioRep.csv"
process_data <- read.csv(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "Ctrl"

# Create display labels for plotting with "-" and "+"
# The names of this vector define the specific order of treatment levels
treatment_labels <- c(
  "Ctrl"      = "-\n-\n-\n-", 
  "PD"        = "+\n-\n-\n-",
  "1ugLPS"    = "-\n+\n-\n-",
  "5ugLPS"    = "-\n-\n+\n-",
  "10ugLPS"   = "-\n-\n-\n+",
  "1ugLPSPD"  = "+\n+\n-\n-",
  "5ugLPSPD"  = "+\n-\n+\n-",
  "10ugLPSPD" = "+\n-\n-\n+"
)

# Extract levels directly from the names of the labels vector
treatment_levels <- names(treatment_labels)

process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)

# Define row labels for left side of strip plots
# Order corresponds to the rows in the strip plot (Top to Bottom)
row_labels <- c(
  "PD", 
  "1 µg/mL LPS", 
  "5 µg/mL LPS", 
  "10 µg/mL LPS"
)

process_data$Well <- factor(process_data$Well)

# CHECK IF BIOLOGICAL_REP COLUMN EXISTS
if ("Biological_Rep" %in% names(process_data)) {
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)
  message("Biological_Rep column found and converted to factor")
} else {
  stop("ERROR: No 'Biological_Rep' column found in CSV. Please check column names above and either:\n",
       "  1. Add a 'Biological_Rep' column to your CSV, or\n",
       "  2. Uncomment one of the options above to create/rename it")
}

# CHECK IF IMAGE COLUMN EXISTS
if ("Image" %in% names(process_data)) {
  process_data$Image <- factor(process_data$Image)
  message("Image column found and converted to factor")
} else {
  warning("WARNING: No 'Image' column found in CSV. Model will proceed without Image random effect.")
}

# Model 1 (raw) - Multiple random effects as separate terms
# This format allows wells and images to vary by treatment and independently accounts for biological replicate variation
RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Treatment:Well) + (1|Image) + (1|Biological_Rep), 
                        data = process_data)

anova(RawFit); summary(RawFit)

# Option 4: Biological_Rep with Image and Well as crossed effects
# RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Biological_Rep) + (1|Image) + (1|Well), 
#                         data = process_data)

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

# Transform: log
process_data$logResults <- log(process_data$Results + 0.5)

# Extract the formula from RawFit and replace the response variable with logResults
raw_formula <- formula(RawFit)
log_formula <- update(raw_formula, logResults ~ .)

# Fit LogFit using the updated formula
LogFit <- lmerTest::lmer(log_formula, data = process_data)

anova(LogFit); summary(LogFit)

# ---- Save ANOVA and summary outputs for RawFit and LogFit ----
save_model_summaries <- function() {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)

  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")

  sink(out_file)
  cat("================================================================\n")
  cat("MIXED EFFECTS MODEL SUMMARIES\n")
  cat("================================================================\n\n")
  
  # Extract and display actual model formulas
  raw_formula <- formula(RawFit)
  log_formula <- formula(LogFit)
  
  cat("=== RawFit Model ===\n")
  cat("Formula: "); print(raw_formula)
  cat("\nThis model accounts for:\n")
  cat("  - Treatment effects (fixed)\n")
  cat("  - Random effects as specified in the formula above\n\n")
  
  cat("ANOVA:\n"); print(anova(RawFit))
  cat("\nSummary:\n"); print(summary(RawFit))
  
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  vc <- as.data.frame(VarCorr(RawFit))
  print(vc)

  cat("\n\n================================================================\n")
  cat("=== LogFit Model ===\n")
  cat("Formula: "); print(log_formula)
  cat("\nThis model accounts for:\n")
  cat("  - Treatment effects (fixed)\n")
  cat("  - Random effects as specified in the formula above\n\n")
  
  cat("ANOVA:\n"); print(anova(LogFit))
  cat("\nSummary:\n"); print(summary(LogFit))
  
  cat("\n\nRandom Effects Variance Components:\n")
  cat("----------------------------------------------------------------\n")
  vc_log <- as.data.frame(VarCorr(LogFit))
  print(vc_log)
  
  cat("\n\n================================================================\n")
  cat("Model Comparison Notes:\n")
  cat("----------------------------------------------------------------\n")
  cat("RawFit uses:  "); cat(deparse(raw_formula), "\n")
  cat("LogFit uses:  "); cat(deparse(log_formula), "\n")
  cat("\nBoth models should have the same random effects structure.\n")
  cat("LogFit applies log transformation: logResults = log(Results + 0.5)\n")
  
  sink()

  message("Saved ANOVA and summary outputs to: ", out_file)
}

# ---- Save diagnostics: ONLY std_resid vs fitted + QQ, plus residual data ----
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
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
  
  message("Saved standardized (Pearson) residual vs fitted and QQ plots for each model to: ", out_dir)
}

# ---------- Helpers: extract p-values and make barplots with stars ----------

extract_treatment_pvals <- function(fit, baseline_name = "Ctrl") {
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

# Updated plotting function with scatter overlay - different symbols per biological replicate
plot_bar_p_values <- function(summ_df, raw_data, response_col, y_label, title_text) {
  # Calculate the maximum height needed (mean + se + some padding for stars)
  max_y <- max(summ_df$mean + summ_df$se, na.rm = TRUE)
  
  # Prepare labels with a leading newline to ensure alignment and spacing
  # This trick forces both the axis text and the annotation to start at the same vertical offset
  # relative to the bottom of the plot panel (y = -Inf)
  aligned_treatment_labels <- setNames(paste0("\n", treatment_labels), names(treatment_labels))
  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  
  ggplot(summ_df, aes(x = Treatment, y = mean)) +
    geom_col(fill = "#717171", alpha = 0.9) +  # Made bars slightly transparent
    geom_jitter(data = raw_data, 
                aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep, color = Biological_Rep), 
                width = 0.2, 
                alpha = 0.6,  # Changed to 0.6 for slight transparency
                size = 1.0) +   # Removed fixed color argument
    scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 7)) +  # Different shapes for up to 8 bio reps
    scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00", "black", "black", "black", "black", "black")) + # Dark Blue, Green, Orange
    scale_x_discrete(labels = aligned_treatment_labels) +  # Apply custom labels with "-" and "+"
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
       geom_text(
      aes(y = mean + se, label = stars), # Position at top of error bar
      vjust = -0.5,                      # Adjust vertical location
      size = 4,                          # Adjust asterick size
      fontface = "bold"                  # Make asterisks bold
    ) +
     scale_y_continuous(expand = expansion(mult = c(0, 0.10))) + # Changed bottom expansion to 0
    labs(
      x = "",  # Remove x-axis title
      y = y_label, 
      title = title_text,
      shape = "Biological Rep",  # Legend title
      color = "Biological Rep"   # Legend title for color (ensures legends merge)
    ) +
    # Add left-side row labels as a SINGLE multi-line annotation
    # y = -Inf places it at the bottom edge of the plot panel
    # vjust =  moves height of labels to align with strip labels
    annotate("text", x = 0.4, y = -Inf, label = combined_row_label, 
             hjust = 1, vjust = 1.1, size = 3, color = "black", lineheight = 0.8) +
    theme_bw() +
    theme(
      # vjust = 1 aligns the top of the text block to the axis tick/panel edge
      # margin(t = 0) ensures it starts exactly where the annotation starts
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, color = "black", lineheight = 0.8, size = 9, margin = margin(t = 0)),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.title.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, color = "black"), # Center the title and make it black
      plot.margin = margin(5, 5, 15, 50),  # Increase left and bottom margins to make room for labels
      legend.position = "right",  # Position legend on the right
      legend.title = element_text(size = 8),  # Make legend title smaller
      legend.text = element_text(size = 7)    # Make legend text smaller too
    ) +
    coord_cartesian(clip = "off")  # Allow drawing outside plot area
}

# ---------- Build and save two barplots: RawFit and LogFit ----------

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base)
dir.create(out_dir, showWarnings = FALSE)

# 1) Raw scale barplot, stars from RawFit fixed effects
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
  title_text   = "Treatment (RawFit)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_RawFit.png"),
  plot     = raw_bar_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 2) Log scale barplot, stars from LogFit fixed effects
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
  y_label      = "log(Results + 0.5) (mean ± SE)",
  title_text   = "Treatment (LogFit)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_LogFit.png"),
  plot     = log_bar_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# ---- Run saves at the end ----
save_diagnostics()
save_model_summaries()
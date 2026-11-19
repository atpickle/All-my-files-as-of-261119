# ============================================================
# SETUP & DATA LOADING
# ============================================================

pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","emmeans","pbkrtest","ggsignif")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# UPDATE FILE NAME Update this path to match Python output filename
input_file  <- '/Users/allisonpickle/Desktop/Gaussian Microglia_Combined_Sigma_RMixedFormatted.csv'
setwd(dirname(input_file))

# Read data (already in long format from Python)
process_data <- read.csv(input_file, stringsAsFactors = FALSE)

cat("\n=== DATA LOADED ===\n")
cat(sprintf("Total rows: %d\n", nrow(process_data)))
cat(sprintf("Columns: %s\n", paste(names(process_data), collapse = ", ")))

# ============================================================
# DEBUG: CHECK DATA STRUCTURE
# ============================================================

cat("\n=== DEBUGGING DATA ===\n")
cat("First 10 rows:\n")
print(head(process_data, 10))

cat("\n\nColumn types:\n")
print(sapply(process_data, class))

cat("\n\nNA counts per column:\n")
print(colSums(is.na(process_data)))

cat("\n\nUnique values per column:\n")
for (col in names(process_data)) {
  cat(sprintf("%s: %d unique values\n", col, length(unique(process_data[[col]]))))
  if (col %in% c("Treatment", "Sigma", "Biological_Rep", "Well", "Image")) {
    cat(sprintf("  Values: %s\n", paste(head(unique(process_data[[col]]), 10), collapse = ", ")))
  }
}

# ============================================================
# DATA PREPARATION
# ============================================================


# UPDATE SET CONTROL GROUP BASED ON WHAT EXCEL SHEET "TREATMENT" COLUMN CONTAINS
# Updated to match actual data values
baseline_treatment <- "Ctrl"
treatment_labels <- c("Ctrl" = "-\n-", "5ugLPS" = "-\n+\n", "5ugLPS" = "-\n+\n")
row_labels <- c("Control", "5ug/mL LPS") # CAN LABEL ANYTHING

# Check if treatment names match
cat("\n=== TREATMENT MATCHING ===\n")
cat("Expected treatments:", paste(names(treatment_labels), collapse = ", "), "\n")
cat("Actual treatments:", paste(unique(process_data$Treatment), collapse = ", "), "\n")

# Factor conversion with checks
process_data$Treatment <- factor(process_data$Treatment, levels = names(treatment_labels))
process_data$Sigma <- factor(process_data$Sigma)
process_data$Well <- factor(process_data$Well)
process_data$Biological_Rep <- factor(process_data$Biological_Rep)
process_data$Image <- factor(process_data$Image)

# Check for NA Results after factorization
cat("\n=== CHECKING FOR ISSUES ===\n")
cat(sprintf("NA in Treatment: %d\n", sum(is.na(process_data$Treatment))))
cat(sprintf("NA in Results: %d\n", sum(is.na(process_data$Results))))
cat(sprintf("Valid rows (non-NA Treatment & Results): %d\n", 
            sum(!is.na(process_data$Treatment) & !is.na(process_data$Results))))

# Remove rows with NA in critical columns
process_data <- process_data[!is.na(process_data$Treatment) & 
                               !is.na(process_data$Results) & 
                               !is.na(process_data$Sigma) &
                               !is.na(process_data$Biological_Rep), ]

cat(sprintf("Rows after removing NAs: %d\n", nrow(process_data)))

if (nrow(process_data) == 0) {
  stop("ERROR: No valid data remaining after removing NAs. Check your CSV file format and column names.")
}

# Print data structure
cat("\n=== DATA STRUCTURE ===\n")
cat(sprintf("Treatment counts:\n"))
print(table(process_data$Treatment, useNA = "ifany"))
cat(sprintf("\nSigma counts:\n"))
print(table(process_data$Sigma, useNA = "ifany"))
cat(sprintf("\nBiological replicates: %s\n", 
            paste(sort(unique(process_data$Biological_Rep)), collapse = ", ")))

# Check biological replicates
n_bioreps <- length(unique(process_data$Biological_Rep))
cat(sprintf("\n=== BIOLOGICAL REPLICATE CHECK ===\n"))
cat(sprintf("Number: %d\n", n_bioreps))

if (n_bioreps < 3) {
  cat("\n⚠️  WARNING: Only", n_bioreps, "biological replicate(s) detected ⚠️\n")
  cat("Simplified model structure will be used. Results are not suitable for biological inference.\n\n")
} else {
  cat(sprintf("\n✓ Sufficient biological replicates (%d ≥ 3)\n\n", n_bioreps))
}

# Log transformation
process_data$logResults <- log(process_data$Results + 0.5)

# ============================================================
# BUILD MODELS - RAW DATA
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("BUILDING INTERACTION MODEL: Treatment × Sigma (RAW DATA)\n")
cat(rep("=", 80), "\n\n")

# Determine random effects structure
has_image <- "Image" %in% names(process_data)
if (n_bioreps >= 3 && has_image) {
  random_effects <- "(1 | Biological_Rep) + (1 | Biological_Rep:Well) + (1 | Well:Image)"
  model_type <- "Full nested with Treatment × Sigma interaction"
} else if (n_bioreps >= 3) {
  random_effects <- "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
  model_type <- "Nested with Treatment × Sigma interaction"
} else if (has_image) {
  random_effects <- "(1 | Well) + (1 | Well:Image)"
  model_type <- "⚠️ SIMPLIFIED with Treatment × Sigma interaction"
} else {
  random_effects <- "(1 | Well)"
  model_type <- "⚠️ SIMPLIFIED with Treatment × Sigma interaction"
}

# Raw data model (REML for estimation)
model_formula <- as.formula(paste("Results ~ Treatment * Sigma +", random_effects))
InteractionModel <- lmerTest::lmer(model_formula, data = process_data)

cat(sprintf("Model: %s\n", model_type))
cat("Formula: "); print(model_formula)
cat("\nANOVA Table:\n"); print(anova(InteractionModel))

# Null and full models for LRT (ML for comparison)
null_formula <- as.formula(paste("Results ~ 1 +", random_effects))
InteractionModel_null <- lme4::lmer(null_formula, data = process_data, REML = FALSE)
InteractionModel_full <- lme4::lmer(model_formula, data = process_data, REML = FALSE)

# ============================================================
# BUILD MODELS - LOG-TRANSFORMED DATA
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("BUILDING INTERACTION MODEL: Treatment × Sigma (LOG-TRANSFORMED)\n")
cat(rep("=", 80), "\n\n")

log_model_formula <- as.formula(paste("logResults ~ Treatment * Sigma +", random_effects))
LogInteractionModel <- lmerTest::lmer(log_model_formula, data = process_data)

cat("Formula: "); print(log_model_formula)
cat("\nANOVA Table:\n"); print(anova(LogInteractionModel))

# Null and full models for LRT
log_null_formula <- as.formula(paste("logResults ~ 1 +", random_effects))
LogInteractionModel_null <- lme4::lmer(log_null_formula, data = process_data, REML = FALSE)
LogInteractionModel_full <- lme4::lmer(log_model_formula, data = process_data, REML = FALSE)

# ============================================================
# ICC CALCULATION
# ============================================================

compute_icc <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  var_random <- sum(vc$vcov[vc$grp != "Residual"])
  var_resid  <- vc$vcov[vc$grp == "Residual"]
  return(var_random / (var_random + var_resid))
}

icc_raw <- compute_icc(InteractionModel)
icc_log <- compute_icc(LogInteractionModel)

cat(sprintf("\nICC (Raw Model): %.4f\n", icc_raw))
cat(sprintf("ICC (Log Model): %.4f\n", icc_log))

# ============================================================
# LIKELIHOOD RATIO TESTS
# ============================================================

lrt_raw <- anova(InteractionModel_null, InteractionModel_full)
lrt_log <- anova(LogInteractionModel_null, LogInteractionModel_full)

cat("\n=== LIKELIHOOD RATIO TEST (Raw Model) ===\n")
print(lrt_raw)

cat("\n=== LIKELIHOOD RATIO TEST (Log Model) ===\n")
print(lrt_log)

# ============================================================
# POST-HOC TESTS
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("POST-HOC COMPARISONS\n")
cat(rep("=", 80), "\n\n")

comparisons <- list(
  sigma = list(formula = ~ Treatment | Sigma, name = "Treatment | Sigma"),
  treatment = list(formula = ~ Treatment, name = "Overall Treatment Effect")
)

emm_results <- list()
for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  emm_results[[name]] <- list(
    emm = emmeans(LogInteractionModel, comp$formula, lmer.df = "satterthwaite"),
    pairs = pairs(emmeans(LogInteractionModel, comp$formula, lmer.df = "satterthwaite"))
  )
  cat(sprintf("\n%s:\n", comp$name))
  print(emm_results[[name]]$pairs)
}

emm_pairs_sigma_df <- as.data.frame(emm_results$sigma$pairs)
emm_pairs_sigma_df$Sigma <- sub("Sigma ", "", emm_pairs_sigma_df$Sigma)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

get_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

save_model_summaries <- function() {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  sink(file.path(out_dir, "Treatment_Sigma_Model_Summary.txt"))
  cat("================================================================\n")
  cat("TREATMENT × SIGMA INTERACTION MODEL SUMMARY\n")
  cat("================================================================\n\n")
  
  cat("=== DATA STRUCTURE ===\n")
  cat(sprintf("Biological replicates: %d (%s)\n", n_bioreps, 
              paste(sort(unique(process_data$Biological_Rep)), collapse = ", ")))
  cat(sprintf("Sigma values: %d (%s)\n\n", length(unique(process_data$Sigma)),
              paste(sort(unique(process_data$Sigma)), collapse = ", ")))
  
  if (n_bioreps < 3) {
    cat("⚠️  WARNING: INSUFFICIENT BIOLOGICAL REPLICATES\n")
    cat("Results not suitable for biological inference.\n\n")
  }
  
  cat("=== RAW DATA MODEL (InteractionModel) ===\n")
  cat("Formula: "); print(model_formula)
  cat(sprintf("Type: %s\n\n", model_type))
  
  cat("ANOVA:\n"); print(anova(InteractionModel))
  cat("\nSummary:\n"); print(summary(InteractionModel))
  cat("\nRandom Effects:\n"); print(as.data.frame(VarCorr(InteractionModel)))
  cat(sprintf("\nICC: %.4f\n", icc_raw))
  
  cat("\n\nLikelihood Ratio Test (Treatment × Sigma Effect):\n")
  cat("----------------------------------------------------------------\n")
  print(lrt_raw)
  
  cat("\n\n=== LOG-TRANSFORMED MODEL (LogInteractionModel) ===\n")
  cat("Formula: "); print(log_model_formula)
  cat("Transform: logResults = log(Results + 0.5)\n\n")
  
  cat("ANOVA:\n"); print(anova(LogInteractionModel))
  cat("\nSummary:\n"); print(summary(LogInteractionModel))
  cat("\nRandom Effects:\n"); print(as.data.frame(VarCorr(LogInteractionModel)))
  cat(sprintf("\nICC: %.4f\n", icc_log))
  
  cat("\n\nLikelihood Ratio Test (Treatment × Sigma Effect):\n")
  cat("----------------------------------------------------------------\n")
  print(lrt_log)
  
  cat("\n\n=== PAIRWISE COMPARISONS (Log Model) ===\n")
  for (name in names(comparisons)) {
    cat(sprintf("\n%s:\n", comparisons[[name]]$name))
    cat("----------------------------------------------------------------\n")
    print(emm_results[[name]]$pairs)
  }
  
  sink()
  message("Saved model summary to: ", file.path(out_dir, "Treatment_Sigma_Model_Summary.txt"))
}

save_diagnostics <- function(save_data = TRUE) {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  models <- list(
    InteractionModel = InteractionModel,
    LogInteractionModel = LogInteractionModel
  )
  
  all_resid <- list()
  
  for (nm in names(models)) {
    fit <- models[[nm]]
    fitted_vals <- fitted(fit)
    std_resid   <- resid(fit, type = "pearson")
    
    png(file.path(out_dir, paste0(nm, "_residuals.png")),
        width = 1200, height = 1000, res = 150)
    plot(fitted_vals, std_resid, xlab = "Fitted Values", 
         ylab = "Standardized (Pearson) Residuals",
         main = paste(nm, "Std Residuals vs Fitted"))
    abline(h = 0, col = "red", lwd = 2)
    dev.off()
    
    png(file.path(out_dir, paste0(nm, "_qq.png")),
        width = 1200, height = 1000, res = 150)
    qq_data <- qqnorm(resid(fit), plot.it = FALSE)
    plot(qq_data$x, qq_data$y, xlab = "Theoretical Quantiles", 
         ylab = "Sample Quantiles",
         main = paste(nm, "QQ Plot"), pch = 16, col = "#0068f9")
    qqline(resid(fit), col = "red", lwd = 2)
    dev.off()
    
    if (save_data) {
      all_resid[[nm]] <- data.frame(
        Model = nm,
        Fitted = fitted_vals,
        StdResidual = std_resid,
        RawResidual = resid(fit)
      )
    }
  }
  
  if (save_data) {
    write.csv(do.call(rbind, all_resid),
              file.path(out_dir, "All_Calculated_Residuals.csv"),
              row.names = FALSE)
  }
  
  message("Saved diagnostic plots to: ", out_dir)
}

# ============================================================
# VISUALIZATION
# ============================================================

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base)
dir.create(out_dir, showWarnings = FALSE)

# Interaction plot
p_sigma <- emmip(LogInteractionModel, Treatment ~ Sigma, CIs = TRUE) +
  scale_color_manual(values = c("#444444", "#FF6B6B"), 
                     labels = c("Control", "LPS 5µg/mL")) +
  labs(title = "Treatment Effect Across Gaussian Sigma Values",
       x = "Gaussian Sigma", y = "log(Results + 0.5)", color = "Treatment") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(out_dir, "Treatment_Sigma_Interaction.png"), 
       p_sigma, width = 8, height = 5, dpi = 300)

# Bar plot function
create_barplot <- function(subset_data, p_val, title_text) {
  summ_df <- subset_data |>
    group_by(Treatment) |>
    summarise(mean = mean(logResults, na.rm = TRUE),
              se = sd(logResults, na.rm = TRUE) / sqrt(sum(!is.na(logResults))),
              .groups = "drop")
  
  stars <- if(length(p_val) > 0) get_stars(p_val) else "ns"
  y_max <- max(summ_df$mean + summ_df$se)
  y_annot <- y_max + 0.15 * (max(subset_data$logResults) - min(subset_data$logResults))
  
  aligned_treatment_labels <- setNames(paste0("\n", treatment_labels), names(treatment_labels))
  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  
  p <- ggplot(summ_df, aes(x = Treatment, y = mean, fill = Treatment)) +
    geom_col(alpha = 0.7) +
    geom_jitter(data = subset_data, 
                aes(x = Treatment, y = logResults, 
                    shape = Biological_Rep, color = Biological_Rep),
                width = 0.2, alpha = 0.6, size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_fill_manual(values = c("#444444", "#FF6B6B")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 7)) +
    scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00", 
                                  "black", "black", "black", "black", "black")) +
    scale_x_discrete(labels = aligned_treatment_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(title = title_text, x = "", y = "log(Results + 0.5) (mean ± SE)",
         shape = "Biological Rep", color = "Biological Rep") +
    annotate("text", x = 0.4, y = -Inf, label = combined_row_label, 
             hjust = 1, vjust = 1.1, size = 3, color = "black", lineheight = 0.8) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.margin = margin(5, 5, 15, 50)
    ) +
    coord_cartesian(clip = "off")
  
  if (stars != "ns") {
    p <- p +
      geom_segment(aes(x = 1, xend = 2, y = y_annot, yend = y_annot), 
                   color = "black", linewidth = 0.5) +
      geom_segment(aes(x = 1, xend = 1, y = y_annot, yend = y_annot - 0.02), 
                   color = "black", linewidth = 0.5) +
      geom_segment(aes(x = 2, xend = 2, y = y_annot, yend = y_annot - 0.02), 
                   color = "black", linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_annot + 0.05, label = stars, 
               size = 7, fontface = "bold")
  }
  return(p)
}

# Bar plots for each sigma
for (sig in levels(process_data$Sigma)) {
  subset_data <- process_data[process_data$Sigma == sig, ]
  p_val <- emm_pairs_sigma_df$p.value[emm_pairs_sigma_df$Sigma == sig]
  p <- create_barplot(subset_data, p_val, paste("Gaussian Sigma", sig))
  ggsave(file.path(out_dir, paste0("Barplot_Sigma_", sig, ".png")),
         p, width = 6, height = 4, dpi = 300)
}

# ============================================================
# RUN SAVES & FINISH
# ============================================================

save_diagnostics()
save_model_summaries()

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")
cat(sprintf("Outputs saved to: %s/\n", out_dir))
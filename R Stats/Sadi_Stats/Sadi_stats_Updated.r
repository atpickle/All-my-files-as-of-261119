# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","emmeans")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Disable pager to prevent "Selection:" prompts
options(pager = "console")

# Suppress singularity and df calculation warnings for large datasets
options(lmerTest.limit = Inf)  # Disable lmerTest df calculation warnings
options(pbkrtest.limit = Inf)  # Disable pbkrtest warnings
emm_options(lmerTest.limit = Inf, pbkrtest.limit = Inf)  # For emmeans

# Working directory
setwd('/Users/allisonpickle/Desktop/Sadi_Data/Ligament Data Raw')

# Data
input_file  <- "Lig data_raw_reform_300C.csv"
process_data <- read.csv(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "as-dealloyed"

# Treatment order and labels - just use the treatment names as they are
treatment_levels <- c(
  "as-dealloyed",
  "300s",
  "360s",
  "600s"
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
                        data = process_data,
                        control = lmerControl(optimizer = "bobyqa"))

# Log-transformed model with same random effects structure
LogFit <- lmerTest::lmer(logResults ~ Treatment + (1|Sample/Image), 
                        data = process_data,
                        control = lmerControl(optimizer = "bobyqa"))

# Check for singular fit
if (isSingular(RawFit)) {
  cat("\n⚠️  WARNING: Singular fit detected in RawFit! Sample variance is zero.\n")
  cat("Refitting with simplified model: Results ~ Treatment + (1|Image)\n\n")
  RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Image), 
                          data = process_data,
                          control = lmerControl(optimizer = "bobyqa"))
}

if (isSingular(LogFit)) {
  cat("\n⚠️  WARNING: Singular fit detected in LogFit! Sample variance is zero.\n")
  cat("Refitting with simplified model: logResults ~ Treatment + (1|Image)\n\n")
  LogFit <- lmerTest::lmer(logResults ~ Treatment + (1|Image), 
                          data = process_data,
                          control = lmerControl(optimizer = "bobyqa"))
}

cat("\n=== RAW FIT MODEL ===\n")
anova(RawFit); summary(RawFit)

cat("\n=== LOG FIT MODEL ===\n")
anova(LogFit); summary(LogFit)

# ---- Sequential Pairwise Comparisons Function ----
perform_sequential_comparisons <- function(fit, model_name) {
  cat("\n=== SEQUENTIAL PAIRWISE COMPARISONS (", model_name, ") ===\n", sep = "")
  cat("Comparing consecutive treatment groups\n\n")
  
  # Get estimated marginal means
  emm <- emmeans(fit, ~ Treatment)
  
  # Use built-in consecutive comparison method
  sequential_results <- pairs(emm, reverse = TRUE)
  
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
  
  # Test significance of random effects using likelihood ratio test
  cat("\n\n================================================================\n")
  cat("TESTING SIGNIFICANCE OF RANDOM EFFECTS (Likelihood Ratio Tests)\n")
  cat("================================================================\n")
  
  # Get the random effects structure
  re_formula <- formula(RawFit)
  
  # Test if we have nested structure
  if (grepl("\\|.*\\/", deparse(re_formula))) {
    cat("\nNested random effects detected. Testing each level:\n\n")
    
    # Fit reduced models
    fit_no_sample <- tryCatch({
      lmerTest::lmer(Results ~ Treatment + (1|Image), 
                     data = process_data,
                     control = lmerControl(optimizer = "bobyqa"))
    }, error = function(e) NULL)
    
    fit_no_image <- tryCatch({
      lmerTest::lmer(Results ~ Treatment + (1|Sample), 
                     data = process_data,
                     control = lmerControl(optimizer = "bobyqa"))
    }, error = function(e) NULL)
    
    # Test Sample variance (compare full vs no sample)
    if (!is.null(fit_no_sample)) {
      lrt_sample <- anova(RawFit, fit_no_sample)
      cat("1. Test for Sample-level variance:\n")
      cat("   H0: Sample variance = 0\n")
      cat("   Chi-square = ", round(lrt_sample$Chisq[2], 4), ", df = ", lrt_sample$`Chi Df`[2], 
          ", p = ", format.pval(lrt_sample$`Pr(>Chisq)`[2], digits = 4), "\n")
      if (!is.na(lrt_sample$`Pr(>Chisq)`[2]) && lrt_sample$`Pr(>Chisq)`[2] < 0.05) {
        cat("   → Sample-to-sample variation IS statistically significant ***\n\n")
      } else {
        cat("   → Sample-to-sample variation is NOT statistically significant\n\n")
      }
    }
    
    # Test Image variance (compare full vs no image)
    if (!is.null(fit_no_image)) {
      lrt_image <- anova(RawFit, fit_no_image)
      cat("2. Test for Image-level variance (within Sample):\n")
      cat("   H0: Image variance = 0\n")
      cat("   Chi-square = ", round(lrt_image$Chisq[2], 4), ", df = ", lrt_image$`Chi Df`[2], 
          ", p = ", format.pval(lrt_image$`Pr(>Chisq)`[2], digits = 4), "\n")
      if (!is.na(lrt_image$`Pr(>Chisq)`[2]) && lrt_image$`Pr(>Chisq)`[2] < 0.05) {
        cat("   → Image-to-image variation IS statistically significant ***\n\n")
      } else {
        cat("   → Image-to-image variation is NOT statistically significant\n\n")
      }
    }
  } else {
    # Simple random effects structure (just Image) - FIX HERE
    cat("\nSimple random effects structure detected. Testing Image variance:\n\n")
    
    # Refit both models with ML (not REML) for likelihood ratio test
    RawFit_ML <- update(RawFit, REML = FALSE)
    fit_fixed <- lm(Results ~ Treatment, data = process_data)
    
    # Calculate likelihood ratio manually
    ll_mixed <- logLik(RawFit_ML)
    ll_fixed <- logLik(fit_fixed)
    chi_sq <- 2 * (as.numeric(ll_mixed) - as.numeric(ll_fixed))
    df_diff <- 1  # One random effect parameter
    p_val <- pchisq(chi_sq, df = df_diff, lower.tail = FALSE)
    
    cat("Test for Image-level variance:\n")
    cat("H0: Image variance = 0\n")
    cat("Chi-square = ", round(chi_sq, 4), ", df = ", df_diff, 
        ", p = ", format.pval(p_val, digits = 4), "\n")
    
    if (!is.na(p_val) && p_val < 0.05) {
      cat("→ Image-to-image variation IS statistically significant ***\n\n")
    } else {
      cat("→ Image-to-image variation is NOT statistically significant\n\n")
    }
  }
  
  # Calculate Intraclass Correlation Coefficients (ICC)
  cat("\n================================================================\n")
  cat("INTRACLASS CORRELATION COEFFICIENTS (ICC)\n")
  cat("================================================================\n")
  cat("ICC quantifies the proportion of variance explained by grouping structure.\n")
  cat("Values closer to 0 indicate little clustering; closer to 1 indicate strong clustering.\n\n")
  
  vc <- as.data.frame(VarCorr(RawFit))
  total_var <- sum(vc$vcov)
  
  for (i in 1:nrow(vc)) {
    if (vc$grp[i] != "Residual") {
      icc <- vc$vcov[i] / total_var
      cat(sprintf("%s ICC = %.6f (%.3f%% of total variance)\n", 
                  vc$grp[i], icc, icc * 100))
      if (icc < 0.01) {
        cat("   → Very weak clustering (negligible - random effect may not be necessary)\n\n")
      } else if (icc < 0.05) {
        cat("   → Weak clustering (minimal impact)\n\n")
      } else if (icc < 0.10) {
        cat("   → Moderate clustering\n\n")
      } else {
        cat("   → Strong clustering (substantial impact)\n\n")
      }
    }
  }
  cat(sprintf("Residual variance = %.4f (%.3f%% of total)\n\n", 
              vc$vcov[nrow(vc)], vc$vcov[nrow(vc)] / total_var * 100))
  
  # Add sequential pairwise comparisons for RawFit
  cat("\n================================================================\n")
  cat("SEQUENTIAL PAIRWISE COMPARISONS (RawFit)\n")
  cat("================================================================\n\n")
  
  emm_raw <- emmeans(RawFit, ~ Treatment)
  n_treatments <- length(treatment_levels)
  sequential_contrasts_raw <- list()
  contrast_names_raw <- c()
  
  for (i in 2:n_treatments) {
    contrast_vec <- rep(0, n_treatments)
    contrast_vec[i - 1] <- 1
    contrast_vec[i] <- -1
    
    contrast_name <- paste(treatment_levels[i - 1], "vs", treatment_levels[i])
    contrast_names_raw <- c(contrast_names_raw, contrast_name)
    sequential_contrasts_raw[[contrast_name]] <- contrast_vec
  }
  
  sequential_results_raw <- contrast(emm_raw, method = sequential_contrasts_raw, adjust = "none")
  
  cat("Comparing consecutive treatment groups:\n")
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
  
  # Test significance of random effects for LogFit
  cat("\n\n================================================================\n")
  cat("TESTING SIGNIFICANCE OF RANDOM EFFECTS (Likelihood Ratio Tests)\n")
  cat("================================================================\n")
  
  # Check LogFit structure
  log_re_formula <- formula(LogFit)
  
  if (grepl("\\|.*\\/", deparse(log_re_formula))) {
    cat("\nNested random effects detected. Testing each level:\n\n")
    
    # Fit reduced models
    fit_log_no_sample <- tryCatch({
      lmerTest::lmer(logResults ~ Treatment + (1|Image), 
                     data = process_data,
                     control = lmerControl(optimizer = "bobyqa"))
    }, error = function(e) NULL)
    
    fit_log_no_image <- tryCatch({
      lmerTest::lmer(logResults ~ Treatment + (1|Sample), 
                     data = process_data,
                     control = lmerControl(optimizer = "bobyqa"))
    }, error = function(e) NULL)
    
    # Test Sample variance
    if (!is.null(fit_log_no_sample)) {
      lrt_log_sample <- anova(LogFit, fit_log_no_sample)
      cat("1. Test for Sample-level variance:\n")
      cat("   H0: Sample variance = 0\n")
      cat("   Chi-square = ", round(lrt_log_sample$Chisq[2], 4), ", df = ", lrt_log_sample$`Chi Df`[2], 
          ", p = ", format.pval(lrt_log_sample$`Pr(>Chisq)`[2], digits = 4), "\n")
      if (!is.na(lrt_log_sample$`Pr(>Chisq)`[2]) && lrt_log_sample$`Pr(>Chisq)`[2] < 0.05) {
        cat("   → Sample-to-sample variation IS statistically significant ***\n\n")
      } else {
        cat("   → Sample-to-sample variation is NOT statistically significant\n\n")
      }
    }
    
    # Test Image variance
    if (!is.null(fit_log_no_image)) {
      lrt_log_image <- anova(LogFit, fit_log_no_image)
      cat("2. Test for Image-level variance (within Sample):\n")
      cat("   H0: Image variance = 0\n")
      cat("   Chi-square = ", round(lrt_log_image$Chisq[2], 4), ", df = ", lrt_log_image$`Chi Df`[2], 
          ", p = ", format.pval(lrt_log_image$`Pr(>Chisq)`[2], digits = 4), "\n")
      if (!is.na(lrt_log_image$`Pr(>Chisq)`[2]) && lrt_log_image$`Pr(>Chisq)`[2] < 0.05) {
        cat("   → Image-to-image variation IS statistically significant ***\n\n")
      } else {
        cat("   → Image-to-image variation is NOT statistically significant\n\n")
      }
    }
  } else {
    # Simple random effects structure - FIX HERE TOO
    cat("\nSimple random effects structure detected. Testing Image variance:\n\n")
    
    # Refit with ML for likelihood ratio test
    LogFit_ML <- update(LogFit, REML = FALSE)
    fit_log_fixed <- lm(logResults ~ Treatment, data = process_data)
    
    # Calculate likelihood ratio manually
    ll_mixed <- logLik(LogFit_ML)
    ll_fixed <- logLik(fit_log_fixed)
    chi_sq <- 2 * (as.numeric(ll_mixed) - as.numeric(ll_fixed))
    df_diff <- 1  # One random effect parameter
    p_val <- pchisq(chi_sq, df = df_diff, lower.tail = FALSE)
    
    cat("Test for Image-level variance:\n")
    cat("H0: Image variance = 0\n")
    cat("Chi-square = ", round(chi_sq, 4), ", df = ", df_diff, 
        ", p = ", format.pval(p_val, digits = 4), "\n")
    
    if (!is.na(p_val) && p_val < 0.05) {
      cat("→ Image-to-image variation IS statistically significant ***\n\n")
    } else {
      cat("→ Image-to-image variation is NOT statistically significant\n\n")
    }
  }
  
  # Calculate ICC for LogFit - ALWAYS RUN THIS
  cat("\n================================================================\n")
  cat("INTRACLASS CORRELATION COEFFICIENTS (ICC)\n")
  cat("================================================================\n")
  cat("ICC quantifies the proportion of variance explained by grouping structure.\n")
  cat("Values closer to 0 indicate little clustering; closer to 1 indicate strong clustering.\n\n")
  
  vc_log <- as.data.frame(VarCorr(LogFit))
  total_var_log <- sum(vc_log$vcov)
  
  for (i in 1:nrow(vc_log)) {
    if (vc_log$grp[i] != "Residual") {
      icc_log <- vc_log$vcov[i] / total_var_log
      cat(sprintf("%s ICC = %.6f (%.3f%% of total variance)\n", 
                  vc_log$grp[i], icc_log, icc_log * 100))
      if (icc_log < 0.01) {
        cat("   → Very weak clustering (negligible - random effect may not be necessary)\n\n")
      } else if (icc_log < 0.05) {
        cat("   → Weak clustering (minimal impact)\n\n")
      } else if (icc_log < 0.10) {
        cat("   → Moderate clustering\n\n")
      } else {
        cat("   → Strong clustering (substantial impact)\n\n")
      }
    }
  }
  cat(sprintf("Residual variance = %.4f (%.3f%% of total)\n\n", 
              vc_log$vcov[nrow(vc_log)], vc_log$vcov[nrow(vc_log)] / total_var_log * 100))
  
  # Add sequential pairwise comparisons for LogFit
  cat("\n================================================================\n")
  cat("SEQUENTIAL PAIRWISE COMPARISONS (LogFit)\n")
  cat("================================================================\n\n")
  
  emm_log <- emmeans(LogFit, ~ Treatment)
  
  sequential_contrasts_log <- list()
  contrast_names_log <- c()
  
  for (i in 2:n_treatments) {
    contrast_vec <- rep(0, n_treatments)
    contrast_vec[i - 1] <- 1
    contrast_vec[i] <- -1
    
    contrast_name <- paste(treatment_levels[i - 1], "vs", treatment_levels[i])
    contrast_names_log <- c(contrast_names_log, contrast_name)
    sequential_contrasts_log[[contrast_name]] <- contrast_vec
  }
  
  sequential_results_log <- contrast(emm_log, method = sequential_contrasts_log, adjust = "none")
  
  cat("Comparing consecutive treatment groups:\n")
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
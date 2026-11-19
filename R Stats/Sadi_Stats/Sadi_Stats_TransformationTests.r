# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","MASS","pbkrtest")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Disable pager to prevent "Selection:" prompts
options(pager = "console")

# Working directory
setwd("/Users/allisonpickle/Desktop/Sadi_Data/Crack Coverage")

# Data
input_file  <- "Crack coverage_225C.csv"
process_data <- read.csv(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "as-dealloyed"

# Treatment order and labels
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
  stop("ERROR: No 'Sample' column found in CSV.")
}

# ---- TRANSFORMATIONS ----

cat("\n=== APPLYING TRANSFORMATIONS ===\n")

# Box-Cox transformation
bc <- boxcox(Results ~ Treatment, data = process_data, plotit = FALSE)
lambda <- bc$x[which.max(bc$y)]
if (lambda == 0) {
  process_data$boxcoxResults <- log(process_data$Results + 0.5)
} else {
  process_data$boxcoxResults <- ((process_data$Results + 0.5)^lambda - 1) / lambda
}
cat("Box-Cox optimal lambda:", round(lambda, 3), "\n")

# Log transformation
process_data$logResults <- log(process_data$Results + 0.5)
cat("Log transformation: log(Results + 0.5)\n")

# Square root transformation
process_data$sqrtResults <- sqrt(process_data$Results + 0.5)
cat("Square root transformation: sqrt(Results + 0.5)\n")

# Cube root transformation
process_data$cbrtResults <- sign(process_data$Results) * abs(process_data$Results)^(1/3)
cat("Cube root transformation: sign(Results) * |Results|^(1/3)\n")

# Inverse transformation
process_data$invResults <- 1 / (process_data$Results + 1)
cat("Inverse transformation: 1 / (Results + 1)\n\n")

# Fit all models (quietly)
cat("Fitting models...\n")
raw_formula <- Results ~ Treatment + (1|Sample) + (1|Image)

RawFit    <- lmerTest::lmer(raw_formula, data = process_data)
BoxCoxFit <- lmerTest::lmer(update(raw_formula, boxcoxResults ~ .), data = process_data)
LogFit    <- lmerTest::lmer(update(raw_formula, logResults ~ .), data = process_data)
SqrtFit   <- lmerTest::lmer(update(raw_formula, sqrtResults ~ .), data = process_data)
CbrtFit   <- lmerTest::lmer(update(raw_formula, cbrtResults ~ .), data = process_data)
InvFit    <- lmerTest::lmer(update(raw_formula, invResults ~ .), data = process_data)

cat("All models fitted successfully!\n")

# ---- Function to compare normality ----
compare_normality <- function() {
  transformations <- list(
    Raw = process_data$Results,
    BoxCox = process_data$boxcoxResults,
    Log = process_data$logResults,
    Sqrt = process_data$sqrtResults,
    CubeRoot = process_data$cbrtResults,
    Inverse = process_data$invResults
  )
  
  results <- data.frame(
    Transformation = names(transformations),
    Shapiro_W = NA,
    Shapiro_p = NA,
    Skewness = NA,
    Kurtosis = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(transformations)) {
    x <- transformations[[i]]
    x <- x[!is.na(x) & !is.infinite(x)]
    
    # Shapiro test (sample if too large)
    if (length(x) > 5000) {
      x_sample <- sample(x, 5000)
      sw <- shapiro.test(x_sample)
      results$Shapiro_W[i] <- sw$statistic
      results$Shapiro_p[i] <- sw$p.value
    } else if (length(x) > 3) {
      sw <- shapiro.test(x)
      results$Shapiro_W[i] <- sw$statistic
      results$Shapiro_p[i] <- sw$p.value
    }
    
    results$Skewness[i] <- mean((x - mean(x))^3) / sd(x)^3
    results$Kurtosis[i] <- mean((x - mean(x))^4) / sd(x)^4 - 3
  }
  
  # Sort by Shapiro p-value (higher is better)
  results <- results[order(-results$Shapiro_p), ]
  
  cat("\n=== NORMALITY COMPARISON (RANKED BY BEST TO WORST) ===\n")
  cat("Higher Shapiro p-value = better normality (want p > 0.05)\n")
  cat("Skewness close to 0 = symmetric distribution\n")
  cat("Kurtosis close to 0 = normal tail behavior\n\n")
  print(results, row.names = FALSE, digits = 4)
  
  return(results)
}

# Run normality comparison
normality_results <- compare_normality()

# ---- Save diagnostics ONLY ----
save_diagnostics <- function(save_data = TRUE) {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("Normal_Diagnostics_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  models <- list(
    RawFit = RawFit,
    BoxCoxFit = BoxCoxFit,
    LogFit = LogFit,
    SqrtFit = SqrtFit,
    CbrtFit = CbrtFit,
    InvFit = InvFit
  )
  
  all_resid <- list()
  
  cat("\n=== GENERATING DIAGNOSTIC PLOTS ===\n")
  
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
         main = paste(nm, "Std Residuals vs Fitted"),
         pch = 16, col = rgb(0, 0.4, 0.9, 0.3))
    abline(h = 0, col = "red", lwd = 2)
    dev.off()
    
    ## 2) QQ plot
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
    
    cat("  -", nm, "plots saved\n")
  }
  
  if (save_data) {
    write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
  
  # Save normality comparison to text file
  sink(file.path(out_dir, "Normality_Comparison.txt"))
  cat("================================================================\n")
  cat("NORMALITY COMPARISON (RANKED BEST TO WORST)\n")
  cat("================================================================\n\n")
  cat("Higher Shapiro p-value = better normality (want p > 0.05)\n")
  cat("Skewness close to 0 = symmetric distribution\n")
  cat("Kurtosis close to 0 = normal tail behavior\n\n")
  print(normality_results, row.names = FALSE, digits = 4)
  cat("\n\n================================================================\n")
  cat("RECOMMENDATION\n")
  cat("================================================================\n\n")
  best_model <- normality_results$Transformation[1]
  cat("Based on normality tests, the best transformation is: ", best_model, "\n\n", sep = "")
  cat("ALWAYS check the diagnostic plots to confirm:\n")
  cat("  1. QQ plot: Points should follow the red line\n")
  cat("  2. Residuals plot: Random scatter around zero\n")
  cat("  3. No patterns or funneling in residuals\n")
  sink()
  
  message("\nSaved diagnostic plots and normality comparison to: ", out_dir)
}

# ---- Run diagnostics ----
save_diagnostics()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Tested transformations: Raw, Box-Cox (λ=", round(lambda, 3), "), Log, Sqrt, CubeRoot, Inverse\n", sep = "")
cat("Best transformation (by normality): ", normality_results$Transformation[1], "\n", sep = "")
cat("\n⚠️  IMPORTANT: Always check the diagnostic plots!\n")
cat("   - QQ plots show normality of residuals\n")
cat("   - Residual plots show homoscedasticity\n\n")
cat("All outputs saved to: MM_Diagnostics_", tools::file_path_sans_ext(basename(input_file)), "\n", sep = "")
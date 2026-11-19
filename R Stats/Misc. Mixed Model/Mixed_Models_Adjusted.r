# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Working directory
setwd("/Users/allisonpickle/Desktop")

# Data
input_file  <- "PD_Microglia Area_RMixedFormated.csv"
process_data <- read.csv(input_file)

# Factors
# Hard‑code the first (baseline) level here:
baseline_treatment <- "Ctrl"

process_data$Treatment <- factor(
  process_data$Treatment,
  levels = c(
    baseline_treatment,
    sort(setdiff(unique(process_data$Treatment), baseline_treatment))
  )
)
process_data$Well <- factor(process_data$Well)

# Model 1 (raw)
RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Treatment:Well), data = process_data)
anova(RawFit); summary(RawFit)

# macOS graphics device (optional)
if (Sys.info()[["sysname"]] == "Darwin" && capabilities("aqua")) {
  if (is.null(dev.list()) || !any(grepl("quartz", names(dev.list())))) quartz()
}

# Helper plotting function: ONLY std residuals vs fitted + QQ
plot_diag <- function(fit, tag) {
  par(mfrow = c(1, 2))

  ## 1) STANDARDIZED (Pearson) residuals vs fitted
  fitted_vals <- fitted(fit)
  std_resid   <- resid(fit, type = "pearson")
  plot(fitted_vals, std_resid,
       xlab = "Fitted Values",
       ylab = "Standardized (Pearson) Residuals",
       main = paste(tag, "Std Residuals vs Fitted"))
  abline(h = 0, col = "red", lwd = 2)

  ## 2) QQ plot
  print(lattice::qqmath(fit, id = 0.05,
                        main = paste(tag, "QQ Plot")))

  par(mfrow = c(1, 1))
}

# Diagnostics for raw model
plot_diag(RawFit, "Raw")

# Transform: log
process_data$logResults <- log(process_data$Results + 0.5)
LogFit <- lmerTest::lmer(logResults ~ Treatment + (1|Treatment:Well), data = process_data)
anova(LogFit); summary(LogFit)
plot_diag(LogFit, "Log")

save_diagnostics()

# ---- Save ANOVA and summary outputs for RawFit and LogFit ----
save_model_summaries <- function() {
  # Build output dir name using input CSV base name (same as diagnostics)
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)

  out_file <- file.path(out_dir, "Model_Summaries.txt")

  sink(out_file)
  cat("=== RawFit (Results ~ Treatment + (1|Treatment:Well)) ===\n\n")
  cat("ANOVA:\n")
  print(anova(RawFit))
  cat("\nSummary:\n")
  print(summary(RawFit))

  cat("\n\n=== LogFit (logResults ~ Treatment + (1|Treatment:Well)) ===\n\n")
  cat("ANOVA:\n")
  print(anova(LogFit))
  cat("\nSummary:\n")
  print(summary(LogFit))
  sink()

  message("Saved ANOVA and summary outputs to: ", out_file)
}

save_model_summaries()

# Save diagnostics: ONLY std_resid vs fitted + QQ, plus residual data
save_diagnostics <- function(save_data = TRUE) {
  # Build output dir name using input CSV base name
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base) # Mixed Models Outputs
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
    
    ## 2) QQ plot
    png(file.path(out_dir, paste0(nm, "_qq.png")),
        width = 1200, height = 1000, res = 150)
    print(lattice::qqmath(fit, id = 0.05,
                          main = paste(nm, "QQ Plot")))
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

save_diagnostics()
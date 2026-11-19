# Packages (quiet)
suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(ggplot2)
  library(lattice)
  library(dplyr)
})

# Silence warnings globally (keeps only explicit messages you print)
options(warn = -1)

# ============================================================
# RUN SETTINGS
# ============================================================

# ---- NEW: point this at a folder containing one or more *_RMixedFormated.csv files ----
input_dir <- '/Users/allisonpickle/Desktop/Neuronal AB Surface Area Coverage'

# Optional: pattern to select just the files you want
file_pattern <- "RMixedFormated.csv$"

# ----------------- FACTORS -----------------
baseline_treatment <- "SD-pup"

treatment_labels <- c(
  "SD-pup" = "+\n-\n-\n-",
  "F344WT" = "+\n+\n-\n-",
  "F344TgAD" = "+\n-\n+\n-",
  "F344TgAD_young" = "+\n-\n-\n+"
)

# ----------------- PLOT LABEL INPUTS (user-defined) -----------------
plot_title_base <- "Neuronal Surface Area Coverage"
x_axis_label    <- ""                # keep "" to match your current style
y_axis_label_expr <- expression(paste("Surface Area Coverage (%)"))

treatment_levels <- names(treatment_labels)
row_labels <- c("SD-pup", "F344WT", "F344TgAD", "F344TgAD_Young")

# ============================================================
# FUNCTIONS (unchanged logic, wrapped so it can run per-file)
# ============================================================

fit_time_model <- function(df, response = c("Results", "logResults"), treatment_levels) {
  response <- match.arg(response)

  if (!("Biological_Rep" %in% names(df))) stop("Biological_Rep column missing.")
  if (!("Well" %in% names(df))) stop("Well column missing.")
  if (!("Treatment" %in% names(df))) stop("Treatment column missing.")
  if (!(response %in% names(df))) stop(paste0(response, " column missing."))

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well           <- droplevels(factor(df$Well))
  df$Treatment      <- droplevels(factor(df$Treatment, levels = treatment_levels))

  has_img <- "Image" %in% names(df)
  if (has_img) df$Image <- droplevels(factor(df$Image))

  n_bio <- nlevels(droplevels(df$Biological_Rep))

  if (n_bio >= 2 && has_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well) + (1 | Well:Image)"
    ))
  } else if (n_bio >= 2) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
    ))
  } else if (has_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Well) + (1 | Well:Image)"
    ))
  } else {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + (1 | Well)"
    ))
  }

  suppressMessages(suppressWarnings(
    lmerTest::lmer(f, data = df)
  ))
}

make_null_full <- function(fit, response_name, df) {
  base_form <- formula(fit)
  rand_terms <- lme4::findbars(base_form)

  make_form <- function(fixed_rhs) {
    if (length(rand_terms) == 0) {
      stats::as.formula(paste0(response_name, " ~ ", fixed_rhs))
    } else {
      rand_part <- paste(
        vapply(rand_terms, function(x) paste0("(", deparse1(x), ")"), character(1)),
        collapse = " + "
      )
      stats::as.formula(paste0(response_name, " ~ ", fixed_rhs, " + ", rand_part))
    }
  }

  null_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(make_form("1"), data = df, REML = FALSE)
  ))
  full_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(make_form("Treatment"), data = df, REML = FALSE)
  ))

  list(null = null_fit, full = full_fit)
}

compute_icc <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  var_random <- sum(vc$vcov[vc$grp != "Residual"])
  var_resid  <- vc$vcov[vc$grp == "Residual"]
  var_random / (var_random + var_resid)
}

extract_treatment_pvals <- function(fit, baseline_name) {
  fe <- coef(summary(fit))
  rn <- rownames(fe)

  is_trt <- grepl("^Treatment", rn)
  if (!any(is_trt)) stop("No Treatment* rows found in fixed effects.")

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

build_bar_summary_with_sig <- function(data, fit, response_col, baseline_name) {
  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  summary_df$Treatment <- factor(summary_df$Treatment, levels = levels(data$Treatment))

  p_df <- extract_treatment_pvals(fit, baseline_name)
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars <- vapply(out$p.value, add_stars, character(1))
  out
}

plot_bar_p_values <- function(
  summ_df,
  raw_data,
  response_col,
  y_label,
  title_text,
  x_label,
  treatment_labels,
  row_labels
) {
  aligned_treatment_labels <- setNames(paste0("\n", treatment_labels), names(treatment_labels))
  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))

  ggplot(summ_df, aes(x = Treatment, y = mean)) +
    geom_col(fill = "#6c6c6c", alpha = 0.7) +
    geom_jitter(
      data = raw_data,
      aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep, color = Biological_Rep),
      width = 0.2, alpha = 0.6, size = 1.0
    ) +
    scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 7)) +
    scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00", "black", "black", "black", "black", "black")) +
    scale_x_discrete(labels = aligned_treatment_labels) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(aes(y = mean + se, label = stars), vjust = -0.5, size = 7, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(x = x_label, y = y_label, title = title_text, shape = "Biological Rep", color = "Biological Rep") +
    annotate("text", x = 0.4, y = -Inf, label = combined_row_label, hjust = 1, vjust = 1.1, size = 3, lineheight = 0.8) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.8, size = 9, margin = margin(t = 0)),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(5, 5, 15, 50),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    ) +
    coord_cartesian(clip = "off")
}

save_diagnostics <- function(out_dir, NestedFit, LogNestedFit, save_data = TRUE) {
  models <- list(NestedFit = NestedFit, LogNestedFit = LogNestedFit)
  all_resid <- list()

  for (nm in names(models)) {
    fit <- models[[nm]]
    fitted_vals <- fitted(fit)
    std_resid   <- resid(fit, type = "pearson")

    png(file.path(out_dir, paste0(nm, "_residuals.png")), width = 1200, height = 1000, res = 150)
    plot(fitted_vals, std_resid, xlab = "Fitted", ylab = "Pearson resid", main = paste(nm, "Residuals vs Fitted"))
    abline(h = 0, col = "red", lwd = 2)
    dev.off()

    png(file.path(out_dir, paste0(nm, "_qq.png")), width = 1200, height = 1000, res = 150)
    qq_data <- qqnorm(resid(fit), plot.it = FALSE)
    plot(qq_data$x, qq_data$y, xlab = "Theoretical", ylab = "Sample", main = paste(nm, "QQ Plot"), pch = 16)
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
    write.csv(do.call(rbind, all_resid), file.path(out_dir, "All_Calculated_Residuals.csv"), row.names = FALSE)
  }
}

save_model_summaries <- function(
  out_dir,
  input_file,
  analysis_data,
  NestedFit,
  LogNestedFit,
  icc_nested,
  icc_log_nested,
  lrt_nested,
  lrt_log_nested
) {
  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")
  sink(out_file)

  cat("MIXED EFFECTS MODEL SUMMARIES\n\n")
  cat("Input: ", input_file, "\n")
  cat("Rows used: ", nrow(analysis_data), "\n")
  cat("Biological reps: ", paste(sort(unique(analysis_data$Biological_Rep)), collapse = ", "), "\n\n")

  cat("NestedFit (raw)\n")
  print(summary(NestedFit))
  cat("\nICC (raw): ", sprintf("%.4f", icc_nested), "\n\n")
  cat("LRT (raw):\n")
  if (is.null(lrt_nested)) cat("Not available (comparison failed).\n") else print(lrt_nested)

  cat("\n\nLogNestedFit\n")
  print(summary(LogNestedFit))
  cat("\nICC (log): ", sprintf("%.4f", icc_log_nested), "\n\n")
  cat("LRT (log):\n")
  if (is.null(lrt_log_nested)) cat("Not available (comparison failed).\n") else print(lrt_log_nested)

  sink()
}

run_one_csv <- function(input_file) {
  # ----------------- INPUT -----------------
  process_data <- read.csv(input_file)

  # --- FORCE expected column names by position: A-F ---
  if (ncol(process_data) < 6) stop("CSV must have at least 6 columns (A-F).")
  names(process_data)[1:6] <- c("Filename", "Treatment", "Biological_Rep", "Well", "Image", "Results")

  # Coerce types (quiet)
  process_data$Treatment      <- as.character(process_data$Treatment)
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)
  process_data$Well           <- factor(process_data$Well)
  process_data$Image          <- factor(process_data$Image)
  process_data$Results        <- suppressWarnings(as.numeric(process_data$Results))

  # ----------------- OUTPUT (same folder as CSV, but per input file) -----------------
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ----------------- FACTORS -----------------
  process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)

  # ----------------- REQUIRED COLUMNS (minimal checks; no extra output) -----------------
  if (!("Biological_Rep" %in% names(process_data))) stop("ERROR: No 'Biological_Rep' column found in CSV.")
  has_image <- "Image" %in% names(process_data)

  # ----------------- BIO REP CHECK (ONLY TERMINAL OUTPUT) -----------------
  n_bioreps <- nlevels(droplevels(process_data$Biological_Rep))
  if (n_bioreps == 1) {
    cat("NOTE: Only 1 biological replicate detected in ", basename(input_file), ". Random BioRep effects not estimable.\n", sep = "")
  }

  # ----------------- COMPLETE CASES (single dataset for all fits) -----------------
  needed <- c("Results", "Treatment", "Well", "Biological_Rep")
  if (has_image) needed <- c(needed, "Image")
  analysis_data <- process_data[stats::complete.cases(process_data[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(analysis_data) == 0) stop("No rows left after complete.cases()")

  # Drop unused levels AFTER filtering
  analysis_data$Treatment      <- droplevels(factor(analysis_data$Treatment, levels = treatment_levels))
  analysis_data$Biological_Rep <- droplevels(factor(analysis_data$Biological_Rep))
  analysis_data$Well           <- droplevels(factor(analysis_data$Well))
  if ("Image" %in% names(analysis_data)) analysis_data$Image <- droplevels(factor(analysis_data$Image))

  # ----------------- MODEL FITS (quiet) -----------------
  NestedFit <- fit_time_model(analysis_data, response = "Results", treatment_levels = treatment_levels)

  analysis_data$logResults <- log(analysis_data$Results + 0.5)
  LogNestedFit <- fit_time_model(analysis_data, response = "logResults", treatment_levels = treatment_levels)

  # ----------------- Null/full for LRT -----------------
  nf_raw <- make_null_full(NestedFit, "Results", analysis_data)
  NestedFit_null <- nf_raw$null
  NestedFit_full <- nf_raw$full

  nf_log <- make_null_full(LogNestedFit, "logResults", analysis_data)
  LogNestedFit_null <- nf_log$null
  LogNestedFit_full <- nf_log$full

  # ----------------- ICC -----------------
  icc_nested <- compute_icc(NestedFit)
  icc_log_nested <- compute_icc(LogNestedFit)

  # ----------------- LRT -----------------
  lrt_nested <- tryCatch(anova(NestedFit_null, NestedFit_full), error = function(e) NULL)
  lrt_log_nested <- tryCatch(anova(LogNestedFit_null, LogNestedFit_full), error = function(e) NULL)

  # ----------------- BARPLOTS -----------------
  nested_bar_summary <- build_bar_summary_with_sig(analysis_data, NestedFit, "Results", baseline_treatment)
  nested_bar_plot <- plot_bar_p_values(
    nested_bar_summary,
    analysis_data,
    "Results",
    y_label    = bquote(.(y_axis_label_expr)~"(mean \u00B1 SE)"),
    title_text = plot_title_base,
    x_label    = x_axis_label,
    treatment_labels = treatment_labels,
    row_labels = row_labels
  )
  ggsave(
    filename = file.path(out_dir, "Barplot_NestedFit.png"),
    plot = nested_bar_plot,
    width = 6, height = 4, dpi = 300
  )

  log_nested_bar_summary <- build_bar_summary_with_sig(analysis_data, LogNestedFit, "logResults", baseline_treatment)
  log_nested_bar_plot <- plot_bar_p_values(
    log_nested_bar_summary,
    analysis_data,
    "logResults",
    y_label    = bquote(log(.(y_axis_label_expr) + 0.5)~"(mean \u00B1 SE)"),
    title_text = plot_title_base,
    x_label    = x_axis_label,
    treatment_labels = treatment_labels,
    row_labels = row_labels
  )
  ggsave(
    filename = file.path(out_dir, "Barplot_LogNestedFit.png"),
    plot = log_nested_bar_plot,
    width = 6, height = 4, dpi = 300
  )

  # ----------------- SAVES -----------------
  save_diagnostics(out_dir, NestedFit, LogNestedFit, save_data = TRUE)
  save_model_summaries(
    out_dir,
    input_file,
    analysis_data,
    NestedFit,
    LogNestedFit,
    icc_nested,
    icc_log_nested,
    lrt_nested,
    lrt_log_nested
  )

  invisible(list(
    input_file = input_file,
    out_dir = out_dir,
    n_rows = nrow(analysis_data)
  ))
}

# ============================================================
# RUN ALL CSVs IN A FOLDER
# ============================================================

csv_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", input_dir)

cat("Found ", length(csv_files), " file(s) in ", input_dir, "\n", sep = "")

results <- lapply(csv_files, function(f) {
  cat("\n--- Running: ", basename(f), " ---\n", sep = "")
  run_one_csv(f)
})
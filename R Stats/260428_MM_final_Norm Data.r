# ============================================================
# Mixed Effects (Normalized): Treatment-only with adaptive nested random effects
# - Uses Results_Norm as the response (col 7)
# - ALSO analyzes Results_log2FC as the response (col 8)
# - Fits:
#     * Results_Norm (raw only)
#     * Results_log2FC (raw only)
# - Diagnostics, ICC, LRT, barplots
# - Runs for ALL matching *_Norm.csv files in input_dir
# ============================================================

suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(ggplot2)
  library(dplyr)
})

options(warn = -1)

# ============================================================
# RUN SETTINGS
# ============================================================
input_dir <- '/Users/allisonpickle/Desktop/Normalized LPS Only/Microglia Morphology by Incubation Time_Split/48hr'
file_pattern <- "\\_Norm\\.csv$"   # only process normalized CSVs

# expected columns are already present; disable positional renaming to avoid duplicates
# Expected schema:
# 1 Filename, 2 Treatment, 3 Biological_Rep, 4 Well, 5 Image, 6 Results, 7 Results_Norm, 8 Results_log2FC
rename_cols_by_position <- NULL

# ----------------- FACTORS / LABELS -----------------
baseline_treatment <- "48hr_Ctrl"

treatment_levels <- c(
  "48hr_Ctrl",
  "48hr_5ng",
  "48hr_100ng",
  "48hr_500ng",
  "48hr_5ug"
)

treatment_labels <- c(
  "48hr_Ctrl"      = "-\n-\n-\n-",
  "48hr_5ng"        = "+\n-\n-\n-",
  "48hr_100ng"      = "-\n+\n-\n-",
  "48hr_500ng"      = "-\n-\n+\n-",
  "48hr_5ug"        = "-\n-\n-\n+"
)

row_labels <- c(
  "5 ng/mL LPS",
  "100 ng/mL LPS",
  "500 ng/mL LPS",
  "5 µg/mL LPS"
)

plot_title_base <- "Treatment"
x_axis_label <- ""
y_axis_label_raw <- "Results_Norm (mean ± SE)"

# log2FC (already log2-transformed; no asinh)
y_axis_label_log2fc_raw <- "log2FC vs 6hr_Ctrl median (mean ± SE)"

# ============================================================
# EDIT: helper for transparent PNG export
# ============================================================
save_png_transparent <- function(filename, plot, width, height, dpi) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "transparent"
  )
}

# ============================================================
# HELPERS
# ============================================================
add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

make_aligned_tick_labels <- function(levels_vec, labels_named_vec) {
  out <- rep("", length(levels_vec))
  names(out) <- levels_vec
  common <- intersect(levels_vec, names(labels_named_vec))
  out[common] <- as.character(labels_named_vec[common])
  out
}

# Fit adaptive random-effects structure for treatment-only
fit_treatment_model <- function(df, response) {
  if (!("Biological_Rep" %in% names(df))) stop("Biological_Rep column missing.")
  if (!("Well" %in% names(df))) stop("Well column missing.")
  if (!("Treatment" %in% names(df))) stop("Treatment column missing.")
  if (!(response %in% names(df))) stop(paste0(response, " column missing."))

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well <- droplevels(factor(df$Well))
  df$Treatment <- droplevels(factor(df$Treatment, levels = treatment_levels))

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

# Build null and full models (same random-effects as fit), ML for LRT
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

# Single ICC (random / (random + resid))
compute_icc <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  var_resid <- attr(lme4::VarCorr(fit), "sc")^2
  var_random <- sum(vc$vcov, na.rm = TRUE)
  var_random / (var_random + var_resid)
}

# baseline-vs-others fixed-effect p-values (assumes relevelled Treatment)
extract_treatment_pvals <- function(fit, levels_vec) {
  fe <- coef(summary(fit))
  rn <- rownames(fe)
  is_trt <- grepl("^Treatment", rn) & !grepl(":", rn)
  if (!any(is_trt)) stop("No Treatment* rows found in fixed effects.")

  trt_rows  <- fe[is_trt, , drop = FALSE]
  trt_names <- sub("^Treatment", "", rownames(trt_rows))

  p_col <- grep("Pr\\(", colnames(trt_rows), value = TRUE)[1]
  if (is.na(p_col)) stop("Could not find p-value column in fixed effects.")

  pvals <- trt_rows[, p_col]
  names(pvals) <- trt_names

  p_df <- data.frame(Treatment = levels_vec, p.value = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(trt_names)) {
    p_df$p.value[p_df$Treatment == trt_names[i]] <- as.numeric(pvals[i])
  }
  p_df
}

build_bar_summary_with_sig <- function(data, fit, response_col) {
  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  summary_df$Treatment <- factor(summary_df$Treatment, levels = levels(data$Treatment))

  p_df <- extract_treatment_pvals(fit, levels(data$Treatment))
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars <- vapply(out$p.value, add_stars, character(1))
  out <- out[!is.na(out$mean), , drop = FALSE]
  out
}

plot_bar_p_values <- function(summ_df, raw_data, response_col, y_label, title_text) {
  raw_data <- raw_data[raw_data$Treatment %in% summ_df$Treatment & !is.na(raw_data$Biological_Rep), , drop = FALSE]
  raw_data$Treatment <- droplevels(raw_data$Treatment)
  raw_data$Biological_Rep <- droplevels(raw_data$Biological_Rep)

  summ_df$Treatment <- factor(as.character(summ_df$Treatment), levels = treatment_levels)
  raw_data$Treatment <- factor(as.character(raw_data$Treatment), levels = treatment_levels)

  aligned_treatment_labels <- make_aligned_tick_labels(treatment_levels, treatment_labels)
  aligned_treatment_labels <- stats::setNames(paste0("\n", aligned_treatment_labels), names(aligned_treatment_labels))

  combined_row_label <- ""
  if (exists("row_labels", inherits = TRUE) && length(row_labels) > 0) {
    combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  }

  hide_legend <- (nlevels(raw_data$Biological_Rep) <= 1)

  p <- ggplot2::ggplot(summ_df, ggplot2::aes(x = Treatment, y = mean)) +
    ggplot2::geom_col(fill = "#6c6c6c", alpha = 0.7) +
    ggplot2::geom_jitter(
      data = raw_data,
      ggplot2::aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep),
      width = 0.2, alpha = 0.6, size = 1.2,
      color = "#1f5fbf"
    ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(y = mean + se, label = stars),
      vjust = -0.5, size = 7, fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(
      title = title_text,
      x = x_axis_label,
      y = y_label,
      shape = "Biological Rep"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      # EDIT: transparent backgrounds for PNG export
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),

      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = -0.2, lineheight = 0.8, size = 9),
      axis.ticks.x = ggplot2::element_line(),
      axis.title.x = ggplot2::element_blank(),
      legend.position = if (hide_legend) "none" else "right",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(5, 5, 45, 55)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!identical(combined_row_label, "")) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0.4, y = -Inf,
      label = combined_row_label,
      hjust = 1, vjust = 1.3,
      size = 3, color = "black",
      lineheight = 0.8
    )
  }

  p
}

# ============================================================
# DIAGNOSTICS
# ============================================================
save_diagnostics <- function(out_dir, fit, save_data = TRUE) {
  nm <- "RawFit"
  all_resid <- list()

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fitted_vals   <- as.numeric(stats::fitted(fit))
  pearson_resid <- as.numeric(stats::resid(fit, type = "pearson"))
  raw_resid     <- as.numeric(stats::resid(fit))

  # ----------------- Residuals vs Fitted -----------------
  res_df <- data.frame(
    FittedValues = fitted_vals,
    PearsonResiduals = pearson_resid
  )

  p_res <- ggplot2::ggplot(res_df, ggplot2::aes(x = FittedValues, y = PearsonResiduals)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linewidth = 0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0(nm, ": Residuals vs Fitted"),
      x = "Fitted Values",
      y = "Standardized (Pearson) Residuals"
    ) +
    ggplot2::theme(
      # EDIT: transparent backgrounds for PNG export
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA)
    )

  save_png_transparent(
    filename = file.path(out_dir, paste0(nm, "_residuals.png")),
    plot = p_res,
    width = 7, height = 5, dpi = 200
  )

  # ----------------- QQ plot -----------------
  qq <- stats::qqnorm(raw_resid, plot.it = FALSE)
  qq_df <- data.frame(TheoreticalQuantiles = qq$x, SampleQuantiles = qq$y)

  qy <- stats::quantile(raw_resid, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  qx <- stats::qnorm(c(0.25, 0.75))
  slope <- (qy[2] - qy[1]) / (qx[2] - qx[1])
  intercept <- qy[1] - slope * qx[1]

  p_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = TheoreticalQuantiles, y = SampleQuantiles)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::geom_abline(slope = slope, intercept = intercept, color = "red", linewidth = 0.7) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0(nm, ": QQ Plot"),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme(
      # EDIT: transparent backgrounds for PNG export
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA)
    )

  save_png_transparent(
    filename = file.path(out_dir, paste0(nm, "_qq.png")),
    plot = p_qq,
    width = 7, height = 5, dpi = 200
  )

  if (save_data) {
    all_resid[[nm]] <- data.frame(
      Model = nm,
      FittedValues = fitted_vals,
      PearsonResiduals = pearson_resid,
      RawResiduals = raw_resid
    )

    utils::write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# SAVE MODEL SUMMARIES (raw only)
# ============================================================
save_model_summaries <- function(
  out_dir,
  input_file,
  analysis_data,
  raw_fit,
  icc_raw,
  lrt_raw
) {
  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")
  sink(out_file)

  cat("MIXED EFFECTS MODEL SUMMARY (Treatment only; normalized response)\n\n")
  cat("Input: ", input_file, "\n", sep = "")
  cat("Rows used: ", nrow(analysis_data), "\n", sep = "")
  cat("Biological reps: ", paste(sort(unique(analysis_data$Biological_Rep)), collapse = ", "), "\n\n", sep = "")

  cat("RAW model (Results_Norm)\n")
  print(formula(raw_fit))
  cat("\n")
  print(summary(raw_fit))
  cat("\nICC (raw): ", sprintf("%.4f", icc_raw), "\n\n", sep = "")

  cat("LRT (raw; ML null vs full)\n")
  if (is.null(lrt_raw)) cat("Not available (comparison failed).\n") else print(lrt_raw)

  sink()
}

# ============================================================
# MAIN PER-FILE RUNNER
# ============================================================
run_one_csv <- function(input_file) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("PROCESSING FILE: ", input_file, "\n", sep = "")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

  # keep original column names; don't auto-mangle
  process_data <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

  # enforce expected schema (only if not already present)
  expected <- c("Filename","Treatment","Biological_Rep","Well","Image","Results","Results_Norm","Results_log2FC")
  if (ncol(process_data) < 8) {
    stop("Expected at least 8 columns. Found: ", ncol(process_data), " in file: ", basename(input_file))
  }
  if (!all(expected %in% names(process_data))) {
    names(process_data)[1:8] <- expected
  }

  # ensure unique column names (dplyr requires this)
  if (any(duplicated(names(process_data)))) {
    dupes <- unique(names(process_data)[duplicated(names(process_data))])
    stop("Duplicate column names in input CSV: ", paste(dupes, collapse = ", "))
  }

  # Optional: rename columns by position (DISABLED above; left here for completeness)
  if (!is.null(rename_cols_by_position)) {
    needed <- max(rename_cols_by_position)
    if (ncol(process_data) < needed) {
      stop(
        "CSV has ", ncol(process_data),
        " columns but rename_cols_by_position needs at least ", needed,
        " columns. File: ", basename(input_file)
      )
    }
    if ("Results_Norm" %in% names(rename_cols_by_position)) {
      colnames(process_data)[rename_cols_by_position[["Results_Norm"]]] <- "Results_Norm"
    }
    if ("Results_log2FC" %in% names(rename_cols_by_position)) {
      colnames(process_data)[rename_cols_by_position[["Results_log2FC"]]] <- "Results_log2FC"
    }
    if (any(duplicated(names(process_data)))) {
      dupes <- unique(names(process_data)[duplicated(names(process_data))])
      stop("Duplicate column names after renaming: ", paste(dupes, collapse = ", "))
    }
  }

  # Required columns (analysis is on LAST TWO: Results_Norm and Results_log2FC)
  if (!("Treatment" %in% names(process_data))) stop("Missing required column: Treatment")
  if (!("Well" %in% names(process_data))) stop("Missing required column: Well")
  if (!("Results_Norm" %in% names(process_data))) stop("Missing required column: Results_Norm")
  if (!("Results_log2FC" %in% names(process_data))) stop("Missing required column: Results_log2FC")
  if (!("Biological_Rep" %in% names(process_data))) stop("Missing required column: Biological_Rep")

  # Types / factors
  process_data$Results_Norm <- suppressWarnings(as.numeric(process_data$Results_Norm))
  process_data$Results_log2FC <- suppressWarnings(as.numeric(process_data$Results_log2FC))
  process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)

  # baseline reference
  if (baseline_treatment %in% levels(process_data$Treatment)) {
    process_data$Treatment <- stats::relevel(process_data$Treatment, ref = baseline_treatment)
  }

  process_data$Well <- factor(process_data$Well)
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)

  has_image <- "Image" %in% names(process_data)
  if (has_image) process_data$Image <- factor(process_data$Image)

  # Filter complete cases used by models:
  needed_norm <- c("Results_Norm", "Treatment", "Well", "Biological_Rep")
  if (has_image) needed_norm <- c(needed_norm, "Image")

  analysis_data_norm <- process_data[stats::complete.cases(process_data[, needed_norm, drop = FALSE]), , drop = FALSE]
  analysis_data_norm <- analysis_data_norm[!is.na(analysis_data_norm$Treatment), , drop = FALSE]
  if (nrow(analysis_data_norm) == 0) stop("No rows left for Results_Norm after filtering complete cases / Treatment levels.")

  analysis_data_norm$Treatment <- droplevels(factor(analysis_data_norm$Treatment, levels = treatment_levels))
  analysis_data_norm$Well <- droplevels(factor(analysis_data_norm$Well))
  analysis_data_norm$Biological_Rep <- droplevels(factor(analysis_data_norm$Biological_Rep))
  if (has_image) analysis_data_norm$Image <- droplevels(factor(analysis_data_norm$Image))

  # Output directory
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  csv_base_clean <- sub("_Norm$", "", csv_base)
  out_dir  <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base_clean, "_Norm"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ============================================================
  # (A) Results_Norm: RAW ONLY
  # ============================================================
  raw_fit <- fit_treatment_model(analysis_data_norm, response = "Results_Norm")

  nf_raw <- make_null_full(raw_fit, "Results_Norm", analysis_data_norm)
  lrt_raw <- tryCatch(anova(nf_raw$null, nf_raw$full), error = function(e) NULL)

  icc_raw <- compute_icc(raw_fit)

  write.csv(
    data.frame(Model = "RawFit_Results_Norm", ICC = icc_raw),
    file.path(out_dir, "ICC_Summary_Results_Norm.csv"),
    row.names = FALSE
  )
  if (!is.null(lrt_raw)) capture.output(lrt_raw, file = file.path(out_dir, "LRT_Raw_Results_Norm_Treatment_vs_Null.txt"))

  summ_raw <- build_bar_summary_with_sig(analysis_data_norm, raw_fit, "Results_Norm")
  p_raw <- plot_bar_p_values(
    summ_df = summ_raw,
    raw_data = analysis_data_norm,
    response_col = "Results_Norm",
    y_label = y_axis_label_raw,
    title_text = paste0(plot_title_base, " (raw; Results_Norm)")
  )
  save_png_transparent(
    filename = file.path(out_dir, "Barplot_RawFit_Results_Norm.png"),
    plot = p_raw,
    width = 6, height = 4, dpi = 300
  )

  save_diagnostics(file.path(out_dir, "Diagnostics_Results_Norm"), raw_fit, save_data = TRUE)
  save_model_summaries(out_dir, input_file, analysis_data_norm, raw_fit, icc_raw, lrt_raw)

  # ============================================================
  # (B) Results_log2FC: RAW ONLY
  # ============================================================
  needed_fc <- c("Results_log2FC", "Treatment", "Well", "Biological_Rep")
  if (has_image) needed_fc <- c(needed_fc, "Image")

  analysis_data_fc <- process_data[stats::complete.cases(process_data[, needed_fc, drop = FALSE]), , drop = FALSE]
  analysis_data_fc <- analysis_data_fc[!is.na(analysis_data_fc$Treatment), , drop = FALSE]

  if (nrow(analysis_data_fc) > 0) {
    analysis_data_fc$Treatment <- droplevels(factor(analysis_data_fc$Treatment, levels = treatment_levels))
    analysis_data_fc$Well <- droplevels(factor(analysis_data_fc$Well))
    analysis_data_fc$Biological_Rep <- droplevels(factor(analysis_data_fc$Biological_Rep))
    if (has_image) analysis_data_fc$Image <- droplevels(factor(analysis_data_fc$Image))

    raw_fit_fc <- fit_treatment_model(analysis_data_fc, response = "Results_log2FC")

    nf_raw_fc <- make_null_full(raw_fit_fc, "Results_log2FC", analysis_data_fc)
    lrt_raw_fc <- tryCatch(anova(nf_raw_fc$null, nf_raw_fc$full), error = function(e) NULL)

    icc_raw_fc <- compute_icc(raw_fit_fc)

    write.csv(
      data.frame(Model = "RawFit_Results_log2FC", ICC = icc_raw_fc),
      file.path(out_dir, "ICC_Summary_Results_log2FC.csv"),
      row.names = FALSE
    )
    if (!is.null(lrt_raw_fc)) capture.output(lrt_raw_fc, file = file.path(out_dir, "LRT_Raw_Results_log2FC_Treatment_vs_Null.txt"))

    summ_raw_fc <- build_bar_summary_with_sig(analysis_data_fc, raw_fit_fc, "Results_log2FC")
    p_raw_fc <- plot_bar_p_values(
      summ_df = summ_raw_fc,
      raw_data = analysis_data_fc,
      response_col = "Results_log2FC",
      y_label = y_axis_label_log2fc_raw,
      title_text = paste0(plot_title_base, " (raw; log2FC)")
    )
    save_png_transparent(
      filename = file.path(out_dir, "Barplot_RawFit_Results_log2FC.png"),
      plot = p_raw_fc,
      width = 6, height = 4, dpi = 300
    )

    save_diagnostics(file.path(out_dir, "Diagnostics_log2FC"), raw_fit_fc, save_data = TRUE)
  } else {
    message("No complete-case rows available for Results_log2FC in: ", basename(input_file), " (skipping log2FC).")
  }

  cat("Done. Outputs saved to: ", out_dir, "\n", sep = "")
  invisible(list(input_file = input_file, out_dir = out_dir, n_rows_norm = nrow(analysis_data_norm)))
}

# ============================================================
# RUN ALL CSVs IN A FOLDER
# ============================================================
csv_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", input_dir)

cat("Found ", length(csv_files), " file(s) in ", input_dir, "\n", sep = "")

results <- lapply(csv_files, function(f) {
  tryCatch(
    run_one_csv(f),
    error = function(e) {
      message("FAILED: ", basename(f), " -> ", conditionMessage(e))
      NULL
    }
  )
})
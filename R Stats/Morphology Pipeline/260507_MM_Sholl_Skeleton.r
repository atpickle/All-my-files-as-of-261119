# ============================================================
# Mixed Effects (Sholl AUC): Treatment × Method with adaptive nested random effects
# + Fits raw + log(Results + 0.5)
# + Diagnostics: residual-vs-fitted + QQ
# + ICC (single ICC = var(random) / (var(random) + var(resid)))
# + LRTs (ML fits):
#     - Treatment (null vs Treatment)
#     - Method (null vs Method)
#     - Interaction (additive vs interaction)
# + Post hoc (requested): within each Method, Ctrl vs 5ugLPS contrast
#     - Uses emmeans (marginal means) with Kenward-Roger df
#     - Holm adjust across the 3 methods
# + Plots (like your Python function):
#     - For each Method, matching raw + natural-log Ctrl vs 5ugLPS plots
#     - Green mean ± SE bars, jittered cells, biological-rep legend, bracket + stars
#     - Uses Mann–Whitney (Wilcoxon rank-sum) within each Method for stars
#     - Log plots summarize log(Results + 0.5), not log of the raw mean
#     - Within each Treatment, raw + log method plots with blue mean ± SE bars,
#       jittered cells and sample counts (no significance tests)
# ============================================================

# ----------------- PACKAGES (quiet) -----------------
suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(ggplot2)
  library(dplyr)
  library(emmeans)  # <-- NEW (post hoc comparisons)
})

options(warn = -1)

# ============================================================
# RUN SETTING
# ============================================================
input_dir <- '/Users/allisonpickle/Desktop/261014_Updated Skeletonization'
file_pattern <- "\\.csv$"

rename_cols_by_position <- NULL

# ============================================================
# FACTORS / LABELS
# ============================================================
baseline_treatment <- "Tri_Ctrl"
baseline_method <- "Skeletonize"

treatment_levels <- c("Tri_Ctrl", "Tri_5ugLPS")
method_levels <- c("Skeletonize", "Thin", "MedialAxis")

treatment_labels <- c(
  "Tri_Ctrl"   = "-\n-",
  "Tri_5ugLPS" = "-\n+"
)

row_labels <- c("5ug/mL LPS")

plot_title_base <- "Treatment × Method"
x_axis_label <- ""
y_axis_label_raw <- "Results (mean ± SE)"
y_axis_label_log <- "log(Results + 0.5) (mean ± SE)"

# ============================================================
# POST HOC SETTINGS (Ctrl vs 5ugLPS within each Method)
# ============================================================
posthoc_adjust <- "holm"        # adjust across methods
alpha <- 0.05
show_ns <- FALSE                # show "n.s." on plots?
use_adjusted_p <- TRUE          # use adjusted p for stars?

p_to_stars <- function(p) {
  if (is.na(p)) return("n.s.")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "n.s."
}

# ============================================================
# HELPERS
# ============================================================
make_aligned_tick_labels <- function(levels_vec, labels_named_vec) {
  out <- rep("", length(levels_vec))
  names(out) <- levels_vec
  common <- intersect(levels_vec, names(labels_named_vec))
  out[common] <- as.character(labels_named_vec[common])
  out
}

add_sig_bracket <- function(p, x1, x2, y, h, label, lw = 0.7) {
  p +
    ggplot2::annotate("segment", x = x1, xend = x1, y = y, yend = y + h, linewidth = lw) +
    ggplot2::annotate("segment", x = x1, xend = x2, y = y + h, yend = y + h, linewidth = lw) +
    ggplot2::annotate("segment", x = x2, xend = x2, y = y + h, yend = y, linewidth = lw) +
    ggplot2::annotate("text", x = (x1 + x2) / 2, y = y + h, label = label, vjust = -0.2, fontface = "bold")
}

# NEW: Mann–Whitney (Wilcoxon rank-sum) Ctrl vs 5ugLPS helper
wilcox_ctrl_vs_lps <- function(df, response_col = "Results") {
  if (!("Treatment" %in% names(df))) stop("wilcox_ctrl_vs_lps(): Treatment column missing.")
  if (!(response_col %in% names(df))) stop(paste0("wilcox_ctrl_vs_lps(): ", response_col, " column missing."))

  d <- df[, c("Treatment", response_col), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d$Treatment <- as.character(d$Treatment)

  x <- d[d$Treatment == "Tri_Ctrl", response_col]
  y <- d[d$Treatment == "Tri_5ugLPS", response_col]

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  if (length(x) < 1 || length(y) < 1) return(NA_real_)

  # exact=FALSE avoids issues with ties / larger N
  stats::wilcox.test(x, y, exact = FALSE)$p.value
}

# ============================================================
# FIT: adaptive random-effects structure for Treatment×Method
# ============================================================
fit_trt_method_model <- function(df, response = c("Results", "logResults")) {
  response <- match.arg(response)

  if (!("Biological_Rep" %in% names(df))) stop("Biological_Rep column missing.")
  if (!("Well" %in% names(df))) stop("Well column missing.")
  if (!("Treatment" %in% names(df))) stop("Treatment column missing.")
  if (!("Method" %in% names(df))) stop("Method column missing.")
  if (!(response %in% names(df))) stop(paste0(response, " column missing."))

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well           <- droplevels(factor(df$Well))
  df$Treatment      <- droplevels(factor(df$Treatment, levels = treatment_levels))
  df$Method         <- droplevels(factor(df$Method, levels = method_levels))

  has_img <- "Image" %in% names(df)
  if (has_img) df$Image <- droplevels(factor(df$Image))

  n_bio <- nlevels(droplevels(df$Biological_Rep))

  if (n_bio >= 2 && has_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Method + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well) + (1 | Well:Image)"
    ))
  } else if (n_bio >= 2) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Method + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
    ))
  } else if (has_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Method + ",
      "(1 | Well) + (1 | Well:Image)"
    ))
  } else {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Method + (1 | Well)"
    ))
  }

  suppressMessages(suppressWarnings(
    lmerTest::lmer(f, data = df)
  ))
}

# ============================================================
# LRTs for Treatment / Method / Interaction (ML fits)
# ============================================================
make_lrt_trt_method <- function(fit, response_name, df) {
  base_form  <- formula(fit)
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

  m0 <- suppressMessages(suppressWarnings(lme4::lmer(make_form("1"),                  data = df, REML = FALSE)))
  mT <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Treatment"),          data = df, REML = FALSE)))
  mM <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Method"),             data = df, REML = FALSE)))
  mA <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Treatment + Method"), data = df, REML = FALSE)))
  mF <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Treatment * Method"), data = df, REML = FALSE)))

  list(
    null = m0,
    trt = mT,
    meth = mM,
    additive = mA,
    full = mF,
    LRT_Treatment = tryCatch(anova(m0, mT), error = function(e) NULL),
    LRT_Method = tryCatch(anova(m0, mM), error = function(e) NULL),
    LRT_Interaction = tryCatch(anova(mA, mF), error = function(e) NULL)
  )
}

# ============================================================
# Single ICC
# ============================================================
compute_icc <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  var_resid <- attr(lme4::VarCorr(fit), "sc")^2
  var_random <- sum(vc$vcov, na.rm = TRUE)
  var_random / (var_random + var_resid)
}

# ============================================================
# Post hoc: Ctrl vs 5ugLPS within each Method (from fitted mixed model)
# ============================================================
posthoc_ctrl_vs_lps_within_method <- function(fit, adjust = "holm") {
  em <- emmeans::emmeans(
    fit,
    specs = ~ Treatment | Method,
    lmer.df = "kenward-roger"
  )

  ct <- emmeans::contrast(em, method = "pairwise", adjust = "none")
  df <- as.data.frame(ct)

  if (!("p.value" %in% names(df))) stop("emmeans output missing p.value column.")
  if (!("Method" %in% names(df))) stop("emmeans output missing Method column.")

  df$P <- df$p.value
  df$P_adj <- stats::p.adjust(df$P, method = adjust)

  df
}

# ============================================================
# Bar summaries + plots
# ============================================================
build_cell_summary <- function(data, response_col) {
  data |>
    dplyr::group_by(Treatment, Method) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      n    = sum(!is.na(.data[[response_col]])),
      .groups = "drop"
    )
}

plot_methods_within_treatment <- function(data, response_col, y_label, title_text) {
  summ <- build_cell_summary(data, response_col)
  summ$Method <- factor(as.character(summ$Method), levels = method_levels)
  data$Method <- factor(as.character(data$Method), levels = method_levels)
  tick_labels <- stats::setNames(
    paste0(as.character(summ$Method), "\n(n=", summ$n, ")"),
    as.character(summ$Method)
  )
  y_min <- min(0, data[[response_col]], summ$mean - summ$se, na.rm = TRUE)
  y_max <- max(0, data[[response_col]], summ$mean + summ$se, na.rm = TRUE)
  span <- max(y_max - y_min, 1e-9)

  ggplot2::ggplot(summ, ggplot2::aes(x = Method, y = mean)) +
    ggplot2::geom_col(fill = "#4C72B0", color = "#4D4D4D", alpha = 0.65, width = 0.8) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.18
    ) +
    ggplot2::geom_jitter(
      data = data,
      ggplot2::aes(x = Method, y = .data[[response_col]]),
      position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 1),
      shape = 21, fill = "black", color = "black", alpha = 0.65, size = 2
    ) +
    ggplot2::scale_x_discrete(labels = tick_labels) +
    ggplot2::scale_y_continuous(
      limits = c(y_min, y_max + 0.05 * span),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(title = title_text, x = "", y = y_label) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
}

plot_bar_trt_method <- function(summ_df, raw_data, response_col, y_label, title_text) {
  raw_data <- raw_data[!is.na(raw_data$Biological_Rep), , drop = FALSE]

  summ_df$Treatment  <- factor(as.character(summ_df$Treatment), levels = treatment_levels)
  summ_df$Method     <- factor(as.character(summ_df$Method), levels = method_levels)
  raw_data$Treatment <- factor(as.character(raw_data$Treatment), levels = treatment_levels)
  raw_data$Method    <- factor(as.character(raw_data$Method), levels = method_levels)

  aligned_treatment_labels <- make_aligned_tick_labels(treatment_levels, treatment_labels)
  aligned_treatment_labels <- stats::setNames(paste0("\n", aligned_treatment_labels), names(aligned_treatment_labels))

  combined_row_label <- ""
  if (exists("row_labels", inherits = TRUE) && length(row_labels) > 0) {
    combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  }

  p <- ggplot2::ggplot(summ_df, ggplot2::aes(x = Treatment, y = mean, fill = Method)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85), alpha = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.85)
    ) +
    ggplot2::geom_jitter(
      data = raw_data,
      ggplot2::aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep, color = Method),
      position = ggplot2::position_jitterdodge(jitter.width = 0.18, dodge.width = 0.85),
      alpha = 0.55,
      size = 1.2
    ) +
    ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::labs(
      title = title_text,
      x = x_axis_label,
      y = y_label,
      fill = "Method",
      color = "Method",
      shape = "Biological Rep"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = -0.2, size = 9, lineheight = 0.8),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(5, 5, 45, 55)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!identical(combined_row_label, "")) {
    p <- p + ggplot2::annotate(
      "text",
      x = 0.4, y = -Inf,
      label = combined_row_label,
      hjust = 1,
      vjust = 1.3,
      size = 3, color = "black",
      lineheight = 0.8
    )
  }

  p
}

# NOTE: This function is currently not used by run_one_csv().
# It was left as-is (but it has brace/loop issues in your current script).
# ============================================================
# DIAGNOSTICS (unchanged)
# ============================================================
save_diagnostics <- function(out_dir, raw_fit, log_fit, save_data = TRUE) {
  models <- list(RawFit = raw_fit, LogFit = log_fit)
  all_resid <- list()

  for (nm in names(models)) {
    fit <- models[[nm]]

    fitted_vals   <- as.numeric(stats::fitted(fit))
    pearson_resid <- as.numeric(stats::resid(fit, type = "pearson"))
    raw_resid     <- as.numeric(stats::resid(fit))

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
      )

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(nm, "_residuals.png")),
      plot = p_res,
      width = 7, height = 5, dpi = 200
    )

    qq <- stats::qqnorm(raw_resid, plot.it = FALSE)
    qq_df <- data.frame(
      TheoreticalQuantiles = qq$x,
      SampleQuantiles = qq$y
    )

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
      )

    ggplot2::ggsave(
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
    }
  }

  if (save_data) {
    utils::write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# SAVE MODEL SUMMARIES (UPDATED to include posthoc table)
# ============================================================
save_model_summaries <- function(
  out_dir,
  input_file,
  analysis_data,
  raw_fit,
  log_fit,
  icc_raw,
  icc_log,
  lrt_raw_list,
  lrt_log_list,
  posthoc_raw,
  posthoc_log
) {
  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")
  sink(out_file)

  cat("MIXED EFFECTS MODEL SUMMARIES (Treatment × Method)\n\n")
  cat("Input: ", input_file, "\n", sep = "")
  cat("Rows used: ", nrow(analysis_data), "\n", sep = "")
  cat("Biological reps: ", paste(sort(unique(analysis_data$Biological_Rep)), collapse = ", "), "\n\n", sep = "")

  cat("RAW model\n")
  print(formula(raw_fit))
  cat("\n")
  print(summary(raw_fit))
  cat("\nICC (raw): ", sprintf("%.4f", icc_raw), "\n\n", sep = "")

  cat("LRTs (raw; ML)\n")
  if (is.null(lrt_raw_list$LRT_Treatment)) cat("  Treatment: Not available\n") else { cat("  Treatment (null vs Treatment)\n"); print(lrt_raw_list$LRT_Treatment) }
  cat("\n")
  if (is.null(lrt_raw_list$LRT_Method)) cat("  Method: Not available\n") else { cat("  Method (null vs Method)\n"); print(lrt_raw_list$LRT_Method) }
  cat("\n")
  if (is.null(lrt_raw_list$LRT_Interaction)) cat("  Interaction: Not available\n") else { cat("  Interaction (additive vs interaction)\n"); print(lrt_raw_list$LRT_Interaction) }
  cat("\n\n")

  cat("POST HOC (raw): Ctrl vs 5ugLPS within each Method (", posthoc_adjust, " adjusted)\n", sep = "")
  print(posthoc_raw)
  cat("\n\n")

  cat("LOG model (log(Results + 0.5))\n")
  print(formula(log_fit))
  cat("\n")
  print(summary(log_fit))
  cat("\nICC (log): ", sprintf("%.4f", icc_log), "\n\n", sep = "")

  cat("LRTs (log; ML)\n")
  if (is.null(lrt_log_list$LRT_Treatment)) cat("  Treatment: Not available\n") else { cat("  Treatment (null vs Treatment)\n"); print(lrt_log_list$LRT_Treatment) }
  cat("\n")
  if (is.null(lrt_log_list$LRT_Method)) cat("  Method: Not available\n") else { cat("  Method (null vs Method)\n"); print(lrt_log_list$LRT_Method) }
  cat("\n")
  if (is.null(lrt_log_list$LRT_Interaction)) cat("  Interaction: Not available\n") else { cat("  Interaction (additive vs interaction)\n"); print(lrt_log_list$LRT_Interaction) }
  cat("\n\n")

  cat("POST HOC (log): Ctrl vs 5ugLPS within each Method (", posthoc_adjust, " adjusted)\n", sep = "")
  print(posthoc_log)
  cat("\n")

  sink()
}

# ============================================================
# MAIN PER-FILE RUNNER (UPDATED: Mann–Whitney asterisks on Plot B)
# ============================================================
run_one_csv <- function(input_file) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("PROCESSING FILE: ", input_file, "\n", sep = "")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

  process_data <- read.csv(input_file, stringsAsFactors = FALSE)

  if (!is.null(rename_cols_by_position)) {
    needed <- max(rename_cols_by_position)
    if (ncol(process_data) < needed) {
      stop(
        "CSV has ", ncol(process_data),
        " columns but rename_cols_by_position needs at least ", needed,
        " columns. File: ", basename(input_file)
      )
    }
    for (nm in names(rename_cols_by_position)) {
      colnames(process_data)[rename_cols_by_position[[nm]]] <- nm
    }
  }

  req <- c("Treatment", "Method", "Well", "Results", "Biological_Rep")
  miss <- setdiff(req, names(process_data))
  if (length(miss) > 0) stop("Missing required column(s): ", paste(miss, collapse = ", "))

  process_data$Results <- suppressWarnings(as.numeric(process_data$Results))
  process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)
  process_data$Method <- factor(process_data$Method, levels = method_levels)

  if (baseline_treatment %in% levels(process_data$Treatment)) {
    process_data$Treatment <- stats::relevel(process_data$Treatment, ref = baseline_treatment)
  }
  if (baseline_method %in% levels(process_data$Method)) {
    process_data$Method <- stats::relevel(process_data$Method, ref = baseline_method)
  }

  process_data$Well <- factor(process_data$Well)
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)

  has_image <- "Image" %in% names(process_data)
  if (has_image) process_data$Image <- factor(process_data$Image)

  needed <- c("Results", "Treatment", "Method", "Well", "Biological_Rep")
  if (has_image) needed <- c(needed, "Image")

  analysis_data <- process_data[stats::complete.cases(process_data[, needed, drop = FALSE]), , drop = FALSE]
  analysis_data <- analysis_data[!is.na(analysis_data$Treatment) & !is.na(analysis_data$Method), , drop = FALSE]
  if (nrow(analysis_data) == 0) stop("No rows left after filtering complete cases / Treatment/Method levels.")

  if (any(!is.finite(analysis_data$Results) | analysis_data$Results <= -0.5)) {
    stop("Results must be finite and greater than -0.5 for log(Results + 0.5).")
  }
  analysis_data$logResults <- log(analysis_data$Results + 0.5)

  analysis_data$Treatment <- droplevels(factor(analysis_data$Treatment, levels = treatment_levels))
  analysis_data$Method <- droplevels(factor(analysis_data$Method, levels = method_levels))
  analysis_data$Well <- droplevels(factor(analysis_data$Well))
  analysis_data$Biological_Rep <- droplevels(factor(analysis_data$Biological_Rep))
  if (has_image) analysis_data$Image <- droplevels(factor(analysis_data$Image))

  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Fit models
  raw_fit <- fit_trt_method_model(analysis_data, response = "Results")
  log_fit <- fit_trt_method_model(analysis_data, response = "logResults")

  # LRTs (ML)
  lrt_raw <- make_lrt_trt_method(raw_fit, "Results", analysis_data)
  lrt_log <- make_lrt_trt_method(log_fit, "logResults", analysis_data)

  if (!is.null(lrt_raw$LRT_Treatment))   capture.output(lrt_raw$LRT_Treatment, file = file.path(out_dir, "LRT_Raw_Treatment_vs_Null.txt"))
  if (!is.null(lrt_raw$LRT_Method))      capture.output(lrt_raw$LRT_Method, file = file.path(out_dir, "LRT_Raw_Method_vs_Null.txt"))
  if (!is.null(lrt_raw$LRT_Interaction)) capture.output(lrt_raw$LRT_Interaction, file = file.path(out_dir, "LRT_Raw_Interaction_vs_Additive.txt"))

  if (!is.null(lrt_log$LRT_Treatment))   capture.output(lrt_log$LRT_Treatment, file = file.path(out_dir, "LRT_Log_Treatment_vs_Null.txt"))
  if (!is.null(lrt_log$LRT_Method))      capture.output(lrt_log$LRT_Method, file = file.path(out_dir, "LRT_Log_Method_vs_Null.txt"))
  if (!is.null(lrt_log$LRT_Interaction)) capture.output(lrt_log$LRT_Interaction, file = file.path(out_dir, "LRT_Log_Interaction_vs_Additive.txt"))

  # ICC
  icc_raw <- compute_icc(raw_fit)
  icc_log <- compute_icc(log_fit)
  write.csv(
    data.frame(Model = c("RawFit", "LogFit"), ICC = c(icc_raw, icc_log)),
    file.path(out_dir, "ICC_Summary.csv"),
    row.names = FALSE
  )

  # POST HOC (requested; still saved to CSV/txt for reporting)
  posthoc_raw <- posthoc_ctrl_vs_lps_within_method(raw_fit, adjust = posthoc_adjust)
  posthoc_log <- posthoc_ctrl_vs_lps_within_method(log_fit, adjust = posthoc_adjust)

  write.csv(posthoc_raw, file.path(out_dir, "PostHoc_Raw_Ctrl_vs_5ugLPS_within_Method.csv"), row.names = FALSE)
  write.csv(posthoc_log, file.path(out_dir, "PostHoc_Log_Ctrl_vs_5ugLPS_within_Method.csv"), row.names = FALSE)

  # Per-method Ctrl vs 5ugLPS with matching raw/log styling.
  for (scale_name in c("Raw", "Log")) {
    response_col <- if (scale_name == "Raw") "Results" else "logResults"
    y_label <- if (scale_name == "Raw") y_axis_label_raw else y_axis_label_log
    for (trt in treatment_levels) {
      sub <- analysis_data[analysis_data$Treatment == trt, , drop = FALSE]
      if (nrow(sub) == 0) next
      treatment_title <- if (trt == "Tri_Ctrl") "Ctrl" else "5 \u00b5g/mL LPS"
      p_methods <- plot_methods_within_treatment(
        data = sub,
        response_col = response_col,
        y_label = if (scale_name == "Raw") "AUC" else "log(AUC + 0.5)",
        title_text = paste0("AUC by Skeleton Method (", treatment_title, ")")
      )
      ggplot2::ggsave(
        file.path(out_dir, paste0("AUC_by_Skeleton_Method_", trt, "_", scale_name, ".png")),
        p_methods, width = 10, height = 6, dpi = 300
      )
    }

    for (m in method_levels) {
      sub <- analysis_data[analysis_data$Method == m, , drop = FALSE]
      if (nrow(sub) == 0) next

      summ <- sub |>
        dplyr::group_by(Treatment) |>
        dplyr::summarise(
          mean = mean(.data[[response_col]], na.rm = TRUE),
          se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
          n    = sum(!is.na(.data[[response_col]])),
          .groups = "drop"
        )
      summ$Treatment <- factor(as.character(summ$Treatment), levels = treatment_levels)

      p <- ggplot2::ggplot(summ, ggplot2::aes(x = Treatment, y = mean)) +
        ggplot2::geom_col(fill = "#55A868", alpha = 0.70, width = 0.65) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.18) +
        ggplot2::geom_jitter(
          data = sub,
          ggplot2::aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep),
          position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 1),
          alpha = 0.65, size = 1.2, color = "black"
        ) +
        ggplot2::labs(
          title = paste0("Ctrl vs 5ugLPS (", m, ")"),
          x = "",
          y = y_label,
          shape = "Biological Rep"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          legend.position = "right"
        )

      # Keep negative log values and the full SE bars visible.
      vals <- sub[[response_col]]
      y_min <- min(0, vals, summ$mean - summ$se, na.rm = TRUE)
      y_max <- max(0, vals, summ$mean + summ$se, na.rm = TRUE)
      span <- max(y_max - y_min, 1e-9)
      y0 <- y_max + 0.06 * span
      h <- 0.05 * span
      y_top <- y0 + h * 1.6

      p <- p + ggplot2::scale_y_continuous(
        limits = c(y_min, y_top),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      )

      # Mann–Whitney stars
      p_mw <- wilcox_ctrl_vs_lps(sub, response_col = response_col)
      lab_mw <- p_to_stars(p_mw)
      is_sig_mw <- is.finite(p_mw) && (p_mw < alpha)

      if (is_sig_mw || show_ns) {
        p <- add_sig_bracket(p, x1 = 1, x2 = 2, y = y0, h = h, label = lab_mw)
      }

      ggplot2::ggsave(
        file.path(out_dir, paste0("Ctrl_vs_5ugLPS_", m, "_", scale_name, ".png")),
        p, width = 5.5, height = 4.3, dpi = 300
      )
    }
  }

  # Diagnostics + summaries
  save_diagnostics(out_dir, raw_fit, log_fit, save_data = TRUE)
  save_model_summaries(
    out_dir, input_file, analysis_data,
    raw_fit, log_fit,
    icc_raw, icc_log,
    lrt_raw, lrt_log,
    posthoc_raw, posthoc_log
  )

  cat("Done. Outputs saved to: ", out_dir, "\n", sep = "")
  invisible(list(input_file = input_file, out_dir = out_dir, n_rows = nrow(analysis_data)))
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
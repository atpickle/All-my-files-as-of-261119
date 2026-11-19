# ============================================================
# Mixed Effects (Sholl AUC): Treatment × RadiusStart with adaptive nested random effects
# Per-radius Ctrl vs LPS plots are saved on raw and natural-log scales.
# Log plots show mean ± SE of log(Results + 0.5), with Mann-Whitney stars.
# Adds separate raw/log R0 vs R9.81 figures for each treatment, using
# paired Wilcoxon signed-rank tests by Key and Holm adjustment across treatments.
# ============================================================

suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(ggplot2)
  library(dplyr)
  library(emmeans)
})

# ----------------------------
# RUN SETTINGS (UPDATE PATH)
# ----------------------------
input_dir <- '/Users/allisonpickle/Desktop/261014_Updated Skeletonization'
file_pattern <- "_RadiusStart_Long\\.csv$"

# If your CSVs do not already have these exact headers, set this mapping.
# rename_cols_by_position <- list(Treatment = 1, RadiusStart = 2, Well = 3, Results = 4, Biological_Rep = 5)
rename_cols_by_position <- NULL

# ----------------------------
# FACTORS / LABELS
# ----------------------------
baseline_treatment <- "Tri_Ctrl"
baseline_radius <- "R0"

treatment_levels <- c("Tri_Ctrl", "Tri_5ugLPS")
radius_levels <- c("R0", "R9.81")   # must match your CSV RadiusStart values

treatment_labels <- c(
  "Tri_Ctrl"   = "-\n-",
  "Tri_5ugLPS" = "-\n+"
)

plot_title_base <- "Treatment × Starting Radius"
x_axis_label <- ""
y_axis_label_raw <- "Results (mean ± SE)"
y_axis_label_log <- "log(Results + 0.5) (mean ± SE)"

# ----------------------------
# POST HOC SETTINGS
# ----------------------------
posthoc_adjust <- "holm"   # adjust across starting radii (2 tests)
alpha <- 0.05
show_ns <- FALSE

compare_paired_radii <- function(data, response_col) {
  if (!("Key" %in% names(data))) stop("Key column required for paired radius comparisons.")
  if (anyNA(data$Key) || any(!nzchar(trimws(as.character(data$Key))))) {
    stop("Every observation needs a nonempty Key for paired radius comparisons.")
  }
  if (anyDuplicated(data[, c("Treatment", "Key", "RadiusStart")])) {
    stop("Duplicate Treatment/Key/RadiusStart rows: cannot uniquely pair cells.")
  }
  rows <- lapply(treatment_levels, function(trt) {
    sub <- data[data$Treatment == trt, , drop = FALSE]
    left <- sub[sub$RadiusStart == radius_levels[1], c("Key", response_col), drop = FALSE]
    right <- sub[sub$RadiusStart == radius_levels[2], c("Key", response_col), drop = FALSE]
    paired <- merge(left, right, by = "Key", suffixes = c("_1", "_2"))
    if (nrow(paired) == 0) stop("No complete radius pairs for ", trt)
    if (nrow(paired) != nrow(left) || nrow(paired) != nrow(right)) {
      stop("Unmatched radius observations for ", trt, "; both radii are required for every Key.")
    }
    x <- paired[[paste0(response_col, "_1")]]
    y <- paired[[paste0(response_col, "_2")]]
    if (any(!is.finite(x) | !is.finite(y))) stop("Paired responses must be finite.")
    p <- if (all(x == y)) 1 else stats::wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value
    data.frame(Treatment = trt, Radius_1 = radius_levels[1], Radius_2 = radius_levels[2],
               Response = response_col, N_pairs = nrow(paired), P = p)
  })
  tab <- do.call(rbind, rows)
  tab$P_adj <- stats::p.adjust(tab$P, method = posthoc_adjust)
  tab$Stars <- vapply(tab$P_adj, p_to_stars, character(1))
  tab
}

plot_paired_radii <- function(data, response_col, y_label, tests, treatment_name) {
  data <- data[data$Treatment == treatment_name, , drop = FALSE]
  if (nrow(data) == 0) stop("No observations for ", treatment_name)
  test <- tests[tests$Treatment == treatment_name, , drop = FALSE]
  if (nrow(test) != 1) stop("Expected one paired-radius test for ", treatment_name)
  summ <- build_cell_summary_radius(data, response_col)
  vals <- data[[response_col]]
  lo <- min(0, vals, summ$mean - summ$se, na.rm = TRUE)
  hi <- max(0, vals, summ$mean + summ$se, na.rm = TRUE)
  span <- max(hi - lo, 1e-9)
  y0 <- hi + 0.06 * span
  h <- 0.05 * span
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = RadiusStart, y = mean)) +
    ggplot2::geom_col(fill = "#55A868", alpha = .70, width = .65) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = .18) +
    ggplot2::geom_jitter(
      data = data,
      ggplot2::aes(y = .data[[response_col]], shape = Biological_Rep),
      position = ggplot2::position_jitter(width = .18, height = 0, seed = 1),
      alpha = .65, size = 1.2, color = "black"
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(lo, y0 + h * 1.6),
      expand = ggplot2::expansion(mult = c(0, .02))
    ) +
    ggplot2::labs(
      title = paste0("Start radius comparison (", treatment_name, ")"),
      x = "", y = y_label, shape = "Biological Rep"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = .5, face = "bold"),
      legend.position = "right"
    )
  if (test$P_adj < alpha || show_ns) {
    p <- add_sig_bracket(p, x1 = 1, x2 = 2, y = y0, h = h, label = test$Stars)
  }
  p
}

# ----------------------------
# Standalone helpers
# ----------------------------
make_aligned_tick_labels <- function(levels_vec, labels_named_vec) {
  out <- stats::setNames(rep("", length(levels_vec)), levels_vec)
  common <- intersect(levels_vec, names(labels_named_vec))
  out[common] <- as.character(labels_named_vec[common])
  out
}

p_to_stars <- function(p) {
  if (!is.finite(p)) return("n.s.")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "n.s."
}

add_sig_bracket <- function(p, x1, x2, y, h, label, lw = 0.7) {
  p +
    ggplot2::annotate("segment", x = x1, xend = x1, y = y, yend = y + h, linewidth = lw) +
    ggplot2::annotate("segment", x = x1, xend = x2, y = y + h, yend = y + h, linewidth = lw) +
    ggplot2::annotate("segment", x = x2, xend = x2, y = y + h, yend = y, linewidth = lw) +
    ggplot2::annotate("text", x = (x1 + x2) / 2, y = y + h, label = label,
                      vjust = -0.2, fontface = "bold")
}

wilcox_ctrl_vs_lps <- function(df, response_col = "Results") {
  x <- df[df$Treatment == "Tri_Ctrl", response_col]
  y <- df[df$Treatment == "Tri_5ugLPS", response_col]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) {
    stop("Both treatments need finite observations for the Mann-Whitney comparison.")
  }
  stats::wilcox.test(x, y, exact = FALSE)$p.value
}

compute_icc <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  var_random <- sum(vc$vcov[vc$grp != "Residual"])
  var_random / (var_random + stats::sigma(fit)^2)
}

save_diagnostics <- function(out_dir, raw_fit, log_fit, save_data = TRUE) {
  models <- list(RawFit = raw_fit, LogFit = log_fit)
  all_resid <- list()
  for (nm in names(models)) {
    fit <- models[[nm]]
    res <- data.frame(
      Model = nm,
      FittedValues = as.numeric(stats::fitted(fit)),
      PearsonResiduals = as.numeric(stats::resid(fit, type = "pearson")),
      RawResiduals = as.numeric(stats::resid(fit))
    )
    p_res <- ggplot2::ggplot(res, ggplot2::aes(x = FittedValues, y = PearsonResiduals)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::geom_hline(yintercept = 0, color = "red", linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = paste0(nm, ": Residuals vs Fitted"),
                    x = "Fitted Values", y = "Pearson Residuals")
    p_qq <- ggplot2::ggplot(res, ggplot2::aes(sample = RawResiduals)) +
      ggplot2::stat_qq(alpha = 0.7, size = 1.2) +
      ggplot2::stat_qq_line(color = "red", linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = paste0(nm, ": QQ Plot"),
                    x = "Theoretical Quantiles", y = "Sample Quantiles")
    ggplot2::ggsave(file.path(out_dir, paste0(nm, "_residuals.png")),
                    p_res, width = 7, height = 5, dpi = 200)
    ggplot2::ggsave(file.path(out_dir, paste0(nm, "_qq.png")),
                    p_qq, width = 7, height = 5, dpi = 200)
    all_resid[[nm]] <- res
  }
  if (save_data) {
    utils::write.csv(do.call(rbind, all_resid),
                      file.path(out_dir, "All_Calculated_Residuals.csv"), row.names = FALSE)
  }
}

save_model_summaries <- function(out_dir, input_file, analysis_data,
                                  raw_fit, log_fit, icc_raw, icc_log,
                                  lrt_raw, lrt_log, posthoc_raw, posthoc_log) {
  con <- file(file.path(out_dir, "Mixed_Model_Summaries.txt"), open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "MIXED EFFECTS MODEL SUMMARIES (Treatment x Starting Radius)",
    paste0("Input: ", input_file),
    paste0("Rows used: ", nrow(analysis_data)),
    paste0("Starting radii: ", paste(radius_levels, collapse = ", ")),
    "Log transform: natural log of Results + 0.5",
    "Post hoc: Ctrl vs 5ugLPS within each starting radius",
    paste0("Post hoc adjustment across radii: ", posthoc_adjust),
    "Plot stars: unadjusted Mann-Whitney tests, not mixed-model post hoc tests.",
    ""
  ), con)
  fits <- list(Raw = raw_fit, Log = log_fit)
  iccs <- c(Raw = icc_raw, Log = icc_log)
  lrts <- list(Raw = lrt_raw, Log = lrt_log)
  posthocs <- list(Raw = posthoc_raw, Log = posthoc_log)
  for (nm in names(fits)) {
    writeLines(paste0(nm, " model"), con)
    capture.output(print(formula(fits[[nm]])), print(summary(fits[[nm]])),
                   file = con, append = TRUE)
    writeLines(paste0("ICC: ", sprintf("%.4f", iccs[[nm]])), con)
    for (test in c("LRT_Treatment", "LRT_RadiusStart", "LRT_Interaction")) {
      writeLines(test, con)
      capture.output(print(lrts[[nm]][[test]]), file = con, append = TRUE)
    }
    writeLines("Mixed-model post hoc comparisons", con)
    capture.output(print(posthocs[[nm]]), file = con, append = TRUE)
    writeLines("", con)
  }
}

# ----------------------------
# FIT: adaptive random-effects structure for Treatment×RadiusStart
# ----------------------------
fit_trt_radius_model <- function(df, response = c("Results", "logResults")) {
  response <- match.arg(response)

  if (!("Biological_Rep" %in% names(df))) stop("Biological_Rep column missing.")
  if (!("Well" %in% names(df))) stop("Well column missing.")
  if (!("Treatment" %in% names(df))) stop("Treatment column missing.")
  if (!("RadiusStart" %in% names(df))) stop("RadiusStart column missing.")
  if (!(response %in% names(df))) stop(paste0(response, " column missing."))

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well           <- droplevels(factor(df$Well))
  df$Treatment      <- droplevels(factor(df$Treatment, levels = treatment_levels))
  df$RadiusStart    <- droplevels(factor(df$RadiusStart, levels = radius_levels))

  has_img <- "Image" %in% names(df)
  if (has_img) df$Image <- droplevels(factor(df$Image))

  n_bio <- nlevels(droplevels(df$Biological_Rep))

  if (n_bio >= 2 && has_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * RadiusStart + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well) + (1 | Well:Image)"
    ))
  } else if (n_bio >= 2) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * RadiusStart + ",
      "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
    ))
  } else if (has_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * RadiusStart + ",
      "(1 | Well) + (1 | Well:Image)"
    ))
  } else {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * RadiusStart + (1 | Well)"
    ))
  }

  suppressMessages(suppressWarnings(
    lmerTest::lmer(f, data = df)
  ))
}

# ----------------------------
# LRTs for Treatment / RadiusStart / Interaction (ML fits)
# ----------------------------
make_lrt_trt_radius <- function(fit, response_name, df) {
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

  m0 <- suppressMessages(suppressWarnings(lme4::lmer(make_form("1"),                       data = df, REML = FALSE)))
  mT <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Treatment"),               data = df, REML = FALSE)))
  mS <- suppressMessages(suppressWarnings(lme4::lmer(make_form("RadiusStart"),                data = df, REML = FALSE)))
  mA <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Treatment + RadiusStart"),    data = df, REML = FALSE)))
  mF <- suppressMessages(suppressWarnings(lme4::lmer(make_form("Treatment * RadiusStart"),    data = df, REML = FALSE)))

  list(
    null = m0,
    trt = mT,
    radius = mS,
    additive = mA,
    full = mF,
    LRT_Treatment = anova(m0, mT),
    LRT_RadiusStart = anova(m0, mS),
    LRT_Interaction = anova(mA, mF)
  )
}

# ----------------------------
# Post hoc: Ctrl vs 5ugLPS within each RadiusStart (emmeans)
# ----------------------------
posthoc_ctrl_vs_lps_within_radius <- function(fit, adjust = "holm") {
  em <- emmeans::emmeans(
    fit,
    specs = ~ Treatment | RadiusStart,
    lmer.df = "kenward-roger"
  )

  ct <- emmeans::contrast(em, method = "pairwise", adjust = "none")
  df <- as.data.frame(ct)

  if (!("p.value" %in% names(df))) stop("emmeans output missing p.value column.")
  if (!("RadiusStart" %in% names(df))) stop("emmeans output missing RadiusStart column.")

  df$P <- df$p.value
  df$P_adj <- stats::p.adjust(df$P, method = adjust)
  df
}

# ----------------------------
# Summaries + Plot A
# ----------------------------
build_cell_summary_radius <- function(data, response_col) {
  data |>
    dplyr::group_by(Treatment, RadiusStart) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      n    = sum(!is.na(.data[[response_col]])),
      .groups = "drop"
    )
}

plot_bar_trt_radius <- function(summ_df, raw_data, response_col, y_label, title_text) {
  raw_data <- raw_data[!is.na(raw_data$Biological_Rep), , drop = FALSE]

  summ_df$Treatment <- factor(as.character(summ_df$Treatment), levels = treatment_levels)
  summ_df$RadiusStart <- factor(as.character(summ_df$RadiusStart), levels = radius_levels)
  raw_data$Treatment <- factor(as.character(raw_data$Treatment), levels = treatment_levels)
  raw_data$RadiusStart <- factor(as.character(raw_data$RadiusStart), levels = radius_levels)

  aligned_treatment_labels <- make_aligned_tick_labels(treatment_levels, treatment_labels)
  aligned_treatment_labels <- stats::setNames(paste0("\n", aligned_treatment_labels), names(aligned_treatment_labels))

  ggplot2::ggplot(summ_df, ggplot2::aes(x = Treatment, y = mean, fill = RadiusStart)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85), alpha = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.2,
      position = ggplot2::position_dodge(width = 0.85)
    ) +
    ggplot2::geom_jitter(
      data = raw_data,
      ggplot2::aes(x = Treatment, y = .data[[response_col]], shape = Biological_Rep, color = RadiusStart),
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
      fill = "Starting Radius",
      color = "Starting Radius",
      shape = "Biological Rep"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = -0.2, size = 9, lineheight = 0.8),
      axis.title.x = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

# ----------------------------
# MAIN (starting-radius version)
# ----------------------------
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

  req <- c("Key", "Treatment", "RadiusStart", "Well", "Results", "Biological_Rep")
  miss <- setdiff(req, names(process_data))
  if (length(miss) > 0) stop("Missing required column(s): ", paste(miss, collapse = ", "))

  process_data$Results <- suppressWarnings(as.numeric(process_data$Results))
  process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)
  unknown_radii <- setdiff(unique(stats::na.omit(process_data$RadiusStart)), radius_levels)
  if (length(unknown_radii) > 0) {
    stop("Unexpected RadiusStart value(s): ", paste(unknown_radii, collapse = ", "),
         ". Expected: ", paste(radius_levels, collapse = ", "))
  }
  process_data$RadiusStart <- factor(process_data$RadiusStart, levels = radius_levels)

  if (baseline_treatment %in% levels(process_data$Treatment)) {
    process_data$Treatment <- stats::relevel(process_data$Treatment, ref = baseline_treatment)
  }
  if (baseline_radius %in% levels(process_data$RadiusStart)) {
    process_data$RadiusStart <- stats::relevel(process_data$RadiusStart, ref = baseline_radius)
  }

  process_data$Well <- factor(process_data$Well)
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)

  has_image <- "Image" %in% names(process_data)
  if (has_image) process_data$Image <- factor(process_data$Image)

  needed <- c("Results", "Treatment", "RadiusStart", "Well", "Biological_Rep")
  if (has_image) needed <- c(needed, "Image")

  analysis_data <- process_data[stats::complete.cases(process_data[, needed, drop = FALSE]), , drop = FALSE]
  analysis_data <- analysis_data[!is.na(analysis_data$Treatment) & !is.na(analysis_data$RadiusStart), , drop = FALSE]
  if (nrow(analysis_data) == 0) stop("No rows left after filtering complete cases / factor levels.")

  if (any(!is.finite(analysis_data$Results) | analysis_data$Results <= -0.5)) {
    stop("Results must be finite and greater than -0.5 for log(Results + 0.5).")
  }
  analysis_data$logResults <- log(analysis_data$Results + 0.5)

  analysis_data$Treatment <- droplevels(factor(analysis_data$Treatment, levels = treatment_levels))
  analysis_data$RadiusStart <- droplevels(factor(analysis_data$RadiusStart, levels = radius_levels))
  analysis_data$Well <- droplevels(factor(analysis_data$Well))
  analysis_data$Biological_Rep <- droplevels(factor(analysis_data$Biological_Rep))
  if (has_image) analysis_data$Image <- droplevels(factor(analysis_data$Image))

  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (scale_name in c("Raw", "Log")) {
    response_col <- if (scale_name == "Raw") "Results" else "logResults"
    tests <- compare_paired_radii(analysis_data, response_col)
    utils::write.csv(
      tests, file.path(out_dir, paste0("PairedWilcoxon_Radii_within_Treatment_", scale_name, ".csv")),
      row.names = FALSE
    )
    for (trt in treatment_levels) {
      p_radii <- plot_paired_radii(
        analysis_data, response_col,
        if (scale_name == "Raw") y_axis_label_raw else y_axis_label_log, tests, trt
      )
      ggplot2::ggsave(
        file.path(out_dir, paste0("Radius_Comparison_", trt, "_", scale_name, ".png")),
        p_radii, width = 5.5, height = 4.3, dpi = 300
      )
    }
  }

  # Fit models
  raw_fit <- fit_trt_radius_model(analysis_data, response = "Results")
  log_fit <- fit_trt_radius_model(analysis_data, response = "logResults")

  # LRTs (ML)
  lrt_raw <- make_lrt_trt_radius(raw_fit, "Results", analysis_data)
  lrt_log <- make_lrt_trt_radius(log_fit, "logResults", analysis_data)

  if (!is.null(lrt_raw$LRT_Treatment))   capture.output(lrt_raw$LRT_Treatment, file = file.path(out_dir, "LRT_Raw_Treatment_vs_Null.txt"))
  capture.output(lrt_raw$LRT_RadiusStart, file = file.path(out_dir, "LRT_Raw_RadiusStart_vs_Null.txt"))
  if (!is.null(lrt_raw$LRT_Interaction)) capture.output(lrt_raw$LRT_Interaction, file = file.path(out_dir, "LRT_Raw_Interaction_vs_Additive.txt"))

  if (!is.null(lrt_log$LRT_Treatment))   capture.output(lrt_log$LRT_Treatment, file = file.path(out_dir, "LRT_Log_Treatment_vs_Null.txt"))
  capture.output(lrt_log$LRT_RadiusStart, file = file.path(out_dir, "LRT_Log_RadiusStart_vs_Null.txt"))
  if (!is.null(lrt_log$LRT_Interaction)) capture.output(lrt_log$LRT_Interaction, file = file.path(out_dir, "LRT_Log_Interaction_vs_Additive.txt"))

  # ICC
  icc_raw <- compute_icc(raw_fit)
  icc_log <- compute_icc(log_fit)
  write.csv(
    data.frame(Model = c("RawFit", "LogFit"), ICC = c(icc_raw, icc_log)),
    file.path(out_dir, "ICC_Summary.csv"),
    row.names = FALSE
  )

  # Post hoc (emmeans)
  posthoc_raw <- posthoc_ctrl_vs_lps_within_radius(raw_fit, adjust = posthoc_adjust)
  posthoc_log <- posthoc_ctrl_vs_lps_within_radius(log_fit, adjust = posthoc_adjust)
  write.csv(posthoc_raw, file.path(out_dir, "PostHoc_Raw_Ctrl_vs_5ugLPS_within_RadiusStart.csv"), row.names = FALSE)
  write.csv(posthoc_log, file.path(out_dir, "PostHoc_Log_Ctrl_vs_5ugLPS_within_RadiusStart.csv"), row.names = FALSE)

  # Plot A: Treatment×RadiusStart
  summ_raw <- build_cell_summary_radius(analysis_data, "Results")
  p_raw <- plot_bar_trt_radius(summ_raw, analysis_data, "Results", y_axis_label_raw, paste0(plot_title_base, " (raw)"))
  ggplot2::ggsave(file.path(out_dir, "Barplot_Raw_Treatment_by_RadiusStart.png"), p_raw, width = 8.5, height = 4.5, dpi = 300)

  summ_log <- build_cell_summary_radius(analysis_data, "logResults")
  p_log <- plot_bar_trt_radius(summ_log, analysis_data, "logResults", y_axis_label_log, paste0(plot_title_base, " (log)"))
  ggplot2::ggsave(file.path(out_dir, "Barplot_Log_Treatment_by_RadiusStart.png"), p_log, width = 8.5, height = 4.5, dpi = 300)

  # Plot B: matching raw/log per-radius Ctrl vs 5ugLPS (Mann–Whitney)
  for (s in radius_levels) {
    sub <- analysis_data[analysis_data$RadiusStart == s, , drop = FALSE]
    if (nrow(sub) == 0) next

    for (scale_name in c("Raw", "Log")) {
    response_col <- if (scale_name == "Raw") "Results" else "logResults"
    y_label <- if (scale_name == "Raw") y_axis_label_raw else y_axis_label_log
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
        title = paste0("Ctrl vs 5ugLPS (", s, ")"),
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

    p_mw <- wilcox_ctrl_vs_lps(sub, response_col = response_col)
    lab_mw <- p_to_stars(p_mw)
    is_sig_mw <- is.finite(p_mw) && (p_mw < alpha)

    if (is_sig_mw || show_ns) {
      p <- add_sig_bracket(p, x1 = 1, x2 = 2, y = y0, h = h, label = lab_mw)
    }

    ggplot2::ggsave(
      file.path(out_dir, paste0("Ctrl_vs_5ugLPS_", s, "_", scale_name, ".png")),
      p, width = 5.5, height = 4.3, dpi = 300
    )
    }
  }

  save_diagnostics(out_dir, raw_fit, log_fit, save_data = TRUE)
  save_model_summaries(out_dir, input_file, analysis_data, raw_fit, log_fit, icc_raw, icc_log, lrt_raw, lrt_log, posthoc_raw, posthoc_log)

  cat("Done. Outputs saved to: ", out_dir, "\n", sep = "")
  invisible(list(input_file = input_file, out_dir = out_dir, n_rows = nrow(analysis_data)))
}

# ============================================================
# RUN ALL CSVs
# ============================================================
csv_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", input_dir)

cat("Found ", length(csv_files), " file(s) in ", input_dir, "\n", sep = "")

results <- lapply(csv_files, function(f) {
  tryCatch(run_one_csv(f), error = function(e) { message("FAILED: ", basename(f), " -> ", conditionMessage(e)); NULL })
})

# ...existing code...
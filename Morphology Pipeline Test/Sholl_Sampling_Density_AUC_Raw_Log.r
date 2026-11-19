# Log Sholl AUC: three paired sampling-density comparisons within each treatment.
# Input CSV: Key, Treatment, StepSize, Results (RAW AUC).
# Additional columns, including Biological_Rep, are ignored.
# Match the same cell Key across S2.5, S5 and Continuous; require complete pairs.
# Two-sided paired Wilcoxon signed-rank tests; Holm across three pairs per treatment.
# Exploratory within one biological replicate: cell pairing does not account for
# additional dependence between different cells from the same image or well.
# Outputs: two log plots, paired-test CSV/report, and raw/log analysis data.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

input_file <- "/Users/allisonpickle/Desktop/Sampling Density/Compiled Step Sizes_Sheet1_SamplingDensity_AUC_Long.csv"
output_dir <- NULL
treatment_levels <- c("Tri_Ctrl", "Tri_5ugLPS")
density_levels <- c("S2.5", "S5", "Continuous")
density_labels <- c(Continuous = "Continuous", S2.5 = "2.5 \u00b5m step", S5 = "5 \u00b5m step")
alpha <- 0.05
posthoc_adjust <- "holm"
show_ns <- FALSE

p_to_stars <- function(p) {
  if (!is.finite(p)) stop("Cannot label a nonfinite p-value.")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "n.s."
}

read_sampling_auc <- function(input_file) {
  if (!file.exists(input_file)) stop("CSV not found: ", input_file)
  dat <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE,
                         colClasses = "character")
  names(dat) <- trimws(names(dat))
  if (anyDuplicated(names(dat))) stop("Duplicate CSV column headers.")
  required <- c("Key", "Treatment", "StepSize", "Results")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing required CSV column(s): ", paste(missing, collapse = ", "))
  dat <- dat[, required]
  for (nm in required) dat[[nm]] <- trimws(dat[[nm]])
  if (!nrow(dat) || anyNA(dat) || any(!nzchar(as.matrix(dat)))) {
    stop("Input must contain complete Key, Treatment, StepSize and raw AUC observations.")
  }
  if (!all(dat$Treatment %in% treatment_levels)) stop("Unexpected Treatment; expected Tri_Ctrl or Tri_5ugLPS.")
  if (!all(dat$StepSize %in% density_levels)) stop("Unexpected StepSize; expected S2.5, S5 or Continuous.")
  dat$Results <- as.numeric(dat$Results)
  if (any(!is.finite(dat$Results) | dat$Results < 0)) stop("Raw AUC must be finite and nonnegative.")
  if (anyDuplicated(dat[, c("Treatment", "Key", "StepSize")])) {
    stop("Duplicate Treatment/Key/StepSize rows; cannot uniquely pair cells.")
  }
  for (trt in treatment_levels) {
    sub <- dat[dat$Treatment == trt, ]
    if (!nrow(sub)) stop("No observations for ", trt)
    for (density in density_levels) {
      if (!setequal(unique(sub$Key), sub$Key[sub$StepSize == density])) {
        stop("Unmatched cells for ", trt, " at ", density, "; every Key needs all three densities.")
      }
    }
  }
  dat$Treatment <- factor(dat$Treatment, levels = treatment_levels)
  dat$StepSize <- factor(dat$StepSize, levels = density_levels)
  dat$CellID <- interaction(dat$Treatment, dat$Key, drop = TRUE)
  dat$log_AUC <- log(dat$Results + 0.5)
  dat
}

paired_density_tests <- function(data, treatment_name) {
  sub <- data[data$Treatment == treatment_name, ]
  pairs <- utils::combn(density_levels, 2, simplify = FALSE)
  tab <- do.call(rbind, lapply(pairs, function(pair) {
    left <- sub[sub$StepSize == pair[1], c("CellID", "log_AUC")]
    right <- sub[sub$StepSize == pair[2], c("CellID", "log_AUC")]
    if (anyDuplicated(left$CellID) || anyDuplicated(right$CellID)) stop("Duplicate paired CellID.")
    matched <- merge(left, right, by = "CellID", suffixes = c("_1", "_2"))
    if (!nrow(matched) || nrow(matched) != nrow(left) || nrow(matched) != nrow(right)) {
      stop("Incomplete cell matching for ", treatment_name, ": ", paste(pair, collapse = " vs "))
    }
    differences <- matched$log_AUC_1 - matched$log_AUC_2
    test <- if (all(differences == 0)) NULL else stats::wilcox.test(
      matched$log_AUC_1, matched$log_AUC_2, paired = TRUE,
      alternative = "two.sided", exact = FALSE
    )
    p <- if (is.null(test)) 1 else test$p.value
    if (!is.finite(p)) stop("Nonfinite paired-test p-value for ", treatment_name)
    data.frame(Treatment = treatment_name, Density_1 = pair[1], Density_2 = pair[2],
               N_pairs = nrow(matched), N_nonzero_differences = sum(differences != 0),
               Median_paired_log_difference = median(differences),
               W = if (is.null(test)) 0 else unname(test$statistic), P = p)
  }))
  tab$P_adj <- stats::p.adjust(tab$P, method = posthoc_adjust)
  tab$Stars_adj <- vapply(tab$P_adj, p_to_stars, character(1))
  tab
}

plot_sampling_auc <- function(data, treatment_name, response_col, tests, y_label) {
  sub <- data[data$Treatment == treatment_name, ]
  summ <- sub |>
    dplyr::group_by(StepSize) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]]),
      se = sd(.data[[response_col]]) / sqrt(dplyr::n()),
      n = dplyr::n(), .groups = "drop"
    )
  lo <- min(0, sub[[response_col]], summ$mean - summ$se, na.rm = TRUE)
  hi <- max(0, sub[[response_col]], summ$mean + summ$se, na.rm = TRUE)
  span <- max(hi - lo, 1e-9)
  tests$P_use <- tests$P_adj
  tests$x1 <- match(tests$Density_1, density_levels)
  tests$x2 <- match(tests$Density_2, density_levels)
  tests <- tests[tests$P_use < alpha | show_ns, ]
  tests <- tests[order(tests$x2 - tests$x1, tests$x1), ]
  y0 <- hi + 0.06 * span
  h <- 0.05 * span
  step <- 0.14 * span
  top <- y0 + max(nrow(tests) - 1, 0) * step + 2 * h
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = StepSize, y = mean)) +
    ggplot2::geom_col(fill = "#4C72B0", alpha = 0.75, width = 0.65) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.18) +
    ggplot2::geom_jitter(
      data = sub, ggplot2::aes(y = .data[[response_col]]),
      position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
      alpha = 0.65, size = 1.2, color = "black"
    ) +
    ggplot2::scale_x_discrete(labels = density_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(lo, top),
                               expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(title = paste0("Sampling density comparison (", treatment_name, ")"),
                  x = "", y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  for (i in seq_len(nrow(tests))) {
    y <- y0 + (i - 1) * step
    p <- p +
      ggplot2::annotate("segment", x = tests$x1[i], xend = tests$x1[i], y = y, yend = y + h, linewidth = 0.7) +
      ggplot2::annotate("segment", x = tests$x1[i], xend = tests$x2[i], y = y + h, yend = y + h, linewidth = 0.7) +
      ggplot2::annotate("segment", x = tests$x2[i], xend = tests$x2[i], y = y + h, yend = y, linewidth = 0.7) +
      ggplot2::annotate("text", x = (tests$x1[i] + tests$x2[i]) / 2, y = y + h,
                        label = p_to_stars(tests$P_use[i]), vjust = -0.2, fontface = "bold")
  }
  p
}

run_sampling_density_auc <- function(input_file, output_dir = NULL) {
  dat <- read_sampling_auc(input_file)
  out_dir <- if (!is.null(output_dir)) output_dir else file.path(dirname(input_file), paste0(
    "Sampling_Density_AUC_Paired_Log_Outputs_", tools::file_path_sans_ext(basename(input_file))))
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory: ", out_dir)
  tables <- lapply(treatment_levels, function(trt) paired_density_tests(dat, trt))
  all_tests <- do.call(rbind, tables)
  utils::write.csv(dat, file.path(out_dir, "Analysis_Data_Raw_Log_AUC.csv"), row.names = FALSE)
  utils::write.csv(all_tests, file.path(out_dir, "PairedWilcoxon_Sampling_Density_Log.csv"), row.names = FALSE)
  writeLines(c(
    "Sampling density comparisons: natural log(AUC + 0.5), within each treatment.",
    "Test: two-sided paired Wilcoxon signed-rank by CellID (Treatment + Key), exact=FALSE.",
    "Holm correction across three comparisons per treatment; stars use adjusted p-values.",
    "All-zero paired differences are assigned p=1; other tests use a normal approximation with continuity correction.",
    "The signed-rank test assumes independent cell pairs and a symmetric paired-difference distribution for a location interpretation.",
    "Cell pairing does not adjust additional clustering among cells from the same image or well.",
    "Exploratory parameter evaluation within one biological replicate; not a biological-replicate treatment test.",
    "Bars show mean +/- SE of individually transformed observations; the test is not specifically a test of means.",
    capture.output(print(all_tests, row.names = FALSE))
  ), file.path(out_dir, "Stat_Results_Sampling_Density_Log.txt"))
  file_tags <- c(Tri_Ctrl = "Ctrl", Tri_5ugLPS = "LPS5ug")
  for (i in seq_along(treatment_levels)) {
    trt <- treatment_levels[i]
    p <- plot_sampling_auc(dat, trt, "log_AUC", tables[[i]], "log(AUC + 0.5) (mean \u00b1 SE)")
    ggplot2::ggsave(file.path(out_dir, paste0(file_tags[[trt]], "_Sampling_Density_AUC_Log.png")),
                    p, width = 7, height = 4.8, dpi = 300)
  }
  message("Outputs saved to: ", out_dir)
  invisible(list(out_dir = out_dir, analysis_data = dat, tests = all_tests))
}

run_sampling_density_auc(input_file, output_dir)

# ============================================================
# Standalone: StepSize compared within each Treatment (Ctrl-only and 5ugLPS-only)
# - For each CSV:
#   - makes 4 plots (raw and natural-log versions):
#       1) Ctrl: S2.5 vs S5
#       2) 5ugLPS: S2.5 vs S5
#   - Log plots summarize log(Results + 0.5), with mean ± SE on that scale.
#   - Saves separate raw/log Mann-Whitney tables and reports.
#   - Mann–Whitney (Wilcoxon rank-sum) for the one step-size pair
#   - (Optional) no p-adjust needed (only 1 test per treatment),
#     but kept as setting for consistency.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# -------------------- SETTINGS --------------------
input_dir <- '/Users/allisonpickle/Desktop/261014_Updated Skeletonization'
file_pattern <- "\\.csv$"

treatment_levels <- c("Tri_Ctrl", "Tri_5ugLPS")
stepsize_levels <- c("S2.5", "S5")  # must match your CSV values in StepSize

alpha <- 0.05
posthoc_adjust <- "holm"   # no real effect with 1 test, but OK to keep
show_ns <- FALSE
use_adjusted_p <- TRUE

y_axis_label_raw <- "Results (mean ± SE)"
y_axis_label_log <- "log(Results + 0.5) (mean ± SE)"

# -------------------- HELPERS --------------------
p_to_stars <- function(p) {
  if (is.na(p)) return("n.s.")
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
    ggplot2::annotate("text", x = (x1 + x2) / 2, y = y + h, label = label, vjust = -0.2, fontface = "bold")
}

wilcox_stepsize_within_treatment <- function(df,
                                            treatment_name,
                                            response_col = "Results",
                                            stepsize_levels_use = stepsize_levels,
                                            adjust = "holm") {
  sub_t <- df[df$Treatment == treatment_name, , drop = FALSE]

  s1 <- stepsize_levels_use[1]
  s2 <- stepsize_levels_use[2]

  x <- sub_t[sub_t$StepSize == s1, response_col]
  y <- sub_t[sub_t$StepSize == s2, response_col]

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  p <- NA_real_
  if (length(x) >= 1 && length(y) >= 1) {
    p <- stats::wilcox.test(x, y, exact = FALSE)$p.value
  }

  res <- data.frame(
    Treatment = treatment_name,
    StepSize_1 = s1,
    StepSize_2 = s2,
    n1 = length(x),
    n2 = length(y),
    P = p,
    stringsAsFactors = FALSE
  )

  res$P_adj <- stats::p.adjust(res$P, method = adjust)
  res$Stars_raw <- p_to_stars(res$P)
  res$Stars_adj <- p_to_stars(res$P_adj)
  res
}

plot_stepsize_within_treatment <- function(df,
                                          treatment_name,
                                          response_col = "Results",
                                          y_label = y_axis_label_raw,
                                          adjust = "holm",
                                          alpha = 0.05,
                                          show_ns = FALSE,
                                          use_adjusted_p = TRUE) {
  sub <- df[df$Treatment == treatment_name, , drop = FALSE]
  if (nrow(sub) == 0) stop("No rows for Treatment = ", treatment_name)

  sub$StepSize <- factor(as.character(sub$StepSize), levels = stepsize_levels)

  summ <- sub |>
    dplyr::group_by(StepSize) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      n    = sum(!is.na(.data[[response_col]])),
      .groups = "drop"
    )
  summ$StepSize <- factor(as.character(summ$StepSize), levels = stepsize_levels)

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = StepSize, y = mean)) +
    ggplot2::geom_col(fill = "#4C72B0", alpha = 0.75, width = 0.65) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0.18) +
    ggplot2::geom_jitter(
      data = sub,
      ggplot2::aes(x = StepSize, y = .data[[response_col]], shape = Biological_Rep),
      position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
      alpha = 0.65, size = 1.2, color = "black"
    ) +
    ggplot2::labs(
      title = paste0("Step size comparison (", treatment_name, ")"),
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

  mw_tab <- wilcox_stepsize_within_treatment(
    df = df,
    treatment_name = treatment_name,
    response_col = response_col,
    stepsize_levels_use = stepsize_levels,
    adjust = adjust
  )

  base_y0 <- y_max + 0.06 * span
  h <- 0.05 * span
  y_top <- base_y0 + h * 2.0

  p <- p + ggplot2::scale_y_continuous(
    limits = c(y_min, y_top),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  )

  p_use <- if (use_adjusted_p) mw_tab$P_adj[1] else mw_tab$P[1]
  lab <- p_to_stars(p_use)
  is_sig <- is.finite(p_use) && (p_use < alpha)

  if (is_sig || show_ns) {
    p <- add_sig_bracket(p, x1 = 1, x2 = 2, y = base_y0, h = h, label = lab)
  }

  list(plot = p, table = mw_tab)
}

write_mw_stepsize_report <- function(mw_df,
                                    out_file,
                                    use_adjusted_p = TRUE,
                                    adjust_name = "holm",
                                    alpha = 0.05) {
  if (is.null(mw_df) || nrow(mw_df) == 0) {
    writeLines("No Mann–Whitney results to report.", con = out_file)
    return(invisible(NULL))
  }

  pick_p <- if (use_adjusted_p) "P_adj" else "P"
  pick_stars <- if (use_adjusted_p) "Stars_adj" else "Stars_raw"

  con <- file(out_file, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("STEP SIZE COMPARISONS WITHIN EACH TREATMENT", con)
  writeLines("Test: Mann–Whitney (Wilcoxon rank-sum), two-sided", con)
  writeLines(paste0("Multiple-comparison adjustment within treatment: ", adjust_name), con)
  writeLines(paste0("Stars computed from ", pick_p, " (alpha = ", alpha, ")"), con)
  writeLines("", con)

  for (trt in unique(mw_df$Treatment)) {
    sub <- mw_df[mw_df$Treatment == trt, , drop = FALSE]
    writeLines(paste0("Treatment: ", trt), con)

    sub_print <- sub[, c("StepSize_1", "StepSize_2", "n1", "n2", "P", "P_adj", pick_p, pick_stars)]
    names(sub_print)[names(sub_print) == pick_p] <- "P_used_for_stars"
    names(sub_print)[names(sub_print) == pick_stars] <- "Stars"

    sub_print$P <- signif(sub_print$P, 4)
    sub_print$P_adj <- signif(sub_print$P_adj, 4)
    sub_print$P_used_for_stars <- signif(sub_print$P_used_for_stars, 4)

    capture.output(print(sub_print, row.names = FALSE), file = con, append = TRUE)
    writeLines("", con)
  }

  invisible(NULL)
}

# -------------------- MAIN --------------------
run_one_csv <- function(input_file) {
  dat <- read.csv(input_file, stringsAsFactors = FALSE)

  req <- c("Treatment", "StepSize", "Well", "Results", "Biological_Rep")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Missing required column(s): ", paste(miss, collapse = ", "))

  dat$Results <- suppressWarnings(as.numeric(dat$Results))
  dat$Treatment <- factor(dat$Treatment, levels = treatment_levels)
  dat$StepSize <- factor(dat$StepSize, levels = stepsize_levels)
  dat$Biological_Rep <- factor(dat$Biological_Rep)

  dat <- dat[stats::complete.cases(dat[, c("Treatment", "StepSize", "Results", "Biological_Rep")]), , drop = FALSE]
  if (nrow(dat) == 0) stop("No complete rows after filtering: ", basename(input_file))

  if (any(!is.finite(dat$Results) | dat$Results <= -0.5)) {
    stop("Results must be finite and greater than -0.5 for log(Results + 0.5).")
  }
  dat$logResults <- log(dat$Results + 0.5)

  base <- tools::file_path_sans_ext(basename(input_file))
  out_dir <- file.path(dirname(input_file), paste0("MW_StepSizeOnly_Outputs_", base))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (scale_name in c("Raw", "Log")) {
  response_col <- if (scale_name == "Raw") "Results" else "logResults"
  y_label <- if (scale_name == "Raw") y_axis_label_raw else y_axis_label_log

  # ---- MW tables ----
  mw_ctrl <- wilcox_stepsize_within_treatment(
    df = dat,
    treatment_name = "Tri_Ctrl",
    response_col = response_col,
    stepsize_levels_use = stepsize_levels,
    adjust = posthoc_adjust
  )

  mw_lps <- wilcox_stepsize_within_treatment(
    df = dat,
    treatment_name = "Tri_5ugLPS",
    response_col = response_col,
    stepsize_levels_use = stepsize_levels,
    adjust = posthoc_adjust
  )

  mw_all <- dplyr::bind_rows(mw_ctrl, mw_lps)

  utils::write.csv(
    mw_all,
    file.path(out_dir, paste0("MannWhitney_StepSizePairwise_within_Treatment_", scale_name, ".csv")),
    row.names = FALSE
  )
  utils::write.csv(mw_ctrl, file.path(out_dir, paste0("MannWhitney_StepSizePairwise_Ctrl_", scale_name, ".csv")), row.names = FALSE)
  utils::write.csv(mw_lps,  file.path(out_dir, paste0("MannWhitney_StepSizePairwise_5ugLPS_", scale_name, ".csv")), row.names = FALSE)

  write_mw_stepsize_report(
    mw_df = mw_all,
    out_file = file.path(out_dir, paste0("Stat_Results_StepSizePairwise_within_Treatment_", scale_name, ".txt")),
    use_adjusted_p = use_adjusted_p,
    adjust_name = posthoc_adjust,
    alpha = alpha
  )

  # ---- Plots ----
  ctrl <- plot_stepsize_within_treatment(
    df = dat,
    treatment_name = "Tri_Ctrl",
    response_col = response_col,
    y_label = y_label,
    adjust = posthoc_adjust,
    alpha = alpha,
    show_ns = show_ns,
    use_adjusted_p = use_adjusted_p
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("Ctrl_StepSize_Comparison_", scale_name, ".png")),
    ctrl$plot, width = 6.2, height = 4.3, dpi = 300
  )

  lps <- plot_stepsize_within_treatment(
    df = dat,
    treatment_name = "Tri_5ugLPS",
    response_col = response_col,
    y_label = y_label,
    adjust = posthoc_adjust,
    alpha = alpha,
    show_ns = show_ns,
    use_adjusted_p = use_adjusted_p
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("LPS5ug_StepSize_Comparison_", scale_name, ".png")),
    lps$plot, width = 6.2, height = 4.3, dpi = 300
  )
  }

  invisible(list(out_dir = out_dir))
}

csv_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", input_dir)

for (f in csv_files) {
  tryCatch(
    run_one_csv(f),
    error = function(e) message("FAILED: ", basename(f), " -> ", conditionMessage(e))
  )
}
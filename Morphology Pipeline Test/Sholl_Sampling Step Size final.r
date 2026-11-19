# ============================================================
# Standalone: StepSize (S2.5, S5, Continuous) compared within each Treatment
# - Input: a single Excel (.xlsx) file, in either format:
#     * Long: Treatment, StepSize, Well, Results
#     * Wide: Row 1 headers = Key, 2.5, 5, Continuous
#       (Key format: Tri_<Treatment>_<...>, e.g. Tri_5ugLPS_04_2)
# - Output:
#   - 2 plots (Ctrl, 5ugLPS), each with all 3 step-size methods side by side
#   - Pairwise Mann-Whitney (Wilcoxon rank-sum) for all 3 pairs per treatment
#   - Holm adjustment across the 3 pairs within each treatment
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

# -------------------- SETTINGS --------------------
input_file  <- '/Users/allisonpickle/Desktop/260501_Morphology Pipeline/Sholl Analysis/Sholl analysis by step size/Sholl Analysis by Step Size.xlsx'
excel_sheet <- 1                       # sheet name or index to read

treatment_levels <- c("Tri_Ctrl", "Tri_5ugLPS")
stepsize_levels  <- c("S2.5", "S5", "Continuous")  # must match StepSize values
stepsize_labels  <- c("S2.5" = "2.5 µm step", "S5" = "5 µm step", "Continuous" = "Continuous")

alpha          <- 0.05
posthoc_adjust <- "holm"
show_ns        <- FALSE   # TRUE = also draw brackets for non-significant pairs
use_adjusted_p <- TRUE

y_axis_label_raw <- "Results (mean ± SE)"

# -------------------- HELPERS --------------------
p_to_stars <- function(p) {
  if (is.na(p)) return("n.s.")
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  "n.s."
}

add_sig_bracket <- function(p, x1, x2, y, h, label, lw = 0.7) {
  p +
    annotate("segment", x = x1, xend = x1, y = y,     yend = y + h, linewidth = lw) +
    annotate("segment", x = x1, xend = x2, y = y + h, yend = y + h, linewidth = lw) +
    annotate("segment", x = x2, xend = x2, y = y + h, yend = y,     linewidth = lw) +
    annotate("text", x = (x1 + x2) / 2, y = y + h, label = label,
             vjust = -0.2, fontface = "bold")
}

# Map Key tokens (Tri_<Treatment>_<...>) to a treatment name in treatment_levels
parse_treatment <- function(t) {
  if (length(t) < 2) return(NA_character_)
  trt <- t[2]
  if (grepl("ctrl|control", trt, ignore.case = TRUE)) return("Tri_Ctrl")
  if (grepl("5ug|LPS", trt, ignore.case = TRUE))      return("Tri_5ugLPS")
  NA_character_
}

# Read the Excel file; returns long-format data (reshapes wide format if needed)
read_wide_or_long <- function(input_file) {
  raw <- as.data.frame(readxl::read_excel(input_file, sheet = excel_sheet),
                       check.names = FALSE)
  names(raw) <- trimws(names(raw))

  # Already long format
  if (all(c("Treatment", "StepSize", "Results") %in% names(raw))) {
    return(raw)
  }

  # Wide format: row 1 headers = Key, 2.5, 5, Continuous
  if (ncol(raw) < 4) {
    stop("Wide format needs 4 columns: Key, 2.5, 5, Continuous")
  }

  expected <- c("Key", "2.5", "5", "Continuous")
  if (all(expected %in% names(raw))) {
    raw <- raw[, expected, drop = FALSE]          # select by header name
  } else {
    warning("Expected headers (", paste(expected, collapse = ", "),
            ") not all found; using first 4 columns by position.")
    raw <- raw[, 1:4, drop = FALSE]
  }
  names(raw) <- c("Key", "S2.5", "S5", "Continuous")

  raw$Key <- as.character(raw$Key)
  raw <- raw[!is.na(raw$Key) & nzchar(raw$Key), , drop = FALSE]

  # Key format: Tri_<Treatment>_<...>
  tokens <- strsplit(raw$Key, "_")
  raw$Treatment <- vapply(tokens, parse_treatment, character(1))
  raw$Well      <- raw$Key

  for (col in c("S2.5", "S5", "Continuous")) {
    raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
  }

  long <- tidyr::pivot_longer(raw, cols = c("S2.5", "S5", "Continuous"),
                              names_to = "StepSize", values_to = "Results")
  as.data.frame(long)
}

# All pairwise Mann-Whitney tests between step-size methods within one treatment
wilcox_stepsize_pairs_within_treatment <- function(df, treatment_name,
                                                   response_col = "Results",
                                                   stepsize_levels_use = stepsize_levels,
                                                   adjust = "holm") {
  sub_t <- df[df$Treatment == treatment_name, , drop = FALSE]
  prs <- utils::combn(stepsize_levels_use, 2, simplify = FALSE)

  res <- do.call(rbind, lapply(prs, function(pr) {
    x <- sub_t[sub_t$StepSize == pr[1], response_col]
    y <- sub_t[sub_t$StepSize == pr[2], response_col]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    p <- NA_real_
    if (length(x) >= 1 && length(y) >= 1) {
      p <- stats::wilcox.test(x, y, exact = FALSE)$p.value
    }
    data.frame(Treatment = treatment_name,
               StepSize_1 = pr[1], StepSize_2 = pr[2],
               n1 = length(x), n2 = length(y), P = p,
               stringsAsFactors = FALSE)
  }))

  res$P_adj     <- stats::p.adjust(res$P, method = adjust)
  res$Stars_raw <- vapply(res$P,     p_to_stars, character(1))
  res$Stars_adj <- vapply(res$P_adj, p_to_stars, character(1))
  res
}

plot_stepsize_within_treatment <- function(df, treatment_name, mw_tab,
                                           response_col = "Results",
                                           y_label = y_axis_label_raw,
                                           alpha = 0.05, show_ns = FALSE,
                                           use_adjusted_p = TRUE) {
  sub <- df[df$Treatment == treatment_name, , drop = FALSE]
  if (nrow(sub) == 0) stop("No rows for Treatment = ", treatment_name)
  sub$StepSize <- factor(as.character(sub$StepSize), levels = stepsize_levels)

  summ <- sub |>
    group_by(StepSize) |>
    summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      n    = sum(!is.na(.data[[response_col]])),
      .groups = "drop"
    )

  p <- ggplot(summ, aes(x = StepSize, y = mean)) +
    geom_col(fill = "#4C72B0", alpha = 0.75, width = 0.65) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.18) +
    geom_jitter(data = sub,
                aes(x = StepSize, y = .data[[response_col]]),
                width = 0.15, alpha = 0.65, size = 1.2, color = "black") +
    scale_x_discrete(labels = stepsize_labels, drop = FALSE) +
    labs(title = paste0("Step size comparison (", treatment_name, ")"),
         x = "", y = y_label) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  vals  <- sub[[response_col]]
  vals  <- vals[is.finite(vals)]
  y_min <- 0
  y_max <- if (length(vals) > 0) max(vals, na.rm = TRUE) else 1
  span  <- max(y_max - y_min, 1e-9)

  # Decide which pairs get brackets
  mw_tab$P_use <- if (use_adjusted_p) mw_tab$P_adj else mw_tab$P
  mw_tab$Label <- vapply(mw_tab$P_use, p_to_stars, character(1))
  mw_tab$Show  <- is.finite(mw_tab$P_use) & (mw_tab$P_use < alpha | show_ns)
  # Order: narrow brackets first so they stack neatly below wider ones
  mw_tab$x1 <- match(mw_tab$StepSize_1, stepsize_levels)
  mw_tab$x2 <- match(mw_tab$StepSize_2, stepsize_levels)
  mw_tab <- mw_tab[order(mw_tab$x2 - mw_tab$x1, mw_tab$x1), , drop = FALSE]
  show_tab <- mw_tab[mw_tab$Show, , drop = FALSE]

  base_y0 <- y_max + 0.06 * span
  h       <- 0.05 * span
  step    <- 0.14 * span
  n_show  <- nrow(show_tab)
  y_top   <- base_y0 + max(n_show - 1, 0) * step + h * 2.0

  p <- p + scale_y_continuous(limits = c(y_min, y_top),
                              expand = expansion(mult = c(0, 0.02)))

  if (n_show > 0) {
    for (i in seq_len(n_show)) {
      p <- add_sig_bracket(p, x1 = show_tab$x1[i], x2 = show_tab$x2[i],
                           y = base_y0 + (i - 1) * step, h = h,
                           label = show_tab$Label[i])
    }
  }
  p
}

write_mw_stepsize_report <- function(mw_df, out_file, use_adjusted_p = TRUE,
                                     adjust_name = "holm", alpha = 0.05) {
  if (is.null(mw_df) || nrow(mw_df) == 0) {
    writeLines("No Mann-Whitney results to report.", con = out_file)
    return(invisible(NULL))
  }
  pick_p     <- if (use_adjusted_p) "P_adj" else "P"
  pick_stars <- if (use_adjusted_p) "Stars_adj" else "Stars_raw"

  con <- file(out_file, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("STEP SIZE COMPARISONS (S2.5 vs S5 vs Continuous) WITHIN EACH TREATMENT", con)
  writeLines("Test: Mann-Whitney (Wilcoxon rank-sum), two-sided, all pairs", con)
  writeLines(paste0("Multiple-comparison adjustment within treatment: ", adjust_name), con)
  writeLines(paste0("Stars computed from ", pick_p, " (alpha = ", alpha, ")"), con)
  writeLines("", con)

  for (trt in unique(mw_df$Treatment)) {
    sub <- mw_df[mw_df$Treatment == trt, , drop = FALSE]
    writeLines(paste0("Treatment: ", trt), con)
    sub_print <- sub[, c("StepSize_1", "StepSize_2", "n1", "n2", "P", "P_adj", pick_stars)]
    names(sub_print)[names(sub_print) == pick_stars] <- "Stars"
    sub_print$P     <- signif(sub_print$P, 4)
    sub_print$P_adj <- signif(sub_print$P_adj, 4)
    capture.output(print(sub_print, row.names = FALSE), file = con, append = TRUE)
    writeLines("", con)
  }
  invisible(NULL)
}

# -------------------- MAIN --------------------
run_one_file <- function(input_file) {
  dat <- read_wide_or_long(input_file)

  req  <- c("Treatment", "StepSize", "Well", "Results")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Missing required column(s): ", paste(miss, collapse = ", "))

  present <- intersect(stepsize_levels, unique(dat$StepSize))
  absent  <- setdiff(stepsize_levels, present)
  if (length(absent) > 0) {
    stop("StepSize value(s) not found: ", paste(absent, collapse = ", "),
         " (found: ", paste(unique(dat$StepSize), collapse = ", "), ")")
  }

  dat$Results   <- suppressWarnings(as.numeric(dat$Results))
  dat$Treatment <- factor(dat$Treatment, levels = treatment_levels)
  dat$StepSize  <- factor(dat$StepSize, levels = stepsize_levels)

  dat <- dat[complete.cases(dat[, c("Treatment", "StepSize", "Results")]), , drop = FALSE]
  if (nrow(dat) == 0) stop("No complete rows after filtering: ", basename(input_file))

  base    <- tools::file_path_sans_ext(basename(input_file))
  out_dir <- file.path(dirname(input_file), paste0("MW_StepSize_vs_Continuous_Outputs_", base))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- MW tables (one per treatment) ----
  mw_list <- lapply(treatment_levels, function(trt) {
    wilcox_stepsize_pairs_within_treatment(dat, trt, "Results", stepsize_levels, posthoc_adjust)
  })
  names(mw_list) <- treatment_levels
  mw_all <- bind_rows(mw_list)

  write.csv(mw_all, file.path(out_dir, "MannWhitney_StepSize_Pairwise_within_Treatment.csv"),
            row.names = FALSE)
  write_mw_stepsize_report(mw_all,
                           file.path(out_dir, "Stat_Results_StepSize_Pairwise_within_Treatment.txt"),
                           use_adjusted_p, posthoc_adjust, alpha)

  # ---- Plots (one per treatment, all 3 methods side by side) ----
  file_tags <- c("Tri_Ctrl" = "Ctrl", "Tri_5ugLPS" = "LPS5ug")
  for (trt in treatment_levels) {
    p <- plot_stepsize_within_treatment(dat, trt, mw_list[[trt]], "Results",
                                        y_axis_label_raw, alpha, show_ns, use_adjusted_p)
    ggsave(file.path(out_dir, paste0(file_tags[[trt]], "_StepSize_vs_Continuous.png")),
           p, width = 7, height = 4.8, dpi = 300)
  }

  invisible(list(out_dir = out_dir))
}

# -------------------- RUN --------------------
if (!file.exists(input_file)) stop("File not found: ", input_file)

tryCatch(
  run_one_file(input_file),
  error = function(e) message("FAILED: ", basename(input_file), " -> ", conditionMessage(e))
)
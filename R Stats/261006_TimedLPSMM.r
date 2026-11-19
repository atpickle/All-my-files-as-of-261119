# ============================================================
# Mixed Effects: Treatment-only with adaptive nested random effects
# FULL VERSION WITH:
# + Explicit nested image variance modeling
# + Variance component outputs
# + Well-mean plotting (avoids pseudoreplication)
# + Raw + log models
# + Diagnostics
# + LRTs
# + Time-point splitting (6hr, 24hr, etc.)
# + Treatment × Time interaction models
# + Interaction line plots
# + Additional _NoLabels barplots retaining x-axis treatment strips
# ============================================================

# ----------------- PACKAGES -----------------
suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(ggplot2)
  library(dplyr)
})

options(warn = -1)

# ============================================================
# SETTINGS
# ============================================================

input_dir <- '/Users/allisonpickle/Desktop/'
file_pattern <- "\\.csv$"

rename_cols_by_position <- c(
  Results = 6
)

baseline_treatment <- "Tri_Ctrl"

treatment_levels <- c(
  "Tri_Ctrl",
  "5ng_LPS",
  "100ng_LPS",
  "500ng_LPS",
  "5ug_LPS"
)

treatment_labels <- c(
  "Tri_Ctrl"  = "-\n-\n-\n-",
  "5ng_LPS"   = "+\n-\n-\n-",
  "100ng_LPS" = "-\n+\n-\n-",
  "500ng_LPS" = "-\n-\n+\n-",
  "5ug_LPS"   = "-\n-\n-\n+"
)

row_labels <- c(
  "5ng/mL LPS",
  "100 ng/mL LPS",
  "500 ng/mL LPS",
  "5 µg/mL LPS"
)

plot_title_base <- "Treatment"
x_axis_label <- ""
y_axis_label_raw <- "Results (mean ± SE)"
y_axis_label_log <- "log(Results + 0.5) (mean ± SE)"

# ============================================================
# HELPERS
# ============================================================

add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  ""
}

make_aligned_tick_labels <- function(levels_vec, labels_named_vec) {
  out <- rep("", length(levels_vec))
  names(out) <- levels_vec
  common <- intersect(levels_vec, names(labels_named_vec))
  out[common] <- as.character(labels_named_vec[common])
  out
}

save_barplot_versions <- function(filename, plot) {
  ggplot2::ggsave(filename, plot, width = 6, height = 4, dpi = 300)

  plot_nolabels <- plot +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  ggplot2::ggsave(
    sub("\\.png$", "_NoLabels.png", filename),
    plot_nolabels,
    width = 6, height = 4, dpi = 300
  )
}

parse_time_treatment_second_element <- function(x) {
  raw <- trimws(as.character(x))

  time_part <- sub("^([^_]+)_.*$", "\\1", raw)
  trt_part  <- sub("^[^_]+_?", "", raw)

  has_us <- grepl("_", raw, fixed = TRUE)
  time_part[!has_us] <- NA_character_
  trt_part[!has_us]  <- NA_character_

  time_part <- trimws(time_part)
  trt_part  <- trimws(trt_part)

  trt_upper <- toupper(trt_part)
  trt_part[trt_upper %in% c("CTRL", "CONTROL")] <- "Tri_Ctrl"
  trt_part[trt_part == ""] <- NA_character_

  list(Time = time_part, Treatment = trt_part)
}

map_short_trt_to_levels <- function(trt) {
  trt <- as.character(trt)
  map <- c(
    "5ng"   = "5ng_LPS",
    "100ng" = "100ng_LPS",
    "500ng" = "500ng_LPS",
    "5ug"   = "5ug_LPS"
  )
  ifelse(trt %in% names(map), unname(map[trt]), trt)
}

# ============================================================
# MODEL FITTING
# ============================================================

fit_treatment_model <- function(df,
                                response = c("Results", "logResults")) {

  response <- match.arg(response)

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well           <- droplevels(factor(df$Well))
  df$Treatment      <- droplevels(factor(df$Treatment, levels = treatment_levels))
  df$BioWell        <- interaction(df$Biological_Rep, df$Well, drop = TRUE)

  has_img <- "Image" %in% names(df)

  if (has_img) {
    df$Image        <- droplevels(factor(df$Image))
    df$BioWellImage <- interaction(df$Biological_Rep, df$Well, df$Image, drop = TRUE)
    df$WellImage    <- interaction(df$Well, df$Image, drop = TRUE)
  }

  n     <- nrow(df)
  n_bio <- nlevels(df$Biological_Rep)

  can_use_img <- FALSE
  if (has_img) {
    can_use_img <-
      (nlevels(df$BioWellImage) < n) &&
      any(table(df$BioWellImage) >= 2)
  }

  if (n_bio >= 2 && n_bio < n && has_img && can_use_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Biological_Rep) + (1 | BioWell) + (1 | BioWellImage)"
    ))
  } else if (n_bio >= 2 && n_bio < n) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Biological_Rep) + (1 | BioWell)"
    ))
  } else if (has_img && can_use_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + ",
      "(1 | Well) + (1 | WellImage)"
    ))
  } else {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment + (1 | Well)"
    ))
  }

  suppressMessages(suppressWarnings(lmerTest::lmer(f, data = df)))
}

# ============================================================
# INTERACTION MODEL FITTING
# ============================================================

fit_interaction_model <- function(df, response = c("Results", "logResults")) {

  response <- match.arg(response)

  df$Biological_Rep <- droplevels(factor(df$Biological_Rep))
  df$Well           <- droplevels(factor(df$Well))
  df$Treatment      <- droplevels(factor(df$Treatment, levels = treatment_levels))
  df$Time           <- droplevels(factor(df$Time))

  df$BioWell <- interaction(df$Biological_Rep, df$Well, drop = TRUE)

  has_img <- "Image" %in% names(df)

  if (has_img) {
    df$Image        <- droplevels(factor(df$Image))
    df$BioWellImage <- interaction(df$Biological_Rep, df$Well, df$Image, drop = TRUE)
    df$WellImage    <- interaction(df$Well, df$Image, drop = TRUE)
  }

  n     <- nrow(df)
  n_bio <- nlevels(df$Biological_Rep)

  can_use_img <- FALSE
  if (has_img) {
    can_use_img <- (nlevels(df$BioWellImage) < n) && any(table(df$BioWellImage) >= 2)
  }

  # Build interaction formula: Treatment * Time + random effects
  if (n_bio >= 2 && n_bio < n && has_img && can_use_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Time + ",
      "(1 | Biological_Rep) + (1 | BioWell) + (1 | BioWellImage)"
    ))
  } else if (n_bio >= 2 && n_bio < n) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Time + ",
      "(1 | Biological_Rep) + (1 | BioWell)"
    ))
  } else if (has_img && can_use_img) {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Time + ",
      "(1 | Well) + (1 | WellImage)"
    ))
  } else {
    f <- stats::as.formula(paste0(
      response, " ~ Treatment * Time + (1 | Well)"
    ))
  }

  suppressMessages(suppressWarnings(lmerTest::lmer(f, data = df)))
}

# ============================================================
# NULL VS FULL MODEL
# ============================================================

make_null_full <- function(fit, response_name, df) {

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

  null_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(make_form("1"), data = df, REML = FALSE)
  ))
  full_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(make_form("Treatment"), data = df, REML = FALSE)
  ))

  list(null = null_fit, full = full_fit)
}

# ============================================================
# VARIANCE COMPONENTS
# ============================================================

compute_variance_components <- function(fit) {
  vc        <- as.data.frame(lme4::VarCorr(fit))
  var_resid <- attr(lme4::VarCorr(fit), "sc")^2

  vc_out <- vc[, c("grp", "vcov")]
  names(vc_out) <- c("Component", "Variance")
  vc_out <- rbind(vc_out, data.frame(Component = "Residual", Variance = var_resid))

  vc_out$Percent <- 100 * vc_out$Variance / sum(vc_out$Variance)
  vc_out
}

# ============================================================
# P-VALUES
# ============================================================

extract_treatment_pvals <- function(fit, levels_vec) {

  fe  <- coef(summary(fit))
  rn  <- rownames(fe)

  is_trt    <- grepl("^Treatment", rn) & !grepl(":", rn)
  trt_rows  <- fe[is_trt, , drop = FALSE]
  trt_names <- sub("^Treatment", "", rownames(trt_rows))

  p_col <- grep("Pr\\(", colnames(trt_rows), value = TRUE)[1]
  pvals <- trt_rows[, p_col]
  names(pvals) <- trt_names

  p_df <- data.frame(Treatment = levels_vec, p.value = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(trt_names)) {
    p_df$p.value[p_df$Treatment == trt_names[i]] <- as.numeric(pvals[i])
  }
  p_df
}

extract_interaction_pvals <- function(fit) {

  fe <- coef(summary(fit))
  rn <- rownames(fe)

  # Main effects
  is_trt  <- grepl("^Treatment", rn) & !grepl(":", rn)
  is_time <- grepl("^Time", rn) & !grepl(":", rn)

  # Interaction terms: Treatment:Time
  is_int <- grepl("^Treatment.*:Time", rn) | grepl("^Time.*:Treatment", rn)

  p_col <- grep("Pr\\(", colnames(fe), value = TRUE)[1]

  list(
    Treatment   = if (any(is_trt))  fe[is_trt, p_col, drop = FALSE] else NULL,
    Time        = if (any(is_time)) fe[is_time, p_col, drop = FALSE] else NULL,
    Interaction = if (any(is_int))  fe[is_int, p_col, drop = FALSE] else NULL
  )
}

# ============================================================
# SUMMARY TABLES WITH TIME-SPECIFIC CONTROL COMPARISONS
# ============================================================

build_bar_summary_with_sig_by_time <- function(data, fit, response_col, time_label) {

  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) /
               sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  summary_df$Treatment <- factor(summary_df$Treatment, levels = levels(data$Treatment))

  # Extract p-values from model (comparing each treatment to baseline)
  p_df <- extract_treatment_pvals(fit, levels(data$Treatment))
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out        <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars  <- vapply(out$p.value, add_stars, character(1))
  out$Time   <- time_label
  out
}

# ============================================================
# INDIVIDUAL PLOTS PER TIME POINT WITH SIGNIFICANCE STARS
# ============================================================

plot_individual_time_points <- function(data, response_col, y_label, title_base) {

  # Compute well-level means for plotting
  plot_data <- data %>%
    group_by(Time, Treatment, Biological_Rep, Well) %>%
    summarise(
      PlotValue = mean(.data[[response_col]], na.rm = TRUE),
      .groups = "drop"
    )

  n_bio_plot  <- dplyr::n_distinct(plot_data$Biological_Rep, na.rm = TRUE)
  show_legend <- n_bio_plot > 1

  aligned_treatment_labels <- make_aligned_tick_labels(treatment_levels, treatment_labels)
  aligned_treatment_labels <- stats::setNames(
    paste0("\n", aligned_treatment_labels),
    names(aligned_treatment_labels)
  )

  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))

  # Split by time and create one plot per time point
  time_points <- sort(unique(as.character(data$Time)))
  
  plots <- lapply(time_points, function(tm) {
    
    # Subset data for this time point
    tm_data <- data[data$Time == tm, , drop = FALSE]
    plot_tm <- plot_data[plot_data$Time == tm, ]
    
    # Fit model for this time point only (Treatment vs. control within time)
    tm_fit <- tryCatch(
      fit_treatment_model(tm_data, response = response_col),
      error = function(e) {
        message("  Model failed for ", tm, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(tm_fit)) {
      return(list(time = tm, plot = NULL))
    }
    
    # Get summary with significance stars (vs. control at this time point)
    summ_tm <- build_bar_summary_with_sig_by_time(tm_data, tm_fit, response_col, tm)
    
    p <- ggplot2::ggplot(summ_tm, ggplot2::aes(x = Treatment, y = mean)) +
      ggplot2::geom_col(fill = "#6c6c6c", alpha = 0.7) +
      ggplot2::geom_jitter(
        data = plot_tm,
        ggplot2::aes(x = Treatment, y = PlotValue, color = Biological_Rep),
        shape = 16,
        width = 0.15,
        alpha = 0.85,
        size = 1
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = mean - se, ymax = mean + se),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = mean + se, label = stars),
        vjust = -0.5, size = 5, fontface = "bold"
      ) +
      ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        title = paste0(title_base, " \u2013 ", tm),
        x = x_axis_label,
        y = y_label,
        color = if (show_legend) "Biological Rep" else NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5,
                                                  vjust = -0.2, lineheight = 0.8, size = 9),
        plot.margin      = ggplot2::margin(5, 5, 45, 55),
        legend.position  = if (show_legend) "right" else "none"
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::annotate(
        "text", x = 0.4, y = -Inf, label = combined_row_label,
        hjust = 1, vjust = 1.3, size = 3, lineheight = 0.8
      )
    
    list(time = tm, plot = p, summary = summ_tm)
  })
  
  # Filter out any failed plots
  plots[!sapply(plots, function(x) is.null(x$plot))]
}

# ============================================================
# PAIRWISE TIME COMPARISONS PER TREATMENT
# ============================================================

compute_consecutive_time_pvals <- function(data, response_col) {
  
  # Order time points chronologically
  time_order <- c("6hr", "24hr", "48hr", "72hr", "96hr")
  available_times <- intersect(time_order, unique(as.character(data$Time)))
  
  if (length(available_times) < 2) {
    return(NULL)  # Need at least 2 time points for comparison
  }
  
  # Results storage
  pval_results <- list()
  
  # Loop through each treatment
  for (trt in treatment_levels) {
    
    trt_data <- data[data$Treatment == trt, ]
    
    if (nrow(trt_data) == 0) next
    
    # Compare consecutive time points
    for (i in 2:length(available_times)) {
      
      time_curr <- available_times[i]
      time_prev <- available_times[i - 1]
      
      # Subset data for these two time points
      compare_data <- trt_data[trt_data$Time %in% c(time_prev, time_curr), ]
      
      if (nrow(compare_data) < 4) next  # Need sufficient data
      
      # Fit mixed model comparing the two time points
      compare_data$Time <- factor(compare_data$Time, levels = c(time_prev, time_curr))
      compare_data$BioWell <- interaction(compare_data$Biological_Rep, compare_data$Well, drop = TRUE)
      
      # Build formula based on available data
      has_img <- "Image" %in% names(compare_data)
      n_bio <- length(unique(compare_data$Biological_Rep))
      
      if (has_img) {
        compare_data$BioWellImage <- interaction(
          compare_data$Biological_Rep, 
          compare_data$Well, 
          compare_data$Image, 
          drop = TRUE
        )
      }
      
      # Simplified model for pairwise comparison
      if (n_bio >= 2 && has_img) {
        f <- stats::as.formula(paste0(
          response_col, " ~ Time + (1 | Biological_Rep) + (1 | BioWell)"
        ))
      } else if (n_bio >= 2) {
        f <- stats::as.formula(paste0(
          response_col, " ~ Time + (1 | Biological_Rep) + (1 | BioWell)"
        ))
      } else {
        f <- stats::as.formula(paste0(
          response_col, " ~ Time + (1 | Well)"
        ))
      }
      
      # Fit model
      fit <- tryCatch(
        suppressMessages(suppressWarnings(lmerTest::lmer(f, data = compare_data))),
        error = function(e) NULL
      )
      
      if (is.null(fit)) next
      
      # Extract p-value for Time effect
      fe <- coef(summary(fit))
      time_row <- grepl("^Time", rownames(fe))
      
      if (any(time_row)) {
        p_col <- grep("Pr\\(", colnames(fe), value = TRUE)[1]
        pval <- as.numeric(fe[time_row, p_col][1])
        
        pval_results[[length(pval_results) + 1]] <- data.frame(
          Treatment = trt,
          Time_From = time_prev,
          Time_To = time_curr,
          p.value = pval,
          stars = add_stars(pval),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(pval_results) == 0) return(NULL)
  
  do.call(rbind, pval_results)
}

# ============================================================
# TREATMENT × TIME INTERACTION LINE PLOT WITH SIGNIFICANCE
# ============================================================

plot_treatment_time_interaction <- function(data, response_col, y_label, title_text) {

  # Compute summary statistics per Treatment × Time
  summ_df <- data %>%
    group_by(Time, Treatment) %>%
    summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  # Order time points chronologically
  time_order <- c("6hr", "24hr", "48hr", "72hr", "96hr")
  available_times <- intersect(time_order, unique(as.character(summ_df$Time)))
  summ_df$Time <- factor(summ_df$Time, levels = available_times)
  
  # Filter out any NA time points
  summ_df <- summ_df[!is.na(summ_df$Time), ]
  
  # Compute consecutive time point p-values
  cat("    Computing consecutive time comparisons per treatment...\n")
  pval_df <- compute_consecutive_time_pvals(data, response_col)
  
  if (!is.null(pval_df)) {
    cat("    Found", nrow(pval_df), "significant comparisons\n")
    print(pval_df[pval_df$stars != "", ])
  }

  # Treatment colors and shapes for clarity
  treatment_colors <- c(
    "Tri_Ctrl"  = "#999999",
    "5ng_LPS"   = "#E69F00",
    "100ng_LPS" = "#56B4E9",
    "500ng_LPS" = "#009E73",
    "5ug_LPS"   = "#D55E00"
  )
  
  treatment_shapes <- c(
    "Tri_Ctrl"  = 15,  # square
    "5ng_LPS"   = 16,  # circle
    "100ng_LPS" = 17,  # triangle
    "500ng_LPS" = 18,  # diamond
    "5ug_LPS"   = 19   # filled circle
  )

  # Treatment labels for legend
  treatment_legend_labels <- c(
    "Tri_Ctrl"  = "Control",
    "5ng_LPS"   = "5 ng/mL LPS",
    "100ng_LPS" = "100 ng/mL LPS",
    "500ng_LPS" = "500 ng/mL LPS",
    "5ug_LPS"   = "5 µg/mL LPS"
  )

  # Base plot
  p <- ggplot2::ggplot(summ_df, ggplot2::aes(x = Time, y = mean, 
                                              color = Treatment, 
                                              group = Treatment,
                                              shape = Treatment)) +
    ggplot2::geom_line(size = 1, alpha = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.2,
      alpha = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = treatment_colors,
      labels = treatment_legend_labels,
      name = "Treatment"
    ) +
    ggplot2::scale_shape_manual(
      values = treatment_shapes,
      labels = treatment_legend_labels,
      name = "Treatment"
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title = title_text,
      x = "Time Point",
      y = y_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5, size = 11),
      axis.text.y      = ggplot2::element_text(size = 11),
      axis.title       = ggplot2::element_text(size = 12, face = "bold"),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(face = "bold"),
      legend.text      = ggplot2::element_text(size = 10),
      plot.margin      = ggplot2::margin(10, 10, 10, 10)
    )
  
  # Add significance stars if available
  if (!is.null(pval_df) && nrow(pval_df) > 0) {
    
    # Filter to only significant comparisons
    sig_df <- pval_df[pval_df$stars != "", ]
    
    if (nrow(sig_df) > 0) {
      
      # Prepare annotation data
      sig_df$Time_From <- factor(sig_df$Time_From, levels = available_times)
      sig_df$Time_To <- factor(sig_df$Time_To, levels = available_times)
      
      # Get x positions for time points (numeric)
      sig_df$x_from <- as.numeric(sig_df$Time_From)
      sig_df$x_to <- as.numeric(sig_df$Time_To)
      sig_df$x_mid <- (sig_df$x_from + sig_df$x_to) / 2
      
      y_range <- diff(range(summ_df$mean, na.rm = TRUE))
      if (!is.finite(y_range) || y_range == 0) y_range <- 1

      # Get y positions: place above the HIGHER of the two time-point means (+ SE),
      # and stack comparisons that share the same x_mid so stars/brackets don't overlap
      sig_df$y_pos <- NA_real_

      for (i in seq_len(nrow(sig_df))) {
        trt <- sig_df$Treatment[i]
        time_to   <- as.character(sig_df$Time_To[i])
        time_from <- as.character(sig_df$Time_From[i])

        row_to   <- summ_df[summ_df$Treatment == trt & as.character(summ_df$Time) == time_to, ]
        row_from <- summ_df[summ_df$Treatment == trt & as.character(summ_df$Time) == time_from, ]

        top_to   <- if (nrow(row_to)   > 0) row_to$mean   + row_to$se   else NA_real_
        top_from <- if (nrow(row_from) > 0) row_from$mean + row_from$se else NA_real_

        top_val <- max(top_to, top_from, na.rm = TRUE)

        if (is.finite(top_val)) {
          sig_df$y_pos[i] <- top_val + 0.08 * y_range
        }
      }

      # Stack overlapping annotations (same x_mid) vertically to avoid collisions.
      # Order from bottom to top by treatment rank: Ctrl, 5ng, 100ng, 500ng, 5ug
      sig_df <- sig_df[order(sig_df$x_mid), ]
      sig_df$trt_rank <- match(as.character(sig_df$Treatment), treatment_levels)

      dup_groups <- split(seq_len(nrow(sig_df)), sig_df$x_mid)
      for (grp in dup_groups) {
        if (length(grp) > 1) {
          ord <- grp[order(sig_df$trt_rank[grp])]
          for (k in seq_along(ord)) {
            # use (k) instead of (k - 1) so even the bottom item gets lifted off its own data
            sig_df$y_pos[ord[k]] <- sig_df$y_pos[ord[k]] + k * 0.22 * y_range
          }
        }
      }

      # Add text annotations for stars (raised slightly above the bracket)
      p <- p + ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(x = x_mid, y = y_pos + 0.05 * y_range, label = stars, color = Treatment),
        inherit.aes = FALSE,
        size = 4,
        fontface = "bold",
        show.legend = FALSE
      )
      
      # Add small brackets to show comparison (anchored at the new y_pos)
      for (i in 1:nrow(sig_df)) {
        bracket_data <- data.frame(
          x = c(sig_df$x_from[i], sig_df$x_from[i], sig_df$x_to[i], sig_df$x_to[i]),
          y = c(sig_df$y_pos[i],
                sig_df$y_pos[i] + 0.025 * y_range,
                sig_df$y_pos[i] + 0.025 * y_range,
                sig_df$y_pos[i]),
          Treatment = sig_df$Treatment[i]
        )
        
        p <- p + ggplot2::geom_path(
          data = bracket_data,
          ggplot2::aes(x = x, y = y, color = Treatment, group = Treatment),
          inherit.aes = FALSE,
          size = 0.5,
          show.legend = FALSE
        )
      }

      # Expand y-axis slightly so stacked stars aren't clipped
      max_y_pos <- max(sig_df$y_pos, na.rm = TRUE) + 0.12 * y_range
      p <- p + ggplot2::coord_cartesian(ylim = c(NA, max_y_pos), clip = "off")
    }
  }
  
  p
}

# ============================================================
# DIAGNOSTICS
# ============================================================

save_diagnostics <- function(out_dir, raw_fit, log_fit) {

  models <- list(RawFit = raw_fit, LogFit = log_fit)

  for (nm in names(models)) {
    fit           <- models[[nm]]
    fitted_vals   <- as.numeric(stats::fitted(fit))
    pearson_resid <- as.numeric(stats::resid(fit, type = "pearson"))
    raw_resid     <- as.numeric(stats::resid(fit))

    res_df <- data.frame(FittedValues = fitted_vals, PearsonResiduals = pearson_resid)

    p_res <- ggplot2::ggplot(res_df, ggplot2::aes(x = FittedValues, y = PearsonResiduals)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, color = "red") +
      ggplot2::theme_bw()

    ggplot2::ggsave(file.path(out_dir, paste0(nm, "_residuals.png")),
                    p_res, width = 7, height = 5, dpi = 200)

    qq    <- stats::qqnorm(raw_resid, plot.it = FALSE)
    qq_df <- data.frame(x = qq$x, y = qq$y)

    p_qq <- ggplot2::ggplot(qq_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::theme_bw()

    ggplot2::ggsave(file.path(out_dir, paste0(nm, "_qq.png")),
                    p_qq, width = 7, height = 5, dpi = 200)
  }
}

# ============================================================
# MODEL SUMMARIES
# ============================================================

save_model_summaries <- function(out_dir, input_file, analysis_data,
                                 raw_fit, log_fit,
                                 vc_raw, vc_log, lrt_raw, lrt_log) {

  out_file <- file.path(out_dir, "Mixed_Model_Summaries.txt")
  sink(out_file)

  cat("MIXED EFFECTS MODEL SUMMARIES\n\n")
  cat("Input: ", input_file, "\n\n", sep = "")

  cat("RAW MODEL\n\n");  print(summary(raw_fit))
  cat("\nVariance Components (raw)\n"); print(vc_raw)
  cat("\nLRT (raw)\n");  print(lrt_raw)

  cat("\n\nLOG MODEL\n\n"); print(summary(log_fit))
  cat("\nVariance Components (log)\n"); print(vc_log)
  cat("\nLRT (log)\n");  print(lrt_log)

  sink()
}

save_interaction_summaries <- function(out_dir, input_file, analysis_data,
                                       raw_fit, log_fit,
                                       vc_raw, vc_log,
                                       pvals_raw, pvals_log) {

  out_file <- file.path(out_dir, "Interaction_Model_Summaries.txt")
  sink(out_file)

  cat("TREATMENT × TIME INTERACTION MODEL SUMMARIES\n\n")
  cat("Input: ", input_file, "\n\n", sep = "")

  cat("RAW MODEL (Treatment * Time)\n\n")
  print(summary(raw_fit))
  cat("\nVariance Components (raw)\n"); print(vc_raw)
  cat("\nP-values (raw)\n")
  cat("  Treatment:\n");   print(pvals_raw$Treatment)
  cat("  Time:\n");        print(pvals_raw$Time)
  cat("  Interaction:\n"); print(pvals_raw$Interaction)

  cat("\n\nLOG MODEL (Treatment * Time)\n\n")
  print(summary(log_fit))
  cat("\nVariance Components (log)\n"); print(vc_log)
  cat("\nP-values (log)\n")
  cat("  Treatment:\n");   print(pvals_log$Treatment)
  cat("  Time:\n");        print(pvals_log$Time)
  cat("  Interaction:\n"); print(pvals_log$Interaction)

  sink()
}


# ============================================================
# PER-TIME-POINT ANALYSIS
# ============================================================

run_one_timepoint <- function(tm_data, time_label, out_dir_base, input_file) {

  cat("  Time point:", time_label, "| n =", nrow(tm_data), "\n")

  tm_data$BioWell <- interaction(tm_data$Biological_Rep, tm_data$Well, drop = TRUE)

  if ("Image" %in% names(tm_data)) {
    tm_data$BioWellImage <- interaction(
      tm_data$Biological_Rep, tm_data$Well, tm_data$Image, drop = TRUE)
    tm_data$WellImage <- interaction(tm_data$Well, tm_data$Image, drop = TRUE)
  }

  out_dir <- file.path(out_dir_base, time_label)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  image_variance <- tm_data %>%
    dplyr::group_by(Treatment, Well) %>%
    dplyr::summarise(
      image_mean = mean(Results, na.rm = TRUE),
      image_sd   = sd(Results, na.rm = TRUE),
      image_cv   = dplyr::if_else(image_mean == 0, NA_real_, 100 * image_sd / image_mean),
      n_images   = dplyr::n(),
      .groups = "drop"
    )

  utils::write.csv(image_variance,
                   file.path(out_dir, "ImageVariation_ByTreatmentWell_Raw.csv"),
                   row.names = FALSE)

  raw_fit <- fit_treatment_model(tm_data, response = "Results")
  log_fit <- fit_treatment_model(tm_data, response = "logResults")

  nf_raw  <- make_null_full(raw_fit, "Results",    tm_data)
  lrt_raw <- tryCatch(anova(nf_raw$null, nf_raw$full),  error = function(e) NULL)

  nf_log  <- make_null_full(log_fit, "logResults", tm_data)
  lrt_log <- tryCatch(anova(nf_log$null, nf_log$full),  error = function(e) NULL)

  vc_raw <- compute_variance_components(raw_fit)
  vc_log <- compute_variance_components(log_fit)

  write.csv(vc_raw, file.path(out_dir, "VarianceComponents_Raw.csv"), row.names = FALSE)
  write.csv(vc_log, file.path(out_dir, "VarianceComponents_Log.csv"), row.names = FALSE)

  # -------------------------------------------------------
  # BAR PLOTS FOR THIS TIME POINT (treatments vs. control)
  # -------------------------------------------------------
  
  # Compute well-level means for plotting
  plot_data_raw <- tm_data %>%
    group_by(Treatment, Biological_Rep, Well) %>%
    summarise(
      PlotValue = mean(Results, na.rm = TRUE),
      .groups = "drop"
    )
  
  plot_data_log <- tm_data %>%
    group_by(Treatment, Biological_Rep, Well) %>%
    summarise(
      PlotValue = mean(logResults, na.rm = TRUE),
      .groups = "drop"
    )
  
  n_bio_plot  <- dplyr::n_distinct(plot_data_raw$Biological_Rep, na.rm = TRUE)
  show_legend <- n_bio_plot > 1
  
  aligned_treatment_labels <- make_aligned_tick_labels(treatment_levels, treatment_labels)
  aligned_treatment_labels <- stats::setNames(
    paste0("\n", aligned_treatment_labels),
    names(aligned_treatment_labels)
  )
  
  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  
  # Raw plot
  summ_raw <- build_bar_summary_with_sig_by_time(tm_data, raw_fit, "Results", time_label)
  
  p_raw <- ggplot2::ggplot(summ_raw, ggplot2::aes(x = Treatment, y = mean)) +
    ggplot2::geom_col(fill = "#6c6c6c", alpha = 0.7) +
    ggplot2::geom_jitter(
      data = plot_data_raw,
      ggplot2::aes(x = Treatment, y = PlotValue, color = Biological_Rep),
      shape = 16,
      width = 0.15,
      alpha = 0.85,
      size = 1
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = mean + se, label = stars),
      vjust = -0.5, size = 5, fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = paste0(plot_title_base, " \u2013 ", time_label, " (raw)"),
      x = x_axis_label,
      y = y_axis_label_raw,
      color = if (show_legend) "Biological Rep" else NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5,
                                                vjust = -0.2, lineheight = 0.8, size = 9),
      plot.margin      = ggplot2::margin(5, 5, 45, 55),
      legend.position  = if (show_legend) "right" else "none"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotate(
      "text", x = 0.4, y = -Inf, label = combined_row_label,
      hjust = 1, vjust = 1.3, size = 3, lineheight = 0.8
    )
  
  save_barplot_versions(file.path(out_dir, "Barplot_RawFit.png"), p_raw)
  
  # Log plot
  summ_log <- build_bar_summary_with_sig_by_time(tm_data, log_fit, "logResults", time_label)
  
  p_log <- ggplot2::ggplot(summ_log, ggplot2::aes(x = Treatment, y = mean)) +
    ggplot2::geom_col(fill = "#6c6c6c", alpha = 0.7) +
    ggplot2::geom_jitter(
      data = plot_data_log,
      ggplot2::aes(x = Treatment, y = PlotValue, color = Biological_Rep),
      shape = 16,
      width = 0.15,
      alpha = 0.85,
      size = 1
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se, ymax = mean + se),
      width = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = mean + se, label = stars),
      vjust = -0.5, size = 5, fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(labels = aligned_treatment_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10))) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = paste0(plot_title_base, " \u2013 ", time_label, " (log)"),
      x = x_axis_label,
      y = y_axis_label_log,
      color = if (show_legend) "Biological Rep" else NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5,
                                                vjust = -0.2, lineheight = 0.8, size = 9),
      plot.margin      = ggplot2::margin(5, 5, 45, 55),
      legend.position  = if (show_legend) "right" else "none"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotate(
      "text", x = 0.4, y = -Inf, label = combined_row_label,
      hjust = 1, vjust = 1.3, size = 3, lineheight = 0.8
    )
  
  save_barplot_versions(file.path(out_dir, "Barplot_LogFit.png"), p_log)

  save_diagnostics(out_dir, raw_fit, log_fit)
  save_model_summaries(out_dir, input_file, tm_data,
                       raw_fit, log_fit,
                       vc_raw, vc_log, lrt_raw, lrt_log)

  cat("    Saved to:", out_dir, "\n")
}


# ============================================================
# MAIN PER-FILE RUNNER
# ============================================================

run_one_csv <- function(input_file) {

  cat("\nProcessing:", basename(input_file), "\n")

  process_data <- read.csv(input_file, stringsAsFactors = FALSE)

  if (ncol(process_data) >= 5) {
    colnames(process_data)[3] <- "Biological_Rep"
    colnames(process_data)[4] <- "Well"
    colnames(process_data)[5] <- "Image"
  } else {
    stop("Expected >= 5 columns. Found: ", ncol(process_data))
  }

  if (!is.null(rename_cols_by_position)) {
    colnames(process_data)[rename_cols_by_position[["Results"]]] <- "Results"
  }

  if (!"Treatment" %in% names(process_data)) {
    stop("No 'Treatment' column found.")
  }

  tt <- parse_time_treatment_second_element(process_data$Treatment)
  process_data$Time      <- tt$Time
  process_data$Treatment <- map_short_trt_to_levels(tt$Treatment)

  process_data$Results <- suppressWarnings(as.numeric(process_data$Results))

  process_data$Treatment <- factor(process_data$Treatment, levels = treatment_levels)
  process_data$Treatment <- stats::relevel(process_data$Treatment, ref = baseline_treatment)

  process_data$Well           <- factor(process_data$Well)
  process_data$Biological_Rep <- factor(process_data$Biological_Rep)

  if ("Image" %in% names(process_data)) {
    process_data$Image <- factor(process_data$Image)
  }

  needed <- c("Results", "Treatment", "Well", "Biological_Rep", "Time")
  if ("Image" %in% names(process_data)) needed <- c(needed, "Image")

  analysis_data <- process_data[complete.cases(process_data[, needed]), ]

  if (nrow(analysis_data) == 0) {
    stop("0 complete-case rows after parsing.")
  }

  analysis_data$logResults <- log(analysis_data$Results + 0.5)

  csv_base     <- tools::file_path_sans_ext(basename(input_file))
  out_dir_base <- file.path(dirname(input_file), paste0("MM_Outputs_", csv_base))
  dir.create(out_dir_base, showWarnings = FALSE, recursive = TRUE)

  # -------------------------------------------------------
  # 1. PER-TIME-POINT ANALYSIS
  # -------------------------------------------------------
  time_points <- sort(unique(as.character(analysis_data$Time)))
  cat("  Time points found:", paste(time_points, collapse = ", "), "\n")

  for (tm in time_points) {

    tm_data <- analysis_data[analysis_data$Time == tm, , drop = FALSE]

    if (nrow(tm_data) < 5) {
      message("  Skipping ", tm, ": too few rows (", nrow(tm_data), ")")
      next
    }

    tryCatch(
      run_one_timepoint(tm_data, tm, out_dir_base, input_file),
      error = function(e) {
        message("  FAILED time point ", tm, ": ", conditionMessage(e))
      }
    )
  }

  # -------------------------------------------------------
  # 2. TREATMENT × TIME INTERACTION (combined analysis)
  # -------------------------------------------------------
  cat("  Running Treatment × Time interaction model...\n")

  analysis_data$BioWell <- interaction(
    analysis_data$Biological_Rep,
    analysis_data$Well,
    drop = TRUE
  )

  if ("Image" %in% names(analysis_data)) {
    analysis_data$BioWellImage <- interaction(
      analysis_data$Biological_Rep,
      analysis_data$Well,
      analysis_data$Image,
      drop = TRUE
    )
    analysis_data$WellImage <- interaction(
      analysis_data$Well,
      analysis_data$Image,
      drop = TRUE
    )
  }

  int_raw <- tryCatch(
    fit_interaction_model(analysis_data, response = "Results"),
    error = function(e) { message("  Interaction raw model failed: ", e$message); NULL }
  )

  int_log <- tryCatch(
    fit_interaction_model(analysis_data, response = "logResults"),
    error = function(e) { message("  Interaction log model failed: ", e$message); NULL }
  )

  # Save interaction model summaries (if both models succeeded)
  if (!is.null(int_raw) && !is.null(int_log)) {

    vc_int_raw <- compute_variance_components(int_raw)
    vc_int_log <- compute_variance_components(int_log)

    write.csv(vc_int_raw, file.path(out_dir_base, "VarianceComponents_Interaction_Raw.csv"), row.names = FALSE)
    write.csv(vc_int_log, file.path(out_dir_base, "VarianceComponents_Interaction_Log.csv"), row.names = FALSE)

    pvals_raw <- extract_interaction_pvals(int_raw)
    pvals_log <- extract_interaction_pvals(int_log)

    save_interaction_summaries(
      out_dir_base, input_file, analysis_data,
      int_raw, int_log,
      vc_int_raw, vc_int_log,
      pvals_raw, pvals_log
    )
  }

  # -------------------------------------------------------
  # 3. INDIVIDUAL TIME POINT PLOTS (always generated)
  # -------------------------------------------------------
  cat("  Generating individual time point plots...\n")
  
  plots_raw <- plot_individual_time_points(
    analysis_data, "Results",
    y_axis_label_raw,
    paste0(plot_title_base, " (raw)")
  )
  
  for (item in plots_raw) {
    save_barplot_versions(
      file.path(out_dir_base, paste0("AllTimeComparison_", item$time, "_Raw.png")),
      item$plot
    )
  }

  plots_log <- plot_individual_time_points(
    analysis_data, "logResults",
    y_axis_label_log,
    paste0(plot_title_base, " (log)")
  )
  
  for (item in plots_log) {
    save_barplot_versions(
      file.path(out_dir_base, paste0("AllTimeComparison_", item$time, "_Log.png")),
      item$plot
    )
  }

  # -------------------------------------------------------
  # 4. TREATMENT × TIME INTERACTION LINE PLOTS (always generated)
  # -------------------------------------------------------
  cat("  Generating Treatment × Time interaction line plots...\n")

  p_interaction_raw <- plot_treatment_time_interaction(
    analysis_data, "Results",
    y_axis_label_raw,
    "Treatment × Time Interaction (raw)"
  )
  ggplot2::ggsave(
    file.path(out_dir_base, "Interaction_TreatmentByTime_Raw.png"),
    p_interaction_raw,
    width = 8, height = 6, dpi = 300
  )

  p_interaction_log <- plot_treatment_time_interaction(
    analysis_data, "logResults",
    y_axis_label_log,
    "Treatment × Time Interaction (log)"
  )
  ggplot2::ggsave(
    file.path(out_dir_base, "Interaction_TreatmentByTime_Log.png"),
    p_interaction_log,
    width = 8, height = 6, dpi = 300
  )

  cat("  All plots saved to:", out_dir_base, "\n")
  cat("Done:", basename(input_file), "\n")
}

# ============================================================
# RUN ALL FILES
# ============================================================

csv_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)

if (length(csv_files) == 0) stop("No CSV files found in: ", input_dir)

# Skip any previously generated output CSVs
csv_files <- csv_files[!grepl("MM_Outputs_", csv_files, fixed = TRUE)]
csv_files <- csv_files[!grepl("^MM_", basename(csv_files))]

cat("Found ", length(csv_files), " file(s)\n", sep = "")

results <- lapply(
  csv_files,
  function(f) {
    tryCatch(
      run_one_csv(f),
      error = function(e) {
        message("FAILED: ", basename(f), " -> ", conditionMessage(e))
        NULL
      }
    )
  }
)

cat("\n=== All files processed ===\n")
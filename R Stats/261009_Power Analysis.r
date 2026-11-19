# ============================================================
# MICROGLIAL MORPHOLOGY POWER ANALYSIS
#
# Estimates power for LPS-versus-control comparisons using:
# - A hierarchical linear mixed-effects model
# - Biological replicate as the independent experimental unit
# - Well and observation-level variability in replicate means
# - Single-step Dunnett adjustment across LPS concentrations
# - Sensitivity analyses for treatment effect and variability
# Main minimum-replicate tables contain one row per LPS concentration
# using the observed effect and variability. Sensitivity scenarios
# are saved separately in MinimumBiologicalReplicates_Sensitivity files.
# Each metric is analyzed both by time point and with all times pooled
# by concentration in Combined_All_Timepoints. The pooled model omits
# Time and does not count repeated times as new biological replicates.
# Power_Analysis_Summaries contains a master minimum-replicate table
# and one summary per morphology metric, excluding sensitivity scenarios.
#
# The simulation works at the biological-replicate treatment-mean
# level. Identifiable lower-level variance components are propagated
# into the variance of each replicate mean; the shared biological-
# replicate random intercept cancels in within-replicate contrasts.
# A separate image component is included only when images have
# repeated observations. Otherwise, image-to-image variation remains
# part of the row-level residual and cannot be separated from cell
# variability unless each input row is an individual microglia.
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required_packages <- c(
  "lmerTest",
  "lme4",
  "emmeans",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following required R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}


# ============================================================
# SETTINGS
# ============================================================

input_dir <- path.expand(
  "/Users/allisonpickle/Desktop/LPS/3 BioReps Analysis/All BioReps Microglia Morphology FINAL"
)

file_pattern <- "\\.csv$"
minimum_biological_replicates <- 3
maximum_biological_replicates <- 50
number_of_simulations <- 1000
dunnett_critical_simulations <- 100000
alpha <- 0.05
target_power <- c(0.80, 0.90)
effect_size_multipliers <- c(0.75, 1.00, 1.25)
variability_multipliers <- c(0.75, 1.00, 1.25)
random_seed <- 261009

baseline_treatment <- "Tri_Ctrl"

treatment_levels <- c(
  "Tri_Ctrl",
  "5ng_LPS",
  "100ng_LPS",
  "500ng_LPS",
  "5ug_LPS"
)

time_order <- c(
  "6hr",
  "24hr",
  "48hr",
  "72hr",
  "96hr"
)

response_columns <- c(
  Raw = "Results",
  Log = "logResults"
)


# ============================================================
# DATA LABEL HELPERS
# ============================================================

parse_time_treatment <- function(x) {
  raw <- trimws(as.character(x))
  has_separator <- grepl("_", raw, fixed = TRUE)

  time <- sub("^([^_]+)_.*$", "\\1", raw)
  treatment <- sub("^[^_]+_?", "", raw)
  time[!has_separator] <- NA_character_
  treatment[!has_separator] <- NA_character_

  treatment_upper <- toupper(trimws(treatment))
  treatment[
    treatment_upper %in% c("CTRL", "CONTROL", "TRI", "TRI_CTRL")
  ] <-
    baseline_treatment

  treatment <- trimws(treatment)
  treatment[treatment == ""] <- NA_character_

  list(
    Time = trimws(time),
    Treatment = treatment
  )
}


map_treatment_names <- function(x) {
  treatment_map <- c(
    "5ng" = "5ng_LPS",
    "100ng" = "100ng_LPS",
    "500ng" = "500ng_LPS",
    "5ug" = "5ug_LPS",
    "5ng_LPS" = "5ng_LPS",
    "100ng_LPS" = "100ng_LPS",
    "500ng_LPS" = "500ng_LPS",
    "5ug_LPS" = "5ug_LPS",
    "Tri_Ctrl" = "Tri_Ctrl"
  )

  x <- trimws(as.character(x))
  matched <- x %in% names(treatment_map)
  x[matched] <- unname(treatment_map[x[matched]])
  x
}


prepare_input_data <- function(input_file) {
  data <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (ncol(data) < 6) {
    stop(
      "Expected at least 6 columns in ",
      basename(input_file),
      "; found ",
      ncol(data),
      "."
    )
  }

  names(data)[3:6] <- c(
    "Biological_Rep",
    "Well",
    "Image",
    "Results"
  )

  if (!"Treatment" %in% names(data)) {
    stop(
      "Treatment column not found in ",
      basename(input_file),
      "."
    )
  }

  parsed <- parse_time_treatment(data$Treatment)
  data$Time <- parsed$Time
  data$Treatment <- map_treatment_names(parsed$Treatment)
  data$Results <- suppressWarnings(as.numeric(data$Results))

  needed_columns <- c(
    "Results",
    "Treatment",
    "Time",
    "Biological_Rep",
    "Well",
    "Image"
  )

  data <- data[
    complete.cases(data[, needed_columns]),
    needed_columns,
    drop = FALSE
  ]

  invalid_treatments <- setdiff(
    unique(data$Treatment),
    treatment_levels
  )
  invalid_times <- setdiff(
    unique(data$Time),
    time_order
  )

  if (length(invalid_treatments) > 0) {
    stop(
      "Unrecognized treatment labels in ",
      basename(input_file),
      ": ",
      paste(invalid_treatments, collapse = ", "),
      "."
    )
  }

  if (length(invalid_times) > 0) {
    stop(
      "Unrecognized time labels in ",
      basename(input_file),
      ": ",
      paste(invalid_times, collapse = ", "),
      "."
    )
  }

  if (nrow(data) == 0) {
    stop("No complete observations in ", basename(input_file), ".")
  }

  if (any(data$Results <= -0.5)) {
    stop(
      "All Results values must be greater than -0.5 to form logResults in ",
      basename(input_file),
      "."
    )
  }

  data$logResults <- log(data$Results + 0.5)
  data$Biological_Rep <- factor(data$Biological_Rep)
  data$Well <- factor(data$Well)
  data$Image <- factor(data$Image)
  data$Treatment <- factor(data$Treatment, levels = treatment_levels)
  data$Time <- factor(data$Time, levels = time_order)
  data$BioWell <- interaction(
    data$Biological_Rep,
    data$Well,
    drop = TRUE
  )
  data$BioWellImage <- interaction(
    data$Biological_Rep,
    data$Well,
    data$Image,
    drop = TRUE
  )

  data
}


# ============================================================
# MIXED MODEL AND DESIGN VARIANCE
# ============================================================

get_random_effects <- function(data) {
  biological_replicates <- nlevels(droplevels(data$Biological_Rep))

  if (biological_replicates < 2) {
    stop("At least 2 observed biological replicates are required.")
  }

  has_repeated_images <- any(table(data$BioWellImage) >= 2)
  has_usable_image_group <-
    nlevels(data$BioWellImage) < nrow(data) &&
    has_repeated_images

  if (has_usable_image_group) {
    return(
      "(1 | Biological_Rep) + (1 | BioWell) + (1 | BioWellImage)"
    )
  }

  "(1 | Biological_Rep) + (1 | BioWell)"
}


fit_morphology_model <- function(data, response) {
  model_data <- droplevels(data)
  random_effects <- get_random_effects(model_data)
  formula <- stats::as.formula(
    paste(response, "~ Treatment +", random_effects)
  )

  lmerTest::lmer(
    formula,
    data = model_data,
    REML = TRUE,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
}


estimate_replicate_mean_variances <- function(data, model) {
  variance_components <- as.data.frame(lme4::VarCorr(model))
  random_variances <- stats::setNames(
    variance_components$vcov[
      variance_components$grp != "Residual"
    ],
    variance_components$grp[
      variance_components$grp != "Residual"
    ]
  )
  random_variances <- random_variances[
    names(random_variances) != "Biological_Rep"
  ]
  residual_variance <- stats::sigma(model)^2

  variance_rows <- list()

  for (treatment in treatment_levels) {
    treatment_data <- data[
      data$Treatment == treatment,
      ,
      drop = FALSE
    ]

    if (nrow(treatment_data) == 0) {
      stop(
        "Treatment ",
        treatment,
        " is absent; all treatments including control are required."
      )
    }

    replicate_variances <- numeric()

    for (biological_replicate in levels(droplevels(data$Biological_Rep))) {
      replicate_data <- treatment_data[
        treatment_data$Biological_Rep == biological_replicate,
        ,
        drop = FALSE
      ]

      if (nrow(replicate_data) == 0) {
        stop(
          "The blocked power analysis requires each biological replicate ",
          "to include every treatment. Missing ",
          treatment,
          " in biological replicate ",
          biological_replicate,
          "."
        )
      }

      number_of_cells <- nrow(replicate_data)
      mean_variance <- residual_variance / number_of_cells

      for (group_name in names(random_variances)) {
        if (!group_name %in% names(replicate_data)) {
          stop(
            "Model variance group ",
            group_name,
            " is not present in the input data."
          )
        }

        group_counts <- table(replicate_data[[group_name]])
        group_weights <- as.numeric(group_counts) / number_of_cells
        mean_variance <- mean_variance +
          sum(group_weights^2) * random_variances[[group_name]]
      }

      replicate_variances <- c(replicate_variances, mean_variance)
    }

    variance_rows[[treatment]] <- data.frame(
      Treatment = treatment,
      MeanVariance = mean(replicate_variances),
      MeanSE = sqrt(mean(replicate_variances)),
      ObservedBiologicalReplicates = length(replicate_variances),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, variance_rows)
}


get_model_treatment_means <- function(model) {
  means <- as.data.frame(
    emmeans::emmeans(model, ~ Treatment)
  )

  required_columns <- c("Treatment", "emmean")
  if (!all(required_columns %in% names(means))) {
    stop("Could not extract treatment means from the fitted mixed model.")
  }

  means <- means[
    match(treatment_levels, as.character(means$Treatment)),
    ,
    drop = FALSE
  ]

  if (anyNA(means$emmean)) {
    stop("Could not estimate a mean for every treatment.")
  }

  means
}


get_dunnett_treatment_effects <- function(model) {
  treatment_emmeans <- emmeans::emmeans(model, ~ Treatment)
  comparisons <- emmeans::contrast(
    treatment_emmeans,
    method = "trt.vs.ctrl",
    ref = match(baseline_treatment, treatment_levels),
    adjust = "dunnett"
  )

  as.data.frame(summary(comparisons))
}


# ============================================================
# DUNNETT-ADJUSTED POWER SIMULATION
# ============================================================

get_dunnett_critical_value <- function(
    mean_variances,
    number_of_replicates,
    alpha
) {
  control_variance <- mean_variances[[baseline_treatment]]
  treatment_variances <- mean_variances[
    setdiff(treatment_levels, baseline_treatment)
  ]
  contrast_variances <- treatment_variances + control_variance
  if (any(!is.finite(contrast_variances)) || any(contrast_variances <= 0)) {
    stop(
      "Dunnett contrasts require finite, positive treatment-control variances."
    )
  }

  contrast_covariance <- matrix(
    control_variance,
    nrow = length(contrast_variances),
    ncol = length(contrast_variances)
  )
  diag(contrast_covariance) <- contrast_variances

  contrast_correlation <- stats::cov2cor(contrast_covariance)
  degrees_of_freedom <- number_of_replicates - 1
  correlated_normals <- matrix(
    stats::rnorm(
      dunnett_critical_simulations * nrow(contrast_correlation)
    ),
    nrow = dunnett_critical_simulations
  ) %*% chol(contrast_correlation)
  multivariate_t_draws <- correlated_normals /
    sqrt(stats::rchisq(
      dunnett_critical_simulations,
      df = degrees_of_freedom
    ) / degrees_of_freedom)
  max_absolute_t <- apply(abs(multivariate_t_draws), 1, max)

  unname(stats::quantile(
    max_absolute_t,
    probs = 1 - alpha,
    names = FALSE
  ))
}


simulate_power_curve <- function(
    treatment_means,
    mean_variances,
    replicate_numbers,
    simulations,
    alpha,
    effect_multiplier,
    variability_multiplier,
    dunnett_critical_values = NULL
) {
  lps_treatments <- setdiff(treatment_levels, baseline_treatment)
  control_mean <- treatment_means[[baseline_treatment]]
  scaled_effects <-
    (treatment_means[lps_treatments] - control_mean) * effect_multiplier
  scaled_mean_variances <-
    mean_variances * variability_multiplier^2

  results <- vector("list", length(replicate_numbers))

  for (index in seq_along(replicate_numbers)) {
    number_of_replicates <- replicate_numbers[[index]]
    critical_value <- if (is.null(dunnett_critical_values)) {
      get_dunnett_critical_value(
        scaled_mean_variances,
        number_of_replicates,
        alpha
      )
    } else {
      dunnett_critical_values[[index]]
    }
    significant_counts <- stats::setNames(
      integer(length(lps_treatments)),
      lps_treatments
    )

    replicate_means <- array(
      stats::rnorm(
        number_of_replicates *
          length(treatment_levels) *
          simulations
      ),
      dim = c(
        number_of_replicates,
        length(treatment_levels),
        simulations
      )
    )
    replicate_means <- sweep(
      replicate_means,
      2,
      sqrt(scaled_mean_variances[treatment_levels]),
      `*`
    )
    simulated_means <- treatment_means[treatment_levels]
    simulated_means[[baseline_treatment]] <- control_mean
    simulated_means[lps_treatments] <-
      control_mean + scaled_effects[lps_treatments]
    replicate_means <- sweep(
      replicate_means,
      2,
      simulated_means,
      `+`
    )

    for (treatment_index in seq_along(lps_treatments)) {
      treatment_column <- match(
        lps_treatments[[treatment_index]],
        treatment_levels
      )
      control_column <- match(baseline_treatment, treatment_levels)
      differences <- matrix(
        replicate_means[, treatment_column, ] -
          replicate_means[, control_column, ],
        nrow = number_of_replicates,
        ncol = simulations
      )
      difference_means <- colMeans(differences)
      sample_variances <- pmax(
        (
          colSums(differences^2) -
            number_of_replicates * difference_means^2
        ) / (number_of_replicates - 1),
        0
      )
      test_statistics <- ifelse(
        sample_variances > 0,
        difference_means /
          sqrt(sample_variances / number_of_replicates),
        ifelse(difference_means == 0, 0, sign(difference_means) * Inf)
      )

      significant_counts[[treatment_index]] <- sum(
        abs(test_statistics) >= critical_value
      )
    }

    estimated_power <- significant_counts / simulations

    results[[index]] <- data.frame(
      BiologicalReplicates = number_of_replicates,
      Treatment = lps_treatments,
      Power = as.numeric(estimated_power),
      MonteCarloSE = sqrt(
        as.numeric(estimated_power) *
          (1 - as.numeric(estimated_power)) /
          simulations
      ),
      DunnettCriticalValue = as.numeric(critical_value),
      Simulations = simulations,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results)
}


run_power_sensitivity <- function(
    treatment_means,
    mean_variances,
    replicate_numbers
) {
  # Multiplying every variance component by the same factor does not
  # change the shared-control contrast correlation or critical value.
  dunnett_critical_values <- vapply(
    replicate_numbers,
    function(number_of_replicates) {
      get_dunnett_critical_value(
        mean_variances,
        number_of_replicates,
        alpha
      )
    },
    numeric(1)
  )

  scenario_results <- list()
  result_index <- 1

  for (effect_multiplier in effect_size_multipliers) {
    for (variability_multiplier in variability_multipliers) {
      scenario_results[[result_index]] <- simulate_power_curve(
        treatment_means = treatment_means,
        mean_variances = mean_variances,
        replicate_numbers = replicate_numbers,
        simulations = number_of_simulations,
        alpha = alpha,
        effect_multiplier = effect_multiplier,
        variability_multiplier = variability_multiplier,
        dunnett_critical_values = dunnett_critical_values
      )
      scenario_results[[result_index]]$EffectSizeMultiplier <-
        effect_multiplier
      scenario_results[[result_index]]$VariabilityMultiplier <-
        variability_multiplier
      result_index <- result_index + 1
    }
  }

  do.call(rbind, scenario_results)
}


summarize_minimum_replicates <- function(power_results) {
  scenario_keys <- unique(
    power_results[
      ,
      c(
        "Treatment",
        "EffectSizeMultiplier",
        "VariabilityMultiplier"
      ),
      drop = FALSE
    ]
  )

  summaries <- vector("list", nrow(scenario_keys))

  for (index in seq_len(nrow(scenario_keys))) {
    scenario <- scenario_keys[index, ]
    matching_rows <- power_results[
      power_results$Treatment == scenario$Treatment &
        power_results$EffectSizeMultiplier ==
          scenario$EffectSizeMultiplier &
        power_results$VariabilityMultiplier ==
          scenario$VariabilityMultiplier,
      ,
      drop = FALSE
    ]

    summary_row <- data.frame(
      Treatment = scenario$Treatment,
      EffectSizeMultiplier = scenario$EffectSizeMultiplier,
      VariabilityMultiplier = scenario$VariabilityMultiplier,
      stringsAsFactors = FALSE
    )

    for (power_target in target_power) {
      target_label <- paste0(
        "MinimumReps_",
        round(power_target * 100)
      )
      qualifying_replicates <- matching_rows$BiologicalReplicates[
        matching_rows$Power >= power_target
      ]
      summary_row[[target_label]] <- if (
        length(qualifying_replicates) > 0
      ) {
        min(qualifying_replicates)
      } else {
        NA_integer_
      }
    }

    summaries[[index]] <- summary_row
  }

  do.call(rbind, summaries)
}


# ============================================================
# POWER OUTPUTS
# ============================================================

save_power_outputs <- function(
    data,
    model,
    response_name,
    time_label,
    output_dir
) {
  treatment_means_table <- get_model_treatment_means(model)
  treatment_means <- stats::setNames(
    treatment_means_table$emmean,
    as.character(treatment_means_table$Treatment)
  )
  design_variances <- estimate_replicate_mean_variances(data, model)
  treatment_effects <- get_dunnett_treatment_effects(model)
  mean_variances <- stats::setNames(
    design_variances$MeanVariance,
    design_variances$Treatment
  )
  replicate_numbers <- seq.int(
    minimum_biological_replicates,
    maximum_biological_replicates
  )

  power_results <- run_power_sensitivity(
    treatment_means = treatment_means,
    mean_variances = mean_variances,
    replicate_numbers = replicate_numbers
  )
  minimum_replicates <- summarize_minimum_replicates(power_results)
  baseline_minimum_replicates <- minimum_replicates[
    minimum_replicates$EffectSizeMultiplier == 1 &
      minimum_replicates$VariabilityMultiplier == 1,
    c(
      "Treatment",
      paste0("MinimumReps_", round(target_power * 100))
    ),
    drop = FALSE
  ]

  write.csv(
    power_results,
    file.path(
      output_dir,
      paste0("PowerCurves_", response_name, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    baseline_minimum_replicates,
    file.path(
      output_dir,
      paste0("MinimumBiologicalReplicates_", response_name, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    minimum_replicates,
    file.path(
      output_dir,
      paste0(
        "MinimumBiologicalReplicates_Sensitivity_",
        response_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  write.csv(
    design_variances,
    file.path(
      output_dir,
      paste0("ReplicateMeanVariance_", response_name, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    treatment_means_table,
    file.path(
      output_dir,
      paste0("EstimatedTreatmentMeans_", response_name, ".csv")
    ),
    row.names = FALSE
  )
  write.csv(
    treatment_effects,
    file.path(
      output_dir,
      paste0("DunnettTreatmentEffects_", response_name, ".csv")
    ),
    row.names = FALSE
  )

  baseline_results <- power_results[
    power_results$EffectSizeMultiplier == 1 &
      power_results$VariabilityMultiplier == 1,
    ,
    drop = FALSE
  ]

  baseline_plot <- ggplot2::ggplot(
    baseline_results,
    ggplot2::aes(
      x = BiologicalReplicates,
      y = Power,
      color = Treatment,
      group = Treatment
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_hline(
      yintercept = target_power,
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1)
    ) +
    ggplot2::scale_x_continuous(
      breaks = replicate_numbers
    ) +
    ggplot2::labs(
      title = paste(
        "Dunnett-adjusted power:",
        time_label,
        response_name
      ),
      subtitle = paste0(
        "Observed effects and variance; family-wise alpha = ",
        alpha
      ),
      x = "Independent biological replicates",
      y = "Power",
      color = "Treatment"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )

  ggplot2::ggsave(
    file.path(
      output_dir,
      paste0("PowerCurve_", response_name, ".png")
    ),
    baseline_plot,
    width = 8,
    height = 5,
    dpi = 300
  )

  sensitivity_plot <- ggplot2::ggplot(
    power_results,
    ggplot2::aes(
      x = BiologicalReplicates,
      y = Power,
      color = factor(EffectSizeMultiplier),
      linetype = factor(VariabilityMultiplier),
      group = interaction(
        EffectSizeMultiplier,
        VariabilityMultiplier
      )
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(
      yintercept = target_power,
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::facet_wrap(~ Treatment) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1)
    ) +
    ggplot2::scale_x_continuous(
      breaks = replicate_numbers
    ) +
    ggplot2::labs(
      title = paste(
        "Power sensitivity:",
        time_label,
        response_name
      ),
      subtitle = paste(
        "Effect-size multiplier shown by color;",
        "variability multiplier shown by line type"
      ),
      x = "Independent biological replicates",
      y = "Power",
      color = "Effect multiplier",
      linetype = "Variability multiplier"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )

  ggplot2::ggsave(
    file.path(
      output_dir,
      paste0("PowerSensitivity_", response_name, ".png")
    ),
    sensitivity_plot,
    width = 10,
    height = 7,
    dpi = 300
  )

  cat(
    "\nMinimum replicate counts for ",
    time_label,
    " (",
    response_name,
    "; effect and variability multipliers = 1):\n",
    sep = ""
  )
  print(baseline_minimum_replicates, row.names = FALSE)
}


write_method_note <- function(data, time_label, output_dir) {
  rows_per_image <- table(data$BioWellImage)
  random_effects <- get_random_effects(data)
  note <- c(
    "Microglial morphology power analysis",
    paste0("Time point: ", time_label),
    paste0("Input observations (Results rows): ", nrow(data)),
    paste0(
      "Independent biological replicates observed: ",
      nlevels(droplevels(data$Biological_Rep))
    ),
    paste0(
      "Rows per biological-replicate/well/image: min=",
      min(rows_per_image),
      ", median=",
      stats::median(as.numeric(rows_per_image)),
      ", max=",
      max(rows_per_image)
    ),
    paste0(
      "Mixed-model random effects: ",
      random_effects
    ),
    paste0(
      "Power uses a blocked biological-replicate treatment-mean ",
      "simulation with single-step Dunnett adjustment (family-wise ",
      "alpha = ",
      alpha,
      ")."
    ),
    paste0(
      "The row-level residual is included in each treatment-mean ",
      "variance. A distinct image variance is estimated only when ",
      "there are repeated observations within images."
    )
  )

  if (time_label == "Combined_All_Timepoints") {
    note <- c(
      note,
      paste0(
        "POOLED ANALYSIS: All observed exposure times are combined by ",
        "concentration. Time is omitted from the fixed-effects model; ",
        "the effect is an observation-weighted pooled effect, not an ",
        "effect guaranteed at each time point."
      ),
      paste0(
        "Biological_Rep IDs are retained across times. The same culture ",
        "observed at multiple times is one independent replicate. This ",
        "requires consistent biological-replicate IDs across times."
      ),
      paste0(
        "Well and image groups are nested within Time and Biological_Rep ",
        "to prevent reused technical IDs at different times from merging. ",
        "This assumes distinct wells at each exposure time."
      ),
      paste0(
        "Future replicates are assumed to have the observed mix of ",
        "times and technical sampling. Time-related variation is not ",
        "separately modeled and can enter the residual/well variance."
      )
    )
  }

  if (all(rows_per_image == 1)) {
    note <- c(
      note,
      paste0(
        "LIMITATION: There is one Results row per image. Image-to-image ",
        "variation is therefore part of the row-level residual; separate ",
        "image-level and individual-microglia variability cannot be ",
        "identified from this dataset."
      )
    )
  } else {
    note <- c(
      note,
      paste0(
        "Interpret the row-level residual according to what each input ",
        "row represents; it is individual-microglia variability only ",
        "when each row is an individual microglia measurement."
      )
    )
  }

  writeLines(
    note,
    file.path(output_dir, "PowerAnalysis_Methods.txt")
  )
}


# ============================================================
# RUN ANALYSIS
# ============================================================

run_analysis_group <- function(time_data, time_label, output_dir) {
  observed_treatments <- unique(as.character(time_data$Treatment))
  missing_treatments <- setdiff(treatment_levels, observed_treatments)

  if (length(missing_treatments) > 0) {
    stop(
      "Analysis group ",
      time_label,
      " is missing treatments: ",
      paste(missing_treatments, collapse = ", "),
      "."
    )
  }

  time_output_dir <- file.path(output_dir, time_label)
  dir.create(
    time_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  write_method_note(time_data, time_label, time_output_dir)

  for (response_name in names(response_columns)) {
    response <- unname(response_columns[[response_name]])
    model <- fit_morphology_model(time_data, response)

    write.csv(
      as.data.frame(lme4::VarCorr(model)),
      file.path(
        time_output_dir,
        paste0("VarianceComponents_", response_name, ".csv")
      ),
      row.names = FALSE
    )

    save_power_outputs(
      data = time_data,
      model = model,
      response_name = response_name,
      time_label = time_label,
      output_dir = time_output_dir
    )
  }
}


run_one_timepoint <- function(data, time_label, output_dir) {
  time_data <- droplevels(data[data$Time == time_label, , drop = FALSE])

  if (nrow(time_data) == 0) {
    stop("No observations for time point ", time_label, ".")
  }

  run_analysis_group(time_data, time_label, output_dir)
}


run_combined_timepoints <- function(data, output_dir) {
  pooled_data <- droplevels(data)
  pooled_data$BioWell <- interaction(
    pooled_data$Biological_Rep,
    pooled_data$Time,
    pooled_data$Well,
    drop = TRUE
  )
  pooled_data$BioWellImage <- interaction(
    pooled_data$Biological_Rep,
    pooled_data$Time,
    pooled_data$Well,
    pooled_data$Image,
    drop = TRUE
  )

  cat("  Combining all time points by concentration\n")
  run_analysis_group(
    pooled_data,
    "Combined_All_Timepoints",
    output_dir
  )
}


run_one_csv <- function(input_file) {
  cat("\nProcessing: ", basename(input_file), "\n", sep = "")

  data <- prepare_input_data(input_file)
  observed_timepoints <- time_order[time_order %in% as.character(data$Time)]

  if (length(observed_timepoints) == 0) {
    stop("No recognized time points in ", basename(input_file), ".")
  }

  output_dir <- file.path(
    dirname(input_file),
    paste0(
      "Power_Analysis_",
      tools::file_path_sans_ext(basename(input_file))
    )
  )
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (time_label in observed_timepoints) {
    cat("  Time point: ", time_label, "\n", sep = "")
    run_one_timepoint(data, time_label, output_dir)
  }

  run_combined_timepoints(data, output_dir)

  cat("  Results saved to: ", output_dir, "\n", sep = "")
  invisible(output_dir)
}


generate_summary_tables <- function(output_dirs) {
  summary_rows <- list()
  power_columns <- paste0("MinimumReps_", round(target_power * 100))
  lps_treatments <- setdiff(treatment_levels, baseline_treatment)

  for (output_dir in output_dirs) {
    morphology <- sub("^Power_Analysis_", "", basename(output_dir))
    groups <- c(time_order, "Combined_All_Timepoints")
    groups <- groups[dir.exists(file.path(output_dir, groups))]

    if (length(groups) == 0) {
      stop("No analysis groups found in: ", output_dir)
    }

    for (group in groups) {
      for (response_name in names(response_columns)) {
        summary_file <- file.path(
          output_dir,
          group,
          paste0("MinimumBiologicalReplicates_", response_name, ".csv")
        )
        if (!file.exists(summary_file)) {
          stop("Missing minimum-replicate results: ", summary_file)
        }

        results <- read.csv(summary_file, stringsAsFactors = FALSE)
        if (
          !all(c("Treatment", power_columns) %in% names(results)) ||
          nrow(results) != length(lps_treatments) ||
          anyDuplicated(results$Treatment) ||
          !setequal(results$Treatment, lps_treatments)
        ) {
          stop("Expected one baseline result per concentration in: ", summary_file)
        }

        results <- results[
          match(lps_treatments, results$Treatment),
          c("Treatment", power_columns),
          drop = FALSE
        ]
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          Morphology = morphology,
          Analysis = if (group == "Combined_All_Timepoints") {
            "Pooled across times"
          } else {
            "Time-specific"
          },
          Time = group,
          ResponseScale = response_name,
          results,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    }
  }

  if (length(summary_rows) == 0) {
    stop("No minimum-replicate results available for summary generation.")
  }

  master_summary <- do.call(rbind, summary_rows)
  summary_dir <- file.path(input_dir, "Power_Analysis_Summaries")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    master_summary,
    file.path(summary_dir, "All_Morphology_MinimumBiologicalReplicates.csv"),
    row.names = FALSE
  )

  for (morphology in unique(master_summary$Morphology)) {
    write.csv(
      master_summary[master_summary$Morphology == morphology, , drop = FALSE],
      file.path(
        summary_dir,
        paste0(morphology, "_MinimumBiologicalReplicates.csv")
      ),
      row.names = FALSE
    )
  }

  writeLines(
    c(
      "Minimum biological replicates: observed effects and variability.",
      "One row per morphology, analysis/time, response scale, and concentration.",
      "MinimumReps_80 and MinimumReps_90 are the estimated counts for 80% and 90% power.",
      paste0(
        "NA means the target was not reached within ",
        minimum_biological_replicates, "-", maximum_biological_replicates,
        " biological replicates; it does not mean zero or missing observations."
      ),
      "Pooled results assume the observed mix of exposure times and are not guarantees at each time.",
      "Raw and Log results are separate analyses; sensitivity scenarios are excluded.",
      "Consult each analysis group's methods and variance outputs for model limitations."
    ),
    file.path(summary_dir, "Summary_Interpretation.txt")
  )
  cat("\nSummary tables saved to: ", summary_dir, "\n", sep = "")
  invisible(master_summary)
}


main <- function() {
  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }

  csv_files <- list.files(
    input_dir,
    pattern = file_pattern,
    full.names = TRUE
  )
  csv_files <- csv_files[
    !grepl(
      "Power_Analysis_",
      basename(csv_files),
      fixed = TRUE
    )
  ]

  if (length(csv_files) == 0) {
    stop("No CSV files found in: ", input_dir)
  }

  set.seed(random_seed)
  cat("Found ", length(csv_files), " CSV file(s).\n", sep = "")

  output_dirs <- character(length(csv_files))
  for (index in seq_along(csv_files)) {
    output_dirs[[index]] <- run_one_csv(csv_files[[index]])
  }

  generate_summary_tables(output_dirs)

  cat("\nPower analyses completed.\n")
}


main()

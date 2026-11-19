# ============================================================
# MIXED MODELS WITH BASELINE COMPARISON
# For CSV with: Filename, Treatment, Threshold_Method, Image_Group, 
#               Biological_Rep, Well, Image, Results
# Includes baseline as reference level for comparisons
# EDIT: all methods are plotted together against Baseline
#       (one plot per treatment) instead of one plot per method
# ============================================================

pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr","emmeans","pbkrtest","ggsignif")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# ============================================================
# UPDATE THIS PATH
# ============================================================
input_file <- '/Users/allisonpickle/Desktop/Compiled Thresholds_Sheet1_RMixedFormated (1).csv'
setwd(dirname(input_file))

# ============================================================
# DATA LOADING
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("LOADING DATA\n")
cat(rep("=", 80), "\n\n")

# Read data
data <- read.csv(input_file, stringsAsFactors = FALSE)

cat(sprintf("✓ Loaded file: %s\n", basename(input_file)))
cat(sprintf("  Total rows: %d\n", nrow(data)))
cat(sprintf("  Columns: %s\n\n", paste(names(data), collapse = ", ")))

# ============================================================
# DEBUG: CHECK DATA STRUCTURE
# ============================================================

cat(rep("=", 80), "\n")
cat("DATA STRUCTURE CHECK\n")
cat(rep("=", 80), "\n\n")

cat("First 10 rows:\n")
print(head(data, 10))

cat("\n\nColumn types:\n")
print(sapply(data, class))

cat("\n\nNA counts per column:\n")
print(colSums(is.na(data)))

cat("\n\nUnique values per column:\n")
for (col in names(data)) {
  cat(sprintf("%s: %d unique values\n", col, length(unique(data[[col]]))))
  if (col %in% c("Treatment", "Threshold_Method", "Biological_Rep", "Well", "Image")) {
    cat(sprintf("  Values: %s\n", paste(head(unique(data[[col]]), 10), collapse = ", ")))
  }
}

# ============================================================
# DATA PREPARATION
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("DATA PREPARATION\n")
cat(rep("=", 80), "\n\n")

treatment_labels <- c("Ctrl" = "-\n-", "5ugLPS" = "-\n+\n")
row_labels <- c("Control", "5ug/mL LPS")

# Check if treatment names match
cat("=== TREATMENT MATCHING ===\n")
cat("Expected treatments:", paste(names(treatment_labels), collapse = ", "), "\n")
cat("Actual treatments:", paste(unique(data$Treatment), collapse = ", "), "\n\n")

# Detect all threshold methods present in data
detected_methods <- unique(data$Threshold_Method)
detected_methods <- detected_methods[!is.na(detected_methods)]

cat("=== THRESHOLD METHODS DETECTED ===\n")
cat("Methods found in data:", paste(detected_methods, collapse = ", "), "\n")

# Check if BASELINE exists in data
has_baseline <- "BASELINE" %in% toupper(detected_methods)

if (has_baseline) {
  # Find the exact case used for baseline
  baseline_name <- detected_methods[toupper(detected_methods) == "BASELINE"]
  cat(sprintf("\n✓ BASELINE method detected: '%s'\n", baseline_name))
  
  # Put BASELINE first, then sort remaining methods alphabetically
  other_methods <- sort(detected_methods[toupper(detected_methods) != "BASELINE"])
  all_methods <- c(baseline_name, other_methods)
  cat("Factor levels (BASELINE as reference):", paste(all_methods, collapse = ", "), "\n\n")
} else {
  cat("\n⚠️  No BASELINE method detected in data\n")
  cat("   Threshold methods will be compared to each other only\n")
  all_methods <- sort(detected_methods)
  cat("Factor levels (alphabetical):", paste(all_methods, collapse = ", "), "\n\n")
}

# Factor conversion
data$Treatment <- factor(data$Treatment, levels = names(treatment_labels))
data$Threshold_Method <- factor(data$Threshold_Method, levels = all_methods)
data$Well <- factor(data$Well)
data$Biological_Rep <- factor(data$Biological_Rep)
if ("Image" %in% names(data)) {
  data$Image <- factor(data$Image)
}

# Check for NA Results after factorization
cat("=== CHECKING FOR ISSUES ===\n")
cat(sprintf("NA in Treatment: %d\n", sum(is.na(data$Treatment))))
cat(sprintf("NA in Threshold_Method: %d\n", sum(is.na(data$Threshold_Method))))
cat(sprintf("NA in Results: %d\n", sum(is.na(data$Results))))
cat(sprintf("Valid rows (non-NA in all critical columns): %d\n", 
            sum(!is.na(data$Treatment) & !is.na(data$Results) & 
                !is.na(data$Threshold_Method) & !is.na(data$Biological_Rep))))

# Remove rows with NA in critical columns
data <- data[!is.na(data$Treatment) & 
             !is.na(data$Results) & 
             !is.na(data$Threshold_Method) &
             !is.na(data$Biological_Rep), ]

cat(sprintf("\nRows after removing NAs: %d\n", nrow(data)))

if (nrow(data) == 0) {
  stop("ERROR: No valid data remaining after removing NAs. Check your CSV file format and column names.")
}

# Print data structure
cat("\n=== DATA STRUCTURE ===\n")
cat("Treatment counts:\n")
print(table(data$Treatment, useNA = "ifany"))
cat("\nThreshold method counts:\n")
print(table(data$Threshold_Method, useNA = "ifany"))
cat(sprintf("\nBiological replicates: %s\n", 
            paste(sort(unique(data$Biological_Rep)), collapse = ", ")))

# Check biological replicates
n_bioreps <- length(unique(data$Biological_Rep))
cat(sprintf("\n=== BIOLOGICAL REPLICATE CHECK ===\n"))
cat(sprintf("Number: %d\n", n_bioreps))

if (n_bioreps < 3) {
  cat("\n⚠️  WARNING: Only", n_bioreps, "biological replicate(s) detected ⚠️\n")
  cat("Simplified model structure will be used. Results are not suitable for biological inference.\n\n")
} else {
  cat(sprintf("\n✓ Sufficient biological replicates (%d ≥ 3)\n\n", n_bioreps))
}

# Log transformation
data$logResults <- log(data$Results + 0.5)

# ============================================================
# BUILD MODELS
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("BUILDING INTERACTION MODEL: Treatment × Threshold Method\n")
cat(sprintf("Methods included: %s\n", paste(all_methods, collapse = ", ")))
if (has_baseline) {
  cat(sprintf("Reference level: %s (all methods compared to this)\n", all_methods[1]))
}
cat(rep("=", 80), "\n\n")

# Determine random effects structure
has_image <- "Image" %in% names(data)
if (n_bioreps >= 3 && has_image) {
  random_effects <- "(1 | Biological_Rep) + (1 | Biological_Rep:Well) + (1 | Well:Image)"
  model_type <- "Full nested with Treatment × Threshold interaction"
} else if (n_bioreps >= 3) {
  random_effects <- "(1 | Biological_Rep) + (1 | Biological_Rep:Well)"
  model_type <- "Nested with Treatment × Threshold interaction"
} else if (has_image) {
  random_effects <- "(1 | Well) + (1 | Well:Image)"
  model_type <- "⚠️ SIMPLIFIED with Treatment × Threshold interaction"
} else {
  random_effects <- "(1 | Well)"
  model_type <- "⚠️ SIMPLIFIED with Treatment × Threshold interaction"
}

cat(sprintf("Model type: %s\n", model_type))
cat(sprintf("Random effects: %s\n\n", random_effects))

# Log-transformed model
log_model_formula <- as.formula(paste("logResults ~ Treatment * Threshold_Method +", random_effects))
LogInteractionModel <- lmerTest::lmer(log_model_formula, data = data)

cat("Formula: "); print(log_model_formula)
cat("\nANOVA Table:\n"); print(anova(LogInteractionModel))

# ============================================================
# POST-HOC TESTS
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("POST-HOC COMPARISONS\n")
cat(rep("=", 80), "\n\n")

# 1. Overall treatment effect comparison
cat("=== Overall Treatment Effect ===\n")
emm_treatment <- emmeans(LogInteractionModel, 
                        pairwise ~ Treatment,
                        lmer.df = "satterthwaite")
print(emm_treatment$contrasts)

# 2. Treatment effect within each threshold method
cat("\n=== Treatment Effect | Threshold Method ===\n")
emm_treatment_by_method <- emmeans(LogInteractionModel,
                                   pairwise ~ Treatment | Threshold_Method,
                                   lmer.df = "satterthwaite")
print(emm_treatment_by_method$contrasts)

emm_pairs_df <- as.data.frame(emm_treatment_by_method$contrasts)

# 3. Compare threshold methods within each treatment
cat("\n=== Threshold Method Comparisons | Treatment ===\n")
emm_method_by_treatment <- emmeans(LogInteractionModel,
                                   pairwise ~ Threshold_Method | Treatment,
                                   lmer.df = "satterthwaite")
print(summary(emm_method_by_treatment$contrasts))

# Store for later use in plots
method_contrasts_df <- as.data.frame(summary(emm_method_by_treatment$contrasts))

# 4. Baseline-specific comparisons (if baseline exists)
if (has_baseline) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("BASELINE COMPARISONS (Dunnett-style contrasts)\n")
  cat(rep("=", 80), "\n\n")
  
  baseline_name <- all_methods[1]
  
  # Compare each method to baseline within each treatment
  cat(sprintf("=== Each Method vs %s | Treatment ===\n", baseline_name))
  emm_vs_baseline <- emmeans(LogInteractionModel,
                             trt.vs.ctrl ~ Threshold_Method | Treatment,
                             ref = 1,  # First level (baseline) as reference
                             lmer.df = "satterthwaite")
  baseline_contrasts <- summary(emm_vs_baseline$contrasts)
  print(baseline_contrasts)
  
  # Overall comparison (collapsing across treatments)
  cat(sprintf("\n=== Each Method vs %s (Overall) ===\n", baseline_name))
  emm_vs_baseline_overall <- emmeans(LogInteractionModel,
                                     trt.vs.ctrl ~ Threshold_Method,
                                     ref = 1,
                                     lmer.df = "satterthwaite")
  baseline_overall <- summary(emm_vs_baseline_overall$contrasts)
  print(baseline_overall)
}

# ============================================================
# VISUALIZATION
# ============================================================

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base, if(has_baseline) "_with_Baseline" else "")
dir.create(out_dir, showWarnings = FALSE)

cat("\n", rep("=", 80), "\n", sep = "")
cat("GENERATING VISUALIZATIONS\n")
cat(rep("=", 80), "\n\n")

# Create dynamic threshold labels - capitalize first letter (except BASELINE stays uppercase)
threshold_labels <- setNames(
  sapply(all_methods, function(x) {
    if (toupper(x) == "BASELINE") {
      "Baseline"
    } else {
      paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
    }
  }),
  all_methods
)

cat("=== PLOT LABELS ===\n")
print(threshold_labels)
cat("\n")

# Interaction plot
p_threshold <- emmip(LogInteractionModel, Treatment ~ Threshold_Method, CIs = TRUE) +
  scale_color_manual(values = c("#444444", "#FF6B6B"), 
                     labels = c("Control", "LPS 5µg/mL")) +
  scale_x_discrete(labels = threshold_labels) +
  labs(title = if(has_baseline) "Treatment Effect: Baseline vs Threshold Methods" else "Treatment Effect Across Threshold Methods",
       x = "Threshold Method", 
       y = "log(Results + 0.5)", 
       color = "Treatment") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "Treatment_Interaction_Plot.png"), 
       p_threshold, width = 10, height = 6, dpi = 300)

cat("✓ Saved: Treatment_Interaction_Plot.png\n")

# Bar plots for each method
get_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

create_barplot <- function(subset_data, p_val, title_text) {
  summ_df <- subset_data |>
    group_by(Treatment) |>
    summarise(mean = mean(logResults, na.rm = TRUE),
              se = sd(logResults, na.rm = TRUE) / sqrt(sum(!is.na(logResults))),
              .groups = "drop")
  
  stars <- if(length(p_val) > 0 && !is.na(p_val)) get_stars(p_val) else "ns"
  y_max <- max(summ_df$mean + summ_df$se)
  y_annot <- y_max + 0.15 * (max(subset_data$logResults) - min(subset_data$logResults))
  
  aligned_treatment_labels <- setNames(paste0("\n", treatment_labels), names(treatment_labels))
  combined_row_label <- paste0("\n", paste(row_labels, collapse = "\n"))
  
  p <- ggplot(summ_df, aes(x = Treatment, y = mean, fill = Treatment)) +
    geom_col(alpha = 0.7) +
    geom_jitter(data = subset_data, 
                aes(x = Treatment, y = logResults, 
                    shape = Biological_Rep, color = Biological_Rep),
                width = 0.2, alpha = 0.6, size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_fill_manual(values = c("#444444", "#FF6B6B")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 7)) +
    scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00", 
                                  "black", "black", "black", "black", "black")) +
    scale_x_discrete(labels = aligned_treatment_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(title = title_text, x = "", y = "log(Results + 0.5) (mean ± SE)",
         shape = "Biological Rep", color = "Biological Rep") +
    annotate("text", x = 0.4, y = -Inf, label = combined_row_label, 
             hjust = 1, vjust = 1.1, size = 3, color = "black", lineheight = 0.8) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.margin = margin(5, 5, 15, 50)
    ) +
    coord_cartesian(clip = "off")
  
  if (stars != "ns") {
    p <- p +
      geom_segment(aes(x = 1, xend = 2, y = y_annot, yend = y_annot), 
                   color = "black", linewidth = 0.5) +
      geom_segment(aes(x = 1, xend = 1, y = y_annot, yend = y_annot - 0.02), 
                   color = "black", linewidth = 0.5) +
      geom_segment(aes(x = 2, xend = 2, y = y_annot, yend = y_annot - 0.02), 
                   color = "black", linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_annot + 0.05, label = stars, 
               size = 7, fontface = "bold")
  }
  return(p)
}

# Bar plots for each threshold method
for (method in levels(data$Threshold_Method)) {
  subset_data <- data[data$Threshold_Method == method, ]
  p_val <- emm_pairs_df$p.value[emm_pairs_df$Threshold_Method == method]
  method_label <- threshold_labels[method]
  p <- create_barplot(subset_data, p_val, paste("Method:", method_label))
  ggsave(file.path(out_dir, paste0("Barplot_", method, ".png")),
         p, width = 6, height = 4, dpi = 300)
  cat(sprintf("✓ Saved: Barplot_%s.png\n", method))
}

# ============================================================
# EDITED: BASELINE vs ALL THRESHOLD METHODS (ONE PLOT PER TREATMENT)
# Order: Baseline, Otsu, Li, Yen, ISODATA, Triangle, Mean, Minimum Error
# ============================================================

if (has_baseline) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("GENERATING COMBINED BASELINE COMPARISON PLOTS\n")
  cat(rep("=", 80), "\n\n")
  
  baseline_name <- all_methods[1]
  baseline_dir <- file.path(out_dir, "Baseline_Comparisons")
  dir.create(baseline_dir, showWarnings = FALSE)
  
  # ---- Desired plotting order (display names) ----
  desired_order <- c("Otsu", "Li", "Yen", "ISODATA", "Triangle", "Mean", "Minimum Error")
  norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
  
  # Match desired names to the method names actually in the data (case/punctuation-insensitive)
  matched <- sapply(desired_order, function(d) {
    hit <- other_methods[norm_name(other_methods) == norm_name(d)]
    if (length(hit) == 0) NA_character_ else hit[1]
  })
  
  missing_desired <- desired_order[is.na(matched)]
  if (length(missing_desired) > 0) {
    cat("  ⚠️  Not found in data (skipped):", paste(missing_desired, collapse = ", "), "\n")
  }
  
  ordered_methods <- unname(matched[!is.na(matched)])
  extra_methods   <- setdiff(other_methods, ordered_methods)  # any methods not in the list go last
  if (length(extra_methods) > 0) {
    cat("  ℹ️  Extra methods appended at end:", paste(extra_methods, collapse = ", "), "\n")
  }
  plot_order <- c(baseline_name, ordered_methods, extra_methods)
  
  # Display labels in plot order
  display_labels <- setNames(threshold_labels[plot_order], plot_order)
  for (i in seq_along(desired_order)) {
    if (!is.na(matched[i])) display_labels[matched[i]] <- desired_order[i]
  }
  
  create_all_methods_vs_baseline_plot <- function(treatment_group, contrast_data) {
    subset_data <- data[data$Treatment == treatment_group &
                        data$Threshold_Method %in% plot_order, ]
    if (nrow(subset_data) == 0) {
      cat(sprintf("  ⚠️  No data for %s\n", treatment_group))
      return(NULL)
    }
    
    subset_data$Threshold_Method <- factor(subset_data$Threshold_Method, levels = plot_order)
    
    summ_df <- subset_data |>
      group_by(Threshold_Method) |>
      summarise(mean = mean(logResults, na.rm = TRUE),
                se = sd(logResults, na.rm = TRUE) / sqrt(sum(!is.na(logResults))),
                .groups = "drop")
    
    # Clean contrast labels (emmeans may wrap names containing spaces in parentheses)
    contrast_clean <- gsub("[()]", "", as.character(contrast_data$contrast))
    
    # Stars for each method vs baseline
    summ_df$stars <- sapply(as.character(summ_df$Threshold_Method), function(m) {
      if (m == baseline_name) return("")
      target <- paste0(m, " - ", baseline_name)
      p <- contrast_data$p.value[contrast_data$Treatment == treatment_group &
                                 contrast_clean == target]
      if (length(p) == 0 || is.na(p[1])) return("ns")
      get_stars(p[1])
    })
    
    # Star position: above the highest of bar+SE or data points for each method
    pt_max <- subset_data |>
      group_by(Threshold_Method) |>
      summarise(pt_max = max(logResults, na.rm = TRUE), .groups = "drop")
    summ_df <- left_join(summ_df, pt_max, by = "Threshold_Method")
    y_top <- max(pmax(summ_df$mean + summ_df$se, summ_df$pt_max))
    y_bot <- min(c(0, subset_data$logResults))
    y_range <- y_top - y_bot
    summ_df$y_star <- pmax(summ_df$mean + summ_df$se, summ_df$pt_max) + 0.04 * y_range
    
    treatment_label <- ifelse(treatment_group == "Ctrl", "Control", "LPS 5µg/mL")
    bar_color <- ifelse(treatment_group == "Ctrl", "#444444", "#FF6B6B")
    
    p <- ggplot(summ_df, aes(x = Threshold_Method, y = mean)) +
      geom_col(fill = bar_color, alpha = 0.7, width = 0.7) +
      geom_jitter(data = subset_data,
                  aes(x = Threshold_Method, y = logResults,
                      shape = Biological_Rep, color = Biological_Rep),
                  width = 0.15, alpha = 0.6, size = 1.8) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
      geom_text(data = summ_df[summ_df$stars != "", ],
                aes(x = Threshold_Method, y = y_star, label = stars),
                size = 5, fontface = "bold", inherit.aes = FALSE) +
      scale_shape_manual(values = c(16, 17, 15, 18, 3, 4, 8, 7)) +
      scale_color_manual(values = c("#08267d", "#228B22", "#FF8C00",
                                    "black", "black", "black", "black", "black")) +
      scale_x_discrete(labels = display_labels) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
      labs(title = sprintf("%s: All Methods vs Baseline", treatment_label),
           x = "Method",
           y = "log(Results + 0.5) (mean ± SE)",
           shape = "Biological Rep",
           color = "Biological Rep") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    return(p)
  }
  
  for (grp in c("Ctrl", "5ugLPS")) {
    p_all <- create_all_methods_vs_baseline_plot(grp, baseline_contrasts)
    if (!is.null(p_all)) {
      fname <- sprintf("Baseline_vs_AllMethods_%s.png", ifelse(grp == "Ctrl", "Control", "LPS"))
      ggsave(file.path(baseline_dir, fname), p_all, width = 10, height = 6, dpi = 300)
      cat(sprintf("✓ Saved: %s\n", fname))
    }
  }
  
  cat(sprintf("\n✓ Combined baseline comparison plots saved to: %s/\n", basename(baseline_dir)))
}

# ============================================================
# SAVE SUMMARY
# ============================================================

cat("\n")
sink(file.path(out_dir, "Analysis_Summary.txt"))
cat("================================================================\n")
cat("THRESHOLD METHODS COMPARISON\n")
if (has_baseline) {
  cat(sprintf("Reference method: %s\n", all_methods[1]))
}
cat(sprintf("Methods analyzed: %s\n", paste(all_methods, collapse = ", ")))
cat(sprintf("Data file: %s\n", basename(input_file)))
cat("================================================================\n\n")

cat("=== TREATMENT EFFECT BY METHOD ===\n")
print(emm_pairs_df)

cat("\n\n=== THRESHOLD METHOD COMPARISONS BY TREATMENT ===\n")
print(summary(emm_method_by_treatment$contrasts))

if (has_baseline) {
  cat("\n\n", rep("=", 64), "\n", sep = "")
  cat(sprintf("BASELINE COMPARISONS (vs %s)\n", all_methods[1]))
  cat(rep("=", 64), "\n\n")
  
  cat("=== By Treatment ===\n")
  print(baseline_contrasts)
  
  cat("\n=== Overall (Across Treatments) ===\n")
  print(baseline_overall)
  
  # EDITED: describes the combined plots
  cat(sprintf("\n\n=== PLOTS GENERATED ===\n"))
  cat("Combined plots (all methods vs Baseline):\n")
  cat("  - Baseline_vs_AllMethods_Control.png\n")
  cat("  - Baseline_vs_AllMethods_LPS.png\n")
}

cat("\n\n=== MODEL SUMMARY ===\n")
print(summary(LogInteractionModel))

cat("\n\n=== ANOVA ===\n")
print(anova(LogInteractionModel))

sink()

cat("\n", rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n")
cat(sprintf("Methods analyzed: %s\n", paste(all_methods, collapse = ", ")))
if (has_baseline) {
  cat(sprintf("✓ Baseline comparisons included (reference: %s)\n", all_methods[1]))
  # EDITED: combined plots instead of per-method plots
  cat(sprintf("✓ Combined comparison plots: %s/Baseline_Comparisons/\n", out_dir))
  cat("  Total plots: 2 (Control and LPS)\n")
}
cat(sprintf("Outputs saved to: %s/\n", out_dir))
cat(rep("=", 80), "\n\n")
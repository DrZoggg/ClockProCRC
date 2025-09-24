

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(tidyr)

#' @title Main Differential Module Analysis Pipeline
#' @description Orchestrates the complete differential module eigengene analysis workflow

execute_differential_module_pipeline <- function(seurat_object,
                                                group_1_identifiers,
                                                group_2_identifiers,
                                                analysis_label = "differential_analysis",
                                                statistical_methodology = "wilcox",
                                                use_harmonized_eigengenes = TRUE,
                                                fold_change_threshold = 0,
                                                minimum_percentage = 0.1,
                                                generate_visualizations = TRUE,
                                                include_comprehensive_metrics = TRUE) {
  
  # Parameter validation and initialization
  validation_result <- validate_differential_analysis_parameters(
    seurat_object = seurat_object,
    group_1_identifiers = group_1_identifiers,
    group_2_identifiers = group_2_identifiers,
    statistical_methodology = statistical_methodology
  )
  
  if (!validation_result$is_valid) {
    stop("Parameter validation failed: ", validation_result$message)
  }
  
  cat("Initiating advanced differential module analysis pipeline...\n")
  cat("Group 1 cells:", length(group_1_identifiers), "\n")
  cat("Group 2 cells:", length(group_2_identifiers), "\n")
  cat("Statistical test:", statistical_methodology, "\n")
  
  # Step 1: Identify differential module eigengenes
  differential_results <- identify_differential_module_eigengenes(
    seurat_object = seurat_object,
    group_1_cells = group_1_identifiers,
    group_2_cells = group_2_identifiers,
    analysis_name = analysis_label,
    statistical_test = statistical_methodology,
    use_harmonized = use_harmonized_eigengenes,
    logfc_threshold = fold_change_threshold,
    min_percentage = minimum_percentage
  )
  
  # Step 2: Generate comprehensive differential analysis report
  if (include_comprehensive_metrics) {
    enriched_results <- enhance_differential_analysis_with_metrics(
      differential_results = differential_results,
      seurat_object = seurat_object,
      group_1_cells = group_1_identifiers,
      group_2_cells = group_2_identifiers,
      analysis_name = analysis_label
    )
    differential_results <- enriched_results$differential_modules
  }
  
  # Step 3: Create visualization
  if (generate_visualizations) {
    volcano_visualization <- create_differential_module_volcano_plot(
      seurat_object = seurat_object,
      differential_modules = differential_results,
      analysis_label = analysis_label,
      display_statistical_cutoffs = TRUE,
      enable_intelligent_labeling = TRUE
    )
    
    print(volcano_visualization)
  }
  
  # Step 4: Generate analysis summary
  analysis_summary <- generate_differential_analysis_summary(
    differential_results = differential_results,
    group_1_size = length(group_1_identifiers),
    group_2_size = length(group_2_identifiers),
    analysis_label = analysis_label
  )
  
  cat("Differential module analysis pipeline completed successfully.\n")
  cat("Significant modules identified:", sum(differential_results$p_val_adj < 0.05, na.rm = TRUE), "\n")
  
  return(list(
    differential_modules = differential_results,
    visualization = if (exists("volcano_visualization")) volcano_visualization else NULL,
    analysis_summary = analysis_summary
  ))
}

#' @title Differential Analysis Parameter Validation Engine
#' @description Comprehensive validation of differential analysis parameters

validate_differential_analysis_parameters <- function(seurat_object, group_1_identifiers,
                                                     group_2_identifiers, statistical_methodology) {
  
  validation_checks <- list()
  
  # Validate Seurat object structure
  if (!inherits(seurat_object, "Seurat")) {
    return(list(is_valid = FALSE, message = "Input must be a valid Seurat object"))
  }
  
  # Validate group 1 barcodes
  missing_group1 <- setdiff(group_1_identifiers, colnames(seurat_object))
  if (length(missing_group1) > 0) {
    validation_checks$group1 <- list(
      is_valid = FALSE,
      message = paste(length(missing_group1), "barcodes in group 1 not found in Seurat object")
    )
  } else {
    validation_checks$group1 <- list(is_valid = TRUE, message = "Group 1 barcodes validated")
  }
  
  # Validate group 2 barcodes
  missing_group2 <- setdiff(group_2_identifiers, colnames(seurat_object))
  if (length(missing_group2) > 0) {
    validation_checks$group2 <- list(
      is_valid = FALSE,
      message = paste(length(missing_group2), "barcodes in group 2 not found in Seurat object")
    )
  } else {
    validation_checks$group2 <- list(is_valid = TRUE, message = "Group 2 barcodes validated")
  }
  
  # Check for barcode overlap
  overlapping_barcodes <- intersect(group_1_identifiers, group_2_identifiers)
  if (length(overlapping_barcodes) > 0) {
    validation_checks$overlap <- list(
      is_valid = FALSE,
      message = paste("Overlapping barcodes detected:", length(overlapping_barcodes))
    )
  } else {
    validation_checks$overlap <- list(is_valid = TRUE, message = "No barcode overlap detected")
  }
  
  # Validate group sizes
  if (length(group_1_identifiers) < 3 || length(group_2_identifiers) < 3) {
    validation_checks$group_sizes <- list(
      is_valid = FALSE,
      message = "Each group must contain at least 3 cells for statistical analysis"
    )
  } else {
    validation_checks$group_sizes <- list(is_valid = TRUE, message = "Group sizes validated")
  }
  
  # Validate statistical methodology
  valid_statistical_tests <- c("wilcox", "t", "bimod", "roc", "negbinom", "poisson", "LR", "MAST")
  if (!statistical_methodology %in% valid_statistical_tests) {
    validation_checks$statistical_test <- list(
      is_valid = FALSE,
      message = paste("Statistical test must be one of:", paste(valid_statistical_tests, collapse = ", "))
    )
  } else {
    validation_checks$statistical_test <- list(is_valid = TRUE, message = "Statistical test validated")
  }
  
  # Validate module eigengenes availability
  if (is.null(seurat_object@misc$harmonized_eigengenes) && is.null(seurat_object@misc$original_eigengenes)) {
    validation_checks$eigengenes <- list(
      is_valid = FALSE,
      message = "Module eigengenes not found. Run module eigengene analysis first."
    )
  } else {
    validation_checks$eigengenes <- list(is_valid = TRUE, message = "Module eigengenes available")
  }
  
  # Compile overall results
  all_valid <- all(sapply(validation_checks, function(x) x$is_valid))
  error_messages <- sapply(validation_checks[!sapply(validation_checks, function(x) x$is_valid)], 
                          function(x) x$message)
  
  return(list(
    is_valid = all_valid,
    message = if (all_valid) "All parameters validated successfully" else paste(error_messages, collapse = "; ")
  ))
}

#' @title Differential Module Eigengene Identifier
#' @description Identifies differentially active modules between groups

identify_differential_module_eigengenes <- function(seurat_object, group_1_cells, group_2_cells,
                                                   analysis_name, statistical_test, use_harmonized,
                                                   logfc_threshold, min_percentage) {
  
  cat("Identifying differential module eigengenes...\n")
  
  # Extract module eigengene data
  eigengene_data <- extract_module_eigengene_matrix(
    seurat_object = seurat_object,
    use_harmonized = use_harmonized
  )
  
  # Prepare data for differential analysis
  prepared_data <- prepare_differential_analysis_data(
    eigengene_matrix = eigengene_data,
    group_1_cells = group_1_cells,
    group_2_cells = group_2_cells
  )
  
  # Perform differential expression analysis
  differential_results <- perform_module_differential_analysis(
    assay_data = prepared_data$assay_object,
    group_1_cells = prepared_data$group_1_cells,
    group_2_cells = prepared_data$group_2_cells,
    statistical_test = statistical_test,
    logfc_threshold = logfc_threshold,
    min_percentage = min_percentage
  )
  
  # Enhance results with additional metrics
  enhanced_results <- enhance_differential_results(
    differential_results = differential_results,
    eigengene_matrix = eigengene_data,
    group_1_cells = group_1_cells,
    group_2_cells = group_2_cells
  )
  
  # Store analysis parameters
  enhanced_results$analysis_parameters <- list(
    statistical_test = statistical_test,
    use_harmonized_eigengenes = use_harmonized,
    group_1_size = length(group_1_cells),
    group_2_size = length(group_2_cells),
    analysis_timestamp = Sys.time()
  )
  
  cat("Differential analysis completed. Modules analyzed:", nrow(enhanced_results), "\n")
  return(enhanced_results)
}

#' @title Module Eigengene Matrix Extractor
#' @description Extracts module eigengene matrix with proper filtering

extract_module_eigengene_matrix <- function(seurat_object, use_harmonized) {
  if (use_harmonized) {
    if (is.null(seurat_object@misc$harmonized_eigengenes)) {
      warning("Harmonized eigengenes not found, using original eigengenes")
      eigengene_matrix <- seurat_object@misc$original_eigengenes
    } else {
      eigengene_matrix <- seurat_object@misc$harmonized_eigengenes
    }
  } else {
    eigengene_matrix <- seurat_object@misc$original_eigengenes
  }
  
  # Remove grey module and ensure non-negative values
  valid_modules <- colnames(eigengene_matrix)[colnames(eigengene_matrix) != "grey"]
  eigengene_matrix <- eigengene_matrix[, valid_modules, drop = FALSE]
  
  # Apply non-negative constraint
  eigengene_matrix[eigengene_matrix < 0] <- 0
  
  return(eigengene_matrix)
}

#' @title Differential Analysis Data Preparer
#' @description Prepares data for differential analysis

prepare_differential_analysis_data <- function(eigengene_matrix, group_1_cells, group_2_cells) {
  # Ensure all cells are present in the eigengene matrix
  available_group1 <- intersect(group_1_cells, rownames(eigengene_matrix))
  available_group2 <- intersect(group_2_cells, rownames(eigengene_matrix))
  
  if (length(available_group1) == 0 || length(available_group2) == 0) {
    stop("No overlapping cells found between groups and eigengene matrix")
  }
  
  # Transpose matrix for Seurat compatibility (genes as rows, cells as columns)
  transposed_matrix <- t(eigengene_matrix)
  
  # Create Seurat assay object
  assay_object <- Seurat::CreateAssayObject(counts = transposed_matrix)
  
  return(list(
    assay_object = assay_object,
    group_1_cells = available_group1,
    group_2_cells = available_group2
  ))
}

#' @title Module Differential Analysis Performer
#' @description Performs statistical testing for module differences

perform_module_differential_analysis <- function(assay_data, group_1_cells, group_2_cells,
                                                statistical_test, logfc_threshold, min_percentage) {
  
  # Use Seurat's FindMarkers function for differential analysis
  differential_results <- Seurat::FindMarkers(
    object = assay_data,
    cells.1 = group_1_cells,
    cells.2 = group_2_cells,
    slot = "counts",
    test.use = statistical_test,
    only.pos = FALSE,
    logfc.threshold = logfc_threshold,
    min.pct = min_percentage,
    verbose = FALSE
  )
  
  # Add module names as a column
  differential_results$module <- rownames(differential_results)
  rownames(differential_results) <- NULL
  
  return(differential_results)
}

#' @title Differential Results Enhancer
#' @description Enhances differential results with additional metrics

enhance_differential_results <- function(differential_results, eigengene_matrix,
                                        group_1_cells, group_2_cells) {
  
  # Calculate additional metrics for each module
  enhanced_results <- differential_results %>%
    rowwise() %>%
    mutate(
      # Calculate mean expression in each group
      mean_expression_group1 = mean(eigengene_matrix[group_1_cells, module], na.rm = TRUE),
      mean_expression_group2 = mean(eigengene_matrix[group_2_cells, module], na.rm = TRUE),
      
      # Calculate fold change manually for verification
      manual_fold_change = mean_expression_group1 - mean_expression_group2,
      
      # Calculate percentage of cells expressing the module
      pct_expressed_group1 = sum(eigengene_matrix[group_1_cells, module] > 0) / length(group_1_cells) * 100,
      pct_expressed_group2 = sum(eigengene_matrix[group_2_cells, module] > 0) / length(group_2_cells) * 100,
      
      # Calculate effect size (Cohen's d)
      sd_group1 = sd(eigengene_matrix[group_1_cells, module], na.rm = TRUE),
      sd_group2 = sd(eigengene_matrix[group_2_cells, module], na.rm = TRUE),
      pooled_sd = sqrt((sd_group1^2 + sd_group2^2) / 2),
      cohens_d = ifelse(pooled_sd > 0, manual_fold_change / pooled_sd, NA),
      
      # Significance categories
      significance_level = case_when(
        p_val_adj < 0.001 ~ "***",
        p_val_adj < 0.01 ~ "**",
        p_val_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    ungroup()
  
  return(enhanced_results)
}

#' @title Differential Analysis Metrics Enhancer
#' @description Adds comprehensive metrics to differential analysis results

enhance_differential_analysis_with_metrics <- function(differential_results, seurat_object,
                                                      group_1_cells, group_2_cells, analysis_name) {
  
  cat("Enhancing differential analysis with comprehensive metrics...\n")
  
  # Retrieve module assignment information
  module_assignments <- seurat_object@misc$module_assignments
  
  # Add module size information
  module_sizes <- module_assignments %>%
    filter(module != "grey") %>%
    group_by(module) %>%
    summarize(module_size = n(), .groups = 'drop')
  
  enhanced_results <- differential_results %>%
    left_join(module_sizes, by = "module") %>%
    mutate(
      # Calculate module conservation score
      conservation_score = -log10(p_val_adj) * abs(avg_log2FC),
      
      # Categorize direction of change
      direction = ifelse(avg_log2FC > 0, "upregulated", "downregulated"),
      
      # Calculate statistical power estimate
      statistical_power = pmin(1, abs(avg_log2FC) * sqrt(pmin(length(group_1_cells), length(group_2_cells)))),
      
      # Quality score combining multiple metrics
      quality_score = (-log10(p_val_adj) + abs(avg_log2FC) + cohens_d) / 3
    )
  
  # Store enhanced results
  seurat_object@misc[[analysis_name]]$enhanced_differential_results <- enhanced_results
  
  return(list(
    differential_modules = enhanced_results,
    seurat_object = seurat_object
  ))
}

#' @title Differential Module Volcano Plot Creator
#' @description Creates volcano plot visualization of differential modules

create_differential_module_volcano_plot <- function(seurat_object, differential_modules,
                                                   analysis_label, display_statistical_cutoffs,
                                                   enable_intelligent_labeling) {
  
  cat("Creating differential module volcano plot...\n")
  
  # Prepare data for plotting
  plot_data <- prepare_volcano_plot_data(
    differential_modules = differential_modules,
    seurat_object = seurat_object
  )
  
  # Create base volcano plot
  volcano_plot <- create_base_volcano_plot(
    plot_data = plot_data,
    display_cutoffs = display_statistical_cutoffs
  )
  
  # Enhance plot with intelligent labeling
  if (enable_intelligent_labeling) {
    volcano_plot <- add_intelligent_plot_labels(
      base_plot = volcano_plot,
      plot_data = plot_data
    )
  }
  
  # Apply professional styling
  finalized_plot <- apply_professional_plot_styling(volcano_plot)
  
  return(finalized_plot)
}

#' @title Volcano Plot Data Preparer
#' @description Prepares data for volcano plot visualization

prepare_volcano_plot_data <- function(differential_modules, seurat_object) {
  # Remove NA values and handle infinite fold changes
  plot_data <- differential_modules %>%
    filter(!is.na(p_val_adj), !is.na(avg_log2FC)) %>%
    mutate(
      # Handle zero p-values by setting to minimum non-zero value
      adjusted_p_value = ifelse(p_val_adj == 0, 
                               min(p_val_adj[p_val_adj > 0], na.rm = TRUE), 
                               p_val_adj),
      
      # Handle infinite fold changes
      finite_log2FC = case_when(
        is.infinite(avg_log2FC) & avg_log2FC > 0 ~ max(avg_log2FC[is.finite(avg_log2FC)], na.rm = TRUE),
        is.infinite(avg_log2FC) & avg_log2FC < 0 ~ min(avg_log2FC[is.finite(avg_log2FC)], na.rm = TRUE),
        TRUE ~ avg_log2FC
      ),
      
      # Calculate transformed p-value
      transformed_pvalue = -log10(adjusted_p_value)
    )
  
  # Add module color information
  module_colors <- seurat_object@misc$module_assignments %>%
    filter(module != "grey") %>%
    select(module, color) %>%
    distinct()
  
  plot_data <- plot_data %>%
    left_join(module_colors, by = "module") %>%
    mutate(
      # Create intelligent labels for significant modules
      display_label = ifelse(p_val_adj < 0.05, module, "")
    )
  
  return(plot_data)
}

#' @title Base Volcano Plot Creator
#' @description Creates the base volcano plot structure

create_base_volcano_plot <- function(plot_data, display_cutoffs) {
  base_plot <- ggplot(plot_data, aes(x = finite_log2FC, y = transformed_pvalue, 
                                    fill = module, color = module)) +
    geom_point(shape = 21, size = 4, color = "black", alpha = 0.8)
  
  # Add statistical cutoffs if requested
  if (display_cutoffs) {
    base_plot <- base_plot +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", alpha = 0.7) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -log10(0.05), ymax = Inf,
               fill = "grey80", alpha = 0.1)
  }
  
  return(base_plot)
}

#' @title Intelligent Plot Label Adder
#' @description Adds intelligent labels to the volcano plot

add_intelligent_plot_labels <- function(base_plot, plot_data) {
  # Identify modules for labeling (most significant or largest fold change)
  labeling_candidates <- plot_data %>%
    filter(p_val_adj < 0.05) %>%
    arrange(p_val_adj, desc(abs(finite_log2FC))) %>%
    head(20)  # Limit number of labels
  
  if (nrow(labeling_candidates) > 0) {
    base_plot <- base_plot +
      ggrepel::geom_text_repel(
        data = labeling_candidates,
        aes(label = module),
        color = "black",
        size = 3,
        max.overlaps = 20,
        min.segment.length = 0.1,
        box.padding = 0.5,
        point.padding = 0.1
      )
  }
  
  return(base_plot)
}

#' @title Professional Plot Styling Applicator
#' @description Applies professional styling to plot

apply_professional_plot_styling <- function(plot_object) {
  styled_plot <- plot_object +
    scale_fill_manual(values = setNames(unique(plot_object$data$color), 
                                       unique(plot_object$data$module))) +
    scale_color_manual(values = setNames(unique(plot_object$data$color), 
                                        unique(plot_object$data$module))) +
    labs(
      x = expression("Log"[2] ~ "Fold Change"),
      y = expression("-Log"[10] ~ "Adjusted P-value"),
      title = "Differential Module Eigengene Analysis",
      subtitle = "Volcano Plot of Module Activity Differences"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(styled_plot)
}

#' @title Differential Analysis Summary Generator
#' @description Generates comprehensive summary of differential analysis

generate_differential_analysis_summary <- function(differential_results, group_1_size,
                                                  group_2_size, analysis_label) {
  
  significant_modules <- differential_results %>% filter(p_val_adj < 0.05)
  
  summary <- list(
    analysis_parameters = differential_results$analysis_parameters,
    summary_statistics = list(
      total_modules_analyzed = nrow(differential_results),
      significant_modules = nrow(significant_modules),
      upregulated_modules = sum(significant_modules$avg_log2FC > 0, na.rm = TRUE),
      downregulated_modules = sum(significant_modules$avg_log2FC < 0, na.rm = TRUE),
      strongest_upregulation = if (nrow(significant_modules) > 0) {
        max(significant_modules$avg_log2FC, na.rm = TRUE)
      } else { NA },
      strongest_downregulation = if (nrow(significant_modules) > 0) {
        min(significant_modules$avg_log2FC, na.rm = TRUE)
      } else { NA }
    ),
    group_information = list(
      group_1_size = group_1_size,
      group_2_size = group_2_size
    )
  )
  
  cat("Differential analysis summary:\n")
  cat("  - Significant modules:", summary$summary_statistics$significant_modules, "\n")
  cat("  - Upregulated:", summary$summary_statistics$upregulated_modules, "\n")
  cat("  - Downregulated:", summary$summary_statistics$downregulated_modules, "\n")
  
  return(summary)
}

# Example usage function
demonstrate_differential_analysis <- function() {
  cat("Differential module analysis pipeline demonstration\n")
  cat("Please use with actual Seurat object containing module eigengenes\n")
}

# Uncomment to test (with actual Seurat object)
# results <- execute_differential_module_pipeline(seurat_object, group1_cells, group2_cells)
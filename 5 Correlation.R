

library(Seurat)
library(Hmisc)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(corrplot)

#' @title Main Module-Trait Correlation Pipeline
#' @description Orchestrates the complete module-trait correlation analysis workflow

execute_module_trait_analysis_pipeline <- function(seurat_object,
                                                  phenotypic_traits = c('CRC_progression'),
                                                  grouping_variable = 'cell_type',
                                                  module_representation = "harmonized_eigengenes",
                                                  correlation_methodology = "pearson",
                                                  subset_variable = NULL,
                                                  subset_values = NULL,
                                                  analysis_identifier = NULL,
                                                  generate_visualizations = TRUE) {
  
  # Parameter validation and initialization
  validation_result <- validate_module_trait_parameters(
    seurat_object = seurat_object,
    phenotypic_traits = phenotypic_traits,
    grouping_variable = grouping_variable,
    module_representation = module_representation,
    correlation_methodology = correlation_methodology
  )
  
  if (!validation_result$is_valid) {
    stop("Parameter validation failed: ", validation_result$message)
  }
  
  analysis_identifier <- analysis_identifier %||% seurat_object@misc$active_analysis
  
  cat("Initiating advanced module-trait correlation analysis pipeline...\n")
  cat("Phenotypic traits:", paste(phenotypic_traits, collapse = ", "), "\n")
  cat("Module representation:", module_representation, "\n")
  
  # Step 1: Prepare data for correlation analysis
  prepared_object <- prepare_module_trait_data(
    seurat_object = seurat_object,
    phenotypic_traits = phenotypic_traits,
    module_representation = module_representation,
    subset_variable = subset_variable,
    subset_values = subset_values,
    analysis_name = analysis_identifier
  )
  
  # Step 2: Perform comprehensive correlation analysis
  correlation_analyzed_object <- perform_comprehensive_correlation_analysis(
    seurat_object = prepared_object,
    phenotypic_traits = phenotypic_traits,
    grouping_variable = grouping_variable,
    correlation_methodology = correlation_methodology,
    analysis_name = analysis_identifier
  )
  
  # Step 3: Generate statistical summaries
  summary_object <- generate_correlation_statistical_summaries(
    seurat_object = correlation_analyzed_object,
    analysis_name = analysis_identifier
  )
  
  # Step 4: Create visualizations
  if (generate_visualizations) {
    visualization_results <- create_module_trait_visualizations(
      seurat_object = summary_object,
      analysis_name = analysis_identifier
    )
    
    # Display correlation heatmap
    if (!is.null(visualization_results$correlation_heatmap)) {
      print(visualization_results$correlation_heatmap)
    }
  }
  
  # Step 5: Generate analysis report
  generate_module_trait_analysis_report(summary_object, analysis_identifier)
  
  cat("Module-trait correlation analysis pipeline completed successfully.\n")
  
  return(list(
    seurat_object = summary_object,
    visualizations = if (exists("visualization_results")) visualization_results else NULL
  ))
}

#' @title Module-Trait Parameter Validation Engine
#' @description Comprehensive validation of module-trait correlation parameters

validate_module_trait_parameters <- function(seurat_object, phenotypic_traits,
                                            grouping_variable, module_representation,
                                            correlation_methodology) {
  
  validation_checks <- list()
  
  # Validate Seurat object structure
  if (!inherits(seurat_object, "Seurat")) {
    return(list(is_valid = FALSE, message = "Input must be a valid Seurat object"))
  }
  
  # Validate phenotypic traits exist in metadata
  missing_traits <- setdiff(phenotypic_traits, colnames(seurat_object@meta.data))
  if (length(missing_traits) > 0) {
    validation_checks$traits <- list(
      is_valid = FALSE,
      message = paste("Traits not found in metadata:", paste(missing_traits, collapse = ", "))
    )
  } else {
    validation_checks$traits <- list(is_valid = TRUE, message = "Phenotypic traits validated")
  }
  
  # Validate trait data types
  valid_trait_types <- c("numeric", "integer", "factor")
  trait_data_types <- sapply(phenotypic_traits, function(trait) {
    class(seurat_object@meta.data[[trait]])
  })
  
  invalid_trait_types <- phenotypic_traits[!trait_data_types %in% valid_trait_types]
  if (length(invalid_trait_types) > 0) {
    validation_checks$trait_types <- list(
      is_valid = FALSE,
      message = paste("Invalid trait types:", paste(invalid_trait_types, collapse = ", "),
                     ". Must be numeric, integer, or factor.")
    )
  } else {
    validation_checks$trait_types <- list(is_valid = TRUE, message = "Trait data types validated")
  }
  
  # Validate module representation option
  valid_representations <- c("harmonized_eigengenes", "original_eigengenes", "module_scores")
  if (!module_representation %in% valid_representations) {
    validation_checks$representation <- list(
      is_valid = FALSE,
      message = paste("Module representation must be one of:", paste(valid_representations, collapse = ", "))
    )
  } else {
    validation_checks$representation <- list(is_valid = TRUE, message = "Module representation validated")
  }
  
  # Validate correlation methodology
  valid_correlation_methods <- c("pearson", "spearman")
  if (!correlation_methodology %in% valid_correlation_methods) {
    validation_checks$correlation_method <- list(
      is_valid = FALSE,
      message = paste("Correlation method must be one of:", paste(valid_correlation_methods, collapse = ", "))
    )
  } else {
    validation_checks$correlation_method <- list(is_valid = TRUE, message = "Correlation method validated")
  }
  
  # Validate grouping variable if provided
  if (!is.null(grouping_variable)) {
    if (!grouping_variable %in% colnames(seurat_object@meta.data)) {
      validation_checks$grouping_var <- list(
        is_valid = FALSE,
        message = paste("Grouping variable", grouping_variable, "not found in metadata")
      )
    } else {
      validation_checks$grouping_var <- list(is_valid = TRUE, message = "Grouping variable validated")
    }
  }
  
  # Validate module assignments exist
  if (is.null(seurat_object@misc$module_assignments)) {
    validation_checks$module_assignments <- list(
      is_valid = FALSE,
      message = "Module assignments not found. Run co-expression network analysis first."
    )
  } else {
    validation_checks$module_assignments <- list(is_valid = TRUE, message = "Module assignments validated")
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

#' @title Module-Trait Data Preparation Engine
#' @description Prepares data for correlation analysis

prepare_module_trait_data <- function(seurat_object, phenotypic_traits, module_representation,
                                     subset_variable, subset_values, analysis_name) {
  
  cat("Preparing data for module-trait correlation analysis...\n")
  
  # Extract module representation data
  module_data <- extract_module_representation_data(
    seurat_object = seurat_object,
    representation_type = module_representation,
    analysis_name = analysis_name
  )
  
  # Process phenotypic traits
  processed_traits <- process_phenotypic_traits(
    seurat_object = seurat_object,
    phenotypic_traits = phenotypic_traits
  )
  
  # Apply subsetting if requested
  if (!is.null(subset_variable)) {
    subset_object <- apply_data_subsetting(
      seurat_object = seurat_object,
      module_data = module_data,
      subset_variable = subset_variable,
      subset_values = subset_values
    )
    seurat_object <- subset_object$seurat_object
    module_data <- subset_object$module_data
  }
  
  # Validate data compatibility
  validate_data_compatibility(module_data, processed_traits$trait_data)
  
  # Store prepared data
  seurat_object@misc[[analysis_name]]$module_trait_data <- list(
    module_representation = module_data,
    phenotypic_traits = processed_traits$trait_data,
    trait_processing_log = processed_traits$processing_log,
    preparation_timestamp = Sys.time()
  )
  
  cat("Data preparation completed. Samples:", nrow(module_data), "Modules:", ncol(module_data), "\n")
  return(seurat_object)
}

#' @title Module Representation Data Extractor
#' @description Extracts module representation data based on specified type

extract_module_representation_data <- function(seurat_object, representation_type, analysis_name) {
  module_data <- switch(representation_type,
    "harmonized_eigengenes" = {
      if (is.null(seurat_object@misc[[analysis_name]]$harmonized_eigengenes)) {
        stop("Harmonized eigengenes not found. Run module eigengene analysis first.")
      }
      seurat_object@misc[[analysis_name]]$harmonized_eigengenes
    },
    "original_eigengenes" = {
      if (is.null(seurat_object@misc[[analysis_name]]$original_eigengenes)) {
        stop("Original eigengenes not found. Run module eigengene analysis first.")
      }
      seurat_object@misc[[analysis_name]]$original_eigengenes
    },
    "module_scores" = {
      if (is.null(seurat_object@misc$module_gene_scores)) {
        stop("Module scores not found. Run module scoring analysis first.")
      }
      # Extract module score columns
      score_columns <- grep("ModuleScore", colnames(seurat_object@meta.data), value = TRUE)
      as.matrix(seurat_object@meta.data[, score_columns, drop = FALSE])
    }
  )
  
  # Filter out grey module if present
  valid_modules <- extract_valid_module_names(seurat_object)
  module_data <- module_data[, colnames(module_data) %in% valid_modules, drop = FALSE]
  
  return(as.data.frame(module_data))
}

#' @title Valid Module Name Extractor
#' @description Extracts valid module names excluding grey module

extract_valid_module_names <- function(seurat_object) {
  module_assignments <- seurat_object@misc$module_assignments
  valid_modules <- unique(module_assignments$module)
  valid_modules <- valid_modules[valid_modules != "grey" & !is.na(valid_modules)]
  return(valid_modules)
}

#' @title Phenotypic Trait Processor
#' @description Processes and validates phenotypic trait data

process_phenotypic_traits <- function(seurat_object, phenotypic_traits) {
  trait_data <- seurat_object@meta.data[, phenotypic_traits, drop = FALSE]
  processing_log <- list()
  
  for (trait_name in phenotypic_traits) {
    trait_vector <- trait_data[[trait_name]]
    
    # Handle factor traits
    if (is.factor(trait_vector)) {
      original_levels <- levels(trait_vector)
      trait_data[[trait_name]] <- as.numeric(trait_vector)
      processing_log[[trait_name]] <- list(
        original_type = "factor",
        converted_to = "numeric",
        levels = original_levels
      )
      warning(sprintf("Converted factor trait '%s' to numeric. Original levels: %s", 
                     trait_name, paste(original_levels, collapse = ", ")))
    }
    
    # Check for missing values
    missing_count <- sum(is.na(trait_data[[trait_name]]))
    if (missing_count > 0) {
      warning(sprintf("Trait '%s' contains %d missing values", trait_name, missing_count))
    }
  }
  
  return(list(trait_data = trait_data, processing_log = processing_log))
}

#' @title Data Subsetting Applicator
#' @description Applies subsetting to data based on specified criteria

apply_data_subsetting <- function(seurat_object, module_data, subset_variable, subset_values) {
  if (!subset_variable %in% colnames(seurat_object@meta.data)) {
    stop("Subset variable '", subset_variable, "' not found in metadata")
  }
  
  subset_indices <- seurat_object@meta.data[[subset_variable]] %in% subset_values
  if (sum(subset_indices) == 0) {
    stop("No cells match the specified subset values")
  }
  
  seurat_object <- seurat_object[, subset_indices]
  module_data <- module_data[subset_indices, , drop = FALSE]
  
  cat("Data subsetting applied. Retained", sum(subset_indices), "cells\n")
  return(list(seurat_object = seurat_object, module_data = module_data))
}

#' @title Data Compatibility Validator
#' @description Validates compatibility between module and trait data

validate_data_compatibility <- function(module_data, trait_data) {
  if (nrow(module_data) != nrow(trait_data)) {
    stop("Module data and trait data have different numbers of samples")
  }
  
  if (nrow(module_data) < 10) {
    warning("Low sample size (n < 10) may affect correlation reliability")
  }
}

#' @title Comprehensive Correlation Analysis Engine
#' @description Performs module-trait correlation analysis

perform_comprehensive_correlation_analysis <- function(seurat_object, phenotypic_traits,
                                                      grouping_variable, correlation_methodology,
                                                      analysis_name) {
  
  cat("Performing comprehensive module-trait correlation analysis...\n")
  
  module_trait_data <- seurat_object@misc[[analysis_name]]$module_trait_data
  module_data <- module_trait_data$module_representation
  trait_data <- module_trait_data$phenotypic_traits
  
  # Perform overall correlation analysis
  overall_results <- compute_correlation_analysis(
    trait_matrix = as.matrix(trait_data),
    module_matrix = as.matrix(module_data),
    correlation_method = correlation_methodology,
    analysis_label = "all_cells"
  )
  
  correlation_results <- list(overall = overall_results)
  
  # Perform group-wise correlations if grouping variable provided
  if (!is.null(grouping_variable)) {
    group_results <- perform_group_wise_correlation_analysis(
      seurat_object = seurat_object,
      trait_data = trait_data,
      module_data = module_data,
      grouping_variable = grouping_variable,
      correlation_methodology = correlation_methodology
    )
    correlation_results <- c(correlation_results, group_results)
  }
  
  # Store correlation results
  seurat_object@misc[[analysis_name]]$correlation_results <- list(
    correlation_matrices = lapply(correlation_results, function(x) x$correlation),
    p_value_matrices = lapply(correlation_results, function(x) x$p_values),
    adjusted_p_value_matrices = lapply(correlation_results, function(x) x$adjusted_p_values),
    analysis_parameters = list(
      correlation_method = correlation_methodology,
      grouping_variable = grouping_variable,
      analysis_timestamp = Sys.time()
    )
  )
  
  cat("Correlation analysis completed. Groups analyzed:", length(correlation_results), "\n")
  return(seurat_object)
}

#' @title Correlation Analysis Computer
#' @description Computes correlations between traits and modules

compute_correlation_analysis <- function(trait_matrix, module_matrix, correlation_method, analysis_label) {
  # Compute correlation matrix using Hmisc
  correlation_result <- Hmisc::rcorr(trait_matrix, module_matrix, type = correlation_method)
  
  # Extract correlation coefficients and p-values
  correlation_matrix <- correlation_result$r[colnames(trait_matrix), colnames(module_matrix), drop = FALSE]
  p_value_matrix <- correlation_result$P[colnames(trait_matrix), colnames(module_matrix), drop = FALSE]
  
  # Adjust p-values for multiple testing
  adjusted_p_values <- adjust_p_values(p_value_matrix)
  
  return(list(
    correlation = correlation_matrix,
    p_values = p_value_matrix,
    adjusted_p_values = adjusted_p_values,
    analysis_label = analysis_label
  ))
}

#' @title P-value Adjustment Applicator
#' @description Applies multiple testing correction to p-values

adjust_p_values <- function(p_value_matrix) {
  # Flatten p-value matrix for adjustment
  p_value_vector <- as.vector(p_value_matrix)
  adjusted_p_vector <- p.adjust(p_value_vector, method = "fdr")
  
  # Reshape back to matrix format
  adjusted_p_matrix <- matrix(adjusted_p_vector, 
                             nrow = nrow(p_value_matrix), 
                             ncol = ncol(p_value_matrix))
  rownames(adjusted_p_matrix) <- rownames(p_value_matrix)
  colnames(adjusted_p_matrix) <- colnames(p_value_matrix)
  
  return(adjusted_p_matrix)
}

#' @title Group-wise Correlation Analysis Performer
#' @description Performs correlation analysis within each group

perform_group_wise_correlation_analysis <- function(seurat_object, trait_data, module_data,
                                                   grouping_variable, correlation_methodology) {
  
  group_levels <- unique(seurat_object@meta.data[[grouping_variable]])
  group_results <- list()
  
  for (group_name in group_levels) {
    group_indices <- seurat_object@meta.data[[grouping_variable]] == group_name
    
    if (sum(group_indices) >= 5) {  # Require minimum group size
      group_trait_data <- trait_data[group_indices, , drop = FALSE]
      group_module_data <- module_data[group_indices, , drop = FALSE]
      
      group_correlation <- compute_correlation_analysis(
        trait_matrix = as.matrix(group_trait_data),
        module_matrix = as.matrix(group_module_data),
        correlation_method = correlation_methodology,
        analysis_label = group_name
      )
      
      group_results[[group_name]] <- group_correlation
    } else {
      warning(sprintf("Group '%s' has insufficient samples (%d), skipping analysis", 
                     group_name, sum(group_indices)))
    }
  }
  
  return(group_results)
}

#' @title Correlation Statistical Summary Generator
#' @description Generates statistical summaries of correlation results

generate_correlation_statistical_summaries <- function(seurat_object, analysis_name) {
  correlation_results <- seurat_object@misc[[analysis_name]]$correlation_results
  
  # Calculate summary statistics
  summary_stats <- calculate_correlation_summary_statistics(correlation_results)
  
  # Identify significant correlations
  significant_findings <- identify_significant_correlations(correlation_results)
  
  # Store summaries
  seurat_object@misc[[analysis_name]]$correlation_summaries <- list(
    summary_statistics = summary_stats,
    significant_correlations = significant_findings,
    summary_timestamp = Sys.time()
  )
  
  cat("Statistical summaries generated. Significant correlations:", 
      nrow(significant_findings$overall), "\n")
  return(seurat_object)
}

#' @title Correlation Summary Statistics Calculator
#' @description Calculates descriptive statistics for correlation results

calculate_correlation_summary_statistics <- function(correlation_results) {
  summary_stats <- list()
  
  for (analysis_group in names(correlation_results$correlation_matrices)) {
    cor_matrix <- correlation_results$correlation_matrices[[analysis_group]]
    p_matrix <- correlation_results$adjusted_p_value_matrices[[analysis_group]]
    
    group_stats <- list(
      mean_correlation = mean(cor_matrix, na.rm = TRUE),
      median_correlation = median(cor_matrix, na.rm = TRUE),
      max_correlation = max(cor_matrix, na.rm = TRUE),
      min_correlation = min(cor_matrix, na.rm = TRUE),
      significant_correlations = sum(p_matrix < 0.05, na.rm = TRUE),
      total_correlations = length(cor_matrix)
    )
    
    summary_stats[[analysis_group]] <- group_stats
  }
  
  return(summary_stats)
}

#' @title Significant Correlation Identifier
#' @description Identifies statistically significant correlations

identify_significant_correlations <- function(correlation_results) {
  significant_findings <- list()
  
  for (analysis_group in names(correlation_results$correlation_matrices)) {
    cor_matrix <- correlation_results$correlation_matrices[[analysis_group]]
    p_matrix <- correlation_results$adjusted_p_value_matrices[[analysis_group]]
    
    # Find significant correlations (FDR < 0.05)
    significant_indices <- which(p_matrix < 0.05, arr.ind = TRUE)
    
    if (nrow(significant_indices) > 0) {
      significant_correlations <- data.frame(
        trait = rownames(cor_matrix)[significant_indices[, 1]],
        module = colnames(cor_matrix)[significant_indices[, 2]],
        correlation = cor_matrix[significant_indices],
        adjusted_p_value = p_matrix[significant_indices],
        stringsAsFactors = FALSE
      )
      
      # Sort by absolute correlation strength
      significant_correlations <- significant_correlations[
        order(-abs(significant_correlations$correlation)),
      ]
    } else {
      significant_correlations <- data.frame(
        trait = character(),
        module = character(),
        correlation = numeric(),
        adjusted_p_value = numeric(),
        stringsAsFactors = FALSE
      )
    }
    
    significant_findings[[analysis_group]] <- significant_correlations
  }
  
  return(significant_findings)
}

#' @title Module-Trait Visualization Generator
#' @description Creates visualizations for module-trait correlations

create_module_trait_visualizations <- function(seurat_object, analysis_name) {
  cat("Generating module-trait correlation visualizations...\n")
  
  visualization_results <- list()
  correlation_results <- seurat_object@misc[[analysis_name]]$correlation_results
  
  # Create correlation heatmap for overall analysis
  if (!is.null(correlation_results$correlation_matrices$overall)) {
    visualization_results$correlation_heatmap <- create_correlation_heatmap(
      correlation_matrix = correlation_results$correlation_matrices$overall,
      p_value_matrix = correlation_results$adjusted_p_value_matrices$overall,
      plot_title = "Module-Trait Correlation Heatmap"
    )
  }
  
  # Create summary bar plot
  visualization_results$summary_plot <- create_correlation_summary_plot(
    correlation_results = correlation_results
  )
  
  return(visualization_results)
}

#' @title Correlation Heatmap Creator
#' @description Creates heatmap visualization of correlation matrix

create_correlation_heatmap <- function(correlation_matrix, p_value_matrix, plot_title) {
  # Prepare annotation for significant correlations
  significance_annotation <- matrix("", nrow = nrow(p_value_matrix), ncol = ncol(p_value_matrix))
  significance_annotation[p_value_matrix < 0.05] <- "*"
  significance_annotation[p_value_matrix < 0.01] <- "**"
  significance_annotation[p_value_matrix < 0.001] <- "***"
  
  # Create heatmap using corrplot
  corrplot::corrplot(
    correlation_matrix,
    method = "color",
    type = "full",
    order = "hclust",
    tl.cex = 0.8,
    tl.col = "black",
    title = plot_title,
    mar = c(0, 0, 2, 0),
    p.mat = p_value_matrix,
    sig.level = 0.05,
    insig = "label_sig",
    pch.cex = 1.2
  )
}

#' @title Correlation Summary Plot Creator
#' @description Creates summary visualization of correlation results

create_correlation_summary_plot <- function(correlation_results) {
  # This function would create additional summary visualizations
  # Implementation depends on specific visualization requirements
  return(NULL)  # Placeholder for actual implementation
}

#' @title Module-Trait Analysis Report Generator
#' @description Generates comprehensive report of analysis results

generate_module_trait_analysis_report <- function(seurat_object, analysis_name) {
  correlation_results <- seurat_object@misc[[analysis_name]]$correlation_results
  correlation_summaries <- seurat_object@misc[[analysis_name]]$correlation_summaries
  
  report <- list(
    analysis_overview = list(
      correlation_method = correlation_results$analysis_parameters$correlation_method,
      groups_analyzed = length(correlation_results$correlation_matrices),
      total_correlations_computed = correlation_summaries$summary_statistics$overall$total_correlations,
      analysis_timestamp = correlation_results$analysis_parameters$analysis_timestamp
    ),
    significant_findings = list(
      overall_significant = nrow(correlation_summaries$significant_correlations$overall),
      strongest_correlation = if (nrow(correlation_summaries$significant_correlations$overall) > 0) {
        correlation_summaries$significant_correlations$overall$correlation[1]
      } else {
        "None"
      }
    )
  )
  
  seurat_object@misc[[analysis_name]]$module_trait_report <- report
  
  cat("Module-trait analysis report:\n")
  cat("  - Correlation method:", report$analysis_overview$correlation_method, "\n")
  cat("  - Significant correlations:", report$significant_findings$overall_significant, "\n")
  cat("  - Strongest correlation:", report$significant_findings$strongest_correlation, "\n")
  
  return(seurat_object)
}

# Example usage function
demonstrate_module_trait_analysis <- function() {
  cat("Module-trait correlation analysis pipeline demonstration\n")
  cat("Please use with actual Seurat object containing module assignments and phenotypic traits\n")
}

# Uncomment to test (with actual Seurat object)
# results <- execute_module_trait_analysis_pipeline(your_seurat_object)
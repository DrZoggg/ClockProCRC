

library(Seurat)
library(harmony)
library(UCell)
library(patchwork)
library(dplyr)
library(purrr)
library(ggplot2)

#' @title Main Module Eigengene Analysis Pipeline
#' @description Orchestrates the complete module eigengene analysis workflow

execute_module_eigengene_pipeline <- function(seurat_object,
                                             harmonization_variable = "Sample",
                                             cell_grouping_variable = "cell_type",
                                             target_cell_group = "SC",
                                             top_gene_count = 25,
                                             scoring_methodology = "UCell",
                                             generate_visualizations = TRUE,
                                             analysis_identifier = NULL) {
  
  # Parameter validation and initialization
  validation_result <- validate_eigengene_parameters(
    seurat_object = seurat_object,
    harmonization_variable = harmonization_variable,
    cell_grouping_variable = cell_grouping_variable,
    target_cell_group = target_cell_group,
    scoring_methodology = scoring_methodology
  )
  
  if (!validation_result$is_valid) {
    stop("Parameter validation failed: ", validation_result$message)
  }
  
  analysis_identifier <- analysis_identifier %||% seurat_object@misc$active_analysis
  
  cat("Initiating advanced module eigengene analysis pipeline...\n")
  cat("Harmonization variable:", harmonization_variable, "\n")
  cat("Target cell group:", target_cell_group, "\n")
  
  # Step 1: Scale data for improved analysis
  scaled_object <- prepare_data_for_eigengene_analysis(seurat_object)
  
  # Step 2: Compute module eigengenes
  eigengene_computed_object <- compute_comprehensive_module_eigengenes(
    seurat_object = scaled_object,
    harmonization_covariate = harmonization_variable,
    analysis_name = analysis_identifier
  )
  
  # Step 3: Extract eigengene matrices
  eigengene_matrices <- extract_eigengene_matrices(eigengene_computed_object, analysis_identifier)
  
  # Step 4: Calculate gene connectivity metrics
  connectivity_analyzed_object <- calculate_gene_module_connectivity(
    seurat_object = eigengene_computed_object,
    grouping_variable = cell_grouping_variable,
    target_group = target_cell_group,
    analysis_name = analysis_identifier
  )
  
  # Step 5: Score cells using top module genes
  scored_object <- compute_cell_module_scores(
    seurat_object = connectivity_analyzed_object,
    top_genes_per_module = top_gene_count,
    scoring_algorithm = scoring_methodology,
    analysis_name = analysis_identifier
  )
  
  # Step 6: Generate visualizations
  if (generate_visualizations) {
    visualization_results <- create_eigengene_visualizations(
      seurat_object = scored_object,
      analysis_name = analysis_identifier
    )
    
    # Display combined feature plots
    if (!is.null(visualization_results$feature_plots)) {
      combined_plot <- arrange_module_plots(visualization_results$feature_plots, columns = 6)
      print(combined_plot)
    }
  }
  
  # Step 7: Generate analysis report
  generate_eigengene_analysis_report(scored_object, analysis_identifier)
  
  cat("Module eigengene analysis pipeline completed successfully.\n")
  
  return(list(
    seurat_object = scored_object,
    eigengene_matrices = eigengene_matrices,
    visualizations = if (exists("visualization_results")) visualization_results else NULL
  ))
}

#' @title Eigengene Parameter Validation Engine
#' @description Comprehensive validation of module eigengene analysis parameters

validate_eigengene_parameters <- function(seurat_object, harmonization_variable,
                                         cell_grouping_variable, target_cell_group,
                                         scoring_methodology) {
  
  validation_checks <- list()
  
  # Validate Seurat object structure
  if (!inherits(seurat_object, "Seurat")) {
    return(list(is_valid = FALSE, message = "Input must be a valid Seurat object"))
  }
  
  # Validate harmonization variable exists
  if (!is.null(harmonization_variable)) {
    if (!harmonization_variable %in% colnames(seurat_object@meta.data)) {
      validation_checks$harmonization_var <- list(
        is_valid = FALSE,
        message = paste("Harmonization variable", harmonization_variable, "not found in metadata")
      )
    } else {
      validation_checks$harmonization_var <- list(is_valid = TRUE, message = "Harmonization variable validated")
    }
  }
  
  # Validate cell grouping variable exists
  if (!cell_grouping_variable %in% colnames(seurat_object@meta.data)) {
    validation_checks$grouping_var <- list(
      is_valid = FALSE,
      message = paste("Cell grouping variable", cell_grouping_variable, "not found in metadata")
    )
  } else {
    validation_checks$grouping_var <- list(is_valid = TRUE, message = "Grouping variable validated")
  }
  
  # Validate target cell group exists
  if (!target_cell_group %in% unique(seurat_object@meta.data[[cell_grouping_variable]])) {
    available_groups <- paste(unique(seurat_object@meta.data[[cell_grouping_variable]]), collapse = ", ")
    validation_checks$target_group <- list(
      is_valid = FALSE,
      message = paste("Target cell group", target_cell_group, "not found. Available:", available_groups)
    )
  } else {
    validation_checks$target_group <- list(is_valid = TRUE, message = "Target cell group validated")
  }
  
  # Validate scoring methodology
  valid_scoring_methods <- c("UCell", "Seurat")
  if (!scoring_methodology %in% valid_scoring_methods) {
    validation_checks$scoring_method <- list(
      is_valid = FALSE,
      message = paste("Scoring method must be one of:", paste(valid_scoring_methods, collapse = ", "))
    )
  } else {
    validation_checks$scoring_method <- list(is_valid = TRUE, message = "Scoring method validated")
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

#' @title Data Preparation for Eigengene Analysis
#' @description Prepares data by scaling and normalization

prepare_data_for_eigengene_analysis <- function(seurat_object) {
  cat("Preparing data for eigengene analysis...\n")
  
  # Identify variable features if not already done
  if (length(VariableFeatures(seurat_object)) == 0) {
    seurat_object <- FindVariableFeatures(seurat_object, verbose = FALSE)
    cat("Variable features identified:", length(VariableFeatures(seurat_object)), "\n")
  }
  
  # Scale data for improved analysis
  variable_features <- VariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = variable_features, verbose = FALSE)
  
  cat("Data scaling completed for", length(variable_features), "variable features\n")
  return(seurat_object)
}

#' @title Comprehensive Module Eigengene Computer
#' @description Computes eigengenes for all modules with optional harmonization

compute_comprehensive_module_eigengenes <- function(seurat_object,
                                                   harmonization_covariate = NULL,
                                                   analysis_name) {
  
  cat("Computing module eigengenes...\n")
  
  # Retrieve module assignments
  module_assignments <- retrieve_module_assignments(seurat_object)
  unique_modules <- extract_unique_modules(module_assignments)
  
  # Initialize storage for eigengene results
  original_eigengenes <- list()
  harmonized_eigengenes <- list()
  
  # Compute eigengenes for each module
  for (module_id in unique_modules) {
    module_result <- compute_single_module_eigengene(
      seurat_object = seurat_object,
      module_identifier = module_id,
      module_assignments = module_assignments,
      harmonization_variable = harmonization_covariate
    )
    
    original_eigengenes[[module_id]] <- module_result$original_eigengene
    if (!is.null(harmonization_covariate)) {
      harmonized_eigengenes[[module_id]] <- module_result$harmonized_eigengene
    }
  }
  
  # Store results in Seurat object
  seurat_object <- store_eigengene_results(
    seurat_object = seurat_object,
    original_eigengenes = original_eigengenes,
    harmonized_eigengenes = if (length(harmonized_eigengenes) > 0) harmonized_eigengenes else NULL,
    analysis_name = analysis_name
  )
  
  cat("Module eigengene computation completed for", length(unique_modules), "modules\n")
  return(seurat_object)
}

#' @title Module Assignments Retriever
#' @description Retrieves module assignments from Seurat object

retrieve_module_assignments <- function(seurat_object) {
  if (is.null(seurat_object@misc$module_assignments)) {
    stop("Module assignments not found. Please run co-expression network analysis first.")
  }
  return(seurat_object@misc$module_assignments)
}

#' @title Unique Module Extractor
#' @description Extracts unique module identifiers

extract_unique_modules <- function(module_assignments) {
  unique_modules <- unique(module_assignments$module)
  unique_modules <- unique_modules[!is.na(unique_modules) & unique_modules != ""]
  return(unique_modules)
}

#' @title Single Module Eigengene Computer
#' @description Computes eigengene for a specific module

compute_single_module_eigengene <- function(seurat_object, module_identifier,
                                           module_assignments, harmonization_variable = NULL) {
  
  # Extract genes belonging to the module
  module_genes <- module_assignments$gene_name[module_assignments$module == module_identifier]
  module_genes <- module_genes[module_genes %in% rownames(seurat_object)]
  
  if (length(module_genes) == 0) {
    warning("No valid genes found for module: ", module_identifier)
    return(NULL)
  }
  
  # Compute module score (eigengene)
  module_score_name <- paste0("Module_", module_identifier)
  seurat_object <- Seurat::AddModuleScore(
    seurat_object,
    features = list(module_genes),
    name = module_identifier,
    search = FALSE
  )
  
  # Extract original eigengene
  original_eigengene <- seurat_object[[paste0(module_identifier, "1")]][, 1]
  names(original_eigengene) <- colnames(seurat_object)
  
  result <- list(original_eigengene = original_eigengene)
  
  # Apply harmonization if requested
  if (!is.null(harmonization_variable)) {
    harmonized_result <- apply_harmonization(
      seurat_object = seurat_object,
      eigengene_vector = original_eigengene,
      harmonization_variable = harmonization_variable,
      module_identifier = module_identifier
    )
    result$harmonized_eigengene <- harmonized_result
  }
  
  return(result)
}

#' @title Eigengene Harmonization Applicator
#' @description Applies harmony batch correction to eigengenes

apply_harmonization <- function(seurat_object, eigengene_vector, harmonization_variable, module_identifier) {
  # Create temporary Seurat object for harmonization
  temp_object <- CreateSeuratObject(
    counts = matrix(eigengene_vector, nrow = 1, dimnames = list("Eigengene", names(eigengene_vector)))
  )
  temp_object@meta.data <- seurat_object@meta.data[colnames(temp_object), , drop = FALSE]
  
  # Run harmony
  temp_object <- RunPCA(temp_object, features = "Eigengene", verbose = FALSE)
  temp_object <- RunHarmony(
    temp_object,
    group.by.vars = harmonization_variable,
    reduction = "pca",
    reduction.save = "harmony",
    verbose = FALSE
  )
  
  # Extract harmonized eigengene
  harmonized_eigengene <- Embeddings(temp_object, reduction = "harmony")[, 1]
  names(harmonized_eigengene) <- colnames(temp_object)
  
  return(harmonized_eigengene)
}

#' @title Eigengene Results Storage Manager
#' @description Stores computed eigengenes in Seurat object

store_eigengene_results <- function(seurat_object, original_eigengenes, harmonized_eigengenes, analysis_name) {
  # Convert lists to matrices
  original_matrix <- do.call(cbind, original_eigengenes)
  colnames(original_matrix) <- names(original_eigengenes)
  
  seurat_object@misc[[analysis_name]]$original_eigengenes <- original_matrix
  
  if (!is.null(harmonized_eigengenes)) {
    harmonized_matrix <- do.call(cbind, harmonized_eigengenes)
    colnames(harmonized_matrix) <- names(harmonized_eigengenes)
    seurat_object@misc[[analysis_name]]$harmonized_eigengenes <- harmonized_matrix
  }
  
  # Record computation parameters
  seurat_object@misc[[analysis_name]]$eigengene_parameters <- list(
    modules_computed = length(original_eigengenes),
    computation_timestamp = Sys.time(),
    harmonization_applied = !is.null(harmonized_eigengenes)
  )
  
  return(seurat_object)
}

#' @title Eigengene Matrix Extractor
#' @description Extracts eigengene matrices from Seurat object

extract_eigengene_matrices <- function(seurat_object, analysis_name) {
  eigengene_data <- list()
  
  if (!is.null(seurat_object@misc[[analysis_name]]$original_eigengenes)) {
    eigengene_data$original <- seurat_object@misc[[analysis_name]]$original_eigengenes
  }
  
  if (!is.null(seurat_object@misc[[analysis_name]]$harmonized_eigengenes)) {
    eigengene_data$harmonized <- seurat_object@misc[[analysis_name]]$harmonized_eigengenes
  }
  
  cat("Eigengene matrices extracted. Original:", !is.null(eigengene_data$original), 
      "Harmonized:", !is.null(eigengene_data$harmonized), "\n")
  
  return(eigengene_data)
}

#' @title Gene Module Connectivity Calculator
#' @description Calculates gene-module connectivity (kME) metrics

calculate_gene_module_connectivity <- function(seurat_object, grouping_variable,
                                              target_group, analysis_name) {
  
  cat("Calculating gene-module connectivity for", target_group, "cells...\n")
  
  # Filter cells by target group
  target_cells <- filter_cells_by_group(seurat_object, grouping_variable, target_group)
  
  if (length(target_cells) == 0) {
    stop("No cells found for target group: ", target_group)
  }
  
  # Extract expression data and eigengenes
  expression_data <- extract_group_expression_data(seurat_object, target_cells)
  eigengene_data <- extract_group_eigengenes(seurat_object, target_cells, analysis_name)
  
  # Calculate connectivity matrix
  connectivity_matrix <- compute_connectivity_matrix(expression_data, eigengene_data)
  
  # Store results
  seurat_object@misc[[analysis_name]]$gene_connectivity <- list(
    connectivity_matrix = connectivity_matrix,
    target_group = target_group,
    calculation_timestamp = Sys.time()
  )
  
  cat("Gene connectivity calculation completed. Matrix dimensions:", dim(connectivity_matrix), "\n")
  return(seurat_object)
}

#' @title Cell Group Filter
#' @description Filters cells by specified group

filter_cells_by_group <- function(seurat_object, grouping_variable, target_group) {
  group_cells <- rownames(seurat_object@meta.data[
    seurat_object@meta.data[[grouping_variable]] == target_group, ])
  return(group_cells)
}

#' @title Group Expression Data Extractor
#' @description Extracts expression data for specific cell group

extract_group_expression_data <- function(seurat_object, target_cells) {
  expression_data <- GetAssayData(seurat_object, slot = 'data')[, target_cells, drop = FALSE]
  return(as.matrix(expression_data))
}

#' @title Group Eigengene Data Extractor
#' @description Extracts eigengene data for specific cell group

extract_group_eigengenes <- function(seurat_object, target_cells, analysis_name) {
  if (is.null(seurat_object@misc[[analysis_name]]$harmonized_eigengenes)) {
    eigengene_data <- seurat_object@misc[[analysis_name]]$original_eigengenes
  } else {
    eigengene_data <- seurat_object@misc[[analysis_name]]$harmonized_eigengenes
  }
  
  return(as.matrix(eigengene_data[target_cells, , drop = FALSE]))
}

#' @title Connectivity Matrix Computer
#' @description Computes correlation-based connectivity matrix

compute_connectivity_matrix <- function(expression_data, eigengene_data) {
  connectivity_matrix <- cor(t(expression_data), eigengene_data, use = "pairwise.complete.obs")
  rownames(connectivity_matrix) <- rownames(expression_data)
  colnames(connectivity_matrix) <- colnames(eigengene_data)
  
  return(connectivity_matrix)
}

#' @title Cell Module Score Computer
#' @description Computes cell scores based on top module genes

compute_cell_module_scores <- function(seurat_object, top_genes_per_module,
                                      scoring_algorithm, analysis_name) {
  
  cat("Computing cell module scores using", scoring_algorithm, "algorithm...\n")
  
  module_assignments <- retrieve_module_assignments(seurat_object)
  unique_modules <- extract_unique_modules(module_assignments)
  
  # Create gene sets for each module
  module_gene_sets <- create_module_gene_sets(module_assignments, unique_modules, top_genes_per_module)
  
  # Compute module scores
  if (scoring_algorithm == "UCell") {
    scored_object <- UCell::AddModuleScore_UCell(
      seurat_object,
      features = module_gene_sets,
      name = "_module_score"
    )
  } else if (scoring_algorithm == "Seurat") {
    scored_object <- Seurat::AddModuleScore(
      seurat_object,
      features = module_gene_sets,
      name = "ModuleScore"
    )
  }
  
  # Store scoring results
  seurat_object@misc[[analysis_name]]$module_scoring <- list(
    scoring_algorithm = scoring_algorithm,
    top_genes_per_module = top_genes_per_module,
    scoring_timestamp = Sys.time()
  )
  
  cat("Cell module scoring completed for", length(unique_modules), "modules\n")
  return(scored_object)
}

#' @title Module Gene Set Creator
#' @description Creates gene sets for each module based on connectivity

create_module_gene_sets <- function(module_assignments, unique_modules, top_gene_count) {
  gene_sets <- list()
  
  for (module_id in unique_modules) {
    module_genes <- module_assignments[module_assignments$module == module_id, ]
    
    # Order genes by connectivity (kME) if available
    if ("kME" %in% colnames(module_genes)) {
      module_genes <- module_genes[order(-module_genes$kME), ]
    }
    
    # Select top genes
    top_genes <- module_genes$gene_name[1:min(top_gene_count, nrow(module_genes))]
    gene_sets[[paste0("Module_", module_id)]] <- top_genes
  }
  
  return(gene_sets)
}

#' @title Eigengene Visualization Generator
#' @description Creates feature plots for module eigengenes

create_eigengene_visualizations <- function(seurat_object, analysis_name) {
  cat("Generating eigengene visualizations...\n")
  
  visualization_results <- list()
  
  # Create feature plots for harmonized eigengenes if available
  if (!is.null(seurat_object@misc[[analysis_name]]$harmonized_eigengenes)) {
    visualization_results$feature_plots <- create_eigengene_feature_plots(
      seurat_object = seurat_object,
      eigengene_matrix = seurat_object@misc[[analysis_name]]$harmonized_eigengenes,
      plot_title_suffix = " (Harmonized)"
    )
  } else if (!is.null(seurat_object@misc[[analysis_name]]$original_eigengenes)) {
    visualization_results$feature_plots <- create_eigengene_feature_plots(
      seurat_object = seurat_object,
      eigengene_matrix = seurat_object@misc[[analysis_name]]$original_eigengenes,
      plot_title_suffix = " (Original)"
    )
  }
  
  return(visualization_results)
}

#' @title Eigengene Feature Plot Creator
#' @description Creates individual feature plots for each eigengene

create_eigengene_feature_plots <- function(seurat_object, eigengene_matrix, plot_title_suffix = "") {
  feature_plots <- list()
  
  for (eigengene_name in colnames(eigengene_matrix)) {
    # Add eigengene to metadata for plotting
    seurat_object[[eigengene_name]] <- eigengene_matrix[, eigengene_name]
    
    # Create feature plot
    feature_plot <- FeaturePlot(
      seurat_object,
      features = eigengene_name,
      order = TRUE,
      raster = TRUE
    ) + 
      ggtitle(paste0(eigengene_name, plot_title_suffix)) +
      theme(plot.title = element_text(face = "bold", size = 10))
    
    feature_plots[[eigengene_name]] <- feature_plot
  }
  
  return(feature_plots)
}

#' @title Module Plot Arranger
#' @description Arranges module plots in a cohesive layout

arrange_module_plots <- function(plot_list, columns = 6) {
  combined_plot <- wrap_plots(plot_list, ncol = columns) +
    plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(face = "bold")))
  
  return(combined_plot)
}

#' @title Eigengene Analysis Report Generator
#' @description Generates comprehensive report of eigengene analysis results

generate_eigengene_analysis_report <- function(seurat_object, analysis_name) {
  eigengene_params <- seurat_object@misc[[analysis_name]]$eigengene_parameters
  connectivity_data <- seurat_object@misc[[analysis_name]]$gene_connectivity
  scoring_data <- seurat_object@misc[[analysis_name]]$module_scoring
  
  report <- list(
    analysis_summary = list(
      modules_analyzed = eigengene_params$modules_computed,
      harmonization_applied = eigengene_params$harmonization_applied,
      target_group = if (!is.null(connectivity_data)) connectivity_data$target_group else "N/A",
      scoring_algorithm = if (!is.null(scoring_data)) scoring_data$scoring_algorithm else "N/A"
    ),
    computational_metrics = list(
      analysis_timestamp = eigengene_params$computation_timestamp,
      connectivity_matrix_dimensions = if (!is.null(connectivity_data)) dim(connectivity_data$connectivity_matrix) else "N/A"
    )
  )
  
  seurat_object@misc[[analysis_name]]$eigengene_report <- report
  
  cat("Eigengene analysis report:\n")
  cat("  - Modules analyzed:", report$analysis_summary$modules_analyzed, "\n")
  cat("  - Harmonization applied:", report$analysis_summary$harmonization_applied, "\n")
  cat("  - Scoring algorithm:", report$analysis_summary$scoring_algorithm, "\n")
  
  return(seurat_object)
}

# Example usage function
demonstrate_eigengene_analysis <- function() {
  cat("Module eigengene analysis pipeline demonstration\n")
  cat("Please use with actual Seurat object containing module assignments\n")
}

# Uncomment to test (with actual Seurat object)
# results <- execute_module_eigengene_pipeline(your_seurat_object)
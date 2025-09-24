

library(Seurat)
library(WGCNA)
library(ggplot2)
library(patchwork)
library(corrplot)
library(igraph)
library(dplyr)
library(purrr)

#' @title Main Co-expression Network Analysis Pipeline
#' @description Orchestrates the complete gene co-expression network workflow

execute_coexpression_analysis_pipeline <- function(seurat_object,
                                                  target_cell_population = "SC",
                                                  grouping_metadata = "cell_type",
                                                  assay_to_use = "RNA",
                                                  expression_data_layer = "data",
                                                  network_architecture = "signed",
                                                  correlation_methodology = "bicor",
                                                  power_evaluation_sequence = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
                                                  analysis_identifier = NULL,
                                                  enable_metacell_processing = TRUE,
                                                  generate_visualizations = TRUE) {
  
  # Parameter validation and initialization
  validation_result <- validate_coexpression_parameters(
    seurat_object = seurat_object,
    target_cell_population = target_cell_population,
    grouping_metadata = grouping_metadata,
    assay_to_use = assay_to_use,
    network_architecture = network_architecture
  )
  
  if (!validation_result$is_valid) {
    stop("Parameter validation failed: ", validation_result$message)
  }
  
  analysis_identifier <- analysis_identifier %||% seurat_object@misc$active_analysis
  
  cat("Initiating advanced co-expression network analysis pipeline...\n")
  cat("Target population:", target_cell_population, "\n")
  cat("Network type:", network_architecture, "\n")
  
  # Step 1: Prepare expression matrix
  expression_prepared_object <- prepare_expression_matrix(
    seurat_object = seurat_object,
    target_population = target_cell_population,
    grouping_variable = grouping_metadata,
    assay_selection = assay_to_use,
    data_layer = expression_data_layer,
    use_metacells = enable_metacell_processing,
    analysis_name = analysis_identifier
  )
  
  # Step 2: Evaluate optimal soft-thresholding powers
  power_optimized_object <- evaluate_optimal_network_powers(
    seurat_object = expression_prepared_object,
    power_sequence = power_evaluation_sequence,
    network_type = network_architecture,
    correlation_function = correlation_methodology,
    analysis_name = analysis_identifier
  )
  
  # Step 3: Generate power analysis visualizations
  if (generate_visualizations) {
    power_visualizations <- create_power_analysis_plots(
      seurat_object = power_optimized_object,
      analysis_name = analysis_identifier
    )
    
    # Display combined plots
    combined_plot <- arrange_power_plots(power_visualizations, columns = 2)
    print(combined_plot)
  }
  
  # Step 4: Construct co-expression network
  network_constructed_object <- construct_gene_coexpression_network(
    seurat_object = power_optimized_object,
    analysis_name = analysis_identifier,
    network_type = network_architecture
  )
  
  # Step 5: Generate network visualizations
  if (generate_visualizations) {
    create_network_dendrogram(
      seurat_object = network_constructed_object,
      plot_title = paste(target_cell_population, "Gene Co-expression Network Dendrogram"),
      analysis_name = analysis_identifier
    )
  }
  
  # Step 6: Extract topological overlap matrix
  topological_matrix <- extract_topological_overlap_matrix(
    seurat_object = network_constructed_object,
    analysis_name = analysis_identifier
  )
  
  # Step 7: Generate analysis report
  generate_network_analysis_report(network_constructed_object, analysis_identifier)
  
  cat("Co-expression network analysis pipeline completed successfully.\n")
  
  return(list(
    seurat_object = network_constructed_object,
    topological_matrix = topological_matrix,
    power_plots = if (exists("power_visualizations")) power_visualizations else NULL
  ))
}

#' @title Co-expression Parameter Validation Engine
#' @description Comprehensive validation of network analysis parameters

validate_coexpression_parameters <- function(seurat_object, target_cell_population,
                                            grouping_metadata, assay_to_use,
                                            network_architecture) {
  
  validation_checks <- list()
  
  # Validate Seurat object structure
  if (!inherits(seurat_object, "Seurat")) {
    return(list(is_valid = FALSE, message = "Input must be a valid Seurat object"))
  }
  
  # Validate grouping variable exists
  if (!grouping_metadata %in% colnames(seurat_object@meta.data)) {
    validation_checks$grouping_var <- list(
      is_valid = FALSE,
      message = paste("Grouping variable", grouping_metadata, "not found in metadata")
    )
  } else {
    validation_checks$grouping_var <- list(is_valid = TRUE, message = "Grouping variable validated")
  }
  
  # Validate target population exists
  if (!target_cell_population %in% unique(seurat_object@meta.data[[grouping_metadata]])) {
    available_populations <- paste(unique(seurat_object@meta.data[[grouping_metadata]]), collapse = ", ")
    validation_checks$target_pop <- list(
      is_valid = FALSE,
      message = paste("Target population", target_cell_population, "not found. Available:", available_populations)
    )
  } else {
    validation_checks$target_pop <- list(is_valid = TRUE, message = "Target population validated")
  }
  
  # Validate assay exists
  assay_to_use <- assay_to_use %||% DefaultAssay(seurat_object)
  if (!assay_to_use %in% names(seurat_object@assays)) {
    validation_checks$assay <- list(
      is_valid = FALSE,
      message = paste("Assay", assay_to_use, "not found in Seurat object")
    )
  } else {
    validation_checks$assay <- list(is_valid = TRUE, message = "Assay validated")
  }
  
  # Validate network type
  valid_network_types <- c("unsigned", "signed", "signed hybrid")
  if (!network_architecture %in% valid_network_types) {
    validation_checks$network_type <- list(
      is_valid = FALSE,
      message = paste("Network type must be one of:", paste(valid_network_types, collapse = ", "))
    )
  } else {
    validation_checks$network_type <- list(is_valid = TRUE, message = "Network type validated")
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

#' @title Expression Matrix Preparation Engine
#' @description Prepares and validates expression data for network analysis

prepare_expression_matrix <- function(seurat_object, target_population, grouping_variable,
                                     assay_selection, data_layer, use_metacells,
                                     analysis_name) {
  
  cat("Preparing expression matrix for", target_population, "cells...\n")
  
  # Retrieve analysis components
  analysis_genes <- retrieve_analysis_genes(seurat_object, analysis_name)
  metacell_reference <- if (use_metacells) retrieve_metacell_reference(seurat_object, analysis_name) else NULL
  
  # Select appropriate data object
  analysis_object <- select_analysis_data_object(seurat_object, metacell_reference, use_metacells)
  
  # Validate assay availability
  assay_selection <- assay_selection %||% DefaultAssay(analysis_object)
  if (!assay_selection %in% names(analysis_object@assays)) {
    stop("Specified assay '", assay_selection, "' not found in the analysis object")
  }
  
  # Filter cells by target population
  filtered_metadata <- filter_cells_by_population(
    metadata_table = analysis_object@meta.data,
    grouping_variable = grouping_variable,
    target_population = target_population
  )
  
  selected_cells <- rownames(filtered_metadata)
  
  # Extract expression data
  expression_data <- extract_expression_data(
    analysis_object = analysis_object,
    assay_name = assay_selection,
    data_slot = data_layer,
    gene_set = analysis_genes,
    cell_set = selected_cells
  )
  
  # Validate gene expression quality
  validated_genes <- validate_expression_quality(expression_data)
  filtered_expression <- expression_data[, validated_genes, drop = FALSE]
  
  cat("Expression matrix prepared. Dimensions:", dim(filtered_expression), "\n")
  
  # Store results in Seurat object
  seurat_object <- store_expression_matrix(
    seurat_object = seurat_object,
    expression_matrix = filtered_expression,
    validated_genes = validated_genes,
    analysis_name = analysis_name
  )
  
  return(seurat_object)
}

#' @title Analysis Data Object Selector
#' @description Chooses between original Seurat object and metacell reference

select_analysis_data_object <- function(seurat_object, metacell_reference, use_metacells) {
  if (use_metacells && !is.null(metacell_reference)) {
    cat("Using metacell reference for analysis\n")
    return(metacell_reference)
  } else {
    cat("Using original Seurat object for analysis\n")
    return(seurat_object)
  }
}

#' @title Cell Population Filter
#' @description Filters metadata to select specific cell population

filter_cells_by_population <- function(metadata_table, grouping_variable, target_population) {
  if (!grouping_variable %in% colnames(metadata_table)) {
    stop("Grouping variable '", grouping_variable, "' not found in metadata")
  }
  
  filtered_metadata <- metadata_table[metadata_table[[grouping_variable]] %in% target_population, , drop = FALSE]
  
  if (nrow(filtered_metadata) == 0) {
    stop("No cells found for target population '", target_population, "'")
  }
  
  cat("Selected", nrow(filtered_metadata), "cells for population", target_population, "\n")
  return(filtered_metadata)
}

#' @title Expression Data Extractor
#' @description Extracts expression data with version compatibility

extract_expression_data <- function(analysis_object, assay_name, data_slot, gene_set, cell_set) {
  # Extract assay data with Seurat version compatibility
  assay_data <- GetAssayData(analysis_object, assay = assay_name, slot = data_slot)
  
  # Ensure genes and cells are available
  available_genes <- intersect(gene_set, rownames(assay_data))
  available_cells <- intersect(cell_set, colnames(assay_data))
  
  if (length(available_genes) == 0) {
    stop("No selected genes found in the expression data")
  }
  
  if (length(available_cells) == 0) {
    stop("No selected cells found in the expression data")
  }
  
  # Extract and transpose expression matrix
  expression_matrix <- t(as.matrix(assay_data[available_genes, available_cells, drop = FALSE]))
  
  return(expression_matrix)
}

#' @title Expression Quality Validator
#' @description Identifies genes with sufficient expression variation

validate_expression_quality <- function(expression_matrix) {
  # Use WGCNA's goodGenes function to identify valid genes
  valid_gene_indices <- WGCNA::goodGenes(expression_matrix, verbose = 0)
  valid_gene_names <- colnames(expression_matrix)[valid_gene_indices]
  
  removed_genes <- ncol(expression_matrix) - length(valid_gene_names)
  if (removed_genes > 0) {
    cat("Removed", removed_genes, "genes with poor expression characteristics\n")
  }
  
  return(valid_gene_names)
}

#' @title Expression Matrix Storage Manager
#' @description Stores prepared expression matrix in Seurat object

store_expression_matrix <- function(seurat_object, expression_matrix, validated_genes, analysis_name) {
  seurat_object@misc[[analysis_name]]$expression_matrix <- expression_matrix
  seurat_object@misc[[analysis_name]]$validated_genes <- validated_genes
  
  # Record preparation parameters
  seurat_object@misc[[analysis_name]]$preparation_parameters <- list(
    matrix_dimensions = dim(expression_matrix),
    preparation_timestamp = Sys.time(),
    genes_retained = length(validated_genes)
  )
  
  return(seurat_object)
}

#' @title Optimal Power Evaluation Engine
#' @description Determines optimal soft-thresholding power for network construction

evaluate_optimal_network_powers <- function(seurat_object, power_sequence, network_type,
                                           correlation_function, analysis_name) {
  
  cat("Evaluating optimal soft-thresholding powers...\n")
  
  expression_matrix <- seurat_object@misc[[analysis_name]]$expression_matrix
  
  # Perform power analysis
  power_analysis_results <- WGCNA::pickSoftThreshold(
    expression_matrix,
    powerVector = power_sequence,
    verbose = 100,
    networkType = network_type,
    corFnc = correlation_function
  )
  
  # Determine optimal power
  optimal_power <- determine_optimal_power(power_analysis_results$fitIndices)
  
  # Store results
  seurat_object@misc[[analysis_name]]$power_analysis <- list(
    fit_indices = power_analysis_results$fitIndices,
    optimal_power = optimal_power,
    analysis_timestamp = Sys.time()
  )
  
  cat("Optimal power determined:", optimal_power, "\n")
  return(seurat_object)
}

#' @title Optimal Power Determiner
#' @description Identifies the best soft-thresholding power

determine_optimal_power <- function(fit_indices) {
  # Find the lowest power that achieves scale-free topology fit R² ≥ 0.8
  satisfactory_powers <- fit_indices$Power[fit_indices$SFT.R.sq >= 0.8]
  
  if (length(satisfactory_powers) > 0) {
    optimal_power <- min(satisfactory_powers)
  } else {
    # If no power reaches 0.8, choose the power with highest R²
    optimal_power <- fit_indices$Power[which.max(fit_indices$SFT.R.sq)]
    warning("No power achieved R² ≥ 0.8. Using power with highest R²: ", optimal_power)
  }
  
  return(optimal_power)
}

#' @title Power Analysis Visualization Generator
#' @description Creates diagnostic plots for power analysis

create_power_analysis_plots <- function(seurat_object, analysis_name) {
  power_data <- seurat_object@misc[[analysis_name]]$power_analysis$fit_indices
  optimal_power <- seurat_object@misc[[analysis_name]]$power_analysis$optimal_power
  
  # Create scale-free topology fit plot
  topology_plot <- ggplot(power_data, aes(x = Power, y = SFT.R.sq)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line(color = "blue", linewidth = 1) +
    geom_vline(xintercept = optimal_power, linetype = "dashed", color = "red", linewidth = 1) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "darkgreen", linewidth = 0.8) +
    labs(
      title = "Scale-Free Topology Model Fit",
      x = "Soft-Thresholding Power",
      y = "Scale-Free Topology Fit (R²)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # Create mean connectivity plot
  connectivity_plot <- ggplot(power_data, aes(x = Power, y = mean.k.)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_line(color = "purple", linewidth = 1) +
    geom_vline(xintercept = optimal_power, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "Mean Gene Connectivity",
      x = "Soft-Thresholding Power",
      y = "Mean Connectivity"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(list(topology_plot = topology_plot, connectivity_plot = connectivity_plot))
}

#' @title Power Plot Arranger
#' @description Arranges power analysis plots in a cohesive layout

arrange_power_plots <- function(plot_list, columns = 2) {
  combined_plot <- wrap_plots(plot_list, ncol = columns) +
    plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(face = "bold")))
  
  return(combined_plot)
}

#' @title Co-expression Network Constructor
#' @description Builds gene co-expression network using WGCNA

construct_gene_coexpression_network <- function(seurat_object, analysis_name, network_type) {
  cat("Constructing gene co-expression network...\n")
  
  expression_matrix <- seurat_object@misc[[analysis_name]]$expression_matrix
  optimal_power <- seurat_object@misc[[analysis_name]]$power_analysis$optimal_power
  
  # Construct network modules
  network_modules <- WGCNA::blockwiseModules(
    expression_matrix,
    power = optimal_power,
    TOMType = network_type,
    numericLabels = TRUE,
    saveTOMs = TRUE,
    saveTOMFileBase = paste0("network_TOM_", analysis_name),
    verbose = 3
  )
  
  # Store network results
  seurat_object@misc[[analysis_name]]$network_construction <- list(
    module_assignment = network_modules,
    construction_power = optimal_power,
    network_type = network_type,
    construction_timestamp = Sys.time()
  )
  
  cat("Network construction completed. Identified", length(unique(network_modules$colors)), "modules\n")
  return(seurat_object)
}

#' @title Network Dendrogram Visualizer
#' @description Creates dendrogram visualization of gene modules

create_network_dendrogram <- function(seurat_object, plot_title, analysis_name) {
  network_data <- seurat_object@misc[[analysis_name]]$network_construction$module_assignment
  
  # Create dendrogram with module colors
  WGCNA::plotDendroAndColors(
    network_data$dendrograms[[1]],
    network_data$colors,
    dendroLabels = FALSE,
    main = plot_title,
    addGuide = TRUE,
    guideHang = 0.05
  )
}

#' @title Topological Overlap Matrix Extractor
#' @description Extracts and returns the topological overlap matrix

extract_topological_overlap_matrix <- function(seurat_object, analysis_name) {
  network_data <- seurat_object@misc[[analysis_name]]$network_construction$module_assignment
  
  # Calculate topological overlap matrix
  topological_matrix <- WGCNA::TOMsimilarityFromExpr(
    network_data$datExpr,
    power = seurat_object@misc[[analysis_name]]$network_construction$construction_power
  )
  
  # Set row and column names
  rownames(topological_matrix) <- colnames(topological_matrix) <- colnames(network_data$datExpr)
  
  cat("Topological overlap matrix extracted. Dimensions:", dim(topological_matrix), "\n")
  return(topological_matrix)
}

#' @title Network Analysis Report Generator
#' @description Generates comprehensive report of network analysis results

generate_network_analysis_report <- function(seurat_object, analysis_name) {
  network_info <- seurat_object@misc[[analysis_name]]$network_construction
  power_info <- seurat_object@misc[[analysis_name]]$power_analysis
  prep_info <- seurat_object@misc[[analysis_name]]$preparation_parameters
  
  report <- list(
    analysis_summary = list(
      genes_analyzed = prep_info$genes_retained,
      optimal_power = power_info$optimal_power,
      modules_identified = length(unique(network_info$module_assignment$colors)),
      network_type = network_info$network_type
    ),
    quality_metrics = list(
      scale_free_fit = max(power_info$fit_indices$SFT.R.sq),
      mean_connectivity = power_info$fit_indices$mean.k.[power_info$fit_indices$Power == power_info$optimal_power]
    )
  )
  
  seurat_object@misc[[analysis_name]]$analysis_report <- report
  
  cat("Network analysis report:\n")
  cat("  - Genes analyzed:", report$analysis_summary$genes_analyzed, "\n")
  cat("  - Optimal power:", report$analysis_summary$optimal_power, "\n")
  cat("  - Modules identified:", report$analysis_summary$modules_identified, "\n")
  cat("  - Scale-free fit R²:", round(report$quality_metrics$scale_free_fit, 3), "\n")
  
  return(seurat_object)
}

# Example usage function
demonstrate_coexpression_analysis <- function() {
  cat("Co-expression network analysis pipeline demonstration\n")
  cat("Please use with actual Seurat object containing expression data\n")
}

# Uncomment to test (with actual Seurat object)
# results <- execute_coexpression_analysis_pipeline(your_seurat_object)
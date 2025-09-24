

library(Seurat)
library(Matrix)
library(FNN)
library(dplyr)
library(purrr)

#' @title Main Metacell Construction Pipeline
#' @description Orchestrates the complete metacell generation workflow

execute_metacell_construction_pipeline <- function(seurat_object,
                                                  grouping_variables = c("cell_type", "Sample"),
                                                  dimensionality_reduction = "harmony",
                                                  neighbor_count = 25,
                                                  maximum_overlap = 10,
                                                  identity_variable = "cell_type",
                                                  assay_selection = NULL,
                                                  data_slot_selection = "counts",
                                                  aggregation_method = "average",
                                                  minimum_group_size = 100,
                                                  target_metacell_count = 1000,
                                                  maximum_iterations = 5000,
                                                  analysis_identifier = NULL,
                                                  enable_verbose_logging = FALSE) {
  
  # Parameter validation and initialization
  validation_result <- validate_metacell_parameters(
    seurat_object = seurat_object,
    grouping_variables = grouping_variables,
    dimensionality_reduction = dimensionality_reduction,
    identity_variable = identity_variable,
    assay_selection = assay_selection,
    data_slot_selection = data_slot_selection
  )
  
  if (!validation_result$is_valid) {
    stop("Parameter validation failed: ", validation_result$message)
  }
  
  analysis_identifier <- analysis_identifier %||% seurat_object@misc$active_analysis
  
  cat("Initiating advanced metacell construction pipeline...\n")
  cat("Grouping variables:", paste(grouping_variables, collapse = ", "), "\n")
  cat("Dimensionality reduction:", dimensionality_reduction, "\n")
  
  # Step 1: Preprocess and group cells
  grouped_object <- preprocess_cell_grouping(
    seurat_object = seurat_object,
    grouping_variables = grouping_variables,
    identity_variable = identity_variable,
    minimum_group_size = minimum_group_size
  )
  
  # Step 2: Construct metacells for each group
  metacell_results <- construct_group_metacells(
    seurat_object = grouped_object$filtered_object,
    cell_groups = grouped_object$valid_groups,
    dimensionality_reduction = dimensionality_reduction,
    neighbor_count = neighbor_count,
    maximum_overlap = maximum_overlap,
    assay_selection = assay_selection,
    data_slot_selection = data_slot_selection,
    aggregation_method = aggregation_method,
    target_metacell_count = target_metacell_count,
    maximum_iterations = maximum_iterations,
    analysis_identifier = analysis_identifier,
    enable_verbose_logging = enable_verbose_logging
  )
  
  # Step 3: Integrate metacell results
  integrated_object <- integrate_metacell_results(
    seurat_object = seurat_object,
    metacell_objects = metacell_results$metacell_list,
    grouping_variables = grouping_variables,
    identity_variable = identity_variable,
    analysis_identifier = analysis_identifier
  )
  
  # Step 4: Normalize metacell expression
  normalized_object <- normalize_metacell_expression(
    seurat_object = integrated_object,
    analysis_identifier = analysis_identifier
  )
  
  # Step 5: Generate construction report
  generate_metacell_construction_report(integrated_object, analysis_identifier)
  
  cat("Metacell construction pipeline completed successfully.\n")
  return(normalized_object)
}

#' @title Metacell Parameter Validation Engine
#' @description Comprehensive validation of all metacell construction parameters

validate_metacell_parameters <- function(seurat_object, grouping_variables, 
                                        dimensionality_reduction, identity_variable,
                                        assay_selection, data_slot_selection) {
  
  validation_checks <- list()
  
  # Validate Seurat object structure
  if (!inherits(seurat_object, "Seurat")) {
    return(list(is_valid = FALSE, message = "Input must be a valid Seurat object"))
  }
  
  # Validate grouping variables exist in metadata
  missing_group_vars <- setdiff(grouping_variables, colnames(seurat_object@meta.data))
  if (length(missing_group_vars) > 0) {
    validation_checks$group_vars <- list(
      is_valid = FALSE,
      message = paste("Missing grouping variables:", paste(missing_group_vars, collapse = ", "))
    )
  } else {
    validation_checks$group_vars <- list(is_valid = TRUE, message = "Grouping variables validated")
  }
  
  # Validate identity variable inclusion
  if (!identity_variable %in% grouping_variables) {
    validation_checks$identity_var <- list(
      is_valid = FALSE,
      message = paste("Identity variable", identity_variable, "must be included in grouping variables")
    )
  } else {
    validation_checks$identity_var <- list(is_valid = TRUE, message = "Identity variable validated")
  }
  
  # Validate dimensionality reduction exists
  if (!dimensionality_reduction %in% names(seurat_object@reductions)) {
    available_reductions <- paste(names(seurat_object@reductions), collapse = ", ")
    validation_checks$reduction <- list(
      is_valid = FALSE,
      message = paste("Reduction", dimensionality_reduction, "not found. Available:", available_reductions)
    )
  } else {
    validation_checks$reduction <- list(is_valid = TRUE, message = "Dimensionality reduction validated")
  }
  
  # Validate assay selection
  assay_selection <- assay_selection %||% DefaultAssay(seurat_object)
  if (!assay_selection %in% names(seurat_object@assays)) {
    validation_checks$assay <- list(
      is_valid = FALSE,
      message = paste("Assay", assay_selection, "not found in Seurat object")
    )
  } else {
    validation_checks$assay <- list(is_valid = TRUE, message = "Assay validated")
  }
  
  # Validate data slot selection
  valid_slots <- c("counts", "data", "scale.data")
  if (!data_slot_selection %in% valid_slots) {
    validation_checks$data_slot <- list(
      is_valid = FALSE,
      message = paste("Data slot must be one of:", paste(valid_slots, collapse = ", "))
    )
  } else {
    validation_checks$data_slot <- list(is_valid = TRUE, message = "Data slot validated")
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

#' @title Cell Grouping Preprocessor
#' @description Organizes cells into groups based on metadata variables

preprocess_cell_grouping <- function(seurat_object, grouping_variables, 
                                    identity_variable, minimum_group_size) {
  
  cat("Preprocessing cell groups...\n")
  
  # Extract and clean metadata
  metadata_df <- seurat_object@meta.data[, grouping_variables, drop = FALSE]
  metadata_df[] <- lapply(metadata_df, as.character)
  
  # Create unique group identifiers
  group_identifiers <- apply(metadata_df, 1, function(x) paste(x, collapse = "##"))
  seurat_object$metacell_group_id <- group_identifiers
  
  # Calculate group sizes and filter
  group_size_table <- table(seurat_object$metacell_group_id)
  valid_group_ids <- names(group_size_table[group_size_table >= minimum_group_size])
  
  if (length(valid_group_ids) == 0) {
    stop("No cell groups meet the minimum size requirement of ", minimum_group_size, " cells")
  }
  
  cat("Identified", length(valid_group_ids), "valid cell groups\n")
  
  # Create subset object for valid groups
  valid_cell_indices <- which(seurat_object$metacell_group_id %in% valid_group_ids)
  filtered_object <- seurat_object[, valid_cell_indices]
  
  return(list(
    filtered_object = filtered_object,
    valid_groups = valid_group_ids,
    group_sizes = group_size_table[valid_group_ids]
  ))
}

#' @title Group-Specific Metacell Constructor
#' @description Builds metacells for each cell group independently

construct_group_metacells <- function(seurat_object, cell_groups,
                                     dimensionality_reduction, neighbor_count,
                                     maximum_overlap, assay_selection, data_slot_selection,
                                     aggregation_method, target_metacell_count,
                                     maximum_iterations, analysis_identifier,
                                     enable_verbose_logging) {
  
  cat("Constructing metacells for", length(cell_groups), "cell groups...\n")
  
  metacell_objects <- list()
  construction_statistics <- list()
  
  for (group_id in cell_groups) {
    if (enable_verbose_logging) {
      cat("Processing group:", group_id, "\n")
    }
    
    # Extract subgroup data
    subgroup_cells <- which(seurat_object$metacell_group_id == group_id)
    subgroup_object <- seurat_object[, subgroup_cells]
    
    # Construct metacells for this group
    metacell_result <- build_single_group_metacells(
      seurat_object = subgroup_object,
      group_identifier = group_id,
      dimensionality_reduction = dimensionality_reduction,
      neighbor_count = neighbor_count,
      maximum_overlap = maximum_overlap,
      assay_selection = assay_selection,
      data_slot_selection = data_slot_selection,
      aggregation_method = aggregation_method,
      target_metacell_count = target_metacell_count,
      maximum_iterations = maximum_iterations
    )
    
    if (!is.null(metacell_result)) {
      metacell_objects[[group_id]] <- metacell_result$metacell_object
      construction_statistics[[group_id]] <- metacell_result$construction_stats
    }
  }
  
  # Remove NULL results
  metacell_objects <- metacell_objects[!sapply(metacell_objects, is.null)]
  
  if (length(metacell_objects) == 0) {
    stop("Metacell construction failed for all groups")
  }
  
  cat("Successfully constructed metacells for", length(metacell_objects), "groups\n")
  
  return(list(
    metacell_list = metacell_objects,
    group_statistics = construction_statistics
  ))
}

#' @title Single Group Metacell Builder
#' @description Constructs metacells for a specific cell group

build_single_group_metacells <- function(seurat_object, group_identifier,
                                        dimensionality_reduction, neighbor_count,
                                        maximum_overlap, assay_selection, data_slot_selection,
                                        aggregation_method, target_metacell_count,
                                        maximum_iterations) {
  
  tryCatch({
    # Extract dimensionality reduction coordinates
    reduction_coordinates <- Embeddings(seurat_object[[dimensionality_reduction]])
    
    # Calculate k-nearest neighbors
    neighbor_graph <- find_k_nearest_neighbors(
      coordinates = reduction_coordinates,
      k_neighbors = neighbor_count
    )
    
    # Perform metacell assignment
    metacell_assignment <- assign_metacell_membership(
      neighbor_graph = neighbor_graph,
      max_overlap = maximum_overlap,
      target_metacell_count = target_metacell_count,
      max_iterations = maximum_iterations
    )
    
    # Aggregate expression data
    aggregated_expression <- aggregate_metacell_expression(
      seurat_object = seurat_object,
      metacell_assignments = metacell_assignment,
      assay_name = assay_selection,
      data_slot = data_slot_selection,
      aggregation_mode = aggregation_method
    )
    
    # Create metacell Seurat object
    metacell_object <- create_metacell_seurat_object(
      expression_matrix = aggregated_expression$expression,
      cell_assignments = aggregated_expression$assignments,
      original_metadata = seurat_object@meta.data,
      group_identifier = group_identifier
    )
    
    # Calculate construction statistics
    stats <- calculate_metacell_statistics(
      original_cell_count = ncol(seurat_object),
      metacell_count = ncol(metacell_object),
      assignment_results = metacell_assignment
    )
    
    return(list(
      metacell_object = metacell_object,
      construction_stats = stats
    ))
    
  }, error = function(e) {
    warning("Metacell construction failed for group ", group_identifier, ": ", e$message)
    return(NULL)
  })
}

#' @title K-Nearest Neighbors Calculator
#' @description Computes nearest neighbors using efficient algorithms

find_k_nearest_neighbors <- function(coordinates, k_neighbors) {
  # Use fast nearest neighbor algorithm
  knn_result <- FNN::get.knn(coordinates, k = k_neighbors)
  
  # Convert to adjacency list format
  neighbor_list <- lapply(1:nrow(coordinates), function(i) {
    knn_result$nn.index[i, ]
  })
  
  return(neighbor_list)
}

#' @title Metacell Membership Assigner
#' @description Assigns cells to metacells with overlap constraints

assign_metacell_membership <- function(neighbor_graph, max_overlap, 
                                      target_metacell_count, max_iterations) {
  
  total_cells <- length(neighbor_graph)
  metacell_assignments <- vector("list", target_metacell_count)
  cell_assignments <- rep(NA, total_cells)
  assigned_cells <- logical(total_cells)
  
  iteration <- 0
  while (sum(assigned_cells) < total_cells && iteration < max_iterations) {
    iteration <- iteration + 1
    
    # Find unassigned cell with most unassigned neighbors
    best_cell <- find_optimal_seed_cell(neighbor_graph, assigned_cells)
    
    if (is.na(best_cell)) break
    
    # Create metacell around seed cell
    metacell_members <- create_metacell_from_seed(
      seed_cell = best_cell,
      neighbor_graph = neighbor_graph,
      assigned_cells = assigned_cells,
      max_overlap = max_overlap
    )
    
    if (length(metacell_members) > 0) {
      metacell_id <- which(sapply(metacell_assignments, is.null))[1]
      if (is.na(metacell_id)) break
      
      metacell_assignments[[metacell_id]] <- metacell_members
      assigned_cells[metacell_members] <- TRUE
    }
  }
  
  # Remove empty metacells
  metacell_assignments <- metacell_assignments[!sapply(metacell_assignments, is.null)]
  
  return(list(
    assignments = metacell_assignments,
    coverage = sum(assigned_cells) / total_cells,
    iterations = iteration
  ))
}

#' @title Optimal Seed Cell Finder
#' @description Identifies the best cell to seed a new metacell

find_optimal_seed_cell <- function(neighbor_graph, assigned_cells) {
  unassigned_cells <- which(!assigned_cells)
  
  if (length(unassigned_cells) == 0) return(NA)
  
  # Score cells by number of unassigned neighbors
  cell_scores <- sapply(unassigned_cells, function(cell) {
    neighbors <- neighbor_graph[[cell]]
    sum(!assigned_cells[neighbors])
  })
  
  best_cell <- unassigned_cells[which.max(cell_scores)]
  return(best_cell)
}

#' @title Metacell Creator from Seed
#' @description Forms a metacell around a seed cell

create_metacell_from_seed <- function(seed_cell, neighbor_graph, assigned_cells, max_overlap) {
  seed_neighbors <- neighbor_graph[[seed_cell]]
  available_neighbors <- seed_neighbors[!assigned_cells[seed_neighbors]]
  
  # Include seed and available neighbors
  metacell_members <- c(seed_cell, available_neighbors)
  
  return(metacell_members)
}

#' @title Metacell Expression Aggregator
#' @description Aggregates expression data for metacell construction

aggregate_metacell_expression <- function(seurat_object, metacell_assignments,
                                         assay_name, data_slot, aggregation_mode) {
  
  # Extract expression data
  expression_data <- GetAssayData(seurat_object, slot = data_slot, assay = assay_name)
  
  # Aggregate expression for each metacell
  aggregated_matrix <- matrix(0, nrow = nrow(expression_data), 
                             ncol = length(metacell_assignments$assignments))
  
  rownames(aggregated_matrix) <- rownames(expression_data)
  colnames(aggregated_matrix) <- paste0("Metacell_", 1:length(metacell_assignments$assignments))
  
  for (i in 1:length(metacell_assignments$assignments)) {
    cell_indices <- metacell_assignments$assignments[[i]]
    
    if (aggregation_mode == "average") {
      aggregated_matrix[, i] <- Matrix::rowMeans(expression_data[, cell_indices, drop = FALSE])
    } else if (aggregation_mode == "sum") {
      aggregated_matrix[, i] <- Matrix::rowSums(expression_data[, cell_indices, drop = FALSE])
    }
  }
  
  return(list(
    expression = aggregated_matrix,
    assignments = metacell_assignments
  ))
}

#' @title Metacell Seurat Object Creator
#' @description Creates a Seurat object from aggregated metacell data

create_metacell_seurat_object <- function(expression_matrix, cell_assignments,
                                         original_metadata, group_identifier) {
  
  # Create new Seurat object
  metacell_object <- CreateSeuratObject(
    counts = expression_matrix,
    project = paste("Metacell", group_identifier)
  )
  
  # Add group information to metadata
  metacell_object$original_group <- group_identifier
  metacell_object$metacell_size <- sapply(cell_assignments$assignments, length)
  
  return(metacell_object)
}

#' @title Metacell Integration Manager
#' @description Integrates metacell objects from different groups

integrate_metacell_results <- function(seurat_object, metacell_objects,
                                      grouping_variables, identity_variable,
                                      analysis_identifier) {
  
  cat("Integrating metacell objects...\n")
  
  if (length(metacell_objects) == 0) {
    stop("No metacell objects to integrate")
  }
  
  # Merge all metacell objects
  combined_metacell <- reduce(metacell_objects, function(x, y) {
    merge(x, y, add.cell.ids = c(deparse(substitute(x)), deparse(substitute(y))))
  })
  
  # Set identity based on specified variable
  Idents(combined_metacell) <- identity_variable
  
  # Store in original Seurat object
  seurat_object@misc[[analysis_identifier]]$metacell_reference <- combined_metacell
  
  # Record integration parameters
  seurat_object@misc[[analysis_identifier]]$integration_parameters <- list(
    integrated_groups = length(metacell_objects),
    total_metacells = ncol(combined_metacell),
    integration_timestamp = Sys.time()
  )
  
  return(seurat_object)
}

#' @title Metacell Expression Normalizer
#' @description Applies normalization to metacell expression data

normalize_metacell_expression <- function(seurat_object, analysis_identifier) {
  metacell_object <- seurat_object@misc[[analysis_identifier]]$metacell_reference
  
  # Apply standard Seurat normalization
  normalized_metacell <- NormalizeData(metacell_object, verbose = FALSE)
  
  # Update stored object
  seurat_object@misc[[analysis_identifier]]$metacell_reference <- normalized_metacell
  
  cat("Metacell expression normalization completed\n")
  return(seurat_object)
}

#' @title Metacell Statistics Calculator
#' @description Computes detailed statistics about metacell construction

calculate_metacell_statistics <- function(original_cell_count, metacell_count, assignment_results) {
  return(list(
    original_cells = original_cell_count,
    resulting_metacells = metacell_count,
    compression_ratio = original_cell_count / metacell_count,
    cell_coverage = assignment_results$coverage,
    iterations_used = assignment_results$iterations
  ))
}

#' @title Metacell Construction Report Generator
#' @description Generates comprehensive QC report for metacell construction

generate_metacell_construction_report <- function(seurat_object, analysis_identifier) {
  metacell_data <- seurat_object@misc[[analysis_identifier]]$metacell_reference
  integration_params <- seurat_object@misc[[analysis_identifier]]$integration_parameters
  
  report <- list(
    construction_summary = list(
      total_metacells_created = integration_params$total_metacells,
      groups_integrated = integration_params$integrated_groups,
      average_metacell_size = mean(metacell_data$metacell_size),
      construction_timestamp = integration_params$integration_timestamp
    ),
    quality_metrics = list(
      median_genes_per_metacell = median(Matrix::colSums(GetAssayData(metacell_data) > 0)),
      total_genes_detected = nrow(metacell_data)
    )
  )
  
  seurat_object@misc[[analysis_identifier]]$construction_report <- report
  
  cat("Metacell construction report:\n")
  cat("  - Total metacells:", report$construction_summary$total_metacells_created, "\n")
  cat("  - Average metacell size:", round(report$construction_summary$average_metacell_size, 1), "\n")
  cat("  - Genes detected:", report$quality_metrics$total_genes_detected, "\n")
  
  return(seurat_object)
}

# Example usage function
demonstrate_metacell_pipeline <- function() {
  # This would use actual Seurat object in practice
  cat("Metacell construction pipeline demonstration\n")
  cat("Please use with actual Seurat object containing cell metadata\n")
}

# Uncomment to test (with actual Seurat object)
# results <- execute_metacell_construction_pipeline(your_seurat_object)
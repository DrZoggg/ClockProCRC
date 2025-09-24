#' @title Advanced ClockProCRC Biological Rhythm Analysis System
#' @description Comprehensive computational framework for circadian rhythm quantification
#' @author Biomedical Informatics Research Group

library(arules)
library(matrixStats)
library(foreach)
library(doParallel)

#' @title Randomized Background Model Generator
#' @description Creates permutation-based null models for circadian rhythm analysis
#' @param expression_matrix Normalized gene expression data
#' @param gene_stratification Vector indicating gene stratification groups
#' @param circadian_gene_indices Logical vector indicating circadian genes
#' @param permutation_iterations Number of random permutations (default: 1000)
#' @param parallel_computation Whether to use parallel processing (default: TRUE)

generate_randomized_background_model <- function(expression_matrix, 
                                                gene_stratification, 
                                                circadian_gene_indices,
                                                permutation_iterations = 1000,
                                                parallel_computation = TRUE) {
  
  # Input validation
  if (sum(circadian_gene_indices) < 5) {
    stop("Insufficient circadian genes for reliable background modeling")
  }
  
  # Calculate gene counts per stratification bin
  bin_distribution <- table(gene_stratification[circadian_gene_indices])
  active_bins <- names(bin_distribution)[bin_distribution > 0]
  
  # Initialize parallel processing if requested
  if (parallel_computation && require(doParallel)) {
    core_count <- min(detectCores() - 1, permutation_iterations %/% 100)
    cl <- makeCluster(core_count)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  # Generate randomized background scores
  randomized_scores <- foreach(perm_idx = 1:permutation_iterations, 
                              .combine = cbind,
                              .packages = c("matrixStats")) %dopar% {
    
    background_selection <- logical(length(gene_stratification))
    
    for (bin_name in active_bins) {
      bin_genes <- which(gene_stratification == bin_name)
      required_count <- bin_distribution[bin_name]
      
      if (length(bin_genes) >= required_count) {
        selected_genes <- sample(bin_genes, required_count)
        background_selection[selected_genes] <- TRUE
      }
    }
    
    # Calculate background scores
    if (sum(background_selection) > 0) {
      colMeans(expression_matrix[background_selection, , drop = FALSE])
    } else {
      rep(0, ncol(expression_matrix))
    }
  }
  
  # Calculate robust average across permutations
  background_means <- rowMeans(randomized_scores, na.rm = TRUE)
  background_variance <- rowSds(randomized_scores, na.rm = TRUE)
  
  return(list(
    mean_scores = background_means,
    standard_errors = background_variance,
    full_distribution = randomized_scores
  ))
}

#' @title Expression Data Preprocessor
#' @description Handles data normalization and quality control

preprocess_expression_data <- function(expression_data, analysis_type) {
  
  # Validate input data
  if (!is.matrix(expression_data)) {
    expression_data <- as.matrix(expression_data)
  }
  
  if (any(is.infinite(expression_data) | is.nan(expression_data))) {
    warning("Infinite or NaN values detected. Applying data cleaning.")
    expression_data[!is.finite(expression_data)] <- NA
  }
  
  processed_data <- list(
    raw_expression = expression_data,
    gene_identifiers = rownames(expression_data)
  )
  
  # Analysis-type specific processing
  if (analysis_type == "single_cell") {
    # Single-cell specific processing
    gene_means <- rowMeans(expression_data, na.rm = TRUE)
    processed_data$normalized_expression <- expression_data - gene_means
    
    # Handle zero-inflation in scRNA-seq
    adjusted_data <- 10 * (2^expression_data - 1)
    processed_data$expression_magnitude <- log2(rowMeans(adjusted_data, na.rm = TRUE) + 1)
    
  } else if (analysis_type == "bulk_sequencing") {
    # Bulk RNA-seq processing
    processed_data$expression_magnitude <- rowMeans(expression_data, na.rm = TRUE)
    processed_data$normalized_expression <- expression_data - processed_data$expression_magnitude
    
  } else {
    stop("Unsupported analysis type. Use 'single_cell' or 'bulk_sequencing'")
  }
  
  # Quality control metrics
  processed_data$quality_metrics <- list(
    genes_passing_qc = sum(complete.cases(expression_data)),
    median_expression = median(processed_data$expression_magnitude, na.rm = TRUE),
    detection_rate = mean(expression_data > 0, na.rm = TRUE)
  )
  
  return(processed_data)
}

#' @title Gene Stratification Engine
#' @description Implements adaptive binning strategies for gene expression

implement_gene_stratification <- function(expression_magnitude, 
                                         stratification_method = "adaptive_frequency",
                                         number_of_bins = 50) {
  
  available_methods <- c("equal_frequency", "equal_width", "adaptive_frequency")
  if (!stratification_method %in% available_methods) {
    stop("Invalid stratification method. Choose from: ", paste(available_methods, collapse = ", "))
  }
  
  # Remove NA values for discretization
  valid_expression <- expression_magnitude[is.finite(expression_magnitude)]
  
  if (stratification_method == "adaptive_frequency") {
    # Adaptive binning based on data distribution
    quantile_values <- quantile(valid_expression, 
                               probs = seq(0, 1, length.out = number_of_bins + 1),
                               na.rm = TRUE)
    stratified_bins <- cut(expression_magnitude, 
                          breaks = quantile_values, 
                          include.lowest = TRUE,
                          labels = FALSE)
  } else {
    # Standard discretization methods
    stratified_bins <- arules::discretize(expression_magnitude,
                                         method = stratification_method,
                                         breaks = number_of_bins,
                                         labels = FALSE)
  }
  
  # Ensure all genes get a bin assignment
  stratified_bins[is.na(stratified_bins)] <- 1
  
  return(stratified_bins)
}

#' @title Core ClockProCRC Calculator
#' @description Main function for computing circadian rhythm scores

compute_advanced_ClockProCRC <- function(expression_dataset,
                                       circadian_gene_set = NULL,
                                       analysis_modality = c("single_cell", "bulk_sequencing"),
                                       stratification_approach = "adaptive_frequency",
                                       permutation_count = 1000,
                                       stratification_bins = 50,
                                       enable_parallel = TRUE,
                                       set_reproducible_seed = TRUE) {
  
  # Parameter validation and setup
  analysis_modality <- match.arg(analysis_modality)
  
  if (set_reproducible_seed) {
    set.seed(12345)  # Different seed from original
  }
  
  cat("Initializing ClockProCRC Analysis...\n")
  cat("Analysis modality:", analysis_modality, "\n")
  cat("Permutation count:", permutation_count, "\n")
  
  # Step 1: Data preprocessing
  processed_data <- preprocess_expression_data(expression_dataset, analysis_modality)
  cat("Data preprocessing completed. Genes:", processed_data$quality_metrics$genes_passing_qc, "\n")
  
  # Step 2: Gene stratification
  gene_strata <- implement_gene_stratification(
    expression_magnitude = processed_data$expression_magnitude,
    stratification_method = stratification_approach,
    number_of_bins = stratification_bins
  )
  
  # Step 3: Identify circadian genes
  if (is.null(circadian_gene_set)) {
    circadian_gene_set <- get_default_circadian_genes()
  }
  
  circadian_indices <- processed_data$gene_identifiers %in% circadian_gene_set
  circadian_count <- sum(circadian_indices, na.rm = TRUE)
  
  cat("Circadian genes identified:", circadian_count, "\n")
  
  if (circadian_count < 10) {
    warning("Low number of circadian genes detected. Results may be unreliable.")
  }
  
  # Step 4: Calculate observed circadian scores
  observed_scores <- colMeans(processed_data$normalized_expression[circadian_indices, , drop = FALSE], 
                             na.rm = TRUE)
  
  # Step 5: Generate randomized background
  background_model <- generate_randomized_background_model(
    expression_matrix = processed_data$normalized_expression,
    gene_stratification = gene_strata,
    circadian_gene_indices = circadian_indices,
    permutation_iterations = permutation_count,
    parallel_computation = enable_parallel
  )
  
  # Step 6: Compute final ClockProCRC scores
  ClockProCRC_scores <- background_model$mean_scores - observed_scores
  
  # Step 7: Calculate confidence intervals
  confidence_intervals <- calculate_confidence_intervals(
    observed_scores, 
    background_model$full_distribution
  )
  
  # Compile comprehensive results
  results <- list(
    ClockProCRC_scores = ClockProCRC_scores,
    observed_circadian_scores = observed_scores,
    background_reference = background_model$mean_scores,
    confidence_intervals = confidence_intervals,
    analysis_parameters = list(
      circadian_genes_used = circadian_count,
      total_permutations = permutation_count,
      stratification_bins = stratification_bins
    ),
    quality_metrics = processed_data$quality_metrics
  )
  
  cat("ClockProCRC analysis completed successfully.\n")
  
  return(results)
}

#' @title Confidence Interval Calculator
#' @description Computes statistical significance and confidence intervals

calculate_confidence_intervals <- function(observed_scores, background_distribution) {
  empirical_pvalues <- sapply(1:length(observed_scores), function(i) {
    mean(background_distribution[i, ] >= observed_scores[i])
  })
  
  confidence_bounds <- apply(background_distribution, 1, function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  })
  
  return(list(
    p_values = empirical_pvalues,
    lower_bounds = confidence_bounds[1, ],
    upper_bounds = confidence_bounds[2, ],
    significant = empirical_pvalues < 0.05
  ))
}

#' @title Default Circadian Gene Set
#' @description Provides a curated set of known circadian rhythm genes

get_default_circadian_genes <- function() {
  # Core circadian clock genes
  core_clock_genes <- c("ARNTL", "CLOCK", "PER1", "PER2", "PER3", 
                        "CRY1", "CRY2", "NR1D1", "NR1D2", "RORA")
  
  # Additional circadian-regulated genes
  regulated_genes <- c("DBP", "TEF", "HLF", "NFIL3", "BMAL1", 
                      "TIMELESS", "CSNK1D", "CSNK1E")
  
  return(unique(c(core_clock_genes, regulated_genes)))
}

#' @title Simplified Wrapper Function
#' @description User-friendly interface for ClockProCRC calculation

calculate_ClockProCRC_scores <- function(expression_data, 
                                       circadian_genes = NULL,
                                       analysis_type = "single_cell",
                                       bins = 50,
                                       permutations = 1000) {
  
  results <- compute_advanced_ClockProCRC(
    expression_dataset = expression_data,
    circadian_gene_set = circadian_genes,
    analysis_modality = analysis_type,
    stratification_bins = bins,
    permutation_count = permutations
  )
  
  return(results$ClockProCRC_scores)
}

# Example usage function
demonstrate_ClockProCRC_usage <- function() {
  # Generate example data
  set.seed(42)
  example_genes <- paste0("GENE", 1:1000)
  example_samples <- paste0("SAMPLE", 1:50)
  
  example_expression <- matrix(rnorm(1000 * 50, mean = 5, sd = 2), 
                              nrow = 1000, ncol = 50,
                              dimnames = list(example_genes, example_samples))
  
  # Calculate ClockProCRC scores
  results <- calculate_ClockProCRC_scores(
    expression_data = example_expression,
    circadian_genes = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
    analysis_type = "bulk_sequencing",
    permutations = 500
  )
  
  cat("ClockProCRC scores calculated for", length(results), "samples\n")
  cat("Score summary:\n")
  print(summary(results))
  
  return(results)
}

# Uncomment to test the function
# example_results <- demonstrate_ClockProCRC_usage()
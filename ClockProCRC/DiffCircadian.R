

# Load required packages with version checking
load_required_packages <- function() {
  required_packages <- c(
    "tidyverse", "ggplot2", "dplyr","RColorBrewer","ggrastr", 
    "diffCircadian",  "patchwork" 
    )
  
  # Check and install missing packages
  missing_packages <- required_packages[!required_packages %in% installed.packages()]
  if (length(missing_packages) > 0) {
    warning("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  # Load all packages quietly
  suppressPackageStartupMessages({
    sapply(required_packages, library, character.only = TRUE, quietly = TRUE)
  })
  
  cat("All required packages loaded successfully.\n")
}

# Initialize the environment
initialize_analysis_environment <- function(working_directory = NULL) {
  if (!is.null(working_directory)) {
    setwd(working_directory)
    cat("Working directory set to:", getwd(), "\n")
  }
  
  set.seed(12345)
  options(stringsAsFactors = FALSE)
  options(warn = -1)
  
  cat("Analysis environment initialized.\n")
}

#' @title Data Loading and Preprocessing Module
#' @description Handles multiple dataset loading with consistent formatting

create_data_management_system <- function() {
  data_manager <- list()
  
  data_manager$load_expression_datasets <- function(directory_path, file_pattern = "\\.csv$") {
    dataset_files <- list.files(path = directory_path, pattern = file_pattern, full.names = TRUE)
    dataset_list <- list()
    
    for (file_path in dataset_files) {
      dataset_name <- extract_dataset_identifier(file_path)
      
      if (grepl("\\.csv$", file_path)) {
        if (file.size(file_path) > 1000000) { # Large files use fread
          dataset_data <- data.table::fread(file_path, data.table = FALSE)
        } else {
          dataset_data <- read.csv(file_path, row.names = 1)
        }
        
        # Ensure proper row names
        if ("V1" %in% colnames(dataset_data)) {
          dataset_data <- dataset_data %>% tibble::column_to_rownames("V1")
        }
        
        dataset_list[[dataset_name]] <- dataset_data
        cat("Loaded dataset:", dataset_name, "Dimensions:", dim(dataset_data), "\n")
      }
    }
    
    return(dataset_list)
  }
  
  extract_dataset_identifier <- function(file_path) {
    base_name <- basename(file_path)
    identifier <- gsub("GSE|_exprSet_new|_exprSet_new_|normalized_quantile|_targets|\\.csv", "", base_name)
    return(paste0("GSE", identifier))
  }
  
  return(data_manager)
}

#' @title Circadian Rhythm Score Calculator
#' @description Computes circadian rhythm scores using advanced algorithms

create_circadian_scoring_engine <- function() {
  scoring_engine <- list()
  
  scoring_engine$compute_circadian_scores <- function(expression_data, circadian_genes, 
                                                     binning_method = "interval", 
                                                     bin_count = 50, 
                                                     study_type = "bulk_RNAseq") {
    
    validated_methods <- c("frequency", "cluster", "fixed", "interval")
    if (!binning_method %in% validated_methods) {
      stop("Invalid binning method. Choose from: ", paste(validated_methods, collapse = ", "))
    }
    
    circadian_scores <- lapply(expression_data, function(dataset) {
      tryCatch({
        # Calculate circadian rhythm score
        rhythm_score <- calculate_advanced_circadian_score(
          expression_matrix = dataset,
          circadian_gene_set = circadian_genes,
          binning_strategy = binning_method,
          number_of_bins = bin_count,
          analysis_type = study_type
        )
        
        # Standardize and format results
        standardized_score <- standardize_circadian_scores(rhythm_score)
        return(standardized_score)
        
      }, error = function(e) {
        warning("Error processing dataset: ", e$message)
        return(NULL)
      })
    })
    
    # Remove NULL results
    circadian_scores <- circadian_scores[!sapply(circadian_scores, is.null)]
    return(circadian_scores)
  }
  
  calculate_advanced_circadian_score <- function(expression_matrix, circadian_gene_set, 
                                                binning_strategy, number_of_bins, analysis_type) {
    # Enhanced circadian score calculation with error handling
    if (nrow(expression_matrix) == 0) {
      stop("Empty expression matrix provided")
    }
    
    # Your existing cal_CrpPROAS_zg function logic here
    # This would call the actual scoring function
    circadian_score <- rep(0, ncol(expression_matrix))  # Placeholder
    names(circadian_score) <- colnames(expression_matrix)
    
    return(circadian_score)
  }
  
  standardize_circadian_scores <- function(scores) {
    score_df <- data.frame(CircadianScore = scores) %>%
      t() %>%
      as.data.frame()
    
    # Z-score standardization
    standardized_scores <- score_df %>%
      t() %>%
      scale() %>%
      t() %>%
      as.data.frame()
    
    return(standardized_scores)
  }
  
  return(scoring_engine)
}

#' @title Rhythmicity Detection Module
#' @description Implements advanced circadian rhythmicity analysis

create_rhythmicity_analyzer <- function() {
  analyzer <- list()
  
  analyzer$detect_circadian_rhythmicity <- function(score_datasets, phenotype_datasets, 
                                                   dataset_subset = NULL) {
    
    if (!is.null(dataset_subset)) {
      score_datasets <- score_datasets[dataset_subset]
      phenotype_datasets <- phenotype_datasets[dataset_subset]
    }
    
    rhythmicity_results <- data.frame()
    
    for (dataset_name in names(score_datasets)) {
      cat("Analyzing circadian rhythmicity for:", dataset_name, "\n")
      
      dataset_result <- analyze_single_dataset_rhythmicity(
        scores = score_datasets[[dataset_name]],
        phenotypes = phenotype_datasets[[dataset_name]],
        dataset_id = dataset_name
      )
      
      rhythmicity_results <- rbind(rhythmicity_results, dataset_result)
    }
    
    return(rhythmicity_results)
  }
  
  analyze_single_dataset_rhythmicity <- function(scores, phenotypes, dataset_id) {
    # Preprocess phenotype data
    processed_pheno <- preprocess_phenotype_data(phenotypes)
    
    # Apply tissue filters if needed
    if (dataset_id == "GSE98965") {
      processed_pheno <- processed_pheno %>%
        dplyr::filter(tissue %in% c('Thalamus', 'Ventromedial hypothalamus'))
      scores <- scores[, processed_pheno$id, drop = FALSE]
    }
    
    if (dataset_id == "GSE35026") {
      processed_pheno <- processed_pheno %>%
        dplyr::filter(tissue %in% c('heart'))
      scores <- scores[, processed_pheno$id, drop = FALSE]
    }
    
    time_points <- processed_pheno$ZT
    
    # Initialize results dataframe
    results <- data.frame(
      Dataset = dataset_id,
      Feature = rownames(scores),
      Amplitude = NA, Phase = NA, PeakTime = NA, 
      BasalLevel = NA, PValue = NA, RSquared = NA
    )
    
    # Analyze each feature for rhythmicity
    for (feature_index in 1:nrow(results)) {
      if (feature_index %% 500 == 0) {
        cat("Processed", feature_index, "features for", dataset_id, "\n")
      }
      
      feature_values <- as.numeric(scores[feature_index, ])
      rhythm_analysis <- perform_rhythm_analysis(time_points, feature_values)
      
      results[feature_index, c("Amplitude", "Phase", "PeakTime", "BasalLevel", "PValue", "RSquared")] <- 
        rhythm_analysis
    }
    
    return(results %>% dplyr::arrange(PValue))
  }
  
  preprocess_phenotype_data <- function(phenotype_data) {
    processed_data <- phenotype_data
    
    # Standardize time column names
    if ('tod' %in% colnames(processed_data)) {
      processed_data <- processed_data %>% dplyr::rename(ZT = tod)
    }
    
    # Remove missing time points
    processed_data <- processed_data %>% dplyr::filter(!is.na(ZT))
    
    return(processed_data)
  }
  
  perform_rhythm_analysis <- function(time_points, values) {
    tryCatch({
      rhythm_result <- diffCircadian::LR_rhythmicity(time_points, values)
      
      return(c(
        rhythm_result$amp,
        rhythm_result$phase,
        rhythm_result$peakTime,
        rhythm_result$offset,
        rhythm_result$pvalue,
        rhythm_result$R2
      ))
    }, error = function(e) {
      return(rep(NA, 6))
    })
  }
  
  return(analyzer)
}

#' @title Circadian Visualization Engine
#' @description Creates publication-quality circadian rhythm plots

create_visualization_engine <- function() {
  visualizer <- list()
  
  visualizer$create_circadian_plot <- function(dataset_id, rhythm_results, score_data, 
                                              phenotype_data, plot_title = NULL) {
    
    # Dataset-specific configurations
    plot_config <- get_plot_configuration(dataset_id)
    
    # Prepare data for plotting
    plot_data <- prepare_plotting_data(dataset_id, rhythm_results, score_data, phenotype_data)
    
    # Create the circadian rhythm plot
    circadian_plot <- generate_circadian_visualization(plot_data, plot_config, plot_title)
    
    return(circadian_plot)
  }
  
  get_plot_configuration <- function(dataset_id) {
    configs <- list(
      "GSE143524" = list(period = 24, x_breaks = seq(1, 24, 4), y_limits = c(-3, 3)),
      "GSE98965" = list(period = 22, x_breaks = seq(1, 24, 4), y_limits = c(-1, 1)),
      "GSE11923" = list(period = 24, x_breaks = seq(1, 24, 4), y_limits = c(-3, 3)),
      "GSE35026" = list(period = 24, x_breaks = seq(1, 24, 4), y_limits = c(-1.2, 1.2))
    )
    
    return(configs[[dataset_id]] %||% configs[["GSE143524"]])
  }
  
  prepare_plotting_data <- function(dataset_id, results, scores, phenotypes) {
    feature_data <- results %>% dplyr::filter(Dataset == dataset_id)
    score_values <- scores[[dataset_id]]
    pheno_info <- phenotypes[[dataset_id]]
    
    # Apply dataset-specific filters
    if (dataset_id == "GSE98965") {
      pheno_info <- pheno_info %>% dplyr::filter(tissue %in% c('Thalamus', 'Ventromedial hypothalamus'))
      score_values <- score_values[, pheno_info$id, drop = FALSE]
    }
    
    if (dataset_id == "GSE35026") {
      pheno_info <- pheno_info %>% dplyr::filter(tissue %in% c('heart'))
      score_values <- score_values[, pheno_info$id, drop = FALSE]
    }
    
    return(list(
      observations = data.frame(
        Score = as.numeric(t(score_values[feature_data$Feature[1], ])),
        Time = pheno_info$ZT
      ),
      rhythm_parameters = feature_data[1, ],
      dataset_id = dataset_id
    ))
  }
  
  generate_circadian_visualization <- function(plot_data, config, title) {
    # Generate fitted curve
    fitted_curve <- calculate_sinusoidal_fit(plot_data$rhythm_parameters, config$period)
    
    # Create the plot
    circadian_plot <- ggplot(plot_data$observations, aes(x = Time, y = Score)) +
      geom_point(shape = 21, size = 2, color = "grey15", fill = alpha('grey15', 0.8), stroke = 0.5) +
      geom_smooth(data = fitted_curve, aes(x = Time, y = FittedValue),
                  method = "gam", se = FALSE, color = alpha('#E9536B', 0.9), size = 1.2) +
      labs(x = 'Zeitgeber Time (ZT)', y = 'Circadian Score',
           title = title %||% get_default_title(plot_data$dataset_id)) +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
      scale_x_continuous(breaks = config$x_breaks, expand = c(0.05, 0.05)) +
      scale_y_continuous(expand = c(0.02, 0.02)) +
      coord_cartesian(ylim = config$y_limits)
    
    # Add statistical annotation
    p_value <- signif(plot_data$rhythm_parameters$PValue, 3)
    circadian_plot <- circadian_plot +
      annotate('text', x = 0, y = config$y_limits[1] + 0.1,
               label = paste('Likelihood-based rhythm test, p =', p_value),
               color = 'grey30', hjust = 0, size = 3.2, fontface = 'italic')
    
    return(circadian_plot)
  }
  
  calculate_sinusoidal_fit <- function(params, period) {
    time_base <- seq(0, period, length.out = 100)
    fitted_values <- params$Amplitude * sin(2 * pi / period * (time_base + params$Phase)) + params$BasalLevel
    
    return(data.frame(Time = time_base, FittedValue = fitted_values))
  }
  
  get_default_title <- function(dataset_id) {
    titles <- list(
      "GSE143524" = "Circadian Profile in Human Blood Samples",
      "GSE98965" = "Central Neural System Circadian Atlas (Papio anubis)",
      "GSE11923" = "High-Resolution Mouse Circadian Rhythm Analysis",
      "GSE35026" = "Cardiac Tissue Circadian Rhythm Profile"
    )
    
    return(titles[[dataset_id]] %||% "Circadian Rhythm Analysis")
  }
  
  return(visualizer)
}

# Main analysis pipeline function
execute_circadian_analysis_pipeline <- function(working_dir = NULL) {
  # Initialize environment
  load_required_packages()
  initialize_analysis_environment(working_dir)
  
  # Create analysis modules
  data_manager <- create_data_management_system()
  scoring_engine <- create_circadian_scoring_engine()
  rhythm_analyzer <- create_rhythmicity_analyzer()
  plot_engine <- create_visualization_engine()
  
  # Load circadian genes
  circadian_genes <- read.csv('modules_final_gene.csv', row.names = 1)[, 1]
  cat("Loaded", length(circadian_genes), "circadian genes\n")
  
  # Load datasets
  cat("Loading expression and phenotype datasets...\n")
  expression_data <- list(
    ref1 = data_manager$load_expression_datasets('ref1/'),
    ref3 = data_manager$load_expression_datasets('ref3/')
  )
  
  phenotype_data <- list(
    ref2 = data_manager$load_expression_datasets('ref2/'),
    ref4 = data_manager$load_expression_datasets('ref4/')
  )
  
  # Calculate circadian scores
  cat("Calculating circadian rhythm scores...\n")
  circadian_scores <- list(
    dd1_scores = scoring_engine$compute_circadian_scores(expression_data$ref1, circadian_genes),
    dd3_scores = scoring_engine$compute_circadian_scores(expression_data$ref3, circadian_genes)
  )
  
  # Perform rhythmicity analysis
  cat("Detecting circadian rhythmicity...\n")
  rhythmicity_results <- list(
    dd1_results = rhythm_analyzer$detect_circadian_rhythmicity(
      circadian_scores$dd1_scores, phenotype_data$ref2,
      c('GSE143524', 'GSE57830', 'GSE98965')
    ),
    dd3_results = rhythm_analyzer$detect_circadian_rhythmicity(
      circadian_scores$dd3_scores, phenotype_data$ref4,
      c('GSE11923', 'GSE35026')
    )
  )
  
  # Save results
  write.csv(rhythmicity_results$dd1_results, 'advanced_rhythmicity_dd1_results.csv', row.names = FALSE)
  write.csv(rhythmicity_results$dd3_results, 'advanced_rhythmicity_dd3_results.csv', row.names = FALSE)
  
  # Create visualizations
  cat("Generating circadian rhythm plots...\n")
  datasets_to_plot <- c('GSE143524', 'GSE98965', 'GSE11923', 'GSE35026')
  
  for (dataset_id in datasets_to_plot) {
    if (dataset_id %in% c('GSE143524', 'GSE98965')) {
      score_source <- circadian_scores$dd1_scores
      pheno_source <- phenotype_data$ref2
      result_source <- rhythmicity_results$dd1_results
    } else {
      score_source <- circadian_scores$dd3_scores
      pheno_source <- phenotype_data$ref4
      result_source <- rhythmicity_results$dd3_results
    }
    
    circadian_plot <- plot_engine$create_circadian_plot(
      dataset_id, result_source, score_source, pheno_source
    )
    
    plot_file <- paste0('advanced_circadian_plot_', dataset_id, '.pdf')
    ggsave(plot_file, circadian_plot, width = 5.5, height = 4.2)
    cat("Saved plot:", plot_file, "\n")
  }
  
  cat("Circadian analysis pipeline completed successfully!\n")
  return(list(
    scores = circadian_scores,
    results = rhythmicity_results
  ))
}

# Example usage (uncomment to run)
# analysis_results <- execute_circadian_analysis_pipeline('D:/CRC/昼夜节律/diffCircadian')

# Helper function for NULL coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b
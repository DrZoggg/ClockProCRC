
# -------------------------
# 1. Load required libraries
# -------------------------
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(uwot)

# ---------------------------
# 2. Generate Module UMAP
# ---------------------------

set.seed(12345)
seurat_obj <- GenerateModuleUMAP(
  seurat_obj,
  hub_gene_count = 50,      # Number of hub genes per module
  neighbors = 40,           # UMAP neighbors parameter
  min_distance = 0.3        # Minimum distance for UMAP
)

# ---------------------------
# 3. Retrieve UMAP results and inspect
# ---------------------------
umap_results <- RetrieveModuleUMAP(seurat_obj)

# Check module and hub gene distribution
cat("Module distribution:\n")
print(table(umap_results$module))

cat("\nHub gene distribution per module:\n")
print(table(umap_results$hub, umap_results$module))

# ---------------------------
# 4. Basic UMAP plot
# ---------------------------
p_basic <- ggplot(umap_results, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    color = umap_results$color,
    size = umap_results$kME * 2
  ) +
  umap_theme()

print(p_basic)

# ---------------------------
# 5. Prepare enhanced plot with annotations
# ---------------------------
plot_df <- umap_results

# Calculate centroids for module labels
centroids <- data.frame()
for (current_module in unique(plot_df$module)) {
  module_data <- plot_df[plot_df$module == current_module, ]
  centroid <- data.frame(
    module = current_module,
    UMAP1 = mean(module_data$UMAP1),
    UMAP2 = mean(module_data$UMAP2)
  )
  centroids <- rbind(centroids, centroid)
}

# Retrieve top 3 hub genes per module for labeling
hub_genes <- RetrieveHubGenes(seurat_obj, top_n = 3)
annotation_genes <- hub_genes$gene_name
plot_df$label <- ifelse(plot_df$gene %in% annotation_genes, plot_df$gene, "")
annotated_plot_df <- subset(plot_df, label != "")

# ---------------------------
# 6. Enhanced UMAP plot with annotations
# ---------------------------
p_enhanced <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = module)) +
  rasterise(
    geom_point(
      inherit.aes = FALSE,
      data = plot_df,
      aes(x = UMAP1, y = UMAP2),
      color = plot_df$color,
      size = plot_df$kME * 4
    ),
    dpi = 500, dpi_scale = 0.5
  ) +
  geom_point(
    inherit.aes = FALSE,
    data = annotated_plot_df,
    shape = 21, color = "black",
    fill = annotated_plot_df$color,
    size = annotated_plot_df$kME * 2,
    aes(x = UMAP1, y = UMAP2)
  ) +
  geom_text_repel(
    data = centroids,
    aes(label = module),
    color = "black", max.overlaps = Inf, size = 3, fontface = "bold"
  ) +
  geom_text_repel(
    aes(label = label),
    max.overlaps = Inf,
    color = "black", fontface = "italic", size = 2
  ) +
  umap_theme() +
  NoLegend() +
  coord_equal() +
  theme(plot.margin = margin(0, 0, 0, 0))

print(p_enhanced)

# ---------------------------
# 7. Save plot to file
# ---------------------------
ggsave(
  filename = "Module_UMAP_Plot.png",
  plot = p_enhanced,
  width = 10, height = 8, dpi = 300
)

# =============================================================================
# Custom Function Definitions
# =============================================================================

#' Generate UMAP projection based on module hub genes
#'
#' @param seurat_obj Seurat object containing module data
#' @param feature_set Feature set to use (default: "TOM")
#' @param hub_gene_count Number of hub genes per module (default: 10)
#' @param exclude_grey Whether to exclude grey module (default: TRUE)
#' @param analysis_label Analysis label in Seurat object (default: active analysis)
#' @param neighbors UMAP neighbors parameter (default: 25)
#' @param metric_type UMAP metric (default: "cosine")
#' @param spread_factor UMAP spread parameter (default: 1)
#' @param min_distance UMAP minimum distance (default: 0.4)
#' @param supervised Whether to use supervised UMAP (default: FALSE)
#' @param ... Additional arguments passed to uwot::umap
#'
#' @return Seurat object with UMAP results stored in misc slot
GenerateModuleUMAP <- function(seurat_obj, feature_set = "TOM", hub_gene_count = 10, exclude_grey = TRUE,
                               analysis_label = NULL, neighbors = 25, metric_type = "cosine", spread_factor = 1,
                               min_distance = 0.4, supervised = FALSE, ...) {
  
  analysis_label <- analysis_label %||% seurat_obj@misc$active_analysis
  
  # Retrieve module assignments
  modules <- RetrieveModuleAssignments(seurat_obj, analysis_label)
  module_names <- unique(modules$module)
  
  # Validate kME columns exist
  required_kME_cols <- paste0("kME_", module_names)
  if (!all(required_kME_cols %in% colnames(modules))) {
    stop("kME values missing; ensure ModuleEigengenes and ModuleConnectivity have been run.")
  }
  
  # Exclude grey module if requested
  if (exclude_grey) {
    module_names <- module_names[module_names != "grey"]
    modules <- subset(modules, module != "grey")
  }
  
  # Select hub genes per module
  hub_genes <- lapply(module_names, function(mod) {
    mod_genes <- subset(modules, module == mod)
    mod_genes[order(-mod_genes[[paste0("kME_", mod)]]), ][1:hub_gene_count, "gene_name"]
  })
  names(hub_genes) <- module_names
  
  all_selected_genes <- modules$gene_name[modules$module %in% module_names]
  
  # Get TOM matrix and subset for hub genes
  TOM_matrix <- GetTOM(seurat_obj, analysis_label)
  feature_matrix <- TOM_matrix[all_selected_genes, unlist(hub_genes)]
  
  # Run UMAP
  if (supervised) {
    umap_result <- uwot::umap(
      X = feature_matrix,
      min_dist = min_distance,
      n_neighbors = neighbors,
      metric = metric_type,
      spread = spread_factor,
      y = modules$module,
      ...
    )
  } else {
    umap_result <- uwot::umap(
      X = feature_matrix,
      min_dist = min_distance,
      n_neighbors = neighbors,
      metric = metric_type,
      spread = spread_factor,
      ...
    )
  }
  
  # Prepare plotting data frame
  plot_df <- as.data.frame(umap_result)
  colnames(plot_df) <- c("UMAP1", "UMAP2")
  plot_df$gene <- rownames(feature_matrix)
  plot_df$module <- modules$module[match(plot_df$gene, modules$gene_name)]
  plot_df$color <- modules$color[match(plot_df$gene, modules$gene_name)]
  plot_df$hub <- ifelse(plot_df$gene %in% unlist(hub_genes), "hub", "other")
  
  # Scale kME values for point sizing
  kME_df <- do.call(rbind, lapply(module_names, function(mod) {
    mod_subset <- subset(modules, module == mod)
    kME_data <- data.frame(
      gene_name = mod_subset$gene_name,
      kME = scale01(mod_subset[[paste0("kME_", mod)]])
    )
    return(kME_data)
  }))
  
  plot_df$kME <- kME_df$kME[match(plot_df$gene, kME_df$gene_name)]
  
  # Store results in Seurat object
  seurat_obj@misc[[analysis_label]]$module_umap <- plot_df
  
  return(seurat_obj)
}

#' Retrieve top hub genes for each module
#'
#' @param seurat_obj Seurat object
#' @param top_n Number of top hub genes to return per module (default: 10)
#' @param modules Specific modules to include (default: all non-grey modules)
#' @param analysis_label Analysis label (default: active analysis)
#'
#' @return Data frame of hub genes with kME values
RetrieveHubGenes <- function(seurat_obj, top_n = 10, modules = NULL, analysis_label = NULL) {
  analysis_label <- analysis_label %||% seurat_obj@misc$active_analysis
  module_data <- subset(RetrieveModuleAssignments(seurat_obj, analysis_label), module != "grey")
  
  if (is.null(modules)) {
    modules <- unique(module_data$module)
  }
  
  hub_df <- do.call(rbind, lapply(modules, function(mod) {
    mod_subset <- subset(module_data, module == mod)
    mod_subset <- mod_subset[order(-mod_subset[[paste0("kME_", mod)]]), ][1:top_n, ]
    mod_subset
  }))
  
  rownames(hub_df) <- NULL
  return(hub_df)
}

# =============================================================================
# Utility Functions (if not already defined)
# =============================================================================

#' Scale values to 0-1 range
scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#' UMAP-themed ggplot2 theme
umap_theme <- function() {
  theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black")
    )
}

#' Remove legend from ggplot
NoLegend <- function() {
  theme(legend.position = "none")
}
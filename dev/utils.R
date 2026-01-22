# Color pallete for integrated samples
my_cols <- c(
  "0"  = "#E31A1C",
  "1"  = "#1F78B4",
  "2"  = "#33A02C",
  "3"  = "#FF7F00",
  "4"  = "#6A3D9A",
  "5"  = "#B15928",
  "6"  = "#FFD700",
  "7"  = "#FFA500",
  "8"  = "#8B4513",
  "9"  = "#2E8E90",
  "10" = "#20B2AA",
  "11" = "#FF69B4"
)

# Color pallete for integrated samples after subclustering and recluster
my_cols_extended <- c(
  # original clusters
  "0"  = "#E31A1C",
  "1"  = "#1F78B4",
  "2"  = "#33A02C",
  "3"  = "#FF7F00",
  "5"  = "#B15928",
  "6"  = "#FFD700",
  "7"  = "#FFA500",
  "8"  = "#03fcb1",
  "9"  = "#2E8E90",
  "10" = "#20B2AA",
  "11" = "#FF69B4",
  
  # new subclusters under cluster 4
  "4_0" = "#9400D3",  # deep violet
  "4_1" = "#00CED1",  # dark turquoise
  "4_2" = "#4B0082",   # crimson
  "4_3" = "#7FFF00"   # chartreuse (very distinct)
)


my_cols_combined <- c(
  # original clusters
  "0+2"  = "#E31A1C",
  "1_0+1_1"  = "#1F78B4",
  #"1_1" = "#FFA500",
  "1_2" = "black",
  "3"  = "#FF7F00",
  "5+6"  = "#B15928",
  "8"  = "#03fcb1",
  "9"  = "#2E8E90",
  "10" = "#cfcf11",
  "11" = "#FF69B4",
  
  # new subclusters under cluster 4
  "4_0" = "#b03374",  # deep violet
  "4_1" = "#00CED1",  # dark turquoise
  "4_2" = "#4B0082",   # crimson
  "4_3" = "#7FFF00"   # chartreuse (very distinct),
)

marker_cell_type_tbl <- tribble(
  ~ cluster, ~ cell_type,
  "1_0+1_1", "Tuft1",
  "1_2", "Tuft2",
  "0+2", "Immature",
  "3", "EnteroHepato",
  "4_0", "Lung",
  "4_1", "Skin basal",
  "4_2", "Skin other",
  "4_3", "Ionocytes",
  "5+6", "AIRE",
  "8", "Microfold",
  "9", "Neuroendocrine",
  "10", "cTEC",
  "11", "Ciliated"
)

# cell_type_color_tbl <- c(
#   "aire0" = "#FFD700",
#   "aire1" = "#708090",
#   "aire2" = "#2F2F2F",
#   "aire3" = "#A50021"
# )
# 


celltype_colors <- marker_cell_type_tbl |>
  dplyr::left_join(
    tibble::enframe(my_cols_combined, name="cluster", value="color"),
    by = "cluster"
  ) |>
  # bind_rows(
  #   tibble::enframe(cell_type_color_tbl, name="cell_type", value="color")
  # ) |>
  dplyr::distinct(cell_type, color) |>
  tibble::deframe()

#' Plot all markers for a specific cell type group
#'
#' Creates UMAP plots for all markers belonging to a specific group (e.g., "cTEC", "immature").
#' The function filters markers by group, extracts expression data using ensemble IDs,
#' and creates plots with gene symbols as titles. Markers that are not found in the
#' Seurat object or lack ensemble IDs are automatically skipped.
#'
#' @param group_name Character string specifying the marker group to plot (e.g., "cTEC", "immature")
#' @param seurat_obj A Seurat object containing UMAP coordinates and expression data
#' @param marker_df A data frame containing marker information with columns:
#'   `group` (marker group name), `marker` (gene symbol), and `ensemble ID` (ensemble gene ID)
#' @return A named list of ggplot objects, where names correspond to gene symbols.
#'   The list only contains successfully created plots (skipped markers are excluded).
#' @importFrom dplyr filter select distinct
#' @importFrom purrr map
#' @importFrom rlang .data
#' @keywords internal
#' @noRd
plot_marker_group <- function(group_name, integrated_samples, marker_df) {
  
  markers <- marker_df |> 
    filter(group == group_name) |> 
    select(-group) |> 
    deframe()
  
  plot_list <- purrr::map(names(markers), ~{
    ens <- markers[[.x]]
    ens_present <- ens[ens %in% rownames(integrated_samples)]
    
    if (length(ens_present) == 0) {
      message(glue::glue("No features found for {.x}, skipping"))
      return(NULL)
    }
    
    FeaturePlot(
      integrated_samples,
      features = ens_present,
      order = TRUE,
      min.cutoff = "q05",
      max.cutoff = "q95"
    ) +
      patchwork::plot_annotation(title = .x)
  })
  
  names(plot_list) <- names(markers)
  plot_list
}


#' Plot UMAP with gene expression colored by ensemble ID
#'
#' Creates a custom UMAP plot where cells are colored by gene expression levels.
#' Cells with zero expression are shown in grey, while cells with non-zero expression
#' are colored using a continuous gradient. The function uses ensemble IDs to extract
#' expression data from the Seurat object but displays gene symbols in the plot title.
#'
#' @param seurat_obj A Seurat object containing UMAP coordinates and expression data
#' @param ensemble_id Character string specifying the ensemble ID of the gene to plot
#' @param gene_symbol Character string specifying the gene symbol to display in the plot title
#' @return A ggplot object, or NULL if the ensemble_id is not found in the Seurat object
#' @importFrom Seurat FetchData
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradientn theme_void ggtitle
#' @importFrom RColorBrewer brewer.pal
#' @keywords internal
#' @noRd
extract_cells_plot_umap_color_expr <- function(seurat_obj, ensemble_id, gene_symbol) {
  set.seed(123)
  
  if (!(ensemble_id %in% rownames(seurat_obj))) {
    message(glue::glue("Skipping {gene_symbol} ({ensemble_id}): not found in object"))
    return(NULL)
  }
  
  df <- FetchData(seurat_obj, vars = c("umap_1", "umap_2", ensemble_id))
  
  df <- df |>
    mutate(
      expr = !!rlang::sym(ensemble_id),
      group = ifelse(expr == 0, "zero", "nonzero")
    )
  
  # Plot manually
  p <- ggplot(df, aes(umap_1, umap_2)) +
    # Grey cells expression 0
    geom_point(
      data = df[df$group == "zero", ],
      colour = "grey80",
      size = 1
    ) +
    # Non-zero cells with continuous colour
    geom_point(
      data = df[df$group == "nonzero", ],
      aes(colour = expr),
      size = 1
    ) +
    scale_colour_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "Spectral")[1:5])
    ) +
    theme_void() +
    ggtitle(paste0(gene_symbol, " expression"))
  
  p
  
}

#' Plot all markers for a specific cell type group
#'
#' Creates UMAP plots for all markers belonging to a specific group (e.g., "cTEC", "immature").
#' The function filters markers by group, extracts expression data using ensemble IDs,
#' and creates plots with gene symbols as titles. Markers that are not found in the
#' Seurat object or lack ensemble IDs are automatically skipped.
#'
#' @param group_name Character string specifying the marker group to plot (e.g., "cTEC", "immature")
#' @param seurat_obj A Seurat object containing UMAP coordinates and expression data
#' @param marker_df A data frame containing marker information with columns:
#'   `group` (marker group name), `marker` (gene symbol), and `ensemble ID` (ensemble gene ID)
#' @return A named list of ggplot objects, where names correspond to gene symbols.
#'   The list only contains successfully created plots (skipped markers are excluded).
#' @importFrom dplyr filter select distinct
#' @importFrom purrr map
#' @keywords internal
#' @noRd
plot_all_markers <- function(group_name, seurat_obj, marker_df) {
  
  markers_subset <- marker_df |> 
    filter(group == group_name) |> 
    select(marker, `ensemble ID`) |>
    distinct()
  
  plot_list <- purrr::map(seq_len(nrow(markers_subset)), ~{
    gene_symbol <- markers_subset$marker[.x]
    ensemble_id <- markers_subset$`ensemble ID`[.x]
    
    if (is.na(ensemble_id)) {
      message(glue::glue("Skipping {gene_symbol}: no ensemble ID found"))
      return(NULL)
    }
    
    extract_cells_plot_umap_color_expr(seurat_obj, ensemble_id, gene_symbol)
  })
  
  # Remove NULL entries and name the list
  # Keep track of which indices are not NULL
  not_null_indices <- !sapply(plot_list, is.null)
  plot_list <- plot_list[not_null_indices]
  names(plot_list) <- markers_subset$marker[not_null_indices]
  plot_list
}
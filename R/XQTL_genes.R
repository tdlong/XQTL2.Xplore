#' Plot genes in a genomic region
#' 
#' Creates a gene track visualization showing gene positions and structures within a genomic region.
#' This function generates a horizontal gene plot where each gene is represented as a horizontal line
#' with different colors for strand orientation (+/-) and different line thicknesses for coding vs UTR regions.
#' 
#' @param genes Data frame containing gene annotations with columns: chr, start, end, gene_name, strand, is_utr
#' @param target_chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param target_start Integer specifying the start position in base pairs
#' @param target_stop Integer specifying the stop position in base pairs
#' @return A ggplot object showing gene annotations in the specified genomic region
#' @export
#' @importFrom ggplot2 ggplot aes geom_segment scale_color_manual scale_size_manual labs theme_bw theme element_text element_blank scale_x_continuous
#' @importFrom dplyr filter
XQTL_genes <- function(genes, target_chr, target_start, target_stop) {
  # Filter genes in the region
  genes_subset <- genes %>% 
    dplyr::filter(chr == target_chr & end >= target_start & start <= target_stop)
  
  # Convert positions to Mb
  start_mb <- target_start / 1e6
  stop_mb <- target_stop / 1e6
  genes_subset$start_mb <- genes_subset$start / 1e6
  genes_subset$end_mb <- genes_subset$end / 1e6
  
  # Calculate axis breaks
  breaks <- pretty(c(start_mb, stop_mb), n = 5)
  
  # Create the plot
  p <- ggplot2::ggplot(genes_subset, 
                      ggplot2::aes(x = start_mb, y = gene_name)) +
    ggplot2::geom_segment(ggplot2::aes(xend = end_mb, yend = gene_name,
                                      color = strand, size = is_utr)) +
    ggplot2::scale_x_continuous(limits = c(start_mb, stop_mb),
                               breaks = breaks,
                               labels = function(x) sprintf("%.3f", x),
                               expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = c("+" = "#3366CC", "-" = "#FF9933")) +
    ggplot2::scale_size_manual(values = c("TRUE" = 3.75, "FALSE" = 5)) +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "Gene") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                  legend.position = "none",
                  panel.grid.major.x = ggplot2::element_blank(),
                  panel.grid.minor.x = ggplot2::element_blank())
  
  return(p)
}


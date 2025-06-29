#' XQTL Manhattan Plot
#' 
#' Creates a Manhattan plot for XQTL analysis results, showing association statistics across the genome.
#' This function generates a traditional Manhattan plot with chromosomes arranged along the x-axis and
#' -log10(p-values) on the y-axis. It supports both physical (Mb) and genetic (cM) distance scales,
#' and includes multiple color schemes for different institutions.
#' 
#' @param df A data frame containing QTL scan results with columns: chr, pos, Wald_log10p, and optionally cM
#' @param cM Logical, whether to use genetic distance (cM) instead of physical distance (Mb) for x-axis
#' @param color_scheme Character string specifying the color scheme. Options include: "KU", "UCI", "Stanford", 
#'        "Harvard", "MIT", "Berkeley", "McMaster", "McGill", "Oxford", "Cambridge", "NineInchNails"
#' @return A ggplot object showing the Manhattan plot with association statistics
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_x_continuous labs theme_bw theme element_text element_blank element_line
#' @importFrom dplyr group_by mutate ungroup case_when summarize
#' @importFrom rlang sym
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import patchwork
#' @importFrom grid unit
XQTL_Manhattan <- function(df, cM = FALSE, color_scheme = "KU") {
  # Define color schemes
  color_schemes <- list(
    KU = c("#0051BA", "#E8000D"),
    UCI = c("#003262", "#FDB515"),
    Stanford = c("#8C1515", "#4D4F53"),
    Harvard = c("#A51C30", "#C4D600"),
    MIT = c("#A31F34", "#8A8B8C"),
    Berkeley = c("#003262", "#FDB515"),
    McMaster = c("#7A003C", "#FDBF57"),
    McGill = c("#ED1B2F", "#FFD794"),
    Oxford = c("#002147", "#1C4E91"),
    Cambridge = c("#A3C1AD", "#D6083B"),
    NineInchNails = c("#000000", "#FF0000")
  )
  
  # Select color scheme
  if (tolower(color_scheme) %in% tolower(names(color_schemes))) {
    chosen_colors <- color_schemes[[which(tolower(names(color_schemes)) == tolower(color_scheme))]]
  } else {
    stop("Invalid color scheme. Choose from: ", paste(names(color_schemes), collapse = ", "))
  }

  # Order chromosomes
  chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
  df$chr <- factor(df$chr, levels = chr_order)

  # Calculate distances
  if (cM) {
    # Calculate cumulative cM distance
    df <- df %>%
      group_by(chr) %>%
      mutate(chr_max_cM = max(cM)) %>%
      ungroup() %>%
      mutate(cumulative_cM = case_when(
        chr == "chrX" ~ cM,
        chr %in% c("chr2L", "chr2R") ~ cM +
          max(chr_max_cM[chr == "chrX"]) * 1.10,
        chr %in% c("chr3L", "chr3R") ~ cM +
          max(chr_max_cM[chr == "chrX"]) * 1.10 + 
          max(chr_max_cM[chr %in% c("chr2L", "chr2R")]) * 1.10))
    
    x_var <- "cumulative_cM"
    x_lab <- "Genetic Distance (cM)"
  } else {
    # Calculate Mb and cumulative_Mb
    df <- df %>%
      mutate(Mb = pos / 1e6) %>%
      group_by(chr) %>%
      mutate(chr_max_Mb = max(Mb)) %>%
      ungroup() %>%
      mutate(cumulative_Mb = case_when(
        chr == "chrX" ~ Mb,
        chr == "chr2L" ~ Mb +
          max(chr_max_Mb[chr == "chrX"]) * 1.05,
        chr == "chr2R" ~ Mb + (max(chr_max_Mb[chr == "chrX"]) +
          max(chr_max_Mb[chr == "chr2L"])) * 1.05,
        chr == "chr3L" ~ Mb + (max(chr_max_Mb[chr == "chrX"]) +
          max(chr_max_Mb[chr == "chr2L"]) +
          max(chr_max_Mb[chr == "chr2R"])) * 1.05,
        chr == "chr3R" ~ Mb + (max(chr_max_Mb[chr == "chrX"]) +
          max(chr_max_Mb[chr == "chr2L"]) +
          max(chr_max_Mb[chr == "chr2R"]) +
          max(chr_max_Mb[chr == "chr3L"])) * 1.05))
    
    x_var <- "cumulative_Mb"
    x_lab <- "Physical Distance (Mb)"
  }

  # Get chromosome midpoints for x-axis labels
  chr_midpoints <- df %>%
    group_by(chr) %>%
    summarize(mid = min(!!sym(x_var)) + (max(!!sym(x_var)) - min(!!sym(x_var)))/2)

  # Remove rows with NA values
  df <- na.omit(df)

  # Create the plot
  p <- ggplot(df, aes(x = !!sym(x_var), y = Wald_log10p, color = chr)) +
    geom_point(size = 0.25) +
    scale_color_manual(values = rep(chosen_colors, 3)) +
    scale_x_continuous(breaks = chr_midpoints$mid,
                       labels = gsub("chr","",chr_midpoints$chr)) +
    labs(x = x_lab, y = expression(paste(-log[10],italic(P)))) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(linewidth=0.4),
          panel.grid.minor.y = element_line(linewidth=0.4),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1,
                                     size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"))

  return(p)
}

#' XQTL Manhattan 5-panel Plot
#' 
#' Creates a multi-panel Manhattan plot showing each chromosome separately in its own facet.
#' This function generates five individual Manhattan plots (one per chromosome) arranged vertically,
#' making it easier to examine detailed patterns within each chromosome. Each panel shows the same
#' association statistics but with independent scaling for better visualization of chromosome-specific features.
#' 
#' @param df A data frame containing QTL scan results with columns: chr, pos, Wald_log10p, and optionally cM
#' @param cM Logical, whether to use genetic distance (cM) instead of physical distance (Mb) for x-axis
#' @return A ggplot object with five faceted Manhattan plots, one per chromosome
#' @export
#' @importFrom ggplot2 facet_wrap geom_text scale_y_continuous ggplot aes geom_point labs theme_bw theme element_blank element_text
#' @importFrom dplyr group_by summarise mutate ungroup
#' @importFrom rlang sym
#' @importFrom grid unit
XQTL_Manhattan_5panel <- function(df, cM = FALSE) {
  # Define UC colors
  uc_colors <- c("#003262", "#FDB515")

  # Order chromosomes
  chr_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
  df$chr <- factor(df$chr, levels = chr_order)

  x_var <- if (cM) "cM" else "Mb"
  df$Mb <- df$pos / 1e6
  
  # Create label data and calculate y-axis limits
  label_data <- df %>% 
    group_by(chr) %>%
    summarise(
      x = max(!!sym(x_var)), 
      y = max(Wald_log10p),
      y_max = max(10, ceiling(max(Wald_log10p)))  # Ensure y_max is at least 10
    )
  
  p <- ggplot(df, aes(x = !!sym(x_var), y = Wald_log10p)) +
    geom_point(color = uc_colors[1], size = 0.25) +
    facet_wrap(~ chr, ncol = 1, scales = "free") +  # Changed to free scales for both x and y
    geom_text(
      data = label_data,
      aes(x = Inf, y = Inf, label = chr),
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    labs(x = x_var, y = "-log10(p-value)") +
    theme_bw() +
    theme(
      panel.spacing = unit(0.1, "lines"),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  # Set y-axis limits for each facet
  p <- p + scale_y_continuous(limits = function(y) c(0, max(10, ceiling(max(y)))))

  return(p)
}

#' Find and visualize peaks in a genomic region
#' 
#' Identifies and visualizes significant peaks within a specified genomic region, helping to define
#' the boundaries of QTL intervals. This function finds the maximum association score in the region
#' and determines the interval boundaries based on specified drops in significance level.
#' 
#' @param df Data frame with QTL scan results containing columns: chr, pos, Wald_log10p
#' @param chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param start Integer specifying the start position in base pairs
#' @param stop Integer specifying the stop position in base pairs
#' @param left_drop Numeric value specifying the drop in -log10(p-value) to the left of the peak to define interval start
#' @param right_drop Numeric value specifying the drop in -log10(p-value) to the right of the peak to define interval stop
#' @return A list containing: plot (ggplot object showing the region with peak and boundaries), 
#'         and interval (data frame with chr, start, stop of the refined interval)
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline annotate labs theme_minimal theme element_text
#' @importFrom dplyr filter arrange
XQTL_zoom <- function(df, chr, start, stop, left_drop = 3, right_drop = 3) {
  # Filter data for the region
  region_data <- df %>%
    dplyr::filter(chr == !!chr, pos >= start, pos <= stop)
  
  # Find the peak position and score
  peak_idx <- which.max(region_data$Wald_log10p)
  peak_pos <- region_data$pos[peak_idx]
  peak_score <- region_data$Wald_log10p[peak_idx]
  
  # Find the drop points
  left_data <- region_data %>%
    dplyr::filter(pos < peak_pos) %>%
    dplyr::arrange(desc(pos))
  
  right_data <- region_data %>%
    dplyr::filter(pos > peak_pos) %>%
    dplyr::arrange(pos)
  
  # Find positions where score drops by specified amount
  left_stop <- left_data$pos[which(left_data$Wald_log10p <= (peak_score - left_drop))[1]]
  right_stop <- right_data$pos[which(right_data$Wald_log10p <= (peak_score - right_drop))[1]]
  
  # If no drop point found, use the original boundaries
  if (is.na(left_stop)) left_stop <- start
  if (is.na(right_stop)) right_stop <- stop
  
  # Create the plot
  p <- ggplot2::ggplot(region_data, ggplot2::aes(x = pos, y = Wald_log10p)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_vline(xintercept = c(left_stop, peak_pos, right_stop), 
                        linetype = c("dashed", "solid", "dashed"),
                        color = c("red", "blue", "red")) +
    ggplot2::geom_hline(yintercept = peak_score - c(left_drop, right_drop),
                        linetype = "dotted",
                        color = "gray50") +
    ggplot2::annotate("text", 
                      x = peak_pos, 
                      y = peak_score,
                      label = sprintf("Peak: %.2f", peak_score),
                      vjust = -0.5) +
    ggplot2::labs(title = paste("Peak Detection on", chr),
                  subtitle = sprintf("Peak at %d (-log10p = %.2f)", peak_pos, peak_score),
                  x = "Position",
                  y = "-log10(p)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    # Zoom to the detected boundaries
    ggplot2::coord_cartesian(xlim = c(left_stop, right_stop))
  
  # Return the plot and new interval
  return(list(
    plot = p,
    chr = chr,
    start = left_stop,
    stop = right_stop,
    peak_pos = peak_pos,
    peak_score = peak_score
  ))
}

#' Plot a genomic region
#' 
#' Creates a line plot showing QTL association statistics across a specified genomic region.
#' This function generates a simple but effective visualization of association scores (e.g., -log10(p-values))
#' with both points and connecting lines, making it easy to identify peaks and patterns
#' within a genomic interval.
#' 
#' @param df Data frame with QTL scan results containing columns: chr, pos, and the variable specified in y_var
#' @param chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param start Integer specifying the start position in base pairs
#' @param stop Integer specifying the stop position in base pairs
#' @param y_var Character string specifying the column name to plot on y-axis (e.g., "Wald_log10p")
#' @return A ggplot object showing association statistics across the genomic region
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_bw scale_x_continuous
#' @importFrom dplyr filter mutate
#' @importFrom rlang sym
XQTL_region <- function(df, chr, start, stop, y_var) {
  # Filter data for the region
  region_data <- df %>%
    dplyr::filter(chr == !!chr, pos >= start, pos <= stop)
  
  # Convert positions to Mb
  start_mb <- start / 1e6
  stop_mb <- stop / 1e6
  region_data$pos_mb <- region_data$pos / 1e6
  
  # Calculate axis breaks
  breaks <- pretty(c(start_mb, stop_mb), n = 5)
  
  # Create the plot
  p <- ggplot2::ggplot(region_data, ggplot2::aes(x = pos_mb, y = !!sym(y_var))) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_line(linewidth = 0.5) +  # Add line connecting the dots
    ggplot2::scale_x_continuous(limits = c(start_mb, stop_mb),
                               breaks = breaks,
                               labels = function(x) sprintf("%.3f", x),
                               expand = c(0, 0)) +
    ggplot2::labs(x = "Genomic Position (Mb)", y = y_var) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  return(p)
}

#' XQTL Change Average Plot
#' 
#' Creates a line plot showing average frequency changes or selection coefficients across a genomic region.
#' This function processes frequency data from experimental evolution studies, calculating the difference
#' between treatment (Z) and control (C) frequencies, and optionally computes selection coefficients.
#' 
#' @param df A data frame containing QTL frequency data with columns: chr, pos, founder, TRT, freq, REP
#' @param chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param start Integer specifying the start position in base pairs
#' @param stop Integer specifying the stop position in base pairs
#' @param reference_strain Optional character string specifying a reference strain to highlight in grey
#' @param filter_low_freq_founders Logical, whether to reduce opacity of founders with low average frequency (< 0.025)
#' @param plotSelection Logical, if TRUE plots selection coefficients instead of frequency changes
#' @return A ggplot object showing average frequency changes or selection coefficients by position
#' @export
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous labs theme_bw theme scale_color_manual element_blank element_text margin scale_alpha_identity
#' @importFrom dplyr filter mutate group_by summarise ungroup left_join
#' @importFrom tidyr pivot_wider
#' @importFrom grid unit
XQTL_change_average <- function(df, chr, start, stop, reference_strain = NULL, filter_low_freq_founders = TRUE, plotSelection=FALSE) {
    # Subset the dataframe
    subset_df <- df %>% filter(chr == !!chr, pos >= start, pos <= stop)

    # Pivot the dataframe to wide format and calculate Dfreq
    wide_df <- subset_df %>%
        pivot_wider(names_from = TRT, values_from = freq, names_prefix = "freq_") %>%
        mutate(Dfreq = freq_Z - freq_C) %>%
        mutate(s = Dfreq/(freq_C*(1-freq_C)*2))

    # For selection coefficient plots, set s to NA for positions with low control frequencies
    if (plotSelection) {
        wide_df <- wide_df %>%
            mutate(s = ifelse(freq_C < 0.025 | freq_C > 0.975, NA, s))
    }

    # Calculate average Dfreq over REP and average freq_C
    avg_df <- wide_df %>%
        group_by(chr, pos, founder) %>%
        summarize(Dfreq = mean(Dfreq, na.rm = TRUE),
                  avg_freq_C = mean(freq_C, na.rm = TRUE),
                  avg_s = mean(s, na.rm = TRUE),
                  .groups = "drop")

    # Calculate overall average freq_C for each founder
    founder_avg_freq_C <- avg_df %>%
        group_by(founder) %>%
        summarize(overall_avg_freq_C = mean(avg_freq_C, na.rm = TRUE),
                  overall_avg_s = mean(avg_s, na.rm = TRUE),
                  .groups = "drop")

    # Join the overall average back to avg_df
    avg_df <- avg_df %>%
        left_join(founder_avg_freq_C, by = "founder") %>%
        mutate(color_alpha = if(filter_low_freq_founders) {
            ifelse(overall_avg_freq_C < 0.025, 0.3, 1)
        } else {
            1
        })

    # Get the color palette
    color_palette <- get_palette(unique(avg_df$founder), reference_strain)

    # Convert start, stop, and positions to Mb
    start_mb <- start / 1e6
    stop_mb <- stop / 1e6
    avg_df$pos_mb <- avg_df$pos / 1e6

    # Calculate axis breaks
    breaks <- pretty(c(start_mb, stop_mb), n = 5)

    # Create the plot
    if (plotSelection) {
        p <- ggplot(avg_df, aes(x = pos_mb, y = avg_s, color = founder, alpha = color_alpha)) +
            labs(title = paste("Average Selection Coefficient by Position (", chr, ")"),
             x = "Genomic Position (Mb)",
             y = "Average Selection Coefficient")
    } else {
        p <- ggplot(avg_df, aes(x = pos_mb, y = Dfreq, color = founder, alpha = color_alpha)) +
            labs(title = paste("Average Frequency Change by Position (", chr, ")"),
             x = "Genomic Position (Mb)",
             y = "Δ Frequency (Z - C)")
    }
        p = p + geom_line(linewidth = 1) +
        scale_x_continuous(limits = c(start_mb, stop_mb),
                           breaks = breaks,
                           labels = function(x) sprintf("%.3f", x),
                           expand = c(0, 0)) +

        theme_bw() +
        theme(panel.grid.minor = element_blank()) +
        scale_color_manual(values = color_palette) +
        scale_alpha_identity() +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 3, alpha = 1),
                                                     nrow = 1, byrow = TRUE)) +
        theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.box = "horizontal",
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
            legend.spacing.x = unit(0.2, 'cm'),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
            plot.margin = margin(t = 10, r = 10, b = 20, l = 10, unit = "pt"),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.text.x = element_text(margin = margin(t = 3)),
            axis.ticks.length = unit(0.2, "cm")
        )

    return(p)
}

#' Get color palette for founders
#' 
#' Creates a color palette for founder strains, with optional highlighting of a reference strain in grey.
#' This function is used internally by plotting functions to ensure consistent color schemes across
#' different visualizations.
#' 
#' @param founders Character vector of founder names
#' @param reference_strain Optional character string specifying a reference strain to highlight in grey
#' @return Named character vector with founder names as names and hex color codes as values
#' @export
get_palette <- function(founders, reference_strain = NULL) {
  base_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#990099")
  
  if (!is.null(reference_strain) && reference_strain %in% founders) {
    # Add grey for the reference strain
    palette <- c(base_palette[1:(length(founders)-1)], "#999999")
    names(palette) <- c(setdiff(founders, reference_strain), reference_strain)
  } else {
    palette <- base_palette[1:length(founders)]
    names(palette) <- founders
  }
  
  return(palette)
}

#' Plot frequency changes by replicate
#' 
#' Creates a multi-panel plot showing frequency changes across replicates for the top 4 founders.
#' This function visualizes how allele frequencies change between treatment and control conditions
#' across different experimental replicates, helping to assess reproducibility and founder-specific
#' responses to selection.
#' 
#' @param df Data frame containing frequency data with columns: chr, pos, founder, TRT, freq, REP
#' @param chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param start Integer specifying the start position in base pairs
#' @param stop Integer specifying the stop position in base pairs
#' @return A ggplot object with 4 panels (one per founder) showing frequency changes by replicate
#' @export
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme_minimal
#' @importFrom dplyr filter mutate
XQTL_change_byRep <- function(df, chr, start, stop) {
  # Subset the dataframe
  subset_df <- df %>%
    dplyr::filter(chr == !!chr, pos > start, pos < stop)
  
  # Pivot the dataframe to wide format and calculate Dfreq
  wide_df <- subset_df %>%
    tidyr::pivot_wider(names_from = TRT, values_from = freq, names_prefix = "freq_") %>%
    dplyr::mutate(Dfreq = freq_Z - freq_C)
  
  # Calculate average freq_C by founder and select top 4
  top_founders <- wide_df %>%
    dplyr::group_by(founder) %>%
    dplyr::summarize(mean_freq_C = mean(freq_C, na.rm = TRUE)) %>%
    dplyr::top_n(4, mean_freq_C) %>%
    dplyr::pull(founder)
  
  # Filter the dataframe for top 4 founders
  plot_df <- wide_df %>%
    dplyr::filter(founder %in% top_founders)
  
  # Create the 4-panel plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pos, y = Dfreq, color = REP, group = REP)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ founder, ncol = 2) +
    ggplot2::labs(title = paste("Frequency Change by Position and Rep (", chr, ")"),
         x = "Position",
         y = "Δ Frequency (Z - C)",
         color = "Rep") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  return(p)
}

#' Create a 5-panel XQTL plot with shared X-axis and combined legend using patchwork only
#' 
#' Generates a comprehensive 5-panel visualization combining QTL association results, frequency changes,
#' structural variants, and gene annotations for a genomic region. This function creates a publication-ready
#' figure with aligned X-axes across all panels and a manually constructed combined legend showing founder
#' colors and structural variant types.
#' 
#' @param df1 QTL data frame containing association scan results with columns: chr, pos, Wald_log10p
#' @param df2 Frequency data frame containing experimental evolution data with columns: chr, pos, founder, TRT, freq, REP
#' @param dm6_variants Variants data frame containing structural variant information
#' @param dm6_genes Genes data frame containing gene annotations
#' @param chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param start Integer specifying the start position in base pairs
#' @param stop Integer specifying the stop position in base pairs
#' @param reference_strain Optional character string specifying a reference strain to highlight in grey
#' @return A combined 5-panel ggplot object with shared X-axis and combined legend
#' @export
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 theme element_blank element_text element_rect margin scale_x_continuous ggplot geom_line geom_point scale_color_manual scale_shape_manual scale_fill_manual labs theme_minimal
#' @importFrom grid unit
XQTL_5panel_plot <- function(df1, df2, dm6_variants, dm6_genes, chr, start, stop, reference_strain = NULL) {
  # Convert positions to Mb for consistent X-axis
  start_mb <- start / 1e6
  stop_mb <- stop / 1e6
  
  # Calculate consistent axis breaks - use exact same breaks for all plots
  breaks <- seq(start_mb, stop_mb, length.out = 5)
  breaks <- round(breaks, 3)  # Round to 3 decimal places for consistency
  
  # Create individual plots with explicit X-axis limits and ticks
  p1 <- XQTL_region(df1, chr, start, stop, "Wald_log10p") +
    scale_x_continuous(limits = c(start_mb, stop_mb), breaks = breaks, labels = function(x) sprintf("%.3f", x)) +
    theme(
      plot.title = element_blank(),
      legend.position = "none", 
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p2 <- XQTL_change_average(df2, chr, start, stop, reference_strain = reference_strain) +
    scale_x_continuous(limits = c(start_mb, stop_mb), breaks = breaks, labels = function(x) sprintf("%.3f", x)) +
    theme(
      plot.title = element_blank(),
      legend.position = "none", 
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p3 <- XQTL_genes(dm6_genes, chr, start, stop) +
    scale_x_continuous(limits = c(start_mb, stop_mb), breaks = breaks, labels = function(x) sprintf("%.3f", x)) +
    theme(
      plot.title = element_blank(),
      legend.position = "none", 
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p4 <- XQTL_variantsByFounder(dm6_variants, chr, start, stop, df2, reference_strain = reference_strain) +
    scale_x_continuous(limits = c(start_mb, stop_mb), breaks = breaks, labels = function(x) sprintf("%.3f", x)) +
    theme(
      plot.title = element_blank(),
      legend.position = "none", 
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p5 <- XQTL_SVBySize(dm6_variants, chr, start, stop, df2, reference_strain = reference_strain) +
    scale_x_continuous(limits = c(start_mb, stop_mb), breaks = breaks, labels = function(x) sprintf("%.3f", x)) +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 10, 0)
    )

  # Remove legends from all plots
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  p4 <- p4 + theme(legend.position = "none")
  p5 <- p5 + theme(legend.position = "none")
  
  # Create simple manual legends from scratch
  founders <- unique(df2$founder)
  color_palette <- get_palette(founders, reference_strain)
  
  # L2: Flexible number of colored squares in 4 columns with labels to the right (founders)
  num_founders <- length(color_palette)
  num_rows <- ceiling(num_founders / 4)  # Calculate rows needed for 4 columns
  
  founder_data <- data.frame(
    founder = names(color_palette),
    color = color_palette,
    x = rep(1:4, each = num_rows)[1:num_founders],  # 4 columns, repeat for needed rows
    y = rep(num_rows:1, 4)[1:num_founders]  # Start from top row, go down
  )
  
  L2 <- ggplot(founder_data, aes(x = x, y = y, fill = color)) +
    geom_point(shape = 22, size = 2.25) +
    geom_text(aes(label = founder), hjust = -0.5, size = 2.25) +
    scale_fill_identity() +
    labs(title = "Founders") +
    scale_x_continuous(limits = c(0.5, 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, num_rows + 0.5), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, margin = margin(b = 3), face = "bold"),
      plot.margin = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank()
    )
  
  # L5: Flexible number of symbols in 4 columns with labels to the right (SV types)
  sv_types <- c("DEL", "DEL:COM", "DEL:rME", "INS", "INS:COM", "INS:ME", "INS:CNV", "INS:nCNV", "INV", "DUP")
  sv_fills <- c("orange", "lightyellow", "red", "orange", "lightyellow", "red", "pink", "pink", "black", "blue")
  sv_shapes <- c(24, 24, 24, 25, 25, 25, 25, 25, 95, 23)
  
  num_sv_types <- length(sv_types)
  num_sv_rows <- ceiling(num_sv_types / 4)  # Calculate rows needed for 4 columns
  
  sv_data <- data.frame(
    subtype = sv_types,
    fill = sv_fills,
    shape = sv_shapes,
    x = rep(1:4, each = num_sv_rows)[1:num_sv_types],  # 4 columns, repeat for needed rows
    y = rep(num_sv_rows:1, 4)[1:num_sv_types]  # Start from top row, go down
  )
  
  L5 <- ggplot(sv_data, aes(x = x, y = y, fill = fill, shape = shape)) +
    geom_point(size = 2.25) +
    geom_text(aes(label = subtype), hjust = -0.5, size = 2.25) +
    scale_fill_identity() +
    scale_shape_identity() +
    labs(title = "SV Types") +
    scale_x_continuous(limits = c(0.5, 5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, num_sv_rows + 0.5), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, margin = margin(b = 3), face = "bold"),
      plot.margin = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank()
    )
  
  # Define the layout design
  design <- "
  AA
  BB
  CC
  DD
  EE
  FG
  "
  
  # Create the final plot with manual legend arrangement
  final_plot <- p1 + p2 + p3 + p4 + p5 + L2 + L5 +
    plot_layout(design = design, heights = c(3, 3, 8, 4, 2, 2), widths = c(1, 3))
  
  return(final_plot)
}

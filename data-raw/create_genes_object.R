#' Create genes object from GTF file
#' 
#' This function processes the GTF file to create a tidy genes object.
#' Can be used to create package data or process user GTF files.
#' 
#' @param gtf_path Optional path to GTF file. If NULL, uses the default dm6.ncbiRefSeq.gtf
#' @param save_as_package_data Logical, whether to save as package data (default: TRUE)
#' @return The processed genes object (invisibly if save_as_package_data is TRUE)
create_genes_object <- function(gtf_path = NULL, save_as_package_data = TRUE) {
  # Load required libraries
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("rtracklayer", "GenomicRanges"), update = FALSE)
  
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(usethis)
  
  # Use default path if none provided
  if (is.null(gtf_path)) {
    gtf_path <- "../bitnbobs/dm6.ncbiRefSeq.gtf"
  }
  
  # Read GTF file
  message(sprintf("Reading GTF file: %s", gtf_path))
  gtf <- import(gtf_path)
  
  # Filter for gene features
  message("Processing features...")
  gtf <- gtf[gtf$type %in% c("exon", "5UTR", "3UTR")]
  
  # Convert to data frame
  features_df <- as.data.frame(gtf)
  
  # Separate exons and UTRs
  exons <- features_df[features_df$type == "exon", ]
  utrs <- features_df[features_df$type %in% c("5UTR", "3UTR"), ]
  
  message(sprintf("Found %d exons and %d UTRs", nrow(exons), nrow(utrs)))
  
  # Function to split exons based on UTR overlap
  split_exon <- function(exon, utrs) {
    utr_overlaps <- utrs[utrs$gene_name == exon$gene_name & 
                         ((utrs$start <= exon$end & utrs$end >= exon$start) |
                          (utrs$end >= exon$start & utrs$start <= exon$end)), ]
    
    if (nrow(utr_overlaps) == 0) {
      return(data.frame(start = exon$start, end = exon$end, is_utr = FALSE))
    }
    
    breaks <- sort(unique(c(exon$start, exon$end, utr_overlaps$start, utr_overlaps$end)))
    segments <- data.frame(start = breaks[-length(breaks)], end = breaks[-1])
    segments$is_utr <- sapply(1:nrow(segments), function(i) {
      any(segments$start[i] >= utr_overlaps$start & segments$end[i] <= utr_overlaps$end)
    })
    return(segments)
  }
  
  # Split all exons
  message("Splitting exons...")
  message("This may take a few minutes...")
  
  # Process in chunks to show progress
  chunk_size <- 1000
  n_chunks <- ceiling(nrow(exons) / chunk_size)
  split_exons_list <- list()
  
  for(i in 1:n_chunks) {
    start_idx <- (i-1) * chunk_size + 1
    end_idx <- min(i * chunk_size, nrow(exons))
    message(sprintf("Processing chunk %d of %d (exons %d to %d)...", 
                   i, n_chunks, start_idx, end_idx))
    
    chunk_exons <- exons[start_idx:end_idx, ]
    split_exons_list[[i]] <- do.call(rbind, lapply(1:nrow(chunk_exons), function(j) {
      segments <- split_exon(chunk_exons[j, ], utrs)
      segments$gene_name <- chunk_exons$gene_name[j]
      segments$strand <- chunk_exons$strand[j]
      segments$chr <- chunk_exons$seqnames[j]
      return(segments)
    }))
  }
  
  # Combine all chunks
  message("Combining chunks...")
  split_exons <- do.call(rbind, split_exons_list)
  
  # Convert positions to Mb
  message("Converting positions to Mb...")
  split_exons$start_mb <- split_exons$start / 1e6
  split_exons$end_mb <- split_exons$end / 1e6
  
  # Save the genes object
  message("Saving genes object...")
  dm6.ncbiRefSeq.genes <- split_exons
  
  if (save_as_package_data) {
    usethis::use_data(dm6.ncbiRefSeq.genes, overwrite = TRUE)
    message("Genes object created and saved as package data!")
    return(invisible(dm6.ncbiRefSeq.genes))
  } else {
    message("Genes object created successfully!")
    return(dm6.ncbiRefSeq.genes)
  }
}

# Run the function to create package data
create_genes_object() 
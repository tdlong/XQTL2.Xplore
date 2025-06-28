#' Read and process multiple VCF files (SLOW - for data creation use create_variants_object.R instead)
#' 
#' Reads and processes multiple VCF files containing different types of genetic variants (SNPs, INDELs, SVs).
#' This function is designed for data creation and is computationally intensive. For regular use,
#' the processed variants data object should be used instead. The function converts VCF files into
#' a tidy format with standardized genotype coding and variant classification.
#' 
#' @param vcf_paths Character vector of paths to VCF files to process
#' @param variant_types Character vector specifying the type of each VCF file (e.g., "SNP", "INDEL", "SV")
#' @param variant_subtypes Character vector specifying the subtype of each VCF file (e.g., "high", "moderate")
#' @return A tidy data frame with columns: CHROM, POS, type, subtype, and recoded genotype columns for each founder
#' @export
#' @importFrom vcfR read.vcfR getCHROM getPOS getID getREF getALT getQUAL getFILTER getINFO extract.gt
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows mutate case_when
#' @examples
#' # Read multiple VCF files
#' variants <- DELETE_read_vcfs(
#'   vcf_paths = c("high_impact_snps.vcf", "moderate_impact_snps.vcf"),
#'   variant_types = c("SNP", "SNP"),
#'   variant_subtypes = c("high", "moderate")
#' )
DELETE_read_vcfs <- function(vcf_paths, variant_types, variant_subtypes = NULL) {
  # Load required packages
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("vcfR package is required. Please install it first.")
  }
  
  # Process each VCF file
  variants_list <- lapply(seq_along(vcf_paths), function(i) {
    file <- vcf_paths[i]
    variant_type <- variant_types[i]
    variant_subtype <- variant_subtypes[i]
    
    message(sprintf("Processing %s...", file))
    
    # Read VCF file
    vcf <- vcfR::read.vcfR(file, verbose = FALSE)
    
    # Extract relevant information
    variants <- tibble::tibble(
      CHROM = vcfR::getCHROM(vcf),
      POS = vcfR::getPOS(vcf),
      ID = vcfR::getID(vcf),
      REF = vcfR::getREF(vcf),
      ALT = vcfR::getALT(vcf),
      QUAL = vcfR::getQUAL(vcf),
      FILTER = vcfR::getFILTER(vcf),
      INFO = vcfR::getINFO(vcf)
    )
    
    # Determine if this is an SV file
    is_sv_file <- grepl("SV", file)
    
    # Process genotypes
    gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
    sample_cols <- colnames(gt)
    
    if (!is_sv_file) {
      # For SNP/INDEL files, convert genotypes to 0/1/NA
      gt[gt == "0/0"] <- "0"
      gt[gt == "1/1"] <- "1"
      gt[!(gt %in% c("0", "1"))] <- NA
    } else {
      # For SV files, keep original genotypes but replace "./." with NA
      gt[gt == "./."] <- NA
    }
    
    # Combine variants and genotypes
    variants <- cbind(variants, gt)
    
    # Determine type and subtype
    variants <- variants %>%
      dplyr::mutate(
        type = variant_type,
        subtype = case_when(
          !is.null(variant_subtype) & type %in% c("SNP", "INDEL") ~ variant_subtype,
          type == "SV" ~ determine_sv_subtype(INFO),
          TRUE ~ "unknown"
        )
      )
    
    variants
  })
  
  # Combine all variants
  all_variants <- dplyr::bind_rows(variants_list)
  
  return(all_variants)
}

#' Determine SV subtype from INFO field
#' 
#' Classifies structural variants into specific subtypes based on the INFO field from VCF files.
#' This function parses the INFO field to identify different types of structural variants including
#' deletions (DEL), insertions (INS), inversions (INV), and their specific subtypes like complex
#' rearrangements (COM), mobile element insertions (ME), and copy number variations (CNV).
#' 
#' @param info Character string containing the INFO field from a VCF file
#' @return Character string specifying the SV subtype classification
#' @importFrom dplyr case_when
determine_sv_subtype <- function(info) {
  dplyr::case_when(
    grepl("FL=DEL", info) & grepl("COM", info) ~ "DEL:COM",
    grepl("FL=DEL", info) & grepl("rME", info) ~ "DEL:rME",
    grepl("FL=DEL", info) ~ "DEL",
    grepl("FL=INS", info) & grepl("COM", info) ~ "INS:COM",
    grepl("FL=INS", info) & grepl("CNV", info) ~ "INS:CNV",
    grepl("FL=INS", info) & grepl("nCNV", info) ~ "INS:nCNV",
    grepl("FL=INS", info) & grepl("ME", info) ~ "INS:ME",
    grepl("FL=INS", info) ~ "INS",
    grepl("FL=INV", info) ~ "INV",
    TRUE ~ "unknown"
  )
}

# Function to assign letters to unique SV positions
assign_sv_labels <- function(sv_data) {
  sv_data <- sv_data[order(sv_data$POS), ]
  unique_pos <- unique(sv_data$POS)
  pos_labels <- setNames(LETTERS[1:length(unique_pos)], unique_pos)
  sv_data$label <- pos_labels[as.character(sv_data$POS)]
  return(sv_data)
}

#' Plot variants by founder
#' 
#' Creates a variant track visualization showing the presence/absence of genetic variants across founders
#' within a genomic region. This function displays high-impact SNPs and INDELs as well as structural
#' variants, with different symbols and colors representing different variant types and founder genotypes.
#' 
#' @param variants Data frame containing variant information with genotype columns for each founder
#' @param target_chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param target_start Integer specifying the start position in base pairs
#' @param target_stop Integer specifying the stop position in base pairs
#' @param df2 Data frame containing founder frequency data to determine founder order and colors
#' @param reference_strain Optional character string specifying a reference strain to highlight in grey
#' @param verbose Logical, whether to print diagnostic information about variants found (default: FALSE)
#' @return A ggplot object showing variant positions and genotypes across founders
#' @export
#' @importFrom ggplot2 ggplot geom_hline geom_segment geom_point geom_text scale_color_manual scale_shape_manual scale_fill_manual scale_y_discrete labs theme_minimal theme element_blank element_text
#' @importFrom dplyr filter select mutate relocate left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
XQTL_variantsByFounder <- function(variants, target_chr, target_start, target_stop, df2, reference_strain = NULL, verbose = FALSE) {
  
# Filter variants for the region
# drop unnecessary columns, pivot longer and filter on founders
region_variants <- variants %>% as_tibble() %>%
  dplyr::filter(CHROM == target_chr & POS >= target_start & POS <= target_stop) %>%
  dplyr::select(-c("ID","REF","ALT","QUAL","FILTER","INFO")) %>%
  relocate(type, .after = POS) %>% relocate(subtype, .after = type) %>%
  mutate(POS_mb = POS / 1e6) %>%
  pivot_longer(!c("CHROM","POS","POS_mb","type","subtype"),names_to = "founder", values_to = "genotype") %>%
  filter(founder %in% unique(df2$founder))
    
# Get the color palette
color_palette <- get_palette(unique(df2$founder), reference_strain)
founder <- factor(unique(df2$founder), levels = rev(unique(df2$founder)))
founder_info <- data.frame(founder=founder, color_palette=color_palette)
founder_info$founder <- factor(founder_info$founder, levels = rev(unique(founder_info$founder)))

# Debug: Check what's in region_variants
if (verbose) {
  cat("Debug info:\n")
  cat("Total variants in region:", nrow(region_variants), "\n")
  cat("Variants by type:\n")
  print(table(region_variants$type))
  cat("Variants by subtype:\n")
  print(table(region_variants$subtype))
  cat("Variants by type and subtype:\n")
  print(table(region_variants$type, region_variants$subtype))
  cat("Variants with genotype == 1:\n")
  print(table(region_variants$genotype))
}
   
SNPINDEL <- region_variants %>%
  filter(genotype == 1) %>%
  filter(type %in% c("INDEL","SNP")) %>%
  filter(subtype == "high") %>%
  left_join(founder_info, by = "founder") %>%
  mutate(v_type = as.factor(paste0(type,"_",subtype))) %>%
  mutate(POS_mb = as.numeric(as.character(POS_mb))) %>%
  filter(!is.na(POS_mb)) %>%
  mutate(founder = factor(founder, levels = levels(founder_info$founder)))

if (verbose) {
  cat("SNPINDEL after filtering:", nrow(SNPINDEL), "\n")
  if(nrow(SNPINDEL) > 0) {
    cat("SNPINDEL types:\n")
    print(table(SNPINDEL$type, SNPINDEL$subtype))
  }
}

# Debug: Let's look at the high impact INDELs specifically
high_indels <- region_variants %>%
  filter(type == "INDEL", subtype == "high")

if (verbose) {
  cat("High impact INDELs found:", nrow(high_indels), "\n")
  if(nrow(high_indels) > 0) {
    cat("High impact INDEL details:\n")
    print(high_indels[, c("CHROM", "POS", "type", "subtype", "founder", "genotype")])
    cat("Genotype distribution for high impact INDELs:\n")
    print(table(high_indels$genotype))
  }
}

# Debug: Let's look at the high impact SNPs and INDELs
high_impact_variants <- region_variants %>%
  filter(type %in% c("INDEL", "SNP"), subtype == "high")

if (verbose) {
  cat("High impact SNPs and INDELs found:", nrow(high_impact_variants), "\n")
  if(nrow(high_impact_variants) > 0) {
    cat("High impact variant details:\n")
    print(high_impact_variants[, c("CHROM", "POS", "type", "subtype", "founder", "genotype")])
    cat("Genotype distribution for high impact variants:\n")
    print(table(high_impact_variants$genotype))
  }
}

SV <- region_variants %>%
  filter(genotype == 1) %>%
  filter(type == "SV") %>%
  left_join(founder_info, by = "founder") %>%
  mutate(POS_mb = as.numeric(as.character(POS_mb))) %>%
  filter(!is.na(POS_mb)) %>%
  mutate(founder = factor(founder, levels = levels(founder_info$founder)))

# Apply the function to your SV data
SV <- assign_sv_labels(SV)

# Create the plot
p <- ggplot() +
  # Add horizontal lines for each founder
  geom_hline(
    data = founder_info,
    aes(yintercept = as.numeric(founder), color = founder),
    size = 0.5
  ) +
  
  # Add SNP and INDEL markers
  geom_segment(
    data = SNPINDEL,
    aes(
      x = POS_mb, 
      xend = POS_mb, 
      y = as.numeric(founder) - 0.2, 
      yend = as.numeric(founder) + 0.2,
      color = v_type
    ),
    size = 0.5
  ) +
  
  # Add SV markers
  geom_point(
    data = SV[!is.na(SV$subtype), ],
    aes(
      x = POS_mb,
      y = as.numeric(founder) + 0.3,
      shape = subtype,
      fill = subtype
    ),
    size = 5,
    stroke = 0.5,
    color = "black"
  ) +
  
  # Add SV labels
  geom_text(
    data = SV[!is.na(SV$subtype), ],
    aes(
      x = POS_mb,
      y = as.numeric(founder) + 0.3,
      label = label
    ),
    size = 2,
    color = "black"
  ) +
  
  # Set color scale for SNPs, INDELs, and founders
  scale_color_manual(
    values = c(
      "SNP_high" = "darkblue",
      "SNP_moderate" = "lightblue",
      "INDEL_high" = "darkgreen",
      "INDEL_moderate" = "lightgreen",
      setNames(founder_info$color_palette, founder_info$founder)
    ),
    name = "Variant Type / Founder"
  ) +
  
  # Set shape and fill scales for SVs
  scale_shape_manual(
    values = c(
      "DEL" = 24, "DEL:COM" = 24, "DEL:rME" = 24,
      "INS" = 25, "INS:COM" = 25, "INS:ME" = 25, 
      "INS:CNV" = 25, "INS:nCNV" = 25,
      "INV" = 95
    ),
    name = "SV Type"
  ) +
  scale_fill_manual(
    values = c(
      "DEL" = "orange", "DEL:COM" = "lightyellow", "DEL:rME" = "red",
      "INS" = "orange", "INS:COM" = "lightyellow", "INS:ME" = "red", 
      "INS:CNV" = "pink", "INS:nCNV" = "pink",
      "INV" = "black"
    ),
    name = "SV Type"
  ) +
  
  # Set y-axis to show founder names
  scale_y_discrete(limits = levels(founder_info$founder)) +
  
  # Set labels and theme
  labs(
    x = "Genomic Position (Mb)",
    y = "Founder"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  )

# Add custom segments for INS:nCNV and INV
for (i in 1:nrow(SV)) {
  sv <- SV[i,]
  if (!is.na(sv$subtype)) {
    y_pos <- as.numeric(sv$founder) + 0.3
    
    if (sv$subtype == "INS:nCNV") {
      p <- p + 
      geom_segment(
        aes(x = sv$POS_mb - 5000/1e6, xend = sv$POS_mb + 5000/1e6, y = y_pos + 0.1, yend = y_pos + 0.1),
        color = "pink", size = 1
      ) +
      geom_segment(
        aes(x = sv$POS_mb - 5000/1e6, xend = sv$POS_mb + 5000/1e6, y = y_pos + 0.15, yend = y_pos + 0.15),
        color = "pink", size = 1
      )
    } else if (sv$subtype == "INV") {
      p <- p + geom_segment(
        aes(x = sv$POS_mb - 10000/1e6, xend = sv$POS_mb + 10000/1e6, y = y_pos, yend = y_pos),
        color = "black", size = 1
      ) +
      geom_segment(
        aes(x = sv$POS_mb - 10000/1e6, xend = sv$POS_mb + 10000/1e6, y = y_pos - 0.2, yend = y_pos + 0.2),
        color = "grey50", size = 0.5
      ) +
      geom_segment(
        aes(x = sv$POS_mb - 10000/1e6, xend = sv$POS_mb + 10000/1e6, y = y_pos + 0.2, yend = y_pos - 0.2),
        color = "grey50", size = 0.5
      )
    }
  }
}

return(p)
}

#############
### now describe the variants
#############

#' Plot structural variants by size
#' 
#' Creates a specialized visualization of structural variants (SVs) showing their relative sizes
#' within a genomic region. This function displays SVs as triangles with heights proportional to
#' their size, with deletions shown below the center line and insertions above. Different SV types
#' are color-coded and labeled for easy identification.
#' 
#' @param variants Data frame containing variant information with REF and ALT columns for size calculation
#' @param target_chr Character string specifying the chromosome to analyze (e.g., "chr2L")
#' @param target_start Integer specifying the start position in base pairs
#' @param target_stop Integer specifying the stop position in base pairs
#' @param df2 Data frame containing founder frequency data to determine founder genotypes
#' @param reference_strain Optional character string specifying a reference strain (not used in this plot)
#' @return A ggplot object showing structural variants with sizes represented by triangle heights
#' @export
#' @importFrom ggplot2 ggplot geom_hline geom_polygon geom_text scale_fill_manual labs theme_minimal theme element_blank element_text geom_segment
#' @importFrom dplyr filter select mutate relocate left_join distinct arrange
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
XQTL_SVBySize <- function(variants, target_chr, target_start, target_stop, df2, reference_strain = NULL) {


region_variants2 <- variants %>% as_tibble() %>%
  dplyr::filter(CHROM == target_chr & POS >= target_start & POS <= target_stop) %>%
  dplyr::select(-c("ID", "QUAL", "FILTER", "INFO")) %>%  # Keep REF and ALT
  mutate(
    REF_length = nchar(REF),
    ALT_length = nchar(ALT),
    POS_mb = POS / 1e6
  ) %>%
  relocate(type, .after = POS) %>% 
  relocate(subtype, .after = type) %>%
  relocate(REF_length, .after = ALT) %>%
  relocate(ALT_length, .after = REF_length) %>%
  relocate(POS_mb, .after = POS) %>%
  pivot_longer(!c("CHROM", "POS", "POS_mb", "type", "subtype", "REF", "ALT", "REF_length", "ALT_length"), 
               names_to = "founder", values_to = "genotype") %>%
  filter(founder %in% unique(df2$founder))

# Get the color palette
color_palette <- get_palette(unique(df2$founder), reference_strain)
founder <- factor(unique(df2$founder), levels = rev(unique(df2$founder)))
founder_info <- data.frame(founder=founder, color_palette=color_palette)
founder_info$founder <- factor(founder_info$founder, levels = rev(unique(founder_info$founder)))

SV2 <- region_variants2 %>%
  filter(genotype == 1) %>%
  filter(type == "SV") %>%
  left_join(founder_info, by = "founder") %>%
  mutate(POS_mb = as.numeric(as.character(POS_mb))) %>%
  filter(!is.na(POS_mb)) %>%
  mutate(founder = factor(founder, levels = levels(founder_info$founder)))
  
# First, let's create our summary table if we haven't already
sv_summary <- SV2 %>%
  assign_sv_labels() %>%
  distinct(CHROM, POS_mb, subtype, REF_length, ALT_length, label) %>%
  arrange(POS_mb) %>%
  mutate(
    y_pos = ifelse(grepl("DEL", subtype), -1, 1),
    size = ifelse(grepl("DEL", subtype), REF_length, ALT_length),
    # Calculate height directly proportional to size
    height = size / max(size) * 0.9  # Adjust the 0.9 factor as needed
  )

# Create triangle coordinates
triangle_data <- sv_summary %>%
  mutate(
    x1 = POS_mb - size/2/1e6,
    x2 = POS_mb,
    x3 = POS_mb + size/2/1e6,
    y1 = y_pos * height,  # Base of triangle
    y2 = 0,               # Point of triangle touches the line
    y3 = y_pos * height   # Base of triangle
  ) %>%
  pivot_longer(
    cols = c(x1, x2, x3, y1, y2, y3),
    names_to = c(".value", "point"),
    names_pattern = "(.)(.)"
  ) %>%
  dplyr::select(-point)

# Get the x-axis limits from the current data
x_min <- min(region_variants2$POS_mb)
x_max <- max(region_variants2$POS_mb)

# Create the plot
p <- ggplot() +
  # Add a central line (light grey, in background)
  geom_hline(yintercept = 0, color = "lightgrey", size = 0.5) +
  
  # Add triangles for SVs
  geom_polygon(data = triangle_data, aes(x = x, y = y, fill = subtype, group = label)) +
  
  # Add labels over the triangles
  geom_text(data = sv_summary, aes(x = POS_mb, y = y_pos * height / 2, label = label), 
            size = 3, color = "black") +
  
  # Set color scale for SVs
  scale_fill_manual(
    values = c(
      "DEL" = "orange", "DEL:COM" = "lightyellow", "DEL:rME" = "red",
      "INS" = "orange", "INS:COM" = "lightyellow", "INS:ME" = "red", 
      "INS:CNV" = "pink", "INS:nCNV" = "pink",
      "INV" = "black"
    ),
    name = "SV Type"
  ) +
  
  # Set labels and theme
  labs(
    x = "Genomic Position (Mb)",
    y = "DEL / INS"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  ) +
  # Set x-axis limits to match the current data
  xlim(x_min, x_max) +
  # Adjust y-axis limits to accommodate all triangles
  ylim(-1, 1)

# Add custom segments for INS:nCNV and INV
for (i in 1:nrow(sv_summary)) {
  sv <- sv_summary[i, ]
  if (sv$subtype == "INS:nCNV") {
    p <- p + 
    geom_segment(
      aes(x = sv$POS_mb - sv$size/2/1e6, xend = sv$POS_mb + sv$size/2/1e6, y = sv$y_pos * sv$height * 0.8, yend = sv$y_pos * sv$height * 0.8),
      color = "pink", size = 1
    ) +
    geom_segment(
      aes(x = sv$POS_mb - sv$size/2/1e6, xend = sv$POS_mb + sv$size/2/1e6, y = sv$y_pos * sv$height * 0.9, yend = sv$y_pos * sv$height * 0.9),
      color = "pink", size = 1
    )
  } else if (sv$subtype == "INV") {
    p <- p + geom_segment(
      aes(x = sv$POS_mb - sv$size/2/1e6, xend = sv$POS_mb + sv$size/2/1e6, y = 0, yend = 0),
      color = "black", size = 1
    ) +
    geom_segment(
      aes(x = sv$POS_mb - sv$size/2/1e6, xend = sv$POS_mb + sv$size/2/1e6, y = -0.1, yend = 0.1),
      color = "grey50", size = 0.5
    ) +
    geom_segment(
      aes(x = sv$POS_mb - sv$size/2/1e6, xend = sv$POS_mb + sv$size/2/1e6, y = 0.1, yend = -0.1),
      color = "grey50", size = 0.5
    )
  }
}

return(p)
}


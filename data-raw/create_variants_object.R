#' Create variants object from VCF files
#' 
#' This function processes the VCF files to create a tidy variants object.
#' Run this function when you need to update the variants data.
#' 
#' @return NULL, but creates/updates the dm6.variants object
process_vcf <- function(file) {
  message(sprintf("Processing %s...", file))
  
  # Read VCF file
  vcf <- read.vcfR(file, verbose = FALSE)
  
  # Extract relevant information
  variants <- tibble(
    CHROM = getCHROM(vcf),
    POS = getPOS(vcf),
    ID = getID(vcf),
    REF = getREF(vcf),
    ALT = getALT(vcf),
    QUAL = getQUAL(vcf),
    FILTER = getFILTER(vcf),
    INFO = getINFO(vcf)
  )
  
  # Determine if this is an SV file
  is_sv_file <- grepl("SV", file)
  
  # Process genotypes
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
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
    mutate(
      type = case_when(
        grepl("snps", file) ~ "SNP",
        grepl("indels", file) ~ "INDEL",
        is_sv_file ~ "SV"
      ),
      subtype = case_when(
        grepl("high", file) & type %in% c("SNP", "INDEL") ~ "high",
        grepl("moderate", file) & type %in% c("SNP", "INDEL") ~ "moderate",
        type == "SV" ~ determine_sv_subtype(INFO),
        TRUE ~ "unknown"
      )
    )
  
  variants
}

determine_sv_subtype <- function(info) {
  case_when(
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

create_variants_object <- function(vcf_files) {
  # Process all VCF files
  variants_list <- lapply(vcf_files, process_vcf)
  
  # Combine all variants
  dm6.variants <- bind_rows(variants_list)
  
  # Save as R data
  usethis::use_data(dm6.variants, overwrite = TRUE)
  
  dm6.variants
}

# List of VCF files to process
vcf_files <- c(
  "../bitnbobs/high_impact_snps.vcf",
  "../bitnbobs/moderate_impact_snps.vcf",
  "../bitnbobs/high_impact_indels.vcf",
  "../bitnbobs/moderate_impact_indels.vcf",
  "../bitnbobs/SV.0328.vcf"
)

# Create the variants object
dm6.variants <- create_variants_object(vcf_files)

# Examine the results
print(dim(dm6.variants))
print(head(dm6.variants))
print(table(dm6.variants$type, dm6.variants$subtype))


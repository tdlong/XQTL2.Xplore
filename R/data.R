#' Processed gene annotations for Drosophila melanogaster (dm6)
#'
#' A dataset containing processed gene annotations from the dm6 NCBI RefSeq GTF file.
#' The data has been processed to include exon and UTR information, with positions
#' converted to megabase (Mb) units for easier plotting.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{seqnames}{Chromosome name (e.g., "chr2L")}
#'   \item{start}{Start position in base pairs}
#'   \item{end}{End position in base pairs}
#'   \item{start_mb}{Start position in megabases}
#'   \item{end_mb}{End position in megabases}
#'   \item{gene_name}{Gene identifier}
#'   \item{strand}{DNA strand ("+" or "-")}
#'   \item{is_utr}{Logical indicating if the region is a UTR}
#' }
#'
#' @source Processed from dm6.ncbiRefSeq.gtf
"dm6.ncbiRefSeq.genes"

#' Processed variant data for Drosophila melanogaster (dm6)
#'
#' A dataset containing processed variant information from multiple VCF files,
#' including SNPs, INDELs, and structural variants. The data has been processed
#' to include type and subtype information for easy filtering and visualization.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{CHROM}{Chromosome name (e.g., "chr2L")}
#'   \item{POS}{Position in base pairs}
#'   \item{type}{Variant type ("SNP", "INDEL", or "SV")}
#'   \item{subtype}{Variant subtype:
#'     \itemize{
#'       \item For SNPs and INDELs: "high" or "moderate" impact
#'       \item For SVs: The structural variant classification from the FL tag
#'     }
#'   }
#'   \item{Additional columns}{Genotype information for each sample}
#' }
#'
#' @source Processed from high_impact_snps.vcf, moderate_impact_snps.vcf,
#' high_impact_indels.vcf, moderate_impact_indels.vcf, and SV.0328.vcf
"dm6.variants"

#' Example QTL scan results from ZINC Hanson experiment
#'
#' This dataset contains QTL scan results for the ZINC Hanson experiment.
#' @format A tibble with 25,977 rows and 8 variables:
#' \describe{
#'   \item{chr}{Chromosome}
#'   \item{pos}{Genomic position (bp)}
#'   \item{Wald_log10p}{-log10(p) from Wald test}
#'   \item{Pseu_log10p}{-log10(p) from Pseudoscan}
#'   \item{Falc_H2}{Falconer's H2}
#'   \item{Cutl_H2}{Cutler's H2}
#'   \item{avg.var}{Average variance}
#'   \item{cM}{Genetic position (cM)}
#' }
"zinc_hanson_pseudoscan"

#' Example founder means from ZINC Hanson experiment
#'
#' This dataset contains founder mean frequencies for the ZINC Hanson experiment.
#' @format A tibble with 3,987,264 rows and 6 variables:
#' \describe{
#'   \item{chr}{Chromosome}
#'   \item{pos}{Genomic position (bp)}
#'   \item{TRT}{Treatment (C or Z)}
#'   \item{REP}{Replicate number}
#'   \item{founder}{Founder strain}
#'   \item{freq}{Allele frequency}
#' }
"zinc_hanson_means" 
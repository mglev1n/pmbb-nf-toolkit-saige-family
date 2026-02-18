#!/usr/bin/env Rscript
#
# identify_gwas_loci.R
#
# Identifies independent genome-wide significant loci from merged GWAS summary
# statistics using distance-based clumping (gwasRtools::get_loci) and annotates
# lead variants with the nearest gene (gwasRtools::get_nearest_gene).
#
# Outputs a CSV with one row per independent locus, to be collected across all
# phenotypes into saige_gwas_loci.csv by the collect_gwas_loci Nextflow process.
#

suppressPackageStartupMessages({
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  if (!requireNamespace("gwasRtools", quietly = TRUE))
    remotes::install_github("lcpilling/gwasRtools")
  if (!requireNamespace("optparse",  quietly = TRUE)) install.packages("optparse")
  if (!requireNamespace("cli",       quietly = TRUE)) install.packages("cli")
  if (!requireNamespace("dplyr",     quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("stringr",   quietly = TRUE)) install.packages("stringr")  
  if (!requireNamespace("vroom",     quietly = TRUE)) install.packages("vroom")
  if (!requireNamespace("tibble",    quietly = TRUE)) install.packages("tibble")
  if (!requireNamespace("scales",    quietly = TRUE)) install.packages("scales")
  if (!requireNamespace("rlang",     quietly = TRUE)) install.packages("rlang")
})

library(optparse)
library(cli)
library(dplyr)
library(vroom)
library(tibble)
library(stringr)

# ---------------------------------------------------------------------------
# extract_loci()
#
# Filters GWAS summary statistics for genome-wide significant SNPs, identifies
# independent loci using gwasRtools::get_loci (distance-based, no LD reference
# required), and annotates lead SNPs with nearest genes via
# gwasRtools::get_nearest_gene.
#
# Parameters
#   df            Data frame of GWAS summary statistics
#   snp_col       Column name for SNP/variant ID  (default: RSID)
#   chr_col       Column name for chromosome       (default: CHR)
#   pos_col       Column name for position         (default: POS_38)
#   maf_col       Column name for effect allele frequency (default: EAF)
#   beta_col      Column name for effect size      (default: B)
#   se_col        Column name for standard error   (default: SE)
#   p_col         Column name for p-value          (default: p_value)
#   p_threshold   Genome-wide significance threshold (default: 5e-8)
#   build         Genome build for gene annotation, 37 or 38 (default: 38)
#   ...           Additional arguments passed to gwasRtools::get_loci
#                 (e.g. distance = 500000)
#
# Returns a tibble of lead SNPs with gene annotation, or an empty tibble.
# ---------------------------------------------------------------------------
extract_loci <- function(df,
                         snp_col        = RSID,
                         chr_col        = CHR,
                         pos_col        = POS_38,
                         maf_col        = EAF,
                         beta_col       = B,
                         se_col         = SE,
                         p_col          = p_value,
                         p_threshold    = 5e-8,
                         locus_distance = 500000,
                         build          = 38) {

  snp_col  <- rlang::enquo(snp_col)
  chr_col  <- rlang::enquo(chr_col)
  pos_col  <- rlang::enquo(pos_col)
  maf_col  <- rlang::enquo(maf_col)
  beta_col <- rlang::enquo(beta_col)
  se_col   <- rlang::enquo(se_col)
  p_col    <- rlang::enquo(p_col)

  snp_col_str  <- rlang::as_name(snp_col)
  chr_col_str  <- rlang::as_name(chr_col)
  pos_col_str  <- rlang::as_name(pos_col)
  maf_col_str  <- rlang::as_name(maf_col)
  beta_col_str <- rlang::as_name(beta_col)
  se_col_str   <- rlang::as_name(se_col)
  p_col_str    <- rlang::as_name(p_col)

  required_cols <- c(snp_col_str, chr_col_str, pos_col_str,
                     maf_col_str, beta_col_str, se_col_str, p_col_str)
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {.field {missing_cols}}")
  }

  initial_rows <- nrow(df)
  cli::cli_alert_info("Starting with {.val {scales::comma(initial_rows)}} rows")

  # Warn on NA values in critical columns
  na_chr <- sum(is.na(df[[chr_col_str]]))
  na_pos <- sum(is.na(df[[pos_col_str]]))
  na_p   <- sum(is.na(df[[p_col_str]]))
  if (na_chr > 0) cli::cli_alert_warning(
    "Found {.val {scales::comma(na_chr)}} rows with missing {.field {chr_col_str}}")
  if (na_pos > 0) cli::cli_alert_warning(
    "Found {.val {scales::comma(na_pos)}} rows with missing {.field {pos_col_str}}")
  if (na_p > 0) cli::cli_alert_warning(
    "Found {.val {scales::comma(na_p)}} rows with missing {.field {p_col_str}}")

  valid_p <- df[[p_col_str]][!is.na(df[[p_col_str]])]
  if (any(valid_p < 0 | valid_p > 1)) {
    cli::cli_alert_warning("Some p-values are outside the expected range [0, 1]")
  }

  # Filter by p-value threshold
  cli::cli_progress_step("Filtering by p-value threshold ({.val {p_threshold}})")
  df_filtered <- dplyr::filter(df, !!p_col < p_threshold)

  after_p_filter <- nrow(df_filtered)
  if (after_p_filter == 0) {
    cli::cli_alert_warning(
      "No genome-wide significant variants found (p < {.val {p_threshold}})")
    return(tibble::tibble())
  }
  cli::cli_alert_info(
    "{.val {scales::comma(after_p_filter)}} variants pass p < {.val {p_threshold}}")

  # Remove missing CHR / POS
  cli::cli_progress_step("Removing rows with missing chromosome or position values")
  df_filtered <- dplyr::filter(df_filtered,
                               !is.na(!!chr_col),
                               !is.na(!!pos_col))

  if (nrow(df_filtered) == 0) {
    cli::cli_abort(
      "No rows remaining after removing missing CHR/POS. Check data quality.")
  }

  # Identify independent loci
  cli::cli_progress_step("Identifying independent loci")
  df_loci <- df_filtered |>
    gwasRtools::get_loci(snp_col     = snp_col_str,
                         chr_col     = chr_col_str,
                         pos_col     = pos_col_str,
                         maf_col     = maf_col_str,
                         beta_col    = beta_col_str,
                         se_col      = se_col_str,
                         p_col       = p_col_str,
                         use_pvalue  = TRUE,
                         n_bases     = locus_distance,
                         p_threshold = p_threshold) |>
    dplyr::filter(lead)

  n_loci <- nrow(df_loci)
  if (n_loci == 0) {
    cli::cli_alert_warning(
      "No independent loci identified. Check data and get_loci parameters.")
    return(df_loci)
  }
  cli::cli_alert_success("Found {.val {n_loci}} independent loci")

  # Annotate with nearest gene
  cli::cli_progress_step(
    "Annotating with nearest genes (build {.val {build}})")
  result <- df_loci |>
    gwasRtools::get_nearest_gene(snp_col = snp_col_str,
                                 chr_col = chr_col_str,
                                 pos_col = pos_col_str,
                                 build   = build)

  cli::cli_alert_success(
    "Complete: {.val {nrow(result)}} loci with gene annotations")

  return(tibble::as_tibble(result))
}

# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--cohort",          type = "character", help = "Cohort name"),
  make_option("--phenotype",       type = "character", help = "Phenotype name"),
  make_option("--sumstats",        type = "character", help = "Path to merged GWAS summary statistics CSV"),
  make_option("--chr_col",         type = "character", default = "CHR",     help = "Chromosome column name [default: %default]"),
  make_option("--pos_col",         type = "character", default = "POS_38",  help = "Position column name [default: %default]"),
  make_option("--snp_col",         type = "character", default = "RSID",    help = "SNP/variant ID column name [default: %default]"),
  make_option("--maf_col",         type = "character", default = "EAF",     help = "Effect allele frequency column name [default: %default]"),
  make_option("--beta_col",        type = "character", default = "B",       help = "Effect size (beta) column name [default: %default]"),
  make_option("--se_col",          type = "character", default = "SE",      help = "Standard error column name [default: %default]"),
  make_option("--p_col",           type = "character", default = "p_value", help = "P-value column name [default: %default]"),
  make_option("--build",           type = "integer",   default = 38L,       help = "Genome build for gene annotation (37 or 38) [default: %default]"),
  make_option("--p_threshold",     type = "double",    default = 5e-8,      help = "Genome-wide significance threshold [default: %default]"),
  make_option("--locus_distance",  type = "integer",   default = 500000L,   help = "Half-window distance (bp) for locus definition [default: %default]")
)

args <- parse_args(OptionParser(
  usage       = "%prog [options]",
  description = "Identify independent GWAS loci using distance-based clumping.",
  option_list = option_list
))

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
if (is.null(args$cohort))     cli::cli_abort("--cohort is required")
if (is.null(args$phenotype))  cli::cli_abort("--phenotype is required")
if (is.null(args$sumstats))   cli::cli_abort("--sumstats is required")

if (!file.exists(args$sumstats))
  cli::cli_abort("--sumstats file not found: {.path {args$sumstats}}")
if (file.size(args$sumstats) == 0)
  cli::cli_abort("--sumstats file is empty: {.path {args$sumstats}}")

if (!args$build %in% c(37L, 38L))
  cli::cli_abort("--build must be 37 or 38 (got: {args$build})")

if (args$p_threshold <= 0 || args$p_threshold >= 1)
  cli::cli_abort("--p_threshold must be between 0 and 1 (got: {args$p_threshold})")

if (args$locus_distance <= 0)
  cli::cli_abort("--locus_distance must be a positive integer (got: {args$locus_distance})")

# ---------------------------------------------------------------------------
# Load and validate data
# ---------------------------------------------------------------------------
cli::cli_h1("SAIGE GWAS Locus Identification")
cli::cli_alert_info("Cohort: {.val {args$cohort}} | Phenotype: {.val {args$phenotype}}")
cli::cli_alert_info("Reading: {.path {args$sumstats}}")

sumstats <- vroom::vroom(args$sumstats, show_col_types = FALSE)

# Confirm all specified column names exist in the data
required_cols <- c(args$snp_col, args$chr_col, args$pos_col,
                   args$maf_col, args$beta_col, args$se_col, args$p_col)
missing_cols  <- setdiff(required_cols, names(sumstats))
if (length(missing_cols) > 0) {
  cli::cli_abort(c(
    "The following columns specified via CLI arguments are not present in {.path {args$sumstats}}:",
    "x" = "{.field {missing_cols}}",
    "i" = "Available columns: {.field {names(sumstats)}}",
    "i" = "Check your gwas_col_names configuration matches the column names in the data."
  ))
}

if (is.character(sumstats[[args$chr_col]]) &&
    any(stringr::str_starts(sumstats[[args$chr_col]], "chr"), na.rm = TRUE)) {
  cli::cli_alert_info(
    "Chromosome column {.field {args$chr_col}} contains 'chr' prefix. ")
    sumstats <- sumstats |>
      dplyr::mutate(!!rlang::sym(args$chr_col) := str_remove(!!rlang::sym(args$chr_col), "chr"))
}

# ---------------------------------------------------------------------------
# Run locus identification
# ---------------------------------------------------------------------------
output_file <- paste0(args$cohort, ".", args$phenotype, ".gwas_loci.csv")

loci <- extract_loci(
  df             = sumstats,
  snp_col        = !!rlang::sym(args$snp_col),
  chr_col        = !!rlang::sym(args$chr_col),
  pos_col        = !!rlang::sym(args$pos_col),
  maf_col        = !!rlang::sym(args$maf_col),
  beta_col       = !!rlang::sym(args$beta_col),
  se_col         = !!rlang::sym(args$se_col),
  p_col          = !!rlang::sym(args$p_col),
  p_threshold    = args$p_threshold,
  locus_distance = args$locus_distance,
  build          = args$build
)

# Add cohort/phenotype columns for traceability in the collected summary table
if (nrow(loci) > 0) {
  loci <- dplyr::mutate(loci,
                        cohort    = args$cohort,
                        phenotype = args$phenotype,
                        .before   = 1)
}

readr::write_csv(loci, output_file)
cli::cli_alert_success("Written: {.path {output_file}} ({nrow(loci)} loci)")

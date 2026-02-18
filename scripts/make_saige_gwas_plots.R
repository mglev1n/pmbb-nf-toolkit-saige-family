#!/usr/bin/env Rscript
#
# make_saige_gwas_plots.R
#
# Generates Manhattan and QQ plots for a single cohort-phenotype GWAS result
# using levinmisc::gg_manhattan_df and levinmisc::gg_qq_df.
#
# Accepts an optional --loci_csv argument (path to per-phenotype loci CSV from
# identify_gwas_loci.R). When provided, up to 2 lead variants per chromosome
# are labelled on the Manhattan plot. When absent or set to the sentinel value
# "NO_FILE", the plot is produced without annotations.
#

suppressPackageStartupMessages({
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  if (!requireNamespace("levinmisc", quietly = TRUE))
    remotes::install_github("mglev1n/levinmisc")
  if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
  if (!requireNamespace("cli",      quietly = TRUE)) install.packages("cli")
  if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("stringr",    quietly = TRUE)) install.packages("stringr")
  if (!requireNamespace("vroom",    quietly = TRUE)) install.packages("vroom")
  if (!requireNamespace("readr",    quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
})

library(optparse)
library(cli)
library(dplyr)
library(vroom)
library(readr)
library(ggplot2)
library(stringr)

# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--cohort",     type = "character", help = "Cohort name"),
  make_option("--phenotype",  type = "character", help = "Phenotype name"),
  make_option("--sumstats",   type = "character", help = "Path to merged GWAS summary statistics CSV"),
  make_option("--phenoTable", type = "character", help = "Path to pheno_summaries.csv"),
  make_option("--chr_col",    type = "character", default = "CHR",     help = "Chromosome column name [default: %default]"),
  make_option("--pos_col",    type = "character", default = "POS_38",  help = "Position column name [default: %default]"),
  make_option("--maf_col",    type = "character", default = "EAF",     help = "Effect allele frequency column name [default: %default]"),
  make_option("--beta_col",   type = "character", default = "B",       help = "Effect size (beta) column name [default: %default]"),
  make_option("--se_col",     type = "character", default = "SE",      help = "Standard error column name [default: %default]"),
  make_option("--p_col",      type = "character", default = "p_value", help = "P-value column name [default: %default]"),
  make_option("--build",      type = "character", default = "hg38",    help = "Genome build for gg_manhattan_df (hg37 or hg38) [default: %default]"),
  make_option("--loci_csv",   type = "character", default = NULL,      help = "Optional path to per-phenotype loci CSV for Manhattan annotations")
)

args <- parse_args(OptionParser(
  usage       = "%prog [options]",
  description = paste("Generate GWAS Manhattan and QQ plots using levinmisc.",
                      "Optionally annotates lead variants from a loci CSV."),
  option_list = option_list
))

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
if (is.null(args$cohort))     cli::cli_abort("--cohort is required")
if (is.null(args$phenotype))  cli::cli_abort("--phenotype is required")
if (is.null(args$sumstats))   cli::cli_abort("--sumstats is required")
if (is.null(args$phenoTable)) cli::cli_abort("--phenoTable is required")

if (!file.exists(args$sumstats))
  cli::cli_abort("--sumstats file not found: {.path {args$sumstats}}")
if (file.size(args$sumstats) == 0)
  cli::cli_abort("--sumstats file is empty: {.path {args$sumstats}}")

if (!file.exists(args$phenoTable))
  cli::cli_abort("--phenoTable file not found: {.path {args$phenoTable}}")

if (!args$build %in% c("hg37", "hg38"))
  cli::cli_abort("--build must be 'hg37' or 'hg38' (got: {.val {args$build}})")

# Sentinel check: treat "NO_FILE" as absent
use_loci <- !is.null(args$loci_csv) &&
            args$loci_csv != "NO_FILE" &&
            file.exists(args$loci_csv) &&
            file.size(args$loci_csv) > 0

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cli::cli_h1("SAIGE GWAS Plots")
cli::cli_alert_info("Cohort: {.val {args$cohort}} | Phenotype: {.val {args$phenotype}}")
cli::cli_alert_info("Reading sumstats: {.path {args$sumstats}}")

sumstats <- vroom::vroom(args$sumstats, show_col_types = FALSE)

# Extract p-values for QQ now, before we filter/modify sumstats
qq_pvals <- sumstats |> dplyr::select(dplyr::all_of(args$p_col))

# Confirm required columns exist
required_cols <- c(args$chr_col, args$pos_col, args$maf_col,
                   args$beta_col, args$se_col, args$p_col)
missing_cols  <- setdiff(required_cols, names(sumstats))
if (length(missing_cols) > 0) {
  cli::cli_abort(c(
    "The following columns specified via CLI arguments are absent from {.path {args$sumstats}}:",
    "x" = "{.field {missing_cols}}",
    "i" = "Available columns: {.field {names(sumstats)}}",
    "i" = "Check your gwas_col_names configuration matches the data."
  ))
}

# ---------------------------------------------------------------------------
# Build plot subtitle from phenotype summary table
# ---------------------------------------------------------------------------
pheno_tbl <- vroom::vroom(args$phenoTable, show_col_types = FALSE)
pheno_row  <- pheno_tbl |>
  dplyr::filter(PHENO == args$phenotype, COHORT == args$cohort)

subtitle <- if (nrow(pheno_row) > 0) {
  if ("Cases" %in% names(pheno_row) && !is.na(pheno_row$Cases[1])) {
    sprintf("Cases = %s  |  Controls = %s",
            scales::comma(pheno_row$Cases[1]),
            scales::comma(pheno_row$Controls[1]))
  } else if ("N" %in% names(pheno_row) && !is.na(pheno_row$N[1])) {
    sprintf("N = %s", scales::comma(pheno_row$N[1]))
  } else {
    NULL
  }
} else {
  cli::cli_alert_warning(
    "No entry found for cohort={.val {args$cohort}}, phenotype={.val {args$phenotype}} in phenoTable")
  NULL
}

# ---------------------------------------------------------------------------
# Build Manhattan input: standard pre-filter
#   - keep p < 0.001 (reduces data volume for plotting)
#   - keep MAF 5–95% OR p < 5e-8 (retain rare significant hits)
# ---------------------------------------------------------------------------
cli::cli_progress_step("Pre-filtering data for Manhattan plot")

p_sym   <- rlang::sym(args$p_col)
maf_sym <- rlang::sym(args$maf_col)

manhattan_df <- sumstats |>
  dplyr::filter(!!p_sym < 0.01) |>
  dplyr::filter(dplyr::between(!!maf_sym, 0.01, 0.99) | !!p_sym < 5e-8)

rm(sumstats)
gc()

if (is.character(manhattan_df[[args$chr_col]]) &&
    any(stringr::str_starts(manhattan_df[[args$chr_col]], "chr"), na.rm = TRUE)) {
  cli::cli_alert_info(
    "Chromosome column {.field {args$chr_col}} contains 'chr' prefix. ")
    manhattan_df <- manhattan_df |>
      dplyr::mutate(!!rlang::sym(args$chr_col) := str_remove(!!rlang::sym(args$chr_col), "chr"))
}

cli::cli_alert_info(
  "{.val {scales::comma(nrow(manhattan_df))}} variants retained for Manhattan plot")

# ---------------------------------------------------------------------------
# Build annotation data frame from loci CSV (optional)
# ---------------------------------------------------------------------------
annotation_df <- NULL

if (use_loci) {
  cli::cli_alert_info("Reading loci for annotations: {.path {args$loci_csv}}")
  loci <- vroom::vroom(args$loci_csv, show_col_types = FALSE)

  if (nrow(loci) > 0) {
    p_sym_loci <- rlang::sym(args$p_col)
    chr_sym    <- rlang::sym(args$chr_col)

    annotation_df <- loci |>
      dplyr::slice_min(!!p_sym_loci, n = 2, with_ties = FALSE, by = !!chr_sym)

    cli::cli_alert_info(
      "Using {.val {nrow(annotation_df)}} lead variants for Manhattan annotations")
  } else {
    cli::cli_alert_warning("Loci CSV is empty — plotting without annotations")
  }
} else {
  cli::cli_alert_info("No loci CSV provided — plotting without lead-variant annotations")
}

# ---------------------------------------------------------------------------
# Output file paths
# ---------------------------------------------------------------------------
out_manhattan <- paste0(args$cohort, ".", args$phenotype, ".manhattan.png")
out_qq        <- paste0(args$cohort, ".", args$phenotype, ".qq.png")
out_manifest  <- paste0(args$cohort, ".", args$phenotype, ".gwas.plots_manifest.csv")

# ---------------------------------------------------------------------------
# Manhattan plot
# ---------------------------------------------------------------------------
cli::cli_progress_step("Generating Manhattan plot")

chr_sym  <- rlang::sym(args$chr_col)
pos_sym  <- rlang::sym(args$pos_col)
beta_sym <- rlang::sym(args$beta_col)
se_sym   <- rlang::sym(args$se_col)
p_sym    <- rlang::sym(args$p_col)

manhattan_plot <- levinmisc::gg_manhattan_df(
  manhattan_df,
  annotation_df = annotation_df,
  chr_col       = !!chr_sym,
  pos_col       = !!pos_sym,
  beta_col      = !!beta_sym,
  se_col        = !!se_sym,
  pval_col      = !!p_sym,
  build         = args$build
)

if (!is.null(subtitle)) {
  manhattan_plot <- manhattan_plot +
    ggplot2::labs(
      title    = paste0("SAIGE GWAS ", args$cohort, ": ",
                        gsub("_", " ", args$phenotype)),
      subtitle = subtitle
    )
}

ggplot2::ggsave(out_manhattan, manhattan_plot,
                width = 14, height = 8, dpi = 300, bg = "white")
cli::cli_alert_success("Saved Manhattan plot: {.path {out_manhattan}}")
rm(manhattan_df, annotation_df, manhattan_plot)
gc()

# ---------------------------------------------------------------------------
# QQ plot
# ---------------------------------------------------------------------------
cli::cli_progress_step("Generating QQ plot")

qq_plot <- levinmisc::gg_qq_df(qq_pvals, pval_col = !!p_sym)

ggplot2::ggsave(out_qq, qq_plot,
                width = 6, height = 6, dpi = 300, bg = "white")
cli::cli_alert_success("Saved QQ plot: {.path {out_qq}}")

# ---------------------------------------------------------------------------
# Plots manifest CSV
# ---------------------------------------------------------------------------
manifest <- tibble::tibble(
  cohort         = args$cohort,
  phenotype      = args$phenotype,
  results_file   = args$sumstats,
  gwas_manhattan = file.path("Plots", out_manhattan),
  gwas_qq        = file.path("Plots", out_qq)
)

readr::write_csv(manifest, out_manifest)
cli::cli_alert_success("Saved manifest: {.path {out_manifest}}")

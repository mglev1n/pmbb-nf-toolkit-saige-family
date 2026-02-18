# SAIGE GWAS Pipeline — Improvement Plan

## Overview

This document describes planned improvements to the `saige_gwas.nf` pipeline, focused on
three areas:

1. **Modularity** — decouple locus identification, plotting, and biofilter annotation into
   independently optional steps; reduce the core pipeline to GWAS execution + summary
   statistic output
2. **Locus identification** — replace the current approach (p-value threshold only, no
   clumping) with distance-based independent locus identification and nearest-gene
   annotation using `gwasRtools::get_loci` and `gwasRtools::get_nearest_gene`
3. **Visualisation** — replace the Python-based `manhattan_plot` module with
   `levinmisc::gg_manhattan_df` and `levinmisc::gg_qq_df`
4. **Input validation** — add upfront parameter and file-existence checks at both the
   Nextflow workflow level and within individual R scripts, so failures are caught early
   with clear error messages

Scope: **`saige_gwas.nf` only**. ExWAS, PheWAS, and Gene PheWAS workflows are out of
scope for this iteration.

---

## Background and Motivation

### Current limitations

| Issue | Detail |
|---|---|
| Locus identification | None. `make_summary_suggestive_gwas` concatenates all variants below `p_cutoff_summarize` and sorts by position. Correlated variants in the same region are all reported. |
| Manhattan/QQ plotting | Uses Python `manhattan_plot` module (`make_saige_gwas_singles_plots.py`). An external dependency not in the container; unmaintained upstream; does not match lab plotting conventions. |
| Modularity | Plotting is always run. Biofilter annotation (`params.annotate`) is the only optional step, but it is also coupled to the plotting path via a hard if/else branch: `annotate=true` forces a different plotting script (`make_saige_gwas_plots_with_annot`). |
| Input validation | `check_input_genetic_data_parameters()` covers genetic data paths and `ftype`, but many other required parameters (e.g. `data_csv`, `cohort_sets`, `p_cutoff_summarize`, column name keys) are not validated upfront. Errors surface late and with cryptic messages. |

### Target state

```
SAIGE_PREPROCESSING → SAIGE_STEP1 → SAIGE_VAR_STEP2 → merge_and_filter  [CORE — always runs]
                                                              │
                    ┌─────────────────┬────────────────────┬─┴──────────────────┐
                    ▼                 ▼                    ▼                    ▼
            identify_loci       make_plots           annotate (biofilter)  [all optional,
          (params.identify_   (params.make_        (params.annotate)        independent]
              loci=true)        plots=true)          (default false)
                    │                 │
         saige_gwas_loci.csv    Manhattan + QQ PNGs
         (replaces                (R/ggplot2)
      saige_gwas_suggestive.csv)
```

The biofilter annotation path is retained as a standalone optional step. It no longer
affects or determines which plotting script is called.

---

## R Packages

Both packages are public GitHub packages and will be installed at runtime inside the
container via `remotes::install_github()`.

| Package | Repository | Purpose |
|---|---|---|
| `gwasRtools` | `lcpilling/gwasRtools` | `get_loci()` for distance-based locus clumping; `get_nearest_gene()` for gene annotation |
| `levinmisc` | `mglev1n/levinmisc` | `gg_manhattan_df()` for Manhattan plots; `gg_qq_df()` for QQ plots |

Version pinning: always-latest-`main` is acceptable for now.

Runtime install pattern (both R scripts):
```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("gwasRtools", quietly = TRUE))
    remotes::install_github("lcpilling/gwasRtools")
if (!requireNamespace("levinmisc", quietly = TRUE))
    remotes::install_github("mglev1n/levinmisc")
```

---

## Column Name Contract

SAIGE raw output columns and their post-`gwas_col_names`-mapping equivalents used by the
R functions are:

| SAIGE raw column | Typical mapped name | R function default arg |
|---|---|---|
| `CHR` | `CHR` | `chr_col = CHR` |
| `POS` | `POS_38` | `pos_col = POS_38` |
| `MarkerID` | `RSID` | `snp_col = RSID` |
| `AF_Allele2` | `EAF` | `maf_col = EAF` |
| `BETA` | `B` | `beta_col = B` |
| `SE` | `SE` | `se_col = SE` |
| `p.value` | `p_value` | `p_col = p_value` |

The Nextflow process shell block will derive the actual column names from
`params.gwas_col_names` and pass them explicitly as CLI arguments to both R scripts.
This eliminates reliance on R function defaults and makes the column mapping traceable
from the pipeline configuration.

Example (in process shell block):
```groovy
def chr_col  = params.gwas_col_names.containsKey('CHR')        ? params.gwas_col_names['CHR']        : 'CHR'
def pos_col  = params.gwas_col_names.containsKey('POS')        ? params.gwas_col_names['POS']        : 'POS_38'
def snp_col  = params.gwas_col_names.containsKey('MarkerID')   ? params.gwas_col_names['MarkerID']   : 'RSID'
def maf_col  = params.gwas_col_names.containsKey('AF_Allele2') ? params.gwas_col_names['AF_Allele2'] : 'EAF'
def beta_col = params.gwas_col_names.containsKey('BETA')       ? params.gwas_col_names['BETA']       : 'B'
def se_col   = params.gwas_col_names.containsKey('SE')         ? params.gwas_col_names['SE']         : 'SE'
def p_col    = params.gwas_col_names.containsKey('p.value')    ? params.gwas_col_names['p.value']    : 'p_value'
```

---

## New Parameters

The following parameters are added. All have sensible defaults so existing configs do not
require modification unless users want to override them.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `params.identify_loci` | Boolean | `true` | Run locus identification and output loci table |
| `params.make_plots` | Boolean | `true` | Generate Manhattan and QQ plots |
| `params.genome_build` | Integer | `38` | Genome build for `get_nearest_gene()` (37 or 38) |
| `params.gwas_locus_distance` | Integer | `500000` | Half-window distance (bp) for `get_loci()` |

`params.annotate` already exists (default `false`) and is unchanged in behaviour; it now
operates as a fully independent optional path.

---

## File Changes

### New files

| File | Description |
|---|---|
| `scripts/identify_gwas_loci.R` | Runs `extract_loci()`, outputs per-phenotype loci CSV |
| `scripts/make_saige_gwas_plots.R` | Runs `gg_manhattan_df()` + `gg_qq_df()`, accepts optional loci CSV for annotations |

### Modified files

| File | Nature of change |
|---|---|
| `processes/saige_postprocessing.nf` | Add `identify_gwas_loci` and `collect_gwas_loci` processes; remove `make_summary_suggestive_gwas` |
| `processes/saige_visualization.nf` | Add `make_saige_gwas_plots_R` process; retain Python processes in file but remove from workflow wiring |
| `processes/saige_helpers.nf` | Add `validate_gwas_params()`; register two new R scripts in `get_script_file_names()` |
| `workflows/saige_gwas.nf` | Call validation first; restructure into independent optional blocks; remove annotate/plot coupling |
| `READMEs/SAIGE_GWAS_docs.md` | Document new parameters, new output files, and modular flags |

---

## Detailed Changes

### 1. `scripts/identify_gwas_loci.R` (new)

**Purpose:** Reads merged GWAS summary statistics, identifies independent loci using
distance-based clumping, and annotates lead variants with nearest genes.

**Command-line interface:**
```
Rscript identify_gwas_loci.R \
  --cohort      <cohort>           \
  --phenotype   <pheno>            \
  --sumstats    <path/to/merged.csv> \
  --chr_col     <col_name>         \
  --pos_col     <col_name>         \
  --snp_col     <col_name>         \
  --maf_col     <col_name>         \
  --beta_col    <col_name>         \
  --se_col      <col_name>         \
  --p_col       <col_name>         \
  --build       <37|38>            \
  --p_threshold <numeric>          \
  --locus_distance <integer>
```

**Script structure:**
```
1. Parse args (optparse or argparse)
2. Install / load gwasRtools, levinmisc, tidyverse, cli
3. VALIDATE INPUTS
   - sumstats file exists and is non-empty
   - build is 37 or 38
   - p_threshold is in (0, 1)
   - locus_distance > 0
4. Read sumstats CSV
5. Validate that all specified column name args exist as columns in the data
6. Call extract_loci(df, snp_col, chr_col, pos_col, maf_col, beta_col, se_col,
                     p_col, p_threshold, build,
                     distance = locus_distance)
7. If 0 rows returned: write an empty CSV with the correct headers and exit 0
   (downstream collector handles empty inputs gracefully)
8. Write {cohort}.{pheno}.gwas_loci.csv
```

**Output:** `{cohort}.{pheno}.gwas_loci.csv` — one row per independent locus lead
variant, with columns from `extract_loci()` plus `gene` and `dist` from
`get_nearest_gene()`.

**`extract_loci()` is defined in the existing lab codebase** (shown in the planning
session). It wraps `gwasRtools::get_loci()` with input validation, NA reporting, and
progress messaging via `cli`. The R script calls it directly.

---

### 2. `scripts/make_saige_gwas_plots.R` (new)

**Purpose:** Generates Manhattan and QQ plots for a single cohort-phenotype combination.
Accepts an optional loci CSV for annotating lead variants on the Manhattan plot.

**Command-line interface:**
```
Rscript make_saige_gwas_plots.R \
  --cohort      <cohort>             \
  --phenotype   <pheno>              \
  --sumstats    <path/to/merged.csv> \
  --phenoTable  <path/to/pheno_summaries.csv> \
  --chr_col     <col_name>           \
  --pos_col     <col_name>           \
  --beta_col    <col_name>           \
  --se_col      <col_name>           \
  --p_col       <col_name>           \
  --maf_col     <col_name>           \
  --build       <hg37|hg38>          \
  [--loci_csv   <path/to/loci.csv>]
```

`--loci_csv` is optional. When absent (or when the sentinel path `NO_FILE` is passed),
the Manhattan plot is produced without lead-variant annotations.

**Script structure:**
```
1. Parse args
2. Install / load levinmisc, tidyverse, cli, scales
3. VALIDATE INPUTS
   - sumstats and phenoTable files exist and are non-empty
   - build is hg37 or hg38
   - if loci_csv provided and is not NO_FILE: check it exists
4. Read sumstats; validate that all column name args exist in the data
5. Read phenoTable; extract N/Cases/Controls for plot subtitle
6. BUILD MANHATTAN DATA
   - Filter: p_col < 0.001 AND (maf_col between 0.05–0.95 OR p_col < 5e-8)
   - This matches the lab's standard pre-filter for gg_manhattan_df
7. BUILD ANNOTATION DF (for lead-variant labels)
   - If loci_csv provided: read and select top 2 lead variants per chromosome
     by p-value (slice_min(p_col, n = 2, with_ties = FALSE, by = chr_col))
   - Otherwise: annotation_df = NULL
8. MANHATTAN PLOT
   levinmisc::gg_manhattan_df(
     filtered_df,
     annotation_df = annotation_df,
     chr_col = chr_col, pos_col = pos_col,
     beta_col = beta_col, se_col = se_col,
     pval_col = p_col, build = build
   )
   Save as {cohort}.{pheno}.manhattan_vertical.png
9. QQ PLOT
   levinmisc::gg_qq_df(full_sumstats_df, pval_col = p_col)
   Save as {cohort}.{pheno}.qq.png
10. MANIFEST
    Write {cohort}.{pheno}.gwas.plots_manifest.csv with columns:
    phenotype, cohort, results_file, gwas_manhattan, gwas_qq
    (matching current manifest format for HTML report compatibility)
```

---

### 3. `processes/saige_postprocessing.nf`

#### Remove: `make_summary_suggestive_gwas`

This process (lines 332–359) concatenates all filtered CSVs across phenotypes and sorts
by CHR/POS. It is replaced by `collect_gwas_loci` below.

#### Add: `identify_gwas_loci`

```nextflow
process identify_gwas_loci {
    publishDir "${launchDir}/${cohort}/Sumstats/"
    label 'safe_to_skip'

    input:
        tuple val(cohort), val(pheno), path(sumstats)
        path(loci_script)

    output:
        tuple val(cohort), val(pheno), path("${cohort}.${pheno}.gwas_loci.csv")

    shell:
        // Derive column names from gwas_col_names param values
        def chr_col  = params.gwas_col_names.containsKey('CHR')        ? params.gwas_col_names['CHR']        : 'CHR'
        def pos_col  = params.gwas_col_names.containsKey('POS')        ? params.gwas_col_names['POS']        : 'POS_38'
        def snp_col  = params.gwas_col_names.containsKey('MarkerID')   ? params.gwas_col_names['MarkerID']   : 'RSID'
        def maf_col  = params.gwas_col_names.containsKey('AF_Allele2') ? params.gwas_col_names['AF_Allele2'] : 'EAF'
        def beta_col = params.gwas_col_names.containsKey('BETA')       ? params.gwas_col_names['BETA']       : 'B'
        def se_col   = params.gwas_col_names.containsKey('SE')         ? params.gwas_col_names['SE']         : 'SE'
        def p_col    = params.gwas_col_names.containsKey('p.value')    ? params.gwas_col_names['p.value']    : 'p_value'
        """
        Rscript ${loci_script} \
          --cohort           ${cohort} \
          --phenotype        ${pheno} \
          --sumstats         ${sumstats} \
          --chr_col          ${chr_col} \
          --pos_col          ${pos_col} \
          --snp_col          ${snp_col} \
          --maf_col          ${maf_col} \
          --beta_col         ${beta_col} \
          --se_col           ${se_col} \
          --p_col            ${p_col} \
          --build            ${params.genome_build} \
          --p_threshold      ${params.p_cutoff_summarize} \
          --locus_distance   ${params.gwas_locus_distance}
        """
    stub:
        """
        touch ${cohort}.${pheno}.gwas_loci.csv
        """
}
```

#### Add: `collect_gwas_loci`

Replaces `make_summary_suggestive_gwas`. Collects all per-phenotype loci CSVs into a
single sorted summary table.

```nextflow
process collect_gwas_loci {
    publishDir "${launchDir}/Summary/"
    label 'safe_to_skip'

    input:
        path(loci_files)

    output:
        path('saige_gwas_loci.csv')

    script:
        """
        #! ${params.my_python}

        import pandas as pd
        import glob

        input_list = [x.strip() for x in "${loci_files}".replace("[","").replace("]","").split()]
        dfs = [pd.read_csv(f) for f in input_list if pd.read_csv(f).shape[0] > 0]

        if dfs:
            pd.concat(dfs).sort_values(by=['CHR', 'POS_38']).to_csv('saige_gwas_loci.csv', index=False)
        else:
            pd.DataFrame().to_csv('saige_gwas_loci.csv', index=False)
        """
    stub:
        '''
        touch saige_gwas_loci.csv
        '''
}
```

**Note on column names in `collect_gwas_loci`:** The CHR/POS_38 sort columns should be
derived from the mapped column names (same pattern as `identify_gwas_loci`) rather than
hard-coded. This is flagged for implementation.

#### Modify: `make_summary_table_with_annot`

Currently outputs `saige_gwas_suggestive.csv` (same filename as
`make_summary_suggestive_gwas`). With the rename to `saige_gwas_loci.csv` for the loci
output, rename this to `saige_gwas_biofilter_annotated.csv` to avoid the clash. The
process logic is otherwise unchanged.

---

### 4. `processes/saige_visualization.nf`

#### Add: `make_saige_gwas_plots_R`

```nextflow
process make_saige_gwas_plots_R {
    publishDir "${launchDir}/Plots/"
    label 'safe_to_skip', 'high_memory_plots'
    memory {
        def fileSizeGb = sumstats.size() / (1024**3)
        def requiredMemGb = Math.max(8, Math.ceil(fileSizeGb * 10))
        return Math.min(requiredMemGb, 256).GB
    }

    input:
        tuple val(cohort), val(pheno), path(sumstats), path(loci_csv)
        path(plot_script)
        path(pheno_table)

    output:
        path "${cohort}.${pheno}.{manhattan_vertical.png,qq.png,gwas.plots_manifest.csv}"

    shell:
        def chr_col  = params.gwas_col_names.containsKey('CHR')        ? params.gwas_col_names['CHR']        : 'CHR'
        def pos_col  = params.gwas_col_names.containsKey('POS')        ? params.gwas_col_names['POS']        : 'POS_38'
        def maf_col  = params.gwas_col_names.containsKey('AF_Allele2') ? params.gwas_col_names['AF_Allele2'] : 'EAF'
        def beta_col = params.gwas_col_names.containsKey('BETA')       ? params.gwas_col_names['BETA']       : 'B'
        def se_col   = params.gwas_col_names.containsKey('SE')         ? params.gwas_col_names['SE']         : 'SE'
        def p_col    = params.gwas_col_names.containsKey('p.value')    ? params.gwas_col_names['p.value']    : 'p_value'
        def build_str = "hg${params.genome_build}"
        def loci_arg  = (loci_csv.name != 'NO_FILE') ? "--loci_csv ${loci_csv}" : ""
        """
        Rscript ${plot_script} \
          --cohort      ${cohort} \
          --phenotype   ${pheno} \
          --sumstats    ${sumstats} \
          --phenoTable  ${pheno_table} \
          --chr_col     ${chr_col} \
          --pos_col     ${pos_col} \
          --maf_col     ${maf_col} \
          --beta_col    ${beta_col} \
          --se_col      ${se_col} \
          --p_col       ${p_col} \
          --build       ${build_str} \
          ${loci_arg}
        """
    stub:
        """
        touch ${cohort}.${pheno}.manhattan_vertical.png
        touch ${cohort}.${pheno}.qq.png
        touch ${cohort}.${pheno}.gwas.plots_manifest.csv
        """
}
```

#### Retain (but unwire): Python plotting processes

`make_saige_gwas_plots` and `make_gwas_plots_with_annot` are kept in the file but are no
longer included in any workflow `include` block. This preserves them for reference during
transition and allows easy rollback if needed.

---

### 5. `processes/saige_helpers.nf`

#### Add: `validate_gwas_params()`

A new validation function called as the first statement inside `workflow SAIGE_GWAS`.
Uses `error()` (not `System.exit(1)` or `println`) for Nextflow-idiomatic failure
reporting.

```groovy
def validate_gwas_params(params) {

    // --- Required: data inputs ---
    if (!params.data_csv)
        error "params.data_csv must be set"
    if (!file(params.data_csv).exists())
        error "params.data_csv file not found: ${params.data_csv}"

    if (!params.cohort_sets)
        error "params.cohort_sets must be set"
    if (!file(params.cohort_sets).exists())
        error "params.cohort_sets file not found: ${params.cohort_sets}"

    // --- Required: at least one phenotype list is non-empty ---
    def has_phenos = (paramToList(params.bin_pheno_list)
                    + paramToList(params.quant_pheno_list)
                    + paramToList(params.survival_pheno_list)).size() > 0
    if (!has_phenos)
        error "At least one of bin_pheno_list, quant_pheno_list, or survival_pheno_list must be non-empty"

    // --- Required: SAIGE scripts ---
    if (!params.step1_script)
        error "params.step1_script must be set"
    if (!file(params.step1_script).exists())
        error "params.step1_script not found: ${params.step1_script}"

    if (!params.step2_script)
        error "params.step2_script must be set"
    if (!file(params.step2_script).exists())
        error "params.step2_script not found: ${params.step2_script}"

    // --- Required: p_cutoff_summarize ---
    if (params.p_cutoff_summarize == null)
        error "params.p_cutoff_summarize must be set"
    if (params.p_cutoff_summarize <= 0 || params.p_cutoff_summarize >= 1)
        error "params.p_cutoff_summarize must be between 0 and 1 (got: ${params.p_cutoff_summarize})"

    // --- Required: gwas_col_names must include a p-value key ---
    if (!params.gwas_col_names || !params.gwas_col_names.containsKey('p.value'))
        error "params.gwas_col_names must include an entry for 'p.value' (the raw SAIGE p-value column)"

    // --- New optional-step params ---
    if (params.genome_build != 37 && params.genome_build != 38)
        error "params.genome_build must be 37 or 38 (got: ${params.genome_build})"

    if (params.gwas_locus_distance <= 0)
        error "params.gwas_locus_distance must be a positive integer (got: ${params.gwas_locus_distance})"

    // --- Conditional: biofilter annotation ---
    if (params.annotate) {
        if (!params.biofilter_loki)
            error "params.annotate is true but params.biofilter_loki is not set"
        if (!file(params.biofilter_loki).exists())
            error "params.biofilter_loki not found: ${params.biofilter_loki}"
        if (!params.biofilter_script)
            error "params.annotate is true but params.biofilter_script is not set"
        if (!file(params.biofilter_script).exists())
            error "params.biofilter_script not found: ${params.biofilter_script}"
    }

    // Genetic data params are already validated by check_input_genetic_data_parameters()
    // which is called separately
}
```

#### Modify: `get_script_file_names()`

Add entries for the two new R scripts:

```groovy
script_names['gwas_loci']  = "${moduleDir}/../scripts/identify_gwas_loci.R"
script_names['gwas_plots_r'] = "${moduleDir}/../scripts/make_saige_gwas_plots.R"
```

The existing `gwas_plots` and `gwas_plots_with_annot` entries are retained in the map
for reference but are no longer referenced by any workflow.

---

### 6. `workflows/saige_gwas.nf`

#### Update `include` block

Remove:
```groovy
include { make_saige_gwas_plots; make_gwas_plots_with_annot } from '../processes/saige_visualization.nf'
```

Add:
```groovy
include { make_saige_gwas_plots_R } from '../processes/saige_visualization.nf'
```

Remove from postprocessing includes:
```groovy
make_summary_suggestive_gwas
```

Add to postprocessing includes:
```groovy
identify_gwas_loci
collect_gwas_loci
```

#### Update `log.info` block

Add new parameters to the startup log:

```groovy
    Postprocessing Options
    ==================================================
    identify_loci           : ${params.identify_loci}
    make_plots              : ${params.make_plots}
    annotate                : ${params.annotate}
    genome_build            : ${params.genome_build}
    gwas_locus_distance     : ${params.gwas_locus_distance}
```

#### Restructure `workflow SAIGE_GWAS`

```groovy
workflow SAIGE_GWAS {
  main:
    // --- VALIDATION (runs before any channels/processes) ---
    validate_gwas_params(params)

    // ... existing chromosome/ftype/script setup unchanged ...

    // --- CORE ---
    // SAIGE_PREPROCESSING, SAIGE_STEP1, SAIGE_VAR_STEP2 unchanged

    (singles_merge_output, filtered_singles_output) =
        merge_and_filter_saige_gwas_output(step2_grouped_output, merge_singles_script)

    // --- OPTIONAL: locus identification ---
    if (params.identify_loci) {
        loci_script   = script_name_dict['gwas_loci']
        loci_per_pheno = identify_gwas_loci(singles_merge_output, loci_script)
        collect_gwas_loci(loci_per_pheno.map { c, p, f -> f }.collect())
    }

    // --- OPTIONAL: plotting ---
    if (params.make_plots) {
        plot_script = script_name_dict['gwas_plots_r']

        // Attach loci CSV to each sumstats tuple, or use sentinel if loci not run
        plot_input = params.identify_loci
            ? singles_merge_output.join(loci_per_pheno, by: [0, 1])
            : singles_merge_output.map { c, p, s -> tuple(c, p, s, file('NO_FILE')) }

        make_saige_gwas_plots_R(plot_input, plot_script, pheno_table)
    }

    // --- OPTIONAL: biofilter annotation (standalone, no plot coupling) ---
    if (params.annotate) {
        biofilter_input = gwas_make_biofilter_positions_input(
            filtered_singles_output.map { c, p, f -> f }.collect()
        )
        BIOFILTER_POSITIONS(Channel.of('GWAS').combine(biofilter_input))
        make_summary_table_with_annot(
            filtered_singles_output.map { c, p, f -> f }.collect(),
            BIOFILTER_POSITIONS.out
        )
    }

    json_params = dump_params_to_json(params, 'saige_gwas')

  emit:
    singles_merge_output
    pheno_table
}
```

---

## Output File Changes

| Old output | New output | Notes |
|---|---|---|
| `Summary/saige_gwas_suggestive.csv` | `Summary/saige_gwas_loci.csv` | One row per independent locus (lead variant only), with nearest gene annotation. Replaces all-variants-below-threshold table. |
| `{cohort}/Sumstats/` | `{cohort}/Sumstats/{cohort}.{pheno}.gwas_loci.csv` | Per-phenotype loci file (new, published here before collection) |
| `Summary/saige_gwas_suggestive.csv` (annotate path) | `Summary/saige_gwas_biofilter_annotated.csv` | Biofilter-annotated summary; renamed to avoid clash with loci output |
| `Plots/*.manhattan_vertical.png` | unchanged | Now generated by R/ggplot2 |
| `Plots/*.qq.png` | unchanged | Now generated by R/ggplot2 |

---

## Input Validation Summary

### Nextflow level (`validate_gwas_params()` in `saige_helpers.nf`)

| Check | Condition triggering error |
|---|---|
| `data_csv` set and exists | param null or file missing |
| `cohort_sets` set and exists | param null or file missing |
| At least one phenotype list non-empty | all three lists empty |
| `step1_script` set and exists | param null or file missing |
| `step2_script` set and exists | param null or file missing |
| `p_cutoff_summarize` in (0, 1) | null or out of range |
| `gwas_col_names` has `p.value` key | missing key |
| `genome_build` is 37 or 38 | any other value |
| `gwas_locus_distance` > 0 | zero or negative |
| If `annotate=true`: `biofilter_loki` + `biofilter_script` set and exist | missing or not found |

Genetic data path validation is already handled by `check_input_genetic_data_parameters()`
and is not duplicated.

### R script level (both `identify_gwas_loci.R` and `make_saige_gwas_plots.R`)

| Check | Action |
|---|---|
| Input file (`--sumstats`) exists and is non-empty | `stop()` with clear message |
| `--phenoTable` exists (plotting script only) | `stop()` |
| `--build` is 37 or 38 | `stop()` |
| All specified column name args are present in the data | `stop()` naming the missing columns |
| `--p_threshold` / `--locus_distance` are valid numerics | `stop()` |
| Zero loci returned (loci script) | Write empty CSV, exit 0, emit `cli_alert_warning` |
| `--loci_csv` is `NO_FILE` or missing (plotting script) | Skip annotation silently, produce plot without labels |

---

## Implementation Order

1. **`scripts/identify_gwas_loci.R`** — write and test locally with example data
2. **`scripts/make_saige_gwas_plots.R`** — write and test locally
3. **`processes/saige_helpers.nf`** — add `validate_gwas_params()` and update
   `get_script_file_names()`
4. **`processes/saige_postprocessing.nf`** — add `identify_gwas_loci`,
   `collect_gwas_loci`; remove `make_summary_suggestive_gwas`; rename
   `make_summary_table_with_annot` output
5. **`processes/saige_visualization.nf`** — add `make_saige_gwas_plots_R`
6. **`workflows/saige_gwas.nf`** — update includes, log block, and workflow body
7. **`READMEs/SAIGE_GWAS_docs.md`** — update parameters table and outputs section
8. **Stub run** — verify with `-stub` flag against existing test config
9. **Test run** — verify with `test_data/GWAS/` data

---

## Out of Scope (this iteration)

- ExWAS, PheWAS, Gene PheWAS workflows
- The annotate plotting path (`make_gwas_plots_with_annot`) — retained in file, not
  wired
- Container rebuild (runtime `remotes::install_github()` is accepted for now)
- LD-based clumping (locus definition is distance-only; no reference panel required)
- Version pinning of `gwasRtools` / `levinmisc` (always-latest-`main`)

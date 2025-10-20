
Documentation for SAIGE ExWAS
=============================

# Module Overview


SAIGE ExWAS is a pipeline for doing whole-exome association study of rare variants and gene burdens with traits using SAIGE software. Please see 
- [Tool Paper Link for Reference](https://www.nature.com/articles/s41588-022-01178-w)
- [Tool Documentation Link for Reference](https://saigegit.github.io/SAIGE-doc/)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_exwas.config)
- [Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)

## Software Requirements


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build saige.sif docker://pennbiobank/saige:latest`

* Docker Command: `docker pull pennbiobank/saige:latest`

* Pull from Google Container Registry: `docker pull gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest`

* Run Command: `nextflow run /path/to/toolkit/module/workflows/saige_exwas.nf`

* Common `nextflow run` flags:

    * `-resume` flag picks up workflow where it left off

    * `-stub` performs a dry run, checks channels without executing code

    * `-profile` selects the compute profiles in nextflow.config

    * `-profile standard` uses the Docker image to execute processes

    * `-profile cluster` uses the Singularity container and submits processes to a queue

    * `-profile all_of_us` uses the Docker image on All of Us Workbench

* More info: [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Detailed Pipeline Steps

## Part I: Setup


1. Start your own tools directory and go there. You may do this in your project analysis directory, but it often makes sense to clone into a general `tools` location

```sh
# Make a directory to clone the pipeline into
TOOLS_DIR="/path/to/tools/directory"
mkdir $TOOLS_DIR
cd $TOOLS_DIR
```

2. Download the source code by cloning from git

```sh
git clone None
cd $TOOLS_DIR/pmbb-nf-toolkit-saige-family
```

3. Build the singularity image
    - you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
    - this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly.


```sh
cd $TOOLS_DIR/pmbb-nf-toolkit-saige-family
singularity build saige.sif docker://pennbiobank/saige:latest
```
## Part II: Configure your run


1. Make a separate analysis/run/working directory.
    - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses.
    - ❗This step is optional, but We Highly recommend making a `tools` directory separate from your `run` directory. We recommend storing the `nextflow.conf` in here as it shouldn't change between runs.


```sh
WDIR="/path/to/analysis/run1"
mkdir -p $WDIR
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/blob/main/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile $MY_PROFILE`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and passing it to the pipeline with `-c /path/to/nextflow.config`.


3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
    - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters.
    - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs).
    - you can compartamentalize your config file as much as you like by passing
    - There are 2 ways to specify the config file during a run:

        - with the `-c` option on the command line: `nextflow run -c SAIGE_FAMILY/configs/saige_exwas.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig SAIGE_FAMILY/configs/saige_exwas.config`

## Part III: Run your analysis


❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-saige-family/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_exwas.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_exwas.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_exwas.config \
   -resume
```
# Pipeline Parameters

## Input Files for SAIGE_ExWAS


* SAIGE Sparse GRM Sample IDs

    * (optional) sample IDs for a sparse relatedness matrix

    * Type: List File

    * Format: txt

    * File Header:


    ```
    1a1
    1a2
    1a3
    1a4
    1a5
    1a6
    1a7
    1a8
    1a9
    1a10
    ```

* Cohort Membership

    * 0/1 table with cohorts as columns and participants as rows - 1 indicates that that row’s participant is a member of that column’s cohort

    * Type: Data Table

    * Format: csv

    * File Header:


    ```
    IID,POP1,POP2,POP3
    1a1,1,0,1
    1a2,1,0,0
    1a3,0,0,0
    1a4,1,0,0
    1a5,1,0,0
    1a6,1,1,0
    1a7,0,0,0
    1a8,1,0,1
    1a9,0,1,1
    
    ```

* SAIGE Step 1 Plink Files

    * a hard-call plink set to use for step 1 (usually also exome or genotype files)

    * Type: Plink Set

    * Format: plink binary

    * File Header:


    ```
    nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.{bed,bim,fam}
    ```

* Phenotypes and Covariates

    * table with participants as rows and all needed phenotypes and covariates as columns

    * Type: Data Table

    * Format: csv

    * File Header:


    ```
    IID,y_quantitative,y_binary,x1,x2,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
    1a1,2.0046544617651,0,1.51178116845085,1,0,0,0,0,0,0,0,0,1,0
    1a2,0.104213400269085,0,0.389843236411431,1,0,0,0,0,0,0,0,0,1,1
    1a3,-0.397498354133647,0,-0.621240580541804,1,0,0,0,0,0,0,0,0,0,1
    1a4,-0.333177899030597,0,-2.2146998871775,1,0,0,0,0,0,0,0,0,1,1
    1a5,1.21333962248852,0,1.12493091814311,1,0,0,0,0,0,0,0,0,1,0
    1a6,-0.275411643032321,0,-0.0449336090152309,1,0,0,0,0,0,0,0,0,1,0
    1a7,0.438532936074923,0,-0.0161902630989461,0,0,0,0,0,0,0,0,0,0,0
    1a8,0.0162938047248591,0,0.943836210685299,0,0,0,0,0,0,0,0,0,1,1
    1a9,0.147167262428064,0,0.821221195098089,1,0,0,0,0,0,0,0,0,1,0
    ```

* SAIGE Group Annotation Files

    * text files formatted like this example from the SAIGE github under extdata/input/group_new_chrposa1a2.txt. They should be chromosome-separated, ending with “1.txt”…”22.txt”. The pipeline will fill in “chr.txt” to the end of the prefix you provide

    * Type: Data Table

    * Format: saige group (txt)

    * File Header:


    ```
    ENSG00000000457 var     1_169853716_C_A 1_169853716_C_T 1_169853717_C_CAGTT
    ENSG00000000457 anno    other_missense  damaging_missense       damaging_missense
    ENSG00000000460 var     1_169795119_C_T 1_169795121_G_C 1_169795123_C_G
    ENSG00000000460 anno    other_missense  other_missense  other_missense
    ```

* Gene Location File

    * CSV file of 

    * Type: Data Table

    * Format: tsv

    * File Header:


    ```
    gene_id chromosome  seq_region_start    seq_region_end  gene_symbol
    GENE1   1   1   90  GS1
    GENE2   2   91  100 GS2
    ```

* SAIGE Step 2 Plink Files

    * This is a set of chromosome-separated hard-call Plink Files for step 2 of SAIGE. The prefix should be indicated such that the chromosome and bed/bim/fam can be appended for each individual file.

    * Type: Plink Set

    * Format: plink binary

    * File Header:


    ```
    genotype_100markers_2chr.chr1.{bed,bim,fam}
    ```

* SAIGE Sparse GRM

    * (optional) a sparse relatedness matrix for faster step 1 computation

    * Type: GRM

    * Format: R sparse matrix

    * File Header:


    ```
    sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
    ```

* Sex Specific Phenotype List

    * A newline-separated list of phenotypes that should be excluded from non-sex-stratified cohorts (e.g., include in AFR_F or AFR_M but exclude from AFR_ALL). Set to 

    * Type: List File

    * Format: txt
## Output Files for SAIGE_ExWAS


* Quantitative Phenotype Violin Plots

    * A violin plot for each quantitative phenotype. One file will be generated for each phenotype and all cohorts will be plotted on each. 

    * Type: Summary Plot

    * Format: png

        * Parallel By: Phenotype

* Singles QQ Plots

    * A QQ Plot of the Null Model vs Log10P results of the analysis for variants.  One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: QQ Plot

    * Format: png

        * Parallel By: Cohort, Phenotype

* Phenotype Summary Table

    * A csv file with summary statistics for the phenotypes. Statistics are computed across all cohorts, phenotypes. For binary phenotypes, it also computes: Total Number, number of cases, Number of Controls, and Prevalence, and summary stats. 

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    COHORT |PHENO     |N  |Controls|Cases|Prevalence         |mean             |std              |min  |25%  |50%  |75%    |max
    AMR_ALL|AAA       |474|464.0   |10.0 |0.02109704641350211|                 |                 |     |     |     |       |
    AMR_ALL|T2D       |546|425.0   |121.0|0.2216117216117216 |                 |                 |     |     |     |       |
    AMR_ALL|BMI_median|531|        |     |                   |29.95419020715631|7.10618710646257 |15.17|24.84|28.69|33.435 |64.465
    AMR_ALL|LDL_median|310|        |     |                   |93.09516129032258|34.95209123947716|9.0  |70.25|90.0 |112.375|291.0
    ```

* Regions Manhattan Plots

    * A dot plot (manhattan plot) of significant gene regions associated with a phenotype. One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: Manhattan Plot

    * Format: png

        * Parallel By: Cohort, Phenotype, Annot Group, MAF

* Regions QQ Plots

    * A QQ Plot of the Null Model vs Log10P results of the analysis for gene regions.  One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: QQ Plot

    * Format: png

        * Parallel By: Cohort, Phenotype, Annot Group, MAF

* Singles Summary Statistics

    * A gzipped, unfiltered TSV (tab-separated) file of the results for the variant (singles) analysis. One file will be created for each unique Cohort, Phenotype, and analysis (regular, cauchy, rare, ultra rare) combination.

    * Type: Summary Statistics

    * Format: tsv.gz

    * File Header:


    ```
    phenotype|chromosome|base_pair_location|variant_id        |other_allele|effect_allele|effect_allele_count|effect_allele_frequency|missing_rate|beta      |standard_error|t_statistic|variance|p_value   |p_value_na|is_spa_test|allele_freq_case|allele_freq
    T2Diab   |21        |41801254          |21_41801254_TCTG_T|TCTG        |T            |277                |0.0046611              |0.0         |-0.099231 |0.167775      |-3.52526   |35.5258 |0.5542179 |0.5542179 |False      |0.00426841      |0.00474126
    T2Diab   |21        |41801360          |21_41801360_C_T   |C           |T            |41                 |0.00068991             |0.0         |-0.864441 |0.633121      |-3.98228   |5.38237 |0.08606924|0.08606924|False      |0.000297796     |0.000769948
    T2Diab   |21        |41801603          |21_41801603_C_T   |C           |T            |24                 |0.00040385             |0.0         |0.322923  |0.570593      |0.991852   |3.07148 |0.5714322 |0.5714322 |False      |0.000496327     |0.000384974
    T2Diab   |21        |41801645          |21_41801645_G_A   |G           |A            |58                 |0.000975971            |0.0         |0.0167811 |0.35132       |0.135962   |8.10206 |0.9619027 |0.9619027 |False      |0.00109192      |0.000952304
    ```

        * Parallel By: Cohort, Phenotype

* Singles Manhattan Plots

    * A dot plot (manhattan plot) of significant variants associated with a phenotype. One plot will be created for each unique combination of phenotype, cohort, annotation group (pLof, etc.), and MAF threshold. 

    * Type: Manhattan Plot

    * Format: png

        * Parallel By: Cohort, Phenotype

* Singles Top Hits Table

    * A FILTERED top hits csv summary file of results including cohort, phenotype, gene, group annotation, p-values, and other counts. One single summary file will be aggregated from all the “top hits” in each “Singles (Variant) Summary Statistics” file. 

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    cohort,phenotype,chromosome,base_pair_location,variant_id,other_allele,effect_allele,effect_allele_count,effect_allele_frequency,missing_rate,beta,standard_error,t_statistic,variance,p_value,p_value_na,is_spa_test,allele_freq_case,allele_freq_ctrl,n_case,n_ctrl,n_case_hom,n_case_het,n_ctrl_hom,n_ctrl_het,n
    EUR_M,BMI_median,1,1203928,1_1203928_G_A,G,A,13,0.000425393,0.0,6.68476,1.49857,2.97666,0.445291,8.167552e-06,,,,,,,,,,,15280.0
    AFR_M,BMI_median,1,6073749,1_6073749_C_T,C,T,14,0.00182292,0.0,8.90012,1.93388,2.37977,0.267386,4.180559e-06,,,,,,,,,,,3840.0
    AFR_F,LDL_median,1,11134343,1_11134343_G_C,G,C,13,0.00127226,0.0,36.8257,8.24963,0.541106,0.0146937,8.047218e-06,,,,,,,,,,,5109.0
    AFR_ALL,LDL_median,1,11828619,1_11828619_C_T,C,T,12,0.000740466,0.0,45.5826,9.2016,0.538359,0.0118106,7.279194e-07,,,,,,,,,,,8103.0
    ```

* Regions Summary Statistics

    * A gzipped, unfiltered TSV (tab-separated) file of the results for the gene (regions) analysis if run. One file will be created for each unique Cohort, Phenotype, and analysis (regular, cauchy, rare, ultra rare) combination.

    * Type: Summary Statistics

    * Format: tsv.gz

    * File Header:


    ```
    phenotype|gene           |annot            |max_maf|p_value           |p_value_burden    |p_value_skat      |beta_burden        |se_burden         |mac   |mac_case|mac_control|rare_var_count|ultrarare_var_count
    T2Diab   |ENSG00000141956|pLoF             |0.0001 |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|pLoF             |0.001  |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|pLoF             |0.01   |0.0479451461682565|0.0479451461682565|0.0479451461682565|0.0588652858997042 |0.0297621953331829|12.0  |5.0     |7.0        |0.0           |9.0
    T2Diab   |ENSG00000141956|damaging_missense|0.0001 |0.464219450219203 |0.464219450219203 |0.464219450219203 |-0.0110759683619445|0.0151328276810456|52.0  |7.0     |45.0       |0.0           |41.0
    ```

        * Parallel By: Cohort, Phenotype

* Binary Phenotype Bar Plots

    * A bar plot for each binary phenotype. One file will be generated for each phenotype and all cohorts will be plotted. 

    * Type: Summary Plot

    * Format: png

        * Parallel By: Phenotype

* Regions Top Hits Table

    * A FILTERED top hits csv summary file of results including cohort, phenotype, gene, group annotation, p-values, and other counts. One single summary file will be aggregated from all the “top hits” in each “Regions Summary Statistics” file. 

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    cohort,phenotype,gene,annot,max_maf,p_value,p_value_burden,p_value_skat,beta_burden,se_burden,mac,mac_case,mac_control,rare_var_count,ultrarare_var_count
    EUR_F,BMI_median,ENSG00000003393,pLoF,0.0001,1.12067607254946e-06,1.12067607254946e-06,1.12067607254946e-06,0.456714431175153,0.0937971711429878,9.0,,,0.0,8.0
    EUR_F,BMI_median,ENSG00000003393,pLoF,0.001,1.12067607254946e-06,1.12067607254946e-06,1.12067607254946e-06,0.456714431175153,0.0937971711429878,9.0,,,0.0,8.0
    EUR_F,BMI_median,ENSG00000003393,pLoF,0.01,1.12067607254946e-06,1.12067607254946e-06,1.12067607254946e-06,0.456714431175153,0.0937971711429878,9.0,,,0.0,8.0
    EUR_M,LDL_median,ENSG00000006530,pLoF,0.01,8.20696106989754e-06,8.20696106989754e-06,8.20696106989754e-06,3.23068908720413,0.724416423633968,3.0,,,0.0,3.0
    ```
## Other Parameters for SAIGE_ExWAS

### Association Test Modeling


* `cont_covars ` (Type: List)

    * Continuous covariates list

* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges
### Post-Processing


* `gene_location_file ` (Type: File Path)

    * This file is used for getting gene-based coordinates for plotting .

    * Corresponding Input File: Gene Location File

        * CSV file of 

        * Type: Data Table

        * Format: tsv

        * File Header:


        ```
        gene_id chromosome  seq_region_start    seq_region_end  gene_symbol
        GENE1   1   1   90  GS1
        GENE2   2   91  100 GS2
        ```

* `region_col_names` (Type: Map (Dictionary))

    * Default SAIGE Region column names mapped to new ones
### Pre-Processing


* `id_col ` (Type: String)

    * ID column label

* `data_csv` (Type: File Path)

    * A csv table with all of the phenotypes and covariates to be tested

    * Corresponding Input File: Phenotypes and Covariates

        * table with participants as rows and all needed phenotypes and covariates as columns

        * Type: Data Table

        * Format: csv

        * File Header:


        ```
        IID,y_quantitative,y_binary,x1,x2,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
        1a1,2.0046544617651,0,1.51178116845085,1,0,0,0,0,0,0,0,0,1,0
        1a2,0.104213400269085,0,0.389843236411431,1,0,0,0,0,0,0,0,0,1,1
        1a3,-0.397498354133647,0,-0.621240580541804,1,0,0,0,0,0,0,0,0,0,1
        1a4,-0.333177899030597,0,-2.2146998871775,1,0,0,0,0,0,0,0,0,1,1
        1a5,1.21333962248852,0,1.12493091814311,1,0,0,0,0,0,0,0,0,1,0
        1a6,-0.275411643032321,0,-0.0449336090152309,1,0,0,0,0,0,0,0,0,1,0
        1a7,0.438532936074923,0,-0.0161902630989461,0,0,0,0,0,0,0,0,0,0,0
        1a8,0.0162938047248591,0,0.943836210685299,0,0,0,0,0,0,0,0,0,1,1
        1a9,0.147167262428064,0,0.821221195098089,1,0,0,0,0,0,0,0,0,1,0
        ```
### QC Options


* `min_maf` (Type: Float)

    * Minimum minor allele frequency for plink QC
### SAIGE Step 1


* `geno` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants, genotype rate  filters out all variants with missing call rates exceeding the provided value 

* `exome_plink_prefix` (Type: Plink Fileset Prefix)

    * Exome plink input files 

* `group_file_prefix` (Type: Chr File Prefix)

    * Has the variant positions for each gene as well as the variant annotation for each variant in the gene in SAIGE format 

    * Corresponding Input File: SAIGE Group Annotation Files

        * text files formatted like this example from the SAIGE github under extdata/input/group_new_chrposa1a2.txt. They should be chromosome-separated, ending with “1.txt”…”22.txt”. The pipeline will fill in “chr.txt” to the end of the prefix you provide

        * Type: Data Table

        * Format: saige group (txt)

        * File Header:


        ```
        ENSG00000000457 var     1_169853716_C_A 1_169853716_C_T 1_169853717_C_CAGTT
        ENSG00000000457 anno    other_missense  damaging_missense       damaging_missense
        ENSG00000000460 var     1_169795119_C_T 1_169795121_G_C 1_169795123_C_G
        ENSG00000000460 anno    other_missense  other_missense  other_missense
        ```

* `hwe` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
### SAIGE Step 2


* `min_mac ` (Type: Float)

    * SAIGE-GENE Step 2 Parameters

* `firth_cutoff` (Type: Float)

    * P-value ()

* `use_firth` (Type: Bool (R: TRUE or FALSE))

    * True to use firth logistic regression
### Workflow


* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified

* `bin_pheno_list` (Type: List)

    * Binary phenotype list

* `quant_pheno_list` (Type: List)

    * Quantitative phenotype list
# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```

params {
    // DATA FILES
    // ----------
    // all default paths are for PMBB WES
    data_csv = "/project/pmbb_codeworks/datasets/CodeWorks_Test_Data/cleaned_test_pheno_covars.csv"
    // data_csv = "/project/pmbb_codeworks/datasets/pmbb_allwas_pheno/common_phecodes_covariate_ALL.csv"
    
    // cohort sets
    cohort_sets = "/project/pmbb_codeworks/datasets/PMBB_Extra/Sample_Lists/Exome_sample_table.csv"
    
    // this is for getting gene-based coordinates for plotting
    gene_location_file = "/project/pmbb_codeworks/datasets/ENSEMBL/homo_sapiens_111_b38.txt"

    // ID column label
    id_col = "PMBB_ID"
    
    // Full list of cohorts (usually ancestry-stratified and/or sex-stratified)
    cohort_list = [
        "PMBB_AMR_ALL", "PMBB_AMR_F","PMBB_AMR_M",
        "PMBB_AFR_ALL", "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_ALL", "PMBB_EUR_F", "PMBB_EUR_M"
        ]
        
    // smaller list of cohorts for testing
    // cohort_list = [
    //     "PMBB_AMR_ALL", "PMBB_AMR_F","PMBB_AMR_M",
    //     "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M"
    //    ]
    
    // subset of cohorts that are female- or male-only which should exclude sex-based covariates
    sex_strat_cohort_list = [
        "PMBB_AMR_F", "PMBB_AMR_M",
        "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_F", "PMBB_EUR_M",
        ]

    // smaller list of sex stratified for testing
    // sex_strat_cohort_list = [
    //     "PMBB_AMR_F","PMBB_AMR_M",
    //     "PMBB_EAS_F", "PMBB_EAS_M",
    //    ]
    
    // Full list of chromosomes
    chromosome_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
    
    // Small list of chormosomes for testing
    // chromosome_list = ["22"]

    // binary and quantitative phenotype [lists] or path to file of newline-separated lists
    // bin_pheno_list = "/project/pmbb_codeworks/datasets/pmbb_allwas_pheno/common_phecodes_list.txt"
    bin_pheno_list = ["T2D", "AAA"]
    quant_pheno_list = ["BMI_median", "LDL_median"]
    // sex_specific pheno file - these will be skipped for _ALL cohorts
    sex_specific_pheno_file = "/project/pmbb_codeworks/datasets/pmbb_allwas_pheno/phecode_Sex_specific.txt"

    // categorical and continuous covariates
    cat_covars = ["SEX"]
    cont_covars = ["DATA_FREEZE_AGE", "Exome_PC1", "Exome_PC2", "Exome_PC3", "Exome_PC4"]
    sex_strat_cat_covars = []
    sex_strat_cont_covars = cont_covars

    // NextFlow, Docker, and Singularity OPTIONS
    // ------------------------------------------
    // default assumes use of the docker container
    my_python = "/opt/conda/bin/python"

    // default paths assume use of the docker container
    step1_script = "/usr/local/bin/step1_fitNULLGLMM.R"
    step2_script = "/usr/local/bin/step2_SPAtests.R"

    // gpu paramater either ON or OFF, need to set config to -c nextflow_gpu.config
    GPU = 'OFF'
    
    // Minimum numbers for filtering cohort-phenotype combinations
    min_bin_cases = 100
    min_quant_n = 1000

    // Config parameters for using precomputed sparse GRM:
    // use_sparse_GRM = true
    // step 1 path should be the small subset of markers used to fit the GRM
    // step1_plink_prefix = "/project/pmbb_codeworks/datasets/PMBB_Extra/SAIGE_Step0_Exome/input/PMBB_exome_random_autosomal_markers"
    // step1_sparse_grm = "/project/pmbb_codeworks/datasets/PMBB_Extra/SAIGE_Step0_Exome/output/PMBB_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
    // step1_sparse_grm_samples = "/project/pmbb_codeworks/datasets/PMBB_Extra/SAIGE_Step0_Exome/output/PMBB_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

    // Config parameters for using real-time FULL GRM:
    use_sparse_GRM = false

    // Genetic Data Inputs:
    // Without a GRM, the exome plink set is used for step 1 because it needs rare variants
    exome_plink_prefix = "/project/pmbb_codeworks/datasets/PMBB-2.0_exome_GL_norm/plink/PMBB-Release-2020-2.0_genetic_exome_GL_norm"
    group_file_prefix = "/project/pmbb_codeworks/datasets/New_VEP_Annotations_2.0/SAIGE_Sets/subset."
    
    // Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
    // Current defaults are recommended by GBMI analysis plan
    maf = 0.01
    geno = 0.01
    hwe = 1E-6
    
    // SAIGE-GENE Step 2 Parameters
    // Current defaults are recommended by BRAVA analysis plan
    min_maf = 0
    min_mac = 0.5
    grouptest_maf = "0.0001,0.001,0.01"
    grouptest_annotation = "pLoF,damaging_missense,other_missense,synonymous,pLoF;damaging_missense,pLoF;damaging_missense;other_missense;synonymous"
    use_firth = "TRUE"
    firth_cutoff = 0.1
    LOCO = "FALSE"

    // P-Value Threshold for Summarizing Results at the End
    p_cutoff_summarize = 0.00001

    // Dictionary (Map) with default SAIGE Region column names mapped to new ones
    regions_col_names = [
        Region: 'gene',
        Group: 'annot',
        max_MAF: 'max_maf',
        Pvalue: 'p_value',
        Pvalue_Burden: 'p_value_burden',
        BETA_Burden: 'beta_burden',
        SE_Burden: 'se_burden',
        Pvalue_SKAT: 'p_value_skat',
        MAC: 'mac',
        MAC_case: 'mac_case',
        MAC_control: 'mac_control',
        Number_rare: 'rare_var_count',
        Number_ultra_rare: 'ultrarare_var_count'
    ]

    // Dictionary (Map) with default SAIGE SingleAssoc column names mapped to new ones
    singles_col_names = [
        CHR: 'chromosome',
        POS: 'base_pair_location',
        MarkerID: 'variant_id',
        Allele1: 'other_allele',
        Allele2: 'effect_allele',
        AC_Allele2: 'effect_allele_count',
        AF_Allele2: 'effect_allele_frequency',
        MissingRate: 'missing_rate',
        BETA: 'beta',
        SE: 'standard_error',
        Tstat: 't_statistic',
        var: 'variance',
        'p.value': 'p_value',
        'p.value.NA': 'p_value_na',
        'Is.SPA': 'is_spa_test',
        AF_case: 'allele_freq_case',
        AF_ctrl: 'allele_freq_ctrl',
        N_case: 'n_case',
        N_ctrl: 'n_ctrl',
        N_case_hom: 'n_case_hom',
        N_case_het: 'n_case_het',
        N_ctrl_hom: 'n_ctrl_hom',
        N_ctrl_het: 'n_ctrl_het',
        N: 'n'
    ]

}

```
## Current `nextflow.config` contents


```
//includeConfig 'configs/saige_exwas.config'
//includeConfig 'configs/saige_gene_phewas.config'
//includeConfig 'configs/saige_variant_phewas.config'

profiles {
    non_docker_dev {
        // run locally without docker
        process.executor = 'local'
    }

    standard {
        // run locally with docker
        process.executor = 'local'
        process.container = 'pennbiobank/saige:latest'
        docker.enabled = true
    }

    cluster {
        // run on LSF cluster
        process.executor = 'lsf'
        process.queue = 'epistasis_normal'
        executor {
            queueSize=500
        }
        process.memory = '15GB'
    	process.container = 'saige.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /project/,/static/'
    }

    all_of_us {
        // CHANGE EVERY TIME! These are specific for each user, see docs
        google.lifeSciences.serviceAccountEmail = 'pet-XXX-@terra.vpc-sc-XXXXXXXX.iam.gserviceaccount.com' // change to user-specific service email
        workDir='gs://fc-secure/path/to/workdir' // change to your user-specific working directory in your workspace bucket
        google.project = 'terra.vpc-sc-XXXXXXXX' // change to your user-specific project ID

        // These should not be changed unless you are an advanced user
        process.container = 'gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest' // GCR SAIGE docker container (static)

        // these are AoU, GCR parameters that should NOT be changed
        process.memory = '15GB' // minimum memory per process (static)
        process.executor = 'google-lifesciences' // AoU uses google-lifesciences (static)
        google.zone = "us-central1-a" // AoU uses central time zone (static)
        google.location = "us-central1"
        google.lifeSciences.debug = true 
        google.lifeSciences.network = "network"
        google.lifeSciences.subnetwork = "subnetwork"
        google.lifeSciences.usePrivateAddress = false
        google.lifeSciences.copyImage = "gcr.io/google.com/cloudsdktool/cloud-sdk:alpine"
        google.enableRequesterPaysBuckets = true
        // google.lifeSciences.bootDiskSize = "20.GB" // probably don't need this
        process{
                withName: 'call_saige_step1_bin' {
                container = '/tnnandi/saige-doe:2'
                }
                }
        process{
                withName: 'call_saige_step1_quant' {
                container = '/tnnandi/saige-doe:2'
                }
        }

    }   
    dnanexus{
        //This is a profile for running on a single machine on dnanexus
        process {
            container = 'pennbiobank/saige:latest'
            executor = 'local'
            memory = '15GB'
        }
        docker {
            enabled = true
        }
        process{
            withName: 'call_saige_step1_bin' {
            container = '/tnnandi/saige-doe:2'
            }
        }
        process{
            withName: 'call_saige_step1_quant' {
            container = '/tnnandi/saige-doe:2'
            }
        }
    }
}

params {
    skip_postprocessing_errors = false
}

process {
    withLabel: safe_to_skip {
        errorStrategy=params.skip_postprocessing_errors ? 'ignore' : 'terminate'
    }
}

```
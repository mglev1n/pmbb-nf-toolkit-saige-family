
Documentation for Single-Variant PheWAS
=======================================

# Module Overview


Single Variant PheWAS is a pipeline that uses SAIGE to perform single variant association tests across all PheCodes for a set of variants specified by the user.

Please see 
- [Tool Paper Link for Reference](https://www.nature.com/articles/s41588-018-0184-y)
- [Tool Documentation Link for Reference](https://saigegit.github.io/SAIGE-doc/)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_variant_phewas.config)
- [Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)

## Software Requirements


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build saige.sif docker://pennbiobank/saige:latest`

* Docker Command: `docker pull pennbiobank/saige:latest`

* Pull from Google Container Registry: `docker pull gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest`

* Run Command: `nextflow run /path/to/toolkit/module/workflows/saige_variant_phewas.nf`

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

        - with the `-c` option on the command line: `nextflow run -c SAIGE_FAMILY/configs/saige_variant_phewas.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig SAIGE_FAMILY/configs/saige_variant_phewas.config`

## Part III: Run your analysis


❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-saige-family/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_variant_phewas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_variant_phewas.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_variant_phewas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_variant_phewas.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_variant_phewas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_variant_phewas.config \
   -resume
```
# Pipeline Parameters

## Input Files for Single-Variant_PheWAS


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

* Quantitative Phenotype List

    * newline-delimited list of quantitative phenotypes to be included

    * Type: List File

    * Format: txt

    * File Header:


    ```
    y_quantitative
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

* Phenotype Descriptions File

    * File containing Phenotype notation,descriptions,and, categories (i.e Phe8.52,Intestinal infection due to C. difficile,Infection)

    * Type: Data Table

    * Format: csv

    * File Header:


    ```
    PHENO,DESCRIPTION,CATEGORY
    x1,description_x1,RED
    x2,description_x2,BLUE
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

* List of SNPs of Interest

    * Tab seperated file with one row for each variant of interest containing the chromosome,start,stop, and rsid in the same format as your genetic data

    * Type: Data Table

    * Format: tsv

    * File Header:


    ```
    1	0	0	rs1
    1	0	1	rs1
    1	0	2	rs2
    1	0	3	rs3
    1	0	5	rs4
    1	0	5	rs5
    1	0	6	rs6
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

* Binary Phenotype List

    * newline-delimited list of binary phenotypes to be tested

    * Type: List File

    * Format: txt

    * File Header:


    ```
    y_binary
    ```

* Sex Specific Phenotype List

    * A newline-separated list of phenotypes that should be excluded from non-sex-stratified cohorts (e.g., include in AFR_F or AFR_M but exclude from AFR_ALL). Set to 

    * Type: List File

    * Format: txt
## Output Files for Single-Variant_PheWAS


* Singles Summary Statistics

    * A gzipped tsv summary file with  phenotype, chromosome, base pair location, variant ID, Other allele, effect allele, effect,beta, statistics, pvalues, and other sumstats. One file will be generated for each combination of cohort, chromosome, and analysis (regular, cauchy, rare, ultra rare). One file will be generated for each cohort, chormosome, and analysis (regular, cauchy, rare, ultra rare) combination. 

    * Type: Summary Statistics

    * Format: tsv.gz

    * File Header:


    ```
    phenotype       chromosome      base_pair_location      variant_id      other_allele    effect_allele   effect_allele_count     effect_allele_frequency missing_rate beta    standard_error  t_statistic     variance        p_value p_value_na      is_spa_test     allele_freq_case        allele_freq_ctrl        n_case       n_ctrl  n_case_hom      n_case_het      n_ctrl_hom      n_ctrl_het
    Phe159.2        4       39199561        4_39199561_G_A  G       A       39      0.000652218     0.0     -1.00478        3.11147 -0.103787       0.103293    0.7467486        0.7467486       False   0.0     0.000653858     75      29823   0       0       0       39
    Phe159.2        4       39214639        4_39214639_A_G  A       G       14      0.000234129     0.0     -1.00338        5.62099 -0.0317569      0.0316501   0.8583262        0.8583262       False   0.0     0.000234718     75      29823   0       0       0       14
    Phe159.2        4       39215918        4_39215918_C_T  C       T       41      0.000685665     0.0     -1.00617        3.49319 -0.0824573      0.0819513   0.7733171        0.7733171       False   0.0     0.000687389     75      29823   0       0       0       41
    ```

        * Parallel By: Cohort, Chromosome

* Singles Summary Statistics Filtered

    * A csv FILTERED summary file with cohort, phenotype, chromosome, base pair location, variant ID, Other allele, effect allele, effect,beta, statistics, pvalues, and other sumstats. One file will be generated for each combination of cohort, chromosome, and analysis (regular, cauchy, rare, ultra rare)

    * Type: Summary Statistics

    * Format: csv

    * File Header:


    ```
    cohort,phenotype,chromosome,base_pair_location,variant_id,other_allele,effect_allele,effect_allele_count,effect_allele_frequency,missing_rate,beta,standard_error,t_statistic,variance,p_value,p_value_na,is_spa_test,allele_freq_case,allele_freq_ctrl,n_case,n_ctrl,n_case_hom,n_case_het,n_ctrl_hom,n_ctrl_het
    EUR_ALL,Phe783.0,4,39274909,4_39274909_C_T,C,T,77,0.00131171,0.0,2.6324,0.606004,5.43601,0.54679,6.999681e-06,1.961538e-13,True,0.015,0.0012178,200,29151,0,6,0,71
    ```

        * Parallel By: Cohort, Chromosome

* Singles Top Hits Table

    * A csv summary file with cohort, phenotype, chromosome, base pair location, variant ID, Other allele, effect allele, effect,beta, statistics, pvalues, and other sumstats

    * Type: Summary Table

    * Format: csv

    * File Header:


    ```
    cohort,phenotype,chromosome,base_pair_location,variant_id,other_allele,effect_allele,effect_allele_count,effect_allele_frequency,missing_rate,beta,standard_error,t_statistic,variance,p_value,p_value_na,is_spa_test,allele_freq_case,allele_freq_ctrl,n_case,n_ctrl,n_case_hom,n_case_het,n_ctrl_hom,n_ctrl_het
    AFR_F,Phe654.0,2,27444495,2_27444495_G_A,G,A,11,0.00080292,0.0,0.481686,0.0759042,-0.430299,1.33634e-15,1.1e-10,0.0E-2147483648,True,0.0,0.000846414,352,6498,0,0,0,11
    AFR_F,Phe643.0,2,27444495,2_27444495_G_A,G,A,11,0.000805153,0.0,12.1317,2.41425,-0.274659,1.47634e-15,2.5e-07,0.0E-2147483648,True,0.0,0.000830189,206,6625,0,0,0,11
    EUR_F,Phe644.0,2,27445431,2_27445431_C_T,C,T,43,0.00161849,0.0,9.01836,1.51694,-0.399219,5.45544e-15,1.4e-09,0.0E-2147483648,True,0.0,0.00163623,144,13140,0,0,0,43
    ```
## Other Parameters for Single-Variant_PheWAS

### Association Test Modeling


* `cont_covars ` (Type: List)

    * Continuous covariates list

* `sex_strat_cont_covars` (Type: List)

    * Continuous covariates for sex stratified cohorts to ensure model converges

* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges
### Post-Processing


* `pheno_descriptions_file` (Type: File Path)

    * file path to phenotype descriptions used for plotting and summary statistics

    * Corresponding Input File: Phenotype Descriptions File

        * File containing Phenotype notation,descriptions,and, categories (i.e Phe8.52,Intestinal infection due to C. difficile,Infection)

        * Type: Data Table

        * Format: csv

        * File Header:


        ```
        PHENO,DESCRIPTION,CATEGORY
        x1,description_x1,RED
        x2,description_x2,BLUE
        ```
### Pre-Processing


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


* `step1_plink_prefix` (Type: Plink Fileset Prefix)

    * Step1 exome plink input fileset  - should be all chromosomes together

    * Corresponding Input File: SAIGE Step 1 Plink Files

        * a hard-call plink set to use for step 1 (usually also exome or genotype files)

        * Type: Plink Set

        * Format: plink binary

        * File Header:


        ```
        nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.{bed,bim,fam}
        ```
### SAIGE Step 2


* `min_mac ` (Type: Float)

    * SAIGE-GENE Step 2 Parameters

* `firth_cutoff` (Type: Float)

    * P-value ()

* `use_firth` (Type: Bool (R: TRUE or FALSE))

    * True to use firth logistic regression

* `snplist` (Type: File Path)

    * file path to list of variants to be analyzed

    * Corresponding Input File: List of SNPs of Interest

        * Tab seperated file with one row for each variant of interest containing the chromosome,start,stop, and rsid in the same format as your genetic data

        * Type: Data Table

        * Format: tsv

        * File Header:


        ```
        1	0	0	rs1
        1	0	1	rs1
        1	0	2	rs2
        1	0	3	rs3
        1	0	5	rs4
        1	0	5	rs5
        1	0	6	rs6
        ```
### Workflow


* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified

* `bin_pheno_list` (Type: List OR File Path)

    * file path to list of binary phenotypes

    * Corresponding Input File: Binary Phenotype List

        * newline-delimited list of binary phenotypes to be tested

        * Type: List File

        * Format: txt

        * File Header:


        ```
        y_binary
        ```
# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```
params {
    // default assumes use of the docker container
    my_python = '/opt/conda/bin/python'
    my_python = '/home/guarelin/mambaforge/envs/py311/bin/python'
    
    // gpu paramater either ON or OFF
    GPU = 'OFF'
    
    data_csv = "/project/pmbb_codeworks/datasets/CodeWorks_Test_Data/cleaned_phewas_pheno_covars.csv"
    cohort_sets = "/project/pmbb_codeworks/datasets/PMBB_Extra/Sample_Lists/Genotype_sample_table.csv"

    // binary and quantitative phenotype lists
    // bin_pheno_list_file = "/project/pmbb_codeworks/datasets/CodeWorks_Test_Data/phecode_list_with_prefix.txt"
    bin_pheno_list = "/project/pmbb_codeworks/projects/SAIGE_FAMILY_TESTING/SAIGE_Gene_PheWAS/test_20_phecodes.txt"
    quant_pheno_list = "/project/pmbb_codeworks/datasets/PMBB_Extra/PheCodes_2.3/lab_list.txt"
    pheno_descriptions_file = "/project/pmbb_codeworks/datasets/Ontology_Things/phecode_descriptions_categories.csv"
    sex_specific_pheno_file = "/project/pmbb_codeworks/datasets/pmbb_allwas_pheno/phecode_Sex_specific.txt"

    //setting file type for step 2 (PLINK/BGEN)
    // ftype = "PLINK"
    ftype = "BGEN"

    // default paths are for PMBB Geno data (PLINK)
    step1_plink_prefix  = "/static/PMBB/PMBB-Release-2020-2.0/Genotype/PMBB-Release-2020-2.0_genetic_genotype"

    // default paths are for PMBB Geno data (BGEN)
    step2_bgen_prefix  = "/project/ritchie_scratch2/PMBB/GSA_V2_45K_NoFilter/merged/bgen/test_subset/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_10Ksubset_"
    samplefile = "/project/ritchie02/projects/gwPheWAS/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_bgen.sample"
    step2_pgen_prefix = "/static/PMBB/PMBB-Release-2020-2.0/Imputed/pgen/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_" 

    // desired snps to extract for step 2 (tsv 4 columns, one snp per line: chr start stop rsid)
    snplist = "/project/pmbb_codeworks/projects/geno_pheno_workbench_dev/SAIGE_FAMILY_GPU/AGMO.tsv"
    // step2_input_prefix = "AGMO" //arbitrary
    step2_plink_prefix = "/project/pmbb_codeworks/projects/geno_pheno_workbench_dev/SAIGE_FAMILY_GPU/AGMO"

    group_file_prefix = "/project/pmbb_codeworks/datasets/BRAVA_Variant_Annots/SAIGE_Sets/subset."

    info = 0.4
    // categorical and continuous covariates
    cat_covars = ["SEX"]
    cont_covars = ["DATA_FREEZE_AGE", "Genotype_PC1","Genotype_PC2","Genotype_PC3",	"Genotype_PC4"]

    // Covariates
    sex_strat_cat_covars = []
    sex_strat_cont_covars = cont_covars

    // P-Value Threshold for Summarizing Results at the End
    p_cutoff_summarize = 0.001

    // ID column label
    id_col = "PMBB_ID"

    min_bin_cases = 1000
    min_quant_n = 5000

    // Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
    // Current defaults are recommended by GBMI analysis plan
    maf = 0.01
    geno = 0.01
    hwe = 1E-6

    // SAIGE-GENE Step 2 Parameters
    // Current defaults are recommended by BRAVA analysis plan
    min_maf = 0
    min_mac = 0.5
    grouptest_maf = "0.01"
    grouptest_annotation = "pLoF,damaging_missense,other_missense,synonymous,pLoF;damaging_missense,pLoF;damaging_missense;other_missense;synonymous"
    use_firth = "TRUE"
    firth_cutoff = 0.1
    LOCO = "FALSE"

    // this is for getting gene-based coordinates for plotting
    // also wrapped in the docker container
    gene_location_file = "/app/NCBI.gene.loc"

    // list of cohorts (usually ancestry-stratified)
    cohort_list = [
        "PMBB_AFR_ALL","PMBB_AFR_F","PMBB_AFR_M",
        "PMBB_AMR_ALL","PMBB_AMR_F","PMBB_AMR_M",
        "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_ALL", "PMBB_EUR_F", "PMBB_EUR_M",
        "PMBB_SAS_ALL", "PMBB_SAS_F", "PMBB_SAS_M",
        ]

    sex_strat_cohort_list = [
        "PMBB_AFR_F","PMBB_AFR_M",
        "PMBB_AMR_M","PMBB_AMR_F",
        "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_F", "PMBB_EUR_M",
        "PMBB_SAS_F", "PMBB_SAS_M"
        ]

    // list of chromosome 
    chromosome_list = ["7"]
   
    // default paths assume use of the docker container
    step1_script = "/usr/local/bin/step1_fitNULLGLMM.R"      
    step2_script = "/usr/local/bin/step2_SPAtests.R"
   
    use_sparse_GRM = false
    // step1_sparse_grm = "/project/pmbb_codeworks/datasets/PMBB_Extra/SAIGE_Step0_Exome/output/PMBB_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
    // step1_sparse_grm_samples = "/project/pmbb_codeworks/datasets/PMBB_Extra/SAIGE_Step0_Exome/output/PMBB_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
   
    // Dictionary (Map) with default SAIGE Region column names mapped to new ones
    region_col_names = [
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
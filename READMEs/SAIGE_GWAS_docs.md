
Documentation for SAIGE GWAS
============================

# Module Overview


SAIGE GWAS is a pipeline for performing genome wide association studies of variants using the R-based SAIGE software. This module has the option of using the biofilter database to provide nearest gene annotation. 

Please see 
- [Tool Paper Link for Reference](https://www.nature.com/articles/s41588-018-0184-y)
- [Tool Documentation Link for Reference](https://saigegit.github.io/SAIGE-doc/)
- [Example Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gwas.config)
- [Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)

## Software Requirements


* [Nextflow version 24.04.3](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity 3.8.3](https://sylabs.io/docs/) OR [Docker 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build saige.sif docker://pennbiobank/saige:latest`

* Docker Command: `docker pull pennbiobank/saige:latest`

* Pull from Google Container Registry: `docker pull gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest`

* Run Command: `nextflow run /path/to/toolkit/module/workflows/saige_gwas.nf`

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

        - with the `-c` option on the command line: `nextflow run -c SAIGE_FAMILY/configs/saige_gwas.config`
        - in the `nextflow.config`: at the top of the file add: `includeConfig SAIGE_FAMILY/configs/saige_gwas.config`

## Part III: Run your analysis


❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly.❗We also HIGHLY recommend doing a TEST run with the included test data in `$TOOLS_DIR/pmbb-nf-toolkit-saige-family/test_data`we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gwas.config \
   -stub

# run an exwas for real
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gwas.config

# resume an exwas run if it was interrupted or ran into an error
nextflow run $TOOLS_DIR/pmbb-nf-toolkit-saige-family/workflows/saige_gwas.nf \
   -profile cluster \
   -c /path/to/nextflow.config \
   -c SAIGE_FAMILY/configs/saige_gwas.config \
   -resume
```
# Pipeline Parameters

## Input Files for SAIGE_GWAS


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

* SAIGE Step 2 VCF Files

    * Type: VCF/VCF Index Set

    * Format: vcf vcf.idx

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

* SAIGE Step 2 BGEN Files

    * Type: BGEN/Sample Set

    * Format: bgen

* Sex Specific Phenotype List

    * A newline-separated list of phenotypes that should be excluded from non-sex-stratified cohorts (e.g., include in AFR_F or AFR_M but exclude from AFR_ALL). Set to 

    * Type: List File

    * Format: txt
## Output Files for SAIGE_GWAS


_No output files defined for this module._
## Other Parameters for SAIGE_GWAS

### Association Test Modeling


* `sex_strat_cont_covars` (Type: List)

    * Continuous covariates for sex stratified cohorts to ensure model converges

* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges
### Post-Processing


* `biofilter_close_dist` (Type: Float)

    * The distance in bp for something to be considered “close” vs “far” with respect to nearest gene annotation. Value is often 5E4
### Pre-Processing


* `id_col ` (Type: String)

    * ID column label

* `min_bin_cases` (Type: Float)

    * For case-control filtering, the minimum number of BINARY phenotype cases you want to keep. Phenotypes with less than this number per cohort will be dropped (Default: 50 if not specified). 

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


* `use_firth` (Type: Bool (R: TRUE or FALSE))

    * True to use firth logistic regression
### Workflow


* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified

* `bin_pheno_list` (Type: List)

    * Binary phenotype list

* `bin_pheno_list` (Type: List OR File Path)

    * file path to list of binary phenotypes
# Configuration and Advanced Workflow Files

## Example Config File Contents (From Path)


```
params {
    // default assumes use of the docker container
    my_python = "/opt/conda/bin/python"
    my_bgenix = "/opt/conda/bin/bgenix"

    //setting file type for step 2 (PLINK/BGEN)
    //ftype = "PLINK"
    ftype = "BGEN"
    GPU="OFF"
    annotate=true   
    
    use_sparse_GRM = false
    step1_script = "/usr/local/bin/step1_fitNULLGLMM.R"
    step2_script = "/usr/local/bin/step2_SPAtests.R"

    data_csv = "/path/to/data/common_ICD_covariate_ALL.csv"

    cohort_sets = "/path/to/data/Imputed_sample_table.csv"

    // default paths are for PMBB Geno data
    step1_plink_prefix  = "/path/to/data/pruned_data"
    step2_plink_prefix = "/path/to/data/pruned_data"
    
    // default paths for Imputed Geno data BGEN
    step2_bgen_prefix  = "/path/to/data/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_chr"
    bgen_samplefile = "/path/to/data/PMBB-Release-2020-2.0_genetic_imputed-topmed-r2_bgen.sample"
    
    
    // categorical and continuous covariates
    cat_covars = ["SEX"]
    cont_covars = ["DATA_FREEZE_AGE", "Genotype_PC1","Genotype_PC2","Genotype_PC3",
                   "Genotype_PC4", "Genotype_PC5","Genotype_PC6","Genotype_PC7",
                   "Genotype_PC8","Genotype_PC9","Genotype_PC10"]

    min_bin_cases = null
    min_quant_n = null

    sex_strat_cat_covars = []
    sex_strat_cont_covars = cont_covars

    // P-Value Threshold for Summarizing Results at the End
    p_cutoff_summarize = 0.00001

    // ID column label
    id_col = "PMBB_ID"

    // Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants
    // Current defaults are recommended by GBMI analysis plan
    maf = 0.01
    geno = 0.01
    hwe = 1E-6
    
    thin_count = ""
    host = "LPC"

   //Step 2 Parameters
    min_maf = 0
    min_mac = 40
    firth_cutoff = 0.1
    LOCO = "TRUE"
    is_imputed_data="TRUE" 
    minInfo=0.3 
    inverseNormalize="TRUE"

 // this is for getting gene-based coordinates for plotting
    // also wrapped in the docker container
    gene_location_file = "/app/NCBI.gene.loc"
         
    cohort_list = [
        "PMBB_AMR_ALL", "PMBB_AMR_F", "PMBB_AMR_M",
        "PMBB_AFR_ALL", "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_ALL", "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_ALL", "PMBB_EUR_F", "PMBB_EUR_M",
        "PMBB_SAS_ALL", "PMBB_SAS_F", "PMBB_SAS_M",
        ]

    // subset of cohorts that are female- or male-only which should exclude sex-based covariates
    sex_strat_cohort_list = [
        "PMBB_AMR_F", "PMBB_AMR_M",
        "PMBB_AFR_F", "PMBB_AFR_M",
        "PMBB_EAS_F", "PMBB_EAS_M",
        "PMBB_EUR_F", "PMBB_EUR_M",
        "PMBB_SAS_F", "PMBB_SAS_M"
        ]

    // binary and quantitative phenotype lists
    bin_pheno_list = "/path/to/data/common_ICD_list.txt"
    quant_pheno_list = []

    sex_specific_pheno_file = "/path/to/data/icd_Sex_specific.txt"
    
    gwas_col_names = [
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
        N_ctrl_het: 'n_ctrl_het'
    ]

    // list of chromosomes
     chromosome_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
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
        process.executor = awsbatch-or-lsf-or-slurm-etc
    }

    standard {
        // run locally with docker
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.container = 'pennbiobank/saige:latest'
        docker.enabled = true
    }

    cluster {
        // run on LSF cluster
        process.executor = awsbatch-or-lsf-or-slurm-etc
        process.queue = 'epistasis_normal'
        executor {
            queueSize=500
        }
        process.memory = '15GB'
    	process.container = 'saige.sif'
        singularity.enabled = true
        singularity.runOptions = '-B /root/,/directory/,/names/'
    }

    all_of_us {
        // CHANGE EVERY TIME! These are specific for each user, see docs
        google.lifeSciences.serviceAccountEmail = service@email.gservicaaccount.com
        workDir = /path/to/workdir/ // can be gs://
        google.project = terra project id

        // These should not be changed unless you are an advanced user
        process.container = 'gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest' // GCR SAIGE docker container (static)

        // these are AoU, GCR parameters that should NOT be changed
        process.memory = '15GB' // minimum memory per process (static)
        process.executor = awsbatch-or-lsf-or-slurm-etc
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
## Advanced Nextflow Users: Take/Emit Info

### Input Channel (take) Description


NONE
### Output Channel (emit) Description


Singles_merge_output- a Channel of two paths to the merged raw and filtered output from SAIGE Step 2
Pheno_table- a Channel for a file that contains sample informations on each phenotype
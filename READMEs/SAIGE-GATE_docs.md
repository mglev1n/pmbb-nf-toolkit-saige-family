
Documentation for SAIGE-GATE
============================

# Module Overview


[Example Module Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gate.config)

[Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)
## Cloning Github Repository


* Command: `git clone https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench.git`

* Navigate to relevant workflow directory...
## Software Requirements


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

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
# Input Files for SAIGE-GATE


* Quantitative Phenotype List

    * newline-delimited list of quantitative phenotypes to be included

    * Type: List File

    * Format: txt

    * File Header:


    ```
    y_quantitative
    ```

* Binary Phenotype List

    * newline-delimited list of binary phenotypes to be tested

    * Type: List File

    * Format: txt

    * File Header:


    ```
    y_binary
    ```

* SAIGE Step 2 BGEN Files

    * Type: BGEN/Sample Set

    * Format: bgen

* SAIGE Step 2 Plink Files

    * This is a set of chromosome-separated hard-call Plink Files for step 2 of SAIGE. The prefix should be indicated such that the chromosome and bed/bim/fam can be appended for each individual file.

    * Type: Plink Set

    * Format: plink binary

    * File Header:


    ```
    genotype_100markers_2chr.chr1.{bed,bim,fam}
    ```

* SAIGE Step 1 Plink Files

    * a hard-call plink set to use for step 1 (usually also exome or genotype files)

    * Type: Plink Set

    * Format: plink binary

    * File Header:


    ```
    nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.{bed,bim,fam}
    ```
# Output Files for SAIGE-GATE


* QQ Plot

    * Type: QQ Plot

        * Parallel By: Cohort, Phenotype

* Manhattan Plots

    * Type: Manhattan Plot

        * Parallel By: Cohort, Phenotype

* Pheno Summary Plots

        * Parallel By: Cohort

* Phenotype Summaries

    * Type: Data Table

        * Parallel By: Cohort

* Summary Suggestive Singles

    * Type: Summary Statistics

        * Parallel By: Cohort, Phenotype
# Parameters for SAIGE-GATE

## Association Test Modeling


* `cont_covars ` (Type: List)

    * Continuous covariates list
## Post-Processing


* `annotate` (Type: Bool (Java: true or false))

    * Whether or not to annotate results with the RSIDs and nearest genes for plotting and summary files.
## Pre-Processing


* `min_survival_cases` (Type: Integer)

    * The min number of cases necessary for time-to-event analysis phenotypes

* `my_bgenix` (Type: String)

    * Path to bgenix executable used for calling bgenix commands within the workflow processes

* `min_quant_n` (Type: Float)

    * For case-control filtering, the minimum number of QUANTITATIVE phenotypes you want to keep. Phenotypes with less than this number per cohort will be dropped (Default: 50 if not specified). 

* `min_bin_cases` (Type: Float)

    * For case-control filtering, the minimum number of BINARY phenotype cases you want to keep. Phenotypes with less than this number per cohort will be dropped (Default: 50 if not specified). 

* `sex_specific_pheno_file` (Type: File Path)

    * A newline-separated list of phenotypes that should only be included in sex-stratified cohorts (e.g., AFR_F but not AFR_ALL).  Can be safely left as null (defaults to an empty List)
## SAIGE Step 1


* `use_sparse_GRM` (Type: Bool (Java: true or false))

    * Whether to use sparse GRM (or full = false)

* `geno` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants, genotype rate  filters out all variants with missing call rates exceeding the provided value 

* `step1_script  ` (Type: File Path)

    * Fits the null logistic/linear mixed model using a full or a sparse genetic relationship matrix (GRM). The GRM estimate the genetic relationship between two individuals over a certain number of SNPs
## SAIGE Step 2


* `ftype` (Type: String)

    * “PLINK” or “BGEN” based on input type for steps 1 & 2 
## Workflow


* `gwas_col_names` (Type: Map (Dictionary))

    * Map of SAIGE Output Column Names to desired Summary Statistics column names

* `event_time_bin` (Type: Integer)

    * 
The width to group event occurrence could be chosen differently based on the nature of your time-to-event data. If the unit of event time is year, you could choose eventTimeBinSize=1; if the unit of event time is month, eventTimeBinSize=1/12 could be a preferred choice.

* `event_time_col` (Type: String)

    * string that defines which column in the phenotype file to use for event time

* `survival_pheno_list` (Type: List OR File Path)

    * List/file with Phenotypes for which to perform time-to-event analysis

* `GPU` (Type: String)

    * whether or not to use GPU processes for SAIGE Step 1

* `my_python` (Type: File Path)

    * Path to the python executable to be used for python scripts - often it comes from the docker/singularity container (/opt/conda/bin/python)

* `quant_pheno_list` (Type: List OR File Path)

    * file path to list of quantitiative phenotypes

    * Corresponding Input File: Quantitative Phenotype List

        * newline-delimited list of quantitative phenotypes to be included

        * Type: List File

        * Format: txt

        * File Header:


        ```
        y_quantitative
        ```

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

    // use only for time-to-event gwas 
    survival_pheno_list = []
    event_time_col = "" 
    event_time_bin =  // event time bin size 
    // (must be in same unit of measure of event time col (months,years,days))
    // (e.g. event_time_col = ???? )
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
## Current Dockerfile for Container/Image


```docker
# build: for python packages, plink, biofilter, and NEAT plots
# FROM continuumio/miniconda3:latest AS build
FROM mambaorg/micromamba:1.5.6 AS build
WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3
# activate micromamba env
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Switch to root user for installation since default micromamba is $MAMBA_USER
USER root
# Install essential build dependencies only
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        wget \
        unzip \
        ca-certificates \
        libtiff5-dev \
    && apt-get clean \
    # cleanup
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    # Configure wget to retry downloads
    && echo "retry_connrefused = on" >> /etc/wgetrc \
    && echo "tries = 5" >> /etc/wgetrc && \
    # Create work directory with proper permissions
    mkdir -p /app && chown -R $MAMBA_USER:$MAMBA_USER /app

# Copy environment.yml file with proper ownership
# COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Switch to mambauser for conda operations
USER $MAMBA_USER

# Install packages from environment.yml
# RUN micromamba install -y -n base -f /tmp/environment.yml && \
#     micromamba clean --all --yes

# Create conda environment with minimal dependencies
RUN micromamba install -y -n base -c conda-forge \
        python=3.12 \
        numpy \
        pandas \
        scipy \
        matplotlib \
        seaborn \
        adjustText \
        apsw \
        sqlite \
        dominate \
        conda-build \
    && micromamba clean --all --yes

# Install bgenix separately since she's fickle
RUN micromamba install -y -n base -c conda-forge \
        bgenix \
    && micromamba clean --all --yes \
    # Verify bgenix is installed
    && micromamba run -n base which bgenix

# Install libtiff separately
RUN micromamba install -y -n base -c conda-forge \
        libtiff \
    && micromamba clean --all --yes

# Set USER to ROOT
USER root

# Handle libtiff separately as it's problematic
RUN export MAMBA_ROOT_PREFIX=/opt/conda && \
    LIBTIFF_PATH=$(find $MAMBA_ROOT_PREFIX -name "libtiff.so.*" -type f | sort -V | tail -n 1) && \
    if [ -n "$LIBTIFF_PATH" ]; then \
        ln -sf $LIBTIFF_PATH $MAMBA_ROOT_PREFIX/lib/libtiff.so.5; \
    else \
        echo "WARNING: Could not find libtiff.so.* files to link" && exit 1; \
    fi


# PLINK: Download and install tools with retries in parallel and keep in separate layer
RUN wget --tries=5 --retry-connrefused https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240526.zip -O plink2.zip & \
    wget --tries=5 --retry-connrefused https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O plink.zip & \
    wait && \
    if [ ! -f plink2.zip ] || [ ! -f plink.zip ]; then \
        echo "ERROR: Failed to download PLINK files" && exit 1; \
    fi && \
    unzip -o plink2.zip && \
    unzip -o plink.zip && \
    if [ ! -f plink ] || [ ! -f plink2 ]; then \
        echo "ERROR: PLINK extraction failed" && exit 1; \
    fi && \
    chmod +x plink plink2 && \
    rm -rf plink2.zip plink.zip

# Install NEAT-plots (shallow clone to speed up)
RUN git clone --depth 1 --single-branch https://github.com/PMBB-Informatics-and-Genomics/NEAT-Plots.git \
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    && rm -rf NEAT-Plots

# Install biofilter with retries
RUN wget --tries=5 --retry-connrefused \
    https://github.com/RitchieLab/biofilter/releases/download/Biofilter-${BIOFILTER_VERSION}/biofilter-${BIOFILTER_VERSION}.tar.gz \
    -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app/ \
    && export MAMBA_ROOT_PREFIX=/opt/conda \
    && $MAMBA_ROOT_PREFIX/bin/python setup.py install \
    && chmod a+rx /app/biofilter.py \
    && rm -rf biofilter.tar.gz

# Create a package with necessary Python files - dynamically get python version
RUN mkdir -p /python-dist/bin && \
    mkdir -p /python-dist/lib && \
    # copy all of conda environment
    cp -a /opt/conda/. /python-dist/ && \
    # Clean up unnecessary files
    rm -rf /python-dist/pkgs/* && \
    rm -rf /python-dist/conda-meta/* && \
    # Determine Python version dynamically
    find /python-dist -name "__pycache__" -type d -exec rm -rf {} +  2>/dev/null || true

# DEV STAGE FOR R AND SAIGE
# ------------------------------------------------------------------------------
# dev: for SAIGE and SAIGE-dependent packages
# dev: build minimal R with SAIGE
# dev: for SAIGE and SAIGE-dependent packages
FROM rocker/tidyverse:4.2.3 AS dev
# FROM rocker/r-ver:4.2.3 AS dev
WORKDIR /app

# Install system dependencies with retry mechanism
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libopenblas-base \
        python3-pip \
        r-cran-devtools \
        git \
        ca-certificates \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libgfortran5 \
        libgomp1 \
        liblapack-dev \
        libopenblas-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --retries 5 cget

# Install R primary dependencies
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'RcppParallel', 'RcppEigen'), repos='https://cloud.r-project.org')" && \
    R -e "install.packages(c('Matrix','methods', 'BH', 'dplyr', 'dbplyr', 'tibble', 'sparsesvd', 'stringr', 'lintools'), repos='https://cloud.r-project.org')" && \
    R -e "install.packages(c('ellipsis', 'vctrs','RSQLite'), repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/BH/BH_1.78.0-0.tar.gz', repos=NULL, type='source')"
# Install R secondary dependencies
RUN R -e "install.packages(c('optparse', 'SPAtest', 'SKAT', 'RhpcBLASctl'), repos='https://cloud.r-project.org')"
# Install development R dependencies 
RUN R -e "install.packages(c('devtools', 'R.utils', 'roxygen2', 'rversions'), repos='https://cloud.r-project.org')"
# Install GitHub R packages separately
RUN R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), devtools::install_github('leeshawn/MetaSKAT'))" 
RUN R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), devtools::install_github('cysouw/qlcMatrix'))"

# Install SAIGE with specific version and error handling
RUN git clone --depth 1 https://github.com/saigegit/SAIGE.git && \
    cd SAIGE && \
    R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), devtools::install_deps('.'))" && \
    R CMD INSTALL . && \
    mkdir -p /usr/local/bin && \
    cp extdata/step1_fitNULLGLMM.R \
       extdata/step2_SPAtests.R \
       extdata/step3_LDmat.R \
       extdata/createSparseGRM.R \
       /usr/local/bin/ && \
    chmod +x /usr/local/bin/*.R && \
    cd .. && \
    rm -rf SAIGE

# MAIN STAGE FINAL
# ------------------------------------------------------------------------------
# main: final image with only necessary packages and scripts
FROM ubuntu:22.04 AS main
# FROM rocker/tidyverse:4.2.3 AS main
# FROM rocker/r-ver:4.2.3 AS main
WORKDIR /app

# Set non-interactive installation
ENV R_VERSION=4.2.3 \
    DEBIAN_FRONTEND=noninteractive

# Install additional system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        libtiff5-dev \
        libgomp1 \
        libopenblas-base \
        wget \
        unzip \
        locales \
        build-essential \
        # add R build dependencies
        gfortran \
        libbz2-dev \
        libpcre2-dev \
        libreadline-dev \
        libx11-dev \
        libxt-dev \
        xorg-dev \
        liblzma-dev \
        zlib1g-dev \
        libcurl4 \
        # Add R package dependencies
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libcairo2-dev \
        libpng-dev \
        libjpeg-dev \
        libicu-dev \
        liblapack-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && locale-gen en_US.UTF-8

# Download and install R
RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
    && tar -xf R-${R_VERSION}.tar.gz \
    && cd R-${R_VERSION} \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

# Copy the entire Python distribution from build stage
COPY --from=build /python-dist /opt/conda

# copy needed binaries from the build stage
COPY --from=build /app/biofilter.py /app/
COPY --from=build /app/manhattan-plot/ /app/manhattan-plot/
COPY --from=build /app/plink /app/plink2 /usr/bin/

# Copy R shared libraries explicitly
COPY --from=dev /usr/lib/R/lib/ /usr/lib/R/lib/
COPY --from=dev /usr/local/lib/R/ /usr/local/lib/R/

# Copy SAIGE and ALL installed libararies from dev stage
# COPY --from=dev /usr/local/lib/R/site-library/ /usr/local/lib/R/site-library/

# # Copy SAIGE and only required R libraries rather than all of them
# COPY --from=dev /usr/local/lib/R/site-library/SAIGE/ /usr/local/lib/R/site-library/SAIGE/
# COPY --from=dev /usr/local/lib/R/site-library/SPAtest/ /usr/local/lib/R/site-library/SPAtest/
# COPY --from=dev /usr/local/lib/R/site-library/SKAT/ /usr/local/lib/R/site-library/SKAT/
# COPY --from=dev /usr/local/lib/R/site-library/MetaSKAT/ /usr/local/lib/R/site-library/MetaSKAT/
# COPY --from=dev /usr/local/lib/R/site-library/Rcpp*/ /usr/local/lib/R/site-library/
# COPY --from=dev /usr/local/lib/R/site-library/Matrix/ /usr/local/lib/R/site-library/Matrix/
# COPY --from=dev /usr/local/lib/R/site-library/BH/ /usr/local/lib/R/site-library/BH/
# COPY --from=dev /usr/local/lib/R/site-library/RhpcBLASctl/ /usr/local/lib/R/site-library/RhpcBLASctl/
# COPY --from=dev /usr/local/lib/R/site-library/data.table/ /usr/local/lib/R/site-library/data.table/
# COPY --from=dev /usr/local/lib/R/site-library/R.utils/ /usr/local/lib/R/site-library/R.utils/
# COPY --from=dev /usr/local/lib/R/site-library/optparse/ /usr/local/lib/R/site-library/optparse/

# COPY SAIGE scripts 
COPY --from=dev /usr/local/bin/step*.R /usr/local/bin/
COPY --from=dev /usr/local/bin/createSparseGRM.R /usr/local/bin/

# Set up environment variables (THESE AREN'T TRANSLATED TO SINGULARITY!!!!)
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8 \
    EDITOR=nano \
    PATH="/usr/bin:/usr/local/bin:/opt/conda/bin:${PATH}" \
    PYTHONPATH="/app/manhattan-plot:${PYTHONPATH:-}" \
    LD_LIBRARY_PATH="/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/lib64:/usr/lib:${LD_LIBRARY_PATH:-}" \
    OMP_NUM_THREADS=1 

USER root

# Create symbolic links for Python to ensure it's in the PATH
RUN ln -sf /opt/conda/bin/python /usr/bin/python && \
    ln -sf /opt/conda/bin/python3 /usr/bin/python3 && \
    ln -sf /opt/conda/bin/bgenix /usr/bin/bgenix && \
    # Make sure ldconfig recognizes our libraries
    ldconfig

# verfy the installation
RUN echo "Checking all required tools:" && \
    which python && python --version || echo "python verification failed" && \
    which Rscript && Rscript --version || echo "Rscript verification failed" && \
    which plink && plink --version || echo "plink verification failed" && \
    which plink2 && plink2 --version || echo "plink2 verification failed" && \
    which bgenix && bgenix -help || echo "bgenix verification failed" && \
    R --vanilla -e "cat('R session information:'); sessionInfo()" || echo "R testing failed" && \
    Rscript -e "library(SAIGE); cat('SAIGE package successfully loaded\n')" || echo "SAIGE package loading failed" && \
    createSparseGRM.R  --help || echo "createSparseGRM.R verification failed" && \
    step1_fitNULLGLMM.R --help || echo "step1_fitNULLGLMM.R verification failed" && \
    step2_SPAtests.R --help || echo "step2_SPAtests.R verification failed" && \
    step3_LDmat.R --help || echo "step3_LDmat.R verification failed"

# Build Command
# docker buildx build -f Dockerfile_ubuntu --platform  linux/amd64 --push -t pennbiobank/saige:latest-ubuntu .
# singularity build saige_ubuntu.sif docker://pennbiobank/saige:latest-ubuntu
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
# Detailed Pipeline Steps


from pathlib import Path

detailed_steps_file = Path("Markdowns/Pipeline_Detailed_Steps.md")

# Write the detailed steps content to a separate file
detailed_steps_file

# Detailed Steps for Runnning One of our Pipelines

Note: test data were obtained from the [SAIGE github repo](https://github.com/saigegit/SAIGE).

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
git clone https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family.git
cd ${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/
```

3. Build the `saige.sif` singularity image
- you may call the image whatever you like, and store it wherever you like. Just make sure you specify the name in `nextflow.conf`
- this does NOT have to be done for every saige-based analysis, but it is good practice to re-build every so often as we update regularly. 

```sh
cd ${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/
singularity build saige.sif docker://pennbiobank/saige:latest
```

## Part II: Configure your run

1. Make a separate analysis/run/working directory.
   - The quickest way to get started, is to run the analysis in the folder the pipeline is run. However, subsequent analyses will over-write results from previous analyses. 
   - ❗This step is optional, but We Highly recommend making a  `tools` directory separate from your `run` directory. The only items that need to be in the run directory are the `nextflow.conf` file and the `${workflow}.conf` file.

```sh
WDIR="/path/to/analysis/run1"
mkdir -p 
cd $WDIR
```

2. Fill out the `nextflow.config` file for your system.
    - See [Nextflow configuration documentation](https://www.nextflow.io/docs/latest/config.html) for information on how to configure this file. An example can be found on our GitHub: [Nextflow Config](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/Example_Configs/nextflow.config).
    - ❗IMPORTANTLY, you must configure a user-defined profile for your run environments (local, docker, saige, cluster, etc.). If multiple profiles are specified, run with a specific profile using `nextflow run -profile ${MY_PROFILE}`.
    - For singularity, The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to` with the location where you built the image above). See [Nextflow Executor Information](https://www.nextflow.io/docs/latest/executor.html) for more details.
    - ⚠️As this file remains mostly unchanged for your system, We recommend storing this file in the `tools/pipeline` directory and symlinking it to your run directory.

3. Create a pipeline-specific `.config` file specifying your run parameters and input files. See Below for workflow-specific parameters and what they mean.
   - Everything in here can be configured in `nextflow.config`, however we find it easier to separate the system-level profiles from the individual run parameters. 
   - Examples can be found in our Pipeline-Specific [Example Config Files](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/Example_Configs/).
   - you can compartamentalize your config file as much as you like by passing 
   - There are 2 ways to specify the config file during a run:
      - with the `-c` option on the command line: `nextflow run -c /path/to/workflow.conf`
      - in the `nextflow.conf`: at the top of the file add: `includeConfig '/path/to/workflow.conf'` 

## Part III: Run your analysis

- ❗We HIGHLY recommend doing a STUB run to test the analysis using the `-stub` flag. This is a dry run to make sure your environment, parameters, and input_files are specified and formatted correctly. 
- ❗We HIGHLY recommend doing a test run with the included test data in `${TOOLS_DIR}/pmbb-nf-toolkit-saige-family/test_data`
- in the `test_data/` directory for each pipeline, we have several pre-configured analyses runs with input data and fully-specified config files.

```sh
# run an exwas stub
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf -stub
# run an exwas for real
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf
# resume an exwas run if it was interrupted or ran into an error
nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -c /path/to/run1/exwas.conf -resume
```

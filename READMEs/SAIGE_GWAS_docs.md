
Documentation for SAIGE GWAS
============================

# Module Overview


SAIGE GWAS is a pipeline for performing genome wide association studies of variants using the R-based SAIGE software. This module has the option of using the biofilter database to provide nearest gene annotation.

[Paper Link for Reference](https://www.nature.com/articles/s41588-018-0184-y)

[Tool Documentation Link](https://saigegit.github.io/SAIGE-doc/)

[Example Module Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gwas.config)

[Example nextflow.config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/nextflow.config)
## Cloning Github Repository:


* Command: `git clone https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench.git`

* Navigate to relevant workflow directory run commands (our pipelines assume all of the nextflow files/scripts are in the current working directory)
## Software Requirements:


* [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)

* [Singularity version 3.8.3](https://sylabs.io/docs/) OR [Docker version 4.30.0](https://docs.docker.com/)
## Commands for Running the Workflow


* Singularity Command: `singularity build saige.sif docker://pennbiobank/saige:latest`

* Docker Command: `docker pull pennbiobank/saige:latest`

* Command to Pull from Google Container Registry: `docker pull gcr.io/verma-pmbb-codeworks-psom-bf87/saige:latest`

* Run Command: `nextflow run /path/to/toolkit/module/workflows/saige_gwas.nf`

* Common `nextflow run` flags:

    * `-resume` flag picks up the workflow where it left off, otherwise, the workflow will rerun from the beginning

    * `-stub` performs a sort of dry run of the whole workflow, checks channels without executing any code

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile` selects the compute profiles we set up in nextflow.config (see nextflow.config file below)

    * `-profile standard` uses the docker image to executes the processes

    * `-profile cluster` uses the singularity container and submits processes to a queue- optimal for HPC or LPC computing systems

    * `-profile all_of_us` uses the docker image to execute pipelines on the All of Us Researcher Workbench

* for more information visit the [Nextflow documentation](https://www.nextflow.io/docs/latest/cli.html)
# Configuration Parameters and Input File Descriptions

## Workflow


* `bin_pheno_list` (Type: List OR File Path)

    * file path to list of binary phenotypes

* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified

* `bin_pheno_list` (Type: List)

    * Binary phenotype list
## Pre-Processing


* `min_quant_n` (Type: Float)

    * For case-control filtering, the minimum number of QUANTITATIVE phenotypes you want to keep. Phenotypes with less than this number per cohort will be dropped (Default: 50 if not specified). 

* `min_bin_cases` (Type: Float)

    * For case-control filtering, the minimum number of BINARY phenotype cases you want to keep. Phenotypes with less than this number per cohort will be dropped (Default: 50 if not specified). 

* `cohort_sets` (Type: File Path)

    * A binary csv table in which the columns are the cohorts and the rows are the individuals. A 1 means that individual is a member of the column’s cohort, and a 0 means they aren’t.

    * Corresponding Input File: Cohort Membership

        * 0/1 table with cohorts as columns and participants as rows - 1 indicates that that row’s participant is a member of that column’s cohort

        * Type: Data Table

        * Format: csv

        * Input File Header:





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

* `data_csv` (Type: File Path)

    * A csv table with all of the phenotypes and covariates to be tested

    * Corresponding Input File: Phenotypes and Covariates

        * table with participants as rows and all needed phenotypes and covariates as columns

        * Type: Data Table

        * Format: csv

        * Input File Header:





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
## Association Test Modeling


* `sex_strat_cont_covars` (Type: List)

    * Continuous covariates for sex stratified cohorts to ensure model converges

* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges
## SAIGE Step 1


* `step1_sparse_grm_samples` (Type: File Path)

    * List of IDs to use in the sparse GRM 

    * Corresponding Input File: SAIGE Sparse GRM Sample IDs

        * (optional) sample IDs for a sparse relatedness matrix

        * Type: List File

        * Format: txt

        * Input File Header:





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

* `step1_sparse_grm` (Type: File Path)

    * File Path to precomputed Sparse GRM

    * Corresponding Input File: SAIGE Sparse GRM

        * (optional) a sparse relatedness matrix for faster step 1 computation

        * Type: GRM

        * Format: R sparse matrix

        * Input File Header:





        ```
        sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
        ```
## SAIGE Step 2


* `ftype` (Type: String)

    * “PLINK” or “BGEN” based on input type for steps 1 & 2 
## Post-Processing


* `biofilter_close_dist` (Type: Float)

    * The distance in bp for something to be considered “close” vs “far” with respect to nearest gene annotation. Value is often 5E4

* `biofilter_script` (Type: File Path)

    * The path to the biofilter script to use. If using the singularity container, should be ‘/app/biofilter.py’

* `biofilter_loki` (Type: File Path)

    * The path to a loki.db file to be used for nearest gene annotation

* `biofilter_build` (Type: String)

    * The build to pass to biofilter - can be 19 or 38

* `annotate` (Type: Bool (Java: true or false))

    * Whether or not to annotate results with the RSIDs and nearest genes for plotting and summary files.
# Output Files from SAIGE_GWAS

# Current Dockerfile for the Container/Image


```docker
# build: for python packages, plink, biofilter, and NEAT plots
FROM continuumio/miniconda3 AS build
WORKDIR /app

# biofilter version argument
ARG BIOFILTER_VERSION=2.4.3

# Combine apt operations and add retries for robustness
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        wget \
        unzip \
        libtiff5-dev \
        libtiff-dev \
        ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # Configure wget to retry downloads
    && echo "retry_connrefused = on" >> /etc/wgetrc \
    && echo "tries = 5" >> /etc/wgetrc

# Install conda packages with explicit channels and retries
RUN conda config --set remote_read_timeout_secs 600 \
    && conda install -y -n base \
        -c conda-forge \
        -c bioconda \
        adjustText \
        apsw \
        sqlite \
        dominate \
        bgenix \
        scipy \
        pandas \
        seaborn \
        matplotlib \
        conda-build \
        numpy \
        libtiff \
    && conda clean --all --yes

# Download and install tools with retries
RUN wget --tries=5 --retry-connrefused https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240526.zip \
    && wget --tries=5 --retry-connrefused https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip \
    && unzip plink2_linux_x86_64_20240526.zip \
    && unzip plink_linux_x86_64_20231211.zip \
    && rm -rf plink2_linux_x86_64_20240526.zip plink_linux_x86_64_20231211.zip

# Handle libtiff symlink
RUN find /opt/conda -name "libtiff.so.*" -type f | sort -V | tail -n 1 | xargs -I {} ln -sf {} /opt/conda/lib/libtiff.so.5 || true

# Install NEAT-plots
RUN git clone --depth 1 https://github.com/PMBB-Informatics-and-Genomics/NEAT-Plots.git \
    && mv NEAT-Plots/manhattan-plot/ /app/ \
    && rm -rf NEAT-Plots

# Install biofilter with retries
RUN wget --tries=5 --retry-connrefused \
    https://github.com/RitchieLab/biofilter/releases/download/Biofilter-${BIOFILTER_VERSION}/biofilter-${BIOFILTER_VERSION}.tar.gz \
    -O biofilter.tar.gz \
    && tar -zxvf biofilter.tar.gz --strip-components=1 -C /app/ \
    && /opt/conda/bin/python setup.py install \
    && chmod a+rx /app/biofilter.py \
    && rm -rf biofilter.tar.gz

# dev: for SAIGE and SAIGE-dependent packages
FROM rocker/tidyverse:4.1.3 AS dev
WORKDIR /tmp

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

# Install R dependencies first
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo'), repos='https://cloud.r-project.org', dependencies=TRUE)" && \
    R -e "install.packages(c('devtools', 'Matrix', 'data.table', 'SPAtest', 'roxygen2', 'rversions'), repos='https://cloud.r-project.org', dependencies=TRUE)" && \
    R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), devtools::install_github('leeshawn/MetaSKAT'))" && \
    R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), devtools::install_github('cysouw/qlcMatrix'))"

# Install SAIGE with improved error handling
RUN git clone --depth 1 https://github.com/saigegit/SAIGE.git && \
    cd SAIGE && \
    R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), devtools::install_deps('.'))" && \
    R CMD INSTALL . && \
    cd .. && \
    mv SAIGE/extdata/step1_fitNULLGLMM.R \
       SAIGE/extdata/step2_SPAtests.R \
       SAIGE/extdata/step3_LDmat.R \
       SAIGE/extdata/createSparseGRM.R \
       /usr/local/bin/ && \
    chmod -R a+x /usr/local/bin/ && \
    rm -R SAIGE

# main: final image with only necessary packages and scripts
FROM ubuntu:20.04 AS main
WORKDIR /app

# Copy files from previous stages
COPY --from=dev /tmp/ /app/
COPY --from=dev /usr/local/ /usr/local/
COPY --from=build /opt/conda/ /opt/conda/
COPY --from=build /app/ /app/

# Install R and required system libraries
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        r-base \
        ca-certificates \
        libtiff5-dev \
        libgomp1 \
        libopenblas-base \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && /opt/conda/bin/conda develop /app/manhattan-plot/ \
    && mv plink2 plink /usr/bin

# Force step_2 to use 1 single thread
ENV OMP_NUM_THREADS=1

USER root
```
# Advanced Nextflow Users: Take/Emit Info

## Input Channel (take) Description


NONE
## Output Channel (emit) Description


Singles_merge_output- a Channel of two paths to the merged raw and filtered output from SAIGE Step 2
Pheno_table- a Channel for a file that contains sample informations on each phenotype
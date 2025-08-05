
Documentation for SAIGE-GATE
============================

# Module Overview


[Example Module Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gate.config)

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

        * Input File Header:





        ```
        y_quantitative
        ```

* `bin_pheno_list` (Type: List OR File Path)

    * file path to list of binary phenotypes

    * Corresponding Input File: Binary Phenotype List

        * newline-delimited list of binary phenotypes to be tested

        * Type: List File

        * Format: txt

        * Input File Header:





        ```
        y_binary
        ```
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
## Association Test Modeling


* `cont_covars ` (Type: List)

    * Continuous covariates list
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
## Post-Processing


* `annotate` (Type: Bool (Java: true or false))

    * Whether or not to annotate results with the RSIDs and nearest genes for plotting and summary files.
# Output Files from SAIGE-GATE


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
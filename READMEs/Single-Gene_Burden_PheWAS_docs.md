
Documentation for Single-Gene Burden PheWAS
===========================================

# Module Overview


Use SAIGE to perform a gene-burden PheWAS for one or more Genes of interest and the rare variants in them

[Paper Link for Reference](https://www.nature.com/articles/s41588-022-01178-w)

[Tool Documentation Link](https://saigegit.github.io/SAIGE-doc/)

[Example Module Config File](https://github.com/PMBB-Informatics-and-Genomics/pmbb-geno-pheno-toolkit/tree/main/Example_Configs/saige_gene_phewas.config)

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

* Run Command: `nextflow run /path/to/toolkit/module/workflows/saige_gene_phewas.nf`

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


* `gene_list_file` (Type: File Path)

    * file path to newline-separated list of ENSEMBL gene ID

    * Corresponding Input File: Gene List

        * newline-delimited list of gene names to be tested

        * Type: List File

        * Format: txt

        * Input File Header:





        ```
        ENSG00000124181
        ENSG00000140463
        ENSG00000138036
        ENSG00000138002
        ENSG00000068885
        ```

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

* `sex_strat_cohort_list` (Type: List)

    * List of cohorts that are sex stratified
## Pre-Processing


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
## QC Options


* `min_maf` (Type: Float)

    * Minimum minor allele frequency for plink QC
## Association Test Modeling


* `cont_covars ` (Type: List)

    * Continuous covariates list

* `sex_strat_cat_covars` (Type: List)

    * Categorical covariates for sex stratified cohorts to ensure model converges
## SAIGE Step 1


* `geno` (Type: Float)

    * Plink parameters for SAIGE Step 1 Input QC which needs a small set of high-quality variants, genotype rate  filters out all variants with missing call rates exceeding the provided value 

* `exome_plink_prefix` (Type: Plink Fileset Prefix)

    * Exome plink input files 

    * Corresponding Input File: SAIGE Exome Plink Files

        * a hard-call plink set of exome data

        * Type: Plink Set

        * Format: plink binary

        * Input File Header:





        ```
        PMBB-Release-2020-2.0_genetic_exome_GL_norm{.bed,.bim,.fam,.log,.pgen,.psam,.pvar}
        ```
## SAIGE Step 2


* `min_mac ` (Type: Float)

    * SAIGE-GENE Step 2 Parameters

* `firth_cutoff` (Type: Float)

    * P-value ()
## Post-Processing


* `pheno_descriptions_file` (Type: File Path)

    * file path to phenotype descriptions used for plotting and summary statistics

    * Corresponding Input File: Phenotype Descriptions File

        * File containing Phenotype notation,descriptions,and, categories (i.e Phe8.52,Intestinal infection due to C. difficile,Infection)

        * Type: Data Table

        * Format: csv

        * Input File Header:





        ```
        PHENO,DESCRIPTION,CATEGORY
        x1,description_x1,RED
        x2,description_x2,BLUE
        ```

* `gene_location_file ` (Type: File Path)

    * This file is used for getting gene-based coordinates for plotting .

    * Corresponding Input File: Gene Location File

        * CSV file of 

        * Type: Data Table

        * Format: tsv

        * Input File Header:





        ```
        gene_id chromosome  seq_region_start    seq_region_end  gene_symbol
        GENE1   1   1   90  GS1
        GENE2   2   91  100 GS2
        ```

* `region_col_names` (Type: Map (Dictionary))

    * Default SAIGE Region column names mapped to new ones
# Output Files from Single-Gene_Burden_PheWAS


* Singles Summary Statistics Filtered

    * A csv FILTERED summary file with cohort, phenotype, chromosome, base pair location, variant ID, Other allele, effect allele, effect,beta, statistics, pvalues, and other sumstats. One file will be generated for each combination of cohort, chromosome, and analysis (regular, cauchy, rare, ultra rare)

    * Type: Summary Statistics

    * Format: csv

    * Parallel By: Cohort, Chromosome

    * Output File Header:





    ```
    cohort,phenotype,chromosome,base_pair_location,variant_id,other_allele,effect_allele,effect_allele_count,effect_allele_frequency,missing_rate,beta,standard_error,t_statistic,variance,p_value,p_value_na,is_spa_test,allele_freq_case,allele_freq_ctrl,n_case,n_ctrl,n_case_hom,n_case_het,n_ctrl_hom,n_ctrl_het
    EUR_ALL,Phe783.0,4,39274909,4_39274909_C_T,C,T,77,0.00131171,0.0,2.6324,0.606004,5.43601,0.54679,6.999681e-06,1.961538e-13,True,0.015,0.0012178,200,29151,0,6,0,71
    ```

* Regions Summary Statistics Filtered

    * A FILTERED top hits csv summary file of results including phenotype, gene, group annotation, p-values, and other counts. One file will be generated for each combination of cohort, chromosome, and analysis (regular, cauchy, rare, ultra rare)

    * Type: Summary Statistics

    * Format: csv

    * Parallel By: Cohort, Chromosome

    * Output File Header:





    ```
    cohort,phenotype,gene,annot,max_maf,p_value,p_value_burden,p_value_skat,beta_burden,se_burden,mac,mac_case,mac_control,rare_var_count,ultrarare_var_count
    EUR_ALL,Phe473.3,ENSG00000157796,pLoF;damaging_missense;other_missense;synonymous,0.001,9.85895741467679e-06,2.30630516571129e-05,3.36463023578037e-05,0.073107224589224,0.0172709262685344,1059.0,17.0,1042.0,18.0,347.0
    ```

* Gene PheWAS QQ Plots

    * Type: QQ Plot

    * Format: png

    * Parallel By: Cohort, Gene, Annot Group, MAF

* Phenotype Summary Table

    * A csv file with summary of all cohort, phenotypes, Total Number, number of cases, Number of Controls, and Prevalence

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    COHORT,PHENO,N,Controls,Cases,Prevalence
    AFR_ALL,Phe10.0,11098,11074.0,24.0,0.002162551811137142
    AFR_ALL,Phe1000.0,10959,10884.0,75.0,0.00684369011771147
    AFR_ALL,Phe1001.0,10818,10637.0,181.0,0.016731373636531707
    ```

* Gene PheWAS Plots

    * A dot plot (manhattan plot) of phenotypes (usually phecodes) rolled into their respective disease categories plotted against their p-values

    * Type: PheWAS Plot

    * Format: png

    * Parallel By: Cohort, Gene, Annot Group, MAF

* Singles Top Hits Table

    * A csv summary file with cohort, phenotype, chromosome, base pair location, variant ID, Other allele, effect allele, effect,beta, statistics, pvalues, and other sumstats. This are all the “filtered” summary files joined together. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    cohort,phenotype,chromosome,base_pair_location,variant_id,other_allele,effect_allele,effect_allele_count,effect_allele_frequency,missing_rate,beta,standard_error,t_statistic,variance,p_value,p_value_na,is_spa_test,allele_freq_case,allele_freq_ctrl,n_case,n_ctrl,n_case_hom,n_case_het,n_ctrl_hom,n_ctrl_het
    AFR_F,Phe654.0,2,27444495,2_27444495_G_A,G,A,11,0.00080292,0.0,0.481686,0.0759042,-0.430299,1.33634e-15,1.1e-10,0.0E-2147483648,True,0.0,0.000846414,352,6498,0,0,0,11
    AFR_F,Phe643.0,2,27444495,2_27444495_G_A,G,A,11,0.000805153,0.0,12.1317,2.41425,-0.274659,1.47634e-15,2.5e-07,0.0E-2147483648,True,0.0,0.000830189,206,6625,0,0,0,11
    EUR_F,Phe644.0,2,27445431,2_27445431_C_T,C,T,43,0.00161849,0.0,9.01836,1.51694,-0.399219,5.45544e-15,1.4e-09,0.0E-2147483648,True,0.0,0.00163623,144,13140,0,0,0,43
    ```

* Regions Top Hits Table

    * A csv summary file with cohort, phenotype, gene, group annotation, max maf, p-value, p-value for burden, p-value SKAT, p-value for any tests defined in config, mac, mac for case and control, rare variant count, and ultra rare variant count. This are all the “filtered” summary files joined together. 

    * Type: Summary Table

    * Format: csv

    * Output File Header:





    ```
    cohort,phenotype,gene,annot,max_maf,p_value,p_value_burden,p_value_skat,beta_burden,se_burden,mac,mac_case,mac_control,rare_var_count,ultrarare_var_count
    EUR_F,Phe654.2,ENSG00000068885,other_missense,0.0001,2.81627590898728e-06,2.81627590898728e-06,2.81627590898728e-06,3.12075302693361,0.666288927246197,60.0,1.0,59.0,0.0,50.0
    EUR_F,Phe665.0,ENSG00000068885,other_missense,0.01,2.31299783622365e-41,7.46670792353132e-05,3.3042826231766494e-42,-0.9560950788528,0.24138315092594,348.0,7.0,341.0,5.0,62.0
    EUR_F,Phe665.0,ENSG00000068885,synonymous,0.0001,1.6691464931375e-07,1.6691464931375e-07,1.6691464931375e-07,-3.00605241896674,0.574457788692838,16.0,0.0,16.0,0.0,13.0
    ```

* Singles Summary Statistics

    * A gzipped tsv summary file with  phenotype, chromosome, base pair location, variant ID, Other allele, effect allele, effect,beta, statistics, pvalues, and other sumstats. One file will be generated for each combination of cohort, chromosome, and analysis (regular, cauchy, rare, ultra rare). One file will be generated for each cohort, chormosome, and analysis (regular, cauchy, rare, ultra rare) combination. 

    * Type: Summary Statistics

    * Format: tsv.gz

    * Parallel By: Cohort, Chromosome

    * Output File Header:





    ```
    phenotype       chromosome      base_pair_location      variant_id      other_allele    effect_allele   effect_allele_count     effect_allele_frequency missing_rate beta    standard_error  t_statistic     variance        p_value p_value_na      is_spa_test     allele_freq_case        allele_freq_ctrl        n_case       n_ctrl  n_case_hom      n_case_het      n_ctrl_hom      n_ctrl_het
    Phe159.2        4       39199561        4_39199561_G_A  G       A       39      0.000652218     0.0     -1.00478        3.11147 -0.103787       0.103293    0.7467486        0.7467486       False   0.0     0.000653858     75      29823   0       0       0       39
    Phe159.2        4       39214639        4_39214639_A_G  A       G       14      0.000234129     0.0     -1.00338        5.62099 -0.0317569      0.0316501   0.8583262        0.8583262       False   0.0     0.000234718     75      29823   0       0       0       14
    Phe159.2        4       39215918        4_39215918_C_T  C       T       41      0.000685665     0.0     -1.00617        3.49319 -0.0824573      0.0819513   0.7733171        0.7733171       False   0.0     0.000687389     75      29823   0       0       0       41
    ```

* Regions Summary Statistics

    * A tsv summary file of results including phenotype, gene, group annotation, p-values, and other counts. One file will be generated for each combination of cohort, chromosome, and analysis (regular, cauchy, rare, ultra rare)

    * Type: Summary Statistics

    * Format: tsv.gz

    * Parallel By: Cohort, Chromosome

    * Output File Header:





    ```
    phenotype       gene    annot   max_maf p_value p_value_burden  p_value_skat    beta_burden     se_burden       mac     mac_case        mac_control     rare_var_count       ultrarare_var_count
    Phe159.2        ENSG00000179029 pLoF    0.0001  0.933222000000003       0.933222000000003       0.933222000000003       -0.0401639507465058     0.479330557312422    2.0     0.0     2.0     0.0     1.0
    Phe159.2        ENSG00000179029 pLoF    0.001   0.933222000000003       0.933222000000003       0.933222000000003       -0.0401639507465058     0.479330557312422    2.0     0.0     2.0     0.0     1.0
    Phe159.2        ENSG00000179029 pLoF    0.01    0.933222000000003       0.933222000000003       0.933222000000003       -0.0401639507465058     0.479330557312422    2.0     0.0     2.0     0.0     1.0
    Phe159.2        ENSG00000179029 damaging_missense       0.0001  0.843473100000001       0.843473100000001       0.843473100000001       -0.0403718054878391 0.20446291771094 13.0    0.0     13.0    0.0     9.
    ```
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
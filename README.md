Documentation for SAIGE Pipelines
=============================

# Overview

SAIGE is a popular tool for performing association tests with multiple types of genetic data, such as:
- Common Variants
- Rare Variants
- Gene Burden

Here are some SAIGE-specific resources:
- [Original SAIGE Paper for Your Reference](https://www.nature.com/articles/s41588-018-0184-y)
- [SAIGE-GENE Paper for Your Reference](https://www.nature.com/articles/s41588-022-01178-w)
- [SAIGE Documentation](https://saigegit.github.io/SAIGE-doc/)

# Software Requirements:
- [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- Some text editor for updating the `nextflow.config` profiles such as `vim` or `nano`
- [Nextflow version 23.04.1.5866](https://www.nextflow.io/docs/latest/cli.html)
- [Singularity version 3.8.3](https://sylabs.io/docs/) OR [Docker version 4.30.0](https://docs.docker.com/)
- [JDK version 11.0.5](https://www.oracle.com/java/technologies/javase/jdk11-archive-downloads.html)

# Steps for Runnning One of the Pipelines with Test Data

Note: test data were obtained from the [SAIGE github repo](https://github.com/saigegit/SAIGE).

## Steps:
1. Start your own project directory and go there
    - `mkdir my_new_saige_project`
    - `cd my_new_saige_project`
2. Build the `saige.sif` singularity image
    - `singularity build saige.sif docker://pennbiobank/saige:latest`
    - you do not necessarily have to do this repeatedly, BUT it's good for versioning to keep it with the project
    - if you want to reuse containers, you could build containers in a more general `containers` directory somewhere
3. Download the source code by cloning from git
    - `git clone https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family.git`
    - You may do this in your project directory, but it often makes sense to clone into a general `tools` location
4. Copy the contents of one of the test data directories
    - `cp -r pmbb-nf-toolkit-saige-family/test_data/ExWAS/test_config_exwas_no_GRM/* .`
5. Fill out the `nextflow.config` file to make sure it matches your system
    - See `nextflow` executor information [here](https://www.nextflow.io/docs/latest/executor.html).
    - From the command line, you can use the text editor `vim` with `vim nextflow.config`
    - The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to` with the location where you built the image in step 2)
6. Do a stub run to test
    - `nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -stub`
    - Note: Replace `/path/to/` with the actual path (full or relative) where you cloned the repository in step 3.
    - ^ make sure you select the correct `profile` you set up
7. Now run the workflow on the actual data!
    - `nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster`
    - Note: Replace `/path/to/` with the actual path (full or relative) where you cloned the repository in step 3.

# Steps for Runnning One of the Pipelines with YOUR Data

## Steps:
1. Start your own project directory and go there
    - `mkdir my_new_saige_project`
    - `cd my_new_saige_project`
2. Build the `saige.sif` singularity image
    - `singularity build saige.sif docker://pennbiobank/saige:latest`
    - you do not necessarily have to do this repeatedly, BUT it can be good for versioning to keep it with the project in case it gets updated later
    - if you want to reuse containers, you could build containers in a more general `containers` directory somewhere
3. Download the source code by cloning from git
    - `git clone https://github.com/PMBB-Informatics-and-Genomics/pmbb-nf-toolkit-saige-family.git`
    - You may do this in your project directory, but it often makes sense to clone into a general `tools` location
4. Set up the appropriate config file for the pipeline you want to run
    - Update any and all input files!
    - Double-check parameters for running the software (QC filters, settings, etc)
    - Refer to the READMEs directory with more info for every parameter and input file
5. Fill out the `nextflow.config` file to make sure it matches your system
    - See `nextflow` executor information [here](https://www.nextflow.io/docs/latest/executor.html).
    - From the command line, you can use the text editor `vim` with `vim nextflow.config`
    - The profile's attribute `process.container` should be set to `'/path/to/saige.sif'` (replace `/path/to/` with the path where you built the image in step 2)
    - Make sure the `nextflow.config` file has an `includeConfig` statement
6. Do a stub run to test
    - `nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster -stub`
    - Note: Replace `/path/to/` with the actual path (full or relative) where you cloned the repository in step 3.
    - ^ make sure you select the correct `profile` you set up
7. Now run the workflow on the actual data!
    - `nextflow run /path/to/pmbb-nf-toolkit-saige-family/workflows/saige_exwas.nf -profile cluster`
    - Note: Replace `/path/to/` with the actual path (full or relative) where you cloned the repository in step 3.

## Notes to keep in mind when applying this process to your data:
- You do not need to copy YOUR data to the Input/ folder (unless you want to). It usually makes sense to add full paths to your `configs/workflow_specific_config` as needed.
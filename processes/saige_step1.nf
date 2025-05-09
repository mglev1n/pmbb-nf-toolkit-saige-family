process plink_qc_for_step1 {
    publishDir "${cohort_dir}/Saige_Step1/"
    machineType 'n2-standard-4'
    memory { params.host == 'AOU' ? '63GB' : '24GB' }
    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)
    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'
    script:
        use_mem = params.host == 'AOU' ? '62000' : '23000'
        """
        plink2 --indep-pairwise 50 5 0.4 \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --memory ${use_mem} \
            --maf ${params.maf} --geno ${params.geno} --not-chr X,Y,MT,XY,PAR1,PAR2 --hwe ${params.hwe} --snps-only \
            --out step1

        shuf -n ${params.thin_count} step1.prune.in >> step1.markerid.list

        #  make input file for step 1 using variant list
        stdbuf -e0 -o0 plink2 --make-bed \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --memory ${use_mem} \
            --extract step1.markerid.list \
            --out step1_plink > plink_qc.log
        """
    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

process plink_qc_for_step1_saige_gene {
    publishDir "${cohort_dir}/Saige_Step1/"
    machineType 'n2-standard-4'
    memory { params.host == 'AOU' ? '63GB' : '24GB' }
    /* commmon variants : ld pruning 50 5 0.4
    / 10-20 : plink --counts --snplist -> 1000
    / 20-430 : plink --counts --snplist -> 1000
    / plink --extract
    */
    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)
    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'
    script:
        use_mem = params.host == 'AOU' ? '62000' : '23000'
        """

        # get counts from whole plink file
        # needs to be plink2 executable to get correct file format out
        plink2 --freq counts \
          --keep ${sample_list} \
          --max-maf 0.01 \
          --mac 10 \
          --bfile ${input_bed.toString().replace('.bed', '')} \
          --out step1_counts

        awk '{ if ((\$5 >= 10 && \$5 <= 20) || (2*\$6-\$5 >= 10 && 2*\$6-\$5 <= 20)) { print \$2} }' ./step1_counts.acount > temp_mac10to20_ids.txt
        awk '{ if ((\$5 >= 20) || (2*\$6-\$5 >= 20)) { print \$2} }' step1_counts.acount > temp_mac20ormore_ids.txt

        num_mac10=\$(< temp_mac10to20_ids.txt wc -l)
        num_mac20=\$(< temp_mac20ormore_ids.txt wc -l)

        if [ "\$num_mac10" -lt 1000 ]; then
            echo "Error: Less than 1000 SNPs with MAC between 10 and 20"
        elif [ "\$num_mac20" -lt 1000 ]; then
            echo "Error: Less than 1000 SNPs with MAC greater than 20"
        else
            # Copy the filtered IDs to the output file
            shuf -n 1000 temp_mac20ormore_ids.txt > ./step1.markerid.list
            shuf -n 1000 temp_mac10to20_ids.txt >> ./step1.markerid.list
            echo "Filtered IDs written to step1.markerid.list"
        fi

        rm temp_mac10to20_ids.txt

        plink2 --indep-pairwise 50 5 0.4 \
            --bfile ${input_bed.toString().replace('.bed', '')} \
            --keep ${sample_list} \
            --maf ${params.maf} --geno ${params.geno} \
            --not-chr X,Y,MT --hwe ${params.hwe} --snps-only \
            --out step1

        # create snp list with variants from all categories
        shuf -n ${params.thin_count} step1.prune.in >> step1.markerid.list

        #  make input file for step 1 using variant list
        stdbuf -e0 -o0 plink2 --make-bed \
          --bfile ${input_bed.toString().replace('.bed', '')} \
          --keep ${sample_list} \
          --not-chr X,Y,MT,XY,PAR1,PAR2 \
          --extract step1.markerid.list \
          --out step1_plink > plink_qc.log

        rm step1_counts.acount
        rm step1.prune.in
        """
    stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

String get_covar_list_args(String cohort, cohort_cat_covars, cohort_cont_covars) {
    String output = ''

    if (cohort_cont_covars.size() > 0 || cohort_cat_covars.size() > 0) {
        output += '--covarColList='
        if (cohort_cont_covars.size() > 0) {
            output += cohort_cont_covars.join(',')
        }
        if (cohort_cat_covars.size() > 0) {
            output += ',' + cohort_cat_covars.join(',') + ' '
        }
        else {
            output += ' '
        }
    }

    if (cohort_cat_covars.size() > 0) {
        output += '--qCovarColList=' + cohort_cat_covars.join(',') + ' '
    }

    return output
}
/*
process call_saige_step1_bin {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    machineType 'n2-standard-4'

    label(params.GPU == 'ON' ? 'gpu_on' : 'gpu_off')

    errorStrategy 'retry'

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        if [ "${params.GPU}"  = "ON" ]; then
            mpirun -n 8 singularity exec --bind /project/:/project/,./:/output,/static/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --nThreads=1 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 Rscript ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --nThreads=29 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi

        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}
*/
process dnanexus_skip_plink_qc_for_step1 {
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort_dir), path(sample_list)

        // parameters
        val plink_prefix

        // new inputs
        tuple path(input_bed), path(input_bim), path(input_fam)

    output:
        // variables
        tuple val(cohort_dir), path('step1_plink.{bed,bim,fam}')
        path 'plink_qc.log'

    script:
    """
    # Check if the required PLINK files exist in the launch directory
    if [[ ! -f "${plink_prefix}.bed" || ! -f "${plink_prefix}.bim" || ! -f "${plink_prefix}.fam" ]]; then
        echo "Error: One or more PLINK files (.bed, .bim, .fam) are missing in the launch directory." >&2
        exit 1
    fi

    # Create symbolic links to mimic output files
    ln -s "${plink_prefix}.bed" step1_plink.bed
    ln -s "${plink_prefix}.bim" step1_plink.bim
    ln -s "${plink_prefix}.fam" step1_plink.fam

    # Create a dummy log file
    echo "PLINK file check completed successfully." > plink_qc.log
    """

stub:
        '''
        touch step1_plink.bed
        touch step1_plink.bim
        touch step1_plink.fam
        touch plink_qc.log
        '''
}

process dnanexus_skip_step1_quant {
    //for this to work, output files must exist in the launch dir
    machineType 'n2-standard-4'
    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    script:
    """
    # Check if the required PLINK files exist in the launch directory
    if [[ ! -f "${pheno}.varianceRatio.txt" || ! -f "${pheno}.rda" || ! -f "${pheno}.log" ]]; then
        echo "Error: One or more PLINK files (.bed, .bim, .fam) are missing in the launch directory." >&2
        exit 1
    fi

    # Create symbolic links to mimic output files
    ln -s "${pheno}.varianceRatio.txt" ${pheno}.varianceRatio.txt
    ln -s "${pheno}.rda" ${pheno}.rda
    ln -s "${pheno}.log" ${pheno}.log

    # Create a dummy log file
    echo "STEP 1 file check completed successfully." > step1.log
    """

stub:
    """
    touch ${pheno}.varianceRatio.txt
    touch ${pheno}.rda
    touch ${pheno}.log
    touch step1.log
    """
}

process dnanexus_skip_step1_bin {
    input:
    tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)

    output:
    tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
    path "${pheno}.log"

    script:
    """
    # Use Nextflow's launchDir
    if [[ ! -f "${launchDir}/${pheno}.varianceRatio.txt" || ! -f "${launchDir}/${pheno}.rda" || ! -f "${launchDir}/${pheno}.log" ]]; then
        echo "Error: One or more files are missing in ${launchDir}" >&2
        ls -la "${launchDir}"
        exit 1
    fi

    ln -sf "${launchDir}/${pheno}.varianceRatio.txt" ./${pheno}.varianceRatio.txt
    ln -sf "${launchDir}/${pheno}.rda" ./${pheno}.rda
    ln -sf "${launchDir}/${pheno}.log" ./${pheno}.log

    echo "STEP 1 file check completed successfully." > step1.log
    """
}

process call_saige_step1_bin {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    errorStrategy 'retry'
    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)

    machineType { (params.GPU == 'ON') ? getGpuMachineTypeChannel(plink_set[1].readLines().size, pheno_file.readLines().size) : 'n2-standard-4' }

    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        if [ "${params.GPU}" = "ON" ] && [ "${params.host}" = "LPC" ]; then
            mpirun -n 8 singularity exec --bind /project/:/project/,./:/output,/static/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --nThreads=1 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}"  = "ON" ] &&  [ "${params.host}" = "DNAnexus" ]; then
            mpirun -n 8 --allow-run-as-root \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno}\
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --nThreads=1 \
            --LOCO=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
         elif [ "${params.GPU}"  = "ON" ] &&  [ "${params.host}" = "AOU" ]; then
            mpirun -n 8 --allow-run-as-root \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno}\
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --nThreads=1 \
            --LOCO=FALSE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 Rscript ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=binary \
            --outputPrefix=${pheno} \
            --LOCO=${params.LOCO} \
            --nThreads=29 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi

        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_quant {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    machineType 'n2-standard-16'
    cpus 15

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        if [ "${params.GPU}"  = "ON" ]; then
            mpirun -n 8 singularity exec --bind /project/:/project/,./:/output,/static/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=quantitative \
            --outputPrefix=${pheno} \
            --nThreads=1 \
            --invNormalize=TRUE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 Rscript ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --traitType=quantitative \
            --outputPrefix=${pheno} \
            --nThreads=29 \
            --LOCO=${params.LOCO} \
            --invNormalize=TRUE \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi

        """
    stub:
        """
        touch ${pheno}.log
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        """
}

process call_saige_step1_survival {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    machineType 'n2-standard-4'
    errorStrategy 'retry'

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        if [ "${params.GPU}" = "ON" ] && [ "${params.host}" = "LPC" ]; then
            mpirun -n 8 singularity exec --bind /project/:/project/,./:/output,/static/ --nv \
            saige-doe-3.sif \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R  \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --outputPrefix=${pheno} \
            --traitType=survival \
            --eventTimeCol=${params.event_time_col} \
            --eventTimeBinSize=${params.event_time_bin} \
            --invNormalize=TRUE \
            --nThreads=29 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}"  = "ON" ] &&  [ "${params.host}" = "DNAnexus" ]; then
            mpirun -n 8 --allow-run-as-root \
            /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --outputPrefix=${pheno} \
            --traitType=survival \
            --eventTimeCol=${params.event_time_col} \
            --eventTimeBinSize=${params.event_time_bin} \
            --invNormalize=TRUE \
            --nThreads=29 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        elif [ "${params.GPU}" = "OFF" ]; then
            echo "${cohort}-${pheno}"
            stdbuf -e0 -o0 Rscript ${params.step1_script} \
            --plinkFile=step1_plink \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno} \
            ${covariate_args} \
            --sampleIDColinphenoFile=${params.id_col} \
            --outputPrefix=${pheno} \
            --traitType=survival \
            --eventTimeCol=${params.event_time_col} \
            --eventTimeBinSize=${params.event_time_bin} \
            --invNormalize=TRUE \
            --LOCO=${params.LOCO} \
            --nThreads=29 \
            --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        fi
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_bin_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    machineType 'n2-standard-4'
    cpus 1

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --plinkFile=${plink_set[0].toString().replace('.bed', '')} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --traitType=binary \
         --outputPrefix=${pheno} \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_quant_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    machineType 'n2-standard-4'
    cpus 1

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --plinkFile=${plink_set[0].toString().replace('.bed', '')} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --outputPrefix=${pheno} \
         --traitType=quantitative \
         --invNormalize=TRUE \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

process call_saige_step1_survival_with_sparse_GRM {
    publishDir "${launchDir}/${cohort}/Saige_Step1/"
    machineType 'n2-standard-4'
    cpus 1

    input:
        // variables
        tuple val(cohort), val(pheno), path(plink_set), path(samples), path(pheno_file)
        tuple path(sparse_grm), path(sparse_grm_samples)
    output:
        // variables
        tuple val(cohort), val(pheno), path("${pheno}.rda"), path("${pheno}.varianceRatio.txt")
        path "${pheno}.log"
    shell:
        covariate_args = get_covar_list_args(cohort,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cat_covars : params.cat_covars,
            params.sex_strat_cohort_list.contains(cohort) ? params.sex_strat_cont_covars : params.cont_covars)
        """
        echo "${cohort}-${pheno}"
        stdbuf -e0 -o0 Rscript ${params.step1_script} \
         --sparseGRMFile=${sparse_grm} \
         --sparseGRMSampleIDFile=${sparse_grm_samples} \
         --useSparseGRMtoFitNULL=TRUE \
         --plinkFile=${plink_set[0].toString().replace('.bed', '')} \
         --phenoFile=${pheno_file} \
         --phenoCol=${pheno} \
         ${covariate_args} \
         --sampleIDColinphenoFile=${params.id_col} \
         --outputPrefix=${pheno} \
         --traitType=survival \
         --eventTimeCol=${params.event_time_col} \
         --eventTimeBinSize=${params.event_time_bin} \
         --invNormalize=TRUE \
         --nThreads=29 \
         --IsOverwriteVarianceRatioFile=TRUE > ${pheno}.log
        """
    stub:
        """
        touch ${pheno}.rda
        touch ${pheno}.varianceRatio.txt
        touch ${pheno}.log
        """
}

include {
    paramToList
    getGpuMachineTypeChannel
} from '../processes/saige_helpers.nf'

workflow SAIGE_STEP1 {
    take:
        cohort_sample_lists
        cohort_pheno_tables
        keep_cohort_bin_pheno_combos
        keep_cohort_quant_pheno_combos
        keep_cohort_survival_pheno_combos
        use_step1_prefix
        IS_GENE
    main:
        bin_pheno_list = paramToList(params.bin_pheno_list)
        bin_pheno = Channel.fromList(bin_pheno_list)
        quant_pheno_list = paramToList(params.quant_pheno_list)
        quant_pheno = Channel.fromList(quant_pheno_list)
        survival_pheno_list = paramToList(params.survival_pheno_list)
        survival_pheno = Channel.fromList(survival_pheno_list)
        plink_suffixes_list = ['.bed', '.bim', '.fam']

        // First decide how to set up plink files - no GRM means we have to run QC. with GRM, we can use them directly
        if (params.use_sparse_GRM) {
            // We can use the input plink files directly - they should be the set of variants used to make the GRM
            cohort_plinkset = cohort_sample_lists.map { cohort, sample_list -> \
                new Tuple(cohort, new Tuple(*plink_suffixes_list.collect { ext -> use_step1_prefix + ext }))
            }
        }
        else {
            // set up inputs for plink_qc
            plink_inputs_tuple = new Tuple(*plink_suffixes_list.collect { ext -> use_step1_prefix + ext })

            if (IS_GENE) {
                // call plink_qc for each cohort (gene version, samples rare variants)
                (cohort_plinkset, logs) = plink_qc_for_step1_saige_gene(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
            }
            else {
                if (params.host == 'DNAnexus') {
                    // Check for existing PLINK files instead of running QC
                    (cohort_plinkset, logs) = dnanexus_skip_plink_qc_for_step1(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
                }
                else {
                    // call plink_qc for each cohort
                    (cohort_plinkset, logs) = plink_qc_for_step1(cohort_sample_lists, use_step1_prefix, plink_inputs_tuple)
                }
            }
        }

        // Add the sample lists and phenotype tables in so we get:
        // (cohort, (plink_trio), sample_list, pheno_table)
        cohort_plinkset_sample_pheno = cohort_plinkset.join(cohort_sample_lists).join(cohort_pheno_tables)

        /*
        Step 1 Input Emission Tuples
        Need: (cohort, phenotype, (plink set), sample_file, pheno_file)
        Have:
        (cohort, (plink_set), sample_file, pheno_file)
        (phenotype)
        Steps:
        1. .combine(bin_pheno or quant_pheno) -> (cohort, (plink set), sample_file, pheno_file, phenotype)
        2. .map() -> rearrange to get cohort, pheno first
        3. .join() -> keep only cohort/pheno combos with enough sample size
        */
        step1_bin_input = cohort_plinkset_sample_pheno \
            .combine(bin_pheno) \
            .map { cohort, plink, sample_file, pheno_file, pheno -> \
                    new Tuple(cohort, pheno, plink, sample_file, pheno_file) } \
            .join(keep_cohort_bin_pheno_combos, by: [0, 1])

        step1_quant_input = cohort_plinkset_sample_pheno \
            .combine(quant_pheno) \
            .map { cohort, plink, sample_file, pheno_file, pheno -> \
                    new Tuple(cohort, pheno, plink, sample_file, pheno_file) } \
            .join(keep_cohort_quant_pheno_combos, by: [0, 1])

        step1_survival_input = cohort_plinkset_sample_pheno \
            .combine(survival_pheno) \
            .map { cohort, plink, sample_file, pheno_file, pheno -> \
                    new Tuple(cohort, pheno, plink, sample_file, pheno_file) } \
            .join(keep_cohort_survival_pheno_combos, by: [0, 1])

        // Now figure out which version of Step 1 We're Calling
        if (params.use_sparse_GRM) {
            sparse_grm_input = new Tuple(params.step1_sparse_grm, params.step1_sparse_grm_samples)

            (step1_bin_output, logs) = call_saige_step1_bin_with_sparse_GRM(step1_bin_input, sparse_grm_input)
            (step1_quant_output, logs) = call_saige_step1_quant_with_sparse_GRM(step1_quant_input, sparse_grm_input)
            (step1_surival_output, logs) = call_saige_step1_survival_with_sparse_GRM(step1_survival_input, sparse_grm_input)
        } else {
            if (params.host == 'DNAnexus') {
                (step1_bin_output, logs) = dnanexus_skip_step1_bin(step1_bin_input)
            }
            //lpc or aou run
            else {
                (step1_bin_output, logs) = call_saige_step1_bin(step1_bin_input)
            }
            //if (params.host == "DNAnexus") {
            //    (step1_bin_output, logs) = dnanexus_skip_step1_quant(step1_quant_input)
            //}
            //lpc or aou run
            //else{
            (step1_quant_output, logs) = call_saige_step1_quant(step1_quant_input)
            //}
            (step1_surival_output, logs) = call_saige_step1_survival(step1_survival_input)
        }
    emit:
        step1_bin_output
        step1_quant_output
        step1_surival_output
}

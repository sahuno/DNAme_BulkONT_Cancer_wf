#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.vcfs = []

process mergeVcfs {
    container '/data1/greenbab/users/ahunos/apps/containers/gatk.sif'

    input:
    path vcfs

    output:
    path 'merged_somatic.vcf.gz'

    script:
    def input_vcfs = vcfs.collect{ "-I ${it}" }.join(' ')
    """
    gatk MergeVcfs \\
        ${input_vcfs} \\
        -O merged_somatic.vcf.gz
    """
}

workflow {
    // Convert params.vcfs to a list if it's not already one.
    def vcfs = params.vcfs instanceof List ? params.vcfs : params.vcfs.tokenize()
    
    if (vcfs.size() < 2) {
        error("You must provide at least two VCF files: --vcfs 'sample1.vcf.gz sample2.vcf.gz'")
    }
    
    mergeVcfs(vcfs)
}



// nextflow run merge_vcfs.nf --vcfs 'sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz' -with-singularity

// nextflow run /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/nextflow/merge_vcfs.nf \
// --vcfs '/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/archived/results/call_snps_indels/D-A-1_4000/snv.vcf.gz /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/archived/results/call_snps_indels/D-A-2_4000/snv.vcf.gz /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/archived/results/call_snps_indels/D-A-3_4000/snv.vcf.gz' \
// -with-singularity /data1/greenbab/users/ahunos/apps/containers/gatk.sif

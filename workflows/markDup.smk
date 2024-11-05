rule mark_duplicates:
    input:
        bams="/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/hac_5mCG_5hmCG/results/mod_bases/044T_v14/044T_v14_modBaseCalls_sorted.bam"
    output:
         markdup_bam="results/mark_duplicates/044T_v14_modBaseCalls_sorted_dup.bam",
         markdup_bam_bai="results/mark_duplicates/044T_v14_modBaseCalls_sorted_dup.bai",
         markdup_stats="results/mark_duplicates/044T_v14_modBaseCalls_marked_dup_metrics.txt"
    log:
        "logs/mark_duplicates/044T_v14.log"
    singularity:
        "/data1/greenbab/users/ahunos/apps/containers/gatk.sif"
    shell:
        """ 
        gatk MarkDuplicates --INPUT {input.bams} --OUTPUT {output.markdup_bam} --METRICS_FILE {output.markdup_stats} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT 2> {log}
        """

# snakemake --use-singularity --cores all
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/markDup.smk --jobs 10 --cores all --use-singularity --singularity-args "--bind /data1/greenbab" -np
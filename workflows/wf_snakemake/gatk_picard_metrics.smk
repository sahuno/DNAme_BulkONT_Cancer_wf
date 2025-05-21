

# Snakefile
import os

# load both config files
configfile: "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/samplingRate4000_n_5000_debub_modbams.yaml"
configfile: "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/picard_config.yaml"
# "config.yaml",

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/gatk_picard_metrics.smk \
# --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" \
# --workflow-profile /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/slurmMinimal \
#   --cores all -np

# outdir= "results"
# extract sample names from the samples mapping
SAMPLES = list(config["samples"].keys())

rule all:
    input:
        expand("{outdir}/metrics/{sample}.picard_alignment_metrics.txt",outdir=config["outdir"],sample=SAMPLES)

rule collect_alignment_metrics:
    """
    Run Picard CollectAlignmentSummaryMetrics on each sorted BAM defined in samples.yaml.
    """
    input:
        bam=lambda wc: config["samples"][wc.sample],
        # bai=lambda wc: config["samples"][wc.sample] + ".bai"
    output:
        metrics="{outdir}/metrics/{sample}.picard_alignment_metrics.txt"
    params:
        reference=config["reference"]
    threads: 1
    resources:
        mem_mb=4000
    log:
        "logs/{outdir}/{sample}.collect_metrics.log"
    singularity:
        config.get("picard_container", None)
    shell:
        r"""
        gatk CollectAlignmentSummaryMetrics \
      -I {input.bam} \
      -O {output.metrics} \
      -R {params.reference} > {log} 2>&1
        """
        # mkdir -p {outdir}/metrics logs

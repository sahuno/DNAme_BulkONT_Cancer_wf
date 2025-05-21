# Snakefile
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/mosdepth_qc.smk \
# --workflow-profile /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/slurmMinimal \
# --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --quiet -R multiqc --forceall -np

ruleorder: multiqc > mosdepth
import os

configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/haplotagged_bam_samples_Rate4000_n_Rate5000_merged.yaml"
SAMPLES = config["samples"]

SINGULARITY_IMAGE_MOSDEPTH = "/data1/greenbab/users/ahunos/apps/containers/mosdepth_multiarch.sif"
SINGULARITY_IMAGE_MULTIQC = "/data1/greenbab/users/ahunos/apps/containers/multiqc_latest.sif"

rule all:
    input:
        expand("results/mosdepth/{sample}.mosdepth.summary.txt", sample=SAMPLES),
        "results/multiQC_mosdepth/multiqc_report.html"

rule mosdepth:
    """
    Run per-sample coverage QC via mosdepth inside your SIF.
    """
    input:
        bam=lambda wc: config["samples"][wc.sample]
        # bam="data/{sample}.bam"
        # bai="data/{sample}.bam.bai"
    output:
        summary="results/mosdepth/{sample}.mosdepth.summary.txt",
        dist="results/mosdepth/{sample}.mosdepth.region.dist.txt",
        perbase="results/mosdepth/{sample}.mosdepth.per-base.bed.gz"
    params:
        by="--by 1000",
        prefix=lambda wc: f"results/mosdepth/{wc.sample}"
    threads: 4
    singularity:
        SINGULARITY_IMAGE_MOSDEPTH
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        mosdepth {params.by} --threads {threads} {params.prefix} {input.bam}
        """

rule multiqc:
    """
    Aggregate all mosdepth summaries with MultiQC.
    """
    output:
        html="results/multiQC_mosdepth/multiqc_report.html"
    singularity:
        SINGULARITY_IMAGE_MULTIQC
    shell:
        """
        mkdir -p $(dirname results/multiQC_mosdepth)
        multiqc results/mosdepth --outdir results/multiQC_mosdepth
        """



    # singularity:
    #     SINGULARITY_IMAGE
# singularity run -B /data1/greenbab /data1/greenbab/users/ahunos/apps/containers/mosdepth_multiarch.sif mosdepth /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/somaticVariations/sampleRate4000/results/longphase_haplotag/D-0-1_4000/D-0-1_4000_haplotagged.bam


# singularity run -B /data1/greenbab /data1/greenbab/users/ahunos/apps/containers/multiqc_latest.sif /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/Alignment_metrics/results/mosdepth --outdir results/multiQC_mosdepth
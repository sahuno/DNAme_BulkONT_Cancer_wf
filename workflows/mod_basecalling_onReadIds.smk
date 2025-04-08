# rm -r .snakemake benchmarks results logs
# author: Samuel ahuno
# purpose: modified basecalling from pod5 to sorted bam file that has secondary alignments; Good for multiple files per sample which needs merging at some point 


# run on slurm cluster as a job
## snakemake -s /data1/greenbab/projects/tcr_ont_DNAme/scripts/workflows/mod_basecalling_onReadIds.smk --workflow-profile /data1/greenbab/projects/tcr_ont_DNAme/scripts/configs/slurmMinimal --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos" --keep-going --forceall -np

# snakemake -s /data1/greenbab/projects/tcr_ont_DNAme/scripts/workflows/mod_basecalling_onReadIds.smk --workflow-profile /data1/greenbab/projects/tcr_ont_DNAme/scripts/configs/slurmMinimal --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos" --keep-going -R mark_duplicates -np


import pandas as pd
import os
from collections import defaultdict

# Define parent project directory
parent_dir = "/data1/greenbab/projects/tcr_ont_DNAme"

# Load config file
configfile: parent_dir + "/scripts/configs/config.yaml"
set_species = config["set_species"]

# Load sample sheet
samples_csv = os.path.join(parent_dir, "scripts/configs/samples_pathsTcells.csv")
samples_df = pd.read_csv(samples_csv, sep='\t')

# Construct unique sample IDs
samples_df["sample_id"] = samples_df.apply(
    lambda row: f"{row['patientID']}_{row['sample']}_{row['conditionTumorNormal']}_{row['readId']}", axis=1
)

# Construct group sample ID for merged output
samples_df["merged_sample"] = samples_df.apply(
    lambda row: f"{row['patientID']}_{row['sample']}_{row['conditionTumorNormal']}", axis=1
)

# Mapping sample_id -> pod5 file path
SAMPLE_MAP = dict(zip(samples_df["sample_id"], samples_df["path"]))
SAMPLE_IDS = list(SAMPLE_MAP.keys())

# Mapping merged_sample -> list of sample_ids
merged_groups = defaultdict(list)
for _, row in samples_df.iterrows():
    merged_groups[row["merged_sample"]].append(row["sample_id"])
MERGED_SAMPLE_IDS = list(merged_groups.keys())

rule all:
    input:
        # Sorted individual BAMs and indexes
        expand("results/sort_mod_bases/{sample_id}.bam", sample_id=SAMPLE_IDS),
        expand("results/sort_mod_bases/{sample_id}.bam.csi", sample_id=SAMPLE_IDS),

        # Merged BAMs
        expand("results/merge_bams/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam", merged_sample=MERGED_SAMPLE_IDS),
        expand("results/merge_bams/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam.bai", merged_sample=MERGED_SAMPLE_IDS),
        expand("results/mark_duplicates/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam", merged_sample=MERGED_SAMPLE_IDS),
        expand("results/mark_duplicates/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bai", merged_sample=MERGED_SAMPLE_IDS),
        expand("results/mark_duplicates/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.txt", merged_sample=MERGED_SAMPLE_IDS)

rule mod_bases:
    input:
        lambda wildcards: SAMPLE_MAP[wildcards.sample_id]
    output:
        mod_calls = temp("results/mod_bases/{sample_id}.bam")
    params:
        methyl_context = "5mCG_5hmCG@latest",
        reference_genome = lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        dirDoradoDownloadedModels="/scratch/greenbab/ahunos/DORADO_CACHE_DIR",
        device = "cuda:all",
        modelQuality = "hac"
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    log:
        "logs/mod_bases/{sample_id}.log"
    benchmark:
        "benchmarks/mod_bases/{sample_id}.benchmark.txt"
    shell:
        """
        dorado basecaller {params.modelQuality},{params.methyl_context} \
            --device {params.device} \
            --models-directory {params.dirDoradoDownloadedModels} \
            --reference {params.reference_genome} \
            --verbose --emit-sam --mm2-opts "-Y" \
            {input} \
            > {output.mod_calls} 2> {log}
        """

rule sort_mod_bases:
    input:
        "results/mod_bases/{sample_id}.bam"
    output:
        mod_calls_sorted_bam = "results/sort_mod_bases/{sample_id}.bam",
        mod_calls_sorted_bam_csi = "results/sort_mod_bases/{sample_id}.bam.csi"
    params:
        samtools_threads = 16
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    log:
        "logs/sort_mod_bases/{sample_id}.log"
    benchmark:
        "benchmarks/sort_mod_bases/{sample_id}.benchmark.txt"
    shell:
        """
        samtools sort --threads {params.samtools_threads} \
            -o {output.mod_calls_sorted_bam} -O BAM --write-index {input} 2> {log}
        """

rule merge_bams:
    input:
        lambda wildcards: expand(
            "results/sort_mod_bases/{sample_id}.bam",
            sample_id=merged_groups[wildcards.merged_sample]
        )
    output:
        merged_bam = "results/merge_bams/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam",
        merged_bai = "results/merge_bams/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam.bai"
    params:
        reference_genome = lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"]
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    threads: 8
    log:
        "logs/merge_bams/{merged_sample}.log"
    benchmark:
        "benchmarks/merge_bams/{merged_sample}.benchmark.txt"
    shell:
        """
        samtools merge --threads {threads} -o {output.merged_bam} {input}
        samtools index -@ {threads} {output.merged_bam} 2> {log}
        """

rule mark_duplicates:
    input:
        bams="results/merge_bams/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam"
    output:
         markdup_bam="results/mark_duplicates/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bam",
         markdup_bam_bai="results/mark_duplicates/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.bai",
         markdup_stats="results/mark_duplicates/{merged_sample}/{merged_sample}_modBaseCalls_dedup_sorted.txt"
    log:
        "logs/mark_duplicates/{merged_sample}/{merged_sample}.log"
    singularity: "/data1/greenbab/users/ahunos/apps/containers/gatk.sif"
    benchmark:
        "benchmarks/mark_duplicates/{merged_sample}/{merged_sample}.benchmark.txt"
    shell:
        """ 
        gatk MarkDuplicates --INPUT {input.bams} --OUTPUT {output.markdup_bam} --METRICS_FILE {output.markdup_stats} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT 2> {log}
        """


#        mkdir -p results/merge_bams/{wildcards.merged_sample}

#clean up; to test script again
# rm -rf .snakemake benchmarks results

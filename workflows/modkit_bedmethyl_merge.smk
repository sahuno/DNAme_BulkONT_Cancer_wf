import yaml

# 
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/modkit_bedmethyl_merge.smk --workflow-profile /data1/greenbab/projects/tcr_ont_DNAme/scripts/configs/slurmMinimal --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos" --keep-going --forceall -np


# Load the sample info from samples.yaml
with open("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/dmrSetDB1i_DMSO/samplesByConditions.yaml") as f:
    samples = yaml.safe_load(f)["samples"]

# Create a dictionary mapping each sample group to its output name
sample_groups = {key: f"{key}.bedmethyl_merged.bed" for key in samples.keys()}

rule all:
    input:
        list(sample_groups.values())

rule bedmethyl_merge:
    input:
        lambda wildcards: samples[wildcards.sample]
    output:
        "{sample}.bedmethyl_merged.bed"
    threads: 8
    shell:
        """
        modkit bedmethyl merge {input} \
            -o {output} \
            -g /data1/greenbab/database/mm10/mm10.fa.fai \
            --threads {threads}
        """

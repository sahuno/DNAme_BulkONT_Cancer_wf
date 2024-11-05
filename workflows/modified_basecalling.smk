#author: Samuel ahuno
#purpose: modified basecalling from pod5 to sorted bam file that has secondary alignments 

#notes on selecting models for dorado modified basecalling
# use sup,5mCG_5hmCG for dna_r9.4.1_e8 & 4kHz Sampling rate
# use sup,5mCG_5hmCG for dna_r10.4.1_e8.2_400bps & 4kHz Sampling rate
# use sup,5mC_5hmC@latest, 6mA@latest for dna_r10.4.1_e8.2_400bps & 5kHz Sampling rate


#to run a single sample without snakemake
# singularity exec -B /data1/greenbab --nv /data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif dorado basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --mm2-opts "-Y" | samtools sort --threads 12 -o ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves_singularity_sorted.bam -O BAM --write-index

#To run snakemake workflow on interactive nodes 
# snakemake -s modified_basecalling.smk  --cores all --jobs unlimited --forcerun --printshellcmds --use-singularity --singularity-args "--bind /data1/greenbab" 

#run on slurm cluster as a job
#1. assume singuilarity arguments (bind paths) are in profile
# snakemake -s modified_basecalling.smk --workflow-profile config/cluster_profiles/slurm --jobs 10 --cores all --use-singularity -np

#2. assume paths are not specified in slurm profile, do it manualy
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/modified_basecalling.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab" --keep-going --forceall -np


#change `parent_dir`to the directory you cloned this github to. 
parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_pod5.yaml"
set_species = config["set_species"]


rule all:
    input: 
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls.bam", samples=config["samples"]),
        expand("results/sort_mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam", samples=config["samples"]),
        expand("results/sort_mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam.csi", samples=config["samples"]),

rule mod_bases:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        # mod_calls="results/mod_bases/{samples}/{samples}_modBaseCalls.bam"
         mod_calls=temp("results/mod_bases/{samples}/{samples}_modBaseCalls.bam")
    params:
        methyl_context="5mCG_5hmCG@latest",
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        device = "cuda:all",
        modelQuality="hac"
    # resources:
    #     mem_mb=32000
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    log:
        "logs/mod_bases/{samples}/{samples}.log"
    benchmark:
        "benchmarks/mod_bases/{samples}/{samples}.benchmark.txt"
    shell:
        """ 
        dorado basecaller {params.modelQuality},{params.methyl_context} --device {params.device} {input} --reference {params.reference_genome} --verbose --emit-sam --mm2-opts "-Y" > {output.mod_calls} 2> {log}
        """

rule sort_mod_bases:
    input:
        "results/mod_bases/{samples}/{samples}_modBaseCalls.bam"
    output:
         mod_calls_sorted_bam="results/sort_mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam",
         mod_calls_sorted_bam_csi="results/sort_mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam.csi"
    params:
        samtools_threads=16
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    log:
        "logs/sort_mod_bases/{samples}/{samples}.log"
    benchmark:
        "benchmarks/sort_mod_bases/{samples}/{samples}.benchmark.txt"
    shell:
        """ 
        samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} -O BAM --write-index {input} 2> {log}
        """


#clean up; to test script again
# rm -rf .snakemake benchmarks results

### Test run without snakemake
# singularity exec --no-home -B /data1/greenbab --nv /data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif basecaller sup,5mC_5hmC@latest,6mA@latest /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5 --device "cuda:all" --reference /data1/greenbab/database/mm10/mm10.fa --verbose --emit-sam --mm2-opts "-Y" #| samtools sort --threads 12 -o ONTwf_5mC5hmC6mA_D01_5000kHz_no_moves_singularity_sorted.bam -O BAM --write-index
# singularity exec --no-home -B /data1/greenbab --nv /data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif bash
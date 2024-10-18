configfile: "config/config.yaml"
configfile: "config/samples.yaml"
#$ conda env export -n pod5 > /lila/data/greenbaum/users/ahunos/apps/dorado_ont_wf/config/pod5.yaml


rule all:
    input: 
        expand("results/pod5/{samples}/{samples}.pod5", samples=config["samples"])

rule pod5:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        "results/{rule}/{samples}/{samples}.pod5"
    log:
      "logs/{rule}/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
        """ 
pod5 convert from_fast5 {input} --output {output} --force-overwrite 2> {log}
        """
# rule mod_bases:
#     input:
#         lambda wildcards: config["samples"][wildcards.samples]
#     output:
#     modified_calls_bam="results/{rule}/{samples}/{samples}_5mCG_5hmCG_calls.bam"
#     modified_calls_sorted_bam="results/{rule}/{samples}/{samples}_5mCG_5hmCG_calls_sorted.bam"
#     modified_calls_sorted_bam_bai="results/{rule}/{samples}/{samples}_5mCG_5hmCG_calls_sorted.bam.bai"
#     CpG_5mC_5hmC_merged_bed="results/{rule}/{samples}/{samples}_5mCG_5hmCG_merged.bed"
#     CpG_5mC_bed="results/{rule}/{samples}/{samples}_5mCG.bed"
#     params:
#         dorado_script=config["dorado"],
#         methyl_context="5mCG_5hmCG",
#         basecall_model_file="/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
#         reference_genome=
#     log:
#         "logs/{rule}/{samples}/{samples}.log"
#     shell:
#         """ 
#         source {params.dorado_script} {input} \
#         {params.methyl_context} \
#         {params.basecall_model_file} {params.reference_genome} \
#         {output.modified_calls_bam} {output.modified_calls_sorted_bam} {output.modified_calls_sorted_bam_bai} \
#         {output.CpG_5mC_5hmC_merged_bed} {output.CpG_5mC_bed} 2> {log}
#         """



# this runs ok
# cd snakemake_template
# snakemake -np #test run with
# snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores
# nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &
# snakemake -np
# sh run_snakefile.sh 
# conda install -c conda-forge mamba
#conda env create -f environment.yml
# /opt/common/CentOS_7/anaconda/anaconda3/bin/conda install -n base -c conda-forge mamba
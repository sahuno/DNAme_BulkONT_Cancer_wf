parent_dir = "/home/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
# configfile: "config/samples_sorted_bams.yaml"
configfile: parent_dir + "config/samples_sorted_bams_iris.yaml"

set_species = "mouse"

rule all:
    input: 
        expand("results/modkit/{samples}/{samples}_modBase_summary.log", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_summary.txt", samples=config["samples"])    

rule modkit:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
         summary_log="results/{rule}/{samples}/{samples}_modBase_summary.log",
         summary_txt="results/{rule}/{samples}/{samples}_modBase_summary.txt"
        #  sample_prob_log="results/{rule}/{samples}/{samples}_modBase_sample_prob.log",
        #  sample_prob_tsv="results/{rule}/{samples}/{samples}_probabilities.tsv",
        #  sample_prob_txt="results/{rule}/{samples}/{samples}_probabilities.txt",
        #  sample_prob_thresh_tsv="results/{rule}/{samples}/{samples}_thresholds.tsv",
        #  modpileup_combined="results/{rule}/{samples}/{samples}_modpileup_combined.bed",
        #  modpileup_5mC="results/{rule}/{samples}/{samples}_modpileup_5mC.bed",
        #  modpileup_5mC_log="results/{rule}/{samples}/{samples}_modpileup_5mC.log",
        #  modpileup_combined_log="results/{rule}/{samples}/{samples}_modpileup_combined.log"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=4,
        modkit_prob_percentiles=0.1,
        outdir= "results/{rule}/{samples}/"
    log:
        "logs/{rule}/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
        """ 
echo ""modkit pileup...""
modkit summary --threads 4 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt}
        """

# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_modkit_test.smk --cores 12 --forcerun --use-conda  -np #dry run with cores

parent_dir = "/home/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_sorted_bams_spect_lilac.yaml"

# configfile: parent_dir + "config/samples_sorted_bams_iris.yaml"
#/home/ahunos/apps/dorado_ont_wf/config/samples_sorted_bams_iris.yaml

set_species = "human"

rule all:
    input: 
        expand("results/modkit/{samples}/{samples}_modBase_summary.log", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_summary.txt", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modpileup_combined.bed", samples=config["samples"]),    
        expand("results/modkit/{samples}/{samples}_modpileup_5mC.bed", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modpileup_5mC.log", samples=config["samples"]),   
        expand("results/modkit/{samples}/{samples}_modpileup_combined.log", samples=config["samples"]),
        expand("results/modkit/{samples}/{samples}_modBase_sample_prob.log", samples=config["samples"]), 
        expand("results/modkit/{samples}/{samples}_probabilities.tsv", samples=config["samples"]), 
        expand("results/modkit/{samples}/{samples}_probabilities.txt", samples=config["samples"]), 
        expand("results/modkit/{samples}/{samples}_thresholds.tsv", samples=config["samples"]) 
 
rule modkit:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
         summary_log="results/{rule}/{samples}/{samples}_modBase_summary.log",
         summary_txt="results/{rule}/{samples}/{samples}_modBase_summary.txt",
         sample_prob_log="results/{rule}/{samples}/{samples}_modBase_sample_prob.log",
         sample_prob_tsv="results/{rule}/{samples}/{samples}_probabilities.tsv",
         sample_prob_txt="results/{rule}/{samples}/{samples}_probabilities.txt",
         sample_prob_thresh_tsv="results/{rule}/{samples}/{samples}_thresholds.tsv",
         modpileup_combined="results/{rule}/{samples}/{samples}_modpileup_combined.bed",
         modpileup_5mC="results/{rule}/{samples}/{samples}_modpileup_5mC.bed",
         modpileup_5mC_log="results/{rule}/{samples}/{samples}_modpileup_5mC.log",
         modpileup_combined_log="results/{rule}/{samples}/{samples}_modpileup_combined.log"

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
modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt} 2> {log}

echo ""modkit pileup...""
# its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log} --only-tabs 2> {log}
echo ""done modkit pileup...""

echo ""modkit combine 5mc 5hmc pileup...""
# its either C or Methylated(i don't care whether it's 5mC or 5hmC)
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log} --only-tabs 2> {log}
echo ""done modkit combine 5mc 5hmc pileup...""


echo ""find the prob filtering threshold...""
modkit sample-probs --threads {params.modkit_threads} \
--percentiles {params.modkit_prob_percentiles} \
--out-dir {params.outdir} \
--prefix {wildcards.samples} --hist \
--log-filepath {output.sample_prob_log} \
{input} 2> {log}

###########5mc & 5hmc, bedgraph outputs ############
modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
--prefix {wildcards.samples} --cpg --combine-mods --ref {params.reference_genome}

#these files will be generated
# prefix_C_CG0_positive.bedgraph
# prefix_C_CG0_negative.bedgraph

###########5mc & 5hmc, bedgraph outputs############
modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
--prefix {wildcards.samples} --cpg --ref {params.reference_genome}

#these files will be generated
# prefix_h_CG0_negative.bedgraph
# prefix_m_CG0_negative.bedgraph
# prefix_h_CG0_positive.bedgraph
# prefix_m_CG0_positive.bedgraph

echo ""done modkit sample-probs...""
        """

# snakemake -s /home/ahunos/apps/dorado_ont_wf/methylation_calling_modkit.smk --cores 12 --forcerun --use-conda  -np #dry run with cores

#run with profile file 
# snakemake -s /home/ahunos/apps/dorado_ont_wf/methylation_calling_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall -np


#run with slurm
# snakemake -s /home/ahunos/apps/dorado_ont_wf/methylation_calling_modkit.smk  --latency-wait 60 --restart-times 2 --keep-going --forceall --use-conda \
# --cluster-config /home/ahunos/apps/dorado_ont_wf/config/cluster_slurm.yaml \
# --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -n {cluster.tasks}" --jobs 10 --cores all
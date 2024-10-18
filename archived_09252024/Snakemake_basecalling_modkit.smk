parent_dir = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
# configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_spectrum009N.yaml"
configfile: "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_split_sample_rates_pod5_triEpi.yaml"

# configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_pod5_dev.yaml"
# configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_pod5_dev.yaml"


#To run with profile file 
#snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs unlimited --cores all 
# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall --latency-wait 60 --restart-times 2 -np


set_species = "mouse"
print(config["samples"])
#less /home/ahunos/apps/dorado_ont_wf/config/samples_pod5_spectrum.yaml

rule all:
    input: 
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam", samples=config["samples"]),
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam.csi", samples=config["samples"]),
        expand("results/dorado_summary/{samples}/{samples}_seq_summary.txt", samples=config["samples"]),
        expand("results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam", samples=config["samples"]),
        expand("results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bai", samples=config["samples"]),
        expand("results/mark_duplicates/{samples}/{samples}_marked_dup_metrics.txt", samples=config["samples"]),
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

rule mod_bases:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
         mod_calls_sorted_bam="results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam",
         mod_calls_sorted_bam_csi="results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam.csi"
    params:
        methyl_context="5mCG_5hmCG",
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        samtools_threads=16,
        device = "cuda:all"
    log:
        "logs/mod_bases/{samples}/{samples}.log"
    benchmark:
        "benchmarks/mod_bases/{samples}/{samples}.benchmark.txt"
    # resources:
    #     nvidia_gpu=1
    shell:
        """ 
        dorado basecaller sup,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log} || dorado basecaller hac,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log}
        """

rule dorado_summary:
    input:
        "results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam"
    output:
         seq_summary="results/dorado_summary/{samples}/{samples}_seq_summary.txt"
    log:
        "logs/dorado_summary/{samples}/{samples}.log"
    benchmark:
        "benchmarks/dorado_summary/{samples}/{samples}.benchmark.txt"
    shell:
        """ 
        dorado summary {input} --verbose > {output.seq_summary} 2> {log}
        """

        
rule mark_duplicates:
    input:
        bams="results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam"
    output:
         markdup_bam="results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam",
         markdup_bam_bai="results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bai",
         markdup_stats="results/mark_duplicates/{samples}/{samples}_marked_dup_metrics.txt"
    log:
        "logs/mark_duplicates/{samples}/{samples}.log"
    benchmark:
        "benchmarks/mark_duplicates/{samples}/{samples}.benchmark.txt"
    shell:
        """ 
        module load java
        gatk MarkDuplicates --INPUT {input.bams} --OUTPUT {output.markdup_bam} --METRICS_FILE {output.markdup_stats} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT 2> {log}
        """

rule modkit:
    input:
        "results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam"
    output:
         summary_log="results/modkit/{samples}/{samples}_modBase_summary.log",
         summary_txt="results/modkit/{samples}/{samples}_modBase_summary.txt",
         sample_prob_log="results/modkit/{samples}/{samples}_modBase_sample_prob.log",
         sample_prob_tsv="results/modkit/{samples}/{samples}_probabilities.tsv",
         sample_prob_txt="results/modkit/{samples}/{samples}_probabilities.txt",
         sample_prob_thresh_tsv="results/modkit/{samples}/{samples}_thresholds.tsv",
         modpileup_combined="results/modkit/{samples}/{samples}_modpileup_combined.bed",
         modpileup_5mC="results/modkit/{samples}/{samples}_modpileup_5mC.bed",
         modpileup_5mC_log="results/modkit/{samples}/{samples}_modpileup_5mC.log",
         modpileup_combined_log="results/modkit/{samples}/{samples}_modpileup_combined.log"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=4,
        modkit_prob_percentiles=0.1,
        outdir= "results/modkit/{samples}/"
    log:
        "logs/modkit/{samples}/{samples}.log"
    benchmark:
        "benchmarks/modkit/{samples}/{samples}.benchmark.txt"
    conda: config["pod5_env"]
    shell:
        """ 
################################################        
echo ""modkit pileup...""
modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt} 2> {log}

################################################
echo ""modkit pileup...""
# its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log} --only-tabs 2> {log}
echo ""done modkit pileup...""

################################################
echo ""modkit combine 5mc 5hmc pileup...""
# its either C or Methylated(i don't care whether it's 5mC or 5hmC)
modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log} --only-tabs 2> {log}
echo ""done modkit combine 5mc 5hmc pileup...""

################################################
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




# echo ""indexing..""
# samtools index -@ {params.samtools_threads} {output.dedup_mod_calls_sorted_bam}

# echo ""delete tmp files..""
# rm {output.mod_calls_bam} {output.mod_calls_sorted_bam} {{output.mod_calls_sorted_bam.bai}}









#run with profile file 
# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall -np

# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --executor slurm --default-resources slurm_partition=componc_gpu slurm_account=greenbab runtime=600 slurm_extra="'--gres gpu:4'" mem_mb_per_cpu=24800 cpus_per_task=24 --jobs 10 --cores all --keep-going --forceall --latency-wait 60 --restart-times 2 -np

###run with slurm
# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk  --latency-wait 60 --restart-times 2 --keep-going --forceall --use-conda \
# --cluster-config /home/ahunos/apps/dorado_ont_wf/config/cluster_slurm.yaml \
# --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -n {cluster.tasks} --gpus 4" --jobs 10 --cores all




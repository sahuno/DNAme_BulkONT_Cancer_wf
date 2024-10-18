parent_dir = "/home/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_pod5_dev.yaml"

set_species = "mouse"

rule all:
    input: 
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam", samples=config["samples"])
        # expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam.csi", samples=config["samples"])
        # expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted_dedup.bam", samples=config["samples"]),
        # expand("results/mod_bases/{samples}/{samples}_seq_summary.txt", samples=config["samples"])


rule mod_bases:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
         mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam",
        #  mod_calls_sorted_bam_csi="results/{rule}/{samples}/{samples}_modBaseCalls_sorted.bam.csi",
        #  dedup_mod_calls_sorted_bam="results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam",
        #  dedup_mod_calls_sorted_bam_bai="results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam.bai",
        #  seq_summary="results/{rule}/{samples}/{samples}_seq_summary.txt"
    params:
        methyl_context="5mCG_5hmCG",
        # basecall_model_file="/lila/data/greenbaum/users/ahunos/refs/dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        samtools_threads=32,
        device = "cuda:all"
    log:
        "logs/{rule}/{samples}/{samples}.log"
    shell:
        """ 
        # Try the first basecalling model
        dorado basecaller sup,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose > {output.mod_calls_sorted_bam} 2> {log}      
        """

# dorado basecaller sup,5mCG_5hmCG@latest {input} --device \"cuda:all\" --emit-sam --reference {params.reference_genome} --verbose | samtools sort -o {output.mod_calls_sorted_bam} --write-index 2> {log}

# dorado basecaller sup,{params.methyl_context}@latest {input} --device \"cuda:all\" --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index 2> {log}

# echo ""seq summary..""
# dorado summary {output.mod_calls_sorted_bam} > {output.seq_summary}
# echo ""marking duplictes..""
# gatk MarkDuplicates \
# --INPUT {output.mod_calls_sorted_bam} \
# --OUTPUT {output.dedup_mod_calls_sorted_bam} \
# --METRICS_FILE marked_dup_metrics.txt

# echo ""indexing..""
# samtools index -@ {params.samtools_threads} {output.dedup_mod_calls_sorted_bam}

# echo ""delete tmp files..""
# rm {output.mod_calls_bam} {output.mod_calls_sorted_bam} {{output.mod_calls_sorted_bam.bai}}







# dorado basecaller hac,5mCG_5hmCG /data1/greenbab/projects/methyl_benchmark_spectrum/data/raw/pod5/044N_v14/pod5/044N_v14.pod5 --device "cuda:all" --emit-sam > 044N_v14.bam

# dorado basecaller hac,5mCG_5hmCG ~/apps/sandbox/modified_basecalling_mninimal_test/044T_v14_minimal.pod5 --device "cuda:all" --emit-sam --reference /data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta -v > 044T_v14_test2.bam

##works well
# # dorado basecaller sup,5mCG_5hmCG@latest ~/apps/sandbox/modified_basecalling_mninimal_test/044T_v14_minimal.pod5 --device "cuda:all" --emit-sam --reference /data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta -v > 044T_v14_sup5mCG_5hmCG.bam | 
# samtools sort 044T_v14_sup5mCG_5hmCG.bam -o 044T_v14_sup5mCG_5hmCG_sorted.bam --write-index
# modkit pileup --threads 4 044T_v14_sup5mCG_5hmCG_sorted.bam 044T_v14_sup5mCG_5hmCG.bed --cpg --ref /data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta --ignore h --log-filepath 044T_v14_sup5mCG_5hmCG.log --only-tabs
# #test basecalling and sorting
# dorado basecaller sup,5mCG_5hmCG@latest ~/apps/sandbox/modified_basecalling_mninimal_test/044T_v14_minimal.pod5 --device "cuda:all" --emit-sam --reference /data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta -v | samtools sort -o 044T_v14_sup5mCG_5hmCG_sorted_combined.bam --write-index







# rule modkit:
#     input:
#         "results/{rule}/{samples}/{samples}_modBaseCalls_sorted_dedup.bam",
#           sample_name=lambda wildcards: wildcards.samples
#     output:
#          summary_log="results/{rule}/{samples}/{samples}_modBase_summary.log",
#          summary_txt="results/{rule}/{samples}/{samples}_modBase_summary.txt",
#          sample_prob_log="results/{rule}/{samples}/{samples}_modBase_sample_prob.log",
#          sample_prob_tsv="results/{rule}/{samples}/{samples}_probabilities.tsv",
#          sample_prob_txt="results/{rule}/{samples}/{samples}_probabilities.txt",
#          sample_prob_thresh_tsv="results/{rule}/{samples}/{samples}_thresholds.tsv",
#          modpileup_combined="results/{rule}/{samples}/{samples}_modpileup_combined.bed",
#          modpileup_5mC="results/{rule}/{samples}/{samples}_modpileup_5mC.bed",
#          modpileup_5mC_log="results/{rule}/{samples}/{samples}_modpileup_5mC.log",
#          modpileup_combined_log="results/{rule}/{samples}/{samples}_modpileup_combined.log"

#     params:
#           reference_genome=config["ref_genome"],
#           modkit_threads=32,
#           modkit_prob_percentiles=0.1,
#           outdir= "results/{rule}/{samples}/"
#     log:
#         "logs/{rule}/{samples}/{samples}.log"
#     shell:
#         """ 
# echo ""modkit summary...""
# modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt}

# echo ""modkit pileup 5mc...""
# # its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
# modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log}

# echo ""modkit pileup 5mc n 5hmC...""
# # its either C or Methylated(i don't care whether it's 5mC or 5hmC)
# modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log}


# echo ""find the prob filtering threshold...""
# modkit sample-probs --threads {params.modkit_threads} \
# --percentiles {params.modkit_prob_percentiles} \
# --out-dir {params.outdir} \
# --prefix {input.sample_name} --hist \
# --log-filepath {output.sample_prob_log} \
# {input}
#         """



###run with slurm
# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk  --latency-wait 60 --restart-times 2 --keep-going --forceall --use-conda \
# --cluster-config /home/ahunos/apps/dorado_ont_wf/config/cluster_slurm.yaml \
# --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -n {cluster.tasks} --gpus 4" --jobs 10 --cores all

# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --workflow-profile /home/ahunos/apps/dorado_ont_wf/config/profile_slurm.yaml --jobs 10 --cores all

# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit_dev.smk --executor slurm --default-resources slurm_partition=componc_gpu slurm_account=greenbab runtime=360 slurm_extra="'--gres gpu:4'" mem_mb_per_cpu=24800 cpus_per_task=24 --jobs 10 --cores all --keep-going --forceall --latency-wait 60 --restart-times 2


# pod5 view /data1/greenbab/projects/methyl_benchmark_spectrum/data/raw/pod5/009T1/results/pod5/009T1/009T1.pod5

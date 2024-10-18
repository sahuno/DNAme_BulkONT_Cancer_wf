#purpose: sample sequenced with differeent flowcells need basecalling done separately and merged. then duplicates marked and modified bases called
import yaml
# sample_yaml_path = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_spectrum009N.yaml" #real samples
sample_yaml_path = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_mergeBam_modkit.yaml" #real
# sample_yaml_path = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_multiple_flowcells.yaml"

parent_dir = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: sample_yaml_path

# configfile: parent_dir + "config/samples_spectrum009N.yaml"
set_species = "human"



# Load samples from YAML file
with open(sample_yaml_path, "r") as f:
    samplesYaml = yaml.safe_load(f)["samples"]

modkit_extensions = {"modkit":["_modBase_summary.log", "_modBase_summary.txt", "_modBase_sample_prob.log","_probabilities.tsv", "_probabilities.txt", "_thresholds.tsv", "_modpileup_combined.bed", "_modpileup_5mC.bed", "_modpileup_5mC.log", "_modpileup_combined.log"]}

#make list to store output files
mod_bases_output_files = []
dorado_summary_output_files = []
merge_bams_output_files = []
mark_duplicates_output_files = []
modkit_output_files = []

mod_bases_logs = []
dorado_summary_logs = []
merge_bams_logs = []
mark_duplicates_logs = []
modkit_logs = []

for sample, files in samplesYaml.items():
    # print(idx, sample, files)
    for i,file in enumerate(files):
        # print(f"results/mod_bases/{sample}/{sample}.{i}.bam")
        # print(f"results/mod_bases/{sample}/{sample}.{i}.bam")
        mod_bases_output_files.append(f"results/mod_bases/{sample}/{sample}.{i}.bam")
        # print(f"results/mod_bases/{sample}/{sample}.{i}.bam.csi")
        mod_bases_logs.append(f"logs/mod_bases/{sample}/{sample}.{i}.log")
        mod_bases_output_files.append(f"results/mod_bases/{sample}/{sample}.{i}.bam.bai")
        dorado_summary_output_files.append(f"results/dorado_summary/{sample}/{sample}.{i}_seq_summary.txt")
        dorado_summary_logs.append(f"logs/dorado_summary/{sample}/{sample}.{i}.log")
    merge_bams_output_files.append(f"results/merge_bams/{sample}/{sample}_modBaseCalls_merged_sorted.bam")
    merge_bams_output_files.append(f"results/merge_bams/{sample}/{sample}_modBaseCalls_merged_sorted.bam.bai")
    merge_bams_logs.append(f"logs/merge_bams/{sample}/{sample}.log")  
    mark_duplicates_output_files.append(f"results/mark_duplicates/{sample}/{sample}_modBaseCalls_merged_sorted_dedup.bam")
    mark_duplicates_output_files.append(f"results/mark_duplicates/{sample}/{sample}_modBaseCalls_merged_sorted_dedup.bai")
    mark_duplicates_output_files.append(f"results/mark_duplicates/{sample}/{sample}_modBaseCalls_merged_sorted_dedup_metrics.txt")
    mark_duplicates_logs.append(f"logs/mark_duplicates/{sample}/{sample}.log")
    for ext in modkit_extensions["modkit"]:
        modkit_output_files.append(f"results/modkit/{sample}/{sample}{ext}")


print(f"config['samples']\n {config["samples"]}")
print(f"\nmod_bases_output_files\n {mod_bases_output_files}")
print(f"\nsamplesKeys\n{samplesYaml.keys()}")



# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/basecalling_mergeBams_mdup_modkit.smk --cores 5 -np

rule all:
    input:
        mod_bases_output_files
        # dorado_summary_output_files,
        # merge_bams_output_files
        # mark_duplicates_output_files
        # modkit_output_files

rule mod_bases:
    input:
        lambda wildcards: expand("{{sample_path}}", sample_path=config['samples'][wildcards.sample])
    output:
         mod_calls_sorted_bam = [file for file in mod_bases_output_files if file.endswith('.bam')],
         mod_calls_sorted_bam_csi = [file for file in mod_bases_output_files if file.endswith('.csi')]
    params:
        methyl_context="5mCG_5hmCG",
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        samtools_threads=16,
        device = "cuda:all"
    log:
        mod_bases_logs
    # benchmark:
    #     "benchmarks/mod_bases/{wildcards.samples}/{wildcards.samples}.benchmark.txt"
    shell:
        """ 
        dorado basecaller sup,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log} || dorado basecaller hac,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log}
        """

# rule dorado_summary:
#     input:
#         [file for file in mod_bases_output_files if file.endswith('.bam')]
#     output:
#          seq_summary=dorado_summary_output_files
#     log:
#         dorado_summary_logs
#     # benchmark:
#         # "benchmarks/dorado_summary/{samples}/{samples}.benchmark.txt"
#     shell:
#         """ 
#         dorado summary {input} --verbose > {output.seq_summary} 2> {log}
#         """

# rule merge_bams:
#     input:
#         bamsAll=[file for file in mod_bases_output_files if file.endswith('.bam')]
#         # bam1=lambda wildcards: config["samples"][wildcards.samples][0]
#     output:
#          merged_bam=[file for file in merge_bams_output_files if file.endswith('.bam')],
#          merged_bai=[file for file in merge_bams_output_files if file.endswith('.bai')]
#     params:
#         reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
#         samtools_threads=16
#     log:
#         merge_bams_logs
#     # benchmark:
#         # "benchmarks/merge_bams/{samples}/{samples}.benchmark.txt"
#     shell:
#         """ 
#         echo "Merging bams for {wildcards.samples}"
#         samtools merge --threads {params.samtools_threads} -o {output.merged_bam} {input.bamsAll} && samtools index -@ 32 {output.merged_bam} 2> {log}
#         """
# #SILENT; [file for file in listfiles if file.endswith('.bam')]
# #echo "first bam file {input.bam1} of sample {wildcards.samples}"
# rule mark_duplicates:
#     input:
#         bams=[file for file in merge_bams_output_files if file.endswith('.bam')]
#     output:
#          markdup_bam=[file for file in mark_duplicates_output_files if file.endswith('.bam')],
#          markdup_bam_bai=[file for file in mark_duplicates_output_files if file.endswith('.bai')],
#          markdup_stats= [file for file in mark_duplicates_output_files if file.endswith('.txt')]
#     log:
#         mark_duplicates_logs
#     # benchmark:
#     #     "benchmarks/mark_duplicates/{samples}/{samples}.benchmark.txt"
#     shell:
#         """ 
#         module load java
#         gatk MarkDuplicates --INPUT {input.bams} --OUTPUT {output.markdup_bam} --METRICS_FILE {output.markdup_stats} --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT 2> {log}
#         """

# rule modkit:
#     input:
#         "results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam"
#     output:
#          summary_log="results/modkit/{samples}/{samples}_modBase_summary.log",
#          summary_txt="results/modkit/{samples}/{samples}_modBase_summary.txt",
#          sample_prob_log="results/modkit/{samples}/{samples}_modBase_sample_prob.log",
#          sample_prob_tsv="results/modkit/{samples}/{samples}_probabilities.tsv",
#          sample_prob_txt="results/modkit/{samples}/{samples}_probabilities.txt",
#          sample_prob_thresh_tsv="results/modkit/{samples}/{samples}_thresholds.tsv",
#          modpileup_combined="results/modkit/{samples}/{samples}_modpileup_combined.bed",
#          modpileup_5mC="results/modkit/{samples}/{samples}_modpileup_5mC.bed",
#          modpileup_5mC_log="results/modkit/{samples}/{samples}_modpileup_5mC.log",
#          modpileup_combined_log="results/modkit/{samples}/{samples}_modpileup_combined.log"
#     params:
#         reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
#         modkit_threads=4,
#         modkit_prob_percentiles=0.1,
#         outdir= "results/modkit/{samples}/"
#     log:
#         "logs/modkit/{samples}/{samples}.log"
#     benchmark:
#         "benchmarks/modkit/{samples}/{samples}.benchmark.txt"
#     conda: config["pod5_env"]
#     shell:
#         """ 
# ################################################        
# echo ""modkit pileup...""
# modkit summary --threads 12 --only-mapped {input} --log-filepath  {output.summary_log} > {output.summary_txt} 2> {log}

# ################################################
# echo ""modkit pileup...""
# # its either C or 5mC(its either unmethylated or 5mC) but what happens to 5hmC when their prob is redistributed? 
# modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_5mC} --cpg --ref {params.reference_genome} --ignore h --log-filepath {output.modpileup_5mC_log} --only-tabs 2> {log}
# echo ""done modkit pileup...""

# ################################################
# echo ""modkit combine 5mc 5hmc pileup...""
# # its either C or Methylated(i don't care whether it's 5mC or 5hmC)
# modkit pileup --threads {params.modkit_threads} {input} {output.modpileup_combined} --cpg --ref {params.reference_genome} --combine-mods --log-filepath {output.modpileup_combined_log} --only-tabs 2> {log}
# echo ""done modkit combine 5mc 5hmc pileup...""

# ################################################
# echo ""find the prob filtering threshold...""
# modkit sample-probs --threads {params.modkit_threads} \
# --percentiles {params.modkit_prob_percentiles} \
# --out-dir {params.outdir} \
# --prefix {wildcards.samples} --hist \
# --log-filepath {output.sample_prob_log} \
# {input} 2> {log}

# ###########5mc & 5hmc, bedgraph outputs ############
# modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
# --prefix {wildcards.samples} --cpg --combine-mods --ref {params.reference_genome}

# #these files will be generated
# # prefix_C_CG0_positive.bedgraph
# # prefix_C_CG0_negative.bedgraph

# ###########5mc & 5hmc, bedgraph outputs############
# modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
# --prefix {wildcards.samples} --cpg --ref {params.reference_genome}

# #these files will be generated
# # prefix_h_CG0_negative.bedgraph
# # prefix_m_CG0_negative.bedgraph
# # prefix_h_CG0_positive.bedgraph
# # prefix_m_CG0_positive.bedgraph

# echo ""done modkit sample-probs...""
#         """


#######
## run code
# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/basecalling_mergeBams_mdup_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall -np
######


####################################################################################
## delete later
####################################################################################
# # Flatten the list of input files and create corresponding output files
# input_files = []
# output_files = []
# rules = ["mod_bases", "dorado_summary", "merge_bams", "mark_duplicates", "modkit"]

# for idx, (sample, files) in enumerate(samples.items()):
#     for file in files:
#       for rule in rules:
#         print(f"results/{rule}/{sample}/{idx}.ext")
#         # input_files.append(file)
#         # output_files.append(file.replace(".pod5", ".out"))  # Adjust extension as needed
######################################################################################
######################################################################################

import yaml
# Load sample data from YAML configuration file
configfile: "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_multiple_flowcells_numbered.yaml"
with open(configfile, "r") as f:
    samples = yaml.safe_load(f)["samples"]

# Rule to generate truncated files
rule all:
    input:
        expand("{sample}_{replicate}_short.pod5", 
               sample=samples.keys(), 
               replicate=samples.values())

rule truncate:
    input:
        "data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/{sample}/{replicate}.pod5"
    output:
        "{sample}_{replicate}_short.pod5"
    shell:
        "samtools -n 20 {input} > {output}"




# """ 
# head -n 20 {input} > {output}
# """

#cat /data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_multiple_flowcells.yamlbash:islogin01:/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config 1029 $ 

        #  markdup_bam="results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bam",
        #  markdup_bam_bai="results/mark_duplicates/{samples}/{samples}_modBaseCalls_sorted_dup.bai",
        #  markdup_stats="results/mark_duplicates/{samples}/{samples}_marked_dup_metrics.txt"

# rule mod_bases:
#     input:
#         lambda wildcards: config["samples"][wildcards.samples]
#     output:

#          mod_calls_sorted_bam = [file for file in mod_bases_output_files if file.endswith('.bam')],
#          mod_calls_sorted_bam_csi = [file for file in mod_bases_output_files if file.endswith('.csi')]
#     params:
#         methyl_context="5mCG_5hmCG",
#         reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
#         samtools_threads=16,
#         device = "cuda:all"
#     log:
#         mod_bases_logs
#     # benchmark:
#     #     "benchmarks/mod_bases/{wildcards.samples}/{wildcards.samples}.benchmark.txt"
#     shell:
#         """ 
#         dorado basecaller sup,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log} || dorado basecaller hac,5mCG_5hmCG@latest {input} --device {params.device} --emit-sam --reference {params.reference_genome} --verbose | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} --write-index  2> {log}
#         """

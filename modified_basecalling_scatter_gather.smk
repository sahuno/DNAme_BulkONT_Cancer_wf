parent_dir = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_split_sample_rates_pod5_triEpi.yaml"

#input sampleName/pod5
# > pod5 repack pod5s/*.pod5 repacked_pods/
# pod5 repack /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5subset/split_by_sample_rate/D-0-1/sample_rate-5000.pod5/*.pod5 /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/raw/TRI_DIVYA_POD5/split_by_sampleRates/results/pod5Repacked/D-0-1/sampleRate5000/

#To run with profile file 
#snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs unlimited --cores all 
# snakemake -s /data1/greenbab/users/ahunos/apps/dorado_ont_wf/Snakemake_basecalling_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall --latency-wait 60 --restart-times 2 -np


set_species = "mouse"
print(config["samples"])
#less /home/ahunos/apps/dorado_ont_wf/config/samples_pod5_spectrum.yaml

rule all:
    input: 
        expand("results/mod_bases/{samples}/{readName}_modBaseCalls_sorted.bam", samples=config["samples"]),
        expand("results/mod_bases/{samples}/{samples}_modBaseCalls_sorted.bam.csi", samples=config["samples"]),

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
    shell:
        """ 
        dorado basecaller sup,5mC_5hmC@latest,6mA@latest {input} --device {params.device} --reference {params.reference_genome} --verbose --emit-sam --mm2-opts "-Y" | samtools sort --threads {params.samtools_threads} -o {output.mod_calls_sorted_bam} -O BAM --write-index  2> {log}
        """

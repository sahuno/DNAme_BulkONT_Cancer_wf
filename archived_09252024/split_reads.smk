#i have 1 big pod5 with numerous reads with sampling rates. this messes up base calling. i want to split the pod5 into multiple pod5s based on the sampling rate.
configfile: "/home/ahunos/apps/dorado_ont_wf/config/config.yaml"
configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_pod5_spectrum.yaml"

rule all:
    input: 
        expand("results/pod5view/{samples}/{samples}_sampleRateSummary.tsv", samples=config["samples"]),
        expand("results/pod5subset/split_by_sample_rate/{samples}/{samples}.done.txt", samples=config["samples"])
rule pod5view:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        out_pod5View="results/{rule}/{samples}/{samples}_sampleRateSummary.tsv"
    log:
      "logs/{rule}/{samples}/{samples}.log"
    conda: config["pod5_env"]
    # params:
        # sample_rate="sample_rate"
    shell:
        "pod5 view -t 4 {input} --include \"read_id, sample_rate\" --output {output.out_pod5View} 2> {log}"


rule pod5subset:
    input:
        pod5 = lambda wildcards: config["samples"][wildcards.samples],
        pod5View="results/pod5view/{samples}/{samples}_sampleRateSummary.tsv"
    output:
        "results/pod5subset/split_by_sample_rate/{samples}/{samples}.done.txt"
    # log:
    #   "logs/{rule}/split_by_sample_rate/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
        "pod5 subset -t 4 --summary {input.pod5View} --columns sample_rate --output results/pod5subset/split_by_sample_rate/{wildcards.samples} {input.pod5} && touch {output}"


# snakemake -s /home/ahunos/apps/dorado_ont_wf/split_reads.smk --cores 12 --forcerun --use-conda  -np #dry run with cores


#run with slurm
# snakemake -s /home/ahunos/apps/dorado_ont_wf/split_reads.smk  --latency-wait 60 --restart-times 2 --keep-going --forceall --use-conda \
# --cluster-config /home/ahunos/apps/dorado_ont_wf/config/cluster_slurm.yaml \
# --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -n {cluster.tasks}" --jobs 10 --cores all
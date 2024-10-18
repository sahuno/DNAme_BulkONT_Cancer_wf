# include: "/home/ahunos/apps/dorado_ont_wf/Snakefile_pod5inspect.smk"
configfile: "/home/ahunos/apps/dorado_ont_wf/config/config.yaml"
configfile: "/home/ahunos/apps/dorado_ont_wf/config/samples_spectrum009N.yaml"
#$ conda env export -n pod5 > /lila/data/greenbaum/users/ahunos/apps/dorado_ont_wf/config/pod5.yaml
#/home/ahunos/apps/dorado_ont_wf/config/samples_fast5_2_pod5.yaml

rule all:
    input: 
        expand("results/pod5/{samples}/{samples}.pod5", samples=config["samples"]),
        expand("results/pod5stats/{samples}/{samples}_pod5_stats.tsv", samples=config["samples"]),
        expand("results/pod5_inspect_debug/{samples}/{samples}_debug_pod5.txt", samples=config["samples"])

rule pod5:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        out_pod5="results/{rule}/{samples}/{samples}.pod5"
    log:
      "logs/{rule}/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
       "POD5_DEBUG=1 pod5 convert fast5 {input} --threads 5 --output {output.out_pod5} --strict --force-overwrite 2> {log}"

rule pod5stats:
    input:
        "results/pod5/{samples}/{samples}.pod5"
    output:
        out_stats="results/{rule}/{samples}/{samples}_pod5_stats.tsv"
    log:
      "logs/{rule}/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
       "POD5_DEBUG=1 pod5 view -t 4 {input} --output {output.out_stats} --force-overwrite 2> {log}"

rule pod5_inspect_debug:
    input:
        "results/pod5/{samples}/{samples}.pod5"
    output:
        "results/pod5_inspect_debug/{samples}/{samples}_debug_pod5.txt"
    log:
      "logs/pod5_inspect_debug/{samples}/{samples}.log"
    conda: config["pod5_env"]
    shell:
       "pod5 inspect debug {input} > {output} 2> {log}"

### add pod5 inspect
#pod5 inspect debug D-S-1.pod5 > D-S-1_pod5_debug.txt



## snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakefile_fast5_to_pod5.smk --workflow-profile /data1/greenbab/users/ahunos/apps/configs/snakemake/slurm --jobs 10 --cores all --use-conda --keep-going --forceall -np --quiet

# this runs ok
# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakefile_fast5_to_pod5.smk --cores 12 --forcerun --use-conda  -np #dry run with cores
#--use-conda 

# snakemake -s /home/ahunos/apps/dorado_ont_wf/Snakefile_fast5_to_pod5.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --use-conda \
# --cluster-config /home/ahunos/apps/dorado_ont_wf/config/cluster_slurm.yaml \
# --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -n {cluster.tasks}" --jobs 10 --cores all


# nohup snakemake -s Snakefile_fast5_pod5.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &
# sh run_snakefile_pod5.sh 

# snakemake -s Snakefile_fast5_pod5.smk -np --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 5:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3
#conda env create -f environment.yml
# /opt/common/CentOS_7/anaconda/anaconda3/bin/conda install -n base -c conda-forge mamba
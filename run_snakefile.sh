#!/bin/bash
#author: samuel ahuno
# purpose: run commands for differnt use cases


#TODO : pull run commands from snakmake files

####################### Archived ####################
# Run snakemake
# snakemake --jobname 's.{jobid}.{rulename}' \
# 	--snakefile $1 \
# 	--use-conda \
# 	--keep-going \
# 	--reason \
# 	--printshellcmds \
# 	--latency-wait 10 \
# 	--rerun-incomplete \
# 	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
# 	-j 5000 \
# 	--cluster-config $2 \
# 	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o out.txt -e err.txt"




# bash /home/ahunos/apps/dorado_ont_wf/run_Snakefile_fast5_pod5.sh /home/ahunos/apps/dorado_ont_wf/Snakefile_fast5_pod5.smk /home/ahunos/apps/dorado_ont_wf/config/cluster_slurm.yaml
# Run snakemake
#argument 1; snakemake file
#argument 2; cluster config file
# snakemake --jobname 's.{jobid}.{rulename}' \
# 	--snakefile $1 \
# 	--use-conda \
# 	--keep-going \
# 	--reason \
# 	--printshellcmds \
# 	--latency-wait 10 \
# 	--rerun-incomplete \
# 	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
# 	-j 5000 \
# 	--cluster-config $2 \
# 	--cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -n {cluster.tasks}"
####################### Archived ####################

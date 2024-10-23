# ONT dorado + modkit snakemake workflow
## date - Sept 3rd 2024
## author - samuel ahuno, mskcc
# DNA Methylation of cancer genomes pipeline Oxford Nanopore

## how to run snakemake workflow
-Add notes


# snakefiles
separate basecalling & fucntional anlysis.
SNV should use the same already workflow as caner genomics 
but improve SNV-based pahsing with methylation 

1. minimal haplotype resolved somatic DNAme + germline DNAme (input:pod5 per file output:phased .modbam)
2. minimal workflow for functional analysis + visualization (modkit dmr reegions, tumor/Normal, bedgrapgs for visualization, read-level methylatyion, stattics of genomewide methylation)
3. minimal workflow for modified basecalling (optional)


future
-batched modified basecalling per sample (for optimized basecalling)

# To run
#sotwares are in singularity container, so can be executed using withut snakemake
singularity exec /data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif dorado --version

#Test run basecalling with singularity container in iris with GPU 
singularity exec --nv /data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif dorado --version

# To run with snakemake
1. install snakemake workflow (if you don't already have it) using the guidelines 
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

2. Run modified basecalling
$ snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/modified_basecalling.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab" --keep-going --forceall -np

NB: pls change path to reference genomes & type of species in `/config/config.yaml` 
b. specifiy the full path to whwere you cloned this repository to
ie. `parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/"`


ONT_tools.sif includes
dorado - modified basecalling (issued by ONT)
modkit - extracting methylation information from Bam & DMR (issued by ONT)
longphase - methylation-based phasing (https://github.com/twolinin/longphase)
pycoQC - interative QC plots (mapping & basecalling quality) ()
samtools - merge, sort and index bam file (Heng Li)
rust/cargo (dependency)

missing
GATK gatkmarkduplicates 

#get gatk docker
`docker pull broadinstitute/gatk`

```
FAQ
i want it get sampling rates of each read?
$pod5 view 

how do you split 1 big pod5 into different sampling rates
$pod5 subset

i want to run modbam for multiple modifcations 5mCG_5hmCG and 6mA
dorado basecaller sup,5mCG_5hmCG@latest,6mA@latest ${O44N_v14Test} --device "cuda:all" --reference $refhg38 --verbose
note if any is not found the pipline would stop

i want summary.tsv file of basecalled bam
dorado summary {input} --verbose > {output.seq_summary} 2> {log}

my bedmethyl is showing coverage instead of methylation
- name your file with `.bedmethyl` extension






### first activate conda
$ conda activate snakemake #we assume there's conda smk environment already installed on your device

### enter the directory of the repository
$cd dorado_ont_wf


## 2 step process
Step 1 (optional). convert .fast5 to .pod5 files
a. edit `samples_fast5_2_pod5.yaml` with paths to fast5
b. lunch program with `$ sh run_Snakefile_fast5_pod5.sh`

Notes on pod5 conversion
i. fast5 to pod5 is a looses conversion. Implies you can delete your original fast5 after conversion without loosing sleep.
ii. `pod5 convert fast5` allows multi-threading with `--threads #` option
iii. pod5 can be installed via `pip install pod5`
iv. You don't need unecessary memory for pod5 conversion. 3GB is enough. You do need lots of cpus though for threading

Step 2. modified base calling with dorado on .pod5 files and subsequent methylation extraction with modkit

## notes on run time
1. pod5 conversion
you dont need lots of memory for pod5 conversion. Max mem=4GB used. Rather increase nthreads=64, nGPUs=4

#merge run scripts
#run snakmake 
sh run_Snakefile_fast5_pod5.sh Snakemake_basecalling_modkit.smk config/cluster_pod5.json


#other snakemake piplines
`methylation_calling_modkit.smk` - methylation calling with modkit only
`split_reads.smk` - split pod files by sample_rate or channels
`Snakemake_basecalling_modkit.smk` - end to end basecalling and methylation calling with modkit



####################################################################################
# note to self
### test run with ``
$snakemake -s Snakefile.smk -np #target filename may change
$snakemake -s Snakefile.smk --cores 12 --forcerun -np #dry run with cores

#run actual pipeline on the cluster
$nohup snakemake -s Snakefile.smk --latency-wait 60 --restart-times 2 --keep-going --forceall --cluster "bsub -J {rule} -R "rusage[mem=32]" -W 1:00 -n 12 -o logs/cluster/{rule}.%J.out -e logs/cluster/{rule}.%J.err" -j 3 &

#alternatievely make a script that run snakmake
$ sh run_snakefile.sh 
$ cat run_snakefile.sh 
#!/bin/bash

# Run snakemake
snakemake --jobname 's.{jobid}.{rulename}' \
	--snakefile Snakefile_agg_stats_ONT.smk \
    --use-conda \
	--keep-going \
	--reason \
	--printshellcmds \
	--latency-wait 10 \
	--rerun-incomplete \
	--stats snakemake_$(date +"%Y%m%d_%H%M%S").stats \
	-j 500 \
	--cluster-config config/cluster.json \
	--cluster "bsub -q {cluster.queue} -n {cluster.threads} -W {cluster.time} -M{cluster.mem} -R\"span[hosts=1] select[mem>{cluster.mem}] rusage[mem={cluster.mem}]\" {cluster.extra} -o out.txt -e err.txt" 


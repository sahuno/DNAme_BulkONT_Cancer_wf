# Snakefile
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/merge_bedmethyl.smk \
# --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" \
# --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs unlimited --cores all -np

import os
import os, glob

# 0) path to your ModKit Singularity image
MODKIT_SIF = "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"

# 1) load treatments from sample.yaml 
configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/treatments_dmr_groups.yaml"
SAMPLES = config["samples"]



# 2) fixed inputs & params
GENOME_CHROM_SIZES = "/data1/greenbab/database/mm10/mm10.sorted.chrom.sizes"
GENOME_FASTA       = "/data1/greenbab/database/mm10/mm10.fa"
LINE1_PROMOTERS    = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/outputs/LINE1_promoters/first_100bp.bed"



# list of treatment groups
SAMPLE_GROUPS = list(SAMPLES.keys())
rule all:
    input:
        expand("results/mergedBedmethyl/{group}/merged_{group}.done.txt",
               group=SAMPLE_GROUPS),
        expand("results/mergedBedmethyl/{group}/merged_{group}.bed",
               group=SAMPLE_GROUPS),
        expand("results/mergedBedmethyl/{group}/merged_{group}.bed.gz",
                group=SAMPLE_GROUPS)
        # expand("results/mergedBedmethyl/{group}/merged_{group}.bw",
        #         group=SAMPLE_GROUPS)

rule run_merge_bedmethyl:
    """
    Run modkit to merge mutiple bed methyls.
    """
    output:
        outdir = directory("results/mergedBedmethyl/{group}"),
        merged_bed= "results/mergedBedmethyl/{group}/merged_{group}.bed",
        merged_bedgz= "results/mergedBedmethyl/{group}/merged_{group}.bed.gz",
        # merged_bedbw= "results/mergedBedmethyl/{group}/merged_{group}.bw",
        done   = touch("results/mergedBedmethyl/{group}/merged_{group}.done.txt")
    params:
        sampleGroup     = lambda wc: " ".join(
                        f"{path} "
                        for name, path in SAMPLES[wc.group].items()
                     ),
        genomeSizes    = GENOME_CHROM_SIZES
    log:
        bedmethyl_log="logs/{group}/mergedBedmethyl_{group}.log"
        # bw_log="logs/{group}/mergedBedmethyl_{group}_bw.log"
    threads: 8
    singularity:
        MODKIT_SIF
    shell:
        r"""
         echo "merging {params.sampleGroup} into {output.merged_bed}"
        
        modkit bedmethyl merge {params.sampleGroup} --out-bed {output.merged_bed} \
          -t {threads} \
          --force \
          --header \
            --genome-sizes {params.genomeSizes} \
          --log-filepath {log.bedmethyl_log}

        touch {output.done}


        sort -k1,1 -k2,2n {output.merged_bed} | bgzip -c > {output.merged_bedgz} 
        
        tabix -p bed {output.merged_bedgz}
        """
#          --header \
        # echo "Converting {output.merged_bedgz} to bedmethyl format"

        # modkit bedmethyl tobigwig --mod-codes m \
        # --suppress-progress \
        # --nthreads {threads} \
        # --sizes {params.genomeSizes} \
        # --log-filepath {log.bw_log} \
        #   {output.merged_bedgz} {output.merged_bedbw}
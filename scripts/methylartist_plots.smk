import os
import pandas as pd
import subprocess


# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/scripts/methylartist_plots.smk \
# --workflow-profile /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/slurmMinimal \
# --cores all --forcerun -np


# ───────────────────────────────────────────────
# USER CONFIGURATION
# ───────────────────────────────────────────────
# Path to your region→ID file
REGION_FILE   = "/data1/greenbab/projects/Sarcoma_DNAme/resources/regions/igv_report.proteins_t2t_chm13v2.txt"
# Genome sizes for bedtools slop/flank
GENOME_SIZES  = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.genome.sizes"
# methylartist references
REF_FA        = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.0.fa"
REF_GTF       = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.0_RefSeq_Liftoff_v5.2.sorted.gtf.gz"
# Define your tumor/normal BAM pairs here
# key = sample name, value = (normal_bam, tumor_bam)
PAIRED_BAMS = {
    "SHAH_H003458": (
        "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/"
        "modBasecalling/results/mark_duplicates/SHAH_H003458_N01_01_WG01_TCDO_SAR_032_FT_N/"
        "SHAH_H003458_N01_01_WG01_TCDO_SAR_032_FT_N_modBaseCalls_dedup_sorted.bam",
        "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/"
        "modBasecalling/results/mark_duplicates/SHAH_H003458_T01_01_WG01_TCDO_SAR_032_PDX_A/"
        "SHAH_H003458_T01_01_WG01_TCDO_SAR_032_PDX_A_modBaseCalls_dedup_sorted.bam"
    ),
    "SHAH_H003459": (
        "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/"
        "modBasecalling/results/mark_duplicates/SHAH_H003459_N01_01_WG01_TCDO_SAR_033_FT_N/"
        "SHAH_H003459_N01_01_WG01_TCDO_SAR_033_FT_N_modBaseCalls_dedup_sorted.bam",
        "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/"
        "modBasecalling/results/mark_duplicates/SHAH_H003459_T01_01_WG01_TCDO_SAR_033_PDX_A/"
        "SHAH_H003459_T01_01_WG01_TCDO_SAR_033_PDX_A_modBaseCalls_dedup_sorted.bam"
    ),
    "SHAH_H003842": (
        "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/"
        "modBasecalling/results/mark_duplicates/SHAH_H003842_N01_02_WG01_TCDO_SAR_061_FT_N/"
        "SHAH_H003842_N01_02_WG01_TCDO_SAR_061_FT_N_modBaseCalls_dedup_sorted.bam",
        "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/"
        "modBasecalling/results/mark_duplicates/SHAH_H003842_T01_01_WG02_TCDO_SAR_061_FT_A/"
        "SHAH_H003842_T01_01_WG02_TCDO_SAR_061_FT_A_modBaseCalls_dedup_sorted.bam"
    )
}
# methylartist settings
MOTIF         = "CG"
MODS          = "m"
PROPORTION    = "1,6,1,3,4"

# Flank & Slop sizes (bp)
FLANK_SIZE    = 1000
SLOP_SIZE     = 1000

# Paths to snakemake bedtools
SNK_BEDTOOLS = "/home/ahunos/miniforge3/envs/snakemake/bin/bedtools"

# Parse region file
if not os.path.exists(REGION_FILE):
    raise FileNotFoundError(f"Region file not found: {REGION_FILE}")
regions = []
with open(REGION_FILE) as fh:
    for line in fh:
        chrom_coords, tid = line.strip().split()
        regions.append((tid, chrom_coords))
REGION_DICT = dict(regions)

# Helper functions

def flank(region):
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    cmd = f"printf '{chrom}\t{start}\t{end}' | {SNK_BEDTOOLS} flank -i stdin -g {GENOME_SIZES} -l {FLANK_SIZE} -r 0"
    out = subprocess.check_output(cmd, shell=True).decode().split()
    return f"{out[1]}-{out[2]}"


def slop(region):
    chrom, coords = region.split(":")
    start, end = coords.split("-")
    cmd = f"printf '{chrom}\t{start}\t{end}' | {SNK_BEDTOOLS} slop -i stdin -g {GENOME_SIZES} -l {SLOP_SIZE} -r 0"
    out = subprocess.check_output(cmd, shell=True).decode().split()
    return f"{out[0]}:{out[1]}-{out[2]}"

# Snakemake rules
rule all:
    input:
        expand(
            "plots/{sample}/{tid}.png",
            sample=list(PAIRED_BAMS.keys()),
            tid=list(REGION_DICT.keys())
        )

rule methylartist_locus:
    output:
        png="plots/{sample}/{tid}.png"
    log:
        "logs/{sample}/{tid}.log"
    params:
        bams=lambda wc: ",".join(PAIRED_BAMS[wc.sample]),
        region=lambda wc: slop(REGION_DICT[wc.tid]),
        flank=lambda wc: flank(REGION_DICT[wc.tid])
    shell:
        r"""
        echo "Running methylartist locus {{wildcards.region}}...\n for {{params.bams}}"
        echo "highlighting region: {{params.flank}}"

        source /home/ahunos/miniforge3/etc/profile.d/conda.sh
        conda activate /home/ahunos/miniforge3/envs/methylartist

        methylartist locus \
          -b {params.bams} \
          -i {params.region} \
          --motif {MOTIF} \
          --mods {MODS} \
          --ref {REF_FA} \
          -g {REF_GTF} \
          -p {PROPORTION} \
          --outfile {output.png} \
        &> {log}
        """
        #   -l {params.flank} \

# methylartist scoredist \
# -b /data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/modBasecalling/results/mark_duplicates/SHAH_H003458_N01_01_WG01_TCDO_SAR_032_FT_N/SHAH_H003458_N01_01_WG01_TCDO_SAR_032_FT_N_modBaseCalls_dedup_sorted.bam,/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/modBasecalling/results/mark_duplicates/SHAH_H003458_T01_01_WG01_TCDO_SAR_032_PDX_A/SHAH_H003458_T01_01_WG01_TCDO_SAR_032_PDX_A_modBaseCalls_dedup_sorted.bam,/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/modBasecalling/results/mark_duplicates/SHAH_H003459_N01_01_WG01_TCDO_SAR_033_FT_N/SHAH_H003459_N01_01_WG01_TCDO_SAR_033_FT_N_modBaseCalls_dedup_sorted.bam,/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/modBasecalling/results/mark_duplicates/SHAH_H003459_T01_01_WG01_TCDO_SAR_033_PDX_A/SHAH_H003459_T01_01_WG01_TCDO_SAR_033_PDX_A_modBaseCalls_dedup_sorted.bam,/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/modBasecalling/results/mark_duplicates/SHAH_H003842_N01_02_WG01_TCDO_SAR_061_FT_N/SHAH_H003842_N01_02_WG01_TCDO_SAR_061_FT_N_modBaseCalls_dedup_sorted.bam,/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/modBasecalling/results/mark_duplicates/SHAH_H003842_T01_01_WG02_TCDO_SAR_061_FT_A/SHAH_H003842_T01_01_WG02_TCDO_SAR_061_FT_A_modBaseCalls_dedup_sorted.bam \
# -m "m" \
# --ref /data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.0.fa \
# --motif CG \
# --outfile global_methylation.png

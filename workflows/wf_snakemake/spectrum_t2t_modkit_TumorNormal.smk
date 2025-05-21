# snakemake -s /data1/greenbab/projects/Sarcoma_DNAme/scripts/workflows/t2t_modkit_TumorNormal.smk \
# --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/modkit_slurmMinimal \
# --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --quiet --forceall -np

# snakemake -s /data1/greenbab/projects/Sarcoma_DNAme/scripts/workflows/t2t_modkit_TumorNormal.smk \
# --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/modkit_slurmMinimal \
# --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" \
# --allowed-rules highQual_modkit_pileupCpGs_sortTabixBed,highQual_modkit_pileupCpGs_BigWigs, modkit_dmr_pairCpGs,modkit_dmr_pair_cross,Unphased_modkit_5mC_pileupBedgraph,Unphased_modkit_5mC_pileup -np


import pandas as pd
import yaml
import os

# ─── 0) CONFIG ────────────────────────────────────────────────────────────────
YAML_CONFIG  = "config.yaml"  # New configuration file for YAML settings

SAMPLES_FILE = "/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/t2t_samples_TumorNormal_clairs.csv"
#uncomment to test on chr19
# SAMPLES_FILE = "/data1/greenbab/projects/Sarcoma_DNAme/data/processed/t2t_analysis/subset_modBam/results/chr19_subset/manifest_chr19_bams.tsv"

#pair_id,tumor,normal
#pairID="pair_044TN"
##T044="/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/hac_5mCG_5hmCG/results/mark_duplicates/044T_v14/044T_v14_modBaseCalls_sorted_dup.bam"
##N044="/data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/sup_5mCG_5hmCG/results/mark_duplicates/044N_v14/044N_v14_modBaseCalls_sorted_dup.bam"



SEED = 1234
OUTDIR       = "results"

# IMAGES
MODKIT_IMG   = "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"  #
ONT_TOOLS_IMG= "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"

# ─── 0) REF ──────────────────────────────────────────────────────────────────
REF          = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.0.fa"
CONTIG    = "chr19"

# RUNNING RESOURCES
THREADS      = 12

# ─── SV CALLING CONFIG ────────────────────────────────────────────────────
SV_OUTDIR          = f"{OUTDIR}/sv"
NORMAL_FILTER_EXPR = "PASS"
TUMOR_FILTER_EXPR  = "PASS"

# Read sample metadata
samples_df = pd.read_csv(SAMPLES_FILE, index_col=0)
SAMPLES    = list(samples_df.index)
VCF_TYPES  = ["tumor", "normal"]

minCov = 5
haplotypes = ["1", "2", "ungrouped"]

MODCODES = ["h", "m","a"]
GENOMESIZES = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.genome.sizes"
GENOMESIZES = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.genome.sizes"



# Pick the correct raw VCF for filtering
def get_raw_vcf(wc):
    base = f"{OUTDIR}/run_clair/{wc.sample}"
    return f"{base}/{wc.sample}.snv.vcf.gz" if wc.type == "tumor" else f"{base}/clair3_normal_germline_output.vcf.gz"


# ─── FINAL OUTPUTS ─────────────────────────────────────────────────────────────
rule all:
    input:        
        # Unphased CpG pileup BED outputs
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased.bed",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_sorted.bed.gz",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_filtered_sorted.bed",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_filtered_sorted.bed.gz",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_filtered.bw",
            type=VCF_TYPES, sample=SAMPLES
        ),

        # Unphased bedGraph outputs
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined.bedgraph",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered.bedgraph",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz.tbi",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz.tbi",
            type=VCF_TYPES, sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/done.{{type}}.{{sample}}.txt",
            type=VCF_TYPES, sample=SAMPLES
        ),
        # unphased DMR results

        # Unphased DMR results (now including raw beds)
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_raw_dmr.bed",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_raw_segmentation.bed",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr.bed.gz",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr.bed.gz.tbi",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_segmentation.bed.gz",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_segmentation.bed.gz.tbi",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_diff.bed",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_diff.bed.gz",
            sample=SAMPLES
        ),
        expand(
            f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_diff.bed.gz.tbi",
            sample=SAMPLES
        )

########################################################
# Unphased ModKit 5mC pileup for Tumor & Normal Cohorts
########################################################
rule Unphased_modkit_5mC_pileup:
    input:
        bam = lambda wc: samples_df.loc[wc.sample, wc.type]
    output:
        raw_bed              = f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased.bed",
        sorted_bed_gz        = f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_sorted.bed.gz",
        filtered_bed         = f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_filtered_sorted.bed",
        filtered_bed_gz      = f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_filtered_sorted.bed.gz",
        bigwig               = f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_filtered.bw"
    params:
        ref         = REF,
        mincov      = minCov,
        genome_sizes= GENOMESIZES
    threads: THREADS
    singularity: MODKIT_IMG
    log:
        f"logs/Unphased_modkit_5mC_pileup/{{type}}/{{sample}}/{{type}}_{{sample}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/Unphased_modkit_5mC_pileup/{wildcards.type}/{wildcards.sample}

        # 1) Raw unphased CpG pileup BED
        modkit pileup \
          --ref {params.ref} \
          --cpg --combine-strands \
          --ignore h --sampling-frac 0.8 --filter-threshold 0.7 \
          --threads {threads} \
          --prefix {wildcards.type}_{wildcards.sample} \
          {input.bam} {output.raw_bed}

        # 2) Sort & compress raw BED
        sort -k1,1 -k2,2n {output.raw_bed} \
          | bgzip -c > {output.sorted_bed_gz}

        # 3) Filter by coverage & sort
        sort -k1,1 -k2,2n {output.raw_bed} \
          | awk -v min_cov={params.mincov} '$5 >= min_cov' \
          > {output.filtered_bed}

        # 4) Compress & index filtered BED
        bgzip -c {output.filtered_bed} > {output.filtered_bed_gz}
        tabix -p bed {output.filtered_bed_gz}
        tabix -p bed {output.sorted_bed_gz}

        # 5) Generate 5mC BigWig from filtered BED
        modkit bedmethyl tobigwig \
          --mod-codes m --suppress-progress \
          --nthreads {threads} \
          --sizes {params.genome_sizes} \
          --log-filepath {log} \
          {output.filtered_bed} {output.bigwig}
        """


rule Unphased_modkit_5mC_pileupBedgraph:
    input:
        bam = lambda wc: samples_df.loc[wc.sample, wc.type]
    output:
        outdir                  = directory(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/"),
        # raw bedGraph
        raw_bg                  = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined.bedgraph",
        # filtered bedGraph
        filt_bg                 = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered.bedgraph",
        # sorted & compressed raw
        raw_bg_sorted_gz        = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz",
        raw_bg_sorted_tbi       = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz.tbi",
        # sorted & compressed filtered
        filt_bg_sorted_gz       = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz",
        filt_bg_sorted_tbi      = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz.tbi",
        # BigWig from filtered
        bigwig                  = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}_unphased_m_CG0_combined_filtered.bw",
        done                    = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/done.{{type}}.{{sample}}.txt"
    params:
        ref         = REF,
        mincov      = minCov,
        seed        = SEED,
        genome_sizes= GENOMESIZES
    threads: THREADS
    singularity: MODKIT_IMG
    log:
        f"logs/Unphased_modkit_5mC_pileupBedgraph/{{type}}/{{sample}}/{{type}}_{{sample}}.log"
    shell:
        r"""
        mkdir -p {OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{wildcards.type}/{wildcards.sample}

        # 1) Generate raw bedGraph
        modkit pileup \
          --ref {params.ref} \
          --cpg --combine-strands \
          --ignore h --bedgraph \
          --filter-threshold 0.7 --sampling-frac 0.8 --seed {params.seed} \
          --threads {threads} \
          --prefix {wildcards.type}_{wildcards.sample}_unphased \
          {input.bam} {output.outdir}

        # 2) Sort & index raw bedGraph
        sort -k1,1 -k2,2n {output.raw_bg} | bgzip -c > {output.raw_bg_sorted_gz}
        tabix -p bed {output.raw_bg_sorted_gz}

        # 3) Filter raw bedGraph by coverage
        awk -v min_cov={params.mincov} '$5 >= min_cov' {output.raw_bg} > {output.filt_bg}

        # 4) Sort & index filtered bedGraph
        sort -k1,1 -k2,2n {output.filt_bg} | bgzip -c > {output.filt_bg_sorted_gz}
        tabix -p bed {output.filt_bg_sorted_gz}

        # 5) Convert filtered bedGraph to BigWig
        /data1/greenbab/users/ahunos/apps/ucsc_tools/bedGraphToBigWig \
          {output.filt_bg_sorted_gz} \
          {params.genome_sizes} \
          {output.bigwig}

        # 6) Mark done
        touch {output.done}
        """

rule modkit_dmr_unphased:
    """
    Call DMRs comparing unphased normal vs tumor 5mC BED outputs,
    then emit raw and post‐processed beds, bgzip/tabix index, and filtered DMRs.
    """
    input:
        normal_bed = f"{OUTDIR}/Unphased_modkit_5mC_pileup/normal/{{sample}}/normal_{{sample}}_unphased_filtered_sorted.bed.gz",
        tumor_bed  = f"{OUTDIR}/Unphased_modkit_5mC_pileup/tumor/{{sample}}/tumor_{{sample}}_unphased_filtered_sorted.bed.gz"
    output:
        # keep raw outputs
        raw_dmr      = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_raw_dmr.bed",
        raw_seg      = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_raw_segmentation.bed",
        # compressed/indexed DMR
        dmr_gz       = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr.bed.gz",
        dmr_tbi      = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr.bed.gz.tbi",
        # compressed/indexed segmentation
        seg_gz       = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_segmentation.bed.gz",
        seg_tbi      = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_segmentation.bed.gz.tbi",
        # filtered DMRs (uncompressed + compressed/indexed)
        dmr_diff     = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_diff.bed",
        dmr_diff_gz  = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_diff.bed.gz",
        dmr_diff_tbi = f"{OUTDIR}/modkit_dmr_unphased/{{sample}}/{{sample}}_unphased_TumorNormal_dmr_diff.bed.gz.tbi"
    params:
        ref           = REF,
        base          = "C",
        mincov        = minCov,
        min_dmr_sites = 3
    threads: THREADS
    singularity: MODKIT_IMG
    shell:
        r"""
        mkdir -p {OUTDIR}/modkit_dmr_unphased/{wildcards.sample}

        # 1) Raw DMR + segmentation
        modkit dmr pair \
          -a {input.normal_bed} \
          -b {input.tumor_bed} \
          -o {output.raw_dmr} \
          --segment {output.raw_seg} \
          --ref {params.ref} \
          --base {params.base} \
          --min-valid-coverage {params.mincov} \
          --threads {threads} \
          --header

        # 2) Compress & index raw DMR
        sort -k1,1 -k2,2n {output.raw_dmr} \
          | bgzip -c > {output.dmr_gz}
        tabix -p bed {output.dmr_gz}

        # 3) Compress & index raw segmentation
        sort -k1,1 -k2,2n {output.raw_seg} \
          | bgzip -c > {output.seg_gz}
        tabix -p bed {output.seg_gz}

        # 4) Filter for “different” DMRs
        awk -F '\t' 'NR>1 && $4=="different" && $6>={params.min_dmr_sites} && ($14>=0.5||$14<=-0.5) && ($15*$16>0)' \
          {output.raw_seg} > {output.dmr_diff}

        # 5) Compress & index filtered DMRs
        sort -k1,1 -k2,2n {output.dmr_diff} \
          | bgzip -c > {output.dmr_diff_gz}
        tabix -p bed {output.dmr_diff_gz}
        """

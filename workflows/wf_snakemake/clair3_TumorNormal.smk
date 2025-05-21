# Snakefile

import pandas as pd

# ─── 0) CONFIG ────────────────────────────────────────────────────────────────
SAMPLES_FILE = "/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/samples_TumorNormal_clairs.csv"       # CSV must have columns: tumor, normal, clairs_model
OUTDIR       = "results"
BIND_DIR     = "/data1/greenbab"
IMAGE        = "/data1/greenbab/users/ahunos/apps/containers/clairs_latest.sif"
CTG_NAMES    = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
REF          = "/data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta"
CONDAPREFIX  = "/opt/conda/envs/clairs"
THREADS      = 12

# ─── 1) READ SAMPLE METADATA ─────────────────────────────────────────────────
samples_df = pd.read_csv(SAMPLES_FILE, index_col = 0)
SAMPLES    = list(samples_df.index)

# ─── 2) RULE: all ──────────────────────────────────────────────────────────────
rule all:
    input:
        expand("results/run_clair/{sample}/{sample}.snp.vcf.gz", sample=SAMPLES),
        expand("results/run_clair/{sample}/done.{sample}.txt", sample=SAMPLES),
        expand("results/filter_pass/{sample}/{sample}.pass.vcf.gz", sample=SAMPLES),
        expand("results/filter_pass/{sample}/{sample}.pass.vcf.gz.tbi", sample=SAMPLES),
        expand('results/longphase_modcall/{sample}/modcall_{sample}.vcf', sample=SAMPLES),
        expand('results/longphase_modcall/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand('results/longphase_coPhase/{sample}/{sample}.vcf', sample=SAMPLES),
        expand('results/longphase_coPhase/{sample}/{sample}_mod.vcf', sample=SAMPLES),
        expand('results/longphase_coPhase/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand('results/longphase_haplotag/{sample}/{sample}_haplotagged.bam', sample=SAMPLES),
        expand('results/longphase_haplotag/{sample}/{sample}_haplotagged.bam.bai', sample=SAMPLES),
        expand('results/longphase_haplotag/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand("results/mutation_NeoantigenInput/{sample}/{sample}_mutation_NeoantigenInput.txt", sample=SAMPLES),
        expand('results/tumor_haplotag/{sample}/{sample}_Tumor_haplotagged.bam', sample=SAMPLES),
        expand('results/tumor_haplotag/{sample}/{sample}_Tumor_haplotagged.bam.bai', sample=SAMPLES),
        expand('results/tumor_haplotag/{sample}/done.{sample}.txt', sample=SAMPLES)
        # expand(f"somaticVarClairs/{OUTDIR}/{sample}", sample=SAMPLES),

# ─── 3) RULE: run_clair ───────────────────────────────────────────────────────
rule run_clair:
    input:
        tumor  = lambda wc: samples_df.loc[wc.sample, "tumor"],
        normal = lambda wc: samples_df.loc[wc.sample, "normal"],
    output:
        outdir=directory("results/run_clair/{sample}"),
        vcf="results/run_clair/{sample}/{sample}.snp.vcf.gz",
        doneflag="results/run_clair/{sample}/done.{sample}.txt"
    threads: THREADS
    resources:
        mem_mb = 768000
    singularity: IMAGE
    params:
        bind         = BIND_DIR,
        image        = IMAGE,
        ctg_names    = CTG_NAMES,
        ref          = REF,
        conda_prefix = CONDAPREFIX,
        # pick up the right Clair3 model per‐sample:
        platform     = lambda wc: samples_df.loc[wc.sample, "clairs_model"]
    log:
        "logs/run_clair/{sample}/{sample}.log"
    shell:
        """
        mkdir -p {output.outdir}
          /opt/bin/run_clairs \
            --tumor_bam_fn   {input.tumor} \
            --normal_bam_fn  {input.normal} \
            --output_prefix {wildcards.sample}.snv \
            --ref_fn         {params.ref} \
            --threads        {threads} \
            --qual           5 \
            --platform       '{params.platform}' \
            --enable_indel_calling \
            --enable_clair3_germline_output \
            --output_dir     {output.outdir} \
            --ctg_name       {params.ctg_names} \
            --conda_prefix   {params.conda_prefix} \
            && touch {output.doneflag} 2> {log}
        """
                        # --print_ref_calls \

            # --ctg_name       {params.ctg_names} \

            # --ctg_name       {params.ctg_names} \
            # --use_longphase_for_intermediate_haplotagging True \

rule filter_pass:
    input:
        vcf = "results/run_clair/{sample}/{sample}.snp.vcf.gz"
    output:
        pass_vcf = "results/filter_pass/{sample}/{sample}.pass.vcf.gz",
        pass_tbi = "results/filter_pass/{sample}/{sample}.pass.vcf.gz.tbi",
        outdir=directory("results/filter_pass/{sample}")
    # singularity: "docker://quay.io/biocontainers/bcftools:1.16--h5a_2"
    shell:
        """
        # extract only PASS records, bgzip & index
        mkdir -p {output.outdir}
        /admin/software/bcftools/bcftools-1.20/bin/bcftools view -f PASS {input.vcf} -Oz -o {output.pass_vcf}
        /data1/greenbab/users/ahunos/apps/htslib-1.20/tabix -p vcf {output.pass_vcf}
        """

rule longphase_modcall:
    input:
        bams=lambda wc: samples_df.loc[wc.sample, "tumor"]
    threads: THREADS
    params:
        reference_genome=REF,
        threads=12,
        out_modcall_suffix='results/longphase_modcall/{sample}/modcall_{sample}'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        done_modcall='results/longphase_modcall/{sample}/done.{sample}.txt',
        out_modcall='results/longphase_modcall/{sample}/modcall_{sample}.vcf'
    log:
      "logs/longphase_modcall/{sample}/{sample}.log"
    shell:
        """
longphase modcall -b {input.bams} -t {threads} -o {params.out_modcall_suffix} -r {params.reference_genome} && touch {output.done_modcall} 2> {log}
        """
# cp results/longphase_modcall/{wildcards.sample}/modcall_{wildcards.sample}.vcf.vcf results/longphase_modcall/{wildcards.sample}/modcall_{wildcards.sample}.vcf


rule longphase_coPhase:
    input:
        bamfile=lambda wc: samples_df.loc[wc.sample, "tumor"],
        modcallfile='results/longphase_modcall/{sample}/modcall_{sample}.vcf',
        snpFile="results/filter_pass/{sample}/{sample}.pass.vcf.gz"
    params:
        reference_genome=REF,
        threads=12,
        out_coPhase_prefix='results/longphase_coPhase/{sample}/{sample}'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        co_phased_mod_vcf='results/longphase_coPhase/{sample}/{sample}_mod.vcf',
        co_phased_vcf='results/longphase_coPhase/{sample}/{sample}.vcf',
        done_longphase_coPhase='results/longphase_coPhase/{sample}/done.{sample}.txt'
    log:
      "logs/longphase_coPhase/{sample}/{sample}.log"
    shell:
        """
longphase phase \
-s {input.snpFile} \
--mod-file {input.modcallfile} \
-b {input.bamfile} \
-r {params.reference_genome} \
-t {params.threads} \
-o {params.out_coPhase_prefix} \
--ont && touch {output.done_longphase_coPhase}
       """

rule tumor_haplotag:
    input:
        bamfile=lambda wc: samples_df.loc[wc.sample, "tumor"],
        svCalls_vcf='results/longphase_coPhase/{sample}/{sample}_mod.vcf',
        tumorPhased_vcf='results/filter_pass/{sample}/{sample}.pass.vcf.gz'
    params:
        reference_genome=REF,
        threads=12,
        out_haplotagged_prefix='results/tumor_haplotag/{sample}/{sample}_Tumor_haplotagged'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        happlotagged_bam='results/tumor_haplotag/{sample}/{sample}_Tumor_haplotagged.bam',
        happlotagged_bam_bai='results/tumor_haplotag/{sample}/{sample}_Tumor_haplotagged.bam.bai',
        done_longphase_haplotag='results/tumor_haplotag/{sample}/done.{sample}.txt'
    log:
      "logs/tumor_haplotag/{sample}/{sample}.log"
    shell:
        """
longphase haplotag \
-s {input.tumorPhased_vcf} \
-r {params.reference_genome} \
-b {input.bamfile} \
-t 8 \
-o {params.out_haplotagged_prefix} && touch {output.done_longphase_haplotag} 2> {log}

samtools index {output.happlotagged_bam}
        """

rule longphase_haplotag:
    input:
        bamfile=lambda wc: samples_df.loc[wc.sample, "tumor"],
        co_phased_mod_vcf='results/longphase_coPhase/{sample}/{sample}_mod.vcf',
        co_phased_vcf='results/longphase_coPhase/{sample}/{sample}.vcf'
    params:
        reference_genome=REF,
        threads=12,
        out_haplotagged_prefix='results/longphase_haplotag/{sample}/{sample}_haplotagged'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        happlotagged_bam='results/longphase_haplotag/{sample}/{sample}_haplotagged.bam',
        happlotagged_bam_bai='results/longphase_haplotag/{sample}/{sample}_haplotagged.bam.bai',
        done_longphase_haplotag='results/longphase_haplotag/{sample}/done.{sample}.txt'
    log:
      "logs/longphase_haplotag/{sample}/{sample}.log"
    shell:
        """
longphase haplotag \
-s {input.co_phased_vcf} \
--mod-file {input.co_phased_mod_vcf} \
-r {params.reference_genome} \
-b {input.bamfile} \
-t 8 \
-o {params.out_haplotagged_prefix} && touch {output.done_longphase_haplotag} 2> {log}

samtools index {output.happlotagged_bam}
        """

rule mutation_NeoantigenInput:
    """
    From the PASS VCF, strip any leading 'chr' from CHROM and emit
    lines like chrom_pos_ref_alt into a plain text file.
    """
    input:
        vcf = "results/longphase_coPhase/{sample}/{sample}.vcf"
    output:
        txt = "results/mutation_NeoantigenInput/{sample}/{sample}_mutation_NeoantigenInput.txt",
        outdir=directory("results/mutation_NeoantigenInput/{sample}")
    params:
        # outdir = lambda wc: os.path.dirname(wc.output.txt)
    # singularity: "docker://quay.io/biocontainers/bcftools:1.16--h5a_2"
    shell:
        r"""
        mkdir -p {output.outdir}
        /admin/software/bcftools/bcftools-1.20/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} \
          | sed 's/^chr//' \
          | awk '{{print $1"_"$2"_"$3"_"$4}}' \
          > {output.txt}
        """

# chrom_pos_ref_alt


                    # --use_gpu \
            # --ctg_name       {params.ctg_names} \

            # --dry_run
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/clair3_TumorNormal.smk \
# --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/slurmMinimal \
# --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --quiet --forceall -np

# --use-singularity --configfile config.yaml --cores 24
# --forceall 
#--rerun-incomplete -
# Snakefile

import pandas as pd

# ─── 0) CONFIG ────────────────────────────────────────────────────────────────
SAMPLES_FILE = "/data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/samples_TumorNormal_clairs.csv"       # CSV must have columns: tumor, normal, clairs_model
OUTDIR       = "results"
BIND_DIR     = "/data1/greenbab"
IMAGE        = "/data1/greenbab/users/ahunos/apps/containers/clairs_latest.sif"
CTG_NAMES    = "chr20"
REF          = "/data1/greenbab/database/hg38/v0/Homo_sapiens_assembly38.fasta"
CONDAPREFIX  = "/opt/conda/envs/clairs"
THREADS      = 12

# ─── 1) READ SAMPLE METADATA ─────────────────────────────────────────────────
samples_df = pd.read_csv(SAMPLES_FILE, index_col = 0)
SAMPLES    = list(samples_df.index)

# ─── 2) RULE: all ──────────────────────────────────────────────────────────────
rule all:
    input:
        expand(f"somaticVarClairs/{OUTDIR}/{{sample}}", sample=SAMPLES)

# ─── 3) RULE: run_clair ───────────────────────────────────────────────────────
rule run_clair:
    input:
        tumor  = lambda wc: samples_df.loc[wc.sample, "tumor"],
        normal = lambda wc: samples_df.loc[wc.sample, "normal"],
    output:
        directory(f"somaticVarClairs/{OUTDIR}/{{sample}}"),
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
    shell:
        """
        mkdir -p {output}
          /opt/bin/run_clairs \
            --tumor_bam_fn   {input.tumor} \
            --normal_bam_fn  {input.normal} \
            --output_prefix {wildcards.sample} \
            --ref_fn         {params.ref} \
            --threads        {threads} \
            --print_ref_calls \
            --qual           5 \
            --platform       '{params.platform}' \
            --enable_indel_calling \
            --enable_clair3_germline_output \
            --use_longphase_for_intermediate_haplotagging True \
            --output_dir     {output} \
            --conda_prefix   {params.conda_prefix} \
        """
                    # --use_gpu \
            # --ctg_name       {params.ctg_names} \

            # --dry_run
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/clair3_TumorNormal.smk \
# --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/slurmMinimal \
# --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --forceall --quiet -np

# --use-singularity --configfile config.yaml --cores 24
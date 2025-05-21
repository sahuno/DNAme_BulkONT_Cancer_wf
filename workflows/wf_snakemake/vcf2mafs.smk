
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/somaticVariations/sampleRate5000/results/snpeff
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/somaticVariations/sampleRate4000/results/snpeff

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/vcf2mafs.smk --cores all --forcerun -np

configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/somaticVariations/samples_snpEff.yaml"
REF = "/data1/greenbab/database/mm10/mm10.fa"

# singularity run /data1/greenbab/users/ahunos/apps/containers/vcf2maf_vep88_1.3.0.sif vcf2maf --help

rule all:
  input:
    expand("vcf2mafs/maf/{sample}.maf", sample=config["samples"]),
    "vcf2mafs/cohort.maf",
    # "vcf2mafs/plots/oncoplot.pdf",
    # "vcf2mafs/driver_genes.tsv"


rule vcf2maf:
  input:
    vcf=lambda wildcards: config["samples"][wildcards.sample],
    ref=REF
  output:
    maf="vcf2mafs/maf/{sample}.maf"
  singularity:
    "/data1/greenbab/users/ahunos/apps/containers/vcf2maf_vep88_1.3.0.sif"    
  shell:
    """
    vcf2maf \
      --input-vcf {input.vcf} \
      --output-maf {output.maf} \
      --ref-fasta {input.ref} \
      --tumor-id {wildcards.sample} \
      --vep-path /usr/local/bin/vep \
      --vep-data /path/to/vep_data
    """


rule merge_maf:
  input:
    maf=expand("vcf2mafs/maf/{sample}.maf", sample=config["samples"])
  output:
    "vcf2mafs/cohort.maf"
  run:
    shell("head -n1 {input.maf[0]} > {output}")
    shell("tail -q -n +2 {input.maf} >> {output}")

# rule make_oncoplot:
#   input:
#     maf="vcf2mafs/cohort.maf"
#   output:
#     pdf="vcf2mafs/plots/oncoplot.pdf"
#   conda:
#     "envs/maftools.yaml"
#   script:
#     "scripts/plot_oncoplot.R"

# rule identify_drivers:
#   input:
#     maf="vcf2mafs/cohort.maf"
#   output:
#     "vcf2mafs/driver_genes.tsv"
#   conda:
#     "envs/dndscv.yaml"
#   script:
#     "scripts/find_drivers.R"

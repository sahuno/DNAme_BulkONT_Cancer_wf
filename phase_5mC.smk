import os
import glob
from datetime import datetime

#purpose: phase snps and 5mc
#remove `--ctg_name="chr19"` options in the `rule call_snps_indels`
# TODO: check why HP tags don't appear in the haplottaged bam file
# Assign the current date to a variable
current_datetime = datetime.now()
formatted_datetime = current_datetime.strftime('%Y%m%d_%H_%M_%S')

parent_dir = "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/"
#set species
set_species = "mouse"

configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_bams_4000sampleRate.yaml" #mouse samples
print(config["samples"]) #sanity check


rule all:
    input:
        expand('results/call_snps_indels/{samples}/snv.vcf.gz', samples=config["samples"]),
        expand('results/call_snps_indels/{samples}/indel.vcf.gz', samples=config["samples"]),
        expand('results/call_snps_indels/{samples}/done.{samples}.txt', samples=config["samples"]),
        expand('results/longphase_modcall/{samples}/modcall_{samples}.vcf', samples=config["samples"]),
        expand('results/longphase_modcall/{samples}/done.{samples}.txt', samples=config["samples"]),
        expand('results/longphase_phase/{samples}/{samples}.vcf', samples=config["samples"]),
        expand('results/longphase_phase/{samples}/{samples}_mod.vcf', samples=config["samples"]),
        expand('results/longphase_phase/{samples}/done.{samples}.txt', samples=config["samples"]),
        expand('results/longphase_haplotag/{samples}/{samples}_haplotagged.bam', samples=config["samples"]),
        expand('results/longphase_haplotag/{samples}/{samples}_haplotagged.bam.bai', samples=config["samples"])

rule call_snps_indels:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        software_dir=config["software_dir"],
        out_dir='results/call_snps_indels/{samples}/',
        threads=12,
        PLATFORM='ont_r10_dorado_sup_4khz'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/clairs-to_latest.sif"
    output:
        done='results/call_snps_indels/{samples}/done.{samples}.txt',
        snps='results/call_snps_indels/{samples}/snv.vcf.gz',
        indels='results/call_snps_indels/{samples}/indel.vcf.gz'
        # sorted_bed='results/sortedBams/{samples}.sorted.bam.bai'
    log:
      "logs/call_snps_indels/{samples}/{samples}.log"
    shell:
        """
/opt/bin/run_clairs_to \
  --tumor_bam_fn {input} \
  --ctg_name="chr19" \
  --ref_fn {params.reference_genome} \
  --threads {params.threads} \
  --platform {params.PLATFORM} \
  --output_dir {params.out_dir} \
  --conda_prefix /opt/micromamba/envs/clairs-to && touch {output.done} 2> {log}
        """
#  /opt/bin/run_clairs_to --help
#  --ctg_name="chr19" \

# set output in params directive
rule longphase_modcall:
    input:
        bams=lambda wildcards: config["samples"][wildcards.samples]
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        software_dir=config["software_dir"],
        threads=12,
        out_modcall_prefix='results/longphase_modcall/{samples}/modcall_{samples}'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        done_modcall='results/longphase_modcall/{samples}/done.{samples}.txt',
        out_modcall='results/longphase_modcall/{samples}/modcall_{samples}.vcf'
    log:
      "logs/longphase_modcall/{samples}/{samples}.log"
    shell:
        """
longphase modcall -b {input.bams} -t 8 -o {output.out_modcall} -r {params.reference_genome} && touch {output.done_modcall} 2> {log}
mv results/longphase_modcall/{wildcards.samples}/modcall_{wildcards.samples}.vcf.vcf results/longphase_modcall/{wildcards.samples}/modcall_{wildcards.samples}.vcf
        """

rule longphase_phase:
    input:
        bamfile=lambda wildcards: config["samples"][wildcards.samples],
        modcallfile='results/longphase_modcall/{samples}/modcall_{samples}.vcf',
        snpFile='results/call_snps_indels/{samples}/snv.vcf.gz'
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        threads=12,
        out_lphase_prefix='results/longphase_phase/{samples}/{samples}'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        co_phased_mod_vcf='results/longphase_phase/{samples}/{samples}_mod.vcf',
        co_phased_vcf='results/longphase_phase/{samples}/{samples}.vcf',
        done_longphase_phase='results/longphase_phase/{samples}/done.{samples}.txt'
    log:
      "logs/longphase_phase/{samples}/{samples}.log"
    shell:
        """
longphase phase \
-s {input.snpFile} \
--mod-file {input.modcallfile} \
-b {input.bamfile} \
-r {params.reference_genome} \
-t {params.threads} \
-o {params.out_lphase_prefix} \
--ont && touch {output.done_longphase_phase}
       """

rule longphase_haplotag:
    input:
        bamfile=lambda wildcards: config["samples"][wildcards.samples],
        snpFile='results/call_snps_indels/{samples}/snv.vcf.gz',
        co_phased_mod_vcf='results/longphase_phase/{samples}/{samples}_mod.vcf'
        # modcall='results/longphase_modcall/{samples}/modcall_{samples}.vcf'
        # co_phased_vcf='results/longphase_phase/{samples}/{samples}.vcf'
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        software_dir=config["software_dir"],
        threads=12,
        out_haplotagged_prefix='results/longphase_haplotag/{samples}/{samples}_haplotagged'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        happlotagged_bam='results/longphase_haplotag/{samples}/{samples}_haplotagged.bam',
        happlotagged_bam_bai='results/longphase_haplotag/{samples}/{samples}_haplotagged.bam.bai'
    log:
      "logs/longphase_haplotag/{samples}/{samples}.log"
    shell:
        """
longphase haplotag \
-s {input.snpFile} \
--mod-file {input.co_phased_mod_vcf} \
-r {params.reference_genome} \
-b {input.bamfile} \
-t 8 \
-o {params.out_haplotagged_prefix} 2> {log}
        """
# --sv-file phased_sv.vcf \ #add sv calls if available

# mv results/longphase_phase/{wildcards.samples}/{wildcards.samples}_mod.vcf.vcf results/longphase_phase/{wildcards.samples}/{wildcards.samples}_mod.vcf
# mv results/longphase_phase/{wildcards.samples}/{wildcards.samples}.vcf.vcf results/longphase_phase/{wildcards.samples}/{wildcards.samples}.vcf
 



#run interactively
# snakemake -s /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/phase_5mC.smk  --cores 12 --forcerun --use-singularity --singularity-args "\"--bind /data1/greenbab\"" -np

#run on slurm
# snakemake -s /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/phase_5mC.smk --workflow-profile /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/phased_modifications/config/slurm --jobs 10 --cores all --use-singularity --singularity-args "'-B "/data1/greenbab"'" --keep-going --forceall -np


# get options
# /data1/greenbab/users/ahunos/apps/containers/clairs-to_latest.sif /opt/bin/run_clairs_to --help

#  singularity shell clair3_latest.sif
# singularity shell clairs_latest.sif
# /opt/bin/run_clair3.sh
# /opt/bin/run_clairs.sh

# find . -name "NEB-WGEMS-7T_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam"
# find . -name "*NEB-WGEMS-7T_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam*"
# mamba create -n methphaser -c bioconda methphaser
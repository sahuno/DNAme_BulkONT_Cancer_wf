#Author: Samuel Ahuno
#purpose: workflow to co-phase 5mc 

import os
import glob
from datetime import datetime

#purpose: phase snps and 5mc
#remove `--ctg_name="chr19"` options in the `rule call_snps_indels`
current_datetime = datetime.now()
formatted_datetime = current_datetime.strftime('%Y%m%d_%H_%M_%S')

parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/"
#set species
set_species = "mouse"

configfile: parent_dir + "config/config.yaml"
# configfile: parent_dir + "config/samples_bams_5000sampleRate.yaml" #mouse samples
configfile: parent_dir + "config/samples_bams_4000sampleRate.yaml" #mouse samples
print(config["samples"]) #sanity check

SAMPLES = list(config["samples"].keys())

motif = "CG"
strands = ["negative", "positive"]
haplotypes = ["1", "2", "ungrouped"]
mods = ["h", "m"]

# pileupRegion: chr19
# pileupBedgraph: True


rule all:
    input:
        expand('results/call_snps_indels/{sample}/{sample}.vcf.gz', sample=SAMPLES),
        expand('results/call_snps_indels/{sample}/indel.vcf.gz', sample=SAMPLES),
        expand('results/call_snps_indels/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand('results/longphase_modcall/{sample}/modcall_{sample}.vcf', sample=SAMPLES),
        expand('results/longphase_modcall/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand('results/longphase_coPhase/{sample}/{sample}.vcf', sample=SAMPLES),
        expand('results/longphase_coPhase/{sample}/{sample}_mod.vcf', sample=SAMPLES),
        expand('results/longphase_coPhase/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand('results/longphase_haplotag/{sample}/{sample}_haplotagged.bam', sample=SAMPLES),
        expand('results/longphase_haplotag/{sample}/{sample}_haplotagged.bam.bai', sample=SAMPLES),
        expand('results/longphase_haplotag/{sample}/done.{sample}.txt', sample=SAMPLES),
        expand("results/filter_pass/{sample}/{sample}.pass.vcf.gz", sample=SAMPLES),
        expand("results/filter_pass/{sample}/{sample}.pass.vcf.gz.tbi", sample=SAMPLES),
        expand("results/mutation_NeoantigenInput/{sample}/{sample}_mutation_NeoantigenInput.txt", sample=SAMPLES),
        expand("results/mPileup_HapModStrandSpecific/{sample}/{sample}.done", sample=SAMPLES),
        # final annotated VCFs and stats under results/snpeff/{sample}/
        expand("results/snpeff/{sample}/{sample}.ann.vcf", sample=SAMPLES),
        expand("results/snpeff/{sample}/{sample}.snpEff_summary.html", sample=SAMPLES)
        # expand("results/mPileupTabixSortBed/{sample}/tabixSorted/{sample}_{haplotypes}.sorted.bed", sample=SAMPLES, haplotypes=haplotypes),
        # expand("results/mPileupTabixSortBed/{sample}/tabixSorted/{sample}_{haplotypes}.sorted.bed.gz", sample=SAMPLES, haplotypes=haplotypes)

# ...exi


rule call_snps_indels:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        out_dir='results/call_snps_indels/{sample}/',
        threads=12,
        PLATFORM='ont_r10_dorado_sup_4khz'
    resources:
        mem_mb = 768000
    singularity: "/data1/greenbab/users/ahunos/apps/containers/clairs-to_latest.sif"
    output:
        done='results/call_snps_indels/{sample}/done.{sample}.txt',
        snps='results/call_snps_indels/{sample}/{sample}.vcf.gz',
        indels='results/call_snps_indels/{sample}/indel.vcf.gz'
        # sorted_bed='results/sortedBams/{sample}.sorted.bam.bai'
    log:
      "logs/call_snps_indels/{sample}/{sample}.log"
    shell:
        """
/opt/bin/run_clairs_to \
  --tumor_bam_fn {input} \
  --ref_fn {params.reference_genome} \
  --threads {params.threads} \
  --use_longphase_for_intermediate_phasing true \
  --use_longphase_for_intermediate_haplotagging true \
  --platform {params.PLATFORM} \
  --output_dir {params.out_dir} \
  --snv_output_prefix {wildcards.sample} \
  --conda_prefix /opt/micromamba/envs/clairs-to \
  && touch {output.done} 2> {log}
        """
    #   --ctg_name="chr19" \
#   --ctg_name="chr19" \

#  /opt/bin/run_clairs_to --help
#  --ctg_name="chr19" # for testing command
# --use_longphase_for_intermediate_haplotagging #until 5hmC & other modifcations are ready

rule filter_pass:
    input:
        vcf = "results/call_snps_indels/{sample}/{sample}.vcf.gz"
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
# set output in params directive
rule longphase_modcall:
    input:
        bams=lambda wildcards: config["samples"][wildcards.sample]
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        threads=12,
        out_modcall_prefix='results/longphase_modcall/{sample}/modcall_{sample}'
    singularity: "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"
    output:
        done_modcall='results/longphase_modcall/{sample}/done.{sample}.txt',
        out_modcall='results/longphase_modcall/{sample}/modcall_{sample}.vcf'
    log:
      "logs/longphase_modcall/{sample}/{sample}.log"
    shell:
        """
longphase modcall -b {input.bams} -t 8 -o {output.out_modcall} -r {params.reference_genome} && touch {output.done_modcall} 2> {log}
mv results/longphase_modcall/{wildcards.sample}/modcall_{wildcards.sample}.vcf.vcf results/longphase_modcall/{wildcards.sample}/modcall_{wildcards.sample}.vcf
        """


rule longphase_coPhase:
    input:
        bamfile=lambda wildcards: config["samples"][wildcards.sample],
        modcallfile='results/longphase_modcall/{sample}/modcall_{sample}.vcf',
        snpFile='results/call_snps_indels/{sample}/{sample}.vcf.gz'
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
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

rule longphase_haplotag:
    input:
        bamfile=lambda wildcards: config["samples"][wildcards.sample],
        co_phased_mod_vcf='results/longphase_coPhase/{sample}/{sample}_mod.vcf',
        co_phased_vcf='results/longphase_coPhase/{sample}/{sample}.vcf'
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
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
        outdir = directory("results/mutation_NeoantigenInput/{sample}")
    # params:
    #     outdir = lambda wc: os.path.dirname(wc.output.txt)
    # singularity: "docker://quay.io/biocontainers/bcftools:1.16--h5a_2"
    shell:
        r"""
        mkdir -p {output.outdir}
        /admin/software/bcftools/bcftools-1.20/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} \
          | sed 's/^chr//' \
          | awk '{{print $1"_"$2"_"$3"_"$4}}' \
          > {output.txt}
        """

rule mPileup_HapModStrandSpecific:
    input:
        hapTagBam='results/longphase_haplotag/{sample}/{sample}_haplotagged.bam'
    output:
        done="results/mPileup_HapModStrandSpecific/{sample}/{sample}.done"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=16,
        modkit_prob_percentiles=0.1,
        partition_tag="HP",
        region_arg = lambda wc: (
            f"--region {config['pileupRegion']}"
            if config.get("pileupRegion") else ""
        ),
        outtype_arg = lambda wc: (
            f"--bedgraph"
            if config.get("pileupBedgraph") else ""
        ),
        outdir = lambda wc: f"results/mPileup_HapModStrandSpecific/{wc.sample}",
        log=lambda wc: f"results/mPileup_HapModStrandSpecific/{wc.sample}/{wc.sample}.log"
    shell:
        """
        mkdir -p {params.outdir}

        ~/.cargo/bin/modkit pileup --threads {params.modkit_threads} \
        {params.outtype_arg} \
        --prefix {wildcards.sample} \
        --motif CG 0 \
        --ref {params.reference_genome} \
        --partition-tag {params.partition_tag} \
        {params.region_arg} \
        --log-filepath {params.log} \
        {input} \
        {params.outdir} && touch {output.done}
        """

rule snpeff:
    """
    Annotate a VCF with SnpEff inside Singularity.
    Outputs go to results/snpeff/{sample}/
    """
    input:
        vcf = "results/longphase_coPhase/{sample}/{sample}.vcf"
    output:
        annotated_vcf="results/snpeff/{sample}/{sample}.ann.vcf",
        stats="results/snpeff/{sample}/{sample}.snpEff_summary.html"
    params:
        img=lambda wc: os.path.join(config["imagesDir"], "snpeff_4.3t.sif"),
        bind_db=config["snpeffDB"],
        genome=config["genome"],
        ud=config["ud"]
    log:
        "logs/snpeff/{sample}.log"
    shell:
        r"""
        # ensure output dir exists
        mkdir -p $(dirname {output.annotated_vcf})

        singularity exec \
           --bind /data1/greenbab \
           --bind {params.bind_db}:/snpEff/data \
           {params.img} \
        java -jar /snpEff/snpEff.jar \
           -v \
           -canon \
           -lof \
           -ud {params.ud} \
           -dataDir /snpEff/data \
           {params.genome} \
           {input.vcf} \
           -o vcf \
           -stats {output.stats} \
        > {output.annotated_vcf} 2> {log}
        """

# rule mPileupTabixSortBed:
#     input:
#         unsorted_bed = "results/mPileup_HapModStrandSpecific/{sample}/{sample}_{haplotypes}.bed"
#         # unsorted_bed = lambda wc: os.path.join("results/mPileup_HapModStrandSpecific",wc.sample, f"{wc.sample}_{wc.haplotypes}.bed")

#     output:
#         sortedBed= "results/mPileupTabixSortBed/{sample}/{sample}_{haplotypes}.sorted.bed",
#         tabixSorted="results/mPileupTabixSortBed/{sample}/sample}_{haplotypes}.sorted.bed.gz"
#     params:
#         outdir = lambda wc: f"results/mPileupTabixSortBed/{wc.sample}"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         /data1/greenbab/users/ahunos/apps/htslib-1.20/bgzip -k {input.unsorted_bed} --output {output.sortedBed}
#         /data1/greenbab/users/ahunos/apps/htslib-1.20/tabix -p bed {output.sortedBed}
#         """





# --sv-file phased_sv.vcf \ #add sv calls if available

# mv results/longphase_coPhase/{wildcards.sample}/{wildcards.sample}_mod.vcf.vcf results/longphase_coPhase/{wildcards.sample}/{wildcards.sample}_mod.vcf
# mv results/longphase_coPhase/{wildcards.sample}/{wildcards.sample}.vcf.vcf results/longphase_coPhase/{wildcards.sample}/{wildcards.sample}.vcf
 


# see clair options 
# singularity run /data1/greenbab/users/ahunos/apps/containers/clairs-to_latest.sif /opt/bin/run_clairs_to --help
# singularity shell /data1/greenbab/users/ahunos/apps/containers/clairs-to_latest.sif /opt/bin/run_clairs_to --help


#D01Bam=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/sandbox/results/call_snps_indels/D-0-1_4000/tmp/phasing_output/phased_bam_output/tumor_chr19.bam
# modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
# --prefix {wildcards.sample} --cpg --combine-mods --ref {params.reference_genome}
# modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
# --prefix {wildcards.sample} --cpg --ref {params.reference_genome}



# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/somaticVar_ClairTumorOnly_revisedSample4000.smk \
# --workflow-profile /data1/greenbab/projects/Sarcoma_DNAme/scripts/configs/slurmMinimal \
# --jobs 10 --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" --keep-going --forceall --quiet -np

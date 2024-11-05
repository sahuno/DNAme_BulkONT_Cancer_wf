parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_haplotaggedBams.yaml"

set_species = "mouse"

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/testHaplotype_modkit.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs unlimited --cores all --use-singularity -np


# Load samples from config
samplesFromConfig = list(config["samples"].keys())
print(f"{samplesFromConfig}")

# Wildcard variables
motifs = ["h", "m", "C"]
Haplotypes = ["1", "2", "ungrouped"]
strands = ["positive", "negative"]
minCoverages = 5

rule all:
    input:
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.bedgraph",
               sample=samplesFromConfig, motif=motifs, strand=strands, Haplotypes=Haplotypes, allow_missing=True),
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_{Haplotypes}_{motif}_CG0_combined.bedgraph",
               sample=samplesFromConfig, motif=motifs, Haplotypes=Haplotypes, allow_missing=True),
        expand("results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_{Haplotypes}_C_CG0_combined.bedgraph",
               sample=samplesFromConfig, Haplotypes=Haplotypes),
        expand("results/filterSortBedgraphs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCoverages}.sorted.bedgraph",
               sample=samplesFromConfig, motif=motifs, strand=strands, Haplotypes=Haplotypes, minCoverages=[minCoverages])

rule ModkitPileup_ModiStrandSpecific:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_h_CG0_negative.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_h_CG0_positive.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_m_CG0_negative.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_m_CG0_positive.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_h_CG0_negative.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_h_CG0_positive.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_m_CG0_negative.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_m_CG0_positive.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_h_CG0_negative.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_h_CG0_positive.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_m_CG0_negative.bedgraph",
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_m_CG0_positive.bedgraph"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=16,
        outdir="results/ModkitPileup_ModiStrandSpecific/{sample}/"
    shell:
        """
        modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
        --prefix {wildcards.sample} --cpg --ref {params.reference_genome} --partition-tag HP
        """

rule ModkitPileup_ModiSpecific:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "results/ModkitPileup_ModiSpecific/{sample}/{sample}_1_h_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiSpecific/{sample}/{sample}_1_m_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiSpecific/{sample}/{sample}_2_h_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiSpecific/{sample}/{sample}_2_m_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiSpecific/{sample}/{sample}_ungrouped_h_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiSpecific/{sample}/{sample}_ungrouped_m_CG0_combined.bedgraph"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=16,
        outdir="results/ModkitPileup_ModiSpecific/{sample}/"
    shell:
        """
        modkit pileup --combine-strands --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
        --prefix {wildcards.sample} --cpg --ref {params.reference_genome} --partition-tag HP
        """

rule ModkitPileup_ModiStrandCollapsed:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_1_C_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_2_C_CG0_combined.bedgraph",
        "results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_ungrouped_C_CG0_combined.bedgraph"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=16,
        outdir="results/ModkitPileup_ModiStrandCollapsed/{sample}/"
    shell:
        """
        modkit pileup --combine-mods --combine-strands --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
        --prefix {wildcards.sample} --cpg --ref {params.reference_genome} --partition-tag HP
        """

rule filterSortBedgraphs:
    input:
        "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.bedgraph",
        #"results/ModkitPileup_ModiSpecific/{sample}/{sample}_{Haplotypes}_{motif}_CG0_combined.bedgraph",
        #"results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_{Haplotypes}_C_CG0_combined.bedgraph"
    output:
        "results/filterSortBedgraphs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCoverages}.sorted.bedgraph"
    params:
        minCov=minCoverages
    shell:
        """
        awk -v min_cov="{params.minCov}" '$5 > min_cov{{print $1, $2, $3, $4}}' "{input}" | sort -k1,1 -k2,2n > {output}
        """

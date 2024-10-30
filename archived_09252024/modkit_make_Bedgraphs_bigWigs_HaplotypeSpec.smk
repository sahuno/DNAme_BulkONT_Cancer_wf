parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_haplotaggedBams.yaml"

set_species = "mouse"

# Make sample list of bedgraphs
motifs = ["h", "m", "C"]
Haplotypes = ["1", "2", "ungrouped"]
strands = ["positive", "negative"]
minCoverages = [5, 10]

samplesFromConfig = list(config["samples"].keys())
print(f"{samplesFromConfig}")

rule all:
    input:
        expand("results/BedGraphs2BigWigs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bw", 
               sample=samplesFromConfig, motif=motifs, strand=strands, minCov=minCoverages, Haplotypes=Haplotypes)

rule modkit:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "results/modkit/{sample}/{sample}_{Haplotypes}_h_CG0_negative.bedgraph",
        "results/modkit/{sample}/{sample}_{Haplotypes}_h_CG0_positive.bedgraph",
        "results/modkit/{sample}/{sample}_{Haplotypes}_m_CG0_positive.bedgraph",
        "results/modkit/{sample}/{sample}_{Haplotypes}_m_CG0_negative.bedgraph",
        "results/modkit/{sample}/{sample}_{Haplotypes}_C_CG0_negative.bedgraph",
        "results/modkit/{sample}/{sample}_{Haplotypes}_C_CG0_positive.bedgraph"
    params:
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        modkit_threads=16,
        modkit_prob_percentiles=0.1,
        outdir="results/modkit/{sample}/"
    shell:
        """
        modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
        --prefix {wildcards.sample} --cpg --combine-mods --ref {params.reference_genome}
        modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
        --prefix {wildcards.sample} --cpg --ref {params.reference_genome}
        """
#  --region chr17
rule filterSortBedgraphs:
    input:
        "results/modkit/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.bedgraph"
    output:
        "results/filterSortBedgraphs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bedgraph"
    wildcard_constraints:
        minCov="|".join(map(str, minCoverages))
    shell:
        """
        awk -v min_cov="{wildcards.minCov}" '$5 > min_cov{{print $1, $2, $3, $4}}' "{input}" | sort -k1,1 -k2,2n > {output}
        """

rule BedGraphs2BigWigs:
    input:
        "results/filterSortBedgraphs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bedgraph"
    output:
        "results/BedGraphs2BigWigs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bw"
    params:
        chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
    shell:
        """
        bedGraphToBigWig {input} {params.chromSizes} {output}
        """

ruleorder: modkit > filterSortBedgraphs > BedGraphs2BigWigs

# sk

# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/modkit_make_Bedgraphs_bigWigs_HaplotypeSpec.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs unlimited --cores all --use-singularity -np

# check if reads are phased?
# samtools view $input.bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^HP:|^PS:/) {print $0; break}}' | less


#template of ooutputs
# -rw-rw-r-- 1 ahunos ahunos 4.8M Oct 29 16:21 D02_1_C_CG0_combined.bedgraph
# -rw-rw-r-- 1 ahunos ahunos 5.4M Oct 29 16:21 D02_2_C_CG0_combined.bedgraph
# -rw-rw-r-- 1 ahunos ahunos  15M Oct 29 16:21 D02_ungrouped_C_CG0_combined.bedgraph

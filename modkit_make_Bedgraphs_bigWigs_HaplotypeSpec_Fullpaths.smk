#samuel ahuno
# how to run
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/modkit_make_Bedgraphs_bigWigs_HaplotypeSpec_Fullpaths.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs unlimited --cores all --use-singularity -np

# to rerun
# rm -rf results .snakemake
# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/modkit_make_Bedgraphs_bigWigs_HaplotypeSpec_Fullpaths.smk

import glob


parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: parent_dir + "config/samples_haplotaggedBams.yaml"

set_species = "mouse"

# Make sample list of bedgraphs
samplesFromConfig = list(config["samples"].keys())
print(f"{samplesFromConfig}")

#wildcards
motifs = ["h", "m", "C"]
motifs_for_splitModNStrand = ["h", "m"]
Haplotypes = ["1", "2", "ungrouped"]
strands = ["positive", "negative"]
# minCoverages = [5, 10]
minCoverages=5

# def find_bedgraph_files():
#     bedgraphs_inputs=glob.glob("results/filteredSortedBedgraphsBigwigs/**/*.bedgraph", recursive=True)
#     bigwig_outputs=[file.replace("filteredSortedBedgraphsBigwigs", "filteredSortedBigWigs").replace(".bedgraph", ".bw") for file in bedgraphs_inputs]
#     return bedgraphs_inputs, bigwig_outputs

# file_pairs = find_bedgraph_files()

##Analysis process 
#1. with strand 5mC/5hmC collapsed, is there difference haploype 1  vrs haplotype 2? 
# expand("results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_1_C_CG0_combined.bedgraph", sample=samplesFromConfig),
# expand("results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_2_C_CG0_combined.bedgraph", sample=samplesFromConfig),

rule all:
    input:
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_h_CG0_negative.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_h_CG0_positive.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_m_CG0_negative.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_1_m_CG0_positive.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_h_CG0_negative.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_h_CG0_positive.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_m_CG0_negative.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_2_m_CG0_positive.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_h_CG0_negative.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_h_CG0_positive.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_m_CG0_negative.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_ungrouped_m_CG0_positive.bedgraph", sample=samplesFromConfig),
        # ModkitPileup_ModiSpecific outputs
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_1_h_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_1_m_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_2_h_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_2_m_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_ungrouped_h_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiSpecific/{sample}/{sample}_ungrouped_m_CG0_combined.bedgraph", sample=samplesFromConfig),
        # ModkitPileup_ModiStrandCollapsed outputs
        expand("results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_1_C_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_2_C_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_ungrouped_C_CG0_combined.bedgraph", sample=samplesFromConfig),
        expand("results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_{strand}.minCov{minCov}.sorted.bedgraph", sample=samplesFromConfig, motifs_for_splitModNStrand=motifs_for_splitModNStrand, strand=strands, minCov=[minCoverages], Haplotypes=Haplotypes, allow_missing=True),
        expand("results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandSpecific/{sample}/bigWig/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_{strand}.minCov{minCov}.sorted.bw", sample=samplesFromConfig, motifs_for_splitModNStrand=motifs_for_splitModNStrand, strand=strands, minCov=[minCoverages], Haplotypes=Haplotypes, allow_missing=True),
        expand("results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiSpecific/{sample}/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_combined.minCov{minCov}.sorted.bedgraph", sample=samplesFromConfig, motifs_for_splitModNStrand=motifs_for_splitModNStrand, minCov=[minCoverages], Haplotypes=Haplotypes, allow_missing=True),
        expand("results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiSpecific/{sample}/bigWig/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_combined.minCov{minCov}.sorted.bw", sample=samplesFromConfig, motifs_for_splitModNStrand=motifs_for_splitModNStrand, minCov=[minCoverages], Haplotypes=Haplotypes, allow_missing=True),
        expand("results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_{Haplotypes}_C_CG0_combined.minCov{minCov}.sorted.bedgraph", sample=samplesFromConfig, minCov=[minCoverages], Haplotypes=Haplotypes, allow_missing=True),
        expand("results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandCollapsed/{sample}/bigWig/{sample}_{Haplotypes}_C_CG0_combined.minCov{minCov}.sorted.bw", sample=samplesFromConfig, minCov=[minCoverages], Haplotypes=Haplotypes, allow_missing=True),
        


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
        modkit_prob_percentiles=0.1,
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
        modkit_prob_percentiles=0.1,
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
        modkit_prob_percentiles=0.1,
        outdir="results/ModkitPileup_ModiStrandCollapsed/{sample}/"
    shell:
        """
        modkit pileup --combine-mods --combine-strands --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
        --prefix {wildcards.sample} --cpg --ref {params.reference_genome} --partition-tag HP
        """

rule filteredSortedBedgraphsBigwigs_ModiStrandSpecific:
    input:
        f"results/ModkitPileup_ModiStrandSpecific/{{sample}}/{{sample}}_{{Haplotypes}}_{{motifs_for_splitModNStrand}}_CG0_{{strand}}.bedgraph"
    output:
        bgs="results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_{strand}.minCov{minCov}.sorted.bedgraph",
        bws="results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandSpecific/{sample}/bigWig/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_{strand}.minCov{minCov}.sorted.bw"
    params:
        minCov=minCoverages,
        chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
    shell:
        """
        awk -v min_cov="{params.minCov}" '$5 > min_cov{{print $1, $2, $3, $4}}' "{input}" | sort -k1,1 -k2,2n > {output.bgs}
        bedGraphToBigWig {output.bgs} {params.chromSizes} {output.bws}
        """

rule filteredSortedBedgraphsBigwigs_ModiSpecific:
    input:
        f"results/ModkitPileup_ModiSpecific/{{sample}}/{{sample}}_{{Haplotypes}}_{{motifs_for_splitModNStrand}}_CG0_combined.bedgraph",

    output:
        bgs="results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiSpecific/{sample}/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_combined.minCov{minCov}.sorted.bedgraph",
        bws="results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiSpecific/{sample}/bigWig/{sample}_{Haplotypes}_{motifs_for_splitModNStrand}_CG0_combined.minCov{minCov}.sorted.bw"
    params:
        minCov=minCoverages,
        chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
        
    shell:
        """
        awk -v min_cov="{params.minCov}" '$5 > min_cov{{print $1, $2, $3, $4}}' "{input}" | sort -k1,1 -k2,2n > {output.bgs}
        bedGraphToBigWig {output.bgs} {params.chromSizes} {output.bws}
        """

rule filteredSortedBedgraphsBigwigs_ModiStrandCollapsed:
    input:
        f"results/ModkitPileup_ModiStrandCollapsed/{{sample}}/{{sample}}_{{Haplotypes}}_C_CG0_combined.bedgraph"
    output:
        bgs="results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_{Haplotypes}_C_CG0_combined.minCov{minCov}.sorted.bedgraph",
        bws="results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiStrandCollapsed/{sample}/bigWig/{sample}_{Haplotypes}_C_CG0_combined.minCov{minCov}.sorted.bw"
    params:
        minCov=minCoverages,
        chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
    shell:
        """
        awk -v min_cov="{params.minCov}" '$5 > min_cov{{print $1, $2, $3, $4}}' "{input}" | sort -k1,1 -k2,2n > {output.bgs}
        bedGraphToBigWig {output.bgs} {params.chromSizes} {output.bws}
        """


# bedGraphToBigWig results/ModkitPileup_ModiSpecific/D-0-1/D-0-1_1_h_CG0_combined.bedgraph /data1/greenbab/database/mm10/mm10.chrom.sizes results/filteredSortedBedgraphsBigwigs/ModkitPileup_ModiSpecific/D-0-1/D-0-1_1_h_CG0_combined.minCov5.sorted.bw

# rule filteredBedGraphs2BigWigs:
#     input:
#         lambda wildcards: wildcards.bedgraph_file
#     output:
#         lambda wildcards: wildcards.bigwig_file
#     params:
#         chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
#     shell:
#         """
#         bedGraphToBigWig {input} {params.chromSizes} {output}
#         """


# bedgraph_inputs, bigwig_outputs = find_bedgraph_files()

# # find_bedgraph_files()[1]
# rule filteredBedGraphs2BigWigs:
#     input:
#         bedgraph_inputs
#     output:
#         bigwig_outputs
#     params:
#         chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
#     shell:
#         """
#         bedGraphToBigWig {input} {params.chromSizes} {output}
#         """



    #     lambda wildcards: bedgraph_inputs
    # output:
    #     lambda wildcards: bigwig_outputs

# results/ModkitPileup_ModiStrandCollapsed/*.bedgraph

# rule filteredBedGraphs2BigWigs:
#     input:
#         "results/filteredBedGraphs2BigWigs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bedgraph"
#     output:
#         "results/filteredBedGraphs2BigWigs/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bw"
#     params:
#         chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
#     shell:
#         """
#         bedGraphToBigWig {input} {params.chromSizes} {output}
#         """



# ruleorder: modkit > filteredSortedBedgraphsBigwigs > BedGraphs2BigWigs



        # if wildcards.motif in ["h", "m"] else
        # f"results/ModkitPileup_ModiSpecific/{wildcards.sample}/{wildcards.sample}_{wildcards.Haplotypes}_C_CG0_combined.bedgraph",
        # "results/ModkitPileup_ModiStrandSpecific/{sample}/{sample}_{Haplotypes}_{motif}_CG0_{strand}.bedgraph",
        # f"results/ModkitPileup_ModiSpecific/{wildcards.sample}/{wildcards.sample}_{wildcards.Haplotypes}_{wildcards.motif}_CG0_combined.bedgraph",
        #  "results/ModkitPileup_ModiStrandCollapsed/{sample}/{sample}_{Haplotypes}_C_CG0_combined.bedgraph"




# # sk

# # snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/modkit_make_Bedgraphs_bigWigs_HaplotypeSpec.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs unlimited --cores all --use-singularity -np
# # check if reads are phased?
# # samtools view $input.bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^HP:|^PS:/) {print $0; break}}' | less


# #template of ooutputs
# #  --combine-mods --combine-strands
# # modkit pileup --threads 16 --bedgraph $bamInprog D02Bam_modkit --prefix D02 --cpg --combine-mods --combine-strands --ref $mm10ref --partition-tag HP
# # D02_1_C_CG0_combined.bedgraph
# # D02_2_C_CG0_combined.bedgraph
# # D02_ungrouped_C_CG0_combined.bedgraph

# # --combine-strands, no --combine-mods 
# # modkit pileup --threads 16 --bedgraph $bamInprog D02_splitMod_modkit --prefix D02 --cpg --combine-strands --ref $mm10ref --partition-tag HP
# # D02_1_h_CG0_combined.bedgraph
# # D02_1_m_CG0_combined.bedgraph
# # D02_2_h_CG0_combined.bedgraph
# # D02_2_m_CG0_combined.bedgraph
# # D02_ungrouped_h_CG0_combined.bedgraph
# # D02_ungrouped_m_CG0_combined.bedgraph

# # no --combine-strands, no --combine-mods 
# # modkit pileup --threads 16 --bedgraph $bamInprog D02_splitMod_Strand_modkit --prefix D02 --cpg --ref $mm10ref --partition-tag HP
# # D02_1_h_CG0_negative.bedgraph
# # D02_1_h_CG0_positive.bedgraph
# # D02_1_m_CG0_negative.bedgraph
# # D02_1_m_CG0_positive.bedgraph
# # D02_2_h_CG0_negative.bedgraph
# # D02_2_h_CG0_positive.bedgraph
# # D02_2_m_CG0_negative.bedgraph
# # D02_2_m_CG0_positive.bedgraph
# # D02_ungrouped_h_CG0_negative.bedgraph
# # D02_ungrouped_h_CG0_positive.bedgraph
# # D02_ungrouped_m_CG0_negative.bedgraph
# # D02_ungrouped_m_CG0_positive.bedgraph




# # pl sprefix with `"results/ModkitPileup_ModiSpecific/{sample}/` like previous ones
# # # D02_1_h_CG0_combined.bedgraph
# # # D02_1_m_CG0_combined.bedgraph
# # # D02_2_h_CG0_combined.bedgraph
# # # D02_2_m_CG0_combined.bedgraph
# # # D02_ungrouped_h_CG0_combined.bedgraph
# # # D02_ungrouped_m_CG0_combined.bedgraph
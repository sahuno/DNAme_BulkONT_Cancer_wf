parent_dir = "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/"
configfile: parent_dir + "config/config.yaml"
configfile: "/data1/greenbab/users/ahunos/apps/dorado_ont_wf/config/samples_009TN_044TN.yaml"

set_species = "human"

# Make sample list of bedgraphs
motifs = ["h", "m", "C"]
strands = ["positive", "negative"]
minCoverages = [5, 10]

samplesFromConfig = list(config["samples"].keys())
print(f"{samplesFromConfig}")

rule all:
    input:
        expand("results/BedGraphs2BigWigs/{sample}/{sample}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bw", 
               sample=samplesFromConfig, motif=motifs, strand=strands, minCov=minCoverages)

rule modkit:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "results/modkit/{sample}/{sample}_h_CG0_negative.bedgraph",
        "results/modkit/{sample}/{sample}_h_CG0_positive.bedgraph",
        "results/modkit/{sample}/{sample}_m_CG0_positive.bedgraph",
        "results/modkit/{sample}/{sample}_m_CG0_negative.bedgraph",
        "results/modkit/{sample}/{sample}_C_CG0_negative.bedgraph",
        "results/modkit/{sample}/{sample}_C_CG0_positive.bedgraph"
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
        "results/modkit/{sample}/{sample}_{motif}_CG0_{strand}.bedgraph"
    output:
        "results/filterSortBedgraphs/{sample}/{sample}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bedgraph"
    wildcard_constraints:
        minCov="|".join(map(str, minCoverages))
    shell:
        """
        awk -v min_cov="{wildcards.minCov}" '$5 > min_cov{{print $1, $2, $3, $4}}' "{input}" | sort -k1,1 -k2,2n > {output}
        """

rule BedGraphs2BigWigs:
    input:
        "results/filterSortBedgraphs/{sample}/{sample}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bedgraph"
    output:
        "results/BedGraphs2BigWigs/{sample}/{sample}_{motif}_CG0_{strand}.minCov{minCov}.sorted.bw"
    params:
        chromSizes=lambda wildcards: config["mm10_chrSizes"] if set_species == "mouse" else config["hg38_chrSizes"]
    shell:
        """
        bedGraphToBigWig {input} {params.chromSizes} {output}
        """

ruleorder: modkit > filterSortBedgraphs > BedGraphs2BigWigs

# sk

#D01Bam=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/sandbox/results/call_snps_indels/D-0-1_4000/tmp/phasing_output/phased_bam_output/tumor_chr19.bam

#check if reads are phased?
# samtools view -h $D01Bam | less
# samtools view $D01Bam | awk '{for(i=12;i<=NF;i++) if($i ~ /^HP:|^PS:/) {print $0; break}}' | less

# mm10ref=/data1/greenbab/database/mm10/mm10.fa


# bamInprog=/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/snps_longphase_modcalls/s4000/sandbox/results/results/longphase_haplotag/D-0-2_4000/D-0-2_4000_haplotagged.bam

# samtools view -h $bamInprog | less
# samtools view $bamInprog | awk '{for(i=12;i<=NF;i++) if($i ~ /^HP:|^PS:/) {print $0; break}}'

### haplotype specific methylation. combine strands and mods
# modkit pileup --threads 16 --bedgraph $bamInprog D02Bam_modkit --prefix D02 --cpg --combine-mods --combine-strands --ref $mm10ref --partition-tag HP
# modkit pileup --threads 16 --bedgraph $bamInprog D02_splitMod_modkit --prefix D02 --cpg --combine-strands --ref $mm10ref --partition-tag HP


#template of ooutputs
# -rw-rw-r-- 1 ahunos ahunos 4.8M Oct 29 16:21 D02_1_C_CG0_combined.bedgraph
# -rw-rw-r-- 1 ahunos ahunos 5.4M Oct 29 16:21 D02_2_C_CG0_combined.bedgraph
# -rw-rw-r-- 1 ahunos ahunos  15M Oct 29 16:21 D02_ungrouped_C_CG0_combined.bedgraph

# modkit pileup --threads {params.modkit_threads} --bedgraph {input} {params.outdir} \
# --prefix {wildcards.sample} --cpg --ref {params.reference_genome}
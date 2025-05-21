# snakemake --snakefile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/workflows/wf_snakemake/modkit_methyl_extraction.smk --use-singularity \
# --workflow-profile /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/slurmMinimal \
# --cores all --use-singularity --singularity-args "--nv -B /data1/greenbab,/scratch/greenbab/ahunos,/data1/shahs3" -R bamCoverage_toBigWig -np



import pandas as pd
import yaml
import os


configfile: "/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/configs/haplotagged_bam_samples_Rate4000_n_Rate5000_merged.yaml" #mouse samples
SAMPLES = list(config["samples"].keys())


print("running modkit pipeline for samples in config file")
print(config["samples"]) #sanity check


# ─── 0) REF ──────────────────────────────────────────────────────────────────
REF          = "/data1/greenbab/database/mm10/mm10.fa"
GENOMESIZES = "/data1/greenbab/database/mm10/mm10.chrom.sizes"
# REF          = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.0.fa"
# GENOMESIZES = "/data1/greenbab/database/T2T_CHM13v2_plusY/chm13v2.genome.sizes"
CONTIG    = "chr19"
# ─── 1) OUTPUT ──────────────────────────────────────────────────────────────────
OUTDIR       = "results"
# IMAGES
MODKIT_IMG   = "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"  # or local: "/data1/greenbab/users/ahunos/apps/containers/modkit_latest.sif"
ONT_TOOLS_IMG= "/data1/greenbab/users/ahunos/apps/containers/ONT_tools.sif"


# RUNNING RESOURCES
THREADS      = 12
THREADS_SV   = 8


minCov = 5
haplotypes = ["1", "2", "ungrouped"]
MODCODES = ["h", "m"]
SEED = 1234



rule all:
    input:
        # expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed.gz",  sample=SAMPLES, haplotypes=haplotypes),
        expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_BigWigs/{{sample}}/{{sample}}_{{haplotypes}}_{{mod}}_sorted.bw",  sample=SAMPLES, haplotypes=haplotypes, mod=MODCODES),
        # modkit pileup CpGs bed (raw)
        expand(f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_1.bed", sample=SAMPLES),
        expand(f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_2.bed", sample=SAMPLES),
        expand(f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_ungrouped.bed", sample=SAMPLES),
        # modkit 5mC pileup CpGs bed (raw)
        expand(f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_1.bed", sample=SAMPLES),
        expand(f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_2.bed", sample=SAMPLES),
        expand(f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_ungrouped.bed", sample=SAMPLES),

        # high‐quality sorted & compressed CpGs beds
        expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotype}}_sorted.bed",sample=SAMPLES, haplotype=haplotypes),
        expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotype}}_sorted.bed.gz", sample=SAMPLES, haplotype=haplotypes),

        # high‐quality sorted & compressed 5mC beds
        expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_sortTabixBed/{{sample}}/{{sample}}_{{haplotype}}_sorted.bed",sample=SAMPLES, haplotype=haplotypes),
        expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_sortTabixBed/{{sample}}/{{sample}}_{{haplotype}}_sorted.bed.gz", sample=SAMPLES, haplotype=haplotypes),
        # expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed.gz",  sample=SAMPLES, haplotypes=haplotypes),
        # expand(f"{OUTDIR}/modkit_pileupCpGsBedgraph/{{sample}}/done.{{sample}}.txt",  sample=SAMPLES),
        # expand(f"{OUTDIR}/modkit_pileupCpGsBedgraph/{{sample}}/done.{{sample}}.txt", sample=SAMPLES),
        expand(f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_BigWigs/{{sample}}/{{sample}}_{{haplotype}}_{{mod}}_sorted.bw", sample=SAMPLES, haplotype=haplotypes, mod=MODCODES),
        # modkit unphased CpGs bed
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased.bed", sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_sorted.bed.gz", sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_filtered_sorted.bed", sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_filtered_sorted.bed.gz", sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_filtered_sorted.bw", sample=SAMPLES),
        # ── outputs from Unphased_modkit_5mC_pileupBedgraph ────────────────
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined.bedgraph",sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered.bedgraph",sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz",sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz.tbi",sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz",sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz.tbi", sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered.bw",sample=SAMPLES),
        expand(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/done.{{sample}}.txt", outdir=OUTDIR, sample=SAMPLES),
        expand(f"{OUTDIR}/bamCoverage_toBigWig/{{sample}}/{{sample}}.bw", sample=SAMPLES),








rule modkit_5mC_5hmC_pileup_Bed:
    input:
        haplotagged_bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        outdir= directory(f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/"),
        mp_phased1_Bed= f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_1.bed",
        mp_phased2_Bed= f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_2.bed",
        mp_phasedUngrouped_Bed= f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_ungrouped.bed",
    singularity: MODKIT_IMG
    params:
        ref= REF,
        mincov= minCov
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{wildcards.sample}
        modkit pileup \
        --ref {params.ref} \
        --cpg --combine-strands \
        --filter-threshold 0.7 \
        --threads {threads} \
        --prefix {wildcards.sample} \
        --partition-tag HP \
        {input.haplotagged_bam} {output.outdir}
        """
        # --region "chr8:50000-90000124" \

rule modkit_5mC_pileup_Bed:
    input:
        haplotagged_bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        outdir= directory(f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/"),
        mp_phased1_Bed= f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_1.bed",
        mp_phased2_Bed= f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_2.bed",
        mp_phasedUngrouped_Bed= f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_ungrouped.bed",
    singularity: MODKIT_IMG
    params:
        ref= REF,
        mincov= minCov
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/modkit_5mC_pileup_Bed/{wildcards.sample}
        modkit pileup \
        --ref {params.ref} \
        --cpg --combine-strands \
        --ignore h \
        --filter-threshold 0.7 \
        --threads {threads} \
        --prefix {wildcards.sample} \
        --partition-tag HP \
        {input.haplotagged_bam} {output.outdir}
        """ 

rule Unphased_modkit_5mC_pileup:
    input:
        haplotagged_bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        mp_unphased_Bed= f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased.bed",
        mp_unphased_sortedBedgz= f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_sorted.bed.gz",
        filtered_mp_unphased_sortedBed= f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_filtered_sorted.bed",
        filtered_mp_unphased_sortedBedgz= f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_filtered_sorted.bed.gz",
        filtered_mp_unphased_bw= f"{OUTDIR}/Unphased_modkit_5mC_pileup/{{sample}}/{{sample}}_unphased_filtered_sorted.bw"

    singularity: MODKIT_IMG
    params:
        ref= REF,
        mincov= minCov,
        genomeSizes = GENOMESIZES,
        setSeed = SEED
    threads: THREADS
    log:
        "logs/Unphased_modkit_5mC_pileup/{sample}/{sample}_unphased.log"
    shell:
        """
        mkdir -p {OUTDIR}/Unphased_modkit_5mC_pileup/{wildcards.sample}

        modkit pileup \
        --ref {params.ref} \
        --cpg --combine-strands \
        --ignore h \
        --sampling-frac 0.8 \
        --filter-threshold 0.7 \
        --threads {threads} \
        --prefix {wildcards.sample} \
        {input.haplotagged_bam} {output.mp_unphased_Bed}

        sort -k1,1 -k2,2n {output.mp_unphased_Bed} | bgzip -c > {output.mp_unphased_sortedBedgz} 
        sort -k1,1 -k2,2n {output.mp_unphased_Bed} | awk -v min_cov={params.mincov} '$5 >= min_cov' > {output.filtered_mp_unphased_sortedBed}
        
        bgzip -c {output.filtered_mp_unphased_sortedBed} > {output.filtered_mp_unphased_sortedBedgz}
        tabix -p bed {output.filtered_mp_unphased_sortedBedgz}
        tabix -p bed {output.mp_unphased_sortedBedgz}
        
        echo "Converting {output.filtered_mp_unphased_sortedBedgz} to bedmethyl format"
        modkit bedmethyl tobigwig --mod-codes m \
        --suppress-progress \
        --nthreads {threads} \
        --sizes {params.genomeSizes} \
        --log-filepath {log} \
          {output.filtered_mp_unphased_sortedBed} {output.filtered_mp_unphased_bw}
        """ 


rule Unphased_modkit_5mC_pileupBedgraph:
    input:
        haplotagged_bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        outdir= directory(f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/"),
        # raw bedGraph
        bedgraph                   = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined.bedgraph",
        # filtered bedGraph
        bedgraph_filtered          = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered.bedgraph",
        # sorted+compressed & index: raw
        bedgraph_sorted            = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz",
        bedgraph_sorted_tbi        = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_sorted.bedgraph.gz.tbi",
        # sorted+compressed & index: filtered
        bedgraph_filt_sorted       = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz",
        bedgraph_filt_sorted_tbi   = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered_sorted.bedgraph.gz.tbi",
        # BigWig from filtered
        bw                         = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/{{sample}}_unphased_m_CG0_combined_filtered.bw",
        # done‐marker
        done                       = f"{OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{{sample}}/done.{{sample}}.txt"
    params:
        ref         = REF,
        mincov      = minCov,
        setSeed     = SEED,
        genomeSizes = GENOMESIZES
    threads: THREADS
    log:
        "logs/Unphased_modkit_5mC_pileupBedgraph/{sample}/{sample}.log"
    singularity: MODKIT_IMG
    shell:
        r"""
        mkdir -p {OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{wildcards.sample}
        # mkdir -p {OUTDIR}/Unphased_modkit_5mC_pileupBedgraph/{wildcards.sample}

        # 1) raw bedGraph
        modkit pileup \
            --ref {params.ref} \
            --cpg --combine-strands \
            --ignore h \
            --bedgraph \
            --filter-threshold 0.7 \
            --sampling-frac 0.8 \
            --seed {params.setSeed} \
            --threads {threads} \
            --prefix {wildcards.sample}_unphased \
            {input.haplotagged_bam} {output.outdir}

        # 2) sort & bgzip + tabix (raw)
        sort -k1,1 -k2,2n {output.bedgraph} \
            | bgzip -c > {output.bedgraph_sorted}
        tabix -p bed {output.bedgraph_sorted}

        # 3) filter by coverage
        awk -v min_cov={params.mincov} '$5 >= min_cov' {output.bedgraph} \
            > {output.bedgraph_filtered}

        # 4) sort & bgzip + tabix (filtered)
        sort -k1,1 -k2,2n {output.bedgraph_filtered} \
            | bgzip -c > {output.bedgraph_filt_sorted}
        tabix -p bed {output.bedgraph_filt_sorted}

        # 5) make BigWig from filtered
        /data1/greenbab/users/ahunos/apps/ucsc_tools/bedGraphToBigWig \
            {output.bedgraph_filt_sorted} \
            {params.genomeSizes} \
            {output.bw}

        # 6) mark done
        touch {output.done}
        """


rule highQual_modkit_pileup_5mC_sortTabixBed:
    input:
        bed= f"{OUTDIR}/modkit_5mC_pileup_Bed/{{sample}}/{{sample}}_{{haplotypes}}.bed"
    output:
        sorted_bed= f"{OUTDIR}/highQual_modkit_pileup_5mC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed",
        sorted_tabixed= f"{OUTDIR}/highQual_modkit_pileup_5mC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed.gz"
    singularity: MODKIT_IMG
    params:
        mincov= minCov
    threads: THREADS
    shell:
        """
        echo "Sorting and filtering {input.bed} for min coverage {params.mincov}"
        sort -k1,1 -k2,2n {input.bed} | awk -v min_cov={params.mincov} '$5 >= min_cov' > {output.sorted_bed} 
        bgzip -c {output.sorted_bed} > {output.sorted_tabixed}
        tabix -p bed {output.sorted_tabixed}
        """       
        # --region "chr19:50000-9000124" \

rule highQual_modkit_pileup_5mC_5hmC_sortTabixBed:
    input:
        bed= f"{OUTDIR}/modkit_5mC_5hmC_pileup_Bed/{{sample}}/{{sample}}_{{haplotypes}}.bed"
    output:
        sorted_bed= f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed",
        sorted_tabixed= f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed.gz"
    singularity: MODKIT_IMG
    params:
        mincov= minCov
    threads: THREADS
    shell:
        """
        echo "Sorting and filtering {input.bed} for min coverage {params.mincov}"
        sort -k1,1 -k2,2n {input.bed} | awk -v min_cov={params.mincov} '$5 >= min_cov' > {output.sorted_bed} 
        bgzip -c {output.sorted_bed} > {output.sorted_tabixed}
        tabix -p bed {output.sorted_tabixed}
        """


rule highQual_modkit_pileup_5mC_5hmC_BigWigs:
    input:
        bed= f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_sortTabixBed/{{sample}}/{{sample}}_{{haplotypes}}_sorted.bed"
    output:
        bedmethyl= f"{OUTDIR}/highQual_modkit_pileup_5mC_5hmC_BigWigs/{{sample}}/{{sample}}_{{haplotypes}}_{{mod}}_sorted.bw"
    singularity: MODKIT_IMG
    params:
        ref= REF,
        genomeSizes = GENOMESIZES
    threads: THREADS
    log:
        "logs/highQual_modkit_pileup_5mC_5hmC_BigWigs/{sample}/{sample}_{haplotypes}_{mod}_sorted.log"
    shell:
        """
        echo "Converting {input.bed} to bedmethyl format"
        modkit bedmethyl tobigwig --mod-codes {wildcards.mod} \
        --suppress-progress \
        --nthreads {threads} \
        --sizes {params.genomeSizes} \
        --log-filepath {log} \
          {input.bed} {output.bedmethyl}
        """
        # --negative-strand-values \


rule bamCoverage_toBigWig:
    input:
        haplotagged_bam=lambda wildcards: config["samples"][wildcards.sample]
    output:
        bamCov_bw= f"{OUTDIR}/bamCoverage_toBigWig/{{sample}}/{{sample}}.bw"
    params:
        ref= REF,
        genomeSizes = GENOMESIZES
    threads: THREADS
    log:
        "logs/bamCoverage_toBigWig/{sample}/{sample}.log"
    shell:
        r"""
        source /home/ahunos/miniforge3/etc/profile.d/conda.sh && conda activate methyl_ONT

        bamCoverage -b {input.haplotagged_bam} \
            -o {output.bamCov_bw} \
            --binSize 10 \
            --normalizeUsing RPKM \
            --smoothLength 50 \
            --extendReads 150 \
            --centerReads -p 6
        """
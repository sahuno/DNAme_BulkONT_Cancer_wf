#!/usr/bin/env python3
"""
plot_methylation_histograms.py

Read multiple genome-wide CpG methylation bedgraph/TSV files and produce a multi-panel
histogram of methylation proportions for each file, with optional chromosome and coverage filtering.

Bins default to 5% increments (0.05) across [0,1] (20 bins) to match published figures.
Coverage filter defaults to >=5 reads per CpG.

Input files must be tab-delimited (bedgraph standard) with columns:
  chrom  start  end  methylation  coverage

Dependencies:
  - polars (>=0.18)
  - matplotlib
"""
import argparse
import math
from pathlib import Path

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator


def parse_args():
    parser = argparse.ArgumentParser(
        description="Multi-panel histograms of genome-wide methylation proportions."
    )
    parser.add_argument(
        "-i", "--inputs", nargs='+', required=True,
        help="Input bedgraph/TSV files (tab-delimited) of CpG methylation."
    )
    parser.add_argument(
        "-o", "--output", default="methylation_histograms.png",
        help="Output image file (PNG, PDF, etc.)."
    )
    parser.add_argument(
        "-b", "--bins", type=int, default=20,
        help="Number of bins for the histogram (default 20 for 5% increments)."
    )
    parser.add_argument(
        "-c", "--cols", type=int, default=2,
        help="Number of columns in the multi-panel grid."
    )
    parser.add_argument(
        "-r", "--chromosomes", nargs='+',
        help="List of chromosomes to include (e.g. chr1 chr2). Default: all."
    )
    parser.add_argument(
        "-m", "--min-coverage", type=int, default=5,
        help="Minimum coverage (reads) per site to include. Default: 5."
    )
    return parser.parse_args()


def load_methylation(path: Path,
                     chroms: list[str] | None = None,
                     min_cov: int = 5) -> pl.Series:
    # Read tab-delimited bedgraph
    df = pl.read_csv(
        path,
        separator='\t',
        has_header=False,
        new_columns=["chrom", "start", "end", "methylation", "coverage"],
        schema_overrides={
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "methylation": pl.Float64,
            "coverage": pl.Int64,
        }
    )
    # Filter by chromosomes
    if chroms:
        df = df.filter(pl.col("chrom").is_in(chroms))
    # Filter by coverage
    if min_cov > 0:
        df = df.filter(pl.col("coverage") >= min_cov)
    return df["methylation"]


def main():
    args = parse_args()
    files = [Path(p) for p in args.inputs]
    n = len(files)

    # Grid layout
    cols = args.cols
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols,
                             figsize=(4 * cols, 3 * rows),
                             sharex=True, sharey=True)
    axes = axes.flatten() if n > 1 else [axes]

    # Plot histograms using NumPy arrays for speed
    for ax, fp in zip(axes, files):
        meth_values = load_methylation(fp,
                                       args.chromosomes,
                                       args.min_coverage).to_numpy()
        counts, bins, patches = ax.hist(
            meth_values,
            bins=args.bins,
            range=(0.0, 1.0),
            alpha=0.7,
            edgecolor='black'
        )
        # Format y-axis in thousands
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{int(y/1000)}"))

        ax.set_title(fp.stem)
        ax.set_xlabel("Methylation Proportion")
        ax.set_ylabel("No. of CpG units (×1,000)")
        ax.set_xlim(0, 1)

    # Hide unused panels
    for ax in axes[n:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()


# 1) Make sure dependencies are installed
# pip install polars matplotlib

# 2) Run the script on two bedgraph files, plotting only chr1 and chr2,
#    arranging panels in 2 columns, using the default 20 bins, and
#    writing to “out.png”:
# python /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/scripts/genome_wide_methylation.py \
#   -i /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-0-1_5000_4000/D-0-1_5000_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-0-2_5000_4000/D-0-2_5000_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-0-3_4000/D-0-3_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-C-1_4000/D-C-1_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-C-2_4000/D-C-2_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-C-3_4000/D-C-3_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-QC-1_4000/D-QC-1_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-QC-2_4000/D-QC-2_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-QC-3_4000/D-QC-3_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-Q-1_5000_4000/D-Q-1_5000_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-Q-2_4000/D-Q-2_4000_unphased_m_CG0_combined.bedgraph \
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-Q-3_4000/D-Q-3_4000_unphased_m_CG0_combined.bedgraph \
#   -c 3 \
#   -b 20 \
#   -m 5 \
#   -r chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM \
#   -o genomeWideDNAme_StandardChrom_QSAT_CKI_mono_combo_minCov5.png

#   -r chr19 \
# chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM

# results/Unphased_modkit_5mC_pileupBedgraph/D-Q-2_4000/D-Q-2_4000_unphased_m_CG0_combined_filtered.bedgraph
# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-Q-2_4000/D-Q-2_4000_unphased_m_CG0_combined_filtered.bedgraph

# /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/*/*.bedgraph
#!/usr/bin/env python3
"""
plot_methylation_bins.py

Read genome-wide CpG methylation bedgraph/TSV files and plot a summary statistic per
variable-size genomic bins, using a flexible chromosome-size file for any reference.

Input methylation files must be tab-delimited with columns:
  chrom  start  end  methylation  coverage

Chromosome sizes must be provided as a two-column TSV (no header):
  chrom\tsize

Dependencies:
  - polars (>=0.18)
  - matplotlib

Usage:
  python plot_methylation_bins.py \
    --chrom-sizes chr_sizes.tsv \
    -i sample1.bedgraph sample2.bedgraph \
    -o methylation_bins.png \
    -r chr1 chr2 \
    -m 10 \
    -b 200000 \
    -s median

Statistic options: mean, median, geomean (geometric mean).
"""
import argparse
from pathlib import Path
import math

import polars as pl
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot methylation statistic per genomic bin with flexible chrom sizes."
    )
    p.add_argument(
        "--chrom-sizes", required=True,
        help="TSV file of chromosome sizes: two columns (chrom, size) with no header."
    )
    p.add_argument(
        "-i", "--inputs", nargs='+', required=True,
        help="Input bedgraph files (tab-delimited) with methylation data."
    )
    p.add_argument(
        "-o", "--output", default="methylation_bins.png",
        help="Output image file (PNG, PDF, etc.)."
    )
    p.add_argument(
        "-r", "--chromosomes", nargs='+',
        help="Chromosomes to include (default: all from chrom-sizes file)."
    )
    p.add_argument(
        "-m", "--min-coverage", type=int, default=5,
        help="Minimum coverage per CpG (default: 5)."
    )
    p.add_argument(
        "-b", "--bin-size", type=int, default=100000,
        help="Bin size in base pairs (default: 100000 = 100kb)."
    )
    p.add_argument(
        "-s", "--statistic", choices=["mean","median","geomean"],
        default="mean",
        help="Statistic to compute per bin: mean, median, geomean."
    )
    return p.parse_args()


def load_chrom_sizes(path: str) -> dict[str,int]:
    df = pl.read_csv(
        path,
        separator='\t',
        has_header=False,
        new_columns=["chrom","size"],
        schema_overrides={"chrom": pl.Utf8, "size": pl.Int64}
    )
    return dict(zip(df["chrom"].to_list(), df["size"].to_list()))


def load_and_filter(path: Path, chroms: list[str] | None, min_cov: int) -> pl.DataFrame:
    df = pl.read_csv(
        path,
        separator='\t',
        has_header=False,
        new_columns=["chrom","start","end","methylation","coverage"],
        schema_overrides={
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "methylation": pl.Float64,
            "coverage": pl.Int64,
        }
    )
    if chroms:
        df = df.filter(pl.col("chrom").is_in(chroms))
    return df.filter(pl.col("coverage") >= min_cov)


def main():
    args = parse_args()
    # Load chromosome sizes from file
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    # Determine which chromosomes to plot
    chrom_list = args.chromosomes or list(chrom_sizes.keys())

    # Compute cumulative offsets
    offsets: dict[str,int] = {}
    cum = 0
    for c in chrom_list:
        size = chrom_sizes.get(c)
        if size is None:
            raise ValueError(f"Chromosome '{c}' not found in sizes file.")
        offsets[c] = cum
        cum += size

    plt.figure(figsize=(12, 4))

    # Plot each sample
    for path in args.inputs:
        df = load_and_filter(Path(path), chrom_list, args.min_coverage)
        df = df.with_columns(
            (pl.col("start") + pl.col("chrom").map_dict(offsets)).alias("cum_pos")
        )
        df = df.with_columns((pl.col("cum_pos") // args.bin_size).alias("bin_idx"))
        # Choose aggregation
        if args.statistic == "mean":
            agg_expr = pl.col("methylation").mean().alias("stat")
        elif args.statistic == "median":
            agg_expr = pl.col("methylation").median().alias("stat")
        else:  # geomean
            agg_expr = (pl.col("methylation").log().mean().exp()).alias("stat")
        summary = (
            df.groupby("bin_idx")
              .agg(agg_expr)
              .sort("bin_idx")
        )
        centers = summary["bin_idx"].to_numpy() * args.bin_size + args.bin_size/2
        plt.plot(centers, summary["stat"].to_numpy(), label=Path(path).stem)

    # Draw chromosome boundaries and set x-ticks
    for c in chrom_list:
        plt.axvline(offsets[c], color='gray', lw=0.5)
    mids = [offsets[c] + chrom_sizes[c]/2 for c in chrom_list]
    plt.xticks(mids, chrom_list, rotation=90, fontsize=8)

    # Labels and legend
    plt.xlabel('Chromosome')
    ylabel_map = {
        'mean': 'Mean 5-mCpG rate',
        'median': 'Median 5-mCpG rate',
        'geomean': 'Geometric mean 5-mCpG rate'
    }
    plt.ylabel(ylabel_map[args.statistic])
    plt.title(f"{ylabel_map[args.statistic]} per {args.bin_size//1000}kb bins")
    plt.ylim(0, 1)
    plt.legend(fontsize=6, ncol=2)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)


if __name__ == '__main__':
    main()


# python /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/scripts/plot_methylation_bins.py \
# --chrom-sizes /data1/greenbab/database/mm10/mm10.sorted_k11_k22n.standard_chrom.sizes \
# -i /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/methylTables/rerun_same_pileup_seed/results/Unphased_modkit_5mC_pileupBedgraph/D-0-1_5000_4000/D-0-1_5000_4000_unphased_m_CG0_combined.bedgraph \
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
#   -r chr1 chr2 chrX chr3 chr4 chr5 chr6 chr7 chr10 chr8 chr14 chr9 chr11 chr13 chr12 chr15 chr16 chr17 chrY chr18 chr19 chrM \
#   -m 5 \
#   -b 100000 \
#   -s median \
#   -o methylation_100kb_median.png


# /data1/greenbab/database/mm10/mm10.sorted_k11_k22n.standard_chrom.sizes
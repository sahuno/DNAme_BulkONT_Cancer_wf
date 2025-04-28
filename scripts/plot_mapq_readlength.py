#!/usr/bin/env python3

# /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/scripts/plot_mapq_readlength.py -o mapq_vs_len_D01_5000_4000.png /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam --hex --log-length

#  --max-reads 50000

import argparse
import sys

import pysam
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot mapping quality vs. read length from one or more BAM files."
    )
    p.add_argument(
        "bams", nargs="+",
        help="Input BAM file(s)."
    )
    p.add_argument(
        "-o", "--output", default="mapq_vs_length.png",
        help="Output image file (png, pdf, etc.)."
    )
    p.add_argument(
        "--hex", action="store_true",
        help="Use hexbin density plot instead of scatter."
    )
    p.add_argument(
        "--max-reads", type=int, default=None,
        help="If set, randomly sample up to this many reads per BAM (for speed)."
    )
    p.add_argument(
        "--log-length", action="store_true",
        help="Plot the log10 of read length instead of the raw read length."
    )
    return p.parse_args()

def collect_metrics(bam_path, max_reads=None):
    """
    Iterate through the BAM, return a DataFrame with columns:
      - mapq
      - length
      - sample (basename of bam)
    """
    sample = bam_path.rsplit("/", 1)[-1]
    bam = pysam.AlignmentFile(bam_path, "rb")
    records = []
    for i, read in enumerate(bam.fetch(until_eof=True)):
        if read.is_unmapped:
            continue
        records.append((read.mapping_quality, read.query_length))
        if max_reads and i >= max_reads:
            break
    bam.close()
    df = pd.DataFrame(records, columns=["mapq", "length"])
    df["sample"] = sample
    return df

def main():
    args = parse_args()

    # Gather data from all BAMs
    dfs = []
    for bam in args.bams:
        try:
            df = collect_metrics(bam, max_reads=args.max_reads)
            dfs.append(df)
        except Exception as e:
            sys.stderr.write(f"Error processing {bam}: {e}\n")
    if not dfs:
        sys.exit("No data collected—check your BAM files.")

    data = pd.concat(dfs, ignore_index=True)

    # Decide on X axis transformation & label
    if args.log_length:
        # avoid log10(0)
        data["length"] = np.log10(data["length"].clip(lower=1))
        xlabel = "Log10(Read length) (bp)"
    else:
        xlabel = "Read length (bp)"

    # 1) Create figure & axes
    fig, ax = plt.subplots(figsize=(8,6))

    # 2) Plot
    if args.hex:
        hb = ax.hexbin(
            data["length"], data["mapq"],
            gridsize=100, mincnt=1
        )
        fig.colorbar(hb, ax=ax, label="Read count")
    else:
        ax.scatter(
            data["length"], data["mapq"],
            s=1, alpha=0.3, rasterized=True
        )

    # 3) Now set labels & title
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Mapping quality (MAPQ)")
    ax.set_title("MAPQ vs. Read length")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved plot to {args.output}")

if __name__ == "__main__":
    main()
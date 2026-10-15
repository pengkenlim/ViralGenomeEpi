#!/usr/bin/env python3
"""Simple pairwise SNP counting for a trimmed MSA FASTA.

This is analogous to module 8 but deliberately excludes gene-aware SNP calls.
Edit the constants at the top before running.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO


INPUT_FASTA = "/media/mngs/48TBRAID5HDD/20261008_230815_Renamed_Out_2/MSA_Output/Consensus_plus_denovo_genomes_aligned_trimmed.fasta"
OUTPUT_DIR = "/media/mngs/48TBRAID5HDD/20261008_230815_Renamed_Out_2/MSA_Output"
DRAW_HISTOGRAM = False


def compute_pairwise_snp_matrix(alignment_file):
    """Compute pairwise SNP distances across an aligned FASTA."""
    alignment = AlignIO.read(alignment_file, "fasta")
    headers = [rec.id for rec in alignment]
    n_seqs = len(headers)

    seq_matrix = np.array([list(str(rec.seq).upper()) for rec in alignment])
    dist_matrix = np.zeros((n_seqs, n_seqs), dtype=int)
    pair_rows = []

    canonical_bases = np.array(["A", "C", "G", "T"])

    for i, j in __import__("itertools").combinations(range(n_seqs), 2):
        seq1 = seq_matrix[i]
        seq2 = seq_matrix[j]
        comparable_mask = np.isin(seq1, canonical_bases) & np.isin(seq2, canonical_bases)
        snp_diff = int(np.sum((seq1 != seq2) & comparable_mask))

        dist_matrix[i, j] = snp_diff
        dist_matrix[j, i] = snp_diff
        pair_rows.append({"Sample_1": headers[i], "Sample_2": headers[j], "SNP_Distance": snp_diff})

    df_matrix = pd.DataFrame(dist_matrix, index=headers, columns=headers)
    df_pairs = pd.DataFrame(pair_rows)
    return df_matrix, df_pairs


def write_outputs(df_matrix, df_pairs, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    matrix_path = os.path.join(output_dir, "pairwise_snp_matrix.csv")
    df_matrix.to_csv(matrix_path)
    print(f"[MATRIX]: wrote {matrix_path}")

    pair_path = os.path.join(output_dir, "pairwise_snp_pairs.csv")
    df_pairs.to_csv(pair_path, index=False)
    print(f"[PAIRS]: wrote {pair_path}")

    if DRAW_HISTOGRAM:
        max_snps = int(df_pairs["SNP_Distance"].max()) if not df_pairs.empty else 10
        bins = np.arange(-0.5, max_snps + 1.5, 1)

        plt.figure(figsize=(8, 5))
        plt.hist(df_pairs["SNP_Distance"], bins=bins, color="#2b5c8f", edgecolor="black", alpha=0.8)
        plt.title("Distribution of Pairwise SNP Differences", fontsize=12, fontweight="bold")
        plt.xlabel("Pairwise SNP Difference", fontsize=10)
        plt.ylabel("Number of Sample Pairs", fontsize=10)
        plt.grid(axis="y", linestyle=":", alpha=0.6)
        plt.tight_layout()

        hist_path = os.path.join(output_dir, "pairwise_snp_histogram.png")
        plt.savefig(hist_path, dpi=300)
        plt.close()
        print(f"[HIST]: wrote {hist_path}")


def main():
    if not os.path.exists(INPUT_FASTA):
        raise FileNotFoundError(f"Input FASTA not found: {INPUT_FASTA}")

    df_matrix, df_pairs = compute_pairwise_snp_matrix(INPUT_FASTA)
    write_outputs(df_matrix, df_pairs, OUTPUT_DIR)

    print("\nDone.")
    print(df_pairs.to_string(index=False))


if __name__ == "__main__":
    main()

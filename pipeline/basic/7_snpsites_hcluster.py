import csv
import glob
import itertools
import os
import subprocess
from Bio import AlignIO, SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
import seaborn as sns

# ==========================================
# User Configuration
# ==========================================
#job_name = "RSV_A"
#accession_list = [
#    "SRR34884131",
#    "SRR34884130",
#    "SRR34884129",
#    "SRR34884128",
#    "SRR34884122",
#    "SRR34884115",
#    "SRR34884239",
#    "SRR34884236",
#    "SRR34884235",
#    "SRR34884234",
#    "SRR34884233",
#    "SRR34884232",
#]

job_name = "RSV_B"
accession_list = [
    "SRR34884179",
    "SRR34884178",
    "SRR34884145",
    "SRR34884143",
    "SRR34884142",
    "SRR34884141",
    "SRR34884140",
    "SRR34884139",
    "SRR34884138",
    "SRR34884137",
    "SRR34884136",
    "SRR34884135",
    "SRR34884134",
    "SRR34884127",
    "SRR34884126",
    "SRR34884125",
    "SRR34884124",
    "SRR34884123",
    "SRR34884120",
    "SRR34884119",
    "SRR34884118",
    "SRR34884117",
    "SRR34884116",
    "SRR34884114",
    "SRR34884241",
    "SRR34884240",
    "SRR34884237",
]

phylogeny_output_dir = (
    "/media/mngs/48TBRAID5HDD/viral_3/phylogeny_results"
)
trimmed_aligned_fasta = os.path.join(
    phylogeny_output_dir, "genomes_aligned_trimmed.fasta"
)
cluster_output_dir = os.path.join(
    phylogeny_output_dir, f"cluster_analysis_{job_name}"
)

# Cluster parameters
snp_cutoff = 3  # Max SNP distance to define outbreak clusters
linkage_method = "average"  # Average linkage clustering as specified in methods
snpsites_bin = "snp-sites"  # Path to snp-sites binary


def subset_alignment_for_accessions(input_fasta, accessions, out_fasta):
    """Subset the alignment to just the requested accessions."""
    records = list(AlignIO.read(input_fasta, "fasta"))
    desired = set(accessions)
    filtered = [rec for rec in records if rec.id in desired]

    missing = sorted(desired - {rec.id for rec in filtered})
    if missing:
        print(
            f"[WARNING]: Missing records in alignment for accessions: {missing}"
        )

    if not filtered:
        raise ValueError(
            f"No records found in {input_fasta} for the requested accession list."
        )

    with open(out_fasta, "w") as handle:
        SeqIO.write(filtered, handle, "fasta")

    print(f"[SUBSET]: Wrote {len(filtered)} sequences to {out_fasta}")
    return out_fasta


def run_snpsites(input_fasta, out_dir):
    """Run snp-sites to extract variant positions into a VCF and SNP alignment."""
    print(
        "[SNP-SITES]: Extracting variant positions from multiple alignment..."
    )
    vcf_out = os.path.join(out_dir, "variants.vcf")
    fasta_out = os.path.join(out_dir, "snp_alignment.fasta")

    # Generate VCF
    cmd_vcf = [snpsites_bin, "-v", "-o", vcf_out, input_fasta]
    subprocess.run(
        cmd_vcf,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )

    # Generate Multi-FASTA of SNP sites
    cmd_fasta = [snpsites_bin, "-m", "-o", fasta_out, input_fasta]
    subprocess.run(
        cmd_fasta,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
    )

    print(
        f"[SNP-SITES]: Variant files generated:\n  - VCF: {vcf_out}\n  - SNP Alignment: {fasta_out}"
    )
    return vcf_out, fasta_out


def compute_pairwise_snp_matrix(alignment_file):
    """Calculates pairwise SNP distances across positions where both samples have canonical ACGT bases."""
    print(
        "[SNP MATRIX]: Reading alignment and calculating pairwise SNP distances..."
    )
    alignment = AlignIO.read(alignment_file, "fasta")
    headers = [rec.id for rec in alignment]
    n_seqs = len(headers)

    # Convert alignment to 2D numpy array of uppercase characters
    seq_matrix = np.array([list(str(rec.seq).upper()) for rec in alignment])

    dist_matrix = np.zeros((n_seqs, n_seqs), dtype=int)
    pair_list = []

    # Strict canonical bases definition
    canonical_bases = np.array(["A", "C", "G", "T"])

    for i, j in itertools.combinations(range(n_seqs), 2):
        seq1 = seq_matrix[i]
        seq2 = seq_matrix[j]

        # A site is comparable ONLY if BOTH samples contain a definitive ACGT base
        is_canonical_1 = np.isin(seq1, canonical_bases)
        is_canonical_2 = np.isin(seq2, canonical_bases)
        comparable_mask = is_canonical_1 & is_canonical_2

        # Count SNP differences ONLY where both bases are canonical AND non-matching
        snp_diff = int(np.sum((seq1 != seq2) & comparable_mask))

        dist_matrix[i, j] = snp_diff
        dist_matrix[j, i] = snp_diff

        pair_list.append(
            {
                "Sample_1": headers[i],
                "Sample_2": headers[j],
                "SNP_Distance": snp_diff,
            }
        )

    df_matrix = pd.DataFrame(dist_matrix, index=headers, columns=headers)
    df_pairs = pd.DataFrame(pair_list)

    return df_matrix, df_pairs


def plot_snp_histogram(df_pairs, out_dir):
    """Generates a histogram distribution of pairwise SNP distances."""
    print("[PLOTTING]: Generating pairwise SNP distance histogram...")
    plt.figure(figsize=(8, 5))

    max_snps = int(df_pairs["SNP_Distance"].max()) if not df_pairs.empty else 10
    bins = np.arange(-0.5, max_snps + 1.5, 1)

    sns.histplot(
        df_pairs["SNP_Distance"],
        bins=bins,
        color="#2b5c8f",
        edgecolor="black",
        alpha=0.8,
    )

    plt.axvline(
        x=snp_cutoff,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Cluster Cutoff ({snp_cutoff} SNPs)",
    )

    plt.title(
        "Distribution of Pairwise SNP Differences",
        fontsize=12,
        fontweight="bold",
    )
    plt.xlabel("Pairwise SNP Difference", fontsize=10)
    plt.ylabel("Number of Sample Pairs", fontsize=10)
    plt.grid(axis="y", linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    fig_path = os.path.join(out_dir, "pairwise_snp_histogram.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[PLOTTING]: Histogram saved to {fig_path}")


def perform_hierarchical_clustering(df_matrix, out_dir):
    """Performs average-linkage hierarchical clustering and plots the dendrogram."""
    print(
        f"[CLUSTERING]: Conducting hierarchical clustering ({linkage_method} linkage)..."
    )

    # Convert square matrix to condensed 1D vector for SciPy
    condensed_dist = squareform(df_matrix.values, force="tovector")

    # Perform linkage
    Z = linkage(condensed_dist, method=linkage_method)

    # Plot Dendrogram
    plt.figure(figsize=(10, 6))
    dendrogram(
        Z,
        labels=df_matrix.index,
        leaf_rotation=90,
        leaf_font_size=10,
        color_threshold=snp_cutoff,
    )

    plt.axhline(
        y=snp_cutoff,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Cutoff = {snp_cutoff} SNPs",
    )

    plt.title(
        f"Hierarchical Clustering Dendrogram ({linkage_method.capitalize()} Linkage)",
        fontsize=12,
        fontweight="bold",
    )
    plt.xlabel("Sample Accession", fontsize=10)
    plt.ylabel("Pairwise SNP Distance", fontsize=10)
    plt.legend(loc="upper right")
    plt.tight_layout()

    dendrogram_path = os.path.join(
        out_dir, "hierarchical_clustering_dendrogram.png"
    )
    plt.savefig(dendrogram_path, dpi=300)
    plt.close()
    print(f"[PLOTTING]: Dendrogram saved to {dendrogram_path}")

    # Assign cluster IDs using fcluster
    cluster_ids = fcluster(Z, t=snp_cutoff, criterion="distance")

    df_clusters = pd.DataFrame(
        {
            "Sample_ID": df_matrix.index,
            "Cluster_ID": [f"Cluster_{cid}" for cid in cluster_ids],
        }
    )

    return df_clusters


def main():
    os.makedirs(cluster_output_dir, exist_ok=True)

    print("\n==========================================")
    print(
        f" Step 7: Pairwise SNP & Hierarchical Clustering Analysis for job '{job_name}'"
    )
    print("==========================================")
    print(f"[JOB]: job_name = {job_name}")
    print(f"[JOB]: accession_list = {accession_list}")
    print(f"[JOB]: cluster_output_dir = {cluster_output_dir}")

    if not os.path.exists(trimmed_aligned_fasta):
        print(
            f"[FATAL]: Input alignment file not found: {trimmed_aligned_fasta}"
        )
        return

    subset_alignment = os.path.join(
        cluster_output_dir, "subset_accessions_alignment.fasta"
    )
    subset_alignment = subset_alignment_for_accessions(
        trimmed_aligned_fasta,
        accession_list,
        subset_alignment,
    )

    # 1. Run snp-sites on the subsetted alignment
    vcf_file, snp_fasta = run_snpsites(subset_alignment, cluster_output_dir)

    # 2. Calculate Pairwise SNP Matrix & Pairwise Differences
    df_matrix, df_pairs = compute_pairwise_snp_matrix(subset_alignment)

    matrix_csv = os.path.join(cluster_output_dir, "pairwise_snp_matrix.csv")
    pairs_csv = os.path.join(cluster_output_dir, "pairwise_snp_pairs.csv")
    df_matrix.to_csv(matrix_csv)
    df_pairs.to_csv(pairs_csv, index=False)
    print(
        f"[SUCCESS]: Saved SNP matrix ({matrix_csv}) and pairs list ({pairs_csv})"
    )

    # 3. Plot Histogram of Pairwise SNP Differences
    plot_snp_histogram(df_pairs, cluster_output_dir)

    # 4. Perform Hierarchical Clustering & Generate Dendrogram
    df_clusters = perform_hierarchical_clustering(
        df_matrix, cluster_output_dir
    )

    # 5. Save Cluster Assignments
    cluster_csv = os.path.join(
        cluster_output_dir, f"clusters_{snp_cutoff}snp_cutoff.csv"
    )
    df_clusters.to_csv(cluster_csv, index=False)

    print("\n==========================================")
    print(" Cluster Summary Results")
    print("==========================================")
    print(df_clusters.to_string(index=False))
    print(f"\n[COMPLETE]: All clustering outputs saved to {cluster_output_dir}")


if __name__ == "__main__":
    main()
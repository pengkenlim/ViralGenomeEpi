import argparse
import importlib.util
import itertools
import json
import os
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import AlignIO, SeqIO
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_config(config_path):
    """Load a Python config file as a Python module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def get_config_value(config, section_name, key, default=None):
    section = getattr(config, section_name, None)
    if isinstance(section, dict):
        return section.get(key, default)
    return default


def resolve_executable_path(config, key):
    value = get_config_value(config, "module_8_params", key, None)
    if value is None:
        raise FileNotFoundError(
            f"Executable '{key}' is not configured. Set an explicit path or command in module_8_params."
        )
    if os.path.exists(value):
        return value
    resolved = shutil.which(value)
    if resolved:
        return resolved
    raise FileNotFoundError(
        f"Configured executable '{value}' for '{key}' was not found. Set an explicit path in module_8_params."
    )


# ---------------------------------------------------------------------------
# Output resolution
# ---------------------------------------------------------------------------

def get_module_7_dir(config):
    analysis_name = (
        get_config_value(config, "module_7_params", "analysis_name", None)
        or get_config_value(config, "module_8_params", "analysis_name", "default")
        or "default"
    )
    return os.path.join(config.cross_module_params["output_dir"], "module_7", analysis_name)


def get_module_8_dir(config):
    analysis_name = (
        get_config_value(config, "module_8_params", "analysis_name", None)
        or get_config_value(config, "module_7_params", "analysis_name", "default")
        or "default"
    )
    return os.path.join(config.cross_module_params["output_dir"], "module_8", analysis_name)


def get_input_samples(config):
    samples = get_config_value(config, "module_8_params", "input_samples", None)
    if samples is None:
        samples = get_config_value(config, "module_7_params", "input_samples", None)
    if samples:
        return list(samples)
    return list(config.cross_module_params.get("samples", []))


def find_trimmed_alignment(config):
    module_7_dir = get_module_7_dir(config)
    if not os.path.isdir(module_7_dir):
        raise FileNotFoundError(f"Module 7 directory not found: {module_7_dir}")

    candidates = [
        os.path.join(module_7_dir, "polished_genomes_aligned_trimmed.fasta"),
        os.path.join(module_7_dir, "genomes_aligned_trimmed.fasta"),
        os.path.join(module_7_dir, "aligned_trimmed.fasta"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path

    matches = []
    for root, _, files in os.walk(module_7_dir):
        for name in files:
            if name.lower().endswith((".fa", ".fasta")) and "trim" in name.lower():
                matches.append(os.path.join(root, name))

    if matches:
        return sorted(matches)[0]

    raise FileNotFoundError(
        f"No trimmed alignment file found in module 7 output: {module_7_dir}"
    )


# ---------------------------------------------------------------------------
# SNP clustering workflow
# ---------------------------------------------------------------------------

def subset_alignment_for_accessions(input_fasta, accessions, out_fasta):
    """Subset an alignment to only the requested sample IDs."""
    records = list(AlignIO.read(input_fasta, "fasta"))
    desired = set(accessions)
    filtered = [rec for rec in records if rec.id in desired]

    missing = sorted(desired - {rec.id for rec in filtered})
    if missing:
        print(f"[WARNING]: Missing records in alignment for accessions: {missing}")

    if not filtered:
        raise ValueError(f"No records found in {input_fasta} for the requested accession list.")

    with open(out_fasta, "w") as handle:
        SeqIO.write(filtered, handle, "fasta")

    print(f"[SUBSET]: Wrote {len(filtered)} sequences to {out_fasta}")
    return out_fasta


def run_snpsites(input_fasta, out_dir, snpsites_bin):
    """Extract variant positions into VCF and SNP-only alignment."""
    print("[SNP-SITES]: Extracting variant positions from the alignment...")
    vcf_out = os.path.join(out_dir, "variants.vcf")
    fasta_out = os.path.join(out_dir, "snp_alignment.fasta")

    cmd_vcf = [snpsites_bin, "-v", "-o", vcf_out, input_fasta]
    subprocess.run(cmd_vcf, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    cmd_fasta = [snpsites_bin, "-m", "-o", fasta_out, input_fasta]
    subprocess.run(cmd_fasta, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    print(
        f"[SNP-SITES]: Variant files generated:\n  - VCF: {vcf_out}\n  - SNP Alignment: {fasta_out}"
    )
    return vcf_out, fasta_out


def compute_pairwise_snp_matrix(alignment_file):
    """Calculate pairwise SNP distances for all sample pairs."""
    print("[SNP MATRIX]: Reading alignment and calculating pairwise SNP distances...")
    alignment = AlignIO.read(alignment_file, "fasta")
    headers = [rec.id for rec in alignment]
    n_seqs = len(headers)

    seq_matrix = np.array([list(str(rec.seq).upper()) for rec in alignment])
    dist_matrix = np.zeros((n_seqs, n_seqs), dtype=int)
    pair_list = []

    canonical_bases = np.array(["A", "C", "G", "T"])

    for i, j in itertools.combinations(range(n_seqs), 2):
        seq1 = seq_matrix[i]
        seq2 = seq_matrix[j]

        is_canonical_1 = np.isin(seq1, canonical_bases)
        is_canonical_2 = np.isin(seq2, canonical_bases)
        comparable_mask = is_canonical_1 & is_canonical_2

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


def plot_snp_histogram(df_pairs, out_dir, snp_cutoff):
    """Generate a histogram of pairwise SNP distances."""
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

    plt.title("Distribution of Pairwise SNP Differences", fontsize=12, fontweight="bold")
    plt.xlabel("Pairwise SNP Difference", fontsize=10)
    plt.ylabel("Number of Sample Pairs", fontsize=10)
    plt.grid(axis="y", linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    fig_path = os.path.join(out_dir, "pairwise_snp_histogram.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[PLOTTING]: Histogram saved to {fig_path}")


def perform_hierarchical_clustering(df_matrix, out_dir, snp_cutoff, linkage_method):
    """Cluster samples using average-linkage distance clustering."""
    print(f"[CLUSTERING]: Conducting hierarchical clustering ({linkage_method} linkage)...")

    condensed_dist = squareform(df_matrix.values, force="tovector")
    Z = linkage(condensed_dist, method=linkage_method)

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

    dendrogram_path = os.path.join(out_dir, "hierarchical_clustering_dendrogram.png")
    plt.savefig(dendrogram_path, dpi=300)
    plt.close()
    print(f"[PLOTTING]: Dendrogram saved to {dendrogram_path}")

    cluster_ids = fcluster(Z, t=snp_cutoff, criterion="distance")
    df_clusters = pd.DataFrame(
        {
            "Sample_ID": df_matrix.index,
            "Cluster_ID": [f"Cluster_{cid}" for cid in cluster_ids],
        }
    )
    return df_clusters


def write_analysis_metadata(output_dir, samples, analysis_name):
    metadata_path = os.path.join(output_dir, "analysis_metadata.json")
    payload = {
        "analysis_name": analysis_name,
        "input_samples": samples,
        "module": "module_8_snpsites_hcluster",
    }
    with open(metadata_path, "w") as handle:
        json.dump(payload, handle, indent=4)
    print(f"[METADATA]: Wrote {metadata_path}")


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run SNP-site extraction and hierarchical clustering for the module 7 alignment."
    )
    parser.add_argument("--config", required=True, help="Path to the config Python file.")
    args = parser.parse_args()

    global config
    config = load_config(args.config)

    input_samples = get_input_samples(config)
    analysis_name = (
        get_config_value(config, "module_8_params", "analysis_name", None)
        or get_config_value(config, "module_7_params", "analysis_name", "default")
        or "default"
    )
    snp_cutoff = float(
        get_config_value(config, "module_8_params", "snp_cutoff", 3)
    )
    linkage_method = (
        get_config_value(config, "module_8_params", "linkage_method", "average")
        or "average"
    )
    snpsites_bin = resolve_executable_path(config, "snpsites_bin_path")

    output_dir = get_module_8_dir(config)
    os.makedirs(output_dir, exist_ok=True)

    print("\n==========================================")
    print(f" Step 8: Pairwise SNP & Hierarchical Clustering Analysis for '{analysis_name}'")
    print("==========================================")
    print(f"[JOB]: analysis_name = {analysis_name}")
    print(f"[JOB]: input_samples = {input_samples}")
    print(f"[JOB]: output_dir = {output_dir}")

    input_alignment = find_trimmed_alignment(config)
    print(f"[INFO]: Using module 7 trimmed alignment: {input_alignment}")

    subset_alignment = os.path.join(output_dir, "subset_accessions_alignment.fasta")
    subset_alignment = subset_alignment_for_accessions(
        input_alignment,
        input_samples,
        subset_alignment,
    )

    vcf_file, snp_fasta = run_snpsites(subset_alignment, output_dir, snpsites_bin)
    del vcf_file, snp_fasta

    df_matrix, df_pairs = compute_pairwise_snp_matrix(subset_alignment)

    matrix_csv = os.path.join(output_dir, "pairwise_snp_matrix.csv")
    pairs_csv = os.path.join(output_dir, "pairwise_snp_pairs.csv")
    df_matrix.to_csv(matrix_csv)
    df_pairs.to_csv(pairs_csv, index=False)
    print(f"[SUCCESS]: Saved SNP matrix ({matrix_csv}) and pairs list ({pairs_csv})")

    plot_snp_histogram(df_pairs, output_dir, int(snp_cutoff))

    df_clusters = perform_hierarchical_clustering(
        df_matrix,
        output_dir,
        int(snp_cutoff),
        linkage_method,
    )

    cluster_csv = os.path.join(output_dir, f"clusters_{int(snp_cutoff)}snp_cutoff.csv")
    df_clusters.to_csv(cluster_csv, index=False)

    write_analysis_metadata(output_dir, input_samples, analysis_name)

    print("\n==========================================")
    print(" Cluster Summary Results")
    print("==========================================")
    print(df_clusters.to_string(index=False))
    print(f"\n[COMPLETE]: All clustering outputs saved to {output_dir}")


if __name__ == "__main__":
    main()

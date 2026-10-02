import argparse
import importlib.util
import json
import os
import shutil
import subprocess

from Bio import AlignIO, SeqIO


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
    value = get_config_value(config, "module_7_params", key, None)
    if value is None:
        raise FileNotFoundError(
            f"Executable '{key}' is not configured. Set an explicit path or command in module_7_params."
        )
    if os.path.exists(value):
        return value
    resolved = shutil.which(value)
    if resolved:
        return resolved
    raise FileNotFoundError(
        f"Configured executable '{value}' for '{key}' was not found. Set an explicit path in module_7_params."
    )


# ---------------------------------------------------------------------------
# Output resolution
# ---------------------------------------------------------------------------

def get_module_6_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_6", sample_name)


def get_module_7_dir(config):
    analysis_name = (
        get_config_value(config, "module_7_params", "analysis_name", "default")
        or "default"
    )
    return os.path.join(config.cross_module_params["output_dir"], "module_7", analysis_name)


def get_input_samples(config):
    samples = get_config_value(config, "module_7_params", "input_samples", None)
    if samples:
        return list(samples)
    return list(config.cross_module_params.get("samples", []))


# ---------------------------------------------------------------------------
# Consensus discovery: locate the actual module 6 outputs we wrote
# ---------------------------------------------------------------------------

def find_polished_fasta(config, sample_name):
    sample_root = get_module_6_dir(config, sample_name)
    if not os.path.isdir(sample_root):
        raise FileNotFoundError(f"Module 6 output directory not found for {sample_name}: {sample_root}")

    candidates = [
        os.path.join(sample_root, f"{sample_name}_consensus_combined.fasta"),
        os.path.join(sample_root, sample_name, f"{sample_name}_final_polished_consensus.fa"),
        os.path.join(sample_root, sample_name, f"{sample_name}_final_polished_consensus.fasta"),
        os.path.join(sample_root, sample_name, f"{sample_name}_polished_consensus.fa"),
        os.path.join(sample_root, sample_name, f"{sample_name}_polished_consensus.fasta"),
    ]

    for path in candidates:
        if os.path.exists(path):
            return path

    # Last-resort search for any FASTA files under the sample directory
    matches = []
    for root, _, files in os.walk(sample_root):
        for name in files:
            if name.lower().endswith((".fa", ".fasta")):
                matches.append(os.path.join(root, name))

    if matches:
        return sorted(matches)[0]

    raise FileNotFoundError(f"No polished consensus FASTA found for {sample_name} in {sample_root}")


# ---------------------------------------------------------------------------
# MAFFT / FastTree workflow
# ---------------------------------------------------------------------------

def run_mafft_alignment(input_fasta, aligned_fasta, threads):
    mafft_bin = resolve_executable_path(config, "mafft_bin_path")
    print(f"[MAFFT]: Running multiple sequence alignment with {threads} threads...")
    cmd = [mafft_bin, "--auto", "--thread", str(threads), input_fasta]
    with open(aligned_fasta, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
    print(f"[MAFFT]: Alignment written to {aligned_fasta}")


def trim_alignment_ends(alignment_file, trimmed_file):
    """Trim terminal columns with gaps or Ns across all sequences."""
    print("[TRIMMING]: Removing terminal 5' and 3' alignment regions containing gaps or missing bases...")
    alignment = AlignIO.read(alignment_file, "fasta")
    num_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()

    start_col = 0
    for col in range(alignment_length):
        column_chars = [alignment[seq_idx, col] for seq_idx in range(num_seqs)]
        if "-" not in column_chars and "N" not in column_chars and "n" not in column_chars:
            start_col = col
            break

    end_col = alignment_length
    for col in range(alignment_length - 1, -1, -1):
        column_chars = [alignment[seq_idx, col] for seq_idx in range(num_seqs)]
        if "-" not in column_chars and "N" not in column_chars and "n" not in column_chars:
            end_col = col + 1
            break

    if start_col >= end_col:
        print("[WARNING]: Trimming boundaries invalid. Keeping the original alignment.")
        AlignIO.write(alignment, trimmed_file, "fasta")
        return

    trimmed_alignment = alignment[:, start_col:end_col]
    AlignIO.write(trimmed_alignment, trimmed_file, "fasta")
    print(
        f"[TRIMMING]: Alignment length reduced from {alignment_length} bp to "
        f"{trimmed_alignment.get_alignment_length()} bp."
    )


def run_fasttree_phylogeny(trimmed_fasta, tree_output_file):
    fasttree_bin = resolve_executable_path(config, "fasttree_bin_path")
    print("[FASTTree]: Building the maximum-likelihood phylogenetic tree using the GTR model...")
    cmd = [fasttree_bin, "-nt", "-gtr", trimmed_fasta]
    with open(tree_output_file, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
    print(f"[SUCCESS]: Tree written to {tree_output_file}")


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def write_analysis_metadata(output_dir, samples, analysis_name):
    metadata_path = os.path.join(output_dir, "analysis_metadata.json")
    payload = {
        "analysis_name": analysis_name,
        "input_samples": samples,
        "module": "module_7_tree_build",
    }
    with open(metadata_path, "w") as handle:
        json.dump(payload, handle, indent=4)
    print(f"[METADATA]: Wrote {metadata_path}")


def main():
    parser = argparse.ArgumentParser(description="Build a phylogeny from the module 6 polished consensus outputs.")
    parser.add_argument("--config", required=True, help="Path to the config Python file.")
    args = parser.parse_args()

    global config
    config = load_config(args.config)

    input_samples = get_input_samples(config)
    analysis_name = (
        get_config_value(config, "module_7_params", "analysis_name", "default")
        or "default"
    )
    threads = int(
        get_config_value(config, "module_7_params", "threads", config.cross_module_params.get("threads", 8))
    )

    output_dir = get_module_7_dir(config)
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nStarting module_7_tree_build...\nConfig: {args.config}\nAnalysis: {analysis_name}\nOutput: {output_dir}\n")

    records = []
    for sample_name in input_samples:
        try:
            fasta_path = find_polished_fasta(config, sample_name)
            sample_records = list(SeqIO.parse(fasta_path, "fasta"))
            if not sample_records:
                print(f"[WARN]: No records found in {fasta_path}; skipping.")
                continue

            for record in sample_records:
                record.id = sample_name
                record.description = ""
                records.append(record)
        except Exception as exc:
            print(f"[SKIP]: Could not include sample {sample_name}: {exc}")

    if not records:
        raise RuntimeError("No valid consensus FASTA records were available for tree building.")

    combined_fasta = os.path.join(output_dir, "polished_genomes_combined.fasta")
    SeqIO.write(records, combined_fasta, "fasta")
    print(f"[INFO]: Collected {len(records)} sequences into {combined_fasta}")

    aligned_fasta = os.path.join(output_dir, "polished_genomes_aligned.fasta")
    run_mafft_alignment(combined_fasta, aligned_fasta, threads)

    trimmed_fasta = os.path.join(output_dir, "polished_genomes_aligned_trimmed.fasta")
    trim_alignment_ends(aligned_fasta, trimmed_fasta)

    tree_output = os.path.join(output_dir, "polished_genomes_tree.nwk")
    run_fasttree_phylogeny(trimmed_fasta, tree_output)

    write_analysis_metadata(output_dir, input_samples, analysis_name)
    print("\n[module_7_tree_build] Finished all samples.")


if __name__ == "__main__":
    main()

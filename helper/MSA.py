#!/usr/bin/env python3
"""Simple MSA + FastTree rerun analogous to module 7, with hardcoded inputs.

Edit the constants below before running.
"""

import os
import shutil
import subprocess

from Bio import AlignIO, SeqIO

INPUT_FASTA = "/media/mngs/48TBRAID5HDD/20261008_230815_Renamed_Out_2/Consensus_plus_denovo_genomes.fasta"
OUTPUT_DIR = "/media/mngs/48TBRAID5HDD/20261008_230815_Renamed_Out_2/MSA_Output"
MAFFT_BIN = "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/mafft"
FASTTREE_BIN = "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/FastTree"
THREADS = 32
REFERENCE_FASTA = None  # e.g. "/PATH/TO/reference.fasta" or None
REFERENCE_ID = "reference"


def resolve_executable(path_or_name):
    if path_or_name and os.path.exists(path_or_name):
        return path_or_name
    resolved = shutil.which(path_or_name) if path_or_name else None
    if resolved:
        return resolved
    raise FileNotFoundError(f"Executable not found: {path_or_name}")


def run_mafft(input_fasta, output_fasta, mafft_bin, threads=32):
    cmd = [mafft_bin, "--auto", "--thread", str(threads), input_fasta]
    with open(output_fasta, "w") as handle:
        subprocess.run(cmd, stdout=handle, stderr=subprocess.DEVNULL, check=True)
    print(f"[MAFFT]: wrote alignment to {output_fasta}")


def trim_alignment_ends(alignment_file, trimmed_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    alignment_length = alignment.get_alignment_length()
    num_seqs = len(alignment)

    start_idx = 0
    for col in range(alignment_length):
        column_chars = [alignment[i, col] for i in range(num_seqs)]
        if "-" not in column_chars and "N" not in column_chars and "n" not in column_chars:
            start_idx = col
            break

    end_idx = alignment_length
    for col in range(alignment_length - 1, -1, -1):
        column_chars = [alignment[i, col] for i in range(num_seqs)]
        if "-" not in column_chars and "N" not in column_chars and "n" not in column_chars:
            end_idx = col + 1
            break

    if start_idx >= end_idx:
        print("[TRIM]: no valid trimming window found; keeping original alignment")
        AlignIO.write(alignment, trimmed_file, "fasta")
        return

    trimmed = alignment[:, start_idx:end_idx]
    AlignIO.write(trimmed, trimmed_file, "fasta")
    print(f"[TRIM]: reduced alignment length from {alignment_length} to {trimmed.get_alignment_length()}")


def run_fasttree(input_fasta, output_tree, fasttree_bin):
    cmd = [fasttree_bin, "-nt", "-gtr", input_fasta]
    with open(output_tree, "w") as handle:
        subprocess.run(cmd, stdout=handle, stderr=subprocess.DEVNULL, check=True)
    print(f"[FASTTREE]: wrote tree to {output_tree}")


def prepare_input_fasta(input_fasta, reference_fasta=None, reference_id="reference", output_dir="."):
    if reference_fasta is None:
        return input_fasta

    records = list(SeqIO.parse(input_fasta, "fasta"))
    ref_records = list(SeqIO.parse(reference_fasta, "fasta"))
    if not ref_records:
        raise ValueError(f"Reference FASTA is empty: {reference_fasta}")

    ref_record = ref_records[0].copy()
    ref_record.id = reference_id
    ref_record.description = ""
    records.append(ref_record)

    temp_input = os.path.join(output_dir, "with_reference.fasta")
    SeqIO.write(records, temp_input, "fasta")
    print(f"[REF]: included reference sequence {reference_id} from {reference_fasta}")
    return temp_input


def main():
    if not os.path.exists(INPUT_FASTA):
        raise FileNotFoundError(f"Input FASTA not found: {INPUT_FASTA}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    mafft_bin = resolve_executable(MAFFT_BIN)
    fasttree_bin = resolve_executable(FASTTREE_BIN)

    working_input = prepare_input_fasta(
        INPUT_FASTA,
        reference_fasta=REFERENCE_FASTA,
        reference_id=REFERENCE_ID,
        output_dir=OUTPUT_DIR,
    )

    base_name = os.path.splitext(os.path.basename(INPUT_FASTA))[0]
    aligned_fasta = os.path.join(OUTPUT_DIR, f"{base_name}_aligned.fasta")
    trimmed_fasta = os.path.join(OUTPUT_DIR, f"{base_name}_aligned_trimmed.fasta")
    tree_file = os.path.join(OUTPUT_DIR, f"{base_name}_tree.nwk")

    run_mafft(working_input, aligned_fasta, mafft_bin, threads=THREADS)
    trim_alignment_ends(aligned_fasta, trimmed_fasta)
    run_fasttree(trimmed_fasta, tree_file, fasttree_bin)

    print("\nDone.")
    print(f"alignment: {aligned_fasta}")
    print(f"trimmed: {trimmed_fasta}")
    print(f"tree: {tree_file}")


if __name__ == "__main__":
    main()

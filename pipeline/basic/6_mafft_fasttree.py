import os
import glob
import subprocess
from Bio import SeqIO
from Bio import AlignIO

# User Configuration
polished_output_dir = "/media/mngs/48TBRAID5HDD/viral_3/ivar_polished_results"
phylogeny_output_dir = "/media/mngs/48TBRAID5HDD/viral_3/phylogeny_results"
accession_list = [
    "SRR34884128",#8
    "SRR34884127",#7
    "SRR34884126",#7
    "SRR34884125",#9
    "SRR34884124",
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
    "SRR34884131",
    "SRR34884130",
    "SRR34884129",
    "SRR34884123",
    "SRR34884122",
    "SRR34884120",
    "SRR34884119",
    "SRR34884118",
    "SRR34884117",
    "SRR34884116",
    "SRR34884115",
    "SRR34884114",
    "SRR34884241",
    "SRR34884240",
    "SRR34884239",
    "SRR34884237",
    "SRR34884236",
    "SRR34884235",
    "SRR34884234",
    "SRR34884233",
    "SRR34884232",
]
threads = 32  # Adjust CPU threads for MAFFT and FastTree

# Binary paths for tools
mafft_bin = "mafft"
fasttree_bin = "FastTree"


def find_polished_fasta(acc):
    """Find the polished consensus FASTA produced by step 5 for one accession."""
    acc_dir = os.path.join(polished_output_dir, acc)
    candidates = [
        os.path.join(acc_dir, f"{acc}_polished_consensus.fa"),
        os.path.join(acc_dir, f"{acc}_polished_consensus.fasta"),
        os.path.join(acc_dir, f"{acc}_consensus.fa"),
    ]
    for path in candidates:
        if os.path.exists(path):
            print(f"[STEP 5 INPUT]: Selected polished FASTA for {acc}: {path}")
            return path

    matches = sorted(glob.glob(os.path.join(acc_dir, "*.fa"))) + sorted(
        glob.glob(os.path.join(acc_dir, "*.fasta"))
    )
    if matches:
        print(f"[STEP 5 INPUT]: Found polished FASTA candidates for {acc}: {matches}")
        return matches[0]

    print(f"[ERROR]: No polished FASTA found for {acc} under {acc_dir}")
    return None


def run_mafft_alignment(input_fasta, aligned_fasta):
    """Aligns genomes using Multiple Alignment using Fast Fourier Transform (MAFFT)."""
    print(f"[MAFFT]: Running multiple sequence alignment with {threads} threads...")
    cmd_mafft = [mafft_bin, "--auto", "--thread", str(threads), input_fasta]
    
    with open(aligned_fasta, "w") as out_f:
        subprocess.run(cmd_mafft, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
    print(f"[MAFFT]: Alignment written to {aligned_fasta}")


def trim_alignment_ends(alignment_file, trimmed_file):
    """Trims 5' and 3' regions of the multiple genome alignment to ensure no missing bases / edge bias."""
    print(f"[TRIMMING]: Removing terminal 5' and 3' regions containing gaps or missing bases...")
    alignment = AlignIO.read(alignment_file, "fasta")
    num_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()

    # Find first column from the left with zero gaps or missing bases across all sequences
    start_col = 0
    for col in range(alignment_length):
        column_chars = [alignment[seq_idx, col] for seq_idx in range(num_seqs)]
        if '-' not in column_chars and 'N' not in column_chars and 'n' not in column_chars:
            start_col = col
            break

    # Find last column from the right with zero gaps or missing bases across all sequences
    end_col = alignment_length
    for col in range(alignment_length - 1, -1, -1):
        column_chars = [alignment[seq_idx, col] for seq_idx in range(num_seqs)]
        if '-' not in column_chars and 'N' not in column_chars and 'n' not in column_chars:
            end_col = col + 1
            break

    if start_col >= end_col:
        print(f"[WARNING]: Trimming boundaries invalid. Keeping original alignment.")
        AlignIO.write(alignment, trimmed_file, "fasta")
        return

    trimmed_alignment = alignment[:, start_col:end_col]
    AlignIO.write(trimmed_alignment, trimmed_file, "fasta")
    print(f"[TRIMMING]: Alignment length reduced from {alignment_length} bp to {trimmed_alignment.get_alignment_length()} bp (retained columns {start_col} to {end_col}).")


def run_fasttree_phylogeny(trimmed_fasta, tree_output_file):
    """Generates a maximum-likelihood phylogenetic tree using FastTree with the GTR model."""
    print(f"[FASTTree]: Building maximum-likelihood phylogenetic tree using GTR model...")
    cmd_fasttree = [fasttree_bin, "-nt", "-gtr", trimmed_fasta]

    with open(tree_output_file, "w") as out_f:
        subprocess.run(cmd_fasttree, stdout=out_f, stderr=subprocess.DEVNULL, check=True)
    print(f"[SUCCESS]: Phylogenetic tree successfully written to {tree_output_file}")


def main():
    os.makedirs(phylogeny_output_dir, exist_ok=True)

    print(f"\n==========================================")
    print(f" Step 5 -> Step 6: Building phylogeny from polished consensus genomes")
    print(f"==========================================")

    polished_records = []
    for acc in accession_list:
        polished_fasta = find_polished_fasta(acc)
        if not polished_fasta:
            print(f"[SKIP]: No polished FASTA found for {acc}; continuing.")
            continue

        records = list(SeqIO.parse(polished_fasta, "fasta"))
        if not records:
            print(f"[WARN]: No records were parsed from {polished_fasta}")
            continue

        for rec in records:
            rec.id = acc
            rec.description = ""
            polished_records.append(rec)

    if not polished_records:
        print(f"[FATAL]: No polished consensus genomes were available for phylogeny. Exiting pipeline.")
        return

    combined_passing_fasta = os.path.join(
        phylogeny_output_dir, "polished_genomes_combined.fasta"
    )
    SeqIO.write(polished_records, combined_passing_fasta, "fasta")
    print(
        f"[INFO]: Collected {len(polished_records)} polished sequences from step 5 into {combined_passing_fasta}"
    )

    print(f"\n==========================================")
    print(f" Step 2: Multiple Sequence Alignment (MAFFT)")
    print(f"==========================================")
    raw_aligned_fasta = os.path.join(phylogeny_output_dir, "genomes_aligned_raw.fasta")
    run_mafft_alignment(combined_passing_fasta, raw_aligned_fasta)

    print(f"\n==========================================")
    print(f" Step 3: Alignment End-Trimming")
    print(f"==========================================")
    trimmed_aligned_fasta = os.path.join(phylogeny_output_dir, "genomes_aligned_trimmed.fasta")
    trim_alignment_ends(raw_aligned_fasta, trimmed_aligned_fasta)

    print(f"\n==========================================")
    print(f" Step 4: Maximum-Likelihood Tree (FastTree)")
    print(f"==========================================")
    tree_file = os.path.join(phylogeny_output_dir, "viral_phylogeny_tree.nwk")
    run_fasttree_phylogeny(trimmed_aligned_fasta, tree_file)

    print(f"\n[COMPLETE]: Phylogeny pipeline finished successfully. Results saved in {phylogeny_output_dir}")


if __name__ == "__main__":
    main()
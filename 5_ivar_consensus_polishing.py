import csv
import glob
import os
import re
import subprocess
from pathlib import Path

# User Configuration
scaffold_dir = "/media/mngs/48TBRAID5HDD/viral_3/scaffolding_and_checkv_results_nucmer_blast"  # Step 4 output directory
validation_reads_dir = "/media/mngs/48TBRAID5HDD/viral_3/kraken2_out"  # Directory containing fastp trimmed reads
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
output_dir = "/media/mngs/48TBRAID5HDD/viral_3/ivar_polished_results"
threads = 64

# Executable Paths / Commands
minimap2_bin = "minimap2"
samtools_bin = "samtools"
ivar_bin = "ivar"
checkv_bin = "checkv"
checkv_db_dir = "/media/mngs/48TBRAID5HDD/viral/checkv_db/checkv-db-v1.5"

# iVar Threshold Parameters
min_qual = 20  # Minimum base quality score to consider (-Q)
min_freq = 0.6  # Minimum allele frequency threshold to call consensus (-t)
min_depth = 10  # Minimum coverage depth to call a base instead of N (-m)


def run_cmd(cmd):
    """Run a shell command and print it."""
    print(f"\n[EXEC]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def clean_scaffold_fasta(scaffold_path, acc_out_dir):
    """Copy scaffold FASTA to output dir, removing malformed non-FASTA lines.

    Accommodates all standard IUPAC nucleotide symbols and drops zero-length sequences.
    """
    cleaned_path = os.path.join(acc_out_dir, f"{Path(scaffold_path).stem}_clean.fasta")
    header_re = re.compile(r"^>.+")
    # Allow standard bases, soft-masked bases, IUPAC ambiguity codes, and gaps
    seq_re = re.compile(r"^[A-Za-z\-]+$")

    kept_headers = 0
    dropped_lines = 0
    dropped_empty_seqs = 0
    current_header = None
    current_seq_lines = []

    def flush_record(header, seq_lines, out_handle):
        nonlocal kept_headers, dropped_empty_seqs
        if not header:
            return
        seq = "".join(seq_lines).strip()
        if not seq:
            dropped_empty_seqs += 1
            return
        out_handle.write(header + "\n")
        # Wrap sequence at 80 chars for clean FASTA formatting
        for i in range(0, len(seq), 80):
            out_handle.write(seq[i:i + 80] + "\n")
        kept_headers += 1

    with open(scaffold_path, "r") as in_f, open(cleaned_path, "w") as out_f:
        for raw_line in in_f:
            line = raw_line.strip()
            if not line:
                continue
            if header_re.match(line):
                flush_record(current_header, current_seq_lines, out_f)
                current_header = line
                current_seq_lines = []
            elif seq_re.match(line):
                current_seq_lines.append(line)
            else:
                dropped_lines += 1
        flush_record(current_header, current_seq_lines, out_f)

    print(
        f"[SCAFFOLD CLEAN]: Kept {kept_headers} headers, dropped {dropped_lines} malformed lines "
        f"and {dropped_empty_seqs} empty sequences."
    )
    return cleaned_path


def find_scaffold_fasta(acc):
    """Locate the scaffold FASTA written by Step 4."""
    acc_dir = os.path.join(scaffold_dir, acc)
    print(f"[SCAFFOLD SEARCH]: Looking for Step 4 scaffold for {acc} under {acc_dir}")

    candidates = [
        os.path.join(acc_dir, f"{acc}_nucmer_scaffold.fasta"),
        os.path.join(acc_dir, f"{acc}_nucmer_1to1_scaffold.fasta"),
        os.path.join(acc_dir, f"{acc}_scaffolded.fasta"),
    ]
    print(f"[SCAFFOLD SEARCH]: Candidate paths: {candidates}")

    for scaffold_path in candidates:
        if os.path.exists(scaffold_path):
            print(f"[SCAFFOLD SEARCH]: Selected scaffold: {scaffold_path}")
            return scaffold_path

    alt_matches = sorted(glob.glob(os.path.join(acc_dir, "*.fasta")))
    if alt_matches:
        print(f"[INFO]: Found scaffold-like FASTA(s) for {acc}: {alt_matches}")
        return alt_matches[0]

    print(
        f"[ERROR]: Scaffold FASTA not found for {acc} under {acc_dir}. "
        f"Checked: {candidates}"
    )
    return None


def find_trimmed_reads(acc):
    """Locate trimmed paired or unpaired reads for alignment."""
    print(f"[READ SEARCH]: Looking for validation reads for {acc} under {validation_reads_dir}/{acc}")
    r1_pattern = os.path.join(
        validation_reads_dir, acc, f"{acc}_dominant_species_1_trimmed.fastq*"
    )
    r2_pattern = os.path.join(
        validation_reads_dir, acc, f"{acc}_dominant_species_2_trimmed.fastq*"
    )
    u_pattern = os.path.join(
        validation_reads_dir,
        acc,
        f"{acc}_dominant_species_unpaired_trimmed.fastq*",
    )

    r1_matches = glob.glob(r1_pattern)
    r2_matches = glob.glob(r2_pattern)
    u_matches = glob.glob(u_pattern)
    print(
        f"[READ SEARCH]: paired R1={r1_matches}, paired R2={r2_matches}, unpaired={u_matches}"
    )

    if r1_matches and r2_matches:
        print(f"[READ SEARCH]: Selected paired reads for {acc}: {r1_matches[0]} and {r2_matches[0]}")
        return ("paired", r1_matches[0], r2_matches[0])
    elif u_matches:
        print(f"[READ SEARCH]: Selected unpaired reads for {acc}: {u_matches[0]}")
        return ("unpaired", u_matches[0], None)
    else:
        print(
            f"[ERROR]: Could not find trimmed FASTQ reads for {acc} in {validation_reads_dir}/{acc}"
        )
        return (None, None, None)


def run_ivar_polishing_pipeline(
    acc, scaffold_fasta, read_type, r1, r2_or_u, acc_out_dir
):
    """Aligns reads to scaffold, sorts BAM, and runs iVar consensus calling."""
    sorted_bam = os.path.join(acc_out_dir, f"{acc}_aligned_sorted.bam")
    consensus_prefix = os.path.join(acc_out_dir, f"{acc}_polished_consensus")
    consensus_fasta = f"{consensus_prefix}.fa"

    # Use the output directory as temp space for samtools sort
    sort_tmp_dir = os.path.join(acc_out_dir, "samtools_sort_tmp")
    os.makedirs(sort_tmp_dir, exist_ok=True)

    # Step A: Map trimmed reads to scaffold using Minimap2 and sort with Samtools
    print(f"[ALIGNMENT]: Mapping {read_type} reads to scaffold for {acc}...")
    print(f"[ALIGNMENT]: Scaffold file: {scaffold_fasta}")
    print(f"[ALIGNMENT]: Read 1: {r1}")
    if read_type == "paired":
        print(f"[ALIGNMENT]: Read 2: {r2_or_u}")

    minimap_cmd = [minimap2_bin, "-ax", "sr", scaffold_fasta]
    if read_type == "paired":
        minimap_cmd.extend([r1, r2_or_u])
    else:
        minimap_cmd.append(r1)

    samtools_view_cmd = [samtools_bin, "view", "-bS", "-"]
    samtools_sort_cmd = [
        samtools_bin,
        "sort",
        "-@",
        str(threads),
        "-T",
        sort_tmp_dir,
        "-o",
        sorted_bam,
        "-",
    ]

    print(f"[ALIGNMENT EXEC]: {' '.join(minimap_cmd)} | {' '.join(samtools_view_cmd)} | {' '.join(samtools_sort_cmd)}")
    
    # DEVNULL on intermediate stderr streams prevents subprocess deadlock on heavy log output
    p1 = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(
        samtools_view_cmd,
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    p3 = subprocess.Popen(
        samtools_sort_cmd,
        stdin=p2.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    p1.stdout.close()
    p2.stdout.close()

    _, p3_err = p3.communicate()
    p2.wait()
    p1.wait()

    if p1.returncode != 0 or p2.returncode != 0 or p3.returncode != 0:
        print(f"[ERROR]: Alignment pipeline failed for {acc}. Sort error: {p3_err.decode(errors='replace')}")
        return None

    if not os.path.exists(sorted_bam) or os.path.getsize(sorted_bam) == 0:
        print(f"[ERROR]: Sorted BAM is missing or empty for {acc}.")
        return None

    print(f"[SUCCESS]: Created sorted BAM file: {sorted_bam}")

    # Index BAM
    print(f"[ALIGNMENT]: Indexing BAM: {sorted_bam}")
    try:
        subprocess.run(
            [samtools_bin, "index", sorted_bam],
            check=True,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        print(f"[ERROR]: samtools index failed for {acc}: {e.stderr.decode(errors='replace')}")
        return None
    print(f"[ALIGNMENT]: BAM index created successfully for {sorted_bam}")

    # Step B: Generate mpileup (depth capped at 1000x) and pipe into iVar consensus
    print(
        f"[IVAR CONSENSUS]: Calling consensus (min_freq={min_freq}, min_qual={min_qual}, min_depth={min_depth})..."
    )

    mpileup_cmd = [
        samtools_bin,
        "mpileup",
        "-aa",
        "-A",
        "-d",
        "1000",
        "-Q",
        str(min_qual),
        sorted_bam,
    ]
    ivar_cmd = [
        ivar_bin,
        "consensus",
        "-p",
        consensus_prefix,
        "-t",
        str(min_freq),
        "-m",
        str(min_depth),
        "-n",
        "N",
    ]

    print(f"[IVAR EXEC]: {' '.join(mpileup_cmd)} | {' '.join(ivar_cmd)}")

    p_mpileup = subprocess.Popen(
        mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    p_ivar = subprocess.Popen(
        ivar_cmd, stdin=p_mpileup.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    p_mpileup.stdout.close()
    _, p_ivar_err = p_ivar.communicate()
    p_mpileup.wait()

    if p_mpileup.returncode != 0 or p_ivar.returncode != 0:
        print(f"[ERROR]: iVar consensus failed for {acc}: {p_ivar_err.decode(errors='replace')}")
        return None

    if not os.path.exists(consensus_fasta):
        print(f"[ERROR]: iVar consensus failed to output {consensus_fasta}.")
        return None

    print(f"[SUCCESS]: Created polished consensus sequence: {consensus_fasta}")
    return consensus_fasta


def count_n_bases(fasta_path):
    """Count the number of N bases in a FASTA sequence using the same logic as Step 4."""
    if not fasta_path or not os.path.exists(fasta_path):
        return 0

    total_n = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            total_n += line.upper().count("N")
    return total_n


def run_checkv_qc(acc, polished_fasta, acc_out_dir):
    """Run CheckV end-to-end on the final polished consensus sequence."""
    print(f"[CHECKV]: Assessing final polished genome completeness for {acc}...")
    checkv_out_dir = os.path.join(acc_out_dir, "checkv_out")
    os.makedirs(checkv_out_dir, exist_ok=True)
    print(f"[CHECKV]: Output directory: {checkv_out_dir}")

    cmd_checkv = [
        checkv_bin,
        "end_to_end",
        polished_fasta,
        checkv_out_dir,
        "-t",
        str(threads),
        "--restart",
    ]
    if checkv_db_dir and os.path.exists(checkv_db_dir):
        cmd_checkv.extend(["-d", checkv_db_dir])

    print(f"[CHECKV EXEC]: {' '.join(cmd_checkv)}")
    subprocess.run(
        cmd_checkv,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    summary_tsv = os.path.join(checkv_out_dir, "quality_summary.tsv")
    if not os.path.exists(summary_tsv):
        print(f"[ERROR]: CheckV did not produce quality_summary.tsv for {acc}.")
        return "N/A", "N/A", "N/A", ""

    with open(summary_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            completeness = (
                row.get("completeness")
                or row.get("checkv_completeness")
                or "N/A"
            )
            quality = row.get("checkv_quality") or row.get("quality") or "N/A"
            kmer_freq = row.get("kmer_freq") or "N/A"
            warnings_str = row.get("warnings") or ""
            print(
                f"[CHECKV RESULTS] {acc}: Completeness={completeness}% | "
                f"Quality={quality} | kmer_freq={kmer_freq}"
            )
            return completeness, quality, kmer_freq, warnings_str

    return "N/A", "N/A", "N/A", ""


def write_qc_report(
    acc,
    scaffold_fasta,
    polished_fasta,
    completeness,
    quality,
    kmer_freq,
    warnings_str,
    acc_out_dir,
):
    """Write a summary report comparing pre- and post-polishing outputs."""
    scaffold_n_count = count_n_bases(scaffold_fasta)
    polished_n_count = count_n_bases(polished_fasta)

    polished_length = 0
    if polished_fasta and os.path.exists(polished_fasta):
        try:
            with open(polished_fasta, "r") as f:
                for line in f:
                    if not line.startswith(">"):
                        polished_length += len(line.strip())
        except Exception as e:
            print(f"[WARNING]: Could not determine polished genome length for {acc}: {e}")

    qc_log = os.path.join(acc_out_dir, f"{acc}_ivar_polished_report.txt")
    with open(qc_log, "w") as out:
        out.write(f"Accession: {acc}\n")
        out.write(f"Input Scaffold FASTA: {scaffold_fasta}\n")
        out.write(f"Input Scaffold N Bases: {scaffold_n_count} bp\n")
        out.write(f"Polished FASTA: {polished_fasta}\n")
        out.write(f"Polished Genome Length: {polished_length} bp\n")
        out.write(f"Polished Genome N Bases: {polished_n_count} bp\n")
        out.write(f"CheckV Completeness: {completeness}%\n")
        out.write(f"CheckV Quality: {quality}\n")
        out.write(f"Kmer Frequency: {kmer_freq}\n")
        out.write(f"CheckV Warnings: {warnings_str}\n")
    print(f"[COMPLETE]: Final report written to {qc_log}")


def process_accession(acc):
    print("\n==========================================")
    print(f" Processing iVar Polishing: {acc}")
    print("==========================================")

    acc_out_dir = os.path.join(output_dir, acc)
    os.makedirs(acc_out_dir, exist_ok=True)
    print(f"[PROCESS]: Beginning polishing workflow for {acc} in {acc_out_dir}")

    scaffold_fasta = find_scaffold_fasta(acc)
    if not scaffold_fasta:
        print(f"[ERROR]: Missing scaffold FASTA for {acc}. Skipping polishing.")
        return

    read_type, r1, r2_or_u = find_trimmed_reads(acc)
    if not read_type:
        print(f"[ERROR]: Missing trimmed FASTQs for {acc}. Skipping polishing.")
        return

    print(f"[PROCESS]: Using scaffold {scaffold_fasta} with read type {read_type}")

    # Clean malformed scaffold FASTA before alignment
    cleaned_scaffold = clean_scaffold_fasta(scaffold_fasta, acc_out_dir)
    print(f"[PROCESS]: Cleaned scaffold saved to {cleaned_scaffold}")

    polished_fasta = run_ivar_polishing_pipeline(
        acc, cleaned_scaffold, read_type, r1, r2_or_u, acc_out_dir
    )
    if not polished_fasta:
        print(f"[ERROR]: Polishing pipeline failed for {acc}. Skipping final CheckV.")
        write_qc_report(
            acc,
            scaffold_fasta,
            None,
            "N/A",
            "N/A",
            "N/A",
            "iVar polishing failed",
            acc_out_dir,
        )
        return

    completeness, quality, kmer_freq, warnings_str = run_checkv_qc(
        acc, polished_fasta, acc_out_dir
    )
    write_qc_report(
        acc,
        scaffold_fasta,
        polished_fasta,
        completeness,
        quality,
        kmer_freq,
        warnings_str,
        acc_out_dir,
    )


def main():
    os.makedirs(output_dir, exist_ok=True)
    for acc in accession_list:
        try:
            process_accession(acc)
        except Exception as e:
            print(f"[ERROR]: iVar polishing pipeline failed for {acc}: {e}")


if __name__ == "__main__":
    main()
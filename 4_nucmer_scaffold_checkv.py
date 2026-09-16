import csv
import glob
import os
import subprocess
from Bio import SeqIO

# User Configuration
spades_dir = "/media/mngs/48TBRAID5HDD/viral_3/spades_out"
ref_selection_dir = "/media/mngs/48TBRAID5HDD/viral_3/automated_reference_selection_blast"
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
output_dir = (
    "/media/mngs/48TBRAID5HDD/viral_3/scaffolding_and_checkv_results_nucmer_blast"
)
threads = 64

nucmer_bin = "nucmer"
delta_filter_bin = "delta-filter"
show_tiling_bin = "show-tiling"
show_coords_bin = "show-coords"
checkv_bin = "checkv"
checkv_db_dir = "/media/mngs/48TBRAID5HDD/viral/checkv_db/checkv-db-v1.5"


def get_fasta_length(fasta_path):
    """Calculates total length of sequences in a FASTA file."""
    total_len = 0
    if not os.path.exists(fasta_path):
        return 0
    with open(fasta_path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                total_len += len(line.strip())
    return total_len


def find_selected_reference(acc):
    """Locate the selected reference genome symlink."""
    acc_dir = os.path.join(ref_selection_dir, acc)
    print(f"[REF LOOKUP]: Searching for selected reference for {acc} in {acc_dir}")
    if not os.path.exists(acc_dir):
        print(f"[REF LOOKUP]: Directory does not exist: {acc_dir}")
        return None

    matches = sorted(
        glob.glob(os.path.join(acc_dir, "selected_reference_*.fasta"))
    )
    if matches:
        print(f"[REF LOOKUP]: Found {len(matches)} selected reference FASTA(s):")
        for m in matches:
            print(f"  - {m}")
        return matches[0]
    print(f"[REF LOOKUP]: No selected reference FASTA found for {acc}")
    return None


def run_nucmer_scaffolding(
    ref_fasta, query_fasta, acc_out_dir, acc, min_id=85.0
):
    """Scaffolds viral contigs while clipping unaligned query overhangs (SPAdes dimers)
    via show-coords strictly at alignment boundaries (zero padding). Preserves internal 
    insertions/deletions and resolves contig overlaps by sequence identity competition.
    """
    prefix = os.path.join(acc_out_dir, f"{acc}_nucmer")
    delta_file = f"{prefix}.delta"
    filtered_delta = f"{prefix}.1to1.delta"
    coords_file = f"{prefix}.coords"
    scaffold_fasta = os.path.join(acc_out_dir, f"{acc}_nucmer_scaffold.fasta")

    # Determine reference length for 3' terminal padding checks
    ref_len = sum(len(r.seq) for r in SeqIO.parse(ref_fasta, "fasta"))

    print(f"[NUCMER]: Starting MUMmer scaffolding for {acc}")

    # 1. Run Nucmer alignment
    cmd_nucmer = [
        nucmer_bin, "--maxmatch", "-l", "15", "-c", "30",
        "-p", prefix, ref_fasta, query_fasta
    ]
    subprocess.run(cmd_nucmer, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Inner helper to filter and parse coordinates for a given identity threshold
    def parse_coords_for_identity(cutoff_id):
        cmd_filter = [delta_filter_bin, "-q", "-l", "30", "-i", str(cutoff_id), delta_file]
        with open(filtered_delta, "w") as out:
            subprocess.run(cmd_filter, check=True, stdout=out, stderr=subprocess.DEVNULL)

        cmd_coords = [show_coords_bin, "-T", "-r", "-l", filtered_delta]
        with open(coords_file, "w") as out:
            subprocess.run(cmd_coords, check=True, stdout=out, stderr=subprocess.DEVNULL)

        alignments = []
        if not os.path.exists(coords_file) or os.path.getsize(coords_file) == 0:
            return alignments

        with open(coords_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 11:
                    continue
                
                # Skip header lines
                if not parts[0].isdigit() or not parts[1].isdigit():
                    continue

                r_start, r_end = int(parts[0]), int(parts[1])
                q_start, q_end = int(parts[2]), int(parts[3])
                
                strand = "+" if q_start <= q_end else "-"
                q_min = min(q_start, q_end)
                q_max = max(q_start, q_end)
                identity = float(parts[6])  # Column 7 in show-coords -T output (% IDY)
                contig_id = parts[-1]      # Last column in show-coords -T output

                alignments.append({
                    "r_start": r_start,
                    "r_end": r_end,
                    "q_min": q_min,
                    "q_max": q_max,
                    "strand": strand,
                    "identity": identity,
                    "contig_id": contig_id
                })
        return alignments

    # 2. Try primary identity threshold (85.0%)
    alignments = parse_coords_for_identity(min_id)

    # 3. Fallback: retry at 70.0% if 0 alignments pass
    if not alignments:
        print(f"[NUCMER WARN - {acc}]: No alignments passed {min_id}% identity. Retrying with fallback cutoff (70.0%)...")
        alignments = parse_coords_for_identity(70.0)

    if not alignments:
        print(f"[NUCMER ERROR - {acc}]: No valid alignment records parsed from {coords_file}")
        return None

    # Sort alignments along reference genome (5' to 3')
    alignments.sort(key=lambda x: x["r_start"])

    # 3.5 PRE-FILTERING: Remove fully nested alignments
    uncontained_alignments = []
    for i, aln in enumerate(alignments):
        is_nested = False
        for j, other in enumerate(alignments):
            if i == j:
                continue
            # Check if 'aln' is completely contained within 'other' on reference coordinates
            if other["r_start"] <= aln["r_start"] and other["r_end"] >= aln["r_end"]:
                # Tie-breaker for identical spans: keep the first encountered record
                if other["r_start"] == aln["r_start"] and other["r_end"] == aln["r_end"]:
                    if j < i:
                        is_nested = True
                        break
                else:
                    is_nested = True
                    print(
                        f"[NUCMER FILTER - {acc}]: Discarding nested contig {aln['contig_id']} "
                        f"({aln['r_start']}-{aln['r_end']}) fully contained within {other['contig_id']} "
                        f"({other['r_start']}-{other['r_end']})"
                    )
                    break
        if not is_nested:
            uncontained_alignments.append(aln)

    alignments = uncontained_alignments

    if not alignments:
        print(f"[NUCMER ERROR - {acc}]: No alignments remained after nested filtering.")
        return None

    # 4. Resolve overlapping contig alignments based on highest sequence identity
    resolved = []
    for aln in alignments:
        while resolved and resolved[-1]["r_end"] >= aln["r_start"]:
            prev = resolved[-1]
            overlap = prev["r_end"] - aln["r_start"] + 1

            if prev["identity"] >= aln["identity"]:
                # Previous contig has higher/equal identity -> Trim 5' ref side of current contig
                aln["r_start"] += overlap
                if aln["strand"] == "+":
                    aln["q_min"] += overlap
                else:  # "-" strand: 5' ref corresponds to q_max
                    aln["q_max"] -= overlap

                # Discard current alignment if fully consumed
                if aln["q_min"] > aln["q_max"] or aln["r_start"] > aln["r_end"]:
                    aln = None
                    break
            else:
                # Current contig has higher identity -> Trim 3' ref side of previous contig
                prev["r_end"] -= overlap
                if prev["strand"] == "+":
                    prev["q_max"] -= overlap
                else:  # "-" strand: 3' ref corresponds to q_min
                    prev["q_min"] += overlap

                # Discard previous alignment if fully consumed and re-check with preceding contig
                if prev["q_min"] > prev["q_max"] or prev["r_start"] > prev["r_end"]:
                    resolved.pop()

        if aln is not None:
            resolved.append(aln)

    if not resolved:
        print(f"[NUCMER ERROR - {acc}]: No alignments remained after overlap resolution.")
        return None

    # 5. Build scaffold sequence using resolved contig coordinates (Zero Padding)
    contigs_dict = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))
    scaffold_chunks = []
    prev_r_end = 0

    for idx, aln in enumerate(resolved):
        cid = aln["contig_id"]
        if cid not in contigs_dict:
            print(f"[NUCMER WARN - {acc}]: Contig {cid} not found in query FASTA.")
            continue

        raw_record = contigs_dict[cid]
        seq_len = len(raw_record.seq)

        # Strict alignment-boundary extraction (0-based conversion, no unaligned padding)
        q_start_idx = max(0, aln["q_min"] - 1)
        q_end_idx = min(seq_len, aln["q_max"])

        # Slice query segment using contig coordinates
        clipped_seq_obj = raw_record.seq[q_start_idx:q_end_idx]

        # Apply reverse complement post-slicing for minus strand contigs
        if aln["strand"] == "-":
            clipped_seq_obj = clipped_seq_obj.reverse_complement()

        clipped_seq = str(clipped_seq_obj)

        # Reference gap N-filling between non-overlapping chunks
        if idx == 0:
            if aln["r_start"] > 1:
                scaffold_chunks.append("N" * (aln["r_start"] - 1))
        else:
            gap_to_prev = aln["r_start"] - prev_r_end - 1
            if gap_to_prev > 0:
                scaffold_chunks.append("N" * gap_to_prev)

        scaffold_chunks.append(clipped_seq)
        prev_r_end = max(prev_r_end, aln["r_end"])

        print(f"[NUCMER - {acc}]: Retained {cid} [{aln['strand']}] (ID: {aln['identity']}%) aligned region ({len(clipped_seq)} bp extracted)")

    # 3' terminal padding
    if prev_r_end < ref_len:
        missing_3prime = ref_len - prev_r_end
        scaffold_chunks.append("N" * missing_3prime)
        print(f"[NUCMER - {acc}]: Appended {missing_3prime} Ns to fill missing 3' terminal region")

    final_sequence = "".join(scaffold_chunks)

    if not final_sequence:
        print(f"[NUCMER ERROR - {acc}]: Scaffold sequence assembly produced 0 bp.")
        return None

    with open(scaffold_fasta, "w") as out_f:
        out_f.write(f">{acc}_nucmer_1to1_scaffold\n{final_sequence}\n")

    print(f"[NUCMER SUCCESS - {acc}]: Scaffold written to {scaffold_fasta} ({len(final_sequence)} bp)")
    return scaffold_fasta


def run_checkv_qc(target_fasta, checkv_out_dir):
    """Run CheckV end-to-end on a FASTA file and return metrics."""
    os.makedirs(checkv_out_dir, exist_ok=True)
    print(f"[CHECKV]: Running CheckV for scaffold: {target_fasta}")
    print(f"[CHECKV]: Output directory: {checkv_out_dir}")

    cmd_checkv = [
        checkv_bin,
        "end_to_end",
        target_fasta,
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
    print(f"[CHECKV]: Looking for summary file at {summary_tsv}")
    if not os.path.exists(summary_tsv):
        print(f"[CHECKV WARNING]: {summary_tsv} was not found.")
        return "N/A", "N/A", "N/A", "quality_summary.tsv missing"

    with open(summary_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        print(f"[CHECKV]: Found {len(rows)} row(s) in quality_summary.tsv")
        for row in rows:
            completeness = (
                row.get("completeness")
                or row.get("checkv_completeness")
                or "N/A"
            )
            quality = row.get("checkv_quality") or row.get("quality") or "N/A"
            kmer_freq = row.get("kmer_freq") or "N/A"
            warnings_str = row.get("warnings") or ""
            print(f"[CHECKV RESULTS]: completeness={completeness}% | quality={quality} | kmer_freq={kmer_freq} | warnings={warnings_str}")
            return completeness, quality, kmer_freq, warnings_str

    print(f"[CHECKV WARNING]: quality_summary.tsv was empty or had no parseable rows.")
    return "N/A", "N/A", "N/A", "No matching record found"


def write_qc_report(
    acc, ref_basename, scaffold_len, ref_len, n_count, checkv_stats, acc_out_dir
):
    """Writes a clean summary report of the scaffold and CheckV metrics."""
    qc_log = os.path.join(acc_out_dir, f"{acc}_nucmer_checkv_report.txt")
    comp, qual, kmer, warn = checkv_stats

    with open(qc_log, "w") as out:
        out.write("==================================================\n")
        out.write(f" MUMmer 1-to-1 Scaffolding QC Report: {acc}\n")
        out.write("==================================================\n")
        out.write(f"Reference Genome:        {ref_basename}\n")
        out.write(f"Reference Genome Length: {ref_len} bp\n")
        out.write(f"Scaffold Length:         {scaffold_len} bp\n")
        out.write(f"Scaffold N Bases:        {n_count} bp\n\n")

        out.write("--- CheckV Metrics ---\n")
        out.write(f"Completeness:  {comp}%\n")
        out.write(f"Quality:       {qual}\n")
        out.write(f"K-mer Freq:    {kmer}\n")
        out.write(f"Warnings:      {warn}\n")

    print(f"[COMPLETE]: QC report written to {qc_log}")


def process_accession(acc):
    print(f"\n==========================================")
    print(f" Processing MUMmer Scaffolding & CheckV: {acc}")
    print(f"==========================================")

    acc_out_dir = os.path.join(output_dir, acc)
    os.makedirs(acc_out_dir, exist_ok=True)
    print(f"[OUTPUT]: Creating per-accession directory: {acc_out_dir}")

    transcripts_fasta = os.path.join(spades_dir, acc, "transcripts.fasta")
    print(f"[INPUT]: Looking for SPAdes transcripts FASTA at {transcripts_fasta}")
    if not os.path.exists(transcripts_fasta):
        print(f"[ERROR]: {transcripts_fasta} not found. Skipping {acc}.")
        return
    print(f"[INPUT]: Found transcripts FASTA ({os.path.getsize(transcripts_fasta)} bytes)")

    ref_fasta = find_selected_reference(acc)
    if not ref_fasta or not os.path.exists(ref_fasta):
        print(f"[ERROR]: Reference FASTA not found for {acc}. Skipping.")
        return

    ref_basename = os.path.basename(ref_fasta)
    ref_len = get_fasta_length(ref_fasta)
    print(f"[INPUT]: Using reference genome {ref_basename} ({ref_len} bp)")

    # 1. Run Nucmer 1-to-1 Tiling Scaffolding
    print(f"[STEP 1]: Running MUMmer-based scaffold construction")
    scaffold_fasta = run_nucmer_scaffolding(
        ref_fasta, transcripts_fasta, acc_out_dir, acc
    )
    if not scaffold_fasta:
        print(f"[ERROR]: MUMmer scaffolding failed for {acc}.")
        return

    # Calculate scaffold length and N base count
    scaffold_record = SeqIO.read(scaffold_fasta, "fasta")
    scaffold_len = len(scaffold_record.seq)
    n_count = scaffold_record.seq.upper().count("N")
    print(f"[STEP 1]: Scaffold generated: {scaffold_fasta} ({scaffold_len} bp, {n_count} Ns)")

    # 2. CheckV evaluation on generated scaffold
    checkv_dir = os.path.join(acc_out_dir, "checkv_nucmer")
    print(f"[STEP 2]: Running CheckV on scaffold in {checkv_dir}")
    checkv_stats = run_checkv_qc(scaffold_fasta, checkv_dir)

    # 3. Write summary report
    write_qc_report(
        acc, ref_basename, scaffold_len, ref_len, n_count, checkv_stats, acc_out_dir
    )


def main():
    os.makedirs(output_dir, exist_ok=True)
    for acc in accession_list:
        try:
            process_accession(acc)
        except Exception as e:
            print(f"[ERROR]: Pipeline failed for {acc}: {e}")


if __name__ == "__main__":
    main()
import argparse
import csv
import glob
import importlib.util
import json
import os
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO


# ---------------------------------------------------------------------------
# iVar Threshold Parameters (Hardcoded for Two-Tiered Rescue Strategy)
# ---------------------------------------------------------------------------
# Tier 1: Strict - High confidence for core genome SNP calling
STRICT_PARAMS = {"min_qual": 20, "min_freq": 0.60, "min_depth": 10}

# Tier 2: Relaxed - Rescue low-coverage terminal UTRs and gaps
# Note: freq=0.51 prevents arbitrary base calls on exact 50/50 mixed-infection ties
RELAXED_PARAMS = {"min_qual": 15, "min_freq": 0.51, "min_depth": 3}


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
    value = get_config_value(config, "module_6_params", key, None)
    if value is None:
        raise FileNotFoundError(
            f"Executable '{key}' is not configured. Set an explicit path or command in module_6_params."
        )
    if os.path.exists(value):
        return value
    resolved = shutil.which(value)
    if resolved:
        return resolved
    raise FileNotFoundError(
        f"Configured executable '{value}' for '{key}' was not found. Set an explicit path in module_6_params."
    )


# ---------------------------------------------------------------------------
# Output path resolution based on the real module contracts
# ---------------------------------------------------------------------------

def get_module_4_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_4", sample_name)

def get_module_5_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_5", sample_name)

def get_module_6_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_6", sample_name)

def get_module_4_summary_path(config, sample_name):
    return os.path.join(get_module_4_dir(config, sample_name), f"{sample_name}_reference_summary.json")


# ---------------------------------------------------------------------------
# Input discovery: read pairs and scaffold inputs
# ---------------------------------------------------------------------------

def locate_fastq_pair(config, sample_name):
    output_dir = config.cross_module_params["output_dir"]
    fastp_dir = os.path.join(output_dir, "module_1", "fastp", sample_name)
    raw_dir = os.path.join(config.cross_module_params["read_dir"], sample_name)
    search_dirs = [p for p in [fastp_dir, raw_dir] if os.path.isdir(p)]

    if not search_dirs: return None, None

    read_1_str = config.cross_module_params.get("read_1_str", "_1.fastq")
    read_2_str = config.cross_module_params.get("read_2_str", "_2.fastq")
    singleton_str = config.cross_module_params.get("singleton_str")

    r1_candidates, r2_candidates = [], []
    for base_dir in search_dirs:
        r1_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}_R1.fastq*")))
        r1_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}{read_1_str}*")))
        r1_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}_1.fastq*")))
        r1_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}*R1*.fastq*")))

        r2_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}_R2.fastq*")))
        r2_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}{read_2_str}*")))
        r2_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}_2.fastq*")))
        r2_candidates.extend(glob.glob(os.path.join(base_dir, f"{sample_name}*R2*.fastq*")))

    r1_matches = sorted(set(r1_candidates))
    r2_matches = sorted(set(r2_candidates))

    if r1_matches and r2_matches: return r1_matches[0], r2_matches[0]
    if r1_matches: return r1_matches[0], None
    return None, None


def find_module_5_scaffolds(config, sample_name, summary):
    sample_dir = get_module_5_dir(config, sample_name)
    if not os.path.isdir(sample_dir):
        raise FileNotFoundError(f"Module 5 output directory not found for {sample_name}: {sample_dir}")

    if summary.get("is_segmented") is True:
        selected_segments = summary.get("selected_segments", [])
        tasks = []
        for seg in selected_segments:
            seg_id = seg.get("segment_seq_id")
            if not seg_id: continue
            patterns = [
                os.path.join(sample_dir, f"{sample_name}_seg_{seg_id}_scaffold.fasta"),
                os.path.join(sample_dir, f"{sample_name}_seg_{seg_id.replace('.', '_')}_scaffold.fasta"),
                os.path.join(sample_dir, f"{sample_name}_seg_*{seg_id}*_scaffold.fasta"),
            ]
            found = None
            for pattern in patterns:
                matches = sorted(glob.glob(pattern))
                if matches:
                    found = matches[0]
                    break
            if found:
                tasks.append({"segment_id": seg_id, "scaffold_fasta": found})

        if tasks: return tasks

        combined_path = os.path.join(sample_dir, f"{sample_name}_segmented_scaffold.fasta")
        if os.path.exists(combined_path):
            return [{"segment_id": "combined", "scaffold_fasta": combined_path}]

        matches = sorted(glob.glob(os.path.join(sample_dir, f"{sample_name}_seg_*_scaffold.fasta")))
        if matches:
            return [{"segment_id": os.path.basename(match).split("_scaffold.fasta")[0].split("_seg_")[-1], "scaffold_fasta": match} for match in matches]

        raise FileNotFoundError(f"No segmented scaffold file found for {sample_name} in {sample_dir}")

    patterns = [
        os.path.join(sample_dir, f"{sample_name}_scaffold.fasta"),
        os.path.join(sample_dir, f"{sample_name}_nucmer_scaffold.fasta"),
        os.path.join(sample_dir, f"{sample_name}_scaffolded.fasta"),
        os.path.join(sample_dir, f"{sample_name}*.fasta"),
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return [{"segment_id": "all", "scaffold_fasta": matches[0]}]

    raise FileNotFoundError(f"No module 5 scaffold FASTA found for {sample_name} under {sample_dir}")


# ---------------------------------------------------------------------------
# Sequence polishing using minimap2 + samtools + iVar (Two-Tiered)
# ---------------------------------------------------------------------------

def run_single_ivar_pass(sorted_bam, prefix, params, samtools_bin, ivar_bin, sample_name):
    """Executes a single samtools mpileup | ivar consensus pipeline."""
    consensus_fasta = f"{prefix}.fa"
    
    mpileup_cmd = [samtools_bin, "mpileup", "-aa", "-A", "-d", "1000", "-Q", str(params["min_qual"]), sorted_bam]
    ivar_cmd = [
        ivar_bin, "consensus", "-p", prefix,
        "-t", str(params["min_freq"]), "-m", str(params["min_depth"]),
        "-q", str(params["min_qual"]), "-n", "N"
    ]

    p_mpileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p_ivar = subprocess.Popen(ivar_cmd, stdin=p_mpileup.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    p_mpileup.stdout.close()
    _, p_ivar_err = p_ivar.communicate()
    p_mpileup.wait()

    if p_mpileup.returncode != 0 or p_ivar.returncode != 0:
        raise RuntimeError(f"iVar pass failed for {sample_name}: {p_ivar_err.decode(errors='replace')}")
        
    if not os.path.exists(consensus_fasta):
        raise RuntimeError(f"iVar failed to produce {consensus_fasta}")
        
    return consensus_fasta


def merge_consensus_fastas(strict_fasta, relaxed_fasta, output_fasta, acc):
    """Merges two consensus FASTAs, preferring the strict call but rescuing Ns with the relaxed call."""
    strict_seq = str(next(SeqIO.parse(strict_fasta, "fasta")).seq).upper()
    relaxed_seq = str(next(SeqIO.parse(relaxed_fasta, "fasta")).seq).upper()
    
    rescued_count = 0
    merged_seq = []
    
    # Handle potential length mismatches gracefully
    min_len = min(len(strict_seq), len(relaxed_seq))
    for i in range(min_len):
        s_base = strict_seq[i]
        r_base = relaxed_seq[i]
        if s_base != "N":
            merged_seq.append(s_base)
        else:
            if r_base != "N":
                merged_seq.append(r_base)
                rescued_count += 1
            else:
                merged_seq.append("N")
                
    # Append any trailing bases from the longer sequence
    if len(strict_seq) > min_len:
        merged_seq.append(strict_seq[min_len:])
    elif len(relaxed_seq) > min_len:
        merged_seq.append(relaxed_seq[min_len:])
                
    final_seq = "".join(merged_seq)
    
    with open(output_fasta, "w") as handle:
        handle.write(f">{acc}_final_polished_consensus\n{final_seq}\n")
        
    print(f"[RESCUE SUCCESS]: Filled {rescued_count} additional N-gaps using relaxed parameters for {acc}.")
    return output_fasta


def analyze_base_transitions(scaffold_fasta, polished_fasta):
    """Compares unpolished scaffold to polished consensus to track N-base transitions."""
    scaffold_seqs = [str(rec.seq).upper() for rec in SeqIO.parse(scaffold_fasta, "fasta")]
    polished_seqs = [str(rec.seq).upper() for rec in SeqIO.parse(polished_fasta, "fasta")]
    
    scaffold_seq = "".join(scaffold_seqs)
    polished_seq = "".join(polished_seqs)
    
    # Truncate to min_len to ensure perfectly consistent math
    min_len = min(len(scaffold_seq), len(polished_seq))
    if min_len == 0:
        return {
            "unpolished_n": 0, "polished_n": 0, "n_to_canonical": 0,
            "canonical_to_n": 0, "net_n_change": 0
        }
        
    scaffold_seq = scaffold_seq[:min_len]
    polished_seq = polished_seq[:min_len]
    
    unpolished_n = scaffold_seq.count('N')
    polished_n = polished_seq.count('N')
    
    n_to_canonical = sum(1 for s, p in zip(scaffold_seq, polished_seq) if s == 'N' and p != 'N')
    canonical_to_n = sum(1 for s, p in zip(scaffold_seq, polished_seq) if s != 'N' and p == 'N')
            
    return {
        "unpolished_n": unpolished_n,
        "polished_n": polished_n,
        "n_to_canonical": n_to_canonical,
        "canonical_to_n": canonical_to_n,
        "net_n_change": polished_n - unpolished_n
    }


def run_ivar_polishing_pipeline(config, sample_name, scaffold_fasta, read_1, read_2, output_dir, label):
    minimap2_bin = resolve_executable_path(config, "minimap2_bin_path")
    samtools_bin = resolve_executable_path(config, "samtools_bin_path")
    ivar_bin = resolve_executable_path(config, "ivar_bin_path")

    threads = int(get_config_value(config, "module_6_params", "threads", get_config_value(config, "module_5_params", "threads", 8)))
    
    sorted_bam = os.path.join(output_dir, f"{label}_aligned_sorted.bam")
    sort_tmp_dir = os.path.join(output_dir, "samtools_sort_tmp")
    os.makedirs(sort_tmp_dir, exist_ok=True)

    # 1. Alignment and Sorting
    minimap_cmd = [minimap2_bin, "-ax", "sr", scaffold_fasta]
    if read_2: minimap_cmd.extend([read_1, read_2])
    else: minimap_cmd.append(read_1)

    samtools_view_cmd = [samtools_bin, "view", "-bS", "-"]
    samtools_sort_cmd = [samtools_bin, "sort", "-@", str(threads), "-T", sort_tmp_dir, "-o", sorted_bam, "-"]

    print(f"[ALIGNMENT]: Mapping reads for {sample_name} -> {scaffold_fasta}")
    p1 = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    p1.stdout.close()
    p2.stdout.close()
    _, p3_err = p3.communicate()
    p2.wait()
    p1.wait()

    if p1.returncode != 0 or p2.returncode != 0 or p3.returncode != 0:
        raise RuntimeError(f"Alignment and sorting failed for {sample_name}: {p3_err.decode(errors='replace')}")

    print(f"[ALIGNMENT]: Indexing BAM: {sorted_bam}")
    subprocess.run([samtools_bin, "index", sorted_bam], check=True, stderr=subprocess.PIPE)

    # 2. Two-Tiered Consensus Calling
    strict_prefix = os.path.join(output_dir, f"{label}_strict_consensus")
    relaxed_prefix = os.path.join(output_dir, f"{label}_relaxed_consensus")
    final_fasta = os.path.join(output_dir, f"{label}_final_polished_consensus.fa")
    
    print(f"[IVAR TIER 1]: Running strict consensus (m={STRICT_PARAMS['min_depth']}, t={STRICT_PARAMS['min_freq']}, q={STRICT_PARAMS['min_qual']})...")
    strict_fasta = run_single_ivar_pass(sorted_bam, strict_prefix, STRICT_PARAMS, samtools_bin, ivar_bin, sample_name)
    
    print(f"[IVAR TIER 2]: Running relaxed rescue consensus (m={RELAXED_PARAMS['min_depth']}, t={RELAXED_PARAMS['min_freq']}, q={RELAXED_PARAMS['min_qual']})...")
    relaxed_fasta = run_single_ivar_pass(sorted_bam, relaxed_prefix, RELAXED_PARAMS, samtools_bin, ivar_bin, sample_name)
    
    # 3. Merge and Rescue Gaps
    print(f"[MERGE]: Combining strict and relaxed consensus to maximize gap-filling...")
    final_fasta = merge_consensus_fastas(strict_fasta, relaxed_fasta, final_fasta, label)
    
    # Cleanup intermediate files
    for f in [strict_fasta, relaxed_fasta, f"{strict_prefix}.qual.txt", f"{relaxed_prefix}.qual.txt"]:
        if os.path.exists(f): os.remove(f)

    return final_fasta


# ---------------------------------------------------------------------------
# CheckV summary parsing
# ---------------------------------------------------------------------------

def clean_checkv_field(value):
    """Aggressively cleans CheckV fields to prevent TSV formatting corruption."""
    if not value: return "N/A"
    return str(value).strip().replace('\n', ' ').replace('\t', ' ').replace('\r', '')

def run_checkv_qc(config, polished_fasta, output_dir):
    checkv_bin = resolve_executable_path(config, "checkv_bin_path")
    checkv_db_dir = get_config_value(config, "module_6_params", "checkv_db_dir", get_config_value(config, "module_5_params", "checkv_db_dir", None))
    threads = int(get_config_value(config, "module_6_params", "threads", get_config_value(config, "module_5_params", "threads", 8)))
    os.makedirs(output_dir, exist_ok=True)

    candidate_threads = [threads] if threads <= 8 else [8, max(1, threads // 2), threads]
    summary_path = os.path.join(output_dir, "quality_summary.tsv")

    for threads_use in candidate_threads:
        cmd = [checkv_bin, "end_to_end", polished_fasta, output_dir, "-t", str(threads_use), "--restart"]
        if checkv_db_dir and os.path.exists(checkv_db_dir):
            cmd.extend(["-d", checkv_db_dir])

        print(f"[CHECKV]: Running CheckV for {polished_fasta} with {threads_use} threads")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0 and os.path.exists(summary_path):
            try:
                with open(summary_path, "r") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    for row in reader:
                        completeness = clean_checkv_field(row.get("completeness") or row.get("checkv_completeness"))
                        quality = clean_checkv_field(row.get("quality") or row.get("checkv_quality"))
                        kmer_freq = clean_checkv_field(row.get("kmer_freq") or row.get("kmer"))
                        warnings = clean_checkv_field(row.get("warnings"))
                        if warnings == "N/A" or warnings == "": warnings = "None"
                        
                        print(f"[CHECKV RESULTS]: Completeness={completeness}% | Quality={quality} | kmer_freq={kmer_freq}")
                        return completeness, quality, kmer_freq, warnings
            except Exception: pass

        if proc.returncode != 0:
            print(f"[CHECKV WARN]: CheckV attempt failed with {threads_use} threads: {proc.stderr.strip()}")

    return "N/A", "N/A", "N/A", "CheckV failed to generate a quality summary" 


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_qc_report(sample_name, scaffold_fasta, polished_fasta, completeness, quality, kmer_freq, warnings, report_path, n_stats):
    with open(report_path, "w") as handle:
        handle.write(f"==================================================\n")
        handle.write(f" iVar Two-Tiered Consensus Polishing Report: {sample_name}\n")
        handle.write(f"==================================================\n")
        handle.write(f"Input scaffold FASTA: {scaffold_fasta}\n")
        handle.write(f"Polished FASTA:       {polished_fasta}\n\n")
        
        handle.write(f"--- N-Base Transition Analysis ---\n")
        handle.write(f"Unpolished Scaffold N-bases: {n_stats['unpolished_n']}\n")
        handle.write(f"Final Polished N-bases:    {n_stats['polished_n']}\n")
        handle.write(f"Net Change in N-bases:     {n_stats['net_n_change']:+d}\n\n")
        
        handle.write(f"--- Gap Filling & Loss Metrics ---\n")
        handle.write(f"Ns rescued to Canonical:   {n_stats['n_to_canonical']} (Gaps successfully filled by reads)\n")
        handle.write(f"Canonical lost to N:       {n_stats['canonical_to_n']} (Coverage dropped below thresholds)\n\n")
        
        handle.write(f"--- CheckV Metrics ---\n")
        handle.write(f"Completeness: {completeness}%\n")
        handle.write(f"Quality:      {quality}\n")
        handle.write(f"K-mer Freq:   {kmer_freq}\n")
        handle.write(f"Warnings:     {warnings}\n")
    print(f"[REPORT]: Wrote {report_path}")
    return report_path


# ---------------------------------------------------------------------------
# Sample processing
# ---------------------------------------------------------------------------

def process_single_sample(config, sample_name):
    summary_path = get_module_4_summary_path(config, sample_name)
    if not os.path.exists(summary_path):
        raise FileNotFoundError(f"Module 4 summary not found for {sample_name}: {summary_path}")

    with open(summary_path, "r") as handle:
        summary = json.load(handle)

    read_1, read_2 = locate_fastq_pair(config, sample_name)
    if not read_1:
        raise FileNotFoundError(f"No trimmed paired reads found for {sample_name} in module 1 fastp output or input read directory.")

    tasks = find_module_5_scaffolds(config, sample_name, summary)
    sample_out_dir = get_module_6_dir(config, sample_name)
    os.makedirs(sample_out_dir, exist_ok=True)

    rows = []
    combined_records = []

    # 1. Polish all segments/scaffolds independently
    for task in tasks:
        label = task["segment_id"]
        scaffold_fasta = task["scaffold_fasta"]
        file_label = f"{sample_name}_seg_{label}" if label != "all" else sample_name
        out_dir = os.path.join(sample_out_dir, file_label)
        os.makedirs(out_dir, exist_ok=True)

        polished_fasta = run_ivar_polishing_pipeline(
            config, sample_name, scaffold_fasta, read_1, read_2, out_dir, file_label
        )
        
        task["polished_fasta"] = polished_fasta
        task["out_dir"] = out_dir
        task["file_label"] = file_label
        
        for record in SeqIO.parse(polished_fasta, "fasta"):
            record.id = f"{sample_name}_{label}"
            record.description = ""
            combined_records.append(record)

    # 2. Write combined polished FASTA
    combined_polished_path = os.path.join(sample_out_dir, f"{sample_name}_consensus_combined.fasta")
    if combined_records:
        SeqIO.write(combined_records, combined_polished_path, "fasta")
        print(f"[COMBINED]: Wrote combined polished consensus FASTA to {combined_polished_path}")

    # 3. Run CheckV on the COMBINED polished genome
    # This solves the "no viral genes detected" issue for segmented viruses
    checkv_dir = os.path.join(sample_out_dir, "checkv_combined")
    completeness, quality, kmer_freq, warnings = run_checkv_qc(config, combined_polished_path, checkv_dir)

    # 4. Generate reports and rows using the combined CheckV stats
    for task in tasks:
        label = task["segment_id"]
        scaffold_fasta = task["scaffold_fasta"]
        polished_fasta = task["polished_fasta"]
        out_dir = task["out_dir"]
        file_label = task["file_label"]
        
        n_stats = analyze_base_transitions(scaffold_fasta, polished_fasta)
        
        report_path = os.path.join(out_dir, f"{file_label}_ivar_polished_report.txt")
        write_qc_report(sample_name, scaffold_fasta, polished_fasta, completeness, quality, kmer_freq, warnings, report_path, n_stats)

        row = {
            "sample_accession": sample_name,
            "is_segmented": bool(summary.get("is_segmented")),
            "segment_id": label,
            "input_scaffold_fasta": scaffold_fasta,
            "polished_fasta": polished_fasta,
            "unpolished_n_bases": n_stats["unpolished_n"],
            "polished_n_bases": n_stats["polished_n"],
            "n_rescued_to_canonical": n_stats["n_to_canonical"],
            "canonical_lost_to_n": n_stats["canonical_to_n"],
            "net_n_change": n_stats["net_n_change"],
            "checkv_completeness": completeness,
            "checkv_quality": quality,
            "checkv_kmer_freq": kmer_freq,
            "checkv_warnings": warnings,
            "report_path": report_path,
        }
        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Aggregated summary writing
# ---------------------------------------------------------------------------

def write_consensus_summary(config, rows):
    module_6_dir = os.path.join(config.cross_module_params["output_dir"], "module_6")
    os.makedirs(module_6_dir, exist_ok=True)
    summary_path = os.path.join(module_6_dir, "Identified_consensus.tsv")

    fieldnames = [
        "sample_accession",
        "is_segmented",
        "segment_id",
        "input_scaffold_fasta",
        "polished_fasta",
        "unpolished_n_bases",
        "polished_n_bases",
        "n_rescued_to_canonical",
        "canonical_lost_to_n",
        "net_n_change",
        "checkv_completeness",
        "checkv_quality",
        "checkv_kmer_freq",
        "checkv_warnings",
        "report_path",
    ]

    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    print(f"[SUMMARY]: Wrote polished consensus table to {summary_path}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Map reads back to scaffold outputs and polish them with iVar consensus calling.")
    parser.add_argument("--config", required=True, help="Path to the pipeline config Python file.")
    args = parser.parse_args()

    config = load_config(args.config)
    all_rows = []

    print(f"\nStarting module_6_consensus_polish...\nConfig file: {args.config}\n")

    for sample_name in config.cross_module_params["samples"]:
        try:
            rows = process_single_sample(config, sample_name)
            all_rows.extend(rows)
        except Exception as exc:
            print(f"[ERROR]: Failed for sample {sample_name}: {exc}")

    write_consensus_summary(config, all_rows)
    print("\n[module_6_consensus_polish] Finished all samples.")


if __name__ == "__main__":
    main()
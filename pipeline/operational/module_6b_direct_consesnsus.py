import argparse
import csv
import importlib.util
import os
import shutil
import subprocess
from collections import Counter
from Bio import SeqIO

"""
Module 6B: Direct reference mapping consensus generation.
Mimics the CZ ID pipeline approach: map reads directly to a reference, 
call consensus with iVar (via samtools mpileup), and evaluate completeness/quality.
"""

def load_config(config_path):
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
    value = get_config_value(config, "module_6b_params", key, None)
    if value is None:
        raise FileNotFoundError(f"Executable '{key}' is not configured in module_6b_params.")
    if os.path.exists(value):
        return value
    resolved = shutil.which(value)
    if resolved:
        return resolved
    raise FileNotFoundError(f"Configured executable '{value}' for '{key}' was not found.")

def get_module_1_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_1", "fastp", sample_name)

def get_module_6b_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_6b", sample_name)

def find_trimmed_reads(sample_dir, sample_name, read_1_str, read_2_str, singleton_str):
    read1 = read2 = singleton = None
    for file_name in sorted(os.listdir(sample_dir)):
        full_path = os.path.join(sample_dir, file_name)
        if not os.path.isfile(full_path) or not file_name.endswith((".fastq", ".fastq.gz")):
            continue
        if read_1_str in file_name:
            read1 = full_path
        elif read_2_str in file_name:
            read2 = full_path
        elif singleton_str is not None and singleton_str in file_name:
            if "_R1" not in file_name and "_R2" not in file_name:
                singleton = full_path

    if read1 is None or read2 is None:
        for file_name in sorted(os.listdir(sample_dir)):
            full_path = os.path.join(sample_dir, file_name)
            if not os.path.isfile(full_path) or not file_name.endswith((".fastq", ".fastq.gz")):
                continue
            if file_name.startswith(f"{sample_name}_R1"):
                read1 = full_path
            elif file_name.startswith(f"{sample_name}_R2"):
                read2 = full_path
            elif file_name.startswith(f"{sample_name}_singleton"):
                singleton = full_path
    return read1, read2, singleton

def run_minimap2_and_ivar(config, sample_name, ref_fasta, read1, read2, singleton, sample_dir):
    minimap2_bin = resolve_executable_path(config, "minimap2_bin_path")
    samtools_bin = resolve_executable_path(config, "samtools_bin_path")
    ivar_bin = resolve_executable_path(config, "ivar_bin_path")
    threads = get_config_value(config, "module_6b_params", "threads", 4)
    
    ivar_min_qual = get_config_value(config, "module_6b_params", "ivar_min_quality", 20)
    ivar_min_freq = get_config_value(config, "module_6b_params", "ivar_min_freq", 0.5)
    ivar_min_depth = get_config_value(config, "module_6b_params", "ivar_min_depth", 10)

    os.makedirs(sample_dir, exist_ok=True)
    prefix = os.path.join(sample_dir, f"{sample_name}_direct_consensus")
    bam_file = f"{prefix}.bam"
    
    print(f"[MODULE 6B]: Mapping reads to reference for {sample_name}...")
    mm2_cmd = f"{minimap2_bin} -a -x sr -t {threads} {ref_fasta} {read1} {read2}"
    if singleton:
        mm2_cmd += f" {singleton}"
    
    sort_cmd = f"{samtools_bin} sort -@ {threads} -o {bam_file}"
    full_cmd = f"{mm2_cmd} | {sort_cmd}"
    subprocess.run(full_cmd, shell=True, check=True, executable='/bin/bash')
    
    print(f"[MODULE 6B]: Indexing BAM for {sample_name}...")
    subprocess.run([samtools_bin, "index", bam_file], check=True)
    
    print(f"[MODULE 6B]: Calling consensus with iVar for {sample_name}...")
    mpileup_cmd = f"{samtools_bin} mpileup -aa -A -d 0 -Q 0 -f {ref_fasta} {bam_file}"
    ivar_cmd = f"{ivar_bin} consensus -p {prefix} -q {ivar_min_qual} -t {ivar_min_freq} -m {ivar_min_depth}"
    full_ivar_cmd = f"{mpileup_cmd} | {ivar_cmd}"
    
    result = subprocess.run(full_ivar_cmd, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
    
    consensus_fasta = f"{prefix}.fasta"
    consensus_fa = f"{prefix}.fa"
    
    # Handle iVar's default .fa extension
    if os.path.exists(consensus_fa):
        print(f"[MODULE 6B]: iVar produced .fa. Renaming to .fasta for consistency.")
        shutil.move(consensus_fa, consensus_fasta)
        return consensus_fasta, bam_file
    elif os.path.exists(consensus_fasta):
        return consensus_fasta, bam_file
    else:
        print(f"[ERROR]: iVar did not produce consensus file.")
        print(f"[DEBUG] iVar stderr:\n{result.stderr}")
        raise RuntimeError(f"iVar failed to produce consensus FASTA for {sample_name}.")

def evaluate_sequence(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return {"length": 0, "n_count": 0, "non_canonical_count": 0, "n_pct": 0.0}
    
    seq = str(records[0].seq).upper()
    length = len(seq)
    counts = Counter(seq)
    
    n_count = counts.get("N", 0)
    canonical = {"A", "C", "G", "T", "N"}
    non_canonical_count = sum(v for k, v in counts.items() if k not in canonical)
    
    n_pct = (n_count / length * 100) if length > 0 else 0.0
    
    return {
        "length": length,
        "n_count": n_count,
        "non_canonical_count": non_canonical_count,
        "n_pct": round(n_pct, 2)
    }

def run_checkv_qc(config, target_fasta, checkv_dir):
    checkv_bin = resolve_executable_path(config, "checkv_bin_path")
    threads = get_config_value(config, "module_6b_params", "threads", 4)
    checkv_db_dir = get_config_value(config, "module_6b_params", "checkv_db_dir", None)
    os.makedirs(checkv_dir, exist_ok=True)
    
    cmd = [checkv_bin, "end_to_end", target_fasta, checkv_dir, "-t", str(threads), "--restart"]
    if checkv_db_dir and os.path.exists(checkv_db_dir):
        cmd.extend(["-d", checkv_db_dir])
        
    print(f"[CHECKV]: Running CheckV for {target_fasta} with {threads} threads")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"[CHECKV WARN]: CheckV failed for {target_fasta}. Stderr: {proc.stderr}")
        return "N/A", "N/A", "N/A", proc.stderr.strip()
        
    summary_tsv = None
    for root, _, files in os.walk(checkv_dir):
        if "quality_summary.tsv" in files:
            summary_tsv = os.path.join(root, "quality_summary.tsv")
            break
            
    if not summary_tsv or not os.path.exists(summary_tsv):
        return "N/A", "N/A", "N/A", "quality_summary.tsv missing"
        
    with open(summary_tsv, "r") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        
    if not rows:
        return "N/A", "N/A", "N/A", "quality_summary.tsv empty"
        
    comp_str = str(rows[0].get("completeness", "N/A")).strip()
    quality = str(rows[0].get("checkv_quality", rows[0].get("quality", "N/A"))).strip()
    kmer_str = str(rows[0].get("kmer_freq", "N/A")).strip()
    warnings = str(rows[0].get("warnings", "None")).strip()
    
    return comp_str, quality, kmer_str, warnings

def process_sample(sample_name, config):
    ref_genomes = get_config_value(config, "module_6b_params", "reference_genomes", {})
    if sample_name not in ref_genomes:
        print(f"[SKIP]: Sample {sample_name} has no configured reference genome in module_6b_params.")
        return None
        
    ref_fasta = ref_genomes[sample_name]
    if not os.path.exists(ref_fasta):
        raise FileNotFoundError(f"Reference genome not found for {sample_name}: {ref_fasta}")
        
    module_1_dir = get_module_1_dir(config, sample_name)
    if not os.path.exists(module_1_dir):
        raise FileNotFoundError(f"Module 1 output directory not found for {sample_name}: {module_1_dir}")
        
    read1, read2, singleton = find_trimmed_reads(
        module_1_dir, sample_name,
        config.cross_module_params["read_1_str"],
        config.cross_module_params["read_2_str"],
        config.cross_module_params.get("singleton_str")
    )
    
    if read1 is None or read2 is None:
        raise FileNotFoundError(f"Trimmed reads not found for {sample_name} in {module_1_dir}")
        
    sample_dir = get_module_6b_dir(config, sample_name)
    
    consensus_fasta, bam_file = run_minimap2_and_ivar(config, sample_name, ref_fasta, read1, read2, singleton, sample_dir)
    seq_stats = evaluate_sequence(consensus_fasta)
    
    checkv_dir = os.path.join(sample_dir, "checkv")
    
    # Prevent CheckV crash on low-coverage (high N%) genomes
    if seq_stats["n_pct"] > 50.0:
        print(f"[CHECKV SKIP]: Skipping CheckV for {sample_name} because genome is {seq_stats['n_pct']}% Ns.")
        comp, quality, kmer, warnings = "N/A", "N/A", "N/A", "Skipped: >50% Ns"
    else:
        comp, quality, kmer, warnings = run_checkv_qc(config, consensus_fasta, checkv_dir)
    
    return {
        "sample_accession": sample_name,
        "reference_fasta": ref_fasta,
        "consensus_fasta": consensus_fasta,
        "length": seq_stats["length"],
        "n_count": seq_stats["n_count"],
        "n_pct": seq_stats["n_pct"],
        "non_canonical_count": seq_stats["non_canonical_count"],
        "checkv_completeness": comp,
        "checkv_quality": quality,
        "checkv_kmer_freq": kmer,
        "checkv_warnings": warnings
    }

def validate_inputs(config):
    ref_genomes = get_config_value(config, "module_6b_params", "reference_genomes", {})
    for sample_name in config.cross_module_params["samples"]:
        if sample_name in ref_genomes:
            ref_fasta = ref_genomes[sample_name]
            if not os.path.exists(ref_fasta):
                raise FileNotFoundError(f"Reference genome not found for {sample_name}: {ref_fasta}")
            print(f"[VALIDATION]: Sample {sample_name} has valid reference genome: {ref_fasta}")
    print("Module 6B inputs validated successfully.")

def main():
    parser = argparse.ArgumentParser(description="Direct reference mapping consensus generation (Module 6B).")
    parser.add_argument("--config", required=True, help="Path to the config Python file.")
    args = parser.parse_args()
    
    config = load_config(args.config)
    print(f"\nStarting module_6b_consensus...\nConfig file location: {args.config}\n")
    
    validate_inputs(config)
    
    summary_rows = []
    for sample_name in config.cross_module_params["samples"]:
        try:
            result = process_sample(sample_name, config)
            if result:
                summary_rows.append(result)
        except Exception as exc:
            print(f"[ERROR]: Failed for sample {sample_name}: {exc}")
            
    if summary_rows:
        module_6b_dir = os.path.join(config.cross_module_params["output_dir"], "module_6b")
        os.makedirs(module_6b_dir, exist_ok=True)
        summary_path = os.path.join(module_6b_dir, "direct_consensus_summary.tsv")
        fieldnames = ["sample_accession", "reference_fasta", "consensus_fasta", "length", 
                      "n_count", "n_pct", "non_canonical_count", "checkv_completeness", 
                      "checkv_quality", "checkv_kmer_freq", "checkv_warnings"]
        with open(summary_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in summary_rows:
                writer.writerow({key: row.get(key, "") for key in fieldnames})
        print(f"[SUMMARY]: Wrote direct consensus summary to {summary_path}")
        
    print("\n[module_6b_consensus] Finished all samples.")

if __name__ == "__main__":
    main()
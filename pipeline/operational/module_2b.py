import argparse
import importlib.util
import os
import subprocess

"""
Module 2B: Reference-guided read extraction.

This module branches from module_1 QC-ed reads. For samples with a configured 
reference genome, it aligns reads permissively using minimap2 and extracts 
the mapped reads into the module_2 sample-specific output directory.
"""

def load_config(config_path):
    """Load a Python config file as a Python module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def find_trimmed_reads(sample_dir, sample_name, read_1_str, read_2_str, singleton_str):
    """Locate paired-end and singleton trimmed FASTQ files for a sample."""
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

    # Fallback to explicit naming if substring matching fails
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


def process_sample(sample_name, config):
    """Align and extract reads for a single sample."""
    ref_genomes = config.module_2b_params.get("reference_genomes", {})
    if sample_name not in ref_genomes:
        print(f"[SKIP]: Sample {sample_name} has no configured reference genome in module_2b_params.")
        return

    ref_fa = ref_genomes[sample_name]
    if not os.path.exists(ref_fa):
        raise FileNotFoundError(f"Reference genome not found for {sample_name}: {ref_fa}")

    module_1_dir = os.path.join(config.cross_module_params["output_dir"], "module_1", "fastp", sample_name)
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

    # SAFEGUARD 1: Check if files are empty (0 bytes). If fastp failed, minimap2 will output nothing.
    if os.path.getsize(read1) == 0 or os.path.getsize(read2) == 0:
        print(f"[WARNING]: Input FASTQ files for {sample_name} are empty (0 bytes). Skipping alignment.")
        return

    # Output to module_2 directory to keep read subsets consolidated
    out_dir = os.path.join(config.cross_module_params["output_dir"], "module_2", sample_name)
    os.makedirs(out_dir, exist_ok=True)

    out_r1 = os.path.join(out_dir, f"{sample_name}_aligned_to_ref_R1.fastq")
    out_r2 = os.path.join(out_dir, f"{sample_name}_aligned_to_ref_R2.fastq")
    out_s = os.path.join(out_dir, f"{sample_name}_aligned_to_ref_singleton.fastq")

    minimap2_bin = config.module_2b_params.get("minimap2_bin_path", "minimap2")
    samtools_bin = config.module_2b_params.get("samtools_bin_path", "samtools")
    threads = config.module_2b_params.get("threads", config.cross_module_params.get("threads", 4))

    # GREEDIER MINIMAP2 SETTINGS:
    # -a: output SAM (required for pipe)
    # -x sr: short read preset
    # -k 11: minimum k-mer size for high sensitivity to divergence
    # -w 3: smaller minimizer window to find more seeds
    # -s 40: lower alignment score threshold to accept weaker/shorter alignments
    # --secondary=no: faster, we only care if the read maps at all
    mm2_cmd = (
        f"{minimap2_bin} -a -x sr -k 11 -w 3 -s 40 --secondary=no -t {threads} "
        f"{ref_fa} {read1} {read2}"
    )
    
    # SAFEGUARD 2: Explicit '-' tells samtools to read from the pipe, preventing parsing ambiguity
    cmd = (
        f"{mm2_cmd} | "
        f"{samtools_bin} view -b -F 4 - | "
        f"{samtools_bin} sort -n -m 1G - | "
        f"{samtools_bin} fastq -1 {out_r1} -2 {out_r2} -s {out_s} -"
    )
    
    print(f"[MODULE 2B]: Aligning and extracting reads for {sample_name} (greedy mode)...")
    print(f"[MODULE 2B]: {cmd}")
    
    # SAFEGUARD 3: Capture stderr so we can see the EXACT minimap2/samtools error if it fails
    result = subprocess.run(cmd, shell=True, executable='/bin/bash', capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"[ERROR]: Minimap2/Samtools pipeline failed for {sample_name}.")
        print(f"--- STDERR OUTPUT ---")
        print(result.stderr.strip())
        print(f"---------------------")
        raise RuntimeError(f"Command failed with return code {result.returncode}")
        
    print(f"[SUCCESS]: Wrote aligned read subsets for {sample_name} to {out_dir}")


def validate_inputs(config):
    """Validate reference genomes before processing."""
    ref_genomes = config.module_2b_params.get("reference_genomes", {})
    for sample_name in config.cross_module_params["samples"]:
        if sample_name in ref_genomes:
            ref_fa = ref_genomes[sample_name]
            if not os.path.exists(ref_fa):
                raise FileNotFoundError(f"Reference genome not found for {sample_name}: {ref_fa}")
            print(f"[VALIDATION]: Sample {sample_name} has valid reference genome: {ref_fa}")
    print("Module 2B inputs validated successfully.")


def main():
    parser = argparse.ArgumentParser(description="Extract reads aligning to a sample-specific reference genome.")
    parser.add_argument("--config", required=True, help="Path to the config.py file")
    args = parser.parse_args()

    config = load_config(args.config)
    print(f"\nStarting module_2b_reference_extract...\nConfig: {args.config}\n")
    
    validate_inputs(config)
    
    for sample_name in config.cross_module_params["samples"]:
        try:
            process_sample(sample_name, config)
        except Exception as exc:
            print(f"[ERROR]: Failed to process sample {sample_name}: {exc}")
            
    print("\n[module_2b_reference_extract] Finished.")


if __name__ == "__main__":
    main()
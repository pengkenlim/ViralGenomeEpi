import argparse
import importlib.util
import os
import subprocess
from pathlib import Path

"""
Module 3: de novo assembly of viral read subsets from module 2.

This module follows the same operational conventions as module_1_qc_trim.py and
module_2_classify_scrub.py:
- config-driven execution with --config,
- validation of configured sample inputs,
- clear per-sample logging,
- subprocess-based execution of the external assembler.

Assemblies are written into nested directories so the output mode is explicit:
    <output_dir>/module_3/viral/<sample>
    <output_dir>/module_3/dominant_virus/<sample>
"""


def load_config(config_path):
    """Load a Python config file as a Python module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def get_reads_subset_mode(config):
    """Normalize a configured reads_subset value into a module-3 mode."""
    subset_value = str(config.module_3_params.get("reads_subset", "viral")).strip().lower()

    aliases = {
        "viral": "viral",
        "viral_plus_unclassified": "viral",
        "viral plus unclassified": "viral",
        "dominant virus only": "dominant_virus",
        "dominant_virus_only": "dominant_virus",
        "dominant virus": "dominant_virus",
    }

    if subset_value not in aliases:
        valid = ", ".join(sorted(aliases))
        raise ValueError(
            f"Unsupported module_3 reads_subset value: '{config.module_3_params.get('reads_subset')}'. "
            f"Expected one of: {valid}"
        )

    return aliases[subset_value]


def locate_subset_reads(sample_dir, sample_name, subset_mode):
    """Locate the paired-end and singleton read files from module 2 for a given subset mode."""
    if subset_mode == "viral":
        read_prefix = "viral_plus_unclassified"
    else:
        read_prefix = "dominant_virus_plus_unclassified"

    read1 = None
    read2 = None
    singleton = None

    for file_name in sorted(os.listdir(sample_dir)):
        full_path = os.path.join(sample_dir, file_name)
        if not os.path.isfile(full_path):
            continue
        if not file_name.endswith((".fastq", ".fastq.gz")):
            continue

        if file_name.startswith(f"{sample_name}_{read_prefix}_R1"):
            read1 = full_path
        elif file_name.startswith(f"{sample_name}_{read_prefix}_R2"):
            read2 = full_path
        elif file_name.startswith(f"{sample_name}_{read_prefix}_singleton"):
            singleton = full_path

    return read1, read2, singleton


def get_thread_count(config):
    """Return the effective thread count for the assembly module."""
    module_threads = config.module_3_params.get("threads")
    if module_threads is not None:
        return int(module_threads)
    return int(config.cross_module_params.get("threads", 4))


def get_memory_limit(config):
    """Return the memory limit for the assembly module."""
    return int(config.module_3_params.get("memory", 500))


def get_assembly_output_dir(config, subset_mode, sample_name):
    """Return the output directory for assembly results for one sample and run mode."""
    mode_dir = "viral" if subset_mode == "viral" else "dominant_virus"
    return os.path.join(config.cross_module_params["output_dir"], "module_3", mode_dir, sample_name)


def run_rnaviralspades(sample_name, read1, read2, singleton, config, subset_mode):
    """Run rnaviralspades.py on the selected read subset for one sample."""
    bin_path = config.module_3_params["rnaviralspades_bin_path"]
    if not os.path.exists(bin_path):
        raise FileNotFoundError(f"rnaviralspades executable not found at {bin_path}")

    output_dir = get_assembly_output_dir(config, subset_mode, sample_name)
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        bin_path,
        "-o",
        str(output_dir),
        "-t",
        str(get_thread_count(config)),
        "-m",
        str(get_memory_limit(config)),
    ]

    if config.module_3_params.get("use_rna_mode", True):
        cmd.append("--rna")

    if read1 is not None and read2 is not None:
        cmd.extend(["-1", str(read1), "-2", str(read2)])
    if singleton is not None:
        cmd.extend(["-s", str(singleton)])

    print(f"[ASSEMBLY]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    # Verify that the expected output files exist.
    scaffolds_file = os.path.join(output_dir, "scaffolds.fasta")
    contigs_file = os.path.join(output_dir, "contigs.fasta")
    if os.path.exists(scaffolds_file):
        print(f"[SUCCESS]: Wrote scaffolds for {sample_name} to {scaffolds_file}")
    elif os.path.exists(contigs_file):
        print(f"[SUCCESS]: Wrote contigs for {sample_name} to {contigs_file}")
    else:
        print(f"[WARNING]: No scaffolds or contigs were produced for {sample_name} in {output_dir}")


def validate_inputs(config):
    """Validate sample directories and required module 2 read subsets before assembly."""
    subset_mode = get_reads_subset_mode(config)
    module_2_root = os.path.join(config.cross_module_params["output_dir"], "module_2")

    for sample_name in config.cross_module_params["samples"]:
        sample_module_2_dir = os.path.join(module_2_root, sample_name)
        if not os.path.exists(sample_module_2_dir):
            raise FileNotFoundError(
                f"Module 2 output directory not found for sample {sample_name}: {sample_module_2_dir}. "
                f"Run module_2_classify_scrub.py before module_3_assembly.py."
            )

        read1, read2, singleton = locate_subset_reads(
            sample_dir=sample_module_2_dir,
            sample_name=sample_name,
            subset_mode=subset_mode,
        )

        if read1 is None or read2 is None:
            raise FileNotFoundError(
                f"Required paired-end read subset for sample {sample_name} not found in {sample_module_2_dir} "
                f"for mode '{subset_mode}'."
            )

        print(
            f"[VALIDATION]: Sample {sample_name} has paired-end read subset '{subset_mode}' "
            f"available in {sample_module_2_dir}."
        )

    print("All module 2 read subsets validated successfully.")


def process_sample(sample_name, config):
    """Run assembly for one sample using the configured subset mode."""
    print(f"\n==========================================")
    print(f" Processing sample: {sample_name}")
    print(f"==========================================")

    subset_mode = get_reads_subset_mode(config)
    sample_module_2_dir = os.path.join(config.cross_module_params["output_dir"], "module_2", sample_name)
    read1, read2, singleton = locate_subset_reads(
        sample_dir=sample_module_2_dir,
        sample_name=sample_name,
        subset_mode=subset_mode,
    )

    if read1 is None or read2 is None:
        raise FileNotFoundError(
            f"Assembly inputs not found for sample {sample_name} in {sample_module_2_dir} "
            f"for subset mode '{subset_mode}'."
        )

    run_rnaviralspades(
        sample_name=sample_name,
        read1=read1,
        read2=read2,
        singleton=singleton,
        config=config,
        subset_mode=subset_mode,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Assemble module 2 viral read subsets with rnaviralspades." 
    )
    parser.add_argument("--config", required=True, help="Path to the config.py file")

    args = parser.parse_args()

    config = load_config(args.config)

    subset_mode = get_reads_subset_mode(config)

    print(f"""
Starting module_3_assembly...

Config file location: {args.config}
Cross-module parameters: {config.cross_module_params}
Module 3 parameters: {config.module_3_params}
Selected subset mode: {subset_mode}
""")

    validate_inputs(config)

    for sample_name in config.cross_module_params["samples"]:
        try:
            process_sample(sample_name, config)
        except Exception as exc:
            print(f"[ERROR]: Failed to assemble sample {sample_name}: {exc}")

    print("\n[module_3_assembly] Finished all samples.")


if __name__ == "__main__":
    main()

import argparse
import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

"""
Module 1: quality control and trimming.

This module is responsible for the first operational step of the pipeline:

1. validate sample inputs from the configured read directory,
2. locate paired-end and singleton read files,
3. run fastp on each sample,
4. write trimmed FASTQ outputs for downstream classification and assembly,
5. aggregate fastp JSON reports for MultiQC.

The script is written in the same operational style as module_2_classify_scrub.py:
- config-driven execution via --config,
- explicit validation of sample directories and read files,
- clear per-sample progress logging,
- subprocess-based command execution with stdout/stderr visible in the terminal.
"""


def load_config(config_path):
    """Load a Python config file as a Python module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def find_sample_read_files(sample_dir, sample_name, read_1_str, read_2_str, singleton_str=None):
    """Locate paired-end and singleton input read files for one sample."""
    read1 = None
    read2 = None
    singleton = None

    for file_name in sorted(os.listdir(sample_dir)):
        full_path = os.path.join(sample_dir, file_name)
        if not os.path.isfile(full_path):
            continue
        if not file_name.endswith((".fastq", ".fastq.gz")):
            continue

        if read_1_str in file_name and "singleton" not in file_name.lower():
            read1 = full_path
        elif read_2_str in file_name and "singleton" not in file_name.lower():
            read2 = full_path
        elif singleton_str is not None and singleton_str in file_name:
            if "_R1" not in file_name and "_R2" not in file_name and "singleton" not in file_name.lower():
                singleton = full_path

    if read1 is None or read2 is None:
        fallback_candidates = sorted(os.listdir(sample_dir))
        for file_name in fallback_candidates:
            full_path = os.path.join(sample_dir, file_name)
            if not os.path.isfile(full_path):
                continue
            if not file_name.endswith((".fastq", ".fastq.gz")):
                continue
            if file_name.startswith(f"{sample_name}_R1"):
                read1 = full_path
            elif file_name.startswith(f"{sample_name}_R2"):
                read2 = full_path
            elif singleton_str is not None and file_name.startswith(f"{sample_name}_singleton"):
                singleton = full_path

    return read1, read2, singleton


def get_thread_count(config):
    """Return the effective thread count for this module from config."""
    module_threads = config.module_1_params.get("threads")
    if module_threads is not None:
        return module_threads
    return config.cross_module_params.get("threads", 4)


def run_fastp(in1, out1, fastp_bin_path, threads, in2=None, out2=None, html_report=None, json_report=None):
    """Run fastp for paired-end or singleton reads."""
    if not os.path.exists(fastp_bin_path):
        raise FileNotFoundError(f"fastp executable not found at {fastp_bin_path}")

    cmd = [
        fastp_bin_path,
        "-i",
        str(in1),
        "-o",
        str(out1),
        "--thread",
        str(threads),
    ]

    if in2 is not None and out2 is not None:
        cmd.extend(["-I", str(in2), "-O", str(out2), "--detect_adapter_for_pe"])

    if html_report is not None:
        cmd.extend(["-h", str(html_report)])
    if json_report is not None:
        cmd.extend(["-j", str(json_report)])

    print(f"[FASTP]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    if not os.path.exists(out1):
        raise FileNotFoundError(f"fastp did not create expected output file: {out1}")
    if in2 is not None and out2 is not None and not os.path.exists(out2):
        raise FileNotFoundError(f"fastp did not create expected output file: {out2}")


def copy_fastp_json_reports(sample_fastp_dir, json_root):
    """Copy fastp JSON reports from each sample directory into the shared MultiQC directory."""
    os.makedirs(json_root, exist_ok=True)
    copied = 0
    for report in sorted(Path(sample_fastp_dir).glob("*.json")):
        shutil.copy2(report, Path(json_root) / report.name)
        copied += 1
    return copied


def run_multiqc(multiqc_bin_path, input_dir, output_dir):
    """Aggregate fastp JSON results into a MultiQC report."""
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        multiqc_bin_path,
        "--outdir",
        output_dir,
        input_dir,
    ]
    print(f"[MULTIQC]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def process_sample(sample_name, config):
    """Run the QC and trimming workflow for a single sample."""
    print(f"\n==========================================")
    print(f" Processing sample: {sample_name}")
    print(f"==========================================")

    sample_input_dir = os.path.join(config.cross_module_params["read_dir"], sample_name)
    if not os.path.exists(sample_input_dir):
        raise FileNotFoundError(
            f"Sample directory for {sample_name} not found: {sample_input_dir}"
        )

    read1, read2, singleton = find_sample_read_files(
        sample_dir=sample_input_dir,
        sample_name=sample_name,
        read_1_str=config.cross_module_params["read_1_str"],
        read_2_str=config.cross_module_params["read_2_str"],
        singleton_str=config.cross_module_params.get("singleton_str"),
    )

    if read1 is None or read2 is None:
        raise FileNotFoundError(
            f"Paired-end read files for sample {sample_name} not found in {sample_input_dir}. "
            f"Expected names containing '{config.cross_module_params['read_1_str']}' and '{config.cross_module_params['read_2_str']}'."
        )

    threads = get_thread_count(config)
    sample_fastp_dir = os.path.join(config.cross_module_params["output_dir"], "module_1", "fastp", sample_name)
    os.makedirs(sample_fastp_dir, exist_ok=True)
    json_root = os.path.join(config.cross_module_params["output_dir"], "module_1", "fastp", "json_files")

    paired_out1 = os.path.join(sample_fastp_dir, f"{sample_name}_R1.fastq.gz")
    paired_out2 = os.path.join(sample_fastp_dir, f"{sample_name}_R2.fastq.gz")
    paired_json = os.path.join(sample_fastp_dir, f"{sample_name}_paired.json")
    paired_html = os.path.join(sample_fastp_dir, f"{sample_name}_paired.html")

    run_fastp(
        in1=read1,
        in2=read2,
        out1=paired_out1,
        out2=paired_out2,
        fastp_bin_path=config.module_1_params["fastp_bin_path"],
        threads=threads,
        html_report=paired_html,
        json_report=paired_json,
    )
    copy_fastp_json_reports(sample_fastp_dir, json_root)
    print(f"[RESULTS]: Wrote paired-end trimmed outputs for {sample_name} to {sample_fastp_dir}")

    if singleton is not None:
        singleton_out = os.path.join(sample_fastp_dir, f"{sample_name}_singleton.fastq.gz")
        singleton_json = os.path.join(sample_fastp_dir, f"{sample_name}_singleton.json")
        singleton_html = os.path.join(sample_fastp_dir, f"{sample_name}_singleton.html")
        run_fastp(
            in1=singleton,
            out1=singleton_out,
            fastp_bin_path=config.module_1_params["fastp_bin_path"],
            threads=threads,
            html_report=singleton_html,
            json_report=singleton_json,
        )
        copy_fastp_json_reports(sample_fastp_dir, json_root)
        print(f"[RESULTS]: Wrote singleton trimmed outputs for {sample_name} to {sample_fastp_dir}")


def validate_inputs(config):
    """Validate configured sample directories and required read files before running fastp."""
    for sample_name in config.cross_module_params["samples"]:
        sample_dir = os.path.join(config.cross_module_params["read_dir"], sample_name)
        if not os.path.exists(sample_dir):
            raise FileNotFoundError(
                f"Sample directory not found for {sample_name}: {sample_dir}"
            )

        read1, read2, singleton = find_sample_read_files(
            sample_dir=sample_dir,
            sample_name=sample_name,
            read_1_str=config.cross_module_params["read_1_str"],
            read_2_str=config.cross_module_params["read_2_str"],
            singleton_str=config.cross_module_params.get("singleton_str"),
        )

        if read1 is None or read2 is None:
            raise FileNotFoundError(
                f"Paired-end reads for sample {sample_name} not found in {sample_dir}. "
                f"Looked for filenames containing '{config.cross_module_params['read_1_str']}' and '{config.cross_module_params['read_2_str']}'."
            )

        print(f"[VALIDATION]: Sample {sample_name} has paired-end reads and optional singleton reads located successfully.")

    print("All sample inputs validated successfully.")


def main():
    parser = argparse.ArgumentParser(
        description="QC and trimming step for paired-end and singleton sequencing reads."
    )
    parser.add_argument("--config", required=True, help="Path to the config.py file")
    args = parser.parse_args()

    config = load_config(args.config)

    print(f"""
Starting module_1_qc_trim...

Config file location: {args.config}
Cross-module parameters: {config.cross_module_params}
Module 1 parameters: {config.module_1_params}
""")

    validate_inputs(config)

    fastp_root = os.path.join(config.cross_module_params["output_dir"], "module_1", "fastp")
    os.makedirs(fastp_root, exist_ok=True)
    json_root = os.path.join(fastp_root, "json_files")
    os.makedirs(json_root, exist_ok=True)

    for sample_name in config.cross_module_params["samples"]:
        try:
            process_sample(sample_name, config)
        except Exception as exc:
            print(f"[ERROR]: Failed to process sample {sample_name}: {exc}")

    multiqc_dir = os.path.join(config.cross_module_params["output_dir"], "module_1", "multiqc")
    json_files = sorted(Path(json_root).glob("*.json"))
    if json_files:
        run_multiqc(
            multiqc_bin_path=config.module_1_params["multi_qc_bin_path"],
            input_dir=json_root,
            output_dir=multiqc_dir,
        )
    else:
        print(f"[WARN]: No fastp JSON reports were found in {json_root}; MultiQC will not run.")

    print("\n[module_1_qc_trim] Finished all samples successfully.")


if __name__ == "__main__":
    main()


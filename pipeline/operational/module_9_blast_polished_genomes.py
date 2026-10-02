import argparse
import csv
import importlib.util
import os
import re
import subprocess


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
    value = get_config_value(config, "module_9_params", key, None)
    if value is None:
        raise FileNotFoundError(
            f"Executable '{key}' is not configured. Set an explicit path or command in module_9_params."
        )
    if os.path.exists(value):
        return value
    resolved = subprocess.run(["which", value], capture_output=True, text=True, check=False)
    if resolved.returncode == 0 and resolved.stdout.strip():
        return resolved.stdout.strip().splitlines()[0]
    raise FileNotFoundError(
        f"Configured executable '{value}' for '{key}' was not found. Set an explicit path in module_9_params."
    )


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def get_module_6_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_6", sample_name)


def get_module_9_dir(config, sample_name=None):
    root = os.path.join(config.cross_module_params["output_dir"], "module_9")
    if sample_name is None:
        return root
    return os.path.join(root, sample_name)


def get_blast_db_path(config):
    blast_db_path = get_config_value(config, "module_9_params", "blast_db_path", None)
    if blast_db_path is None:
        raise FileNotFoundError("module_9_params['blast_db_path'] must be set to a valid BLAST database prefix.")
    return blast_db_path


def get_blastn_bin_path(config):
    return resolve_executable_path(config, "blastn_bin_path")


def parse_taxids(raw_value):
    if raw_value is None or raw_value == "":
        return []
    if isinstance(raw_value, (int, float)):
        return [int(raw_value)]
    values = []
    for chunk in re.split(r"[;,|]", str(raw_value)):
        chunk = chunk.strip()
        if not chunk:
            continue
        match = re.search(r"(\d+)", chunk)
        if match:
            values.append(int(match.group(1)))
    return values


def extract_blast_taxonomy(staxids, sscinames):
    taxids = parse_taxids(staxids)
    subject_taxid = taxids[0] if taxids else None

    species_sci_name = ""
    if sscinames:
        species_sci_name = str(sscinames).split(";")[0].strip()
        if not species_sci_name:
            species_sci_name = str(sscinames).split("|")[0].strip()

    return {
        "subject_accession": "",
        "subject_taxid": subject_taxid,
        "subject_name": species_sci_name,
        "species_taxid": subject_taxid,
        "species_sci_name": species_sci_name,
    }


def find_polished_genome(config, sample_name):
    """Prefer the final combined polished genome written by module 6."""
    sample_root = get_module_6_dir(config, sample_name)
    if not os.path.isdir(sample_root):
        raise FileNotFoundError(f"Module 6 output directory not found for {sample_name}: {sample_root}")

    preferred_names = [
        f"{sample_name}_consensus_combined.fasta",
        f"{sample_name}_combined_consensus.fasta",
        f"{sample_name}_final_polished_consensus.fa",
        f"{sample_name}_polished_consensus.fasta",
    ]
    for filename in preferred_names:
        candidate = os.path.join(sample_root, filename)
        if os.path.exists(candidate):
            return candidate

    matches = []
    for root, _, files in os.walk(sample_root):
        for filename in files:
            lower_name = filename.lower()
            if lower_name.endswith((".fa", ".fasta")) and (
                "consensus" in lower_name or "polished" in lower_name or "genome" in lower_name
            ):
                matches.append(os.path.join(root, filename))

    if matches:
        return sorted(matches)[0]

    raise FileNotFoundError(
        f"No polished genome FASTA found for sample {sample_name} under {sample_root}. "
        "Expected output from module 6 such as <sample>_consensus_combined.fasta."
    )


# ---------------------------------------------------------------------------
# Database validation
# ---------------------------------------------------------------------------

def validate_blast_db(blast_db_path):
    if not blast_db_path:
        raise FileNotFoundError("No BLAST database path was provided for module 9.")

    db_suffixes = [
        ".nhr", ".nin", ".nsq", ".ntf", ".nto",
        ".nal", ".nog", ".nos", ".not", ".ndb",
    ]
    if os.path.isdir(blast_db_path):
        files = os.listdir(blast_db_path)
        if files:
            return
    if any(os.path.exists(f"{blast_db_path}{suffix}") for suffix in db_suffixes):
        return
    raise FileNotFoundError(f"BLAST database not found at {blast_db_path}. Set module_9_params['blast_db_path'] to the database prefix.")


# ---------------------------------------------------------------------------
# BLAST execution
# ---------------------------------------------------------------------------

def run_blastn_for_sample(config, sample_name, query_fasta, output_tsv):
    blastn_bin = get_blastn_bin_path(config)
    blast_db_path = get_blast_db_path(config)
    validate_blast_db(blast_db_path)

    evalue = float(get_config_value(config, "module_9_params", "evalue", 1e-6))
    max_target_seqs = int(get_config_value(config, "module_9_params", "max_target_seqs", 20))
    threads = int(get_config_value(config, "module_9_params", "threads", 8))
    outfmt = get_config_value(
        config,
        "module_9_params",
        "outfmt",
        "6 qseqid sseqid pident length qlen slen bitscore evalue staxids sscinames stitle",
    )

    command = [
        blastn_bin,
        "-task",
        "megablast",
        "-num_threads",
        str(threads),
        "-query",
        query_fasta,
        "-db",
        blast_db_path,
        "-out",
        output_tsv,
        "-outfmt",
        outfmt,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
    ]

    print(f"[BLASTN]: Running module 9 BLAST for {sample_name} -> {blast_db_path}")
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return output_tsv


# ---------------------------------------------------------------------------
# Summary writing
# ---------------------------------------------------------------------------

def parse_blast_report(config, blast_report_path):
    if not os.path.exists(blast_report_path) or os.path.getsize(blast_report_path) == 0:
        return []

    rows = []
    with open(blast_report_path, "r") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < 11:
                continue
            qseqid, sseqid, pident, length, qlen, slen, bitscore, evalue, staxids, sscinames, stitle = fields[:11]
            taxonomy = extract_blast_taxonomy(staxids, sscinames)
            rows.append(
                {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "pident": float(pident),
                    "length": int(length),
                    "qlen": int(qlen),
                    "slen": int(slen),
                    "bitscore": float(bitscore),
                    "evalue": float(evalue),
                    "staxids": staxids,
                    "sscinames": sscinames,
                    "stitle": stitle,
                    "subject_accession": sseqid,
                    "subject_taxid": taxonomy["subject_taxid"],
                    "subject_name": taxonomy["subject_name"],
                    "species_taxid": taxonomy["species_taxid"],
                    "species_sci_name": taxonomy["species_sci_name"],
                }
            )
    return rows


def write_blast_summary(config, rows):
    module_9_root = get_module_9_dir(config)
    os.makedirs(module_9_root, exist_ok=True)
    summary_path = os.path.join(module_9_root, "Identified_blast_hits.tsv")

    fieldnames = [
        "sample_name",
        "query_fasta",
        "qseqid",
        "sseqid",
        "subject_accession",
        "subject_taxid",
        "subject_name",
        "species_taxid",
        "species_sci_name",
        "pident",
        "alignment_length",
        "query_length",
        "subject_length",
        "bitscore",
        "evalue",
        "staxids",
        "sscinames",
        "stitle",
    ]

    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    print(f"[SUMMARY]: Wrote module 9 BLAST summary to {summary_path}")
    return summary_path


def process_sample(config, sample_name):
    query_fasta = find_polished_genome(config, sample_name)
    sample_out_dir = get_module_9_dir(config, sample_name)
    os.makedirs(sample_out_dir, exist_ok=True)

    blast_tsv = os.path.join(sample_out_dir, f"{sample_name}_module_9_blastn.tsv")
    run_blastn_for_sample(config, sample_name, query_fasta, blast_tsv)

    hits = parse_blast_report(config, blast_tsv)
    summarized_rows = []
    for hit in hits:
        summarized_rows.append(
            {
                "sample_name": sample_name,
                "query_fasta": query_fasta,
                "qseqid": hit["qseqid"],
                "sseqid": hit["sseqid"],
                "subject_accession": hit.get("subject_accession", hit["sseqid"]),
                "subject_taxid": hit.get("subject_taxid"),
                "subject_name": hit.get("subject_name"),
                "species_taxid": hit.get("species_taxid"),
                "species_sci_name": hit.get("species_sci_name"),
                "pident": hit["pident"],
                "alignment_length": hit["length"],
                "query_length": hit["qlen"],
                "subject_length": hit["slen"],
                "bitscore": hit["bitscore"],
                "evalue": hit["evalue"],
                "staxids": hit["staxids"],
                "sscinames": hit["sscinames"],
                "stitle": hit["stitle"],
            }
        )

    if not summarized_rows:
        summarized_rows.append(
            {
                "sample_name": sample_name,
                "query_fasta": query_fasta,
                "qseqid": "",
                "sseqid": "",
                "subject_accession": "",
                "subject_taxid": "",
                "subject_name": "",
                "species_taxid": "",
                "species_sci_name": "",
                "pident": "",
                "alignment_length": "",
                "query_length": "",
                "subject_length": "",
                "bitscore": "",
                "evalue": "",
                "staxids": "",
                "sscinames": "",
                "stitle": "No significant BLAST hits",
            }
        )

    print(f"[RESULT]: {sample_name} produced {len(summarized_rows)} BLAST record(s) in {blast_tsv}")
    return summarized_rows


def main():
    parser = argparse.ArgumentParser(description="BLAST polished genomes from module 6 against a user-specified reference database.")
    parser.add_argument("--config", required=True, help="Path to the pipeline config Python file.")
    args = parser.parse_args()

    config = load_config(args.config)

    sample_list = get_config_value(config, "module_9_params", "input_samples", None)
    if sample_list is None:
        sample_list = config.cross_module_params["samples"]

    all_rows = []
    print(f"\nStarting module_9_blast_polished_genomes...\nConfig file: {args.config}\n")

    for sample_name in sample_list:
        try:
            all_rows.extend(process_sample(config, sample_name))
        except Exception as exc:
            print(f"[ERROR]: Failed for sample {sample_name}: {exc}")

    write_blast_summary(config, all_rows)
    print("\n[module_9_blast_polished_genomes] Finished all samples.")


if __name__ == "__main__":
    main()

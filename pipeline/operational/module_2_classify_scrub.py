import argparse
import gzip
import importlib.util
import os
import subprocess
from collections import defaultdict
from pathlib import Path

"""
Module 2: Kraken2 taxonomic classification and read-level separation.

This module is designed to run after the QC/trim step in module 1. It assumes
that each sample has been trimmed and stored in a sample-specific output
subdirectory, typically under:

    <output_dir>/module_1/fastp/<sample_name>/

For each sample, the script:

1. finds paired-end and singleton trimmed FASTQ files,
2. runs Kraken2 classification with k2 classify,
3. parses the classification output,
4. identifies reads assigned to the configured viral taxonomic scope,
5. splits reads into two non-overlapping groups:
      a) non-viral reads
      b) viral + unclassified reads
6. writes out the filtered read sets for downstream use.

The design matches the operational style of module_1_qc_trim.py, including:
- argparse-driven config loading,
- explicit validation of sample and file existence,
- standard printed progress logs,
- per-sample processing,
- command execution through subprocess.run(..., check=True).
"""


def load_config(config_path):
    """Load a Python config file as a module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def normalize_read_id(read_id):
    """Normalize a read identifier so paired-end mate suffixes /1 and /2 match."""
    if read_id is None:
        return None
    return read_id.strip().rstrip("/1").rstrip("/2")


def open_fq(filepath, mode):
    """Open FASTQ files as text, including gzipped inputs."""
    return gzip.open(filepath, mode) if str(filepath).endswith(".gz") else open(filepath, mode)


def parse_kraken_report(report_path):
    """Parse a Kraken2 report file into a list of taxonomic entries.

    Each entry contains the taxon metadata used to determine whether a read
    assignment falls within the configured viral scope.
    """
    entries = []
    with open(report_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            name = parts[5]
            depth = len(name) - len(name.lstrip())
            entries.append(
                {
                    "percent": parts[0].strip(),
                    "reads_covered": int(parts[1].strip()),
                    "reads_assigned": int(parts[2].strip()),
                    "rank": parts[3].strip(),
                    "taxid": parts[4].strip(),
                    "name": name.strip(),
                    "depth": depth,
                }
            )
    return entries


def read_taxonomy_nodes(nodes_path):
    """Load a taxonomy tree from NCBI nodes.dmp into parent->children mappings."""
    if not nodes_path or not os.path.exists(nodes_path):
        return {}, {}

    children_by_parent = defaultdict(set)
    parent_by_child = {}
    with open(nodes_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or "\t|" not in line:
                continue
            fields = [part.strip() for part in line.split("\t|")]
            if len(fields) < 2:
                continue
            try:
                taxid = int(fields[0])
                parent_taxid = int(fields[1])
            except ValueError:
                continue
            parent_by_child[taxid] = parent_taxid
            children_by_parent[parent_taxid].add(taxid)
    return parent_by_child, children_by_parent


def collect_taxonomy_descendants(root_taxid, parent_by_child, children_by_parent):
    """Return the root taxid plus every descendant taxid beneath it."""
    if root_taxid is None:
        return set()

    descendants = {str(root_taxid)}
    queue = [int(root_taxid)]
    seen = set()

    while queue:
        current = queue.pop(0)
        if current in seen:
            continue
        seen.add(current)
        for child in sorted(children_by_parent.get(current, set())):
            child_str = str(child)
            if child_str not in descendants:
                descendants.add(child_str)
                queue.append(child)

    return descendants


def get_taxids_in_scope(report_paths, viral_scope_taxid, nodes_path=None):
    """Return the exact descendant taxid set for the configured scope.

    Prefer the NCBI taxonomy dump when available; this avoids relying on the
    report indentation depth as a proxy for true taxonomic ancestry.
    """
    scope_taxids = set()

    if nodes_path and os.path.exists(nodes_path):
        parent_by_child, children_by_parent = read_taxonomy_nodes(nodes_path)
        scope_taxids = collect_taxonomy_descendants(viral_scope_taxid, parent_by_child, children_by_parent)
        if scope_taxids:
            return scope_taxids

    for report_path in report_paths:
        if not os.path.exists(report_path):
            continue

        entries = parse_kraken_report(report_path)
        viral_idx = next(
            (i for i, entry in enumerate(entries) if entry["taxid"] == str(viral_scope_taxid)),
            None,
        )
        if viral_idx is None:
            continue

        viral_depth = entries[viral_idx]["depth"]
        scope_taxids.add(str(viral_scope_taxid))

        for entry in entries:
            if entry["depth"] > viral_depth and entry["taxid"] not in ("", "0"):
                scope_taxids.add(str(entry["taxid"]))

    return scope_taxids


def get_taxonomy_nodes_path(config):
    """Resolve the taxonomy nodes dump path from module config if available."""
    if hasattr(config, "module_2_params") and config.module_2_params.get("nodes_path"):
        return config.module_2_params.get("nodes_path")
    if hasattr(config, "module_4_params") and config.module_4_params.get("nodes_path"):
        return config.module_4_params.get("nodes_path")
    return None


def parse_kraken_classification(classification_path):
    """Parse a Kraken2 classification output file into read -> {status, taxid}."""
    assignments = {}
    with open(classification_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            status, read_id, taxid = parts[0], parts[1], parts[2]
            assignments[normalize_read_id(read_id)] = {
                "status": status,
                "taxid": str(taxid).strip(),
            }
    return assignments


def split_reads_by_viral_status(classification_assignments, viral_taxids, include_unclassified=True):
    """Split read IDs into viral + optional-unclassified vs non-viral reads.

    The broad bucket used downstream should include all reads in the configured viral
    scope plus any unclassified reads. This is the set used to make the filtered
    FASTQ outputs, and it should be a superset of the strict viral-only bucket.
    """
    viral_reads = set()
    non_viral = set()

    for read_id, info in classification_assignments.items():
        if read_id is None:
            continue
        status = info["status"]
        taxid = info["taxid"]

        if status == "U":
            if include_unclassified:
                viral_reads.add(read_id)
            else:
                non_viral.add(read_id)
        elif taxid in viral_taxids:
            viral_reads.add(read_id)
        else:
            non_viral.add(read_id)

    return non_viral, viral_reads


def split_reads_by_viral_only(classification_assignments, viral_taxids):
    """Split read IDs into viral-only vs non-viral-or-unclassified groups."""
    viral_only = set()
    non_viral_or_unclassified = set()

    for read_id, info in classification_assignments.items():
        if read_id is None:
            continue
        status = info["status"]
        taxid = info["taxid"]

        if status != "U" and taxid in viral_taxids:
            viral_only.add(read_id)
        else:
            non_viral_or_unclassified.add(read_id)

    return non_viral_or_unclassified, viral_only


def find_dominant_viral_species(report_paths, viral_scope_taxid):
    """Identify the dominant viral species from Kraken2 report files.

    This follows the same logic used in the original basic pipeline: find the
    virus root taxid, sum species-level reads beneath it, and select the species
    with the largest lineage read count.
    """
    species_counts = defaultdict(int)
    species_names = {}

    for report_path in report_paths:
        if not os.path.exists(report_path):
            continue
        entries = parse_kraken_report(report_path)
        viral_idx = next(
            (i for i, entry in enumerate(entries) if entry["taxid"] == str(viral_scope_taxid)),
            None,
        )
        if viral_idx is None:
            continue

        viral_depth = entries[viral_idx]["depth"]
        i = viral_idx + 1
        while i < len(entries) and entries[i]["depth"] > viral_depth:
            entry = entries[i]
            if entry["rank"] == "S" and entry["reads_covered"] > 0:
                species_counts[str(entry["taxid"])] += entry["reads_covered"]
                species_names[str(entry["taxid"])] = entry["name"]
            i += 1

    if not species_counts:
        return None

    dominant_taxid = max(species_counts, key=lambda taxid: (species_counts[taxid], taxid))
    return {
        "taxid": dominant_taxid,
        "name": species_names.get(dominant_taxid, "Unknown"),
        "reads": species_counts[dominant_taxid],
    }


def get_allowed_taxids(report_path, viral_scope_taxid, dominant_taxid):
    """Return the taxids belonging to the dominant viral lineage and its descendants."""
    entries = parse_kraken_report(report_path)
    viral_idx = next(
        (i for i, entry in enumerate(entries) if entry["taxid"] == str(viral_scope_taxid)),
        None,
    )
    dom_idx = next(
        (i for i, entry in enumerate(entries) if entry["taxid"] == str(dominant_taxid)),
        None,
    )
    if viral_idx is None or dom_idx is None:
        return set()

    viral_depth = entries[viral_idx]["depth"]
    dom_depth = entries[dom_idx]["depth"]
    allowed = set()

    current_depth = dom_depth
    for i in range(dom_idx, viral_idx - 1, -1):
        entry = entries[i]
        if entry["depth"] < current_depth:
            if entry["depth"] > viral_depth:
                allowed.add(str(entry["taxid"]))
            current_depth = entry["depth"]
            if current_depth <= viral_depth:
                break

    i = dom_idx + 1
    while i < len(entries) and entries[i]["depth"] > dom_depth:
        allowed.add(str(entries[i]["taxid"]))
        i += 1

    return allowed


def summarize_classification_assignments(assignments, viral_taxids, dominant_taxids, dominant_species_name):
    """Summarize unique read assignments into the requested taxonomic categories."""
    total_reads = len(assignments)
    human = 0
    viral = 0
    dominant_species = 0
    unclassified = 0
    non_viral_non_human = 0

    for info in assignments.values():
        taxid = str(info.get("taxid", "")).strip()
        status = info.get("status", "")

        if status == "U":
            unclassified += 1
        elif taxid == "9606":
            human += 1
        else:
            if taxid in viral_taxids:
                viral += 1
                if taxid in dominant_taxids or dominant_taxids is None:
                    dominant_species += 1
            else:
                non_viral_non_human += 1

    return {
        "total_reads": total_reads,
        "Human": human,
        "Viral": viral,
        "Dominant viral species": dominant_species,
        "Non-viral / Non-human": non_viral_non_human,
        "Unclassified": unclassified,
        "dominant_species_name": dominant_species_name,
    }


def dominant_virus_plus_unclassified_ids(assignments, dominant_taxids, viral_taxids=None):
    """Return IDs for unclassified reads plus reads in the dominant viral lineage.

    The dominant lineage is treated as a subset of the configured viral scope, and
    the returned set is therefore constrained to that exact viral subtree to avoid
    unbounded expansion beyond the intended scope.
    """
    if viral_taxids is not None:
        dominant_taxids = set(dominant_taxids) & set(viral_taxids)

    keep = set()
    for read_id, info in assignments.items():
        if read_id is None:
            continue
        if info.get("status") == "U":
            keep.add(read_id)
        elif str(info.get("taxid", "")).strip() in dominant_taxids:
            keep.add(read_id)
    return keep


def format_count(value, total):
    """Format count as integer + percentage of the total for this read type."""
    if total == 0:
        return "0 (0.0%)"
    percent = (float(value) / float(total)) * 100.0
    return f"{value} ({percent:.1f}%)"


def write_summary_report(report_rows, output_path):
    """Write a TSV summary report of paired and singleton read counts by taxonomic category."""
    with open(output_path, "w") as handle:
        handle.write(
            "Sample\tReadType\tIdentity of dominant virus species\tDominant virus species reads\tAll viral reads\tHuman reads\tNon-human + Non-viral reads\tUnclassified reads\n"
        )
        for row in report_rows:
            handle.write(
                f"{row['sample']}\t{row['read_type']}\t"
                f"{row['dominant_species_identity']}\t{row['Dominant viral species']}\t"
                f"{row['Viral']}\t{row['Human']}\t{row['Non-viral / Non-human']}\t{row['Unclassified']}\n"
            )


def write_filtered_fastq(in_fastq, out_fastq, keep_ids):
    """Write a FASTQ file containing only reads whose normalized IDs are in keep_ids."""
    count = 0
    with open_fq(in_fastq, "rt") as infile, open(out_fastq, "w") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break
            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()
            if not seq:
                break

            read_id = normalize_read_id(header[1:].split()[0])
            if read_id in keep_ids:
                outfile.write(f"{header}{seq}{plus}{qual}")
                count += 1
    return count


def write_filtered_paired_fastq(r1_in, r2_in, r1_out, r2_out, keep_ids):
    """Filter a paired-end FASTQ pair using a shared set of retained read IDs."""
    count = 0
    with (
        open_fq(r1_in, "rt") as in1,
        open_fq(r2_in, "rt") as in2,
        open(r1_out, "w") as out1,
        open(r2_out, "w") as out2,
    ):
        while True:
            h1 = in1.readline()
            h2 = in2.readline()
            if not h1 or not h2:
                break
            s1 = in1.readline(); s2 = in2.readline()
            p1 = in1.readline(); p2 = in2.readline()
            q1 = in1.readline(); q2 = in2.readline()
            if not s1 or not s2:
                break

            read_id_1 = normalize_read_id(h1[1:].split()[0])
            read_id_2 = normalize_read_id(h2[1:].split()[0])

            if read_id_1 in keep_ids and read_id_2 in keep_ids:
                out1.write(f"{h1}{s1}{p1}{q1}")
                out2.write(f"{h2}{s2}{p2}{q2}")
                count += 1
    return count


def run_kraken2_classify(input_reads, output_file, report_file, db_path, k2_bin_path, threads, paired=False):
    """Run Kraken2 classification with k2 classify on one sample read set."""
    cmd = [
        k2_bin_path,
        "classify",
        "--db",
        db_path,
        "--threads",
        str(threads),
        "--output",
        str(output_file),
        "--report",
        str(report_file),
        "--use-daemon"
    ]

    if paired:
        cmd.extend(["--paired", str(input_reads[0]), str(input_reads[1])])
    else:
        cmd.append(str(input_reads[0]))

    print(f"[KRAKEN2]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def find_trimmed_reads(sample_dir, sample_name, read_1_str, read_2_str, singleton_str):
    """Locate paired-end and singleton trimmed FASTQ files for a sample."""
    read1 = None
    read2 = None
    singleton = None

    for file_name in sorted(os.listdir(sample_dir)):
        full_path = os.path.join(sample_dir, file_name)
        if not os.path.isfile(full_path):
            continue

        if read_1_str in file_name and file_name.endswith((".fastq", ".fastq.gz")):
            read1 = full_path
        elif read_2_str in file_name and file_name.endswith((".fastq", ".fastq.gz")):
            read2 = full_path
        elif singleton_str is not None and singleton_str in file_name and file_name.endswith((".fastq", ".fastq.gz")):
            if "_R1" not in file_name and "_R2" not in file_name and "_singleton" not in file_name:
                singleton = full_path

    if read1 is None or read2 is None:
        # allow more explicit module_1 naming fallbacks
        for file_name in sorted(os.listdir(sample_dir)):
            full_path = os.path.join(sample_dir, file_name)
            if not os.path.isfile(full_path):
                continue
            if file_name.startswith(f"{sample_name}_R1") and file_name.endswith((".fastq", ".fastq.gz")):
                read1 = full_path
            elif file_name.startswith(f"{sample_name}_R2") and file_name.endswith((".fastq", ".fastq.gz")):
                read2 = full_path
            elif file_name.startswith(f"{sample_name}_singleton") and file_name.endswith((".fastq", ".fastq.gz")):
                singleton = full_path

    return read1, read2, singleton


def process_sample(sample_name, config):
    """Process one sample through Kraken2 classification and taxon-based read separation."""
    print(f"\n==========================================")
    print(f" Processing sample: {sample_name}")
    print(f"==========================================")

    module_1_dir = os.path.join(
        config.cross_module_params["output_dir"],
        "module_1",
        "fastp",
        sample_name,
    )
    if not os.path.exists(module_1_dir):
        raise FileNotFoundError(
            f"Trimmed output directory for sample {sample_name} not found: {module_1_dir}"
        )

    read1, read2, singleton = find_trimmed_reads(
        sample_dir=module_1_dir,
        sample_name=sample_name,
        read_1_str=config.cross_module_params["read_1_str"],
        read_2_str=config.cross_module_params["read_2_str"],
        singleton_str=config.cross_module_params.get("singleton_str"),
    )

    if read1 is None or read2 is None:
        raise FileNotFoundError(
            f"Paired-end trimmed read files not found for sample {sample_name} in {module_1_dir}. "
            f"Expected files containing '{config.cross_module_params['read_1_str']}' and '{config.cross_module_params['read_2_str']}'."
        )

    sample_out_dir = os.path.join(
        config.cross_module_params["output_dir"],
        "module_2",
        sample_name,
    )
    os.makedirs(sample_out_dir, exist_ok=True)

    paired_classification_out = os.path.join(sample_out_dir, f"{sample_name}_paired.kraken")
    paired_report_out = os.path.join(sample_out_dir, f"{sample_name}_paired.report")
    run_kraken2_classify(
        input_reads=[read1, read2],
        output_file=paired_classification_out,
        report_file=paired_report_out,
        db_path=config.module_2_params["kraken2_db_path"],
        k2_bin_path=config.module_2_params["Kraken2_k2_bin_path"],
        threads=config.module_2_params.get(
            "threads",
            config.cross_module_params.get("threads", 4),
        ),
        paired=True,
    )

    paired_assignments = parse_kraken_classification(paired_classification_out)
    nodes_path = get_taxonomy_nodes_path(config)
    paired_viral_taxids = get_taxids_in_scope(
        [paired_report_out],
        config.module_2_params["scope_to_keep_taxId"],
        nodes_path,
    )
    non_viral_paired, viral_plus_unclassified_paired = split_reads_by_viral_status(
        paired_assignments,
        paired_viral_taxids,
        include_unclassified=True,
    )
    _, viral_only_paired = split_reads_by_viral_only(paired_assignments, paired_viral_taxids)

    non_viral_r1 = os.path.join(sample_out_dir, f"{sample_name}_non_viral_R1.fastq")
    non_viral_r2 = os.path.join(sample_out_dir, f"{sample_name}_non_viral_R2.fastq")
    viral_r1 = os.path.join(sample_out_dir, f"{sample_name}_viral_plus_unclassified_R1.fastq")
    viral_r2 = os.path.join(sample_out_dir, f"{sample_name}_viral_plus_unclassified_R2.fastq")
    viral_only_r1 = os.path.join(sample_out_dir, f"{sample_name}_viral_R1.fastq")
    viral_only_r2 = os.path.join(sample_out_dir, f"{sample_name}_viral_R2.fastq")

    paired_non_viral_count = write_filtered_paired_fastq(
        r1_in=read1,
        r2_in=read2,
        r1_out=non_viral_r1,
        r2_out=non_viral_r2,
        keep_ids=non_viral_paired,
    )
    paired_viral_count = write_filtered_paired_fastq(
        r1_in=read1,
        r2_in=read2,
        r1_out=viral_r1,
        r2_out=viral_r2,
        keep_ids=viral_plus_unclassified_paired,
    )
    paired_viral_only_count = write_filtered_paired_fastq(
        r1_in=read1,
        r2_in=read2,
        r1_out=viral_only_r1,
        r2_out=viral_only_r2,
        keep_ids=viral_only_paired,
    )

    print(
        f"[PAIRED GROUPS]: non_viral={paired_non_viral_count} pairs; "
        f"viral_plus_unclassified={paired_viral_count} pairs; "
        f"viral_only={paired_viral_only_count} pairs"
    )

    singleton_report_out = None
    singleton_assignments = {}
    singleton_taxids = set()
    singleton_summary = None
    if singleton is not None:
        singleton_classification_out = os.path.join(sample_out_dir, f"{sample_name}_singleton.kraken")
        singleton_report_out = os.path.join(sample_out_dir, f"{sample_name}_singleton.report")
        run_kraken2_classify(
            input_reads=[singleton],
            output_file=singleton_classification_out,
            report_file=singleton_report_out,
            db_path=config.module_2_params["kraken2_db_path"],
            k2_bin_path=config.module_2_params["Kraken2_k2_bin_path"],
            threads=config.module_2_params.get(
                "threads",
                config.cross_module_params.get("threads", 4),
            ),
            paired=False,
        )

        singleton_assignments = parse_kraken_classification(singleton_classification_out)
        singleton_taxids = get_taxids_in_scope(
            [singleton_report_out],
            config.module_2_params["scope_to_keep_taxId"],
            nodes_path,
        )
        non_viral_singleton, viral_plus_unclassified_singleton = split_reads_by_viral_status(
            singleton_assignments,
            singleton_taxids,
            include_unclassified=True,
        )
        _, viral_only_singleton = split_reads_by_viral_only(singleton_assignments, singleton_taxids)

        non_viral_singleton_out = os.path.join(sample_out_dir, f"{sample_name}_non_viral_singleton.fastq")
        viral_singleton_out = os.path.join(sample_out_dir, f"{sample_name}_viral_plus_unclassified_singleton.fastq")
        viral_only_singleton_out = os.path.join(sample_out_dir, f"{sample_name}_viral_singleton.fastq")

        singleton_non_viral_count = write_filtered_fastq(
            singleton,
            non_viral_singleton_out,
            non_viral_singleton,
        )
        singleton_viral_count = write_filtered_fastq(
            singleton,
            viral_singleton_out,
            viral_plus_unclassified_singleton,
        )
        singleton_viral_only_count = write_filtered_fastq(
            singleton,
            viral_only_singleton_out,
            viral_only_singleton,
        )

        print(
            f"[SINGLETON GROUPS]: non_viral={singleton_non_viral_count} reads; "
            f"viral_plus_unclassified={singleton_viral_count} reads; "
            f"viral_only={singleton_viral_only_count} reads"
        )

    report_paths = [paired_report_out]
    if singleton_report_out is not None:
        report_paths.append(singleton_report_out)
    dominant_species = find_dominant_viral_species(report_paths, config.module_2_params["scope_to_keep_taxId"])
    dominant_species_name = dominant_species["name"] if dominant_species else "None"
    dominant_taxid = dominant_species["taxid"] if dominant_species else None
    dominant_taxids = set()
    if dominant_taxid is not None:
        dominant_taxids = get_allowed_taxids(report_paths[0], config.module_2_params["scope_to_keep_taxId"], dominant_taxid)

    dominant_paired_ids = (
        dominant_virus_plus_unclassified_ids(paired_assignments, dominant_taxids, paired_viral_taxids)
        if dominant_taxids
        else set()
    )
    dominant_paired_out_r1 = os.path.join(sample_out_dir, f"{sample_name}_dominant_virus_plus_unclassified_R1.fastq")
    dominant_paired_out_r2 = os.path.join(sample_out_dir, f"{sample_name}_dominant_virus_plus_unclassified_R2.fastq")
    dominant_paired_count = write_filtered_paired_fastq(
        r1_in=read1,
        r2_in=read2,
        r1_out=dominant_paired_out_r1,
        r2_out=dominant_paired_out_r2,
        keep_ids=dominant_paired_ids,
    )
    print(f"[PAIRED DOMINANT VIRUS+UNCLASSIFIED]: {dominant_paired_count} pairs")

    dominant_singleton_ids = set()
    if singleton is not None and dominant_taxids:
        dominant_singleton_ids = dominant_virus_plus_unclassified_ids(
            singleton_assignments,
            dominant_taxids,
            singleton_taxids,
        )
    dominant_singleton_out = os.path.join(sample_out_dir, f"{sample_name}_dominant_virus_plus_unclassified_singleton.fastq")
    dominant_singleton_count = write_filtered_fastq(
        singleton,
        dominant_singleton_out,
        dominant_singleton_ids,
    ) if singleton is not None else 0
    if singleton is not None:
        print(f"[SINGLETON DOMINANT VIRUS+UNCLASSIFIED]: {dominant_singleton_count} reads")

    paired_summary = summarize_classification_assignments(
        paired_assignments,
        paired_viral_taxids,
        dominant_taxids,
        dominant_species_name,
    )
    paired_row = {
        "sample": sample_name,
        "read_type": "paired-end",
        "Human": format_count(paired_summary["Human"], paired_summary["total_reads"]),
        "Viral": format_count(paired_summary["Viral"], paired_summary["total_reads"]),
        "Dominant viral species": format_count(paired_summary["Dominant viral species"], paired_summary["total_reads"]),
        "Non-viral / Non-human": format_count(paired_summary["Non-viral / Non-human"], paired_summary["total_reads"]),
        "Unclassified": format_count(paired_summary["Unclassified"], paired_summary["total_reads"]),
        "dominant_species_identity": f"{dominant_species_name} (taxid {dominant_taxid})" if dominant_species else "No viral species detected",
    }

    singleton_row = None
    if singleton is not None:
        singleton_summary = summarize_classification_assignments(
            singleton_assignments,
            singleton_taxids,
            dominant_taxids,
            dominant_species_name,
        )
        singleton_row = {
            "sample": sample_name,
            "read_type": "singleton",
            "Human": format_count(singleton_summary["Human"], singleton_summary["total_reads"]),
            "Viral": format_count(singleton_summary["Viral"], singleton_summary["total_reads"]),
            "Dominant viral species": format_count(singleton_summary["Dominant viral species"], singleton_summary["total_reads"]),
            "Non-viral / Non-human": format_count(singleton_summary["Non-viral / Non-human"], singleton_summary["total_reads"]),
            "Unclassified": format_count(singleton_summary["Unclassified"], singleton_summary["total_reads"]),
            "dominant_species_identity": f"{dominant_species_name} (taxid {dominant_taxid})" if dominant_species else "No viral species detected",
        }
    else:
        singleton_row = {
            "sample": sample_name,
            "read_type": "singleton",
            "Human": "0 (0.0%)",
            "Viral": "0 (0.0%)",
            "Dominant viral species": "0 (0.0%)",
            "Non-viral / Non-human": "0 (0.0%)",
            "Unclassified": "0 (0.0%)",
            "dominant_species_identity": "No singleton reads",
        }

    print(f"[RESULTS]: Read groups written to {sample_out_dir}")
    print(f"[DOMINANT SPECIES]: {dominant_species_name} (taxid {dominant_taxid})")
    return [paired_row, singleton_row]


def validate_inputs(config):
    """Validate sample names and expected trimmed read directories before analysis."""
    sample_names = config.cross_module_params["samples"]
    trimmed_root = os.path.join(config.cross_module_params["output_dir"], "module_1", "fastp")

    for sample_name in sample_names:
        sample_dir = os.path.join(trimmed_root, sample_name)
        if not os.path.exists(sample_dir):
            raise FileNotFoundError(
                f"Trimmed reads directory not found for sample {sample_name}: {sample_dir}. "
                f"Run module_1_qc_trim.py before module_2_classify_scrub.py."
            )

        read1, read2, singleton = find_trimmed_reads(
            sample_dir=sample_dir,
            sample_name=sample_name,
            read_1_str=config.cross_module_params["read_1_str"],
            read_2_str=config.cross_module_params["read_2_str"],
            singleton_str=config.cross_module_params.get("singleton_str"),
        )

        if read1 is None or read2 is None:
            raise FileNotFoundError(
                f"Trimmed paired reads not found for sample {sample_name} in {sample_dir}."
            )

    print("All sample directories and trimmed read inputs validated successfully.")


def main():
    parser = argparse.ArgumentParser(
        description="Classify trimmed reads with Kraken2 and split them into non-viral vs viral+unclassified groups."
    )
    parser.add_argument("--config", required=True, help="Path to the config.py file")
    args = parser.parse_args()

    config = load_config(args.config)

    print(f"""
Starting module_2_classify_scrub...

Config file location: {args.config}
Cross-module parameters: {config.cross_module_params}
Module 2 parameters: {config.module_2_params}
""")

    validate_inputs(config)

    output_root = os.path.join(config.cross_module_params["output_dir"], "module_2")
    os.makedirs(output_root, exist_ok=True)

    report_rows = []
    for sample_name in config.cross_module_params["samples"]:
        try:
            report_rows.extend(process_sample(sample_name, config))
        except Exception as exc:
            print(f"[ERROR]: Failed to process sample {sample_name}: {exc}")

    report_path = os.path.join(output_root, "read_classification_report.tsv")
    if os.path.exists(os.path.join(output_root, "report.txt")):
        os.remove(os.path.join(output_root, "report.txt"))
    write_summary_report(report_rows, report_path)
    print(f"[REPORT]: Wrote module 2 summary report to {report_path}")

    # stop daemon
    cmd = [
        config.module_2_params["Kraken2_k2_bin_path"],
        "clean",
        "--stop-daemon",
    ]
    try:
        subprocess.run(cmd, check=True)
    except Exception:
        pass

    print("\n[module_2_classify_scrub] Finished all samples.")


if __name__ == "__main__":
    main()

import argparse
import csv
import importlib.util
import json
import os
import shutil
import subprocess
from typing import List, Optional

from Bio import SeqIO


def load_config(config_path):
    """Load a Python config file as a Python module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")
    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def get_output_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_5", sample_name)

def get_module_3_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_3",  sample_name)

def get_module_4_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_4", sample_name)

def find_query_assembly(config, sample_name):
    """Return the assembled contigs/scaffolds FASTA for a sample."""
    sample_dir = get_module_3_dir(config, sample_name)

    path = os.path.join(sample_dir, "combined_transcripts.fasta")
    if os.path.exists(path):
        return path
    raise FileNotFoundError(f"No assembly FASTA found for sample {sample_name} in {sample_dir}.")

def get_reference_summary_path(config, sample_name):
    return os.path.join(get_module_4_dir(config, sample_name), f"{sample_name}_reference_summary.json")

def get_selected_reference_fasta(config, sample_name):
    """Return the selected reference FASTA for a sample from module 4 output."""
    summary_path = get_reference_summary_path(config, sample_name)
    if not os.path.exists(summary_path):
        raise FileNotFoundError(f"Module 4 summary not found for sample {sample_name}: {summary_path}")
    with open(summary_path, "r") as handle: summary = json.load(handle)
    if summary.get("is_segmented") is True or summary.get("selected_segments"): return None

    explicit_accession = get_explicit_ref_genome(config, sample_name)
    if explicit_accession:
        genome_accession = explicit_accession
    else:
        genome_accession = summary.get("genome_accession") or summary.get("best_reference_accession")
    if not genome_accession: raise ValueError(f"No genome accession recorded for {sample_name}.")

    module_4_ref_dir = os.path.join(config.cross_module_params["output_dir"], "module_4", "ref_genomes")
    target = os.path.join(module_4_ref_dir, f"{genome_accession}.fasta")
    if os.path.exists(target): return target

    symlink = os.path.join(get_module_4_dir(config, sample_name), f"selected_reference_{genome_accession}.fasta")
    if os.path.exists(symlink): return symlink
    raise FileNotFoundError(f"Reference FASTA not found for {sample_name}: {genome_accession}.")

def get_fasta_length(fasta_path):
    total_length = 0
    if not os.path.exists(fasta_path): return 0
    with open(fasta_path, "r") as handle:
        for line in handle:
            if not line.startswith(">"): total_length += len(line.strip())
    return total_length

def get_config_value(config, key, default=None):
    return config.module_5_params.get(key, default)


def get_explicit_ref_genome(config, sample_name):
    """Return an explicitly configured reference accession for a sample, if any."""
    raw = get_config_value(config, "explicit_ref_genome", {})
    if not isinstance(raw, dict):
        return None
    value = raw.get(sample_name)
    if value is None:
        return None
    return str(value).strip() or None


def get_excluded_contigs(config, sample_name):
    """Return contig IDs to ignore during scaffolding for a specific sample."""
    raw = get_config_value(config, "exclude_contigs", {})
    if not isinstance(raw, dict):
        return []
    sample_list = raw.get(sample_name, [])
    if sample_list is None:
        return []
    if isinstance(sample_list, str):
        sample_list = [sample_list]
    return [str(item) for item in sample_list]


def prepare_query_fasta_for_scaffolding(config, sample_name, query_fasta):
    """Filter out excluded contigs from the assembly used for scaffolding."""
    excluded = get_excluded_contigs(config, sample_name)
    if not excluded:
        return query_fasta

    sample_dir = get_output_dir(config, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    filtered_path = os.path.join(sample_dir, f"{sample_name}_filtered_scaffold_contigs.fasta")
    kept_records = []
    for record in SeqIO.parse(query_fasta, "fasta"):
        if record.id in excluded:
            print(f"[EXCLUDE CONTIG - {sample_name}]: Ignoring contig '{record.id}' from scaffolding.")
            continue
        kept_records.append(record)

    if not kept_records:
        raise ValueError(
            f"All contigs were excluded for sample {sample_name} in module_5 exclude_contigs: {excluded}"
        )

    with open(filtered_path, "w") as handle:
        SeqIO.write(kept_records, handle, "fasta")
    return filtered_path


def resolve_executable_path(config, key):
    configured = get_config_value(config, key, None)
    if configured is None:
        raise FileNotFoundError(
            f"Executable '{key}' is not configured. Set an explicit path or command in module_5_params."
        )
    if os.path.exists(configured):
        return configured
    resolved = shutil.which(configured)
    if resolved:
        return resolved
    raise FileNotFoundError(
        f"Configured executable '{configured}' for '{key}' was not found. Set an explicit path in module_5_params."
    )


# =====================================================================
# MINIMAP2 EXECUTION & PARSING
# =====================================================================

def parse_paf_alignments(paf_file, min_identity=70.0):
    """Parses PAF file into a list of alignment dictionaries. Converts to 1-based inclusive for stitching."""
    alignments = []
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return alignments

    with open(paf_file, "r") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 12: continue
            
            try:
                matches = int(parts[9])
                block_len = int(parts[10])
                if block_len == 0: continue
                identity = (matches / block_len) * 100.0
            except ValueError: continue

            if identity < min_identity: continue

            alignments.append({
                "contig_id": parts[0],
                "r_start": int(parts[7]) + 1,
                "r_end": int(parts[8]),
                "q_min": int(parts[2]) + 1,
                "q_max": int(parts[3]),
                "strand": parts[4],
                "identity": identity,
                "ref_cov_raw": (int(parts[8]) - int(parts[7]))
            })
    return alignments

def calculate_total_reference_coverage(alignments, ref_length):
    """Calculates the merged reference coverage from a list of alignments."""
    if not alignments or ref_length == 0: return 0.0
    intervals = sorted([(a['r_start'], a['r_end']) for a in alignments])
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1: merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else: merged.append((start, end))
    covered_bases = sum(end - start + 1 for start, end in merged)
    return covered_bases / ref_length


# =====================================================================
# ACCEPTANCE THRESHOLDS & ALIGNMENT MODES
# =====================================================================

# Alignment modes, ordered from most to least stringent.
ALIGNMENT_MODES = [
    {"name": "conserved",    "args": ["-c", "-x", "asm5"],         "min_identity": 70.0},
    {"name": "divergent_70", "args": ["-c", "-k", "14", "-w", "5"], "min_identity": 70.0},
    {"name": "divergent_60", "args": ["-c", "-k", "14", "-w", "5"], "min_identity": 60.0},
]

# Slice: a single contiguous block must cover >=98% of the reference;
# the resulting scaffold may contain at most 5% Ns (mostly terminal padding).
SLICE_MIN_REF_COVERAGE = 0.98
SLICE_MAX_N_PCT = 5.0

# Stitch: merged blocks must cover >=98% of the reference and the
# resulting scaffold may contain at most 2% Ns (stricter: joins are inferred).
STITCH_MIN_REF_COVERAGE = 0.98
STITCH_MAX_N_PCT = 2.0


# =====================================================================
# STRATEGY 1: BEST SLICE EXTRACTION (With Terminal Gap-Filling)
# =====================================================================

def try_best_slice_extraction(paf_file, query_fasta, ref_length, sample_name, sample_dir, suffix="", min_identity=70.0):
    alignments = parse_paf_alignments(paf_file, min_identity=min_identity)
    if not alignments: return None

    best_aln = max(alignments, key=lambda x: x['ref_cov_raw'])
    max_ref_cov = best_aln['ref_cov_raw'] / ref_length

    # Trigger fallback if coverage is below the slice threshold
    if max_ref_cov < SLICE_MIN_REF_COVERAGE:
        return None

    print(f"[SLICE MODE - {sample_name}{suffix}]: Single contig {best_aln['contig_id']} covers {max_ref_cov*100:.1f}% of reference. Extracting slice...")

    q_start_0based = best_aln['q_min'] - 1
    q_end_0based = best_aln['q_max']
    
    contig_seq = next((rec.seq for rec in SeqIO.parse(query_fasta, "fasta") if rec.id == best_aln['contig_id']), None)
    if contig_seq is None: raise RuntimeError(f"Contig {best_aln['contig_id']} missing from assembly.")
    
    sliced_seq = contig_seq[q_start_0based:q_end_0based]
    if best_aln['strand'] == "-": sliced_seq = sliced_seq.reverse_complement()
        
    # GAP FILL TERMINAL ENDS (Matches stitch mode behavior)
    left_gap = "N" * (best_aln['r_start'] - 1)
    right_gap = "N" * (ref_length - best_aln['r_end'])
    final_seq = left_gap + str(sliced_seq) + right_gap
    
    scaffold_fasta = os.path.join(sample_dir, f"{sample_name}{suffix}_scaffold.fasta")
    with open(scaffold_fasta, "w") as handle:
        handle.write(f">{sample_name}{suffix}_scaffold\n{final_seq}\n")
        
    print(f"[SLICE SUCCESS - {sample_name}{suffix}]: Extracted {len(sliced_seq)} bp slice (Total with terminal gaps: {len(final_seq)} bp).")
    return {
        "scaffold_fasta": scaffold_fasta,
        "scaffold_type": f"Sliced ({best_aln['contig_id']})"
    }


# =====================================================================
# STRATEGY 2: MULTI-CONTIG STITCHING
# =====================================================================

def write_stitch_contig_locations(sample_name, sample_dir, suffix, contig_locations):
    """Write a TSV summarising the reference intervals contributed by each stitched contig."""
    summary_path = os.path.join(sample_dir, f"{sample_name}{suffix}_stitch_contig_locations.tsv")
    fieldnames = [
        "order", "contig_id", "reference_start", "reference_end", "query_start", "query_end",
        "strand", "final_scaffold_start", "final_scaffold_end", "length"
    ]
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for entry in contig_locations:
            writer.writerow({key: entry.get(key, "") for key in fieldnames})
    return summary_path


def stitch_contigs_from_paf(paf_file, query_fasta, ref_length, sample_name, sample_dir, suffix="", min_identity=70.0):
    alignments = parse_paf_alignments(paf_file, min_identity=min_identity)
    if not alignments:
        return None  # Return None to allow fallback loop to continue

    alignments.sort(key=lambda x: x['r_start'])

    uncontained = []
    for i, aln in enumerate(alignments):
        is_nested = False
        for j, other in enumerate(alignments):
            if i == j: continue
            if other['r_start'] <= aln['r_start'] and other['r_end'] >= aln['r_end']:
                if other['r_start'] == aln['r_start'] and other['r_end'] == aln['r_end']:
                    if j < i: is_nested = True; break
                else: is_nested = True; break
        if not is_nested: uncontained.append(aln)

    if not uncontained:
        return None

    resolved = []
    for aln in uncontained:
        current = aln
        while resolved and resolved[-1]['r_end'] >= current['r_start']:
            prev = resolved[-1]
            overlap = prev['r_end'] - current['r_start'] + 1

            if prev['identity'] >= current['identity']:
                current['r_start'] += overlap
                if current['strand'] == "+": current['q_min'] += overlap
                else: current['q_max'] -= overlap
                if current['q_min'] > current['q_max'] or current['r_start'] > current['r_end']:
                    current = None; break
            else:
                prev['r_end'] -= overlap
                if prev['strand'] == "+": prev['q_max'] -= overlap
                else: prev['q_min'] += overlap
                if prev['q_min'] > prev['q_max'] or prev['r_start'] > prev['r_end']:
                    resolved.pop()

        if current is not None: resolved.append(current)

    if not resolved:
        return None

    contigs_by_id = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))
    scaffold_chunks = []
    contig_locations = []
    prev_r_end = 0
    used_contigs = []
    final_offset = 1

    for idx, aln in enumerate(resolved):
        contig_id = aln['contig_id']
        if contig_id not in contigs_by_id: continue
        used_contigs.append(contig_id)

        contig_record = contigs_by_id[contig_id]
        seq_len = len(contig_record.seq)

        q_start_idx = max(0, aln['q_min'] - 1)
        q_end_idx = min(seq_len, aln['q_max'])

        clipped_seq = str(contig_record.seq[q_start_idx:q_end_idx])
        if aln['strand'] == "-": clipped_seq = str(contig_record.seq[q_start_idx:q_end_idx].reverse_complement())

        if idx == 0 and aln['r_start'] > 1:
            scaffold_chunks.append("N" * (aln['r_start'] - 1))
            final_offset += aln['r_start'] - 1
        else:
            gap_to_prev = aln['r_start'] - prev_r_end - 1
            if gap_to_prev > 0:
                scaffold_chunks.append("N" * gap_to_prev)
                final_offset += gap_to_prev

        scaffold_chunks.append(clipped_seq)
        contig_locations.append({
            "order": len(contig_locations) + 1,
            "contig_id": contig_id,
            "reference_start": aln['r_start'],
            "reference_end": aln['r_end'],
            "query_start": aln['q_min'],
            "query_end": aln['q_max'],
            "strand": aln['strand'],
            "final_scaffold_start": final_offset,
            "final_scaffold_end": final_offset + len(clipped_seq) - 1,
            "length": len(clipped_seq),
        })
        final_offset += len(clipped_seq)
        prev_r_end = max(prev_r_end, aln['r_end'])

    if prev_r_end < ref_length:
        scaffold_chunks.append("N" * (ref_length - prev_r_end))

    final_sequence = "".join(scaffold_chunks)
    if not final_sequence:
        return None

    scaffold_fasta = os.path.join(sample_dir, f"{sample_name}{suffix}_scaffold.fasta")
    with open(scaffold_fasta, "w") as handle:
        handle.write(f">{sample_name}{suffix}_scaffold\n{final_sequence}\n")

    location_summary_path = os.path.join(sample_dir, f"{sample_name}{suffix}_stitch_contig_locations.tsv")
    print(f"[STITCH SUCCESS - {sample_name}{suffix}]: Stitched {len(resolved)} contigs into {len(final_sequence)} bp scaffold.")
    return {
        "scaffold_fasta": scaffold_fasta,
        "scaffold_type": f"Stitched ({', '.join(used_contigs)})",
        "contig_locations": contig_locations,
        "contig_location_summary": location_summary_path,
    }


# =====================================================================
# ORCHESTRATOR (6-Tier Fallback Loop)
# =====================================================================

def build_tier_order(strategy):
    """Return the evaluation order as a list of (mode_dict, technique) tuples.

    technique_first (default):
        conserved-slice -> divergent_70-slice -> divergent_60-slice
        -> conserved-stitch -> divergent_70-stitch -> divergent_60-stitch
        (Prioritises a single contiguous alignment block over inferred joins.)

    stringency_first (legacy):
        conserved-slice -> conserved-stitch -> divergent_70-slice
        -> divergent_70-stitch -> divergent_60-slice -> divergent_60-stitch
        (Prioritises alignment stringency over scaffolding technique.)
    """
    if strategy == "stringency_first":
        return [(mode, tech) for mode in ALIGNMENT_MODES for tech in ("slice", "stitch")]
    # Default: technique_first
    return [(mode, tech) for tech in ("slice", "stitch") for mode in ALIGNMENT_MODES]


def run_minimap2_scaffolding(config, ref_fasta, query_fasta, sample_name, sample_dir, segment_id=None):
    """Orchestrates the 6-tier fallback loop.

    The tier order is controlled by module_5_params["fallback_strategy"]:
      - "technique_first"  (default): slice at all stringencies, then stitch at all stringencies.
      - "stringency_first" (legacy):  slice+stitch per stringency, from most to least stringent.

    Minimap2 PAF outputs are cached per alignment mode (lazily), so each mode
    is aligned at most once regardless of tier order.
    """
    minimap2_bin = resolve_executable_path(config, "minimap2_bin_path")
    suffix = f"_seg_{segment_id}" if segment_id else ""
    ref_length = sum(len(record.seq) for record in SeqIO.parse(ref_fasta, "fasta"))

    strategy = get_config_value(config, "fallback_strategy", "technique_first")
    tier_order = build_tier_order(strategy)
    print(f"[SCAFFOLDING - {sample_name}{suffix}]: Strategy = {strategy}. Tier order: "
          + " -> ".join(f"{m['name']}-{t}" for m, t in tier_order))

    paf_cache = {}
    accepted_paf = None

    def get_paf(mode):
        """Lazily run minimap2 for an alignment mode and cache the PAF path."""
        if mode["name"] not in paf_cache:
            paf_file = os.path.join(sample_dir, f"{sample_name}_minimap2{suffix}_{mode['name']}.paf")
            print(f"[MINIMAP2 - {sample_name}{suffix}]: Running in {mode['name']} mode...")
            subprocess.run(
                [minimap2_bin] + mode["args"] + ["-o", paf_file, ref_fasta, query_fasta],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            paf_cache[mode["name"]] = paf_file
        return paf_cache[mode["name"]]

    try:
        for tier_idx, (mode, technique) in enumerate(tier_order, start=1):
            mode_name = mode["name"]
            paf_file = get_paf(mode)

            if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
                print(f"[TIER {tier_idx}/6 - {mode_name}-{technique}]: No alignments found. Skipping.")
                continue

            if technique == "slice":
                result = try_best_slice_extraction(
                    paf_file, query_fasta, ref_length, sample_name, sample_dir, suffix,
                    min_identity=mode["min_identity"]
                )
                if not result:
                    print(f"[TIER {tier_idx}/6 - {mode_name}-slice]: No single contig covers "
                          f">= {SLICE_MIN_REF_COVERAGE*100:.0f}% of reference. Skipping.")
                    continue

                seq = str(next(SeqIO.parse(result["scaffold_fasta"], "fasta")).seq)
                n_pct = (seq.upper().count('N') / len(seq) * 100) if seq else 100.0

                if n_pct <= SLICE_MAX_N_PCT:
                    result["alignment_mode"] = mode_name
                    print(f"[TIER {tier_idx}/6 - {mode_name}-slice]: ACCEPTED ({n_pct:.1f}% Ns).")
                    accepted_paf = paf_file
                    if "contig_locations" in result and result["contig_locations"]:
                        result["contig_location_summary"] = write_stitch_contig_locations(
                            sample_name, sample_dir, suffix, result["contig_locations"]
                        )
                        print(f"[STITCH LOCATIONS - {sample_name}{suffix}]: Wrote contig contribution map to {result['contig_location_summary']}")
                    return result
                else:
                    print(f"[TIER {tier_idx}/6 - {mode_name}-slice]: REJECTED ({n_pct:.1f}% Ns > {SLICE_MAX_N_PCT}%). Discarding.")
                    os.remove(result["scaffold_fasta"])
                    if result.get("contig_location_summary") and os.path.exists(result["contig_location_summary"]):
                        os.remove(result["contig_location_summary"])

            else:  # stitch
                result = stitch_contigs_from_paf(
                    paf_file, query_fasta, ref_length, sample_name, sample_dir, suffix,
                    min_identity=mode["min_identity"]
                )
                if not result:
                    print(f"[TIER {tier_idx}/6 - {mode_name}-stitch]: No stitchable alignments. Skipping.")
                    continue

                seq = str(next(SeqIO.parse(result["scaffold_fasta"], "fasta")).seq)
                n_pct = (seq.upper().count('N') / len(seq) * 100) if seq else 100.0
                alignments = parse_paf_alignments(paf_file, min_identity=mode["min_identity"])
                coverage = calculate_total_reference_coverage(alignments, ref_length)

                if coverage >= STITCH_MIN_REF_COVERAGE and n_pct <= STITCH_MAX_N_PCT:
                    result["alignment_mode"] = mode_name
                    print(f"[TIER {tier_idx}/6 - {mode_name}-stitch]: ACCEPTED "
                          f"({coverage*100:.1f}% coverage, {n_pct:.1f}% Ns).")
                    accepted_paf = paf_file
                    if "contig_locations" in result and result["contig_locations"]:
                        result["contig_location_summary"] = write_stitch_contig_locations(
                            sample_name, sample_dir, suffix, result["contig_locations"]
                        )
                        print(f"[STITCH LOCATIONS - {sample_name}{suffix}]: Wrote contig contribution map to {result['contig_location_summary']}")
                    return result
                else:
                    print(f"[TIER {tier_idx}/6 - {mode_name}-stitch]: REJECTED "
                          f"({coverage*100:.1f}% coverage, {n_pct:.1f}% Ns). Discarding.")
                    os.remove(result["scaffold_fasta"])
                    if result.get("contig_location_summary") and os.path.exists(result["contig_location_summary"]):
                        os.remove(result["contig_location_summary"])

        raise RuntimeError(f"All 6 alignment/scaffolding strategies failed for {sample_name}{suffix}.")
    finally:
        # Keep the PAF for the tier that produced the accepted scaffold; clean up the others.
        for cached_paf in paf_cache.values():
            if cached_paf != accepted_paf and os.path.exists(cached_paf):
                os.remove(cached_paf)


def run_checkv_qc(config, target_fasta, checkv_dir):
    """Run CheckV on the scaffolded FASTA and return aggregated completeness/quality metrics."""
    checkv_bin = resolve_executable_path(config, "checkv_bin_path")
    base_threads = int(get_config_value(config, "threads", 8))
    checkv_db_dir = get_config_value(config, "checkv_db_dir", None)
    os.makedirs(checkv_dir, exist_ok=True)
    
    candidate_threads = [base_threads] if base_threads <= 8 else [8, max(1, base_threads // 2), base_threads]
    last_error = None
    
    for threads in candidate_threads:
        cmd = [checkv_bin, "end_to_end", target_fasta, checkv_dir, "-t", str(threads), "--restart"]
        if checkv_db_dir and os.path.exists(checkv_db_dir): 
            cmd.extend(["-d", checkv_db_dir])
        
        print(f"[CHECKV]: Running CheckV for {target_fasta} with {threads} threads")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        lower_output = (proc.stdout + proc.stderr).lower()
        
        is_failure = proc.returncode != 0 or "error:" in lower_output or "hmmsearch tasks failed" in lower_output
        if not is_failure: 
            break
        
        last_error = (proc.stdout + proc.stderr).strip()
        if "hmmsearch tasks failed" in lower_output: 
            return "N/A", "N/A", "N/A", last_error
        if threads != candidate_threads[-1]: 
            print(f"[CHECKV WARN]: Retrying with fewer threads...")
    else: 
        raise RuntimeError(f"CheckV failed.\nError: {last_error}")

    summary_tsv = None
    for root, _, files in os.walk(checkv_dir):
        if "quality_summary.tsv" in files:
            summary_tsv = os.path.join(root, "quality_summary.tsv")
            break
            
    if not summary_tsv: 
        return "N/A", "N/A", "N/A", "quality_summary.tsv missing"

    with open(summary_tsv, "r") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        
    if not rows: 
        return "N/A", "N/A", "N/A", "quality_summary.tsv empty"

    # --- BULLETPROOF AGGREGATION ---
    completions = []
    for r in rows:
        comp_str = str(r.get("completeness", "")).strip()
        if comp_str and comp_str not in ["N/A", "NA", ""]:
            try:
                completions.append(float(comp_str))
            except ValueError:
                pass
                
    avg_comp = sum(completions) / len(completions) if completions else "N/A"
    
    quality_rank = {"High-quality": 4, "Medium-quality": 3, "Low-quality": 2, "Not-determined": 1}
    valid_qualities = []
    for r in rows:
        q = str(r.get("checkv_quality", r.get("quality", ""))).strip()
        if q in quality_rank:
            valid_qualities.append(q)
    best_quality = max(valid_qualities, key=lambda q: quality_rank[q]) if valid_qualities else "Not-determined"
    
    kmer_freqs = []
    for r in rows:
        kf_str = str(r.get("kmer_freq", "")).strip()
        if kf_str and kf_str not in ["N/A", "NA", ""]:
            try:
                kmer_freqs.append(float(kf_str))
            except ValueError:
                pass
                
    avg_kmer = sum(kmer_freqs) / len(kmer_freqs) if kmer_freqs else "N/A"
    
    warning_list = []
    for r in rows:
        w = str(r.get("warnings", "")).strip()
        if w and w not in ["N/A", "NA", "None", ""]:
            warning_list.append(w)
    warnings = "; ".join(warning_list) if warning_list else "None"
    
    comp_out = str(round(avg_comp, 2)) if avg_comp != "N/A" else "N/A"
    kmer_out = str(round(avg_kmer, 4)) if avg_kmer != "N/A" else "N/A"
    
    return comp_out, best_quality, kmer_out, warnings


def write_qc_report(config, sample_name, reference_fasta, scaffold_fasta, scaffold_len, reference_len, checkv_stats, sample_dir, is_concatenated=False, scaffold_type="Unknown", alignment_mode="Unknown", n_info="N/A", contig_location_summary=None):
    completeness, quality, kmer, warnings = checkv_stats
    report_path = os.path.join(sample_dir, f"{sample_name}_scaffold_checkv_report.txt")
    
    with open(report_path, "w") as handle:
        handle.write(f"==================================================\n")
        handle.write(f" Minimap2 Scaffold QC Report: {sample_name}\n")
        handle.write(f"==================================================\n")
        handle.write(f"Reference FASTA:      {reference_fasta}\n")
        handle.write(f"Reference Length:     {reference_len} bp\n")
        handle.write(f"Scaffold FASTA:       {scaffold_fasta}\n")
        handle.write(f"Scaffold Length:      {scaffold_len} bp\n\n")
        handle.write(f"--- Assembly & Alignment Metrics ---\n")
        handle.write(f"Scaffold Type:      {scaffold_type}\n")
        handle.write(f"Alignment Mode:     {alignment_mode.capitalize()}\n")
        handle.write(f"Gap-filled Bases:   {n_info}\n")
        if contig_location_summary and os.path.exists(contig_location_summary):
            handle.write(f"Contig Locations:    {contig_location_summary}\n")
        handle.write(f"\n")
        handle.write(f"--- CheckV Metrics ---\n")
        if is_concatenated:
            handle.write(f"Completeness:       {completeness}% (Estimated from concatenated segments)\n")
        else:
            handle.write(f"Completeness:       {completeness}%\n")
        handle.write(f"Quality:            {quality}\n")
        handle.write(f"K-mer Freq:         {kmer}\n")
        handle.write(f"Warnings:           {warnings}\n")
        
        if is_concatenated:
            handle.write(f"\n[NOTE]: Individual segment CheckV metrics are available in:\n")
            handle.write(f"        {os.path.join(sample_dir, 'checkv_segments/quality_summary.tsv')}\n")
        
    print(f"[QC REPORT]: Wrote {report_path}")
    return report_path


# =====================================================================
# SAMPLE PROCESSORS
# =====================================================================

def process_non_segmented_sample(sample_name, config, summary, query_fasta):
    sample_dir = get_output_dir(config, sample_name)
    reference_fasta = get_selected_reference_fasta(config, sample_name)
    
    print(f"\n==========================================\n Processing non-segmented sample: {sample_name}\n==========================================")
    
    scaffold_query_fasta = prepare_query_fasta_for_scaffolding(config, sample_name, query_fasta)
    scaffold_result = run_minimap2_scaffolding(config, reference_fasta, scaffold_query_fasta, sample_name, sample_dir)
    scaffold_fasta = scaffold_result["scaffold_fasta"]
    
    # Calculate gap-filled bases
    scaffold_seq = str(next(SeqIO.parse(scaffold_fasta, "fasta")).seq)
    scaffold_len = len(scaffold_seq)
    n_count = scaffold_seq.upper().count('N')
    n_pct = (n_count / scaffold_len * 100) if scaffold_len > 0 else 0
    n_info = f"{n_count} ({n_pct:.1f}%)"
    
    reference_len = get_fasta_length(reference_fasta)

    checkv_stats = run_checkv_qc(config, scaffold_fasta, os.path.join(sample_dir, "checkv"))
    write_qc_report(config, sample_name, reference_fasta, scaffold_fasta, scaffold_len, reference_len, checkv_stats, sample_dir, 
                    scaffold_type=scaffold_result["scaffold_type"], 
                    alignment_mode=scaffold_result["alignment_mode"],
                    n_info=n_info,
                    contig_location_summary=scaffold_result.get("contig_location_summary"))

    return {
        "sample_accession": sample_name, "reference_fasta": reference_fasta, "scaffold_fasta": scaffold_fasta,
        "reference_length": reference_len, "scaffold_length": scaffold_len,
        "checkv_completeness": checkv_stats[0], "checkv_quality": checkv_stats[1],
        "checkv_kmer_freq": checkv_stats[2], "checkv_warnings": checkv_stats[3],
        "scaffold_type": scaffold_result["scaffold_type"],
        "alignment_mode": scaffold_result["alignment_mode"],
        "gap_filled_bases": n_info,
        "contig_location_summary": scaffold_result.get("contig_location_summary", "")
    }


def process_segmented_sample(sample_name, config, summary, query_fasta):
    sample_dir = get_output_dir(config, sample_name)
    selected_segments = summary.get("selected_segments", [])
    
    if not selected_segments:
        print(f"[ERROR]: Sample {sample_name} is marked as segmented but no segments found in summary.")
        return None

    print(f"\n==========================================\n Processing SEGMENTED sample: {sample_name} ({len(selected_segments)} segments)\n==========================================")
    
    final_scaffold_records = []
    total_ref_len = 0
    segment_details = []
    total_n_count = 0
    total_scaffold_len = 0
    
    for seg in selected_segments:
        seg_id = seg["segment_seq_id"]
        ref_fasta = seg["fasta_path"]
        
        print(f"\n[SEGMENT PROCESSING - {sample_name}]: Processing segment {seg_id}...")
        scaffold_query_fasta = prepare_query_fasta_for_scaffolding(config, sample_name, query_fasta)
        scaffold_result = run_minimap2_scaffolding(config, ref_fasta, scaffold_query_fasta, sample_name, sample_dir, segment_id=seg_id)
        
        for record in SeqIO.parse(scaffold_result["scaffold_fasta"], "fasta"):
            record.id = f"{sample_name}_{seg_id}"
            record.description = f"segment={seg_id}"
            final_scaffold_records.append(record)
            
            seq = str(record.seq).upper()
            total_scaffold_len += len(seq)
            total_n_count += seq.count('N')
            total_ref_len += get_fasta_length(ref_fasta)
            
        segment_details.append(f"{seg_id}: {scaffold_result['scaffold_type']} ({scaffold_result['alignment_mode']})")
        if scaffold_result.get("contig_location_summary"):
            print(f"[SEGMENT LOCATIONS - {sample_name}]: {scaffold_result['contig_location_summary']}")

    # 1. Write multi-FASTA for downstream tools (iVar, MAFFT) and individual segment CheckV
    combined_scaffold_fasta = os.path.join(sample_dir, f"{sample_name}_segmented_scaffold.fasta")
    SeqIO.write(final_scaffold_records, combined_scaffold_fasta, "fasta")
    print(f"[SUCCESS - {sample_name}]: Wrote combined segmented scaffold to {combined_scaffold_fasta}")
    
    # 2. Run CheckV on multi-FASTA to generate individual segment metrics for user inspection
    checkv_dir_multi = os.path.join(sample_dir, "checkv_segments")
    print(f"[CHECKV]: Running CheckV on individual segments for detailed inspection...")
    run_checkv_qc(config, combined_scaffold_fasta, checkv_dir_multi)

    # 3. Create a concatenated sequence for overall completeness estimation
    spacer = "NNNNN"
    concatenated_seq = spacer.join([str(rec.seq) for rec in final_scaffold_records])
    concatenated_fasta = os.path.join(sample_dir, f"{sample_name}_concatenated_scaffold.fasta")
    with open(concatenated_fasta, "w") as handle:
        handle.write(f">{sample_name}_concatenated_genome\n{concatenated_seq}\n")
        
    # 4. Run CheckV on the concatenated sequence for the main report
    checkv_dir_concat = os.path.join(sample_dir, "checkv_concatenated")
    print(f"[CHECKV]: Running CheckV on concatenated genome for overall completeness estimation...")
    checkv_stats = run_checkv_qc(config, concatenated_fasta, checkv_dir_concat)
    
    # Calculate gap-filled bases (based on actual assembled segments, excluding artificial spacers)
    n_pct = (total_n_count / total_scaffold_len * 100) if total_scaffold_len > 0 else 0
    n_info = f"{total_n_count} ({n_pct:.1f}%)"
    segment_summary_str = "; ".join(segment_details)
    
    write_qc_report(config, sample_name, "Multi-segment Reference (see JSON)", combined_scaffold_fasta, total_scaffold_len, total_ref_len, checkv_stats, sample_dir, 
                    is_concatenated=True, scaffold_type=segment_summary_str, alignment_mode="Mixed (per segment)", n_info=n_info)

    return {
        "sample_accession": sample_name, 
        "reference_fasta": "Multi-segment (see summary JSON)", 
        "scaffold_fasta": combined_scaffold_fasta,
        "reference_length": total_ref_len, 
        "scaffold_length": total_scaffold_len,
        "checkv_completeness": checkv_stats[0], 
        "checkv_quality": checkv_stats[1],
        "checkv_kmer_freq": checkv_stats[2], 
        "checkv_warnings": checkv_stats[3],
        "scaffold_type": segment_summary_str,
        "alignment_mode": "Mixed (per segment)",
        "gap_filled_bases": n_info,
        "contig_location_summary": ""
    }


def process_sample(sample_name, config):
    sample_dir = get_output_dir(config, sample_name)
    os.makedirs(sample_dir, exist_ok=True)
    summary_path = get_reference_summary_path(config, sample_name)
    
    with open(summary_path, "r") as handle: 
        summary = json.load(handle)
    
    query_fasta = find_query_assembly(config, sample_name)
    
    if summary.get("is_segmented") is True and summary.get("selected_segments"):
        return process_segmented_sample(sample_name, config, summary, query_fasta)
    else:
        return process_non_segmented_sample(sample_name, config, summary, query_fasta)


def validate_inputs(config):
    for sample_name in config.cross_module_params["samples"]:
        summary_path = get_reference_summary_path(config, sample_name)
        if not os.path.exists(summary_path): 
            raise FileNotFoundError(f"Module 4 summary missing for {sample_name}")
        with open(summary_path, "r") as handle: 
            summary = json.load(handle)
            
        find_query_assembly(config, sample_name)
        
        if summary.get("is_segmented") is True and summary.get("selected_segments"):
            for seg in summary["selected_segments"]:
                if not os.path.exists(seg["fasta_path"]):
                    raise FileNotFoundError(f"Segment FASTA missing for {sample_name}: {seg['fasta_path']}")
            print(f"[VALIDATION]: Sample {sample_name} has valid reference segments identified in module 4.")
        else:
            get_selected_reference_fasta(config, sample_name)
            print(f"[VALIDATION]: Sample {sample_name} has valid reference genome identified in module 4.")

    print("All module 4 selections validated successfully.")


def main():
    parser = argparse.ArgumentParser(description="Minimap2 hybrid scaffold (Slice/Stitch) and CheckV QC.")
    parser.add_argument("--config", required=True, help="Path to the config Python file.")
    args = parser.parse_args()
    config = load_config(args.config)
    print(f"\nStarting module_5_scaffold_qc...\nConfig file location: {args.config}\n")
    
    validate_inputs(config)
    summary_rows = []
    for sample_name in config.cross_module_params["samples"]:
        try:
            result = process_sample(sample_name, config)
            if result: summary_rows.append(result)
        except Exception as exc: 
            print(f"[ERROR]: Failed for sample {sample_name}: {exc}")

    if summary_rows:
        module_5_dir = os.path.join(config.cross_module_params["output_dir"], "module_5")
        os.makedirs(module_5_dir, exist_ok=True)
        summary_path = os.path.join(module_5_dir, "Identified_scaffolds.tsv")
        fieldnames = ["sample_accession", "reference_fasta", "scaffold_fasta", "reference_length", 
                      "scaffold_length", "checkv_completeness", "checkv_quality", "checkv_kmer_freq", "checkv_warnings",
                      "scaffold_type", "alignment_mode", "gap_filled_bases", "contig_location_summary"]
        with open(summary_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in summary_rows: writer.writerow({key: row.get(key, "") for key in fieldnames})
        print(f"[SUMMARY]: Wrote scaffold QC table to {summary_path}")
    print("\n[module_5_scaffold_qc] Finished all samples.")

if __name__ == "__main__":
    main()
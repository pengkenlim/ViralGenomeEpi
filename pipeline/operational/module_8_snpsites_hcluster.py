import argparse
import gzip
import importlib.util
import io
import itertools
import json
import os
import re
import shutil
import subprocess
import urllib.request
import zipfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import AlignIO, SeqIO
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform


GLOBAL_COL = "Global SNP Difference"
BASE_CODES = {"A": 0, "C": 1, "G": 2, "T": 3}
LOW_COV_FRAC = 0.90


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
    value = get_config_value(config, "module_8_params", key, None)
    if value is None:
        raise FileNotFoundError(
            f"Executable '{key}' is not configured. Set an explicit path or command in module_8_params."
        )
    if os.path.exists(value):
        return value
    resolved = shutil.which(value)
    if resolved:
        return resolved
    raise FileNotFoundError(
        f"Configured executable '{value}' for '{key}' was not found. Set an explicit path in module_8_params."
    )


def get_msa_reference(config):
    """Return module 7's configured MSA reference accession, or None."""
    for key in ("msa_reference", "msa_ref", "MSA_reference"):
        value = get_config_value(config, "module_7_params", key, None)
        if value:
            return str(value).strip()
    return None


# ---------------------------------------------------------------------------
# Output resolution
# ---------------------------------------------------------------------------

def get_module_7_dir(config):
    analysis_name = (
        get_config_value(config, "module_7_params", "analysis_name", None)
        or get_config_value(config, "module_8_params", "analysis_name", "default")
        or "default"
    )
    return os.path.join(config.cross_module_params["output_dir"], "module_7", analysis_name)


def get_module_8_dir(config):
    analysis_name = (
        get_config_value(config, "module_8_params", "analysis_name", None)
        or get_config_value(config, "module_7_params", "analysis_name", "default")
        or "default"
    )
    return os.path.join(config.cross_module_params["output_dir"], "module_8", analysis_name)


def get_input_samples(config):
    samples = get_config_value(config, "module_8_params", "input_samples", None)
    if samples is None:
        samples = get_config_value(config, "module_7_params", "input_samples", None)
    if samples:
        return list(samples)
    return list(config.cross_module_params.get("samples", []))


def find_trimmed_alignment(config):
    module_7_dir = get_module_7_dir(config)
    if not os.path.isdir(module_7_dir):
        raise FileNotFoundError(f"Module 7 directory not found: {module_7_dir}")

    candidates = [
        os.path.join(module_7_dir, "polished_genomes_aligned_trimmed.fasta"),
        os.path.join(module_7_dir, "genomes_aligned_trimmed.fasta"),
        os.path.join(module_7_dir, "aligned_trimmed.fasta"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path

    matches = []
    for root, _, files in os.walk(module_7_dir):
        for name in files:
            if name.lower().endswith((".fa", ".fasta")) and "trim" in name.lower():
                matches.append(os.path.join(root, name))

    if matches:
        return sorted(matches)[0]

    raise FileNotFoundError(
        f"No trimmed alignment file found in module 7 output: {module_7_dir}"
    )


def find_ref_aligned_alignment(config):
    path = os.path.join(get_module_7_dir(config), "polished_genomes_ref_aligned.fasta")
    return path if os.path.isfile(path) else None


# ---------------------------------------------------------------------------
# Reference GFF acquisition (cached; three fallbacks)
# ---------------------------------------------------------------------------

def _extract_gff_from_zip_bytes(data, dest_path):
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        members = [m for m in zf.namelist() if m.endswith(".gff")]
        if not members:
            raise RuntimeError("ZIP contains no .gff member.")
        preferred = [m for m in members if m.endswith("genomic.gff")]
        pick = (preferred or members)[0]
        with zf.open(pick) as src, open(dest_path, "wb") as out:
            shutil.copyfileobj(src, out)
    return dest_path


def _download_gff_via_cli(config, accession, zip_path):
    binary = None
    for section in ("module_8_params", "module_7_params"):
        value = get_config_value(config, section, "datasets_bin_path", None)
        if value:
            binary = value
            break
    if binary is None:
        binary = shutil.which("datasets")
    if not binary:
        raise FileNotFoundError("datasets CLI not configured/found.")
    if not os.path.exists(binary):
        resolved = shutil.which(binary)
        if not resolved:
            raise FileNotFoundError(f"datasets CLI '{binary}' not found.")
        binary = resolved
    cmd = [binary, "download", "genome", "accession", accession, "--include", "gff", "--filename", zip_path]
    subprocess.run(cmd, check=True, capture_output=True)


def _download_gff_via_api(accession):
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{accession}/download?include_annotation_type=GENOME_GFF&filename={accession}_annotation.zip"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "ViralGenomeEpi-module8"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        return resp.read()


def _download_gff_via_ftp(accession):
    prefix, rest = accession.split("_", 1)
    digits = rest.split(".")[0].zfill(9)
    grouped = "/".join(digits[i:i + 3] for i in range(0, 9, 3))
    list_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{grouped}/"
    req = urllib.request.Request(list_url, headers={"User-Agent": "ViralGenomeEpi-module8"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        html = resp.read().decode("utf-8", errors="replace")
    dirs = re.findall(r'href="(' + re.escape(accession) + r'[^"/]*/)"', html)
    if not dirs:
        raise RuntimeError(f"No assembly directory matched {accession} at {list_url}")
    dirname = dirs[0].rstrip("/")
    gff_url = f"{list_url}{dirname}/{dirname}_genomic.gff.gz"
    req = urllib.request.Request(gff_url, headers={"User-Agent": "ViralGenomeEpi-module8"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        raw = resp.read()
    try:
        return gzip.decompress(raw)
    except OSError:
        return raw  # already uncompressed


def download_reference_gff(config, accession, output_dir):
    """Fetch the RefSeq GFF for the accession; cache under reference_annotation/."""
    cache_dir = os.path.join(output_dir, "reference_annotation")
    os.makedirs(cache_dir, exist_ok=True)
    cache_path = os.path.join(cache_dir, f"{accession}.gff")
    if os.path.isfile(cache_path):
        print(f"[GFF]: Using cached annotation {cache_path}")
        return cache_path

    zip_path = cache_path + ".zip"
    errors = []

    # 1) NCBI Datasets CLI
    try:
        _download_gff_via_cli(config, accession, zip_path)
        with open(zip_path, "rb") as fh:
            data = fh.read()
        _extract_gff_from_zip_bytes(data, cache_path)
        if os.path.exists(zip_path):
            os.remove(zip_path)
        print(f"[GFF]: Downloaded via datasets CLI -> {cache_path}")
        return cache_path
    except Exception as exc:
        errors.append(f"datasets CLI: {exc}")

    # 2) Datasets REST API
    try:
        data = _download_gff_via_api(accession)
        _extract_gff_from_zip_bytes(data, cache_path)
        print(f"[GFF]: Downloaded via Datasets REST API -> {cache_path}")
        return cache_path
    except Exception as exc:
        errors.append(f"Datasets API: {exc}")

    # 3) NCBI FTP (all-genomes path)
    try:
        raw = _download_gff_via_ftp(accession)
        with open(cache_path, "wb") as fh:
            fh.write(raw)
        print(f"[GFF]: Downloaded via NCBI FTP -> {cache_path}")
        return cache_path
    except Exception as exc:
        errors.append(f"FTP: {exc}")

    if os.path.exists(zip_path):
        os.remove(zip_path)
    raise RuntimeError("All GFF download methods failed: " + " | ".join(errors))


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def parse_gff_genes(gff_path):
    """Parse 'gene' features from a GFF; identifier from Name= (fallback locus_tag, ID)."""
    genes = []
    seen = set()
    with open(gff_path, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            attrs = fields[8]
            name = None
            for key in ("Name", "locus_tag", "ID"):
                match = re.search(rf"{key}=([^;\n]+)", attrs)
                if match:
                    name = match.group(1).strip()
                    break
            if not name:
                continue
            if name in seen:
                print(f"[GFF WARNING]: Duplicate gene name '{name}'; keeping first occurrence.")
                continue
            seen.add(name)
            genes.append({"name": name, "start": start, "end": end, "strand": fields[6]})
    genes.sort(key=lambda g: g["start"])
    return genes


# ---------------------------------------------------------------------------
# Gene-aware SNP counting on the reference-included (untrimmed) alignment
# ---------------------------------------------------------------------------

def write_gene_msas(records, char_arrays, genes, gene_cols_masks, output_dir):
    """Slice the untrimmed ref-included alignment into per-gene MSAs (no trimming)."""
    print("[GENE MSA]: Writing untrimmed per-gene alignment slices...")
    files = []
    used = set()
    for gene, cols in zip(genes, gene_cols_masks):
        strand_suffix = "plus" if gene["strand"] == "+" else "minus"
        safe_name = re.sub(r"[^\w.\-]", "_", gene["name"])
        fname = f"genome_ref_{safe_name}_{gene['start']}_{gene['end']}_{strand_suffix}.fasta"
        if fname in used:
            fname = f"genome_ref_{safe_name}_{gene['start']}_{gene['end']}_{strand_suffix}_dup.fasta"
        used.add(fname)
        path = os.path.join(output_dir, fname)
        with open(path, "w") as handle:
            for rec, arr in zip(records, char_arrays):
                sub = "".join(arr[cols].tolist())
                handle.write(f">{rec.id}\n{sub}\n")
        files.append(fname)
    print(f"[GENE MSA]: Wrote {len(files)} gene-specific MSAs to {output_dir}")
    return files


def build_gene_analysis(config, output_dir, input_samples):
    """Return a gene-analysis context dict, or None if gene-aware mode is unavailable."""
    msa_ref = get_msa_reference(config)
    if not msa_ref:
        print("[GENE-AWARE]: No 'msa_reference' in module_7_params; skipping gene-aware analysis.")
        return None
    if not re.match(r"^(GCF|GCA)_\d+(\.\d+)?$", msa_ref):
        print(f"[GENE-AWARE]: 'msa_reference' ('{msa_ref}') is not a GCF/GCA accession; skipping.")
        return None

    ref_aligned = find_ref_aligned_alignment(config)
    if ref_aligned is None:
        print("[GENE-AWARE]: polished_genomes_ref_aligned.fasta not found in module 7; skipping.")
        return None

    try:
        gff_path = download_reference_gff(config, msa_ref, output_dir)
    except Exception as exc:
        print(f"[GENE-AWARE WARNING]: GFF download failed ({exc}); skipping gene-aware analysis.")
        return None

    genes = parse_gff_genes(gff_path)
    if not genes:
        print("[GENE-AWARE WARNING]: No 'gene' features parsed from GFF; skipping.")
        return None

    records = list(AlignIO.read(ref_aligned, "fasta"))
    ref_id = msa_ref.replace(" ", "_")
    ref_idx = next((i for i, r in enumerate(records) if r.id == ref_id), None)
    if ref_idx is None:
        ref_idx = next((i for i, r in enumerate(records) if msa_ref in r.id), None)
    if ref_idx is None:
        print("[GENE-AWARE WARNING]: Reference row not found in ref-included alignment; skipping.")
        return None
    ref_id = records[ref_idx].id

    ids_all = [r.id for r in records]
    wanted = [s for s in input_samples if s in ids_all]
    missing = sorted(set(input_samples) - set(wanted))
    if missing:
        print(f"[GENE-AWARE WARNING]: Samples absent from ref-included alignment: {missing}")
    if not wanted:
        print("[GENE-AWARE WARNING]: No input samples present in ref-included alignment; skipping.")
        return None

    ids = wanted + [ref_id]
    idx_of = {rid: i for i, rid in enumerate(ids_all)}
    sel_idx = [idx_of[rid] for rid in ids]

    # Encode alignment: 0-3 = A/C/G/T, -2 = gap, -1 = N/ambiguous
    # FIX: Added .upper() to ensure lowercase bases are correctly recognized
    char_arrays = [np.array(list(str(r.seq).upper())) for r in records]
    n_all, L = len(records), len(char_arrays[0])
    M = np.full((n_all, L), -1, dtype=np.int8)
    for i, arr in enumerate(char_arrays):
        M[i][arr == "-"] = -2
        for ch, code in BASE_CODES.items():
            M[i][arr == ch] = code
    valid_all = M >= 0

    # Reference coordinate lift: column -> 1-based reference position (0 = ref gap column)
    ref_arr = char_arrays[ref_idx]
    nongap = ref_arr != "-"
    ref_pos = np.zeros(L, dtype=np.int64)
    ref_pos[nongap] = np.cumsum(nongap)[nongap]

    # Per-gene column masks, pairwise counts, and long-format rows
    gene_counts = {}   # (id_a, id_b) -> {gene_name: (snps, compared)}
    long_rows = []
    gene_cols_masks = []
    kept_genes = []
    for gene in genes:
        cols = np.where((ref_pos >= gene["start"]) & (ref_pos <= gene["end"]))[0]
        if cols.size == 0:
            print(f"[GENE-AWARE WARNING]: Gene '{gene['name']}' has no aligned columns; skipping.")
            continue
        kept_genes.append(gene)
        gene_cols_masks.append(cols)

        Mk = M[:, cols]
        Vk = valid_all[:, cols]
        for a, b in itertools.combinations(range(len(ids)), 2):
            ia, ib = sel_idx[a], sel_idx[b]
            cmp_mask = Vk[ia] & Vk[ib]
            snps = int(np.sum((Mk[ia] != Mk[ib]) & cmp_mask))
            compared = int(np.sum(cmp_mask))
            gene_counts.setdefault((ids[a], ids[b]), {})[gene["name"]] = (snps, compared)
            gene_counts.setdefault((ids[b], ids[a]), {})[gene["name"]] = (snps, compared)
            source = "sample-reference" if ref_id in (ids[a], ids[b]) else "sample-sample"
            long_rows.append({
                "Sample_1": ids[a],
                "Sample_2": ids[b],
                "Gene": gene["name"],
                "SNP_Distance": snps,
                "Compared_Bases": compared,
                "Gene_Aligned_Columns": int(cols.size),
                "Low_Coverage_Flag": bool(compared < LOW_COV_FRAC * cols.size),
                "Pair_Source": source,
            })

    # Whole-genome globals from the ref-included alignment (for sample-reference pairs)
    ref_global = {}
    for a in range(len(wanted)):
        ia, ib = sel_idx[a], sel_idx[ref_idx]
        cmp_mask = valid_all[ia] & valid_all[ib]
        ref_global[(wanted[a], ref_id)] = int(np.sum((M[ia] != M[ib]) & cmp_mask))

    # --- CONSISTENCY GUARD ---
    # Check if coordinate mapping is working: if global SNPs > 0 but gene SNPs sum to 0, warn.
    for pair, counts in gene_counts.items():
        if pair in ref_global:
            global_snps = ref_global[pair]
            gene_sum = sum(c[0] for c in counts.values())
            if global_snps > 10 and gene_sum == 0:
                print(f"[WARNING] Consistency check failed for pair {pair}: "
                      f"Global SNPs = {global_snps}, but sum of gene SNPs = {gene_sum}. "
                      f"Check reference coordinate mapping or GFF intervals.")

    # Untrimmed per-gene MSA slices (all sequences in the ref-included alignment)
    gene_msa_files = write_gene_msas(records, char_arrays, kept_genes, gene_cols_masks, output_dir)

    return {
        "msa_reference": msa_ref,
        "ref_id": ref_id,
        "gff_path": gff_path,
        "alignment_path": ref_aligned,
        "genes": kept_genes,
        "gene_counts": gene_counts,
        "ref_global": ref_global,
        "long_rows": long_rows,
        "gene_msa_files": gene_msa_files,
        "wanted_samples": wanted,
    }


# ---------------------------------------------------------------------------
# SNP clustering workflow (legacy core, unchanged semantics)
# ---------------------------------------------------------------------------

def subset_alignment_for_accessions(input_fasta, accessions, out_fasta):
    """Subset an alignment to only the requested sample IDs."""
    records = list(AlignIO.read(input_fasta, "fasta"))
    desired = set(accessions)
    filtered = [rec for rec in records if rec.id in desired]

    missing = sorted(desired - {rec.id for rec in filtered})
    if missing:
        print(f"[WARNING]: Missing records in alignment for accessions: {missing}")

    if not filtered:
        raise ValueError(f"No records found in {input_fasta} for the requested accession list.")

    with open(out_fasta, "w") as handle:
        SeqIO.write(filtered, handle, "fasta")

    print(f"[SUBSET]: Wrote {len(filtered)} sequences to {out_fasta}")
    return out_fasta


def run_snpsites(input_fasta, out_dir, snpsites_bin):
    """Extract variant positions into VCF and SNP-only alignment."""
    print("[SNP-SITES]: Extracting variant positions from the alignment...")
    vcf_out = os.path.join(out_dir, "variants.vcf")
    fasta_out = os.path.join(out_dir, "snp_alignment.fasta")

    cmd_vcf = [snpsites_bin, "-v", "-o", vcf_out, input_fasta]
    subprocess.run(cmd_vcf, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    cmd_fasta = [snpsites_bin, "-m", "-o", fasta_out, input_fasta]
    subprocess.run(cmd_fasta, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    print(
        f"[SNP-SITES]: Variant files generated:\n  - VCF: {vcf_out}\n  - SNP Alignment: {fasta_out}"
    )
    return vcf_out, fasta_out


def compute_pairwise_snp_matrix(alignment_file):
    """Calculate pairwise SNP distances for all sample pairs."""
    print("[SNP MATRIX]: Reading alignment and calculating pairwise SNP distances...")
    alignment = AlignIO.read(alignment_file, "fasta")
    headers = [rec.id for rec in alignment]
    n_seqs = len(headers)

    seq_matrix = np.array([list(str(rec.seq).upper()) for rec in alignment])
    dist_matrix = np.zeros((n_seqs, n_seqs), dtype=int)
    pair_list = []

    canonical_bases = np.array(["A", "C", "G", "T"])

    for i, j in itertools.combinations(range(n_seqs), 2):
        seq1 = seq_matrix[i]
        seq2 = seq_matrix[j]

        is_canonical_1 = np.isin(seq1, canonical_bases)
        is_canonical_2 = np.isin(seq2, canonical_bases)
        comparable_mask = is_canonical_1 & is_canonical_2

        snp_diff = int(np.sum((seq1 != seq2) & comparable_mask))

        dist_matrix[i, j] = snp_diff
        dist_matrix[j, i] = snp_diff

        pair_list.append(
            {
                "Sample_1": headers[i],
                "Sample_2": headers[j],
                "SNP_Distance": snp_diff,
            }
        )

    df_matrix = pd.DataFrame(dist_matrix, index=headers, columns=headers)
    df_pairs = pd.DataFrame(pair_list)

    return df_matrix, df_pairs


def plot_snp_histogram(df_pairs, out_dir, snp_cutoff):
    """Generate a histogram of pairwise SNP distances."""
    print("[PLOTTING]: Generating pairwise SNP distance histogram...")
    plt.figure(figsize=(8, 5))

    max_snps = int(df_pairs[GLOBAL_COL].max()) if not df_pairs.empty else 10
    bins = np.arange(-0.5, max_snps + 1.5, 1)

    sns.histplot(
        df_pairs[GLOBAL_COL],
        bins=bins,
        color="#2b5c8f",
        edgecolor="black",
        alpha=0.8,
    )

    plt.axvline(
        x=snp_cutoff,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Cluster Cutoff ({snp_cutoff} SNPs)",
    )

    plt.title("Distribution of Pairwise SNP Differences", fontsize=12, fontweight="bold")
    plt.xlabel("Pairwise SNP Difference", fontsize=10)
    plt.ylabel("Number of Sample Pairs", fontsize=10)
    plt.grid(axis="y", linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    fig_path = os.path.join(out_dir, "pairwise_snp_histogram.png")
    plt.savefig(fig_path, dpi=300)
    plt.close()
    print(f"[PLOTTING]: Histogram saved to {fig_path}")


def perform_hierarchical_clustering(df_matrix, out_dir, snp_cutoff, linkage_method):
    """Cluster samples using average-linkage distance clustering."""
    print(f"[CLUSTERING]: Conducting hierarchical clustering ({linkage_method} linkage)...")

    condensed_dist = squareform(df_matrix.values, force="tovector")
    Z = linkage(condensed_dist, method=linkage_method)

    plt.figure(figsize=(10, 6))
    dendrogram(
        Z,
        labels=df_matrix.index,
        leaf_rotation=90,
        leaf_font_size=10,
        color_threshold=snp_cutoff,
    )

    plt.axhline(
        y=snp_cutoff,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Cutoff = {snp_cutoff} SNPs",
    )

    plt.title(
        f"Hierarchical Clustering Dendrogram ({linkage_method.capitalize()} Linkage)",
        fontsize=12,
        fontweight="bold",
    )
    plt.xlabel("Sample Accession", fontsize=10)
    plt.ylabel("Pairwise SNP Distance", fontsize=10)
    plt.legend(loc="upper right")
    plt.tight_layout()

    dendrogram_path = os.path.join(out_dir, "hierarchical_clustering_dendrogram.png")
    plt.savefig(dendrogram_path, dpi=300)
    plt.close()
    print(f"[PLOTTING]: Dendrogram saved to {dendrogram_path}")

    cluster_ids = fcluster(Z, t=snp_cutoff, criterion="distance")
    df_clusters = pd.DataFrame(
        {
            "Sample_ID": df_matrix.index,
            "Cluster_ID": [f"Cluster_{cid}" for cid in cluster_ids],
        }
    )
    return df_clusters


def write_analysis_metadata(output_dir, samples, analysis_name, gene_context=None, fallback_reason=None):
    metadata_path = os.path.join(output_dir, "analysis_metadata.json")
    payload = {
        "analysis_name": analysis_name,
        "input_samples": samples,
        "module": "module_8_snpsites_hcluster",
        "global_column": GLOBAL_COL,
        "global_distance_source": "module 7 trimmed samples-only alignment (sample-sample pairs)",
    }
    if gene_context:
        payload.update({
            "gene_aware_analysis": True,
            "msa_reference": gene_context["msa_reference"],
            "reference_id_in_alignment": gene_context["ref_id"],
            "gff_path": gene_context["gff_path"],
            "gene_counting_alignment": gene_context["alignment_path"],
            "gene_names": [g["name"] for g in gene_context["genes"]],
            "gene_intervals": {g["name"]: [g["start"], g["end"], g["strand"]] for g in gene_context["genes"]},
            "reference_pairs_included": True,
            "reference_pair_global_source": "ref-included untrimmed alignment, all comparable columns",
            "gene_msa_files": gene_context["gene_msa_files"],
        })
    else:
        payload.update({"gene_aware_analysis": False, "fallback_reason": fallback_reason})
    with open(metadata_path, "w") as handle:
        json.dump(payload, handle, indent=4)
    print(f"[METADATA]: Wrote {metadata_path}")


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run SNP-site extraction and hierarchical clustering for the module 7 alignment."
    )
    parser.add_argument("--config", required=True, help="Path to the config Python file.")
    args = parser.parse_args()

    config = load_config(args.config)

    input_samples = get_input_samples(config)
    analysis_name = (
        get_config_value(config, "module_8_params", "analysis_name", None)
        or get_config_value(config, "module_7_params", "analysis_name", "default")
        or "default"
    )
    snp_cutoff = float(
        get_config_value(config, "module_8_params", "snp_cutoff", 3)
    )
    linkage_method = (
        get_config_value(config, "module_8_params", "linkage_method", "average")
        or "average"
    )
    snpsites_bin = resolve_executable_path(config, "snpsites_bin_path")

    output_dir = get_module_8_dir(config)
    os.makedirs(output_dir, exist_ok=True)

    print("\n==========================================")
    print(f" Step 8: Pairwise SNP & Hierarchical Clustering Analysis for '{analysis_name}'")
    print("==========================================")
    print(f"[JOB]: analysis_name = {analysis_name}")
    print(f"[JOB]: input_samples = {input_samples}")
    print(f"[JOB]: output_dir = {output_dir}")

    input_alignment = find_trimmed_alignment(config)
    print(f"[INFO]: Using module 7 trimmed alignment: {input_alignment}")

    subset_alignment = os.path.join(output_dir, "subset_accessions_alignment.fasta")
    subset_alignment = subset_alignment_for_accessions(
        input_alignment,
        input_samples,
        subset_alignment,
    )

    vcf_file, snp_fasta = run_snpsites(subset_alignment, output_dir, snpsites_bin)
    del vcf_file, snp_fasta

    df_matrix, df_pairs = compute_pairwise_snp_matrix(subset_alignment)
    df_pairs = df_pairs.rename(columns={"SNP_Distance": GLOBAL_COL})

    matrix_csv = os.path.join(output_dir, "pairwise_snp_matrix.csv")
    df_matrix.to_csv(matrix_csv)

    # Histogram + clustering stay samples-only (reference excluded by design)
    plot_snp_histogram(df_pairs, output_dir, int(snp_cutoff))

    df_clusters = perform_hierarchical_clustering(
        df_matrix,
        output_dir,
        int(snp_cutoff),
        linkage_method,
    )

    cluster_csv = os.path.join(output_dir, f"clusters_{int(snp_cutoff)}snp_cutoff.csv")
    df_clusters.to_csv(cluster_csv, index=False)

    # ------------------------------------------------------------------
    # Gene-aware extension (optional; falls back silently to legacy mode)
    # ------------------------------------------------------------------
    gene_context = build_gene_analysis(config, output_dir, input_samples)
    fallback_reason = None
    if gene_context is None:
        fallback_reason = "msa_reference undefined or gene-aware prerequisites unavailable"

    if gene_context:
        genes = gene_context["genes"]
        gene_cols = [f"{g['name']} SNP Difference" for g in genes]

        # Build numeric (float) gene columns in one pass; np.nan marks pairs that
        # could not be compared in the ref-included alignment.
        col_values = {col: [] for col in gene_cols}
        for _, row in df_pairs.iterrows():
            counts = gene_context["gene_counts"].get((row["Sample_1"], row["Sample_2"]))
            for gene in genes:
                col = f"{gene['name']} SNP Difference"
                if counts and gene["name"] in counts:
                    col_values[col].append(float(counts[gene["name"]][0]))
                else:
                    col_values[col].append(np.nan)
        for col in gene_cols:
            df_pairs[col] = col_values[col]

        ref_id = gene_context["ref_id"]
        new_rows = []
        for sample in gene_context["wanted_samples"]:
            row = {
                "Sample_1": sample,
                "Sample_2": ref_id,
                GLOBAL_COL: gene_context["ref_global"][(sample, ref_id)],
            }
            counts = gene_context["gene_counts"].get((sample, ref_id), {})
            for gene in genes:
                value = counts.get(gene["name"])
                row[f"{gene['name']} SNP Difference"] = float(value[0]) if value else np.nan
            new_rows.append(row)
        if new_rows:
            df_pairs = pd.concat([df_pairs, pd.DataFrame(new_rows)], ignore_index=True)

        df_pairs = df_pairs[["Sample_1", "Sample_2", GLOBAL_COL] + gene_cols]

        # Enforce integer output: nullable Int64 writes real counts without a
        # decimal point while still representing missing pairs as NA.
        snp_cols = [GLOBAL_COL] + gene_cols
        df_pairs[snp_cols] = df_pairs[snp_cols].astype("Int64")

        df_long = pd.DataFrame(gene_context["long_rows"])
        long_csv = os.path.join(output_dir, "pairwise_snp_pairs_per_gene_long.csv")
        df_long.to_csv(long_csv, index=False)
        print(f"[SUCCESS]: Saved per-gene long-format distances to {long_csv}")
    else:
        df_pairs = df_pairs[["Sample_1", "Sample_2", GLOBAL_COL]]
        df_pairs[GLOBAL_COL] = df_pairs[GLOBAL_COL].astype("Int64")

    pairs_csv = os.path.join(output_dir, "pairwise_snp_pairs.csv")
    df_pairs.to_csv(pairs_csv, index=False, na_rep="NA")
    print(f"[SUCCESS]: Saved SNP matrix ({matrix_csv}) and pairs list ({pairs_csv})")

    write_analysis_metadata(output_dir, input_samples, analysis_name, gene_context, fallback_reason)

    print("\n==========================================")
    print(" Cluster Summary Results")
    print("==========================================")
    print(df_clusters.to_string(index=False))
    print(f"\n[COMPLETE]: All clustering outputs saved to {output_dir}")


if __name__ == "__main__":
    main()
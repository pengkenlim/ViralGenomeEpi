# Module 6B — Direct Reference-Mapping Consensus

**Script:** [`pipeline/operational/module_6b_direct_consesnsus.py`](../../pipeline/operational/module_6b_direct_consesnsus.py)

## Purpose

Module 6B generates consensus genomes **without *de novo* assembly**, mimicking the CZ ID-style approach: trimmed reads are mapped directly onto a configured reference genome, and the consensus is called from the pileup. This pathway is the fallback of choice for low-titer or highly divergent samples whose assemblies fail in Modules 3–5, and for amplicon data where a reference is already known. It additionally computes detailed depth/coverage metrics and runs CheckV when the genome is complete enough to support it.

## Input Requirements

- Trimmed paired-end FASTQs from Module 1 (`<output_dir>/module_1/fastp/<sample>/`); singletons are appended to the alignment when present.
- A per-sample reference map in `module_6b_params["reference_genomes"]` — typically pointing at Module 4's downloaded genomes (`<output_dir>/module_4/ref_genomes/<GCF>.fasta`). **Samples absent from this map are skipped.**

```python
module_6b_params = {
    "threads": 8,
    "minimap2_bin_path": "minimap2",
    "samtools_bin_path": "samtools",
    "ivar_bin_path": "ivar",
    "checkv_bin_path": "checkv",
    "checkv_db_dir": "/path/to/your/checkv/db",   # optional
    "ivar_min_quality": 20,
    "ivar_min_freq": 0.5,
    "ivar_min_depth": 10,
    "reference_genomes": {
        "SRR34884179": "/abs/path/module_4/ref_genomes/GCF_000855545.1.fasta",
    },
}
```

## Core Algorithm / Logic

1. **Mapping** — `minimap2 -a -x sr <ref> <R1> <R2> [singleton] | samtools sort` → coordinate-sorted BAM, then `samtools index`.
2. **Consensus calling** — the classic mpileup-to-iVar pipeline:

```bash
samtools mpileup -aa -A -d 0 -Q 0 -f <ref.fa> <sample>_direct_consensus.bam \
  | ivar consensus -p <sample>_direct_consensus -q 20 -t 0.5 -m 10
```

   `-aa -A -d 0 -Q 0` emit all positions (including zero-coverage and orphan reads) with no depth or base-quality caps, leaving thresholding entirely to iVar: positions below `ivar_min_depth` become `N`; elsewhere the majority allele (≥ `ivar_min_freq`, quality ≥ `ivar_min_quality`) is called. iVar's `.fa` output is normalized to `.fasta` for consistency.
3. **Sequence evaluation** — `length`, `n_count`, `n_pct`, and `non_canonical_count` (bases outside A/C/G/T/N) of the consensus.
4. **Depth metrics** — computed from `samtools depth -a` over the sorted BAM:

| Column | Definition |
|---|---|
| `coverage_1x_pct` | % of reference bases with depth ≥ 1 |
| `coverage_10x_pct` | % of reference bases with depth ≥ 10 (the `ivar_min_depth` default) |
| `mean_depth_overall` | Mean depth over **all** reference bases (zeros included) |
| `mean_depth_covered` | Mean depth over covered bases only |
| `min_depth_covered` | Shallowest covered position |
| `max_depth` | Deepest position |

5. **CheckV gating** — CheckV's gene-annotation stage (HMMER `hmmsearch`) can crash on sequences that are mostly `N`. The module therefore **skips CheckV when `n_pct > 50`** and records `checkv_warnings = "Skipped: >50% Ns"`; otherwise CheckV runs and `checkv_completeness` / `checkv_quality` / `checkv_kmer_freq` / `checkv_warnings` are captured. A CheckV failure degrades to `N/A` fields plus the captured stderr, never a crash.

## Visualization (IGV)

The indexed BAM is designed for direct inspection in the **Integrative Genomics Viewer (IGV)**:

```bash
# Files (per sample), under <output_dir>/module_6b/<sample>/
#   <sample>_direct_consensus.bam
#   <sample>_direct_consensus.bam.bai
```

1. Open IGV and load the reference genome used for mapping: *Genomes → Load Genome from File…* → the sample's `reference_fasta`.
2. Load the alignment: *File → Load from File…* → `<sample>_direct_consensus.bam` (the `.bai` index must sit next to it, which the module guarantees).
3. Inspect the **coverage track** for dropouts: regions with zero or &lt;10× coverage correspond to `N` stretches in `<sample>_direct_consensus.fasta` and to depressed `coverage_10x_pct` in the summary. Comparing dropout coordinates against `checkv_warnings` helps distinguish genuine low-titer samples from mapping to a mismatched reference.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_6b_params["reference_genomes"]` | `{}` | **Required.** `{sample: /path/to/reference.fasta}`; unlisted samples are skipped |
| `module_6b_params["minimap2_bin_path"]` | — | Minimap2 executable (explicit path or `PATH`-resolvable name) |
| `module_6b_params["samtools_bin_path"]` | — | samtools executable |
| `module_6b_params["ivar_bin_path"]` | — | iVar executable |
| `module_6b_params["checkv_bin_path"]` | — | CheckV executable |
| `module_6b_params["checkv_db_dir"]` | `None` | CheckV DB (optional; used only if the directory exists) |
| `module_6b_params["ivar_min_quality"]` | `20` | iVar `-q` minimum base quality |
| `module_6b_params["ivar_min_freq"]` | `0.5` | iVar `-t` minimum allele frequency |
| `module_6b_params["ivar_min_depth"]` | `10` | iVar `-m` minimum depth to call a base (below → `N`) |
| `module_6b_params["threads"]` | `4` | Mapping threads |

## Output Files

Under `<output_dir>/module_6b/`:

| Path | Description |
|---|---|
| `<sample>/<sample>_direct_consensus.fasta` | Consensus genome (Ns at low-depth positions); input to Module 7 when `consensus_source = "module_6b"` |
| `<sample>/<sample>_direct_consensus.bam` / `.bam.bai` | Sorted + indexed alignment for IGV inspection |
| `<sample>/checkv/` | CheckV output (when run) |
| `direct_consensus_summary.tsv` | Cohort summary: `sample_accession`, `reference_fasta`, `consensus_fasta`, `length`, `n_count`, `n_pct`, `non_canonical_count`, `coverage_1x_pct`, `coverage_10x_pct`, `mean_depth_overall`, `mean_depth_covered`, `min_depth_covered`, `max_depth`, `checkv_completeness`, `checkv_quality`, `checkv_kmer_freq`, `checkv_warnings` |

## Error Handling &amp; Edge Cases

- **Sample not in `reference_genomes`** → `[SKIP]` log line; no output.
- **Missing reference / Module 1 directory / trimmed reads** → `FileNotFoundError` at validation or per-sample processing.
- **Unresolvable executable** → `FileNotFoundError` naming the config key (`resolve_executable_path` also searches `PATH`).
- **iVar produces no consensus** → `RuntimeError` with iVar's captured stderr.
- **`n_pct > 50`** → CheckV skipped by design (HMMER crash prevention); treated as a quality signal, not an error.
- **Empty depth output** → all depth metrics reported as `0` rather than raising.

## Usage

```bash
python pipeline/operational/module_6b_direct_consesnsus.py --config my_config.py
```

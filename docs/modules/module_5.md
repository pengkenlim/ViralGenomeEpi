# Module 5 — Reference-Guided Scaffolding &amp; Quality Control

**Script:** [`pipeline/operational/module_5_scaffold_qc.py`](../../pipeline/operational/module_5_scaffold_qc.py)

## Purpose

Module 5 converts a sample's fragmented *de novo* contigs into a full-length, reference-ordered scaffold. It aligns the assembly to the Module 4 reference genome with Minimap2 and applies a **6-tier hybrid fallback** that combines single-contig **slicing** and multi-contig **stitching** across three alignment stringency modes, accepting the first reconstruction that meets strict coverage and N-content thresholds. Every accepted scaffold is then quality-assessed with **CheckV**. Segmented viruses are handled per segment.

## Input Requirements

- `<output_dir>/module_3/<sample>/combined_transcripts.fasta` — assembled contigs (written/verified by Module 4).
- `<output_dir>/module_4/<sample>/<sample>_reference_summary.json` — the reference selection record, from which the reference FASTA in `module_4/ref_genomes/` is resolved.
- For segmented samples, the summary's `selected_segments[*].fasta_path` entries must exist.
- **Overrides** (checked first, bypassing Module 4 resolution):
  - `module_5_params["override_references"][sample] = "/abs/path/to/ref.fasta"` — full FASTA-path override (validated for existence), and/or
  - `module_5_params["explicit_ref_genome"][sample] = "NC_045512.2"` — accession-level override, and
  - `module_5_params["exclude_contigs"][sample] = ["contig_7", ...]` — contig IDs removed before alignment/scaffolding.

```python
module_5_params = {
    "minimap2_bin_path": "/PATH/TO/minimap2",
    "checkv_bin_path": "checkv",
    "checkv_db_dir": "/PATH/TO/checkv-db-v1.5",
    "threads": 64,
    "fallback_strategy": "technique_first",     # default; or "stringency_first"
    "exclude_contigs": {"sampleA": ["contigname_1", "contigname_2"]},
    "explicit_ref_genome": {"sampleA": "NC_045512.2"},
}
```

## Core Algorithm / Logic

### Slice vs. stitch

- **Slice** — if a *single* contig aligns across **≥ 98% of the reference**, the aligned block is excised from the contig (reverse-complemented when on the minus strand) and padded with `N`s at the unaligned reference termini, yielding a reference-length scaffold. Accepted when the result contains **≤ 5% Ns**.
- **Stitch** — when no single contig suffices, non-redundant alignments are sorted by reference coordinate and merged into one scaffold, with `N` runs filling inter-contig gaps and terminal ends. Accepted when total reference coverage is **≥ 98%** and the result contains **≤ 2% Ns** (stricter, because stitching introduces more junctions).

### The 6-tier fallback loop

Three Minimap2 alignment modes are crossed with the two techniques:

| Mode | Minimap2 arguments | Minimum alignment identity |
|---|---|---|
| `conserved` | `-c -x asm5` | 70% |
| `divergent_70` | `-c -k 14 -w 5` | 70% |
| `divergent_60` | `-c -k 14 -w 5` | 60% |

`fallback_strategy` controls the traversal order of the resulting six `(mode, technique)` tiers:

- `technique_first` (**default**): slice@conserved → slice@divergent_70 → slice@divergent_60 → stitch@conserved → stitch@divergent_70 → stitch@divergent_60. Prefers single-contig scaffolds (fewer junctions) before relaxing alignment stringency.
- `stringency_first`: slice → stitch per mode before moving to the next mode. Prefers the highest-confidence alignments even if stitching is required.

The **first tier whose result passes its thresholds wins** and is recorded as `alignment_mode` in `Identified_scaffolds.tsv`. Rejected candidate scaffolds are deleted immediately to keep the output directory unambiguous, and PAF files are cached per mode (each of the three modes runs Minimap2 at most once) and cleaned up on exit.

### Identity-weighted overlap trimming (stitch mode)

Adjacent contig alignments frequently overlap on the reference; naively concatenating them would **duplicate sequence** at the junctions. Module 5 resolves overlaps by trimming the *lower-identity* alignment:

1. Alignments are filtered to the mode's identity threshold and **contained (nested) intervals are discarded** — a contig block fully inside another contributes nothing new.
2. Walking left-to-right along the reference, when the current alignment overlaps the previously accepted one, the alignment with the **lower percent identity is trimmed by the exact overlap length** (adjusting both reference and query coordinates, strand-aware); if trimming empties it, it is dropped. Equal identity keeps the earlier (leftmost) alignment and trims the later one.
3. Only then are the clipped contig subsequences joined with gap-filling `N`s, so every reference base is contributed by exactly one contig — the highest-identity donor available.

For stitched scaffolds, a per-contig coordinate table (`<sample>_stitch_contig_locations.tsv`) records which reference interval each contig contributed, enabling full provenance tracing.

### CheckV QC

Each accepted scaffold is assessed with `checkv end_to_end ... --restart`. To survive `hmmsearch` crashes at high parallelism, the module retries with reduced thread counts (`8 → threads/2 → threads`). Aggregated metrics — mean `completeness`, best `checkv_quality` (ranked High &gt; Medium &gt; Low &gt; Not-determined), mean `kmer_freq`, and concatenated `warnings` — feed both the per-sample report and the cohort TSV. **Segmented samples** are additionally evaluated both per segment (`checkv_segments/`) and as an `NNNNN`-joined concatenated genome (`checkv_concatenated/`) for an overall completeness estimate.

### Failure logging

When all six tiers fail, the module does not raise: it records every tier's failure reason (with observed coverage and N% where available) into the `failed_attempt_metrics` and `failure_reason` columns of `Identified_scaffolds.tsv`, marks `scaffold_fasta` as `FAILED`, and logs a consolidated per-sample warning — the effective `_scaffold_failure` record for the sample.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_5_params["minimap2_bin_path"]` | — | Minimap2 executable (≥ 2.x) |
| `module_5_params["checkv_bin_path"]` | `"checkv"` | CheckV executable |
| `module_5_params["checkv_db_dir"]` | — | CheckV database directory |
| `module_5_params["fallback_strategy"]` | `"technique_first"` | Tier traversal: `technique_first` or `stringency_first` |
| `module_5_params["exclude_contigs"]` | `{}` | Per-sample contig IDs dropped before scaffolding |
| `module_5_params["explicit_ref_genome"]` | `{}` | Per-sample accession override for the scaffold reference |
| `module_5_params["override_references"]` | `{}` | Per-sample absolute reference FASTA path override |
| `module_5_params["threads"]` | `8` | Mapping/QC threads |

Thresholds are constants in the script: `SLICE_MIN_REF_COVERAGE = 0.98`, `SLICE_MAX_N_PCT = 5.0`, `STITCH_MIN_REF_COVERAGE = 0.98`, `STITCH_MAX_N_PCT = 2.0`.

## Output Files

Under `<output_dir>/module_5/`:

| Path | Description |
|---|---|
| `<sample>/<sample>_scaffold.fasta` | Non-segmented scaffold (reference length, N-padded) |
| `<sample>/<sample>_seg_<id>_scaffold.fasta` | Per-segment scaffolds (segmented viruses) |
| `<sample>/<sample>_segmented_scaffold.fasta` | Multi-record FASTA of all segments |
| `<sample>/<sample>_concatenated_scaffold.fasta` | `NNNNN`-joined concatenated genome (completeness estimation) |
| `<sample>/<sample>_stitch_contig_locations.tsv` | Reference coordinates contributed by each stitched contig |
| `<sample>/<sample>_scaffold_checkv_report.txt` | Human-readable QC report (type, mode, N stats, CheckV metrics) |
| `<sample>/checkv/`, `checkv_segments/`, `checkv_concatenated/` | Raw CheckV outputs |
| `Identified_scaffolds.tsv` | Cohort summary: `sample_accession`, `reference_fasta`, `scaffold_fasta`, `reference_length`, `scaffold_length`, `checkv_completeness`, `checkv_quality`, `checkv_kmer_freq`, `checkv_warnings`, `scaffold_type`, `alignment_mode`, `gap_filled_bases`, `failed_attempt_metrics`, `failure_reason` |

## Error Handling &amp; Edge Cases

- **All tiers fail** → sample recorded as `FAILED` with full per-tier reasons in `Identified_scaffolds.tsv`; the batch continues.
- **Missing Module 4 summary / reference / assembly** → `FileNotFoundError` at validation (override samples skip Module 4 validation but their FASTA path is checked).
- **Segmented sample with missing segment FASTA** → validation error naming the missing path.
- **CheckV `hmmsearch` failure** → automatic retry with fewer threads; persistent failure raises only after all candidates are exhausted, or returns `N/A` fields when HMMER itself is the failure point.
- **Override samples are assumed non-segmented** — do not use `override_references` for multi-segment genomes.

## Usage

```bash
python pipeline/operational/module_5_scaffold_qc.py --config my_config.py
```

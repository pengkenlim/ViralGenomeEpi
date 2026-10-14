# Module 6 — Consensus Polishing (Two-Tier iVar)

**Script:** [`pipeline/operational/module_6_consensus_polish.py`](../../pipeline/operational/module_6_consensus_polish.py)

## Purpose

Module 6 turns Module 5 scaffolds into final polished consensus genomes. Trimmed reads are re-mapped to each scaffold, and iVar calls the consensus in **two passes** — a strict pass for high-confidence SNP positions and a relaxed pass to rescue low-coverage termini and remaining N-gaps — which are then merged. Each polished genome is QC'd with CheckV.

## Input Requirements

- Module 5 scaffolds: `<output_dir>/module_5/<sample>/<sample>_scaffold.fasta` (or per-segment `<sample>_seg_<id>_scaffold.fasta` files, resolved via the Module 4 summary's segmentation status).
- Trimmed paired reads from Module 1 (`<output_dir>/module_1/fastp/<sample>/`); falls back to the raw input directory when necessary.
- Executables: `minimap2`, `samtools`, `ivar`, `checkv` (+ CheckV DB).

## Core Algorithm / Logic

1. **Re-mapping** — `minimap2 -ax sr <scaffold> <R1> <R2> | samtools view -bS | samtools sort`, followed by `samtools index`.
2. **Two-tier iVar consensus** (`samtools mpileup` → `ivar consensus`):

| Tier | Min quality (`-q`) | Min frequency (`-t`) | Min depth (`-m`) | Role |
|---|---|---|---|---|
| Strict | 20 | 0.60 | 10 | High-confidence core-genome SNP calls |
| Relaxed | 15 | 0.51 | 3 | Rescue of low-coverage terminal UTRs and gap regions |

   The relaxed frequency of 0.51 is deliberate: it prevents arbitrary base calls on exact 50/50 allelic ties (e.g. mixed infections), which would be irreproducible.
3. **Merging** — the strict consensus is the backbone; positions where the strict pass emitted `N` but the relaxed pass produced a canonical base are rescued from the relaxed sequence. Per-base transition statistics (`unpolished_n`, `polished_n`, `n_to_canonical`, `canonical_to_n`, `net_n_change`) quantify gap-filling versus coverage loss.
4. **Segmentation handling** — segments are polished independently, then combined into `<sample>_consensus_combined.fasta`; CheckV runs on the *combined* genome so that segmented viruses are not penalized as "no viral genes detected".
5. **CheckV QC** — same crash-resilient retry ladder as Module 5 (`8 → threads/2 → threads`), reporting `completeness`, `quality`, `kmer_freq`, `warnings`.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_6_params["minimap2_bin_path"]` | — | Minimap2 executable |
| `module_6_params["samtools_bin_path"]` | — | samtools executable |
| `module_6_params["ivar_bin_path"]` | — | iVar executable |
| `module_6_params["checkv_bin_path"]` | `"checkv"` | CheckV executable |
| `module_6_params["checkv_db_dir"]` | falls back to Module 5 value | CheckV database directory |
| `module_6_params["threads"]` | falls back to Module 5/`cross_module_params` | Mapping and QC threads |

The iVar thresholds themselves are hardcoded constants (`STRICT_PARAMS`, `RELAXED_PARAMS`) to keep the two-tier rescue strategy consistent across runs.

## Output Files

Under `<output_dir>/module_6/`:

| Path | Description |
|---|---|
| `<sample>/<file_label>_aligned_sorted.bam` (+ `.bai`) | Sorted, indexed read-to-scaffold alignment |
| `<sample>/<file_label>_final_polished_consensus.fa` | Merged strict+relaxed polished consensus (per scaffold/segment) |
| `<sample>/<file_label>_ivar_polished_report.txt` | Human-readable report: N-transition analysis + CheckV metrics |
| `<sample>/<sample>_consensus_combined.fasta` | **Primary deliverable** — final polished consensus genome(s); input to Module 7 when `consensus_source = "module_6"` |
| `<sample>/checkv_combined/` | CheckV output for the combined genome |
| `Identified_consensus.tsv` | Cohort summary: `sample_accession`, `is_segmented`, `segment_id`, `input_scaffold_fasta`, `polished_fasta`, `unpolished_n_bases`, `polished_n_bases`, `n_rescued_to_canonical`, `canonical_lost_to_n`, `net_n_change`, `checkv_completeness`, `checkv_quality`, `checkv_kmer_freq`, `checkv_warnings`, `report_path` |

(`file_label` is `<sample>` for non-segmented genomes, `<sample>_seg_<id>` per segment otherwise.)

## Error Handling &amp; Edge Cases

- **Missing Module 4 summary or trimmed reads** → `FileNotFoundError`.
- **Alignment/sort pipeline failure** → `RuntimeError` with the captured `samtools sort` stderr.
- **iVar fails to emit a consensus** → `RuntimeError` naming the sample.
- **CheckV field sanitization** — all CheckV text fields are stripped of embedded newlines/tabs so `Identified_consensus.tsv` can never be corrupted by malformed tool output.
- **High-N genomes** are *not* skipped here (unlike Module 6B); if CheckV itself fails, fields are recorded as `N/A` with a warning rather than crashing the module.

## Usage

```bash
python pipeline/operational/module_6_consensus_polish.py --config my_config.py
```

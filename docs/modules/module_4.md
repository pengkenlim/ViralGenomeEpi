# Module 4 — Reference Genome Selection &amp; Download

**Script:** [`pipeline/operational/module_4_download_ref_genome.py`](../../pipeline/operational/module_4_download_ref_genome.py)

## Purpose

Module 4 identifies the most appropriate reference genome for each sample's *de novo* assembly and downloads it from NCBI. It BLASTs the assembled transcripts against a local viral RefSeq database, clusters redundant hits, aggregates evidence per genome, resolves full taxonomy, detects segmented viruses, and fetches the chosen genome FASTA — producing the reference upon which Module 5 scaffolding and Module 6B direct mapping depend.

## Input Requirements

- `<output_dir>/module_3/<sample>/combined_transcripts.fasta` — the combined assembly FASTA. If absent, the module itself merges any `transcripts.fasta` found across Module 3 modes (`viral`, `viral_over_assembly`, `dominant_virus`, `dominant_virus_over_assembly`, `aligned_to_ref`, `aligned_to_ref_over_assembly`), deduplicating repeated FASTA identifiers with numeric suffixes.
- Local reference data:
  - Viral RefSeq **BLAST database** (`blast_db_path`, validated via the presence of `<db>.nin`),
  - **`seq2genome_map.json`** mapping sequence accessions → genome/assembly accessions (built by [`setup/prepare_database.py`](../../setup/prepare_database.py)),
  - **`assembly_data_report.jsonl`** (NCBI assembly metadata),
  - NCBI taxonomy dumps `nodes.dmp` / `names.dmp`.
- Network access to the **NCBI Datasets API** for genome download (existing local copies are reused).

## Core Algorithm / Logic

1. **BLAST search** — `blastn` of `combined_transcripts.fasta` against the viral RefSeq DB (`-max_target_seqs 20 -evalue 1e-5 -max_hsps 1`, tabular output `qseqid sseqid pident length evalue bitscore`).
2. **Hit filtering** — only hits with **alignment length ≥ 500 bp and identity ≥ 70%** contribute.
3. **Best-reference selection (non-segmented path)** — bitscores of passing hits are aggregated per genome accession (via `seq2genome_map.json`, with basename fallback for version mismatches); the genome with the highest total bitscore wins.
4. **Segment detection** — the candidate genome is downloaded and inspected: **>1 sequence ⇒ segmented virus**. For segmented genomes, per-segment representatives are chosen by **cd-hit-est clustering (80% identity)** of the passing BLAST hits: within each homology cluster, the member with the highest aggregated bitscore becomes the segment representative, and each is mapped to its genome accession and segment sequence ID.
5. **Taxonomy resolution** — the organism taxid from `assembly_data_report.jsonl` is walked up `nodes.dmp` to the species rank; scientific/common names come from `names.dmp`.
6. **Download &amp; symlinking** — genome FASTAs are stored centrally in `module_4/ref_genomes/<GCF_accession>.fasta`; per-sample `selected_reference_<acc>.fasta` symlinks point into that store. For segmented samples, per-segment FASTAs are extracted.
7. **Reporting** — a per-sample JSON summary (`<sample>_reference_summary.json`) records the genome accession, segmentation status, selected segments, and taxonomy; a cohort table `Identified_species.tsv` aggregates the results.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_4_params["blastn_bin_path"]` | `"blastn"` | BLASTN executable |
| `module_4_params["blastdbcmd_bin_path"]` | `"blastdbcmd"` | Sequence extraction from the local BLAST DB |
| `module_4_params["cd_hit_est_bin_path"]` | `"cd-hit-est"` | Clustering tool for segment resolution |
| `module_4_params["blast_db_path"]` | local fallback path | Viral RefSeq BLAST DB prefix — **set explicitly** |
| `module_4_params["seq2genome_map_path"]` | local fallback path | Sequence→genome accession map — **set explicitly** |
| `module_4_params["assembly_data_report_path"]` | local fallback path | NCBI assembly metadata JSONL — **set explicitly** |
| `module_4_params["nodes_path"]` / `["names_path"]` | local fallback paths | NCBI taxonomy dumps — **set explicitly** |
| `module_4_params["threads"]` | `8` | BLAST threads |

> **Warning:** several defaults are hardcoded to the original authors' machine (`/home/mngs/...`). Always configure these paths explicitly in your config file.

## Output Files

Under `<output_dir>/module_4/`:

| Path | Description |
|---|---|
| `<sample>/<sample>_blastn_viral_refseq.tsv` | Raw BLAST hits of the assembly |
| `<sample>/<sample>_reference_summary.json` | Selection record: `genome_accession` / `best_reference_accession`, `is_segmented`, `selected_segments` (with `segment_seq_id` and `fasta_path`), taxonomy fields |
| `<sample>/selected_reference_<acc>.fasta` | Symlink to the chosen reference |
| `ref_genomes/<GCF_accession>.fasta` | Central store of downloaded genomes (shared across samples and reused by Modules 5, 6B, 7) |
| `Identified_species.tsv` | Cohort-level table of selected references, species names, and segmentation status |

## Error Handling &amp; Edge Cases

- **Missing BLAST DB (`.nin`)** or missing `combined_transcripts.fasta` → `FileNotFoundError` before any search.
- **No passing BLAST hits** (≥500 bp, ≥70% identity) → `ValueError("No high-quality matches found...")`; the sample cannot proceed to scaffolding — review assembly quality or supply `module_5_params["explicit_ref_genome"]`.
- **Ambiguous accession mapping** — falls back to basename (`NC_xxxxxx`) matching within `seq2genome_map.json`; raises `KeyError` if still unresolved.
- **Download failure / tiny ZIP** → logged, temp files cleaned, and the function returns `False` rather than leaving corrupt references.
- **Segmented samples missing segment FASTAs** are caught at Module 5 validation.

## Usage

```bash
python pipeline/operational/module_4_download_ref_genome.py --config my_config.py
```

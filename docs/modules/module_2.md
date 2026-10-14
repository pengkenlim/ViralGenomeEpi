# Module 2 — Taxonomic Classification &amp; Read Scrubbing

**Script:** [`pipeline/operational/module_2_classify_scrub.py`](../../pipeline/operational/module_2_classify_scrub.py)

## Purpose

Module 2 taxonomically classifies the trimmed reads from Module 1 with **Kraken2** (`k2 classify`) and splits them into analytically useful, non-overlapping read sets: non-viral, viral-only, viral + unclassified, and dominant-virus + unclassified. This both scrubs host (human) reads out of the analysis and concentrates downstream assembly on the viral fraction of the metagenome.

## Input Requirements

- Trimmed FASTQs from Module 1 at `<output_dir>/module_1/fastp/<sample>/` (paired-end required; singletons processed separately when present).
- A Kraken2 database (`kraken2_db_path`).
- NCBI taxonomy dumps `nodes.dmp` (and optionally `names.dmp`) to resolve the exact descendant taxid set of the configured scope; when `nodes.dmp` is unavailable, the module falls back to inferring scope membership from Kraken2 report indentation depth.

## Core Algorithm / Logic

1. **Classification** — paired mates are classified jointly (`k2 classify --paired`), singletons separately; both produce `.kraken` (per-read assignment) and `.report` (taxon summary) files.
2. **Scope resolution** — `scope_to_keep_taxId` (default **10239**, *Viruses*) is expanded to the complete set of descendant taxids using `nodes.dmp`, so membership testing reflects true taxonomy rather than report indentation.
3. **Read-set partitioning** —
   - `viral_plus_unclassified`: reads classified within the viral scope **plus all unclassified (`U`) reads** (broad, assembly-oriented set).
   - `viral`: reads classified within the viral scope only (strict set).
   - `non_viral`: everything else (including human, taxid 9606).
   - `dominant_virus_plus_unclassified`: unclassified reads plus reads in the **dominant viral species lineage** (the species-rank taxon with the largest covered-read count beneath the viral root), intersected with the viral scope to prevent unbounded lineage expansion.
4. **Paired-end read ID normalization** — Kraken2 reports paired-end assignments using the read's shared identifier, while FASTQ headers may carry `/1` and `/2` mate suffixes. `normalize_read_id()` strips these suffixes (`read_id.strip().rstrip("/1").rstrip("/2")`) before any dictionary insertion or set-membership test. This is critical: without normalization, the two mates of a pair would occupy *different* dictionary keys, causing one mate's classification to **overwrite** the other's and corrupting both the taxonomic counts and the paired FASTQ filtering (a pair is written only when *both* normalized IDs are in the keep set). Summary counts in `read_classification_report.tsv` therefore reflect unique read pairs, not doubled mates.
5. **Reporting** — per-sample counts are emitted for human, all-viral, dominant-viral-species, non-viral/non-human, and unclassified reads.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_2_params["Kraken2_k2_bin_path"]` | — | Path to the Kraken2 `k2` executable |
| `module_2_params["kraken2_db_path"]` | — | Kraken2 database directory |
| `module_2_params["scope_to_keep_taxId"]` | `10239` | Root taxid of the taxonomic scope to retain (10239 = Viruses) |
| `module_2_params["nodes_path"]` | — | NCBI `nodes.dmp` for exact descendant resolution |
| `module_2_params["names_path"]` | — | NCBI `names.dmp` (taxon name lookup) |
| `module_2_params["threads"]` | falls back to `cross_module_params["threads"]` | Classification threads |

## Output Files

Under `<output_dir>/module_2/`:

| Path | Description |
|---|---|
| `<sample>/<sample>_paired.kraken`, `<sample>_paired.report` | Kraken2 per-read assignments and summary (paired) |
| `<sample>/<sample>_singleton.kraken`, `<sample>_singleton.report` | Same for singleton reads |
| `<sample>/<sample>_non_viral_{R1,R2,singleton}.fastq` | Scrubbed (non-viral) read sets |
| `<sample>/<sample>_viral_{R1,R2,singleton}.fastq` | Strict viral-only read sets |
| `<sample>/<sample>_viral_plus_unclassified_{R1,R2,singleton}.fastq` | Viral + unclassified read sets (default Module 3 input) |
| `<sample>/<sample>_dominant_virus_plus_unclassified_{R1,R2,singleton}.fastq` | Dominant-species + unclassified sets |
| `read_classification_report.tsv` | Cohort summary: `Sample`, `ReadType`, `Identity of dominant virus species`, `Dominant virus species reads`, `All viral reads`, `Human reads`, `Non-human + Non-viral reads`, `Unclassified reads` |

## Error Handling &amp; Edge Cases

- **Missing Module 1 output or mates** → `FileNotFoundError` with the expected suffixes printed.
- **No viral scope in a report** → scope falls back to report-depth inference across available reports; if the viral root taxid is absent everywhere, the scope set is empty and all reads fall into the non-viral bucket.
- **No species-level reads** → no dominant species is selected; dominant-virus outputs are skipped for that sample.
- **Kraken2 failure** → `subprocess.run(..., check=True)` raises with the exact command echoed to the log.

## Usage

```bash
python pipeline/operational/module_2_classify_scrub.py --config my_config.py
```

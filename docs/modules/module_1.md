# Module 1 — Quality Control &amp; Adapter Trimming

**Script:** [`pipeline/operational/module_1_qc_trim.py`](../../pipeline/operational/module_1_qc_trim.py)

## Purpose

Module 1 is the entry point of the operational pipeline. It performs adapter trimming and quality filtering of raw paired-end reads (and optional singleton reads) with **fastp**, and aggregates all per-sample fastp JSON reports into a single **MultiQC** report for cohort-level QC review.

## Input Requirements

- Raw FASTQ files organized as **one sub-directory per sample** under `cross_module_params["read_dir"]`:

```text
<read_dir>/
└── <sample_name>/
    ├── <sample_name>_1.fastq(.gz)      # matched by read_1_str
    ├── <sample_name>_2.fastq(.gz)      # matched by read_2_str
    └── <sample_name>.fastq(.gz)        # optional singleton, matched by singleton_str
```

- File detection is suffix-driven via `read_1_str` / `read_2_str` / `singleton_str` from `cross_module_params`. Files whose names contain `singleton` are excluded from the paired-end match; if substring matching fails, the module falls back to the explicit prefixes `<sample>_R1`, `<sample>_R2`, and `<sample>_singleton`.
- Both paired-end mates are mandatory; singletons are processed only when `singleton_str` is not `None`.

## Core Algorithm / Logic

1. **Input validation** — before any computation, every configured sample directory and read file is checked for existence (`validate_inputs`).
2. **Paired-end fastp run** — each sample is processed with:

```bash
fastp -i <R1> -I <R2> -o <sample>_R1.fastq.gz -O <sample>_R2.fastq.gz \
      --detect_adapter_for_pe --thread <threads> \
      -h <sample>_paired.html -j <sample>_paired.json
```

   `--detect_adapter_for_pe` enables automatic paired-end adapter detection.
3. **Singleton fastp run** (if present) — same tool, single-end mode.
4. **Report aggregation** — all per-sample fastp JSON files are copied into `module_1/fastp/json_files/` and MultiQC is run over that directory to produce the cohort report.
5. **Output verification** — the module raises if fastp fails to create an expected output file.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `cross_module_params["read_dir"]` | — | Directory of per-sample raw FASTQ folders |
| `cross_module_params["read_1_str"]` | `"_1.fastq"` | Suffix identifying mate 1 |
| `cross_module_params["read_2_str"]` | `"_2.fastq"` | Suffix identifying mate 2 |
| `cross_module_params["singleton_str"]` | `".fastq"` | Suffix for singleton reads; `None` disables singleton handling |
| `module_1_params["fastp_bin_path"]` | — | Path to the fastp executable |
| `module_1_params["multi_qc_bin_path"]` | `"multiqc"` | MultiQC command (pip-installed into the environment) |
| `module_1_params["threads"]` | falls back to `cross_module_params["threads"]` | fastp thread count |

## Output Files

All outputs land under `<output_dir>/module_1/`:

| Path | Description |
|---|---|
| `fastp/<sample>/<sample>_R1.fastq.gz`, `..._R2.fastq.gz` | Trimmed paired-end reads (primary input to Modules 2, 2B, 6, 6B) |
| `fastp/<sample>/<sample>_singleton.fastq.gz` | Trimmed singleton reads (if present) |
| `fastp/<sample>/<sample>_paired.{html,json}` | Per-sample fastp reports (paired) |
| `fastp/<sample>/<sample>_singleton.{html,json}` | Per-sample fastp reports (singleton) |
| `fastp/json_files/` | Aggregated JSON reports (MultiQC input) |
| `multiqc/` | Cohort-level MultiQC report |

## Error Handling &amp; Edge Cases

- **Missing sample directory or paired files** → `FileNotFoundError` before any processing begins (fail-fast validation).
- **fastp produces no output** → explicit `FileNotFoundError` after each run, preventing silent propagation of empty files downstream.
- **0-byte trimmed outputs** are not re-checked here, but Module 2B explicitly guards against them (skips alignment with a warning).

## Usage

```bash
python pipeline/operational/module_1_qc_trim.py --config my_config.py
```

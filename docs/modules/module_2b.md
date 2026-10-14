# Module 2B — Greedy Reference-Guided Read Extraction

**Script:** [`pipeline/operational/module_2b.py`](../../pipeline/operational/module_2b.py)

## Purpose

Module 2B is the entry point of the **direct-mapping pathway**, designed for low-titer or highly divergent samples in which taxonomic classification recovers too few reads for a successful *de novo* assembly. For each sample with a configured reference genome, it aligns the Module 1 trimmed reads to that reference with a deliberately permissive ("greedy") Minimap2 configuration and extracts only the mapped reads, producing a concentrated read set for Module 3 (`reads_subset = "aligned_to_ref"`) or for direct inspection.

## Input Requirements

- Trimmed paired-end FASTQs from Module 1 at `<output_dir>/module_1/fastp/<sample>/` (singletons are also extracted when present).
- A per-sample reference genome map in `module_2b_params["reference_genomes"]`:

```python
module_2b_params = {
    "minimap2_bin_path": "minimap2",
    "samtools_bin_path": "samtools",
    "threads": 8,
    "reference_genomes": {
        "SRR34884179": "/abs/path/to/module_4/ref_genomes/GCF_000855545.1.fasta",
        # Samples not listed here are SKIPPED.
    },
}
```

References are typically the genomes selected by Module 4 (e.g. `<output_dir>/module_4/ref_genomes/<accession>.fasta`), but any FASTA may be supplied.

## Core Algorithm / Logic

For each configured sample, the module runs a single piped command:

```bash
minimap2 -a -x sr -k 11 -w 3 -s 40 --secondary=no -t <threads> <ref.fa> <R1> <R2> \
  | samtools view -b -F 4 - \
  | samtools sort -n -m 1G - \
  | samtools fastq -1 <out_R1> -2 <out_R2> -s <out_singleton> -
```

The Minimap2 parameters are tuned well beyond the `sr` preset defaults to **rescue divergent and low-abundance viral reads** that a standard alignment would reject:

| Flag | Meaning | Why it matters here |
|---|---|---|
| `-x sr` | Short-read preset | Baseline for paired-end Illumina data |
| `-k 11` | K-mer size 11 (preset default: 21) | Seeds matches in highly mutated/divergent genomes |
| `-w 3` | Minimizer window 3 (preset default: 11) | Retains far more seed candidates, increasing sensitivity |
| `-s 40` | Minimum chaining score 40 (preset default: 100+) | Accepts weaker/shorter alignments typical of low-titer, partially degraded viral reads |
| `--secondary=no` | Drop secondary alignments | Only "does the read map?" matters; also faster |
| `-F 4` (samtools) | Discard unmapped | Extraction of mapped reads only |

`samtools sort -n` (name-sort) precedes `samtools fastq` so that R1/R2 pairs are correctly paired in the output files.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_2b_params["reference_genomes"]` | `{}` | **Required.** `{sample: /path/to/reference.fasta}`; unlisted samples are skipped |
| `module_2b_params["minimap2_bin_path"]` | `"minimap2"` | Minimap2 executable |
| `module_2b_params["samtools_bin_path"]` | `"samtools"` | samtools executable |
| `module_2b_params["threads"]` | falls back to `cross_module_params["threads"]` | Alignment threads |

## Output Files

Written to the Module 2 sample directory so that all read subsets stay consolidated:

| Path | Description |
|---|---|
| `<output_dir>/module_2/<sample>/<sample>_aligned_to_ref_R1.fastq` | Mapped mate-1 reads |
| `<output_dir>/module_2/<sample>/<sample>_aligned_to_ref_R2.fastq` | Mapped mate-2 reads |
| `<output_dir>/module_2/<sample>/<sample>_aligned_to_ref_singleton.fastq` | Mapped singletons |

These files are consumed by Module 3 when `module_3_params["reads_subset"] = "aligned_to_ref"`.

## Error Handling &amp; Edge Cases

- **Sample absent from `reference_genomes`** → logged `[SKIP]`; the module continues with the next sample.
- **Missing reference FASTA** → `FileNotFoundError` at pre-flight validation (`validate_inputs`) before any alignment runs.
- **Empty (0-byte) input FASTQs** — if fastp produced empty outputs (e.g. a failed lane), the module logs a warning and skips the sample instead of letting Minimap2 emit nothing downstream.
- **Piped-command failures** — `stderr` is captured and printed verbatim on non-zero exit so the exact Minimap2/samtools error is visible; a `RuntimeError` carries the return code.
- **Per-sample isolation** — an exception in one sample is logged (`[ERROR]`) without aborting the batch.

## Usage

```bash
python pipeline/operational/module_2b.py --config my_config.py
```

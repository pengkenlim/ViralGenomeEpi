# Module 3 — *De Novo* Assembly with rnaviralSPAdes

**Script:** [`pipeline/operational/module_3_assembly.py`](../../pipeline/operational/module_3_assembly.py)

## Purpose

Module 3 assembles the viral read subsets produced by Module 2 (or Module 2B) into contigs using **rnaviralSPAdes**, an RNA-virus–aware assembler from the SPAdes family. The resulting contigs are the raw material from which Module 4 selects a reference genome and Module 5 builds full-length scaffolds.

## Input Requirements

- Read subsets under `<output_dir>/module_2/<sample>/`, selected by `module_3_params["reads_subset"]`:

| `reads_subset` value | Files consumed |
|---|---|
| `viral` (default) | `<sample>_viral_plus_unclassified_{R1,R2,singleton}.fastq` |
| `dominant_virus` | `<sample>_dominant_virus_plus_unclassified_{R1,R2,singleton}.fastq` |
| `aligned_to_ref` | `<sample>_aligned_to_ref_{R1,R2,singleton}.fastq` (Module 2B output) |

  Aliases such as `viral_plus_unclassified`, `dominant virus only`, and `aligned to ref` are normalized to the three canonical modes.
- Paired-end reads (R1 + R2) are mandatory; singletons are passed to SPAdes when present.

## Core Algorithm / Logic

1. **Subset resolution and validation** — the configured `reads_subset` mode is normalized and every sample's input files are verified to exist before any assembly starts.
2. **Assembly** — per sample:

```bash
rnaviralspades.py -o <output_dir> -t <threads> -m <memory> \
    [--rna] [-k 21,33,55,77,99,127] \
    -1 <R1> -2 <R2> [-s <singleton>]
```

   - `--rna` (enabled by default via `use_rna_mode`) invokes RNA-seq assembly behavior appropriate for RNA viral genomes.
   - `over_assembly = True` widens the k-mer sweep (`-k 21,33,55,77,99,127`) to improve recovery of low-coverage or heterogeneous genomes at the cost of runtime; the output directory is suffixed `_over_assembly` so runs never clobber each other.
3. **Output verification** — the module logs success if `scaffolds.fasta` or `contigs.fasta` appears, and warns otherwise (assembly can legitimately fail on very sparse read sets).

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_3_params["rnaviralspades_bin_path"]` | — | Path to `rnaviralspades.py` (SPAdes ≥ 4.x) |
| `module_3_params["reads_subset"]` | `"viral"` | Which Module 2/2B read set to assemble |
| `module_3_params["use_rna_mode"]` | `True` | Pass `--rna` to the assembler |
| `module_3_params["over_assembly"]` | `False` | Enable extended k-mer sweep (`-k 21,33,55,77,99,127`) |
| `module_3_params["memory"]` | `500` | Memory limit passed to SPAdes `-m` (GB in the template; note SPAdes interprets `-m` as GB) |
| `module_3_params["threads"]` | falls back to `cross_module_params["threads"]` | Assembly threads |

## Output Files

Under `<output_dir>/module_3/`:

| Path | Description |
|---|---|
| `<mode>/<sample>/` (mode ∈ `viral`, `dominant_virus`, `aligned_to_ref`, plus `_over_assembly` variants) | Full SPAdes output: `transcripts.fasta`, `contigs.fasta`, `scaffolds.fasta`, graphs, logs |
| `<sample>/combined_transcripts.fasta` | Combined assembly FASTA written by Module 4 (deduplicated headers) — the actual input to reference selection and Module 5 scaffolding |

## Error Handling &amp; Edge Cases

- **Missing Module 2 directory or subset reads** → `FileNotFoundError` at validation, with a hint to run Module 2 (and optionally Module 2B) first.
- **Missing assembler binary** → `FileNotFoundError` before launching.
- **Empty/invalid `reads_subset` string** → `ValueError` listing the accepted aliases.
- **Assembly failure on one sample** (e.g. insufficient viral reads) → logged as `[ERROR]`; the batch continues. Downstream modules (4, 5) will skip samples lacking `combined_transcripts.fasta`.
- **Low-titer samples that fail here** are candidates for the direct-mapping pathway (Modules 2B → 6B), which does not require assembly.

## Usage

```bash
python pipeline/operational/module_3_assembly.py --config my_config.py
```

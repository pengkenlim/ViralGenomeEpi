# Module 7 — Multiple Sequence Alignment &amp; Phylogeny

**Script:** [`pipeline/operational/module_7_tree_build.py`](../../pipeline/operational/module_7_tree_build.py)

## Purpose

Module 7 collects the per-sample consensus genomes, builds a multiple sequence alignment (MSA) with **MAFFT**, trims the ragged alignment termini, and infers an approximate maximum-likelihood phylogeny with **FastTree**. It can run on any subset of samples under a named analysis, and can optionally include a reference genome as an outgroup — which also unlocks the gene-aware SNP analysis in Module 8.

## Input Requirements

- Consensus FASTAs, resolved per sample according to the **`consensus_source` toggle** in `module_7_params`:

| `consensus_source` | Directory searched | Preferred filenames |
|---|---|---|
| `"module_6"` (default) | `<output_dir>/module_6/<sample>/` | `<sample>_consensus_combined.fasta`, then `.../final_polished_consensus.fa(fasta)`, `.../polished_consensus.fa(fasta)` |
| `"module_6b"` | `<output_dir>/module_6b/<sample>/` | `<sample>_direct_consensus.fasta` (`.fa`) |

  If none of the preferred names exists, the module falls back to the first FASTA found recursively under the sample's consensus directory.
- Optional `msa_reference`: either a filesystem path to a reference FASTA or an accession resolvable against `<output_dir>/module_4/ref_genomes/` (`.fasta` / `.fa` / `.fna`).

## Core Algorithm / Logic

1. **Sample collection** — each input sample contributes exactly one record to the MSA, re-headed to the sample accession. Samples without a locatable consensus are logged and skipped; the run aborts only if *no* valid records remain.
2. **Workflow 1 (always): samples-only phylogeny** —
   - MAFFT `--auto --thread N` on the combined FASTA,
   - **terminal end-trimming**: 5′/3′ alignment columns are removed inward until the first column containing no gap (`-`) or `N`/`n` in *any* sequence, eliminating the ragged low-confidence ends characteristic of viral consensus genomes,
   - FastTree `-nt -gtr` (GTR nucleotide model) on the trimmed MSA.
3. **Workflow 2 (optional): reference-included analysis** — when `msa_reference` is set, the reference record (ID normalized with spaces → underscores) is appended and the full combine → align → trim → tree pipeline is repeated with `ref_`-prefixed filenames. The reference acts as an **outgroup for rooting** the tree, and the *untrimmed* reference-included alignment (`polished_genomes_ref_aligned.fasta`) is what Module 8 slices per gene for its gene-aware SNP counting.
4. **Metadata** — `analysis_metadata.json` records `analysis_name`, `input_samples`, `consensus_source`, `msa_reference`, and whether the reference-included tree was built.

The `consensus_source` toggle is what allows the **direct-mapping pathway** (Module 6B) to feed the same downstream phylogenetics: set `"module_6b"` and Module 7 transparently substitutes `module_6b/<sample>/<sample>_direct_consensus.fasta` for the Module 6 outputs.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_7_params["mafft_bin_path"]` | — | MAFFT executable |
| `module_7_params["fasttree_bin_path"]` | — | FastTree executable |
| `module_7_params["consensus_source"]` | `"module_6"` | `"module_6"` (de novo polished) or `"module_6b"` (direct mapping) |
| `module_7_params["msa_reference"]` | `None` | Reference accession or FASTA path used as outgroup; also the annotation source for Module 8 gene-aware analysis |
| `module_7_params["input_samples"]` | all `cross_module_params["samples"]` | Sample subset for this analysis |
| `module_7_params["analysis_name"]` | `"default"` | Names the output subdirectory, enabling multiple independent analyses |
| `module_7_params["threads"]` | falls back to `cross_module_params["threads"]` | MAFFT threads |

## Output Files

Under `<output_dir>/module_7/<analysis_name>/`:

| Path | Description |
|---|---|
| `polished_genomes_combined.fasta` | Concatenated input consensus genomes |
| `polished_genomes_aligned.fasta` | MAFFT MSA (samples only) |
| `polished_genomes_aligned_trimmed.fasta` | End-trimmed MSA — **the input to Module 8's SNP pipeline** |
| `polished_genomes_tree.nwk` | FastTree phylogeny (samples only) |
| `polished_genomes_ref_combined.fasta` | Input genomes + reference (workflow 2) |
| `polished_genomes_ref_aligned.fasta` | Reference-included MSA, **untrimmed** — basis of Module 8 per-gene slicing |
| `polished_genomes_aligned_ref_trimmed.fasta` | Trimmed reference-included MSA |
| `polished_genomes_ref_tree.nwk` | Outgroup-included phylogeny |
| `analysis_metadata.json` | Provenance record for the analysis |

## Error Handling &amp; Edge Cases

- **Missing consensus directory / FASTA for a sample** → `[SKIP]` with the reason; other samples proceed.
- **No valid consensus records at all** → `RuntimeError` (nothing to align).
- **Unresolved `msa_reference`** → `FileNotFoundError` (checked as file path, then in `module_4/ref_genomes/`).
- **Degenerate trimming** — if trimming would remove the entire alignment (no fully gap/N-free column exists), the untrimmed alignment is kept and a warning is printed.
- **Unconfigured executables** → `FileNotFoundError` naming the missing config key; `PATH`-resolvable names are accepted.

## Usage

```bash
python pipeline/operational/module_7_tree_build.py --config my_config.py
```

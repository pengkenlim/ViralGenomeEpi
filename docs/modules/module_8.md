# Module 8 — Gene-Aware SNP Distances &amp; Transmission Clustering

**Script:** [`pipeline/operational/module_8_snpsites_hcluster.py`](../../pipeline/operational/module_8_snpsites_hcluster.py)

## Purpose

Module 8 converts the Module 7 alignment into epidemiological evidence: it extracts variant sites with **snp-sites**, computes **pairwise SNP distances**, performs **hierarchical clustering** at a configurable SNP cutoff, and — when a reference was included in Module 7 — decomposes distances **per gene** against the reference annotation. The result is a set of transmission-cluster assignments with full pairwise distance provenance.

## Input Requirements

- Module 7 outputs under `<output_dir>/module_7/<analysis_name>/`:
  - `polished_genomes_aligned_trimmed.fasta` (trimmed, samples-only MSA) — the basis for the global SNP matrix and clustering. Fallback names `genomes_aligned_trimmed.fasta` / `aligned_trimmed.fasta` are also recognized.
  - `polished_genomes_ref_aligned.fasta` (untrimmed, reference-included MSA) — required for the optional gene-aware analysis.
- For gene-aware mode, `module_7_params["msa_reference"]` must be a **GCF_/GCA_ accession** so that the matching RefSeq GFF can be fetched.
- `snp-sites` executable; `datasets` CLI or network access for GFF download (three fallbacks: NCBI Datasets CLI → Datasets REST API → NCBI FTP), cached under `module_8/<analysis_name>/reference_annotation/`.

## Core Algorithm / Logic

### Gap-stripped pairwise SNP counting

For every sample pair, the distance is computed over the trimmed MSA using **only columns where both sequences carry canonical bases (A/C/G/T)**:

```text
SNP_Distance = Σ columns where base_i ≠ base_j AND base_i, base_j ∈ {A, C, G, T}
```

Columns containing a gap (`-`) or ambiguous/`N` base in *either* sequence are excluded from both the numerator and the comparable denominator. This **gap-stripping** is essential: a single 1-bp assembly indel shifts every downstream column in a gapped alignment, and naive character-by-character comparison would count each shifted position as a difference — turning one artifactual indel into hundreds of false SNPs (a failure mode routinely seen in genomes flagged "bad" by Nextclade QC). Gap-stripped counting confines each true difference to exactly one counted SNP, preventing assembly artifacts from splitting genuine transmission clusters.

The same canonical-base rule is applied by `compute_pairwise_snp_matrix` to build the square `pairwise_snp_matrix.csv` and the long-format `pairwise_snp_pairs.csv` (`Global SNP Difference` column).

### Gene-aware SNP decomposition (optional)

When `msa_reference` is configured and `polished_genomes_ref_aligned.fasta` exists:

1. The reference **GFF** is downloaded and its `gene` features parsed (name from `Name=`, then `locus_tag`, then `ID`).
2. Alignment columns are lifted to reference coordinates via the reference row's non-gap positions; each gene's interval defines a column mask.
3. For every pair (sample–sample **and** sample–reference), per-gene SNP counts and compared-base counts are computed with the same gap-stripped rule, restricted to the gene's columns.
4. **Low-coverage flagging** — a pair/gene row is marked `Low_Coverage_Flag = True` when `Compared_Bases < 90%` of the gene's aligned columns (`LOW_COV_FRAC = 0.90`), warning that the per-gene distance is based on limited overlap and should not be over-interpreted.
5. **Consistency guard** — if a sample–reference pair has >10 global SNPs but zero summed gene SNPs, a warning is emitted (indicates a coordinate-mapping problem).
6. Per-gene, **untrimmed MSA slices** (`genome_ref_<gene>_<start>_<end>_<strand>.fasta`) are written for visual inspection.
7. A built-in consistency check plus integer (`Int64`, `NA` for uncomputable pairs) output formatting keeps the wide table analysis-ready.

Gene-aware mode degrades gracefully: without `msa_reference`, without the ref-included alignment, or if the GFF download fails, Module 8 runs in legacy samples-only mode and records the reason in `analysis_metadata.json` (`gene_aware_analysis: false`).

### Clustering

- SciPy hierarchical clustering (`scipy.cluster.hierarchy.linkage`) on the square distance matrix with the configured `linkage_method` (template default `single`).
- Flat clusters are cut at `snp_cutoff` (default **3** SNPs) via `fcluster(..., criterion="distance")` and written as `clusters_<N>snp_cutoff.csv`.
- A SNP-distance histogram and a dendrogram (both with the cutoff line drawn when `draw_cutoff_line = True`) are rendered at 300 dpi.

### QC gating recommendation

**Exclude samples with Nextclade "bad" QC (frameshifts, pervasive ambiguous calls) from `input_samples` before clustering.** Gap-stripping protects *retained* pairs from artifact inflation, but frameshifted genomes can still distort alignment structure globally; the cleanest practice is to gate on Nextclade QC upstream and let Module 8 cluster only high-confidence genomes.

## Key Configuration Parameters

| Parameter (dict) | Default | Description |
|---|---|---|
| `module_8_params["snpsites_bin_path"]` | — | snp-sites executable |
| `module_8_params["snp_cutoff"]` | `3` | Maximum pairwise SNP distance defining a transmission cluster |
| `module_8_params["linkage_method"]` | `"single"` (template) | SciPy linkage: `single`, `average`, `complete`, … |
| `module_8_params["draw_cutoff_line"]` | `True` | Draw the cutoff marker on histogram and dendrogram |
| `module_8_params["input_samples"]` | module 7's, else all samples | Samples included in the clustering run |
| `module_8_params["analysis_name"]` | module 7's, else `"default"` | Output subdirectory under `module_8/` |
| `module_8_params["datasets_bin_path"]` | `PATH` lookup | Optional explicit NCBI `datasets` CLI for GFF download |
| `module_7_params["msa_reference"]` | `None` | GCF/GCA accession enabling the gene-aware analysis |

## Output Files

Under `<output_dir>/module_8/<analysis_name>/`:

| Path | Description |
|---|---|
| `subset_accessions_alignment.fasta` | Trimmed MSA subset to `input_samples` |
| `variants.vcf` | snp-sites variant calls |
| `snp_alignment.fasta` | SNP-only alignment |
| `pairwise_snp_matrix.csv` | Square SNP distance matrix (gap-stripped) |
| `pairwise_snp_pairs.csv` | Long-format pairs: `Sample_1`, `Sample_2`, `Global SNP Difference`, plus per-gene `<gene> SNP Difference` columns and sample–reference rows when gene-aware mode is active |
| `pairwise_snp_pairs_per_gene_long.csv` | Tidy per-gene table: `Sample_1`, `Sample_2`, `Gene`, `SNP_Distance`, `Compared_Bases`, `Gene_Aligned_Columns`, `Low_Coverage_Flag`, `Pair_Source` (`sample-sample` / `sample-reference`) |
| `pairwise_snp_histogram.png` | Distribution of pairwise SNP distances with cutoff line |
| `hierarchical_clustering_dendrogram.png` | Clustering dendrogram colored/cut at `snp_cutoff` |
| `clusters_<N>snp_cutoff.csv` | `Sample_ID` → `Cluster_ID` assignments |
| `genome_ref_<gene>_*.fasta` | Untrimmed per-gene MSA slices (gene-aware mode) |
| `reference_annotation/<acc>.gff` | Cached reference GFF |
| `analysis_metadata.json` | Provenance: samples, gene-aware status, GFF path, gene intervals, fallback reason |

## Error Handling &amp; Edge Cases

- **Missing Module 7 directory or trimmed alignment** → `FileNotFoundError` (with candidate-name search).
- **Requested samples missing from the alignment** → warned and dropped; a `ValueError` is raised only if *none* are found.
- **GFF unavailable / unparseable / no `gene` features** → gene-aware mode skipped with a logged warning; legacy analysis proceeds.
- **Reference row absent from the ref-included alignment** → gene-aware mode skipped.
- **All three GFF download routes fail** → `RuntimeError` listing each route's error (datasets CLI, REST API, FTP).
- **Duplicate gene names in the GFF** → first occurrence kept, warning logged.

## Usage

```bash
python pipeline/operational/module_8_snpsites_hcluster.py --config my_config.py
```

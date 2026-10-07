# Operational pipeline

Config-driven, module-based workflow for ViralGenomeEpi (modules 1–9). It supersedes the legacy one-off scripts in [../basic](../basic). Every module takes the same single argument, `--config <path>`, and writes into `<output_dir>/module_N/`.

## Configuration

All run settings live in a Python config file; [sample_config.py](sample_config.py) is the template (the repo root also has per-dataset examples such as `config_RSV_A.py`). It is loaded as a Python module and defines:

- `cross_module_params`: `read_dir`, `samples`, `output_dir`, `read_1_str` / `read_2_str` / `singleton_str` (FASTQ suffixes; `singleton_str=None` ignores singletons), default `threads`.
- `module_1_params` … `module_9_params`: per-module tool paths, thread counts and thresholds.
- `module_7_params`, `module_8_params`, `module_9_params` accept `input_samples` to run on a subset; modules 7 and 8 also use `analysis_name` to create `module_7/<analysis_name>` and `module_8/<analysis_name>`.

Placeholders such as `/PATH/TO/...` and the empty `output_dir` must be filled in. `read_dir` is expected to contain one sub-directory per sample holding its FASTQ files.

### Parameters per module

| Module | Key parameters |
|---|---|
| 1 | `fastp_bin_path`, `multi_qc_bin_path`, `threads` |
| 2 | `Kraken2_k2_bin_path`, `kraken2_db_path`, `scope_to_keep_taxId` (10239 = Viruses), `threads` |
| 3 | `rnaviralspades_bin_path`, `reads_subset` (`viral` or `dominant_virus`), `use_rna_mode`, `over_assembly`, `memory` (passed to SPAdes `-m`, GB), `threads` |
| 4 | `blastn_bin_path`, `blastdbcmd_bin_path`, `cd_hit_est_bin_path`, `blast_db_path`, `seq2genome_map_path`, `assembly_data_report_path`, `nodes_path`, `names_path`, `threads` |
| 5 | `checkv_bin_path`, `checkv_db_dir`, `minimap2_bin_path`, optional `fallback_strategy` (`technique_first` default, or `stringency_first`), `threads` |
| 6 | `minimap2_bin_path`, `samtools_bin_path`, `ivar_bin_path`, `checkv_bin_path`, `checkv_db_dir`, `threads` |
| 7 | `mafft_bin_path`, `fasttree_bin_path`, `input_samples`, `analysis_name`, `threads` |
| 8 | `snpsites_bin_path`, `snp_cutoff` (default 3), `linkage_method` (SciPy; `single` in the template), `input_samples`, `analysis_name`, `threads` |
| 9 | `blastn_bin_path`, `blast_db_path` (user-defined DB prefix), `evalue`, `max_target_seqs`, `outfmt`, `input_samples`, `threads` |

## Modules

| # | Script | What it does | Main outputs (under `<output_dir>/`) |
|---|---|---|---|
| 1 | [module_1_qc_trim.py](module_1_qc_trim.py) | Adapter trimming and quality filtering with fastp; aggregates fastp JSON reports with MultiQC | `module_1/fastp/<sample>/`, `module_1/fastp/json_files/`, `module_1/multiqc/` |
| 2 | [module_2_classify_scrub.py](module_2_classify_scrub.py) | Kraken2 classification (`k2 classify`); splits reads into non-viral, viral, viral + unclassified, and dominant-virus + unclassified sets | `module_2/<sample>/<sample>_{non_viral,viral,viral_plus_unclassified,dominant_virus_plus_unclassified}_{R1,R2,singleton}.fastq`, `module_2/read_classification_report.tsv` |
| 3 | [module_3_assembly.py](module_3_assembly.py) | De novo assembly with rnaviralSPAdes on the chosen read subset | `module_3/viral/<sample>/` or `module_3/dominant_virus/<sample>/` (`_over_assembly` suffix when `over_assembly` is on) |
| 4 | [module_4_download_ref_genome.py](module_4_download_ref_genome.py) | BLASTs transcripts against the viral RefSeq DB, clusters hits (cd-hit-est, 80% identity), picks the best reference by aggregated bitscore (hits ≥500 bp and ≥70% identity), resolves taxonomy, detects segmented viruses, and downloads the genome from the NCBI Datasets API | `module_4/<sample>/<sample>_reference_summary.json`, `selected_reference_<acc>.fasta`, `module_4/ref_genomes/`, `module_4/Identified_species.tsv` |
| 5 | [module_5_scaffold_qc.py](module_5_scaffold_qc.py) | Reference-guided scaffolding with minimap2 using a 6-tier fallback (see below), then CheckV QC; handles segmented viruses per segment | `module_5/<sample>/<sample>[_seg_<id>]_scaffold.fasta`, `<sample>_scaffold_checkv_report.txt`, `module_5/Identified_scaffolds.tsv` |
| 6 | [module_6_consensus_polish.py](module_6_consensus_polish.py) | Maps reads back to the scaffold (minimap2 + samtools) and calls consensus with iVar in a strict then relaxed pass, then CheckV QC | `module_6/<sample>/<sample>_consensus_combined.fasta`, `*_aligned_sorted.bam`, `*_ivar_polished_report.txt`, `module_6/Identified_consensus.tsv` |
| 7 | [module_7_tree_build.py](module_7_tree_build.py) | Combines polished genomes, aligns with MAFFT, trims terminal gap/N columns, builds a FastTree tree | `module_7/<analysis_name>/polished_genomes_{combined,aligned,aligned_trimmed}.fasta`, `polished_genomes_tree.nwk`, `analysis_metadata.json` |
| 8 | [module_8_snpsites_hcluster.py](module_8_snpsites_hcluster.py) | Subsets the trimmed alignment, extracts SNPs with snp-sites, computes pairwise SNP distances, hierarchical clustering at `snp_cutoff` | `module_8/<analysis_name>/variants.vcf`, `snp_alignment.fasta`, `pairwise_snp_matrix.csv`, `pairwise_snp_pairs.csv`, `pairwise_snp_histogram.png`, `hierarchical_clustering_dendrogram.png`, `analysis_metadata.json` |
| 9 | [module_9_blast_polished_genomes.py](module_9_blast_polished_genomes.py) | Optional: BLASTs polished consensus genomes (module 6) against a user-defined DB (e.g. NCBI core nt) to confirm strain/taxon; expects `staxids` and `sscinames` in `outfmt` | `module_9/<sample>/<sample>_module_9_blastn.tsv`, `module_9/Identified_blast_hits.tsv` |

Modules 7–9 are run on selected samples, so they can be repeated with different `input_samples` / `analysis_name`.

### Module 5 scaffolding tiers

Three minimap2 modes (`conserved`: `-c -x asm5`, ≥70% identity; `divergent_70`: `-c -k 14 -w 5`, ≥70%; `divergent_60`: same, ≥60%) are combined with two techniques:

- **slice**: one contig block covering ≥98% of the reference, scaffold ≤5% Ns
- **stitch**: merged blocks covering ≥98% of the reference, scaffold ≤2% Ns

`technique_first` (default) tries slice at all three modes, then stitch at all three. `stringency_first` tries slice then stitch per mode. The first accepted tier wins and is recorded as `alignment_mode` in `Identified_scaffolds.tsv`.

### Module 6 iVar passes

Strict: min quality 20, min frequency 0.60, min depth 10. Relaxed (rescues low-coverage termini): min quality 15, min frequency 0.51, min depth 3.

## Running

```bash
CFG=pipeline/operational/sample_config.py   # or your own config
for m in 1_qc_trim 2_classify_scrub 3_assembly 4_download_ref_genome \
         5_scaffold_qc 6_consensus_polish 7_tree_build 8_snpsites_hcluster \
         9_blast_polished_genomes; do
  python pipeline/operational/module_$m.py --config $CFG || break
done
```

Module 9 is optional.

## Requirements

- Python 3 with Biopython, SciPy, NumPy, pandas, matplotlib, seaborn (as imported by modules 7–8), plus `multiqc` and `checkv` (pip-installable)
- External tools: fastp, Kraken2 (`k2`) and DB, rnaviralSPAdes (SPAdes ≥ 4), BLAST+ (`blastn`, `blastdbcmd`), cd-hit-est, minimap2, samtools, iVar, MAFFT, FastTree, snp-sites, CheckV DB
- Module 4 reference data: viral RefSeq BLAST DB, `seq2genome_map.json` (built by [../../setup/prepare_database.py](../../setup/prepare_database.py)), `assembly_data_report.jsonl`, NCBI taxonomy `nodes.dmp` / `names.dmp`; network access to the NCBI Datasets API
- Module 4 falls back to hardcoded default paths under `/home/mngs/...` if the matching config keys are missing, so set them explicitly.

## Tests

```bash
pytest tests/
```

Currently covers module 9 ([../../tests](../../tests)).

## See also

- [../../README.md](../../README.md): repository overview
- [../basic](../basic): legacy pipeline scripts

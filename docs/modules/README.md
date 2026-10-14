# ViralGenomeEpi — Module Documentation

Detailed, per-module documentation for the operational pipeline ([`pipeline/operational/`](../../pipeline/operational)). Each page covers purpose, exact input requirements, core algorithm, key configuration parameters, output files, and error handling. For a high-level overview, quick-start instructions, and the pipeline architecture diagram, see the [root README](../../README.md).

## Main pathway (Modules 1 → 8)

| # | Module | Script | Document |
|---|---|---|---|
| 1 | QC &amp; adapter trimming (fastp / MultiQC) | `module_1_qc_trim.py` | [module_1.md](module_1.md) |
| 2 | Kraken2 classification &amp; read scrubbing | `module_2_classify_scrub.py` | [module_2.md](module_2.md) |
| 3 | *De novo* assembly (rnaviralSPAdes) | `module_3_assembly.py` | [module_3.md](module_3.md) |
| 4 | Reference selection &amp; download (BLAST + cd-hit-est + NCBI) | `module_4_download_ref_genome.py` | [module_4.md](module_4.md) |
| 5 | 6-tier slice/stitch scaffolding + CheckV | `module_5_scaffold_qc.py` | [module_5.md](module_5.md) |
| 6 | Two-tier iVar consensus polishing | `module_6_consensus_polish.py` | [module_6.md](module_6.md) |
| 7 | MAFFT alignment &amp; FastTree phylogeny | `module_7_tree_build.py` | [module_7.md](module_7.md) |
| 8 | Gap-stripped, gene-aware SNP clustering | `module_8_snpsites_hcluster.py` | [module_8.md](module_8.md) |

## Alternative direct-mapping pathway

| # | Module | Script | Document |
|---|---|---|---|
| 2B | Greedy Minimap2 read extraction (`-x sr -k 11 -w 3 -s 40`) | `module_2b.py` | [module_2b.md](module_2b.md) |
| 6B | Direct reference-mapping consensus (`samtools mpileup \| ivar consensus`) | `module_6b_direct_consesnsus.py` | [module_6b.md](module_6b.md) |

The direct-mapping pathway bypasses assembly for low-titer/divergent samples: **Module 1 → Module 2B → (optional Module 3 with `reads_subset = "aligned_to_ref"`) → Module 6B**, then rejoins the main flow at Module 7 by setting `module_7_params["consensus_source"] = "module_6b"`.

## Optional module

Module 9 (`module_9_blast_polished_genomes.py`) BLASTs the polished consensus genomes against a user-defined database for strain/taxon confirmation. It is optional and configured via `module_9_params` (`blastn_bin_path`, `blast_db_path`, `evalue`, `max_target_seqs`, `outfmt`, `input_samples`); outputs land under `<output_dir>/module_9/` (`<sample>_module_9_blastn.tsv`, `Identified_blast_hits.tsv`).

## Conventions shared by all modules

- Every module takes exactly one argument: `--config <path-to-config.py>`.
- Every module writes under `<output_dir>/module_<N>/` (Module 2B writes into `module_2/`; Module 6B into `module_6b/`).
- Per-module settings live in `module_<N>_params` dictionaries; shared settings in `cross_module_params` (see the [annotated template](../../pipeline/operational/sample_config.py)).
- Modules 7 and 8 support `input_samples` + `analysis_name` for reproducible, named re-analysis of sample subsets.

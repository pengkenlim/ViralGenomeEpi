# ViralGenomeEpi

Workflow for viral genome recovery and epidemiological analysis from metagenomic sequencing data. Adapted from https://pubmed.ncbi.nlm.nih.gov/40313286/; the pipeline has been validated on that paper's sequencing data by reproducing the reported viral linkages.

## Status

Research pipeline, not a turnkey package. The **operational pipeline** in [pipeline/operational](pipeline/operational) is config-driven (`--config`) and is the version to use. The original hardcoded scripts are kept in [pipeline/basic](pipeline/basic) for reference only.

## Operational workflow

| # | Module | Purpose |
|---|---|---|
| 1 | `module_1_qc_trim.py` | fastp trimming/QC, MultiQC report |
| 2 | `module_2_classify_scrub.py` | Kraken2 classification; split into non-viral, viral, viral + unclassified and dominant-virus read sets |
| 3 | `module_3_assembly.py` | rnaviralSPAdes assembly of the selected read subset |
| 4 | `module_4_download_ref_genome.py` | BLAST + cd-hit-est reference selection, taxonomy resolution, segmented-virus detection, NCBI genome download |
| 5 | `module_5_scaffold_qc.py` | minimap2 reference-guided scaffolding (6-tier slice/stitch fallback), optional contig exclusion before scaffolding, and CheckV QC |
| 6 | `module_6_consensus_polish.py` | read re-mapping and iVar consensus (strict then relaxed pass), CheckV QC |
| 7 | `module_7_tree_build.py` | MAFFT alignment, end trimming, FastTree phylogeny |
| 8 | `module_8_snpsites_hcluster.py` | snp-sites, pairwise SNP distances, hierarchical clustering at a SNP cutoff |
| 9 | `module_9_blast_polished_genomes.py` | optional BLAST of polished genomes against a user-defined DB for strain/taxon confirmation |

Modules 7–9 can be re-run on sample subsets, with an `analysis_name` separating outputs. Full details (parameters, outputs, tiers and thresholds) are in [pipeline/operational/README.md](pipeline/operational/README.md).

## Repository layout

- [pipeline/operational/](pipeline/operational): modules 1–9 and [sample_config.py](pipeline/operational/sample_config.py), the config template
- [pipeline/basic/](pipeline/basic): legacy scripts (Kraken2 scrub, SPAdes, Mash/BLAST reference selection, nucmer + CheckV, iVar, MAFFT/FastTree, snp-sites clustering)
- [setup/prepare_database.py](setup/prepare_database.py): builds `seq2genome_map.json` and a combined FASTA from a dehydrated NCBI viral genome download (paths inside are hardcoded; edit before use)
- [tests/](tests): pytest tests (currently module 9)
- `config_*.py` (RSV A/B, HMPV A2/B1-2, RhV A/B): example per-dataset configs from past runs
- `dependencies/` and `test_data/output/` are git-ignored local directories for tools, databases and test outputs

## Quick start

1. Copy [sample_config.py](pipeline/operational/sample_config.py) and set `read_dir`, `samples`, `output_dir`, tool and database paths.
2. Run the modules in order:

```bash
CFG=my_config.py
for m in 1_qc_trim 2_classify_scrub 3_assembly 4_download_ref_genome \
         5_scaffold_qc 6_consensus_polish 7_tree_build 8_snpsites_hcluster; do
  python pipeline/operational/module_$m.py --config $CFG || break
done
python pipeline/operational/module_9_blast_polished_genomes.py --config $CFG   # optional
```

Test one sample first before running a full batch.

## Inputs and dependencies

- Raw FASTQ reads: one directory per sample under `read_dir`, paired-end with optional singletons
- Databases: Kraken2 DB, viral RefSeq BLAST DB with `seq2genome_map.json`, `assembly_data_report.jsonl`, NCBI taxonomy dumps (`nodes.dmp`, `names.dmp`), CheckV DB; optionally a BLAST DB for module 9
- Tools: fastp, MultiQC, Kraken2, rnaviralSPAdes, BLAST+, cd-hit-est, minimap2, samtools, iVar, CheckV, MAFFT, FastTree, snp-sites
- Python: Biopython, SciPy, NumPy, pandas, matplotlib, seaborn
- Network access to the NCBI Datasets API for reference download (module 4)

## Caveats

- Tool paths are set in the config and some defaults in module 4 point to local paths; set them explicitly.
- Limited error handling; failed samples or ambiguous reference choices may need manual inspection.
- Not yet packaged (conda/Docker/Snakemake) and few automated tests.

## License

See [LICENSE](LICENSE).

# Operational pipeline

This directory contains the current module-based workflow for ViralGenomeEpi. It replaces the legacy one-off scripts with a clearer run order, a shared configuration file, and explicit module outputs.

## Purpose

The operational pipeline is designed to make the workflow more reproducible and easier to adapt across datasets. Each module handles one stage of the analysis and writes outputs into a predictable location under the project output directory.

## Main configuration

The shared run settings live in [sample_config.py](sample_config.py). That file defines:

- the input samples
- the raw read directory and file suffixes
- the shared output root
- the tool paths for each module
- module-specific thresholds and analysis names

The config is intended to be edited before running the workflow on a new system or sequencing project.

## Module order

The operational workflow is organized in this order:

1. Module 1: QC trimming and read filtering
2. Module 2: taxonomic classification and viral read selection
3. Module 3: assembly of viral reads
4. Module 4: reference genome selection and download
5. Module 5: scaffold QC and coverage validation
6. Module 6: read mapping and consensus polishing
7. Module 7: alignment and phylogenetic tree construction
8. Module 8: SNP extraction and clustering

## Operational scripts

- [module_1_qc_trim.py](module_1_qc_trim.py): trims adapters and quality-filters reads
- [module_2_classify_scrub.py](module_2_classify_scrub.py): classifies reads with Kraken2 and keeps viral reads
- [module_3_assembly.py](module_3_assembly.py): assembles viral contigs/scaffolds
- [module_4_download_ref_genome.py](module_4_download_ref_genome.py): chooses the best reference genome for each sample and downloads it
- [module_5_scafold_qc.py](module_5_scafold_qc.py): maps scaffolds to the reference and runs QC checks
- [module_6_consensus_polish.py](module_6_consensus_polish.py): maps reads back to the scaffold and calls the polished consensus
- [module_7_tree_build.py](module_7_tree_build.py): aligns consensus sequences and builds a phylogeny
- [module_8_snpsites_hcluster.py](module_8_snpsites_hcluster.py): extracts SNPs and clusters sequences by genetic distance

## Expected output layout

Outputs are organized under the root output directory, usually something like:

- output/module_1
- output/module_2
- output/module_3
- output/module_4
- output/module_5
- output/module_6
- output/module_7
- output/module_8

For each analysis, module 7 and module 8 use an `analysis_name` value to create a subdirectory such as `module_7/RSV_group1` or `module_8/RSV_group1`.

## How to run the workflow

Activate the project environment and then execute the modules in order with the config file:

```bash
source pipeline/operational/.venv/bin/activate
python pipeline/operational/module_1_qc_trim.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_2_classify_scrub.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_3_assembly.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_4_download_ref_genome.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_5_scafold_qc.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_6_consensus_polish.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_7_tree_build.py --config pipeline/operational/sample_config.py
python pipeline/operational/module_8_snpsites_hcluster.py --config pipeline/operational/sample_config.py
```

## Notes on configuration

Before running the workflow:

- review all binary paths in [sample_config.py](sample_config.py)
- confirm the output directory is writable
- ensure the sample list matches your dataset
- verify that the reference database and CheckV DB paths exist
- update the `analysis_name` when running a distinct analysis group

## Relationship to the legacy pipeline

The legacy scripts under [../basic](../basic) are still useful for historical context and debugging. The operational pipeline is the preferred implementation for running the project in a structured way.

## See also

- [../../README.md](../../README.md) for the repository overview
- [../basic/README.md](../basic/README.md) for the legacy pipeline guide

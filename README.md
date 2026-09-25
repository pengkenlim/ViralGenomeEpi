# ViralGenomeEpi

Prototype workflow for viral genome recovery and epidemiological analysis from metagenomic sequencing data. Pipeline adapted from https://pubmed.ncbi.nlm.nih.gov/40313286/. The pipeline has been validated with sequencing data from said paper by reproducing viral linkages. 

## Status

This repository is a research/development pipeline, not a production-ready turnkey package. The scripts are written as workflow steps for a specific computational environment and are designed to be adapted and debugged as needed. Several values are hardcoded to local paths and project-specific accession lists. More work needs to be done.

## Overview

In essence, the pipeline combines:

- read-level taxonomic filtering with Kraken2
- viral read enrichment and cleaning
- de novo assembly with SPAdes/rSPAdes
- reference selection via Mash and/or BLAST-taxonomy logic
- reference-guided scaffolding with MUMmer
- quality check with CheckV
- consensus polishing with iVar
- multiple sequence alignment and phylogenetic tree building with MAFFT/FastTree
- SNP-based clustering using snp-sites and hierarchical clustering

This is intended to support a batch-style metagenomic viral genomics project, especially for epidemiological comparison across multiple accessions.

## Repository layout

- `1_run_kraken2_scrub_human.py` — viral read enrichment, Kraken2 classification, read filtering, and output FASTQ generation
- `2_default_spades_assembly.py` — de novo viral assembly from enriched reads
- `3_mash_refgenome_dowload.py` — reference selection using Mash distance and NCBI genome retrieval
- `3_blast_taxid_refgenome_download.py` — alternative reference selection using BLAST taxonomy-weighting and NCBI representative genome retrieval
- `4_nucmer_scaffold_checkv.py` — reference-guided contig scaffolding, overlap resolution, and CheckV evaluation
- `5_ivar_consensus_polishing.py` — read mapping and polished consensus generation via iVar
- `6_mafft_fasttree.py` — MAFFT alignment and FastTree phylogeny
- `7_snpsites_hcluster.py` — SNP extraction and clustering analyses

## Pipeline order

The intended workflow is sequential:

1. Run Kraken2-based viral read extraction and filtering.
2. Assemble the enriched reads with SPAdes.
3. Choose a suitable reference genome for each sample.
4. Scaffold assembled contigs against the selected reference.
5. Validate assembly quality with CheckV.
6. Polish the scaffold by mapping reads and calling a consensus.
7. Build a whole-genome alignment and phylogenetic tree.
8. Extract SNPs and cluster sequences for epidemiological comparison.

## Typical inputs

This repo assumes the following inputs are available:

- raw paired-end or single-end FASTQ reads
- a Kraken2 database with a viral taxonomy scope
- a prebuilt Mash viral reference database (for the Mash workflow)
- a local BLAST nt database (for BLAST-based reference selection)
- external bioinformatics dependencies installed and callable in PATH or with explicit executable paths in the scripts

## External dependencies

The scripts rely on the following tools, typically installed in a bioinformatics environment:

- Kraken2
- fastp
- SPAdes or rnaviralspades
- Mash
- BLAST+ (`blastn`)
- NCBI Entrez / NCBI datasets access (for genome retrieval)
- MUMmer (`nucmer`, `show-coords`, `delta-filter`, `show-tiling`)
- CheckV
- minimap2
- samtools
- iVar
- MAFFT
- FastTree
- snp-sites
- Biopython

## Important configuration notes

Before running the scripts in a new environment:

- update all hardcoded paths in each script
- confirm the accession list matches the samples you want to process
- verify that the Kraken2 database path and taxonomy scope are appropriate
- confirm that input read directories and output directories exist and are writable
- confirm that binaries are available in PATH or replace the hardcoded executable paths

At present, the scripts are not designed around CLI arguments or configuration files; they are stepwise scripts focused on a specific local workflow.

## Running the workflow

The scripts are designed to be executed in order, usually as direct Python commands:

```bash
python 1_run_kraken2_scrub_human.py
python 2_default_spades_assembly.py
python 3_mash_refgenome_dowload.py
# or
python 3_blast_taxid_refgenome_download.py
python 4_nucmer_scaffold_checkv.py
python 5_ivar_consensus_polishing.py
python 6_mafft_fasttree.py
python 7_snpsites_hcluster.py
```

Because the scripts use hardcoded directories and per-sample accession lists, it is best practice to test one accession first before scaling to the full batch.

## Data expectations and outputs

The outputs are organized by accession directory under the configured output base directories. Typical outputs include:

- cleaned taxonomically filtered read FASTQs
- assembly outputs from SPAdes
- candidate reference genomes
- aligned and scaffolded sequences
- CheckV QC reports
- polished consensus FASTAs
- aligned genomes and phylogenetic tree files
- SNP alignment and clustering outputs

The exact file names are defined in the scripts and may vary slightly depending on which branch of the workflow is used.

## Current caveats

This is still a developmental workflow and not a polished production pipeline. Some expected issues include:

- hardcoded local file paths
- script-specific assumptions about directory naming
- no formal CLI interface
- limited error handling and validation in some steps
- dependence on a particular project layout and sequencing setup
- manual debugging may be required for failed samples or reference selection inconsistencies

## Recommended workflow for future use

For ongoing development, I am planning to:

1. convert each script to a common configuration file or CLI options
2. replace absolute paths with environment variables or config-driven arguments
3. add standard logging and reproducible run metadata
4. add a dry-run mode and sample-level skip logic
5. validate one accession end-to-end before larger batch runs
6. add tests for file discovery, assembly QC, and reference selection logic
7. package the project into a reproducible environment (conda, Docker, or Nextflow/Snakemake)

## Minimal handoff guidance

Immediate tasks:

- the path definitions at the top of each script
- the accession list in each script
- whether the input directory structure still matches the expected naming
- whether the external tools are installed and available
- whether the reference database paths are still valid

The scripts are best treated as a research pipeline template rather than a finished toolchain.

## License

See the repository LICENSE file for licensing details.

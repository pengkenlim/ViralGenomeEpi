"""Operational pipeline configuration.

This file defines the shared settings used by the scripts in pipeline/operational.
The config is imported as a Python module, so each dictionary is read directly by
individual modules.

Key conventions:
- cross_module_params contains values shared across modules.
- module_X_params contains parameters that are specific to each analysis stage.
- The "samples" list is the main set of samples processed by the pipeline.
- Modules 7 and 8 can override that list with their own input_samples when a
  smaller subset should be analyzed as a focused phylogeny/clustering run.
- analysis_name groups related outputs under a subdirectory under output_dir, e.g.
  module_7/RSV_group1 or module_8/RSV_group1.

Edit paths carefully when moving the project to a new machine or changing the
reference databases or software installation locations.
"""

# ---------------------------------------------------------------------------
# Shared settings used across modules
# ---------------------------------------------------------------------------
cross_module_params = {
    "read_dir": "/PATH/TO/REPO/test_data/input",
    # Folder containing sample directories and raw FASTQ files.
    "samples": [
        "SRR34884124",
        "SRR34884125",
        "SRR34884126",
        "SRR34884127",
        "SRR34884128",
        "SRR34884276",
        "SRR34884296",
        "SRR34884297",
        "SRR34884298",
    ],
    # The full sample set for the default pipeline run.
    "output_dir": "",
    # Root directory where module_* outputs are stored.
    "read_1_str": "_1.fastq",
    # File suffix used to detect paired-end read 1 FASTQs.
    "read_2_str": "_2.fastq",
    # File suffix used to detect paired-end read 2 FASTQs.
    "singleton_str": ".fastq",
    # Suffix used to detect singleton reads that do not match read_1_str/read_2_str.
    # Set to None if singleton reads should be ignored.
    "threads": 8,
    # Default thread count used by modules that do not override their own value.
}


# ---------------------------------------------------------------------------
# Module 1: QC and adapter trimming
# ---------------------------------------------------------------------------
module_1_params = {
    "threads": 8,
    # Number of threads for fastp. If set to None, the cross-module value is used.
    "fastp_bin_path": "/PATH/TO/fastp",
    # Path to the fastp executable.
    "multi_qc_bin_path": "multiqc", #installed by pip
    # Path or command name for MultiQC (installed in the Python environment).
}


# ---------------------------------------------------------------------------
# Module 2: Taxonomic classification and read filtering
# ---------------------------------------------------------------------------
module_2_params = {
    "threads": 64,
    # Higher thread count for classification/processing.
    "Kraken2_k2_bin_path": "/PATH/TO/kraken2/k2",
    # Path to the Kraken2 k-mer database executable.
    "kraken2_db_path": "/PATH/TO/k2_standard_20260626/",
    # Directory containing the Kraken2 reference database.
    "scope_to_keep_taxId": 10239,
    # Viral taxonomy ID kept during the classification step (e.g., viruses).
}


# ---------------------------------------------------------------------------
# Module 3: Assembly of viral reads into contigs/scaffolds
# ---------------------------------------------------------------------------
module_3_params = {
    "threads": 64,
    # CPU threads for assembly.
    "memory": 500,
    # Memory allocation in MB for the assembler.
    "rnaviralspades_bin_path": "/PATH/TO/SPAdes-4.3.0-Linux/bin/rnaviralspades.py",
    # Full path to the rnaviralSPAdes script for RNA virus assembly.
    "reads_subset": "viral",
    # Which reads to assemble: usually the viral subset or dominant-virus subset.
    "use_rna_mode": True,
    # Enable RNA-mode assembly behavior in the assembler.
}


# ---------------------------------------------------------------------------
# Module 4: Reference selection and metadata enrichment
# ---------------------------------------------------------------------------
module_4_params = {
    "threads": 32,
    # Threads for BLAST and related database queries.
    "blastn_bin_path": "/PATH/TO/blastn",
    # Executable for nucleotide BLAST.
    "blastdbcmd_bin_path": "/PATH/TO/blastdbcmd",
    # Tool used to fetch sequence records from the local BLAST DB.
    "cd_hit_est_bin_path": "/PATH/TO/cd-hit-est",
    # Sequence clustering tool used during transcript/reference matching.
    "blast_db_path": "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/viral_refseq_blast_db",
    # Local BLAST database of viral reference sequences.
    "seq2genome_map_path": "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/seq2genome_map.json",
    # Mapping from assembled transcript identifiers to full genome accessions.
    "assembly_data_report_path": "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/assembly_data_report.jsonl",
    # Metadata file containing assembly-level genome information.
    "nodes_path": "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/nodes.dmp",
    # NCBI taxonomy nodes table used to map tax IDs to parent/lineage information.
    "names_path": "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/names.dmp",
    # NCBI taxonomy names table used to resolve scientific names and common names.
}


# ---------------------------------------------------------------------------
# Module 5: Scaffold QC and reference coverage validation
# ---------------------------------------------------------------------------
module_5_params = {
    "threads": 64,
    # Threads for mapping / QC workflow.
    "checkv_bin_path": "checkv", # installed by pip
    # CheckV executable used for contamination and completeness assessment.
    "checkv_db_dir": "/PATH/TO/checkv-db-v1.5",
    # Directory containing the CheckV reference database.
    "minimap2_bin_path": "/PATH/TO/minimap2-2.31_x64-linux/minimap2",
    # Minimap2 executable used for read mapping and scaffold validation.
}


# ---------------------------------------------------------------------------
# Module 6: Reference-based polishing and consensus generation
# ---------------------------------------------------------------------------
module_6_params = {
    "threads": 64,
    # Threads for read mapping and iVar polishing.
    "minimap2_bin_path": "/PATH/TO/minimap2-2.31_x64-linux/minimap2",
    # Minimap2 alignment executable.
    "samtools_bin_path": "/PATH/TO/samtools",
    # Samtools executable for BAM indexing and manipulation.
    "ivar_bin_path": "/PATH/TO/ivar",
    # iVar executable used to generate consensus sequences from mapped reads.
    "checkv_bin_path": "checkv", #installed by pip
    # CheckV executable for QC checks on assembled consensus sequences.
    "checkv_db_dir": "/PATH/TO/checkv-db-v1.5",
    # CheckV database directory.
}


# ---------------------------------------------------------------------------
# Module 7: Multiple sequence alignment and phylogeny building
# ---------------------------------------------------------------------------
module_7_params = {
    "threads": 32,
    # Threads for MAFFT alignment.
    "mafft_bin_path": "/PATH/TO/mafft",
    # Path to the MAFFT alignment executable.
    "fasttree_bin_path": "/PATH/TO/FastTree",
    # Path to FastTree used to build the approximate maximum-likelihood tree.
    "input_samples": [
        "SRR34884124",
        "SRR34884125",
        "SRR34884126",
        "SRR34884127",
        "SRR34884128",
    ],
    # Subset of samples to include in this independent phylogeny analysis.
    "analysis_name": "RSV_group1",
    # Output directory name under module_7, e.g. output/module_7/RSV_group1.
}


# ---------------------------------------------------------------------------
# Module 8: SNP-site extraction and clustering analysis
# ---------------------------------------------------------------------------
module_8_params = {
    "threads": 32,
    # Threads for SNP analysis and clustering tasks.
    "snpsites_bin_path": "/PATH/TO/snp-sites",
    # Executable name or full path for the snp-sites binary.
    "snp_cutoff": 3,
    # Maximum pairwise SNP distance used to define the cluster threshold.
    "linkage_method": "single",
    # Hierarchical clustering method used by SciPy (e.g. single, average, complete).
    "input_samples": [
        "SRR34884124",
        "SRR34884125",
        "SRR34884126",
        "SRR34884127",
        "SRR34884128",
    ],
    # Subset of samples to include in this clustering run.
    "analysis_name": "RSV_group1",
    # Output directory name under module_8, e.g. output/module_8/RSV_group1.
}
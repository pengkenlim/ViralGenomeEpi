import os
import subprocess
from pathlib import Path

# User Configuration
spades_path = "/media/mngs/48TBRAID5HDD/viral/SPAdes-4.3.0-Linux/bin/rnaviralspades.py"  # or full path e.g., "/usr/local/bin/spades.py"
read_dir = "/media/mngs/48TBRAID5HDD/viral_3/kraken2_out" # Path where scrubbed fastq files are saved

accession_list = [
    "SRR34884128",#8
    "SRR34884127",#7
    "SRR34884126",#7
    "SRR34884125",#9
    "SRR34884124",
    "SRR34884179",
    "SRR34884178",
    "SRR34884145",
    "SRR34884143",
    "SRR34884142",
    "SRR34884141",
    "SRR34884140",
    "SRR34884139",
    "SRR34884138",
    "SRR34884137",
    "SRR34884136",
    "SRR34884135",
    "SRR34884134",
    "SRR34884131",
    "SRR34884130",
    "SRR34884129",
    "SRR34884123",
    "SRR34884122",
    "SRR34884120",
    "SRR34884119",
    "SRR34884118",
    "SRR34884117",
    "SRR34884116",
    "SRR34884115",
    "SRR34884114",
    "SRR34884241",
    "SRR34884240",
    "SRR34884239",
    "SRR34884237",
    "SRR34884236",
    "SRR34884235",
    "SRR34884234",
    "SRR34884233",
    "SRR34884232",
]
accession_list = ["SRR34884135"]
output_dir = "/media/mngs/48TBRAID5HDD/viral_3/spades_out"
threads = 16  # Adjust CPU threads for SPAdes
memory_gb = 500  # Memory limit in GB
use_rna_mode = True  # Set to True for rSPAdes (--rna flag)



def run_cmd(cmd):
    """Utility to run shell commands safely."""
    print(f"\n[EXEC]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def process_accession(acc):
    print(f"\n==========================================")
    print(f" Running SPAdes Assembly for: {acc}")
    print(f"==========================================")

    # Output directory for this specific sample
    sample_out_dir = Path(output_dir) / acc
    sample_out_dir.mkdir(parents=True, exist_ok=True)

    # Locate input FASTQ files (supports both gzipped and plain FASTQs)
    acc_dir = Path(read_dir) / acc if (Path(read_dir) / acc).is_dir() else Path(read_dir)

    # Use dominant-species reads produced by step 1
    r1 = list(acc_dir.glob(f"{acc}_dominant_species_1_trimmed.fastq*"))
    r2 = list(acc_dir.glob(f"{acc}_dominant_species_2_trimmed.fastq*"))
    u = list(acc_dir.glob(f"{acc}_dominant_species_unpaired_trimmed.fastq*"))

    # Build SPAdes command
    cmd = [
        spades_path,
        "-o",
        str(sample_out_dir),
        "-t",
        str(threads),
        "-m",
        str(memory_gb),
    ]

    # Add RNA mode if specified
    if use_rna_mode:
        cmd.append("--rna")

    has_inputs = False

    # Attach paired-end inputs
    if r1 and r2:
        cmd.extend(["-1", str(r1[0]), "-2", str(r2[0])])
        has_inputs = True
        print(f"[INPUT]: Found paired reads: {r1[0].name}, {r2[0].name}")

    # Attach unpaired input
    if u:
        cmd.extend(["-s", str(u[0])])
        has_inputs = True
        print(f"[INPUT]: Found unpaired reads: {u[0].name}")

    if not has_inputs:
        print(f"[SKIP]: No matching FASTQ files found for {acc} in {acc_dir}")
        return

    # Execute SPAdes
    run_cmd(cmd)
    
    # Summary report on key output
    scaffolds_file = sample_out_dir / "scaffolds.fasta"
    contigs_file = sample_out_dir / "contigs.fasta"
    
    if scaffolds_file.exists():
        print(f"[SUCCESS]: Scaffolds generated at: {scaffolds_file}")
    elif contigs_file.exists():
        print(f"[SUCCESS]: Contigs generated at: {contigs_file}")


def main():
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for acc in accession_list:
        try:
            process_accession(acc)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR]: SPAdes assembly failed for {acc}: {e}")
        except Exception as e:
            print(f"[ERROR]: Unexpected error for {acc}: {e}")


if __name__ == "__main__":
    main()
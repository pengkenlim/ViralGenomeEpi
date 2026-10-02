import os
import glob
import subprocess
import re
import zipfile
import shutil
import json
from pathlib import Path

# =============================================================================
# User Configuration
# =============================================================================

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
enriched_reads_dir = "/media/mngs/48TBRAID5HDD/viral_3/kraken2_out"  # Output directory from step 1 (dominant-species FASTQs)
ref_database_msh = "/media/mngs/48TBRAID5HDD/viral2/refseq.genomes.k21s1000.msh" # Pre-sketched Mash database of viral reference genomes
ref_genome_dir = "/media/mngs/48TBRAID5HDD/viral_3/ref_genome"          # Central directory to store downloaded reference genomes
output_dir = "/media/mngs/48TBRAID5HDD/viral_3/automated_reference_selection"
mash_path = "/media/mngs/48TBRAID5HDD/viral2/mash-Linux64-v2.3/mash"  # Full path to Mash binary; update if not in $PATH
datasets_bin_path = "" # Optional: Supply full path to NCBI 'datasets' binary if needed (e.g., "/path/to/datasets")
threads = 32

# =============================================================================
# Functions
# =============================================================================

def locate_enriched_read_files(acc):
    """Locates paired (R1, R2) and unpaired dominant-species FASTQ files for an accession."""
    acc_dir = os.path.join(enriched_reads_dir, acc)

    r1_files = sorted(glob.glob(os.path.join(acc_dir, f"{acc}_dominant_species_1.fastq*")))
    r2_files = sorted(glob.glob(os.path.join(acc_dir, f"{acc}_dominant_species_2.fastq*")))
    unpaired_files = sorted(glob.glob(os.path.join(acc_dir, f"{acc}_dominant_species_unpaired.fastq*")))

    return r1_files, r2_files, unpaired_files

def extract_clean_accession(ref_header_string):
    """Extracts a standard RefSeq/GenBank assembly accession (e.g., GCF_000855545.1) from a messy header."""
    match = re.search(r'(GC[AF]_\d+\.\d+)', ref_header_string)
    if match:
        return match.group(1)
    return ref_header_string.split(".")[0] + "." + ref_header_string.split(".")[1]

def get_fasta_stats(fasta_path):
    """Calculates number of sequences and total base length from a FASTA file."""
    num_seqs = 0
    total_bases = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                num_seqs += 1
            else:
                total_bases += len(line.strip())
    return num_seqs, total_bases

def fetch_assembly_metadata(assembly_accession):
    """Queries NCBI datasets API report to get taxid and scientific name."""
    try:
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/dataset_report"
        cmd = ["curl", "-s", api_url]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        data = json.loads(result.stdout)
        
        reports = data.get("reports", [])
        if reports:
            organism = reports[0].get("organism", {})
            tax_id = organism.get("tax_id", "N/A")
            sci_name = organism.get("sci_name", "N/A")
            return tax_id, sci_name
    except Exception:
        pass
    return "N/A", "N/A"

def download_ncbi_genome(assembly_accession, target_fasta_path):
    """Downloads a genome fasta using NCBI's 'datasets' tool (via custom path or PATH) or API fallback."""
    if os.path.exists(target_fasta_path) and os.path.getsize(target_fasta_path) > 0:
        print(f"  -> Reference already present locally: {target_fasta_path} (Skipping download)")
        return True

    print(f"  -> Downloading reference {assembly_accession} via NCBI Datasets...")
    os.makedirs(ref_genome_dir, exist_ok=True)
    
    zip_output = os.path.join(ref_genome_dir, f"{assembly_accession}.zip")
    extract_dir = os.path.join(ref_genome_dir, f"{assembly_accession}_extracted")
    
    try:
        bin_to_use = None
        if datasets_bin_path and os.path.exists(datasets_bin_path):
            bin_to_use = datasets_bin_path
        else:
            bin_to_use = shutil.which("datasets")
            
        if bin_to_use:
            cmd = [bin_to_use, "download", "genome", "accession", assembly_accession, "--include", "genome", "--filename", zip_output]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        else:
            api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/download?include_annotation_type=GENOME_FASTA"
            curl_cmd = ["curl", "-s", api_url, "-o", zip_output]
            subprocess.run(curl_cmd, check=True)
            
        if not os.path.exists(zip_output) or os.path.getsize(zip_output) < 100:
            print(f"[ERROR]: Downloaded package for {assembly_accession} is empty or missing.")
            return False
            
        with zipfile.ZipFile(zip_output, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
            
        fna_files = []
        for root, dirs, files in os.walk(extract_dir):
            for file in files:
                if file.endswith(".fna") or file.endswith(".fa") or file.endswith(".fasta"):
                    fna_files.append(os.path.join(root, file))
                    
        if not fna_files:
            print(f"[ERROR]: Could not find genomic FASTA inside downloaded package for {assembly_accession}")
            if os.path.exists(zip_output): os.remove(zip_output)
            if os.path.exists(extract_dir): shutil.rmtree(extract_dir)
            return False
            
        shutil.copy(fna_files[0], target_fasta_path)
        
        if os.path.exists(zip_output): os.remove(zip_output)
        if os.path.exists(extract_dir): shutil.rmtree(extract_dir)
        
        print(f"  -> Successfully downloaded and configured reference {assembly_accession}.")
        return True
        
    except Exception as e:
        print(f"[ERROR]: NCBI Datasets download failed for {assembly_accession}: {e}")
        if os.path.exists(zip_output): os.remove(zip_output)
        if os.path.exists(extract_dir): shutil.rmtree(extract_dir)
        return False

def select_best_reference(acc, acc_out_dir):
    """Sketches reads, queries Mash database, identifies best reference, and handles download & summary reporting."""
    print(f"\n[MASH]: Processing reference selection for {acc}...")
    
    r1_files, r2_files, unpaired_files = locate_enriched_read_files(acc)
    all_fastqs = []
    for r1, r2 in zip(r1_files, r2_files):
        all_fastqs.extend([r1, r2])
    all_fastqs.extend(unpaired_files)
    
    if not all_fastqs:
        print(f"[ERROR]: No enriched FASTQ files found for {acc} in {enriched_reads_dir}")
        return None

    sample_msh = os.path.join(acc_out_dir, f"{acc}_reads.msh")
    cmd_sketch = [mash_path, "sketch", "-m", "2", "-r", "-c", "200", "-o", sample_msh] + all_fastqs
    subprocess.run(cmd_sketch, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    dist_tsv = os.path.join(acc_out_dir, f"{acc}_mash_distances.tsv")
    cmd_dist = [mash_path, "dist", ref_database_msh, sample_msh]
    
    with open(dist_tsv, "w") as dist_f:
        subprocess.run(cmd_dist, stdout=dist_f, stderr=subprocess.DEVNULL, check=True)

    best_ref_id = None
    min_dist = 1.0
    
    with open(dist_tsv, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            ref_name = parts[0]
            distance = float(parts[2])

            if distance < min_dist:
                min_dist = distance
                best_ref_id = ref_name

    if not best_ref_id:
        print(f"[ERROR]: Could not determine best reference match for {acc}.")
        return None

    ref_basename = Path(best_ref_id).name
    assembly_accession = extract_clean_accession(ref_basename)
    
    print(f"[MASH SUCCESS]: Closest reference found: '{ref_basename}' (Clean Accession: {assembly_accession}) with Mash distance {min_dist:.4f}")

    central_ref_fasta = os.path.join(ref_genome_dir, f"{assembly_accession}.fasta")
    
    success = download_ncbi_genome(assembly_accession, central_ref_fasta)
    if not success:
        return None

    # Compute overview metrics requested
    tax_id, sci_name = fetch_assembly_metadata(assembly_accession)
    num_seqs, total_bases = get_fasta_stats(central_ref_fasta)

    # Save summary report to JSON
    summary_data = {
        "sample_accession": acc,
        "best_reference_accession": assembly_accession,
        "scientific_name": sci_name,
        "tax_id": tax_id,
        "pairwise_mash_distance": min_dist,
        "number_of_sequences": num_seqs,
        "total_bases_length": total_bases
    }
    
    summary_json_path = os.path.join(acc_out_dir, f"{acc}_reference_summary.json")
    with open(summary_json_path, "w") as json_f:
        json.dump(summary_data, json_f, indent=4)
        
    print("  -> Reference Summary Recorded:")
    print(f"     • Scientific Name: {sci_name} (TaxID: {tax_id})")
    print(f"     • Pairwise Distance: {min_dist:.5f}")
    print(f"     • Sequences: {num_seqs} | Total Bases: {total_bases:,} bp")
    print(f"     • Saved to: {summary_json_path}")

    sample_ref_symlink = os.path.join(acc_out_dir, f"selected_reference_{assembly_accession}.fasta")
    if os.path.lexists(sample_ref_symlink):
        os.remove(sample_ref_symlink)
    os.symlink(central_ref_fasta, sample_ref_symlink)

    return sample_ref_symlink

def main():
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(ref_genome_dir, exist_ok=True)
    for acc in accession_list:
        acc_out_dir = os.path.join(output_dir, acc)
        os.makedirs(acc_out_dir, exist_ok=True)
        try:
            selected_ref = select_best_reference(acc, acc_out_dir)
            if selected_ref:
                print(f"[READY]: Reference prepared for RagTag -> {selected_ref}\n")
        except Exception as e:
            print(f"[ERROR]: Reference selection failed for {acc}: {e}")

if __name__ == "__main__":
    main()
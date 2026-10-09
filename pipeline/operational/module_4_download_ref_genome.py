import argparse
import csv
import importlib.util
import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path


def load_config(config_path):
    """Load a Python config file as a Python module object."""
    spec = importlib.util.spec_from_file_location("config", config_path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load config from: {config_path}")

    config = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config)
    return config


def get_output_dir(config, sample_name):
    return os.path.join(config.cross_module_params["output_dir"], "module_4", sample_name)

def get_ref_genome_dir(config):
    return os.path.join(config.cross_module_params["output_dir"], "module_4", "ref_genomes")

def _normalize_fasta_header(header, occurrence, total_occurrences):
    """Append a numerical suffix to duplicate FASTA identifiers.

    For repeated headers, the first occurrence becomes <id>_1, the second
    <id>_2, and so on. A single occurrence is left unchanged.
    """
    header = header.strip()
    if total_occurrences <= 1:
        return header

    if " " in header:
        id_part, suffix_part = header.split(None, 1)
        return f"{id_part}_{occurrence} {suffix_part}"
    return f"{header}_{occurrence}"


def get_module_3_dir(config, sample_name):
    fasta_records = []
    seen_counts = {}

    for subdir in ["viral", "viral_over_assembly", "dominant_virus", "dominant_virus_over_assembly", "aligned_to_ref_over_assembly", "aligned_to_ref"]:
        potential_fasta = os.path.join(config.cross_module_params["output_dir"], "module_3", subdir, sample_name, "transcripts.fasta")
        if not os.path.exists(potential_fasta):
            continue

        with open(potential_fasta, "r") as handle:
            fasta_text = handle.read().strip()

        if not fasta_text:
            continue

        current_records = []
        parts = fasta_text.split(">")
        for part in parts:
            if not part.strip():
                continue
            lines = part.strip().splitlines()
            if not lines:
                continue
            header = lines[0].strip()
            seq = "".join(line.strip() for line in lines[1:])
            current_records.append((header, seq))
            header_id = header.split()[0]
            seen_counts[header_id] = seen_counts.get(header_id, 0) + 1

        fasta_records.extend(current_records)

    if fasta_records:
        combined_fasta_path = os.path.join(config.cross_module_params["output_dir"], "module_3", sample_name, "combined_transcripts.fasta")
        os.makedirs(os.path.dirname(combined_fasta_path), exist_ok=True)

        occurrence_counts = {}
        with open(combined_fasta_path, "w") as handle:
            for header, seq in fasta_records:
                header_id = header.split()[0]
                occurrence_counts[header_id] = occurrence_counts.get(header_id, 0) + 1
                final_header = _normalize_fasta_header(header, occurrence_counts[header_id], seen_counts.get(header_id, 1))
                handle.write(f">{final_header}\n{seq}\n")

    return os.path.join(config.cross_module_params["output_dir"], "module_3", sample_name)

def locate_transcripts_file(config, sample_name):
    #locate and copy the assembled transcript FASTA from module 3
    sample_dir = get_module_3_dir(config, sample_name)
    path = os.path.join(sample_dir, "combined_transcripts.fasta")
    if os.path.exists(path):
        return path
    raise FileNotFoundError(f"No assembled transcript FASTA found for sample {sample_name} in {sample_dir}.")

def get_blast_db_path(config):
    return config.module_4_params.get("blast_db_path", "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/viral_refseq_blast_db")

def get_blastdbcmd_bin_path(config):
    return config.module_4_params.get("blastdbcmd_bin_path", "blastdbcmd")

def get_cd_hit_est_bin_path(config):
    return config.module_4_params.get("cd_hit_est_bin_path", "cd-hit-est")

def get_genome_map_path(config):
    return config.module_4_params.get("seq2genome_map_path", "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/seq2genome_map.json")

def get_assembly_report_path(config):
    return config.module_4_params.get("assembly_data_report_path", "/home/mngs/GitHubRepos/ViralGenomeEpi/dependencies/assembly_data_report.jsonl")

def get_nodes_path(config):
    return config.module_4_params.get("nodes_path", "/media/mngs/48TBRAID5HDD/Temp_output_ken/taxdb_data/nodes.dmp")

def get_names_path(config):
    return config.module_4_params.get("names_path", "/media/mngs/48TBRAID5HDD/Temp_output_ken/taxdb_data/names.dmp")


# --- Taxonomy Parsing Functions ---
def read_taxonomy_nodes(nodes_path):
    nodes = {}
    if not nodes_path or not os.path.exists(nodes_path): return nodes
    with open(nodes_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or "\t|" not in line: continue
            fields = [part.strip() for part in line.split("\t|")]
            if len(fields) < 3: continue
            try:
                taxid_int = int(fields[0]); parent_int = int(fields[1])
                nodes[taxid_int] = {"parent_taxid": parent_int, "rank": fields[2]}
            except ValueError: continue
    return nodes

def read_taxonomy_names(names_path):
    names = {}
    if not names_path or not os.path.exists(names_path): return names
    with open(names_path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or "\t|" not in line: continue
            fields = [part.strip() for part in line.split("\t|")]
            if len(fields) < 3: continue
            try: taxid_int = int(fields[0])
            except ValueError: continue
            if taxid_int not in names: names[taxid_int] = {}
            if fields[3] if len(fields) > 3 else "" == "scientific name": names[taxid_int]["scientific_name"] = fields[1]
            if fields[3] if len(fields) > 3 else "" == "genbank common name": names[taxid_int]["common_name"] = fields[1]
            if "canonical_name" not in names[taxid_int]: names[taxid_int]["canonical_name"] = fields[1]
    return names

def find_species_taxid_for_organism_taxid(organism_taxid, nodes_path, names_path):
    if organism_taxid is None: return None, None
    nodes = read_taxonomy_nodes(nodes_path); names = read_taxonomy_names(names_path)
    current_taxid = int(organism_taxid); seen = set()
    while current_taxid not in seen:
        seen.add(current_taxid); node = nodes.get(current_taxid)
        if node and node.get("rank") == "species":
            tax_name = names.get(current_taxid, {}).get("scientific_name") or names.get(current_taxid, {}).get("canonical_name")
            return current_taxid, tax_name
        if node is None: break
        current_taxid = int(node.get("parent_taxid", -1))
        if current_taxid <= 0: break
    return None, None

def normalize_species_metadata(genome_accession, report_path, nodes_path=None, names_path=None):
    organism_taxid, organism_name = None, None
    if os.path.exists(report_path):
        with open(report_path, "r") as handle:
            for line in handle:
                line = line.strip()
                if not line: continue
                try: record = json.loads(line)
                except json.JSONDecodeError: continue
                accession = record.get("currentAccession") or record.get("accession")
                if accession and accession.lower() == str(genome_accession).lower():
                    organism = record.get("organism", {})
                    organism_taxid = organism.get("taxId"); organism_name = organism.get("organismName")
                    break
    if organism_taxid is not None:
        species_taxid, species_sci_name = find_species_taxid_for_organism_taxid(organism_taxid, nodes_path, names_path)
        if species_taxid is not None and species_sci_name is not None:
            return organism_taxid, organism_name, species_taxid, species_sci_name
    return organism_taxid, organism_name, None, None


# --- Utility Functions ---
def get_fasta_stats(fasta_path):
    num_seqs = 0; total_bases = 0
    with open(fasta_path, "r") as handle:
        for line in handle:
            if line.startswith(">"): num_seqs += 1
            else: total_bases += len(line.strip())
    return num_seqs, total_bases

def detect_segmented_virus(num_seqs):
    return num_seqs > 1

def extract_clean_accession(ref_header_string):
    cleaned = str(ref_header_string).strip()
    match = re.search(r"(NC_[0-9]+\.[0-9]+|NG_[0-9]+\.[0-9]+|NT_[0-9]+\.[0-9]+)", cleaned)
    return match.group(1) if match else cleaned

def find_genome_accession(ref_accession, seq2genome_map_path):
    if not os.path.exists(seq2genome_map_path): raise FileNotFoundError(f"Genome mapping file not found: {seq2genome_map_path}")
    with open(seq2genome_map_path, "r") as handle: genome_map = json.load(handle)
    if ref_accession in genome_map: return genome_map[ref_accession]
    ref_basename = ref_accession.split(".")[0]
    for key, value in genome_map.items():
        if key.split(".")[0] == ref_basename: return value
    raise KeyError(f"No genome accession found for reference accession {ref_accession}")


def download_ncbi_genome(assembly_accession, target_fasta_path, config):
    if os.path.exists(target_fasta_path) and os.path.getsize(target_fasta_path) > 0:
        print(f"  -> Reference already present locally: {target_fasta_path} (skipping download)"); return True
    ref_genome_dir = get_ref_genome_dir(config); os.makedirs(ref_genome_dir, exist_ok=True)
    zip_output = os.path.join(ref_genome_dir, f"{assembly_accession}.zip"); extract_dir = os.path.join(ref_genome_dir, f"{assembly_accession}_extracted")
    try:
        api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/download?include_annotation_type=GENOME_FASTA"
        subprocess.run(["curl", "-fsSL", api_url, "-o", zip_output], check=True)
        if not os.path.exists(zip_output) or os.path.getsize(zip_output) < 100: return False
        with zipfile.ZipFile(zip_output, "r") as zip_ref: zip_ref.extractall(extract_dir)
        fasta_files = [os.path.join(root, filename) for root, _, files in os.walk(extract_dir) for filename in files if filename.endswith((".fna", ".fa", ".fasta"))]
        if not fasta_files: return False
        fasta_files.sort(); shutil.copy(fasta_files[0], target_fasta_path)
        if os.path.exists(zip_output): os.remove(zip_output)
        if os.path.exists(extract_dir): shutil.rmtree(extract_dir)
        return True
    except Exception as exc:
        print(f"[ERROR]: NCBI download failed for {assembly_accession}: {exc}")
        if os.path.exists(zip_output): os.remove(zip_output)
        if os.path.exists(extract_dir): shutil.rmtree(extract_dir)
        return False


# --- BLAST Hit Extraction and CD-HIT Clustering ---

def extract_blast_hit_fastas(config, blast_tsv_path, blast_db_path, output_fasta_path):
    """Extract the actual FASTA sequences of the passing BLAST hits from the local DB."""
    sseqids = set()
    with open(blast_tsv_path, "r") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                try:
                    pident = float(parts[2]); length = int(parts[3])
                    if length >= 500 and pident >= 70.0:
                        # Keep ORIGINAL ID for blastdbcmd lookup
                        sseqids.add(parts[1].strip())
                except ValueError: continue

    if not sseqids:
        raise ValueError("No valid BLAST hits to extract.")

    ids_file = output_fasta_path + ".ids"
    with open(ids_file, "w") as f:
        f.write("\n".join(sseqids))

    blastdbcmd_bin = get_blastdbcmd_bin_path(config)
    cmd = [blastdbcmd_bin, "-db", blast_db_path, "-entry_batch", ids_file, "-out", output_fasta_path]
    
    # UNMUTED for debugging: If blastdbcmd fails (e.g., missing -parse_seqids), you will see the error now.
    subprocess.run(cmd, check=True)

    if os.path.exists(ids_file): os.remove(ids_file)
    return output_fasta_path


def run_cd_hit_clustering(config, input_fasta, output_fasta, identity=0.80):
    """Cluster the extracted BLAST hits using CD-HIT."""
    cd_hit_est_bin = get_cd_hit_est_bin_path(config)
    cmd = [cd_hit_est_bin, "-i", input_fasta, "-o", output_fasta, "-c", str(identity), "-n", "8", "-d", "0"]
    
    # UNMUTED for debugging: If cd-hit-est fails or is not found, you will see the error now.
    subprocess.run(cmd, check=True)
    
    return output_fasta + ".clstr"


def parse_cd_hit_clusters_and_select_best(clstr_file, blast_tsv_path, seq2genome_map_path):
    """
    Parse the CD-HIT .clstr file. For each cluster, sum the bitscores of its members 
    and select the single best representative (highest bitscore).
    """
    genome_map = {}
    if seq2genome_map_path and os.path.exists(seq2genome_map_path):
        with open(seq2genome_map_path, "r") as handle: 
            genome_map = json.load(handle)

    # 1. Parse BLAST TSV to get bitscores per sseqid
    sseqid_bitscores = {}
    with open(blast_tsv_path, "r") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                try:
                    raw_sseqid = parts[1].strip()
                    pident = float(parts[2])
                    length = int(parts[3])
                    bitscore = float(parts[5])
                    
                    # Clean the ID to match the genome_map.json format
                    clean_sseqid = raw_sseqid.split("|")[1]
                    
                    if length >= 500 and pident >= 70.0:
                        sseqid_bitscores[clean_sseqid] = sseqid_bitscores.get(clean_sseqid, 0.0) + bitscore
                except ValueError:
                    continue

    # 2. Parse .clstr file to get clusters
    clusters = []
    current_cluster = []
    
    if not os.path.exists(clstr_file):
        raise FileNotFoundError(f"CD-HIT cluster file not found: {clstr_file}. CD-HIT may have failed or the input FASTA was empty.")

    with open(clstr_file, "r") as handle:
        for line in handle:
            if line.startswith(">Cluster"):
                if current_cluster: 
                    clusters.append(current_cluster)
                current_cluster = []
            else:
                # FIXED REGEX: Captures the ID strictly BEFORE the '...' that CD-HIT appends.
                # It allows dots inside the accession (like NC_007366.1) but stops at '...'
                match = re.search(r">([^\s]+?)(?:\.\.\.|\s|$)", line)
                if match:
                    raw_seq_id = match.group(1)
                    # Apply the exact same cleaning to ensure it matches the BLAST dictionary keys
                    clean_seq_id = raw_seq_id.split("|")[1]
                    current_cluster.append(clean_seq_id)
                    
    if current_cluster:
        clusters.append(current_cluster)

    # 3. Select best representative per cluster
    selected_segments = []
    for cluster in clusters:
        # Filter out any seq_ids that didn't make it into sseqid_bitscores
        valid_cluster = [seq for seq in cluster if seq in sseqid_bitscores]
        if not valid_cluster:
            # Optional: Print debug info if a cluster fails to map, so you know why
            # print(f"[DEBUG] Cluster failed to map. Keys in dict: {list(sseqid_bitscores.keys())[:3]}... Cluster IDs: {cluster[:3]}")
            continue
            
        best_sseqid = max(valid_cluster, key=lambda x: sseqid_bitscores.get(x, 0.0))
        best_score = sseqid_bitscores.get(best_sseqid, 0.0)
        
        # Map back to genome accession
        genome_id = genome_map.get(best_sseqid)
        if genome_id is None:
            ref_basename = best_sseqid.split(".")[0]
            for key, value in genome_map.items():
                if key.split(".")[0] == ref_basename:
                    genome_id = value
                    break
                    
        selected_segments.append({
            "segment_seq_id": best_sseqid,
            "genome_accession": genome_id or "UNKNOWN",
            "bitscore": best_score
        })
        
    selected_segments.sort(key=lambda x: x["bitscore"], reverse=True)
    print(f"[CD-HIT CLUSTERING]: Identified {len(clusters)} distinct segment homology groups.")
    print(f"[CD-HIT CLUSTERING]: Successfully mapped {len(selected_segments)} segments to reference genomes.")
    return selected_segments


def extract_specific_sequence(genome_fasta_path, target_seq_id, output_fasta_path):
    extracted = False
    with open(genome_fasta_path, "r") as in_handle, open(output_fasta_path, "w") as out_handle:
        write_mode = False
        for line in in_handle:
            if line.startswith(">"):
                if target_seq_id in line: write_mode = True; out_handle.write(line)
                else: write_mode = False
            elif write_mode: out_handle.write(line); extracted = True
    if not extracted: raise RuntimeError(f"Could not find sequence {target_seq_id} in {genome_fasta_path}.")
    return output_fasta_path


# --- Main Pipeline Logic ---

def select_best_reference_by_aggregation(blast_tsv_path, seq2genome_map_path=None):
    genome_scores = {}; genome_map = {}
    if seq2genome_map_path and os.path.exists(seq2genome_map_path):
        with open(seq2genome_map_path, "r") as handle: genome_map = json.load(handle)
    with open(blast_tsv_path, "r") as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            try:
                raw_ref_id = parts[1].strip()
                clean_ref_id = raw_ref_id.replace("ref", "").replace("|", "")
                pident = float(parts[2]); length = int(parts[3]); bitscore = float(parts[5])
            except ValueError: continue
            if length < 500 or pident < 70.0: continue
            
            genome_id = clean_ref_id
            if seq2genome_map_path:
                genome_id = genome_map.get(clean_ref_id)
                if genome_id is None:
                    ref_basename = clean_ref_id.split(".")[0]
                    for key, value in genome_map.items():
                        if key.split(".")[0] == ref_basename: genome_id = value; break
                if genome_id is None: genome_id = clean_ref_id
                
            genome_scores[genome_id] = genome_scores.get(genome_id, 0.0) + bitscore
            
    if not genome_scores: raise ValueError("No high-quality matches found in BLAST results.")
    best_genome_id = max(genome_scores, key=genome_scores.get)
    print(f"[SELECTION]: Best aggregate genome is {best_genome_id} with {genome_scores[best_genome_id]:.2f} total bitscore.")
    return best_genome_id


def select_best_reference_for_sample(sample_name, config):
    output_dir = get_output_dir(config, sample_name); os.makedirs(output_dir, exist_ok=True)
    transcripts_path = locate_transcripts_file(config, sample_name)
    blast_db_path = get_blast_db_path(config)
    if not os.path.exists(f"{blast_db_path}.nin"): raise FileNotFoundError(f"BLAST database not found at {blast_db_path}")

    blast_tsv_path = os.path.join(output_dir, f"{sample_name}_blastn_viral_refseq.tsv")
    blastn_bin = config.module_4_params.get("blastn_bin_path", "blastn")
    threads = int(config.module_4_params.get("threads", 8))

    print(f"\n[BLASTN]: Processing reference selection for {sample_name}...")
    cmd = [blastn_bin, "-query", transcripts_path, "-db", blast_db_path, "-out", blast_tsv_path,
           "-outfmt", "6 qseqid sseqid pident length evalue bitscore", "-max_target_seqs", "20",
           "-evalue", "1e-5", "-num_threads", str(threads), "-max_hsps", "1"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    genome_map_path = get_genome_map_path(config)
    genome_accession = select_best_reference_by_aggregation(blast_tsv_path, genome_map_path)
    
    temp_fasta = os.path.join(get_ref_genome_dir(config), f"temp_{genome_accession}.fasta")
    download_ncbi_genome(genome_accession, temp_fasta, config)
    num_seqs, total_bases = get_fasta_stats(temp_fasta)
    if os.path.exists(temp_fasta): os.remove(temp_fasta)

    is_segmented = detect_segmented_virus(num_seqs)
    assembly_report_path = get_assembly_report_path(config)
    nodes_path = get_nodes_path(config); names_path = get_names_path(config)
    organism_taxid, organism_name, species_taxid, species_sci_name = normalize_species_metadata(
        genome_accession, assembly_report_path, nodes_path=nodes_path, names_path=names_path)

    summary_data = {
        "sample_accession": sample_name, "best_reference_accession": genome_accession,
        "genome_accession": genome_accession, "organism_taxid": organism_taxid,
        "organism_name": organism_name, "species_taxid": species_taxid,
        "species_name": species_sci_name, "species_sci_name": species_sci_name,
        "number_of_sequences": num_seqs, "total_bases_length": total_bases,
        "blast_db": blast_db_path, "is_segmented": is_segmented
    }

    if is_segmented:
        print(f"  -> [SEGMENTED DETECTED]: Reference has {num_seqs} sequences. Switching to CD-HIT homology clustering.")
        
        # 1. Extract FASTAs of passing BLAST hits
        hits_fasta = os.path.join(output_dir, f"{sample_name}_blast_hits.fasta")
        extract_blast_hit_fastas(config, blast_tsv_path, blast_db_path, hits_fasta)

        # 2. Run CD-HIT
        clustered_fasta = os.path.join(output_dir, f"{sample_name}_clustered_hits.fasta")
        clstr_file = run_cd_hit_clustering(config, hits_fasta, clustered_fasta, identity=0.80)
        
        # 3. Select best representative per cluster
        selected_segments = parse_cd_hit_clusters_and_select_best(clstr_file, blast_tsv_path, genome_map_path)
        segment_fasta_paths = []
        
        for seg in selected_segments:
            seq_id = seg["segment_seq_id"]; seg_genome = seg["genome_accession"]
            print(f"  -> Selecting segment: {seq_id} (Parent Genome: {seg_genome})")
            target_fasta = os.path.join(get_ref_genome_dir(config), f"{seq_id}.fasta")
            parent_genome_fasta = os.path.join(get_ref_genome_dir(config), f"{seg_genome}.fasta")
            download_ncbi_genome(seg_genome, parent_genome_fasta, config)
            extract_specific_sequence(parent_genome_fasta, seq_id, target_fasta)
            segment_fasta_paths.append({
                "segment_seq_id": seq_id, "genome_accession": seg_genome,
                "fasta_path": target_fasta, "bitscore": seg["bitscore"]
            })
            
        summary_data["selected_segments"] = segment_fasta_paths
        # Cleanup temp files
        #for f in [hits_fasta, clustered_fasta, clstr_file]:
        #    if os.path.exists(f): os.remove(f)
            
    else:
        target_fasta = os.path.join(get_ref_genome_dir(config), f"{genome_accession}.fasta")
        if not download_ncbi_genome(genome_accession, target_fasta, config):
            raise RuntimeError(f"Failed to download genome FASTA for sample {sample_name}.")
        summary_data["selected_segments"] = []

    summary_json = os.path.join(output_dir, f"{sample_name}_reference_summary.json")
    with open(summary_json, "w") as handle: json.dump(summary_data, handle, indent=4)

    primary_symlink = os.path.join(output_dir, f"selected_reference_{genome_accession}.fasta")
    if os.path.lexists(primary_symlink): os.remove(primary_symlink)
    if is_segmented and summary_data["selected_segments"]:
        os.symlink(summary_data["selected_segments"][0]["fasta_path"], primary_symlink)
    else:
        os.symlink(os.path.join(get_ref_genome_dir(config), f"{genome_accession}.fasta"), primary_symlink)

    print(f"     • Selected primary genome accession: {genome_accession}")
    print(f"     • Is Segmented: {is_segmented}")
    print(f"     • Saved to: {summary_json}")
    return primary_symlink


def write_identified_species_summary(config, summary_rows):
    module_4_dir = os.path.join(config.cross_module_params["output_dir"], "module_4")
    os.makedirs(module_4_dir, exist_ok=True)
    output_path = os.path.join(module_4_dir, "Identified_species.tsv")
    fieldnames = ["sample_accession", "best_reference_accession", "genome_accession", "organism_taxid",
                  "organism_name", "species_taxid", "species_name", "species_sci_name",
                  "number_of_sequences", "total_bases_length", "blast_db", "is_segmented", "segment_count"]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in summary_rows:
            row_to_write = row.copy()
            segs = row_to_write.get("selected_segments", [])
            row_to_write["segment_count"] = len(segs) if isinstance(segs, list) else 0
            if "selected_segments" in row_to_write: del row_to_write["selected_segments"]
            writer.writerow({key: row_to_write.get(key, "") for key in fieldnames})
    print(f"[SUMMARY]: Wrote identified species table to {output_path}")


def validate_inputs(config):
    for sample_name in config.cross_module_params["samples"]:
        locate_transcripts_file(config, sample_name)
        print(f"[VALIDATION]: Sample {sample_name} has assembled transcripts available for module 4.")
    print("All module 3 transcript assemblies validated successfully.")


def main():
    parser = argparse.ArgumentParser(description="Select sample-specific reference genomes from assembled transcripts using BLAST.")
    parser.add_argument("--config", required=True, help="Path to the config Python file.")
    args = parser.parse_args()
    config = load_config(args.config)
    print(f"\nStarting module_4_download_ref_genome...\nConfig file location: {args.config}\n")
    validate_inputs(config)
    summary_rows = []
    for sample_name in config.cross_module_params["samples"]:
        try:
            selected_ref = select_best_reference_for_sample(sample_name, config)
            summary_json = os.path.join(get_output_dir(config, sample_name), f"{sample_name}_reference_summary.json")
            if os.path.exists(summary_json):
                with open(summary_json, "r") as handle: summary_rows.append(json.load(handle))
            print(f"[READY]: Reference prepared for sample {sample_name} -> {selected_ref}\n")
        except Exception as exc: 
            print(f"[ERROR]: Failed for sample {sample_name}: {exc}")
    write_identified_species_summary(config, summary_rows)
    print("\n[module_4_download_ref_genome] Finished all samples.")


if __name__ == "__main__":
    main()
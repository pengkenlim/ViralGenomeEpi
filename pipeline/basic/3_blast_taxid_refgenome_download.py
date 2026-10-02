import csv
import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import Entrez, SeqIO

# =============================================================================
# User Configuration
# =============================================================================

Entrez.email = "pengkenlim.sbs@gmail.com"  # REQUIRED by NCBI API
blastn_path = "blastn"  # or full path to binary, e.g. "/path/to/ncbi-blast-16.0.0+/bin/blastn"
read_dir ="/media/mngs/48TBRAID5HDD/viral_3/spades_out" 
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
local_nt_db = "/media/mngs/48TBRAID5HDD/NHS_RMg_platform/db/blastdb/core_nt"
ref_genome_dir = "/media/mngs/48TBRAID5HDD/viral_3/ref_genome"
output_dir = "/media/mngs/48TBRAID5HDD/viral_3/automated_reference_selection_blast"
threads = 96

# =============================================================================
# Helper functions
# =============================================================================


def extract_longest_sequence(transcripts_fasta, output_single_fasta):
    """Extract the longest transcript from a transcript FASTA."""
    longest_seq = None
    max_len = 0

    for record in SeqIO.parse(transcripts_fasta, "fasta"):
        if len(record.seq) > max_len:
            max_len = len(record.seq)
            longest_seq = record

    if longest_seq:
        SeqIO.write(longest_seq, output_single_fasta, "fasta")
        return longest_seq.id, max_len
    return None, 0


def get_fasta_stats(fasta_path):
    """Return sequence count and total bases for a FASTA."""
    num_seqs = 0
    total_bases = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                num_seqs += 1
            else:
                total_bases += len(line.strip())
    return num_seqs, total_bases


def download_ref_genome_for_taxid(species_taxid, target_fasta_path):
    """Download a representative genome FASTA for the BLAST target TaxID."""
    os.makedirs(ref_genome_dir, exist_ok=True)

    if os.path.exists(target_fasta_path) and os.path.getsize(target_fasta_path) > 0:
        print(f"[CACHE]: Reference genome for TaxID {species_taxid} already exists.")
        return True

    print(f"[DOWNLOAD]: Fetching representative genome for TaxID {species_taxid}...")

    query = f"txid{species_taxid}[Organism:exp] AND refseq[filter]"
    try:
        handle = Entrez.esearch(db="nuccore", term=query, retmax=5)
        search_results = Entrez.read(handle)
        handle.close()
        id_list = search_results.get("IdList", [])

        if not id_list:
            query_fallback = f"txid{species_taxid}[Organism:exp]"
            handle = Entrez.esearch(db="nuccore", term=query_fallback, retmax=5)
            search_results = Entrez.read(handle)
            handle.close()
            id_list = search_results.get("IdList", [])

        if not id_list:
            print(f"[ERROR]: No genome sequences found on NCBI for TaxID {species_taxid}.")
            return False

        target_id = id_list[0]
        fetch_handle = Entrez.efetch(
            db="nuccore", id=target_id, rettype="fasta", retmode="text"
        )
        seq_data = fetch_handle.read()
        fetch_handle.close()

        with open(target_fasta_path, "w") as out:
            out.write(seq_data)

        print(f"[SUCCESS]: Saved reference genome to {target_fasta_path}")
        return True
    except Exception as e:
        print(f"[ERROR]: Download failed for TaxID {species_taxid}: {e}")
        return False


def select_best_reference(acc, acc_out_dir):
    """BLAST the longest transcript, pick the strongest target TaxID, and download that TaxID's representative genome."""
    print(f"\n[BLAST-TAXID]: Processing reference selection for {acc}...")

    transcripts_file = os.path.join(read_dir, acc, "transcripts.fasta")
    if not os.path.exists(transcripts_file):
        print(f"[ERROR]: File not found: {transcripts_file}. Skipping.")
        return None

    temp_query_fasta = os.path.join(acc_out_dir, f"{acc}_longest_transcript.fasta")
    seq_id, seq_len_bp = extract_longest_sequence(transcripts_file, temp_query_fasta)
    if not seq_id:
        print(f"[ERROR]: No valid sequences found in {transcripts_file}.")
        return None

    seq_len_kb = seq_len_bp / 1000.0
    print(f"[LONGEST SEQUENCE]: Selected '{seq_id}'")
    print(f"[SEQUENCE SIZE]: {seq_len_bp} bp ({seq_len_kb:.2f} kb)")
    print(f"[EXEC]: Launching megablast for {acc} ({seq_len_kb:.2f} kb contig)...")

    blast_out = os.path.join(acc_out_dir, f"{acc}_blast_hits.tsv")
    blast_cmd = [
        blastn_path,
        "-task",
        "megablast",
        "-query",
        temp_query_fasta,
        "-db",
        local_nt_db,
        "-outfmt",
        "6 qseqid sseqid pident length staxids sscinames bitscore evalue",
        "-max_target_seqs",
        "10",
        "-max_hsps",
        "1",
        "-num_threads",
        str(threads),
        "-out",
        blast_out,
    ]
    subprocess.run(blast_cmd, check=True)

    taxid_votes = defaultdict(float)
    taxid_hit_map = defaultdict(list)

    with open(blast_out, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            (
                qseqid,
                sseqid,
                pident,
                length,
                staxid_raw,
                ssciname,
                bitscore,
                evalue,
            ) = row

            staxid = staxid_raw.split(";")[0].strip()
            if not staxid:
                continue

            bitscore_value = float(bitscore)
            taxid_votes[staxid] += bitscore_value
            taxid_hit_map[staxid].append(
                {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "pident": pident,
                    "length": length,
                    "staxid": staxid,
                    "ssciname": ssciname,
                    "bitscore": bitscore_value,
                    "evalue": evalue,
                }
            )

    if os.path.exists(temp_query_fasta):
        os.remove(temp_query_fasta)

    if not taxid_votes:
        print(f"[WARNING]: No BLAST hits generated for {acc}.")
        return None

    selected_taxid = max(taxid_votes, key=lambda k: taxid_votes[k])
    best_hit = max(taxid_hit_map[selected_taxid], key=lambda x: x["bitscore"])
    selected_species = best_hit.get("ssciname", "Unknown")

    print(
        f"[BLAST TARGET TAXID]: {selected_taxid} ({selected_species}) "
        f"selected from the strongest weighted BLAST target hits"
    )
    print(
        f"[TOP HIT]: {best_hit['sseqid']} | identity={best_hit['pident']}% | "
        f"bitscore={best_hit['bitscore']} | evalue={best_hit['evalue']}"
    )

    central_ref_fasta = os.path.join(ref_genome_dir, f"taxid_{selected_taxid}.fasta")
    success = download_ref_genome_for_taxid(selected_taxid, central_ref_fasta)
    if not success:
        return None

    num_seqs, total_bases = get_fasta_stats(central_ref_fasta)

    summary_data = {
        "sample_accession": acc,
        "best_reference_accession": selected_taxid,
        "scientific_name": selected_species,
        "tax_id": selected_taxid,
        "pairwise_mash_distance": None,
        "number_of_sequences": num_seqs,
        "total_bases_length": total_bases,
        "selected_by": "blast_target_taxid",
    }
    summary_json_path = os.path.join(acc_out_dir, f"{acc}_reference_summary.json")
    with open(summary_json_path, "w") as json_f:
        json.dump(summary_data, json_f, indent=4)

    print("  -> Reference Summary Recorded:")
    print(f"     • Scientific Name: {selected_species} (TaxID: {selected_taxid})")
    print(f"     • Sequences: {num_seqs} | Total Bases: {total_bases:,} bp")
    print(f"     • Saved to: {summary_json_path}")

    sample_ref_symlink = os.path.join(acc_out_dir, f"selected_reference_taxid_{selected_taxid}.fasta")
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

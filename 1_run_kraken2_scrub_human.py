import gzip
import os
import subprocess
from collections import defaultdict
from pathlib import Path

# User Configuration
k2_path = "/home/mngs/GitHubRepos/PrepareMG/dependencies/kraken2/k2"
kraken2_db_path = (
    "/media/mngs/48TBRAID5HDD/AMR/kraken2_dbs/k2_standard_20260626"
)
use_k2_daemon = True
read_dir = "/media/mngs/48TBRAID5HDD/viral/SRR_fastq_data"
viral_scope_taxid = (
    "10239"  # NCBI TaxID for Viruses; species must fall within this taxonomic scope
)
fastp_path = "/media/mngs/48TBRAID5HDD/viral2/fastp"  # Path to fastp executable
fastp_timeout = 1800  # Max seconds to wait for a single fastp run before aborting
fastp_retries = 3  # Retry fastp after a timeout/crash before failing the accession
output_dir = "/media/mngs/48TBRAID5HDD/viral_3/kraken2_out"

threads = 32  # Adjust CPU threads for Kraken2

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

def open_fq(filepath, mode):
    return (
        gzip.open(filepath, mode)
        if str(filepath).endswith(".gz")
        else open(filepath, mode)
    )


def parse_kraken_report(report_path):
    """Parse a Kraken2 report file into a list of entry dicts."""
    entries = []
    with open(report_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            name = parts[5]
            depth = len(name) - len(name.lstrip())
            entries.append(
                {
                    "percent": parts[0].strip(),
                    "reads_covered": int(parts[1].strip()),
                    "reads_assigned": int(parts[2].strip()),
                    "rank": parts[3].strip(),
                    "taxid": parts[4].strip(),
                    "name": name.strip(),
                    "depth": depth,
                }
            )
    return entries


def find_dominant_viral_species(report_paths, viral_scope_taxid):
    """Identify the viral species with the most lineage reads summed across reports."""
    species_counts = defaultdict(int)
    species_names = {}

    for report_path in report_paths:
        entries = parse_kraken_report(report_path)
        viral_idx = next(
            (i for i, e in enumerate(entries) if e["taxid"] == viral_scope_taxid),
            None,
        )
        if viral_idx is None:
            continue
        viral_depth = entries[viral_idx]["depth"]
        i = viral_idx + 1
        while i < len(entries) and entries[i]["depth"] > viral_depth:
            e = entries[i]
            if e["rank"] == "S" and e["reads_covered"] > 0:
                species_counts[e["taxid"]] += e["reads_covered"]
                species_names[e["taxid"]] = e["name"]
            i += 1

    if not species_counts:
        return None

    dom_taxid = max(species_counts, key=lambda t: (species_counts[t], t))
    return {
        "taxid": dom_taxid,
        "name": species_names[dom_taxid],
        "reads": species_counts[dom_taxid],
    }


def get_allowed_taxids(report_path, viral_scope_taxid, dominant_taxid):
    """Return taxids on the path from (but not including) viral scope to the
    dominant species, plus all descendants of the dominant species.
    """
    entries = parse_kraken_report(report_path)
    viral_idx = next(
        (i for i, e in enumerate(entries) if e["taxid"] == viral_scope_taxid),
        None,
    )
    dom_idx = next(
        (i for i, e in enumerate(entries) if e["taxid"] == dominant_taxid), None
    )
    if viral_idx is None or dom_idx is None:
        return set()

    viral_depth = entries[viral_idx]["depth"]
    dom_depth = entries[dom_idx]["depth"]

    allowed = set()

    # Ancestors of the dominant species below the viral scope (inclusive of dominant)
    current_depth = dom_depth
    for i in range(dom_idx, viral_idx - 1, -1):
        e = entries[i]
        if e["depth"] < current_depth:
            if e["depth"] > viral_depth:
                allowed.add(e["taxid"])
            current_depth = e["depth"]
            if current_depth <= viral_depth:
                break

    # Descendants of the dominant species
    i = dom_idx + 1
    while i < len(entries) and entries[i]["depth"] > dom_depth:
        allowed.add(entries[i]["taxid"])
        i += 1

    return allowed


def extract_dominant_species_reads(
    kraken_out_file,
    fastq1_in,
    fastq2_in=None,
    fastq1_out=None,
    fastq2_out=None,
    allowed_taxids=None,
):
    """Extract reads that are unclassified or assigned to the dominant species lineage.

    Handles both paired-end and single-end FASTQs. Input FASTQs may be gzipped.
    """
    if allowed_taxids is None:
        allowed_taxids = set()

    print(f"[DOMINANT SPECIES EXTRACTION]: Parsing {kraken_out_file}...")

    keep_read_ids = set()
    with open(kraken_out_file, "r") as kfile:
        for line in kfile:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                status, read_id, taxid = parts[0], parts[1], parts[2]
                if status == "U" or taxid in allowed_taxids:
                    keep_read_ids.add(read_id)

    print(
        f"[DOMINANT SPECIES EXTRACTION]: Keeping {len(keep_read_ids)} read IDs."
    )

    def filter_single_fastq(in_path, out_path):
        count_kept = 0
        with open_fq(in_path, "rt") as infile, open(out_path, "w") as outfile:
            while True:
                header = infile.readline()
                if not header:
                    break
                seq = infile.readline()
                strand = infile.readline()
                qual = infile.readline()

                read_id = header[1:].split()[0].rstrip("/1").rstrip("/2")

                if read_id in keep_read_ids:
                    outfile.write(f"{header}{seq}{strand}{qual}")
                    count_kept += 1
        return count_kept

    if fastq2_in and fastq2_out:
        count_kept = 0
        with (
            open_fq(fastq1_in, "rt") as in1,
            open_fq(fastq2_in, "rt") as in2,
            open(fastq1_out, "w") as out1,
            open(fastq2_out, "w") as out2,
        ):
            while True:
                h1, s1, st1, q1 = (
                    in1.readline(),
                    in1.readline(),
                    in1.readline(),
                    in1.readline(),
                )
                h2, s2, st2, q2 = (
                    in2.readline(),
                    in2.readline(),
                    in2.readline(),
                    in2.readline(),
                )
                if not h1 or not h2:
                    break

                read_id = h1[1:].split()[0].rstrip("/1").rstrip("/2")

                if read_id in keep_read_ids:
                    out1.write(f"{h1}{s1}{st1}{q1}")
                    out2.write(f"{h2}{s2}{st2}{q2}")
                    count_kept += 1
        print(f"[DOMINANT SPECIES EXTRACTION]: Kept {count_kept} read pairs.")
    else:
        count_kept = filter_single_fastq(fastq1_in, fastq1_out)
        print(
            f"[DOMINANT SPECIES EXTRACTION]: Kept {count_kept} unpaired reads."
        )


def extract_non_human_reads(
    kraken_out_file, fastq1_in, fastq2_in=None, fastq1_out=None, fastq2_out=None
):
    """Parses Kraken2 output and extracts reads NOT assigned to human (TaxID 9606).

    Handles both single-end/unpaired and paired-end FASTQs natively.
    """
    print(f"[CUSTOM EXTRACTION]: Parsing {kraken_out_file}...")

    # Step 1: Collect Read IDs to EXCLUDE (Human reads)
    human_read_ids = set()
    with open(kraken_out_file, "r") as kfile:
        for line in kfile:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                status, read_id, taxid = parts[0], parts[1], parts[2]
                # If classified as Human (9606), mark for removal
                if status == "C" and taxid == "9606":
                    human_read_ids.add(read_id)

    print(
        f"[CUSTOM EXTRACTION]: Found {len(human_read_ids)} human reads to remove."
    )

    # Step 2: Stream through FASTQs and filter out human reads
    def filter_single_fastq(in_path, out_path):
        count_kept = 0
        with (
            open_fq(in_path, "rt") as infile,
            open(out_path, "w") as outfile,
        ):
            while True:
                header = infile.readline()
                if not header:
                    break
                seq = infile.readline()
                strand = infile.readline()
                qual = infile.readline()

                # Extract read ID (everything after @ up to first space or slash)
                read_id = header[1:].split()[0].rstrip("/1").rstrip("/2")

                if read_id not in human_read_ids:
                    outfile.write(f"{header}{seq}{strand}{qual}")
                    count_kept += 1
        return count_kept

    if fastq2_in and fastq2_out:
        # Paired-end filtering
        count_kept = 0
        with (
            open_fq(fastq1_in, "rt") as in1,
            open_fq(fastq2_in, "rt") as in2,
            open(fastq1_out, "w") as out1,
            open(fastq2_out, "w") as out2,
        ):

            while True:
                h1, s1, st1, q1 = (
                    in1.readline(),
                    in1.readline(),
                    in1.readline(),
                    in1.readline(),
                )
                h2, s2, st2, q2 = (
                    in2.readline(),
                    in2.readline(),
                    in2.readline(),
                    in2.readline(),
                )
                if not h1 or not h2:
                    break

                read_id = h1[1:].split()[0].rstrip("/1").rstrip("/2")

                if read_id not in human_read_ids:
                    out1.write(f"{h1}{s1}{st1}{q1}")
                    out2.write(f"{h2}{s2}{st2}{q2}")
                    count_kept += 1
        print(f"[CUSTOM EXTRACTION]: Kept {count_kept} non-human read pairs.")
    else:
        # Single-end/Unpaired filtering
        count_kept = filter_single_fastq(fastq1_in, fastq1_out)
        print(
            f"[CUSTOM EXTRACTION]: Kept {count_kept} non-human unpaired reads."
        )


def run_kraken_classify(input_reads, output_file, report_file, paired=False):
    """Run Kraken2 classification using the k2 classify subcommand."""
    cmd = [
        k2_path,
        "classify",
        "--db",
        kraken2_db_path,
        "--threads",
        str(threads),
    ]

    if use_k2_daemon:
        cmd.append("--use-daemon")

    if paired:
        cmd.extend(["--paired", str(input_reads[0]), str(input_reads[1])])
    else:
        cmd.append(str(input_reads[0]))

    cmd.extend(["--output", str(output_file), "--report", str(report_file)])
    print(f"[K2 CLASSIFY]: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def run_fastp_trimming(
    in1, in2=None, out1=None, out2=None, html_report=None, json_report=None
):
    """Trims FASTQ inputs using fastp with timeout/retry protection."""
    if not os.path.exists(fastp_path):
        raise FileNotFoundError(f"fastp executable not found at {fastp_path}")

    fastp_cmd = [
        fastp_path,
        "-i",
        str(in1),
        "-o",
        str(out1),
        "-w",
        "8",
    ]

    if in2 and out2:
        fastp_cmd.extend(["-I", str(in2), "-O", str(out2)])

    if html_report:
        fastp_cmd.extend(["-h", str(html_report)])
    if json_report:
        fastp_cmd.extend(["-j", str(json_report)])

    for attempt in range(1, fastp_retries + 1):
        print(
            f"[FASTP TRIMMING]: Attempt {attempt}/{fastp_retries} for {in1}"
        )

        for p in [out1, out2]:
            if p and os.path.exists(p):
                os.remove(p)

        try:
            result = subprocess.run(
                fastp_cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=fastp_timeout,
            )
            print(f"[FASTP TRIMMING]: Completed successfully on attempt {attempt}.")
            if out1 and not os.path.exists(out1):
                raise FileNotFoundError(f"fastp did not create output file: {out1}")
            if out2 and not os.path.exists(out2):
                raise FileNotFoundError(f"fastp did not create output file: {out2}")
            return
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
            stderr = ""
            if hasattr(e, "stderr") and e.stderr:
                stderr = e.stderr.decode(errors="replace") if isinstance(e.stderr, bytes) else str(e.stderr)
            elif hasattr(e, "stdout") and e.stdout:
                stderr = e.stdout.decode(errors="replace") if isinstance(e.stdout, bytes) else str(e.stdout)

            print(f"[FASTP TRIMMING]: Attempt {attempt} failed: {e}")
            if stderr:
                print(stderr.strip())

            if attempt < fastp_retries:
                print(f"[FASTP TRIMMING]: Retrying fastp after failure...\n")
                continue

            raise RuntimeError(
                f"fastp failed after {fastp_retries} attempts for {in1}"
            ) from e


def process_accession(acc):
    print(f"\n==========================================")
    print(f" Processing Accession: {acc}")
    print(f"==========================================")

    out_path = Path(output_dir) / acc
    out_path.mkdir(parents=True, exist_ok=True)

    r1 = Path(read_dir) / f"{acc}_1.fastq.gz"
    r2 = Path(read_dir) / f"{acc}_2.fastq.gz"
    u = Path(read_dir) / f"{acc}.fastq.gz"

    # --- 1. PROCESS PAIRED READS ---
    if r1.exists() and r2.exists():
        paired_kraken_out = out_path / f"{acc}_paired.kraken"
        paired_report_out = out_path / f"{acc}_paired.report"

        run_kraken_classify(
            input_reads=[r1, r2],
            output_file=paired_kraken_out,
            report_file=paired_report_out,
            paired=True,
        )

        # Extract Non-Human Paired Reads
        clean_r1 = out_path / f"{acc}_scrubbed_1.fastq"
        clean_r2 = out_path / f"{acc}_scrubbed_2.fastq"
        extract_non_human_reads(
            paired_kraken_out,
            r1,
            fastq2_in=r2,
            fastq1_out=clean_r1,
            fastq2_out=clean_r2,
        )

    # --- 2. PROCESS UNPAIRED READS ---
    if u.exists():
        unpaired_kraken_out = out_path / f"{acc}_unpaired.kraken"
        unpaired_report_out = out_path / f"{acc}_unpaired.report"

        run_kraken_classify(
            input_reads=[u],
            output_file=unpaired_kraken_out,
            report_file=unpaired_report_out,
            paired=False,
        )

        # Extract Non-Human Unpaired Reads
        clean_u = out_path / f"{acc}_scrubbed_unpaired.fastq"
        extract_non_human_reads(unpaired_kraken_out, u, fastq1_out=clean_u)

    # --- 3. DETERMINE DOMINANT VIRAL SPECIES & EXTRACT ITS READS ---
    report_paths = []
    if r1.exists() and r2.exists():
        report_paths.append(paired_report_out)
    if u.exists():
        report_paths.append(unpaired_report_out)

    dominant_species = find_dominant_viral_species(
        report_paths, viral_scope_taxid
    )

    if dominant_species:
        print(
            f"\n[Dominant viral species]: {dominant_species['name']} "
            f"(TaxID {dominant_species['taxid']}, "
            f"{dominant_species['reads']} lineage reads)"
        )

        # Record the dominant species for this accession
        record_file = out_path / f"{acc}_dominant_viral_species.txt"
        with open(record_file, "w") as rf:
            rf.write(f"Accession: {acc}\n")
            rf.write(f"Dominant Viral Species: {dominant_species['name']}\n")
            rf.write(f"TaxID: {dominant_species['taxid']}\n")
            rf.write(f"Lineage Reads Covered: {dominant_species['reads']}\n")
        print(f"[INFO]: Dominant species recorded in {record_file}")

        # Use one report to derive the allowed taxid lineage
        allowed_taxids = get_allowed_taxids(
            report_paths[0], viral_scope_taxid, dominant_species["taxid"]
        )

        # Extract dominant-species subset from the already scrubbed FASTQs
        if r1.exists() and r2.exists():
            dom_r1 = out_path / f"{acc}_dominant_species_1.fastq"
            dom_r2 = out_path / f"{acc}_dominant_species_2.fastq"
            extract_dominant_species_reads(
                paired_kraken_out,
                clean_r1,
                fastq2_in=clean_r2,
                fastq1_out=dom_r1,
                fastq2_out=dom_r2,
                allowed_taxids=allowed_taxids,
            )

            # --- 4A. FASTP TRIMMING FOR PAIRED READS ---
            trimmed_r1 = out_path / f"{acc}_dominant_species_1_trimmed.fastq"
            trimmed_r2 = out_path / f"{acc}_dominant_species_2_trimmed.fastq"
            fastp_html = out_path / f"{acc}_fastp_paired.html"
            fastp_json = out_path / f"{acc}_fastp_paired.json"

            run_fastp_trimming(
                in1=dom_r1,
                in2=dom_r2,
                out1=trimmed_r1,
                out2=trimmed_r2,
                html_report=fastp_html,
                json_report=fastp_json,
            )

        if u.exists():
            dom_u = out_path / f"{acc}_dominant_species_unpaired.fastq"
            extract_dominant_species_reads(
                unpaired_kraken_out,
                clean_u,
                fastq1_out=dom_u,
                allowed_taxids=allowed_taxids,
            )

            # --- 4B. FASTP TRIMMING FOR UNPAIRED READS ---
            trimmed_u = (
                out_path / f"{acc}_dominant_species_unpaired_trimmed.fastq"
            )
            fastp_html_u = out_path / f"{acc}_fastp_unpaired.html"
            fastp_json_u = out_path / f"{acc}_fastp_unpaired.json"

            run_fastp_trimming(
                in1=dom_u,
                out1=trimmed_u,
                html_report=fastp_html_u,
                json_report=fastp_json_u,
            )
    else:
        print(
            f"[WARNING]: No viral species found within scope "
            f"{viral_scope_taxid} for {acc}; skipping dominant-species extraction."
        )


def main():
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for acc in accession_list:
        try:
            process_accession(acc)
        except Exception as e:
            print(f"[ERROR]: Processing failed for {acc}: {e}")


if __name__ == "__main__":
    main()
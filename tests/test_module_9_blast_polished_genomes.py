import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "pipeline" / "operational" / "module_9_blast_polished_genomes.py"


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Expected module at {MODULE_PATH} before implementation")
    spec = importlib.util.spec_from_file_location("module_9_blast_polished_genomes", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_find_polished_genome_prefers_combined_consensus(tmp_path):
    module = load_module()
    out_dir = tmp_path / "output"
    sample_dir = out_dir / "module_6" / "sampleA"
    sample_dir.mkdir(parents=True)
    combined = sample_dir / "sampleA_consensus_combined.fasta"
    combined.write_text(">seq\nACGT\n")

    class Config:
        cross_module_params = {"output_dir": str(out_dir)}

    assert module.find_polished_genome(Config(), "sampleA") == str(combined)


def test_get_blast_db_path_uses_module_9_value():
    module = load_module()

    class Config:
        cross_module_params = {"output_dir": "/tmp/output"}
        module_9_params = {"blast_db_path": "/custom/blast_db"}

    assert module.get_blast_db_path(Config()) == "/custom/blast_db"


def test_run_blastn_for_sample_uses_configured_threads(monkeypatch, tmp_path):
    module = load_module()
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        if cmd[:2] == ["which", "blastn"]:
            return type("Result", (), {"returncode": 0, "stdout": "/usr/bin/blastn\n", "stderr": ""})()
        return type("Result", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    db_dir = tmp_path / "blast_db"
    db_dir.mkdir()
    (db_dir / "blast_db.nin").write_text("stub")

    class Config:
        cross_module_params = {"output_dir": str(tmp_path / "output")}
        module_9_params = {
            "blastn_bin_path": "blastn",
            "blast_db_path": str(db_dir / "blast_db"),
            "evalue": 1e-6,
            "max_target_seqs": 20,
            "outfmt": "6 qseqid sseqid pident length qlen slen bitscore evalue staxids sscinames stitle",
            "threads": 12,
        }

    module.run_blastn_for_sample(Config(), "sampleA", str(tmp_path / "query.fa"), str(tmp_path / "out.tsv"))
    blast_cmd = calls[-1]
    assert "-num_threads" in blast_cmd
    assert blast_cmd[blast_cmd.index("-num_threads") + 1] == "12"


def test_extract_blast_taxonomy_uses_staxids_and_sscinames():
    module = load_module()

    result = module.extract_blast_taxonomy("1000;2000", "Human orthopoxvirus;Other virus")
    assert result["subject_taxid"] == 1000
    assert result["species_taxid"] == 1000
    assert result["species_sci_name"] == "Human orthopoxvirus"

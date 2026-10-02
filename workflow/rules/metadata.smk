rule record_effective_config:
    input:
        project_path(config["project"]["config_file"])
    output:
        output_path("effective_config")
    params:
        serialized=yaml.safe_dump(config, sort_keys=False)
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text(params.serialized)


rule record_input_manifest:
    input:
        raw=input_path("fasta"),
        card=[str(RGI_DATABASE / "card.json")] if RGI_ENABLED else [],
        swiss=[database_path("swissprot_diamond")] if SWISS_ENABLED else [],
        code=[str(p) for p in sorted(WORKFLOW_DIR.rglob("*")) if p.is_file() and p.suffix in {".py", ".smk"}] + [str(WORKFLOW_DIR / "Snakefile"), str(REPO_ROOT / "run.py"), str(REPO_ROOT / "scripts/setup_databases.py"), str(REPO_ROOT / "environment.yml")]
    output:
        output_path("input_manifest")
    run:
        import hashlib
        import csv
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["role", "path", "bytes", "sha256"])
            records = [("input_fasta", input.raw)]
            records += [("card_json", p) for p in input.card]
            records += [("swissprot_diamond", p) for p in input.swiss]
            records += [("pipeline_code", p) for p in input.code]
            for role, name in records:
                path = Path(name)
                digest = hashlib.sha256()
                with path.open("rb") as source:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        digest.update(block)
                writer.writerow([role, str(path.resolve()), path.stat().st_size, digest.hexdigest()])


rule record_software_versions:
    output:
        output_path("software_versions")
    params:
        executables={key: tool(key) for key in ["python", "snakemake", "mmseqs", "prodigal"] + (["diamond"] if SWISS_ENABLED or RGI_ENABLED else []) + (["rgi"] if RGI_ENABLED else [])},
        configuration=yaml.safe_dump(config, sort_keys=True)
    run:
        import subprocess
        import csv
        import json
        import datetime
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["component", "executable", "version_or_detail"])
            writer.writerow(["arg-annotator", str(REPO_ROOT), (REPO_ROOT / "VERSION").read_text().strip()])
            writer.writerow(["recorded_utc", "", datetime.datetime.now(datetime.timezone.utc).isoformat()])
            for name, executable in params.executables.items():
                flag = "version" if name == "mmseqs" else "-v" if name == "prodigal" else "--version"
                command = [executable, "main", "--version"] if name == "rgi" else [executable, flag]
                result = subprocess.run(command, capture_output=True, text=True, timeout=30)
                if result.returncode:
                    raise RuntimeError(f"Version query failed: {command}: {result.stderr}")
                detail = " | ".join((result.stdout + result.stderr).strip().splitlines())
                writer.writerow([name, executable, detail])
            if RGI_ENABLED:
                card = json.loads((RGI_DATABASE / "card.json").read_text())
                writer.writerow(["CARD", str(RGI_DATABASE / "card.json"), card.get("_version", "unknown; see input SHA256")])
            writer.writerow(["SwissProt", database_path("swissprot_diamond") if SWISS_ENABLED else "disabled", "See database_manifest.json and input_manifest.tsv for release provenance and SHA256"])

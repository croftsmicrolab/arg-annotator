#!/usr/bin/env python3
"""Prepare external CARD and Swiss-Prot databases; never bundle them in Git."""
from __future__ import annotations

import argparse
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parent.parent
CARD_URL = "https://card.mcmaster.ca/latest/data"
SWISS_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "arg-annotator/0.1.0"})
    part = destination.with_suffix(destination.suffix + ".part")
    print(f"Downloading {url}", flush=True)
    try:
        with urllib.request.urlopen(request, timeout=180) as source, part.open("wb") as target:
            metadata = {"url": url, "resolved_url": source.url, "last_modified": source.headers.get("Last-Modified", "")}
            shutil.copyfileobj(source, target)
        part.replace(destination)
    finally:
        part.unlink(missing_ok=True)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Download current CARD and Swiss-Prot releases")
    parser.add_argument("--card-json", type=Path, help="Already downloaded/extracted CARD card.json")
    parser.add_argument("--swissprot-fasta", type=Path, help="Reviewed Swiss-Prot FASTA or FASTA.gz")
    parser.add_argument("--swissprot-release", default="", help="Release label from your downloaded release notes")
    parser.add_argument("--database-dir", type=Path, default=ROOT / "databases")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Explicitly replace already prepared databases in this directory")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not args.download and not (args.card_json and args.swissprot_fasta):
        parser.error("Use --download, or supply both --card-json and --swissprot-fasta")
    for key, path in [("CARD", args.card_json), ("Swiss-Prot", args.swissprot_fasta)]:
        if path and not path.is_file():
            parser.error(f"{key} source does not exist: {path}")
    for name in ["rgi", "diamond"]:
        if not shutil.which(name):
            parser.error(f"{name} not found; activate the Conda environment first")
    db = args.database_dir.resolve()
    card_dir, swiss_dir, workdir = db / "card", db / "swissprot", db / "rgi"
    for directory in (card_dir, swiss_dir, workdir):
        directory.mkdir(parents=True, exist_ok=True)
    card_target = workdir / "localDB/card.json"
    swiss_target = swiss_dir / "uniprot_sprot.dmnd"
    manifest_path = db / "database_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    # Fail rather than overwrite a database whose provenance is unknown.
    for key, target in [("CARD", card_target), ("SwissProt", swiss_target)]:
        if target.exists() and not args.force and key not in manifest:
            parser.error(f"{target} already exists without a setup manifest. Use another directory or --force explicitly.")
    with tempfile.TemporaryDirectory(prefix="setup_", dir=db) as temporary:
        temporary = Path(temporary)
        if not card_target.exists() or args.force:
            if args.card_json:
                card_source = args.card_json.resolve()
                card_meta = {"source": str(card_source)}
            else:
                archive = temporary / "card.tar.bz2"
                card_meta = download(CARD_URL, archive)
                card_meta["download_sha256"] = sha256(archive)
                with tarfile.open(archive, "r:*") as package:
                    matches = [member for member in package.getmembers() if member.isfile() and Path(member.name).name == "card.json"]
                    if len(matches) != 1:
                        raise RuntimeError("Expected one card.json in the CARD archive")
                    card_source = temporary / "card.json"
                    with package.extractfile(matches[0]) as source, card_source.open("wb") as target:
                        shutil.copyfileobj(source, target)
            card = json.loads(card_source.read_text())
            if "_version" not in card:
                raise ValueError("CARD JSON lacks _version; provide the curated CARD card.json")
            if card_source == card_target.resolve():
                raise ValueError("Provide an external CARD source, not the destination database itself")
            # Rebuilding avoids retaining cached reference indexes from an older CARD release.
            if (workdir / "localDB").exists() and not args.force:
                raise ValueError("Existing incomplete localDB; choose another database directory or explicitly use --force")
            if (workdir / "localDB").exists():
                shutil.rmtree(workdir / "localDB")
            subprocess.run(["rgi", "load", "--card_json", str(card_source), "--local"], cwd=workdir, check=True)
            if not card_target.is_file():
                raise RuntimeError("RGI did not create localDB/card.json")
            manifest["CARD"] = {**card_meta, "version": card["_version"], "sha256": sha256(card_target), "prepared_path": str(card_target)}
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        else:
            print("Reusing prepared CARD database")
        if not swiss_target.exists() or args.force:
            if args.swissprot_fasta:
                swiss_source = args.swissprot_fasta.resolve()
                swiss_meta = {"source": str(swiss_source)}
            else:
                swiss_source = temporary / "uniprot_sprot.fasta.gz"
                swiss_meta = download(SWISS_URL, swiss_source)
            swiss_meta["source_sha256"] = sha256(swiss_source)
            if swiss_source.suffix == ".gz":
                plain = temporary / "swissprot.fasta"
                with gzip.open(swiss_source, "rb") as source, plain.open("wb") as target:
                    shutil.copyfileobj(source, target)
            else:
                plain = swiss_source
            # Build in temporary space and replace the finished index atomically.
            prefix = temporary / "uniprot_sprot"
            subprocess.run(["diamond", "makedb", "--in", str(plain), "--db", str(prefix), "--threads", str(args.threads)], check=True)
            prefix.with_suffix(".dmnd").replace(swiss_target)
            manifest["SwissProt"] = {**swiss_meta, "release": args.swissprot_release or "unspecified; identify from source checksum/date", "sha256": sha256(swiss_target), "prepared_path": str(swiss_target)}
        else:
            print("Reusing prepared Swiss-Prot database")
    manifest["prepared_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    manifest["commands"] = {name: subprocess.run([name, "main", "--version"] if name == "rgi" else [name, "--version"], capture_output=True, text=True, check=True).stdout.strip() for name in ["rgi", "diamond"]}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Databases ready. Provenance: {manifest_path}")
    if db != ROOT / "databases":
        print("For this custom directory, update the three databases paths in config/config.yaml.")


if __name__ == "__main__":
    main()

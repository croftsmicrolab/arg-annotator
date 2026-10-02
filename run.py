#!/usr/bin/env python3
"""Run arg-annotator from an activated Conda environment."""
import argparse
import copy
import os
from pathlib import Path
import shutil
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parent

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=ROOT / "config/config.yaml")
    p.add_argument("--input", type=Path, help="Processed nucleotide insert FASTA (not FASTQ).")
    p.add_argument("--sample", help="Sample label for documentation and plots.")
    p.add_argument("--output", type=Path, help="Independent results directory; relative to repository.")
    p.add_argument("--cores", type=int, help="Total CPU cap; defaults to resources.default_threads (8).")
    p.add_argument("--top-percent", type=float, help="Top percent, e.g. 5 or 10. Activates top_fraction mode.")
    p.add_argument("--min-support", type=int, help="Keep every RepORF with support >= this number (at least 2).")
    p.add_argument("--skip-annotation", action="store_true", help="Run clustering, selection, and plots without databases.")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--unlock", action="store_true", help="Unlock an interrupted run after confirming it has stopped.")
    p.add_argument("--keep-going", action="store_true")
    a = p.parse_args()
    if a.cores is not None and a.cores < 1:
        p.error("--cores must be positive")
    if a.top_percent is not None and a.min_support is not None:
        p.error("Choose either --top-percent or --min-support")
    if a.top_percent is not None and not 0 < a.top_percent <= 100:
        p.error("--top-percent must be in (0, 100]")
    if a.min_support is not None and a.min_support < 2:
        p.error("--min-support must be at least 2")
    config_path = a.config.expanduser().resolve()
    c = copy.deepcopy(yaml.safe_load(config_path.read_text()))
    cores = a.cores if a.cores is not None else int(c["resources"].get("default_threads", 8))
    if cores < 1:
        p.error("resources.default_threads must be positive")
    c["project"]["root"] = str(ROOT)
    c["project"]["config_file"] = str(config_path)
    if a.input:
        c["inputs"]["fasta"] = str(a.input.expanduser().resolve())
    if a.sample:
        c["project"]["sample"] = a.sample
    if a.top_percent is not None:
        c["selection"].update(mode="top_fraction", nominal_fraction=a.top_percent / 100)
    if a.min_support is not None:
        c["selection"].update(mode="minimum_support", minimum_support=a.min_support)
    if a.skip_annotation:
        c["annotation"]["rgi"]["enabled"] = False
        c["annotation"]["swissprot"]["enabled"] = False
    out = a.output or Path(c["outputs"]["root"])
    if not out.is_absolute():
        out = ROOT / out
    out = out.resolve()
    raw = Path(c["inputs"]["fasta"])
    if not raw.is_absolute():
        raw = ROOT / raw
    if not raw.is_file():
        p.error(f"Input FASTA does not exist: {raw}")
    if out == ROOT or raw.is_relative_to(out):
        p.error("Results must be separate from the repository root and input FASTA")
    snakemake = shutil.which(c.get("tools", {}).get("snakemake", "snakemake"))
    if not snakemake:
        p.error("snakemake not found; activate the environment from environment.yml")
    out.mkdir(parents=True, exist_ok=True)
    c.setdefault("run", {})["output_root"] = str(out)
    c["run"].update(cores=cores, scheduler="greedy")
    # A separate working directory also isolates Snakemake state for concurrent samples.
    work = out / ".workflow"
    work.mkdir(exist_ok=True)
    runtime = work / "runtime_config.yaml"
    text = yaml.safe_dump(c, sort_keys=False)
    if not runtime.exists() or runtime.read_text() != text:
        runtime.write_text(text)
    cmd = [snakemake, "--snakefile", str(ROOT / "workflow/Snakefile"), "--configfile", str(runtime),
           "--cores", str(cores), "--scheduler", "greedy", "--rerun-incomplete", "--printshellcmds"]
    for enabled, flag in [(a.dry_run, "--dry-run"), (a.unlock, "--unlock"), (a.keep_going, "--keep-going")]:
        if enabled:
            cmd.append(flag)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLCONFIGDIR"] = str(out / "logs/matplotlib")
    env["XDG_CACHE_HOME"] = str(out / "logs/cache")
    env["TMPDIR"] = str(out / ".tmp")
    Path(env["TMPDIR"]).mkdir(exist_ok=True)
    print(f"arg-annotator | sample={c['project']['sample']} | results={out}", flush=True)
    return subprocess.call(cmd, cwd=work, env=env)

if __name__ == "__main__":
    sys.exit(main())

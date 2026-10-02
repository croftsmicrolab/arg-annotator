#!/usr/bin/env python3
"""Isolate RGI's shell-sensitive paths and verify its actual outputs."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--query", required=True, type=Path)
    p.add_argument("--database", required=True, type=Path)
    p.add_argument("--txt", required=True, type=Path)
    p.add_argument("--json", required=True, type=Path)
    p.add_argument("--rgi", default="rgi")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--enabled", type=int, choices=[0, 1], default=1)
    p.add_argument("--include-loose", type=int, choices=[0, 1], default=1)
    p.add_argument("--low-quality", type=int, choices=[0, 1], default=1)
    p.add_argument("--clean", type=int, choices=[0, 1], default=1)
    a = p.parse_args()
    for target in [a.txt, a.json]:
        target.parent.mkdir(parents=True, exist_ok=True)
    if not a.enabled or a.query.stat().st_size == 0:
        a.txt.write_text("")
        a.json.write_text("{}\n")
        print("RGI skipped: disabled or no candidates")
        return
    # RGI 6.0.8 builds internal shell command strings without quoting paths.
    # Its localDB path is derived from cwd, so a temporary symlink provides a
    # path without spaces without copying or changing the source CARD JSON.
    base = Path(tempfile.gettempdir())
    if any(char.isspace() for char in str(base)):
        base = Path("/tmp").resolve()
    with tempfile.TemporaryDirectory(prefix="arg_rgi_", dir=base) as directory:
        work = Path(directory)
        (work / "localDB").symlink_to(a.database.resolve(), target_is_directory=True)
        query = work / "query.faa"
        shutil.copyfile(a.query, query)
        prefix = work / "result"
        command = [a.rgi, "main", "-i", str(query), "-o", str(prefix),
                   "-t", "protein", "-a", "DIAMOND", "-n", str(a.threads), "--local"]
        for enabled, flag in [(a.include_loose, "--include_loose"), (a.low_quality, "--low_quality"), (a.clean, "--clean")]:
            if enabled:
                command.append(flag)
        env = os.environ.copy()
        env["TMPDIR"] = str(work)
        subprocess.run(command, cwd=work, env=env, check=True)
        for suffix, target in [(".txt", a.txt), (".json", a.json)]:
            source = prefix.with_suffix(suffix)
            if not source.is_file():
                raise RuntimeError(f"RGI exited without {suffix} output; inspect the RGI logs")
            shutil.copyfile(source, target)


if __name__ == "__main__":
    main()

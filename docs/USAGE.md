# Usage

## Input and upstream processing

Input is plain, uncompressed nucleotide FASTA of already processed inserts or
CCS sequences. This workflow does not perform basecalling, FASTQ QC, assembly,
vector trimming or library demultiplexing. Preserve one record per read when
read-count support is intended. Already deduplicated FASTA loses the original
weights; abundance encoded in a FASTA header is not interpreted.

By default only A/C/G/T/N are accepted, empty sequences fail, sequences become
uppercase, and unsafe/duplicate IDs become unique IDs. Original headers and
ID mappings are preserved in `00_input/`. Run separately for different selections
or libraries unless pooled read-count support is deliberately intended.

## Configuration and overrides

`run.py` loads `config/config.yaml`. Use `--config path/to/custom.yaml` for a
complete custom configuration. CLI input, sample, output and selection overrides
take precedence. Input paths on the CLI are resolved from the invoking directory;
relative paths in the YAML and `--output` are resolved from the repository root.

```bash
python run.py --input /data/inserts.fasta --sample my_library --output results/my_library --cores 8
python run.py --help
```

Use `--top-percent 5` to select the percentage mode, or `--min-support 2` for the
support mode. These options are mutually exclusive. YAML equivalent:

```yaml
selection:
  mode: top_fraction   # alternatively minimum_support
  nominal_fraction: 0.05
  minimum_support: 2   # direct cutoff, or support-one fallback cutoff
```

Other selection metadata in the default configuration documents fixed implementation
behavior; rounding is ceiling, tie order is support/member-count/ID, and support-one
records are always exported. Editing those documentation fields does not change
the algorithm. Thresholds must be integers ≥2 and fractions must be in (0,1].

Clustering coverage and identity values are MMseqs2 parameters, not percent values.
For protein annotation use RGI `input_type: protein`, aligner `DIAMOND`, and
`local: true`; the prepared RGI database is under `databases/rgi/localDB`.
Swiss-Prot sensitivity can be `fast`, `sensitive`, `more-sensitive`,
`very-sensitive` or `ultra-sensitive`. Its input is predicted protein, so the
command is `blastp`. `blastx` would apply to nucleotide queries and is not used here.

## Execution, resuming and diagnostics

Activate the environment before running. `--cores` caps total concurrent CPU
threads; the resource settings in the YAML cap individual tools. Always try
`--dry-run` before a larger run. Snakemake metadata and runtime settings are
stored inside your output's `.workflow/`, so no existing project metadata is used.
Rerun the same command to resume incomplete work. Use separate output directories
when comparing configurations. `--unlock` is for a stale Snakemake lock after
confirming no active process is using that output directory.

A normal run requires prepared databases even for a dry-run; `--skip-annotation`
allows database-free clustering and plots. Empty primary selections are valid:
candidate files are empty/header-only, annotations are skipped and validation
still runs. An input with no predictable ORFs cannot proceed to ORF clustering;
inspect Prodigal's logs and input length/content. Invalid/missing parent mappings
cause workflow failure rather than a valid final result.

MMseqs intermediate databases and temporary clustering files live inside the
chosen output directory and can be large. Preserve scientific outputs before
removing an entire output directory. The workflow does not automatically delete
your results or modify the input/reference source files. RGI can create cached
reference indexes in its configured localDB; keep that directory writable and
use a separate database working directory when changing CARD versions.

## Plots and table interpretation

Both curves sort positive counts descending. They use a linear rank x-axis,
logarithmic count y-axis, blue line and dotted grids. Plotting uses Matplotlib,
so the runtime does not need R.
The weighted curve reads `reporf_support`, not `orf_cluster_member_count`.
Ties in the plotted rank TSV use ID ordering; the candidate-selection table
also uses ORF member count as its second key. Tied values have identical heights.

Each candidate has one combined annotation row. Raw hit files preserve multiple
hits; the combined table reports best evidence and hit counts. Low confidence,
partial ORF and no-hit rows are retained. Inspect RGI model scores/coverage and
Swiss-Prot coverage alongside the entire insert before assigning function.

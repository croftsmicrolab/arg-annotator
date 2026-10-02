# Computational methods

Processed nucleotide insert sequences are normalized and clustered using
MMseqs2 with minimum sequence identity 0.99, coverage 0.99, coverage mode 0
and cluster mode 0. ORFs are predicted from representative inserts using
Prodigal in metagenomic mode. Predicted **nucleotide** ORFs are clustered with
minimum sequence identity 0.99, coverage 0.80, coverage mode 1 and cluster mode 2.
MMseqs2's coverage mode and cluster mode are explicitly retained because
changing them can alter membership and representative identities.

For each representative ORF cluster R, abundance-weighted support is:

```text
support(R) = sum(size(C) for C in unique_parent_insert_clusters(R))
```

Here `size(C)` counts original input FASTA records in insert cluster C. A parent
cluster contributes once to a given ORF cluster even when multiple ORFs from that
parent occur in it. Clusters are ranked by support descending, ORF member count
descending and representative ID ascending.

In top-fraction mode the nominal number selected is `ceil(N × f)`, with default
`f=0.05`. If the nominal boundary has support 1, selection switches to all clusters
with support at least the configured cutoff (default 2). In alternative minimum-
support mode, every cluster meeting that cutoff is selected. Support-one clusters
are exported separately. Top-fraction selection does not expand all boundary ties
when boundary support exceeds 1; therefore software/representative changes can
affect exactly which tied records are retained.

Selected representative proteins are annotated with RGI against a locally loaded
CARD reference, using protein input, the DIAMOND aligner, inclusion of Loose hits,
low-quality handling and clean output. Trailing translation stop markers are removed
from the RGI query FASTA. DIAMOND blastp searches the reviewed Swiss-Prot reference
in sensitive mode with E-value threshold 1e-3 and at most 25 target sequences.
These settings are configurable and their effective values are saved for each run.

Annotation results are merged with support and ORF coordinates/completeness.
Candidate rows are retained irrespective of annotation confidence, including
partial ORFs and no-hit candidates. RGI hit categories and Swiss-Prot identity,
query/subject coverage, E-values and scores support evidence assessment. The
software's confidence classes are descriptive heuristics rather than experimental
resistance validation; weak/fragmentary acetyltransferase similarities should not
be interpreted as validated ARGs.

Insert-cluster sizes and weighted RepORF supports are plotted as descending rank
curves with a linear rank axis and logarithmic count axis. Structural validation
checks complete memberships, matching nucleotide/protein IDs, unique parent
mapping, support recalculation, selection and annotation retention. Effective
configuration, software/reference versions and SHA256 manifests accompany outputs.

Use this paragraph as a starting point for manuscript Methods, replacing defaults
with the effective parameters and actual database releases of the analysis reported.

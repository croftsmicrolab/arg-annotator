Place processed nucleotide insert/CCS sequences in `sequences.fasta`, or use
`python run.py --input /path/to/sequences.fasta`. Input must be plain FASTA,
not FASTQ or paired short reads. Each FASTA record contributes one count.
Perform read QC, adapter/vector removal and insert extraction upstream.
User sequences are ignored by Git.

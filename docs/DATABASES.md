# External database preparation

Activate the environment, then choose one method:

```bash
python scripts/setup_databases.py --download
# Or use exact locally downloaded releases:
python scripts/setup_databases.py --card-json /data/CARD/card.json \
  --swissprot-fasta /data/SwissProt/uniprot_sprot.fasta.gz \
  --swissprot-release YOUR_RELEASE
```

`--download` streams the current CARD archive and reviewed Swiss-Prot FASTA
from their official servers. Only `card.json` is read from the CARD archive;
RGI loads it with `rgi load --card_json ... --local`, run inside `databases/rgi/`.
DIAMOND builds `databases/swissprot/uniprot_sprot.dmnd` from Swiss-Prot protein
FASTA. Download archives/decompressed temporary files are removed afterward.
References and indexes remain available for subsequent runs and are ignored by Git.

Prepared databases are reused by default. `--force` explicitly rebuilds them;
RGI cached indexes are cleared to prevent reuse of a different CARD version.
Never point this option at a protected scientific database directory. To keep
multiple releases, use `--database-dir /data/references/release_label` and update
all three `databases` paths in your YAML. RGI's localDB must sit directly under
its configured working directory, as required by RGI local database resolution.

The setup manifest `databases/database_manifest.json` records source path/URL,
download metadata, CARD version, optional Swiss-Prot release label, checksum and
preparation time. Every run also hashes the actual CARD JSON and Swiss-Prot DMND.
For latest Swiss-Prot, obtain the release label from the official release notes
and pass `--swissprot-release`; a checksum/date without that label is explicitly
reported as unspecified rather than inventing a version.

Software version pins do not pin a reference release. Reproducing manuscript
results requires retaining the exact input and references used originally.
CARD and UniProt have their own data-use/citation terms; cite those resources
alongside this software.

Official sources:

- [CARD downloads](https://card.mcmaster.ca/download)
- [CARD current data archive](https://card.mcmaster.ca/latest/data)
- [RGI local database loading documentation](https://github.com/arpcard/rgi/blob/master/docs/rgi_load.rst)
- [UniProt download documentation](https://www.uniprot.org/help/downloads)
- [Reviewed Swiss-Prot FASTA](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)

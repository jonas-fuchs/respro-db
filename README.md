# ResPro Database Conversion Pattern

This repository is intended as a general home for ResPro-compatible, auto-updated
resistance databases. The HSV Daehne/Jaki conversion is the first implemented
example, and more source databases are expected to be added.

## Repository purpose

- fetch curated upstream source data (Zenodo, release assets, or similar)
- transform source data into ResPro-compatible TSV artifacts
- track non-migrated source rows with explicit reasons
- auto-update outputs via GitHub Actions and open update PRs only when a new
  upstream version is detected

## Recommended structure per conversion

Each source database should use its own folder under databases:

```text
databases/<source_name>/
  scripts/
    convert.py
    requirements.txt
  output/
    rules.tsv
    formula-rules.tsv        # only when combination rules are present
    metadata.json
    non-migrated-rules.txt
```

And one workflow per source in .github/workflows:

```text
.github/workflows/<source_name>-autobump.yml
```
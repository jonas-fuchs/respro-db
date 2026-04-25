# ResPro Database Conversion Pattern

This repository is intended as a general home for ResPro-compatible, auto-updated antiviral resistance databases.

## Repository purpose

- fetch curated upstream source data (Zenodo, Github, release assets, or similar)
- transform source data into ResPro-compatible TSV artifacts
- track non-migrated source rows with explicit reasons
- auto-update outputs via GitHub Actions and open update PRs only when a new upstream version is detected

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

## You would like support for your database?

Either open a PR here and add the required files or open an issue which database you would like to have supported. We will do our best to make it possible. Importantly, we will need some kind of API to that database to be able to check for updates in regular intervals and autoupdate the respective files in the repo.

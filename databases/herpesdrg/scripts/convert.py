#!/usr/bin/env python3
"""
Convert HerpesDRG TSV into ResPro-compatible TSV artifacts.

Strategy for handling multiple IC50 values and phenotypes from the same publication:
- When a mutation-antiviral pair appears in multiple sources (by publication),
  the data is aggregated:
  1. IC50 values: Calculate median. Add comment about count only if N>1.
  2. Phenotypes: Detect conflicts (resistant vs sensitive -> contradictory).
  3. Publications: Join with commas.

Outputs:
- rules.tsv (required, contains aggregated atomic rules)
- formula-rules.tsv (optional, only when grouped co-mutation rows are emitted)
- metadata.json (required)
- non-migrated-rules.txt (audit trail)
"""

import argparse
import csv
import hashlib
import json
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

SOURCE_URL = "https://raw.githubusercontent.com/ojcharles/herpesdrg-db/main/herpesdrg-db.tsv"

RULES_COLUMNS = [
    "gene",
    "reference_identifier",
    "position",
    "reference",
    "mutation",
    "antiviral",
    "group_id",
    "member_id",
    "phenotype",
    "fold_ic50",
    "publication",
    "source",
    "comment",
]

FORMULA_COLUMNS = [
    "group_id",
    "antiviral",
    "expression",
    "phenotype",
    "fold_ic50",
    "publication",
    "source",
    "comment",
]

NON_MIGRATED_COLUMNS = [
    "reason",
    "mutation_id",
    "virus",
    "gene",
    "aa_change",
    "co_gene",
    "co_aa",
    "status",
    "note",
    "details",
]

ANTIVIRAL_COLUMNS = [
    "Ganciclovir",
    "Aciclovir",
    "Cidofovir",
    "Foscarnet",
    "Brincidofovir",
    "Letermovir",
    "Brivudine",
    "Penciclovir",
    "Tomeglovir",
    "Maribavir",
    "Cyclopropavir",
    "Amenamevir",
    "Pritelivir",
]

REQUIRED_INPUT_COLUMNS = {
    "mutation_id",
    "virus",
    "gene",
    "aa_change",
    "ref_title",
    "ref_link",
    "ref_doi",
    "co_gene",
    "co_aa",
    "note",
    "status",
    *ANTIVIRAL_COLUMNS,
}

# User-provided fixed mappings.
REFERENCE_BY_VIRUS = {
    "hcmv": "NC_006273",
    "vzv": "NC_001348",
    "adeno": "AC_000008.1",
    "hsv1": "NC_001806",
    "hsv2": "NC_001798",
    "hhv6b": "MF511171",
}


def norm(v: object) -> str:
    if v is None:
        return ""
    return str(v).strip()


def normalize_virus(raw_virus: str) -> str:
    return raw_virus.strip().lower()


def should_exclude(status: str) -> tuple[bool, str]:
    if status.strip().upper() != "A":
        return True, "status_not_active"
    return False, ""


def parse_mutation(aa_change: str) -> tuple[str, str, int]:
    text = aa_change.strip().replace(" ", "")

    insertion = re.fullmatch(r"([A-Za-z])(\d+)(?:insert|ins)([A-Za-z]+)", text, flags=re.IGNORECASE)
    if insertion:
        ref = insertion.group(1).upper()
        pos = int(insertion.group(2))
        inserted = insertion.group(3).upper()
        return ref, f"{ref}{pos}{ref}{inserted}", pos

    sub = re.fullmatch(r"([A-Za-z*])(\d+)([A-Za-z*])", text)
    if sub:
        ref = sub.group(1).upper()
        pos = int(sub.group(2))
        alt = sub.group(3).upper()
        return ref, "*" if alt == "*" else alt, pos

    stop_word = re.fullmatch(r"([A-Za-z])(\d+)(stop)", text, flags=re.IGNORECASE)
    if stop_word:
        return stop_word.group(1).upper(), "*", int(stop_word.group(2))

    fs = re.fullmatch(r"([A-Za-z])(\d+)(frameshift\*?|fs.*)", text, flags=re.IGNORECASE)
    if fs:
        ref = fs.group(1).upper()
        pos = int(fs.group(2))
        return ref, f"{ref}{pos}fsX", pos

    pref_del = re.fullmatch(r"del([A-Za-z])(\d+)", text, flags=re.IGNORECASE)
    if pref_del:
        ref = pref_del.group(1).upper()
        pos = int(pref_del.group(2))
        return ref, f"{ref}{pos}del", pos

    simple_del = re.fullmatch(r"([A-Za-z])(\d+)del", text, flags=re.IGNORECASE)
    if simple_del:
        ref = simple_del.group(1).upper()
        pos = int(simple_del.group(2))
        return ref, f"{ref}{pos}del", pos

    range_del = re.fullmatch(r"([A-Za-z])?(\d+)-(\d+)del", text, flags=re.IGNORECASE)
    if range_del:
        ref = (range_del.group(1) or "").upper()
        start = int(range_del.group(2))
        end = int(range_del.group(3))
        mut = f"del{start}_{end}" if not ref else f"{ref}{start}-{end}del"
        return ref, mut, start

    alt_pref_del = re.fullmatch(r"([A-Za-z])del(\d+)", text, flags=re.IGNORECASE)
    if alt_pref_del:
        ref = alt_pref_del.group(1).upper()
        pos = int(alt_pref_del.group(2))
        return ref, f"{ref}{pos}del", pos

    numeric_only = re.fullmatch(r"\d+", text)
    if numeric_only:
        raise ValueError("numeric_only_mutation_not_supported")

    raise ValueError("unparseable_mutation")


def parse_fold_and_phenotype(value: str) -> tuple[str, str]:
    token = value.strip()
    if not token:
        return "", ""

    lowered = token.lower()
    if lowered == "resistant":
        return "", "resistant"
    if lowered == "polymorphism":
        return "", "sensitive"

    try:
        numeric = float(token)
    except ValueError:
        return "", ""

    if numeric < 0:
        return "", ""
    return f"{numeric:g}", ""


def tsv_from_rows(rows: list[dict], columns: list[str]) -> str:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(norm(row.get(col, "")) for col in columns))
    return "\n".join(lines) + "\n"


def checksum(content: str) -> str:
    return "sha256:" + hashlib.sha256(content.encode("utf-8")).hexdigest()


def non_migrated_text(rows: list[dict]) -> str:
    lines = [
        "# Non-migrated HerpesDRG rows",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
        "\t".join(NON_MIGRATED_COLUMNS),
    ]
    for row in sorted(
        rows,
        key=lambda r: (norm(r.get("reason")), norm(r.get("virus")), norm(r.get("gene")), norm(r.get("aa_change"))),
    ):
        lines.append("\t".join(norm(row.get(col, "")) for col in NON_MIGRATED_COLUMNS))
    return "\n".join(lines) + "\n"


def add_non_migrated(bucket: list[dict], source_row: dict, reason: str, details: str = "") -> None:
    bucket.append(
        {
            "reason": reason,
            "mutation_id": norm(source_row.get("mutation_id")),
            "virus": norm(source_row.get("virus")),
            "gene": norm(source_row.get("gene")),
            "aa_change": norm(source_row.get("aa_change")),
            "co_gene": norm(source_row.get("co_gene")),
            "co_aa": norm(source_row.get("co_aa")),
            "status": norm(source_row.get("status")),
            "note": norm(source_row.get("note")),
            "details": details,
        }
    )


def build_publication(doi: str, link: str) -> str:
    if doi:
        return doi if doi.lower().startswith("doi:") else f"doi:{doi}"
    return ""


def split_combo_mutations(aa_change: str) -> list[str]:
    return [part.strip() for part in aa_change.split(";") if part.strip()]


def make_member_id(gene: str, position: int, mutation: str, idx: int) -> str:
    token = re.sub(r"[^A-Za-z0-9_]+", "_", mutation).strip("_") or "mut"
    return f"{gene}_{position}_{token}_{idx}"


def convert(source_rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    rules_rows = []
    formula_rows = []
    non_migrated = []
    # Dictionary to aggregate single-mutation rules: (gene, ref_id, position, mutation, antiviral) -> aggregation
    aggregated_rules = {}
    group_counter = 0

    def next_group_id() -> str:
        nonlocal group_counter
        group_counter += 1
        return f"G{group_counter:05d}"

    for row in source_rows:
        virus_raw = norm(row.get("virus"))
        gene = norm(row.get("gene"))
        aa_change = norm(row.get("aa_change"))
        note = norm(row.get("note"))
        status = norm(row.get("status"))

        virus_key = normalize_virus(virus_raw)
        if not virus_key:
            add_non_migrated(non_migrated, row, "unknown_virus", "No supported virus mapping")
            continue

        if virus_key not in REFERENCE_BY_VIRUS:
            add_non_migrated(non_migrated, row, "missing_reference_mapping", "Supported virus but no fixed reference")
            continue

        exclude, reason = should_exclude(status)
        if exclude:
            add_non_migrated(non_migrated, row, reason)
            continue

        if not gene:
            add_non_migrated(non_migrated, row, "missing_gene")
            continue
        if not aa_change:
            add_non_migrated(non_migrated, row, "missing_aa_change")
            continue

        mutation_parts = split_combo_mutations(aa_change)
        parsed_mutations = []
        for part in mutation_parts:
            try:
                ref_aa, mutation_token, pos = parse_mutation(part)
            except ValueError as exc:
                add_non_migrated(non_migrated, row, str(exc), details=f"failed_part={part}")
                parsed_mutations = []
                break
            parsed_mutations.append((ref_aa, mutation_token, pos))
        if not parsed_mutations:
            continue

        publication = build_publication(norm(row.get("ref_doi")), norm(row.get("ref_link")))
        source_label = "HerpesDRG"

        has_antiviral_data = False
        emitted_for_row = 0
        for antiviral_column in ANTIVIRAL_COLUMNS:
            raw_value = norm(row.get(antiviral_column))
            if not raw_value:
                continue

            fold_ic50, phenotype = parse_fold_and_phenotype(raw_value)
            if not fold_ic50 and not phenotype:
                continue

            has_antiviral_data = True
            antiviral_norm = antiviral_column.lower()

            if len(parsed_mutations) == 1:
                ref_aa, mutation_token, pos = parsed_mutations[0]
                ref_id = REFERENCE_BY_VIRUS[virus_key]

                # Aggregation key: (gene, ref_id, position, mutation, antiviral) without publication
                aggr_key = (gene, ref_id, pos, mutation_token, antiviral_norm)
                
                if aggr_key not in aggregated_rules:
                    aggregated_rules[aggr_key] = {
                        "gene": gene,
                        "reference_identifier": ref_id,
                        "position": pos,
                        "reference": ref_aa,
                        "mutation": mutation_token,
                        "antiviral": antiviral_norm,
                        "group_id": "",
                        "member_id": "",
                        "source": source_label,
                        "publications": [],
                        "phenotypes": [],
                        "ic50_values": [],
                    }
                
                # Aggregate data from this row
                if publication:
                    aggregated_rules[aggr_key]["publications"].append(publication)
                if phenotype:
                    aggregated_rules[aggr_key]["phenotypes"].append(phenotype)
                if fold_ic50:
                    try:
                        aggregated_rules[aggr_key]["ic50_values"].append(float(fold_ic50))
                    except ValueError:
                        pass
                emitted_for_row += 1
                continue

            group_id = next_group_id()
            member_ids = []
            for idx, (ref_aa, mutation_token, pos) in enumerate(parsed_mutations, start=1):
                member_id = make_member_id(gene, pos, mutation_token, idx)
                member_ids.append(member_id)
                rules_rows.append(
                    {
                        "gene": gene,
                        "reference_identifier": REFERENCE_BY_VIRUS[virus_key],
                        "position": pos,
                        "reference": ref_aa,
                        "mutation": mutation_token,
                        "antiviral": antiviral_norm,
                        "group_id": group_id,
                        "member_id": member_id,
                        "phenotype": "",
                        "fold_ic50": "",
                        "publication": publication,
                        "source": source_label,
                        "comment": note,
                    }
                )
            formula_rows.append(
                {
                    "group_id": group_id,
                    "antiviral": antiviral_norm,
                    "expression": "(" + " AND ".join(member_ids) + ")",
                    "phenotype": phenotype,
                    "fold_ic50": fold_ic50,
                    "publication": publication,
                    "source": source_label,
                    "comment": note,
                }
            )
            emitted_for_row += 1

        if not has_antiviral_data:
            add_non_migrated(non_migrated, row, "no_antiviral_signal")

    # Finalize aggregated single-mutation rules
    for aggr_key, aggr_data in aggregated_rules.items():
        # Calculate median IC50 if present
        ic50_values = aggr_data["ic50_values"]
        if ic50_values:
            ic50_values.sort()
            n = len(ic50_values)
            if n % 2 == 1:
                median_ic50 = ic50_values[n // 2]
            else:
                median_ic50 = (ic50_values[n // 2 - 1] + ic50_values[n // 2]) / 2.0
            fold_ic50 = f"{median_ic50:g}" if median_ic50 else ""
        else:
            fold_ic50 = ""

        # Determine phenotype: check for conflicts
        phenotypes = aggr_data["phenotypes"]
        phenotype = ""
        if phenotypes:
            unique_phenotypes = set(phenotypes)
            # Check for contradiction: resistant vs sensitive
            if "resistant" in unique_phenotypes and "sensitive" in unique_phenotypes:
                phenotype = "contradictory"
            else:
                # Use the first phenotype (or most common if needed)
                phenotype = phenotypes[0]

        # Join publications with comma
        publication = ",".join(p for p in aggr_data["publications"] if p)

        # Add comment about IC50 count if multiple values
        comment = ""
        if len(ic50_values) > 1:
            comment = f"{len(ic50_values)} IC50 values - displayed is median"

        rule = {
            "gene": aggr_data["gene"],
            "reference_identifier": aggr_data["reference_identifier"],
            "position": aggr_data["position"],
            "reference": aggr_data["reference"],
            "mutation": aggr_data["mutation"],
            "antiviral": aggr_data["antiviral"],
            "group_id": "",
            "member_id": "",
            "phenotype": phenotype,
            "fold_ic50": fold_ic50,
            "publication": publication,
            "source": aggr_data["source"],
            "comment": comment,
        }
        rules_rows.append(rule)

    rules_rows.sort(
        key=lambda r: (
            norm(r["reference_identifier"]),
            norm(r["gene"]),
            int(r["position"]),
            norm(r["mutation"]),
            norm(r["antiviral"]),
            norm(r["source"]),
        )
    )
    formula_rows.sort(key=lambda r: (norm(r["group_id"]), norm(r["antiviral"])))
    return rules_rows, formula_rows, non_migrated


def fetch_source_commit_date(source_url: str) -> str:
    """Query the GitHub API for the latest commit date of the source file.

    Returns the committer date as a YYYY-MM-DD string, or today's date as fallback.
    """
    # Parse owner/repo/path from a raw.githubusercontent.com URL
    # e.g. https://raw.githubusercontent.com/owner/repo/branch/path
    import re as _re

    match = _re.match(
        r"https://raw\.githubusercontent\.com/([^/]+)/([^/]+)/([^/]+)/(.+)",
        source_url,
    )
    if not match:
        print(
            f"WARNING: Cannot parse GitHub raw URL to fetch commit date: {source_url}",
            file=sys.stderr,
        )
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")

    owner, repo, branch, path = match.groups()
    api_url = (
        f"https://api.github.com/repos/{owner}/{repo}/commits"
        f"?path={path}&sha={branch}&per_page=1"
    )
    try:
        with urllib.request.urlopen(api_url) as response:
            commits = json.loads(response.read().decode("utf-8"))
        commit_date = commits[0]["commit"]["committer"]["date"][:10]
        print(f"Source commit date: {commit_date}", file=sys.stderr)
        return commit_date
    except Exception as exc:
        print(
            f"WARNING: Could not fetch commit date from GitHub API ({exc}); using today.",
            file=sys.stderr,
        )
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def download_source(source_url: str) -> str:
    with urllib.request.urlopen(source_url) as response:
        return response.read().decode("utf-8")


def validate_input_header(reader: csv.DictReader) -> None:
    header = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_COLUMNS - header)
    if missing:
        raise ValueError(f"Missing required source columns: {', '.join(missing)}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert HerpesDRG TSV to ResPro TSV artifacts")
    parser.add_argument("--source-url", default=SOURCE_URL, help="Upstream TSV URL")
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent.parent / "output"),
        help="Output directory for generated artifacts",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    maintainer_update = fetch_source_commit_date(args.source_url)
    print(f"Downloading {args.source_url} …", file=sys.stderr)
    source_text = download_source(args.source_url)

    reader = csv.DictReader(source_text.splitlines(), delimiter="\t")
    validate_input_header(reader)
    source_rows = list(reader)

    rules_rows, formula_rows, non_migrated_rows = convert(source_rows)

    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    rules_path = out_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")

    formula_path = out_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = tsv_from_rows(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
        print(f"Written {formula_path} ({len(formula_rows)} rows).", file=sys.stderr)
    else:
        formula_path.unlink(missing_ok=True)

    metadata = {
        "maintainers": ["Oscar Charles"],
        "contact": "oscar.charles.18@ucl.ac.uk",
        "publication_pmid": "39192205",
        "website": "https://github.com/ojcharles/herpesdrg-db",
        "description": "Comprehensive resource for human herpesvirus antiviral drug resistance genotyping.",
        "maintainer_update": maintainer_update,
        "license": "MIT",
        "tsv_checksum": checksum(rules_content),
    }
    metadata_path = out_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    non_migrated_content = non_migrated_text(non_migrated_rows)
    non_migrated_path = out_dir / "non-migrated-rules.txt"
    non_migrated_path.write_text(non_migrated_content, encoding="utf-8")

    print(f"Written {rules_path} ({len(rules_rows)} rows).", file=sys.stderr)
    print(f"Written {metadata_path}.", file=sys.stderr)
    print(
        f"Written {non_migrated_path} ({len(non_migrated_rows)} aggregated entries).",
        file=sys.stderr,
    )
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

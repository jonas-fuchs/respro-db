#!/usr/bin/env python3
"""
Convert HerpesDRG TSV into ResPro-compatible TSV artifacts.

Outputs:
- rules.tsv (required)
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
    "adeno_ad5": "AC_000008.1",
    "hsv1": "NC_001806",
    "hsv2": "NC_001798",
    "hhv6b": "MF511171",
}

VIRUS_ALIASES = {
    "hcmv": "hcmv",
    "human cytomegalovirus": "hcmv",
    "cmv": "hcmv",
    "vzv": "vzv",
    "varicella zoster virus": "vzv",
    "adeno": "adeno_ad5",
    "adeno ad5": "adeno_ad5",
    "adeno ad 5": "adeno_ad5",
    "adenovirus 5": "adeno_ad5",
    "adeno5": "adeno_ad5",
    "hsv1": "hsv1",
    "herpes simplex virus 1": "hsv1",
    "hsv-1": "hsv1",
    "hsv 1": "hsv1",
    "hsv2": "hsv2",
    "herpes simplex virus 2": "hsv2",
    "hsv-2": "hsv2",
    "hsv 2": "hsv2",
    "hhv6b": "hhv6b",
    "hhv-6b": "hhv6b",
    "hhv 6b": "hhv6b",
    "hhv6-b": "hhv6b",
}

EXCLUSION_NOTE_PATTERNS = [
    r"\banecdotal\b",
    r"cannot be inferred",
    r"cannot infer",
    r"ref wrong",
    r"review value does not match",
    r"unable to substantiate",
    r"unknown co[- ]?mut",
    r"co-occur",
    r"not resistance data",
    r"no data",
    r"ambig",
]


def norm(v: object) -> str:
    if v is None:
        return ""
    return str(v).strip()


def normalize_virus(raw_virus: str) -> str:
    folded = re.sub(r"\s+", " ", raw_virus.strip().lower())
    if folded in VIRUS_ALIASES:
        return VIRUS_ALIASES[folded]
    compact = re.sub(r"[^a-z0-9]+", "", folded)
    if compact in {"hcmv", "humancytomegalovirus", "cmv"}:
        return "hcmv"
    if compact in {"vzv", "varicellazostervirus"}:
        return "vzv"
    if compact in {"adeno", "adenoad5", "adeno5", "adenovirus5", "adenoadenoad5"}:
        return "adeno_ad5"
    if compact in {"hsv1", "herpessimplexvirus1", "hsvtype1"}:
        return "hsv1"
    if compact in {"hsv2", "herpessimplexvirus2", "hsvtype2"}:
        return "hsv2"
    if compact in {"hhv6b", "hhv6", "hhv6bvirus"}:
        return "hhv6b"
    return ""


def should_exclude(status: str, note: str) -> tuple[bool, str]:
    if status.strip().upper() != "A":
        return True, "status_not_active"
    note_l = note.lower()
    for pattern in EXCLUSION_NOTE_PATTERNS:
        if re.search(pattern, note_l):
            return True, "flagged_note"
    return False, ""


def parse_mutation(aa_change: str) -> tuple[str, str, int]:
    text = aa_change.strip().replace(" ", "")

    sub = re.fullmatch(r"([A-Za-z*])(\d+)([A-Za-z*])", text)
    if sub:
        ref = sub.group(1).upper()
        pos = int(sub.group(2))
        alt = sub.group(3).upper()
        return ref, "*" if alt == "*" else alt, pos

    stop_word = re.fullmatch(r"([A-Za-z])(\d+)(stop)", text, flags=re.IGNORECASE)
    if stop_word:
        return stop_word.group(1).upper(), "*", int(stop_word.group(2))

    fs = re.fullmatch(r"([A-Za-z])(\d+)(frameshift|fs.*)", text, flags=re.IGNORECASE)
    if fs:
        return fs.group(1).upper(), "fs", int(fs.group(2))

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

    pref_del = re.fullmatch(r"([A-Za-z])del(\d+)", text, flags=re.IGNORECASE)
    if pref_del:
        ref = pref_del.group(1).upper()
        pos = int(pref_del.group(2))
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
    if lowered in {"resistant", "resistance", "r"}:
        return "", "resistant"
    if lowered in {"sensitive", "susceptible", "s"}:
        return "", "sensitive"
    if lowered in {"intermediate", "i"}:
        return "", "intermediate"
    if lowered == "polymorphism":
        return "", "unknown"

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
    parts = []
    if doi:
        cleaned = doi if doi.lower().startswith("doi:") else f"doi:{doi}"
        parts.append(cleaned)
    if link:
        parts.append(link)
    return ",".join(parts)


def convert(source_rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    rules_rows = []
    formula_rows = []
    non_migrated = []
    seen_rule_keys = set()

    for row in source_rows:
        mutation_id = norm(row.get("mutation_id"))
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

        exclude, reason = should_exclude(status, note)
        if exclude:
            add_non_migrated(non_migrated, row, reason)
            continue

        if not gene:
            add_non_migrated(non_migrated, row, "missing_gene")
            continue
        if not aa_change:
            add_non_migrated(non_migrated, row, "missing_aa_change")
            continue

        try:
            ref_aa, mutation_token, pos = parse_mutation(aa_change)
        except ValueError as exc:
            add_non_migrated(non_migrated, row, str(exc))
            continue

        publication = build_publication(norm(row.get("ref_doi")), norm(row.get("ref_link")))
        source_label = f"HerpesDRG mutation_id={mutation_id}"

        emitted_for_row = 0
        for antiviral_column in ANTIVIRAL_COLUMNS:
            raw_value = norm(row.get(antiviral_column))
            if not raw_value:
                continue

            fold_ic50, phenotype = parse_fold_and_phenotype(raw_value)
            rule = {
                "gene": gene,
                "reference_identifier": REFERENCE_BY_VIRUS[virus_key],
                "position": pos,
                "reference": ref_aa,
                "mutation": mutation_token,
                "antiviral": antiviral_column.lower(),
                "group_id": "",
                "member_id": "",
                "phenotype": phenotype,
                "fold_ic50": fold_ic50,
                "publication": publication,
                "source": source_label,
                "comment": note,
            }

            dedupe_key = (
                rule["gene"],
                rule["reference_identifier"],
                rule["position"],
                rule["mutation"],
                rule["antiviral"],
                rule["publication"],
            )
            if dedupe_key in seen_rule_keys:
                continue
            seen_rule_keys.add(dedupe_key)
            rules_rows.append(rule)
            emitted_for_row += 1

        if emitted_for_row == 0:
            add_non_migrated(non_migrated, row, "no_antiviral_signal")

        # Conservative handling: complex co-mutation context is excluded from formula output.
        if norm(row.get("co_aa")) or norm(row.get("co_gene")):
            add_non_migrated(
                non_migrated,
                row,
                "co_mutation_context_excluded",
                "Co-mutation rows are logged conservatively; no formula group emitted.",
            )

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

    print(f"Downloading source TSV: {args.source_url}", file=sys.stderr)
    source_text = download_source(args.source_url)

    reader = csv.DictReader(source_text.splitlines(), delimiter="\t")
    validate_input_header(reader)
    source_rows = list(reader)

    rules_rows, formula_rows, non_migrated_rows = convert(source_rows)

    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    (out_dir / "rules.tsv").write_text(rules_content, encoding="utf-8")

    formula_path = out_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = tsv_from_rows(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
    else:
        formula_path.unlink(missing_ok=True)

    metadata = {
        "maintainers": ["HerpesDRG Team"],
        "contact": "https://github.com/ojcharles/herpesdrg-db",
        "publication_pmid": "",
        "website": "https://github.com/ojcharles/herpesdrg-db",
        "description": "Curated multi-virus antiviral resistance rules converted from HerpesDRG TSV.",
        "maintainer_update": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "license": "See source repository",
        "tsv_checksum": checksum(rules_content),
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    non_migrated_content = non_migrated_text(non_migrated_rows)
    (out_dir / "non-migrated-rules.txt").write_text(non_migrated_content, encoding="utf-8")

    print(f"rules.tsv rows: {len(rules_rows)}", file=sys.stderr)
    print(f"formula-rules.tsv rows: {len(formula_rows)}", file=sys.stderr)
    print(f"non-migrated entries: {len(non_migrated_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()

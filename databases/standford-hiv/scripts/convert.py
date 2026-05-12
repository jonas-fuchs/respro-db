#!/usr/bin/env python3
"""
Convert Stanford HIVDB ASI XML into ResPro-compatible TSV artifacts.

Outputs:
- rules.tsv
- formula-rules.tsv, only when grouped logical rules are emitted
- metadata.json
- non-migrated-rules.txt

Important mapping:
- Atomic score rule -> rules.tsv row with antiviral + score.
- Boolean rule -> atomic member rows in rules.tsv + expression row in formula-rules.tsv.
- Equal-score MAX(...) -> OR formula.
- Mixed-score MAX(...) -> non-migrated audit entry to avoid double-counting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_SOURCE_URL = (
    "https://raw.githubusercontent.com/hivdb/hivfacts/main/data/algorithms/"
    "HIVDB_latest.xml"
)

DEFAULT_MUTATION_TYPES_URL = (
    "https://raw.githubusercontent.com/hivdb/hivfacts/main/data/"
    "mutation-type-pairs_hiv1.csv"
)

REFERENCE_IDENTIFIER = "MN919177" # HIV-1 subtype B reference sequence used by Stanford HIVDB ASI

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
    "score",
    "ic50",
    "publication",
    "source",
    "comment",
]

FORMULA_COLUMNS = [
    "group_id",
    "antiviral",
    "expression",
    "phenotype",
    "score",
    "ic50",
    "publication",
    "source",
    "comment",
]

NON_MIGRATED_COLUMNS = [
    "reason",
    "drug",
    "gene",
    "score",
    "raw_rule",
    "details",
]

SOURCE_LABEL = "Stanford HIVDB"
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
RESERVED_EXPR_WORDS = {"AND", "OR", "NOT", "XOR"}

GENE_ORDER = {"CA": 0, "PR": 1, "RT": 2, "IN": 3}

# Stanford HIVDB uses short drug codes. Use canonical full names for better
# downstream metadata enrichment (e.g., PubChem lookup).
DRUG_NAME_MAP = {
    "3TC": "lamivudine",
    "ABC": "abacavir",
    "ATV/R": "atazanavir",
    "AZT": "zidovudine",
    "BIC": "bictegravir",
    "CAB": "cabotegravir",
    "D4T": "stavudine",
    "DDI": "didanosine",
    "DOR": "doravirine",
    "DPV": "dapivirine",
    "DRV/R": "darunavir",
    "DTG": "dolutegravir",
    "EFV": "efavirenz",
    "ETR": "etravirine",
    "EVG": "elvitegravir",
    "FPV/R": "fosamprenavir",
    "FTC": "emtricitabine",
    "IDV/R": "indinavir",
    "ISL": "islatravir",
    "LEN": "lenacapavir",
    "LPV/R": "lopinavir",
    "NFV": "nelfinavir",
    "NVP": "nevirapine",
    "RAL": "raltegravir",
    "RPV": "rilpivirine",
    "SQV/R": "saquinavir",
    "TDF": "tenofovir disoproxil",
    "TPV/R": "tipranavir",
}


@dataclass(frozen=True)
class SourceInfo:
    requested_url: str
    xml_url: str
    pointer_filename: str
    source_version: str
    source_date: str
    commit_date: str
    source_sha256: str


@dataclass(frozen=True)
class RuleItem:
    lhs: str
    score: str
    raw: str
    max_scope: str


@dataclass(frozen=True)
class MutationTerm:
    position: int
    reference_hint: str
    mutation: str


@dataclass(frozen=True)
class ParsedExpression:
    expression: str
    members: tuple[dict, ...]
    is_formula: bool


def norm(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def http_get_bytes(url: str) -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "respro-stanford-hivdb-converter/1.0",
            "Accept": "application/vnd.github+json, text/plain, */*",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def http_get_text(url: str) -> str:
    return http_get_bytes(url).decode("utf-8-sig")


def checksum_text(text: str) -> str:
    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def checksum_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def parse_finite_score(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if parsed != parsed or parsed in (float("inf"), float("-inf")):
        return None
    return parsed


def version_from_hivdb_filename(filename: str) -> str:
    match = re.fullmatch(r"HIVDB_(.+)\.xml", filename)
    if not match:
        raise ValueError(f"Could not extract HIVDB version from filename: {filename}")
    return match.group(1)


def versions_url_for(source_url: str) -> str:
    return urllib.parse.urljoin(source_url, "versions.json")


def resolve_source(source_url: str) -> SourceInfo:
    source_url = source_url.strip()

    if source_url.endswith("HIVDB_latest.xml"):
        pointer_filename = http_get_text(source_url).strip()
        if not re.fullmatch(r"HIVDB_.+\.xml", pointer_filename):
            raise ValueError(
                f"HIVDB_latest.xml returned unexpected content: {pointer_filename!r}"
            )
        xml_url = urllib.parse.urljoin(source_url, pointer_filename)
    else:
        pointer_filename = ""
        xml_url = source_url

    filename = Path(urllib.parse.urlparse(xml_url).path).name
    source_version = version_from_hivdb_filename(filename)

    versions = json.loads(http_get_text(versions_url_for(source_url)))
    source_date = ""
    for version, date, virus in versions.get("HIVDB", []):
        if version == source_version and virus == "HIV1":
            source_date = date
            break

    if not source_date:
        raise ValueError(f"Could not find HIVDB {source_version} HIV1 in versions.json")

    xml_bytes = http_get_bytes(xml_url)
    source_sha256 = checksum_bytes(xml_bytes)

    return SourceInfo(
        requested_url=source_url,
        xml_url=xml_url,
        pointer_filename=pointer_filename,
        source_version=source_version,
        source_date=source_date,
        commit_date=fetch_source_commit_date(xml_url),
        source_sha256=source_sha256,
    )


def fetch_source_commit_date(raw_url: str) -> str:
    match = re.match(
        r"https://raw\.githubusercontent\.com/([^/]+)/([^/]+)/([^/]+)/(.+)",
        raw_url,
    )
    if not match:
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")

    owner, repo, branch, path = match.groups()
    api_url = (
        f"https://api.github.com/repos/{owner}/{repo}/commits"
        f"?path={urllib.parse.quote(path)}&sha={urllib.parse.quote(branch)}&per_page=1"
    )

    try:
        commits = json.loads(http_get_text(api_url))
        return commits[0]["commit"]["committer"]["date"][:10]
    except Exception as exc:
        eprint(f"WARNING: Could not fetch GitHub commit date: {exc}")
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")


def text_of(parent: ET.Element, child_name: str, default: str = "") -> str:
    child = parent.find(child_name)
    if child is None or child.text is None:
        return default
    return child.text.strip()


def split_csv_text(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def build_drug_to_gene(root: ET.Element) -> dict[str, str]:
    definitions = root.find("DEFINITIONS")
    if definitions is None:
        raise ValueError("XML has no DEFINITIONS element")

    gene_to_classes: dict[str, list[str]] = {}
    class_to_drugs: dict[str, list[str]] = {}

    for gene_def in definitions.findall("GENE_DEFINITION"):
        gene = text_of(gene_def, "NAME")
        classes = split_csv_text(text_of(gene_def, "DRUGCLASSLIST"))
        if gene:
            gene_to_classes[gene] = classes

    for drug_class in definitions.findall("DRUGCLASS"):
        class_name = text_of(drug_class, "NAME")
        drugs = split_csv_text(text_of(drug_class, "DRUGLIST"))
        if class_name:
            class_to_drugs[class_name] = drugs

    drug_to_gene: dict[str, str] = {}
    for gene, classes in gene_to_classes.items():
        for class_name in classes:
            for drug in class_to_drugs.get(class_name, []):
                drug_to_gene[drug] = gene

    if not drug_to_gene:
        raise ValueError("Could not map drugs to genes from DEFINITIONS")

    return drug_to_gene


def strip_outer_parens(text: str) -> str:
    text = text.strip()
    if not (text.startswith("(") and text.endswith(")")):
        return text

    depth = 0
    for idx, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0 and idx != len(text) - 1:
                return text

    return text[1:-1].strip()


def split_top_level(text: str, sep: str = ",") -> list[str]:
    parts: list[str] = []
    start = 0
    paren = brace = bracket = 0

    for idx, char in enumerate(text):
        if char == "(":
            paren += 1
        elif char == ")":
            paren = max(0, paren - 1)
        elif char == "{":
            brace += 1
        elif char == "}":
            brace = max(0, brace - 1)
        elif char == "[":
            bracket += 1
        elif char == "]":
            bracket = max(0, bracket - 1)
        elif char == sep and paren == brace == bracket == 0:
            part = text[start:idx].strip()
            if part:
                parts.append(part)
            start = idx + 1

    tail = text[start:].strip()
    if tail:
        parts.append(tail)

    return parts


def split_score_item(part: str) -> tuple[str, str] | None:
    if "=>" not in part:
        return None
    lhs, score = part.split("=>", 1)
    return strip_outer_parens(lhs.strip()), score.strip()


def normalize_antiviral(drug: str) -> str:
    token = norm(drug).upper()
    if token in DRUG_NAME_MAP:
        return DRUG_NAME_MAP[token]
    return norm(drug).lower()


def extract_rule_items(condition: str) -> list[RuleItem]:
    condition = re.sub(r"\s+", " ", condition or "").strip()
    if not condition:
        return []

    if condition.upper().startswith("SCORE FROM"):
        condition = condition[len("SCORE FROM") :].strip()

    condition = strip_outer_parens(condition)
    items: list[RuleItem] = []
    max_counter = 0

    def visit(fragment: str, max_scope: str = "") -> None:
        nonlocal max_counter
        fragment = strip_outer_parens(fragment.strip())

        for part in split_top_level(fragment):
            part = part.strip()
            if not part:
                continue

            if part.upper().startswith("MAX"):
                max_counter += 1
                inner = part[3:].strip()
                inner = strip_outer_parens(inner)
                visit(inner, max_scope=f"MAX{max_counter:05d}")
                continue

            split_item = split_score_item(part)
            if split_item is None:
                continue

            lhs, score = split_item
            items.append(RuleItem(lhs=lhs, score=score, raw=part, max_scope=max_scope))

    visit(condition)
    return items


def expand_mutation_group(raw_group: str) -> list[str]:
    token = raw_group.strip()
    lower = token.lower()

    # Whole-token shorthands.
    if lower in {"insertion", "insert", "ins"}:
        return ["_"]
    if lower in {"deletion", "del"}:
        return ["-"]

    out = []
    for char in token:
        if char in {"_", "-", "*"}:
            out.append(char)
            continue

        if char.isalpha():
            # Stanford ASI may use lowercase d/i as deletion/insertion shorthand.
            if char == "d":
                out.append("-")
                continue
            if char == "i":
                out.append("_")
                continue

            aa = char.upper()
            if aa in CANONICAL_AA:
                out.append(aa)

    return out


TERM_PAREN_RE = re.compile(
    r"(?<![A-Za-z0-9])(?P<pos>\d{1,4})\(\s*(?P<not>NOT\s+)?(?P<aas>[A-Za-z*_\-/]+)\s*\)",
    re.IGNORECASE,
)

TERM_SIMPLE_RE = re.compile(
    r"(?<![A-Za-z0-9])(?P<ref>[A-Z*]?)(?P<pos>\d{1,4})(?P<aas>[A-Za-z*_\-/]+)(?![A-Za-z0-9])",
    re.IGNORECASE,
)


def parse_term_text(text: str) -> tuple[list[MutationTerm], bool]:
    text = text.strip()

    paren = TERM_PAREN_RE.fullmatch(text)
    if paren:
        position = int(paren.group("pos"))
        negative = bool(paren.group("not"))
        aas = expand_mutation_group(paren.group("aas"))
        return [
            MutationTerm(position=position, reference_hint="", mutation=aa)
            for aa in aas
        ], negative

    simple = TERM_SIMPLE_RE.fullmatch(text)
    if simple:
        reference_hint = (simple.group("ref") or "").upper()
        position = int(simple.group("pos"))
        aas_raw = simple.group("aas")

        aas = expand_mutation_group(aas_raw)

        terms = []
        for aa in aas:
            if reference_hint and aa == reference_hint:
                continue
            terms.append(
                MutationTerm(
                    position=position,
                    reference_hint=reference_hint,
                    mutation=aa,
                )
            )
        return terms, False

    return [], False


COMMENT_MUT_RE = re.compile(r"\b([A-Z])(\d{1,4})([A-Z*_\-/]+)\b")


def infer_reference_hints_from_comments(root: ET.Element) -> dict[tuple[str, int, str], str]:
    hints: dict[tuple[str, int, str], str] = {}

    for comment in root.findall(".//COMMENT_STRING"):
        comment_id = comment.attrib.get("id", "")
        gene_match = re.match(r"^(CA|PR|RT|IN)", comment_id)
        if not gene_match:
            continue

        gene = gene_match.group(1)
        text = text_of(comment, "TEXT")
        searchable = f"{comment_id} {text}"

        for match in COMMENT_MUT_RE.finditer(searchable):
            reference = match.group(1)
            position = int(match.group(2))
            for mutation in expand_mutation_group(match.group(3)):
                if mutation != reference:
                    hints[(gene, position, mutation)] = reference

    return hints


def infer_reference_by_position_from_mutation_types(
    mutation_types_text: str,
) -> dict[tuple[str, int], str]:
    reader = csv.DictReader(io.StringIO(mutation_types_text))
    if not reader.fieldnames:
        return {}

    required = {"gene", "position", "aas"}
    if not required.issubset(set(reader.fieldnames)):
        return {}

    seen: dict[tuple[str, int], set[str]] = defaultdict(set)

    for row in reader:
        gene = norm(row.get("gene"))
        position_text = norm(row.get("position"))
        aas = norm(row.get("aas"))

        if not gene or not position_text.isdigit():
            continue

        key = (gene, int(position_text))
        for aa in expand_mutation_group(aas):
            if aa in CANONICAL_AA:
                seen[key].add(aa)

    inferred = {}
    for key, observed in seen.items():
        missing = CANONICAL_AA - observed
        if len(missing) == 1:
            inferred[key] = next(iter(missing))

    return inferred


def choose_reference(
    gene: str,
    term: MutationTerm,
    by_mutation: dict[tuple[str, int, str], str],
    by_position: dict[tuple[str, int], str],
) -> str:
    if term.reference_hint and term.reference_hint in CANONICAL_AA:
        return term.reference_hint

    if (gene, term.position, term.mutation) in by_mutation:
        return by_mutation[(gene, term.position, term.mutation)]

    if (gene, term.position) in by_position:
        return by_position[(gene, term.position)]

    return ""


def make_group_id(drug: str, gene: str, raw_rule: str, score: str) -> str:
    seed = f"{drug}|{gene}|{raw_rule}|{score}".encode("utf-8")
    return "G" + hashlib.sha256(seed).hexdigest()[:12].upper()


def make_member_id(group_id: str, index: int) -> str:
    member_id = f"{group_id}_M{index:03d}"
    if member_id.upper() in RESERVED_EXPR_WORDS:
        member_id = f"{group_id}_MEMBER_{index:03d}"
    return member_id


def make_shared_member_id(gene: str, position: int, reference: str, mutation: str) -> str:
    seed = f"{gene}|{position}|{reference}|{mutation}".encode("utf-8")
    member_id = "M" + hashlib.sha256(seed).hexdigest()[:14].upper()
    if member_id.upper() in RESERVED_EXPR_WORDS:
        member_id = "M_" + member_id
    return member_id


def normalize_logic_chunk(chunk: str) -> str:
    chunk = chunk.replace("&&", " AND ")
    chunk = chunk.replace("||", " OR ")
    chunk = chunk.replace("+", " AND ")
    chunk = re.sub(r"\bAND\b", " AND ", chunk, flags=re.IGNORECASE)
    chunk = re.sub(r"\bOR\b", " OR ", chunk, flags=re.IGNORECASE)
    chunk = re.sub(r"\bNOT\b", " NOT ", chunk, flags=re.IGNORECASE)
    chunk = re.sub(r"\bXOR\b", " XOR ", chunk, flags=re.IGNORECASE)
    chunk = re.sub(r"\s+", " ", chunk)
    return chunk.strip()


def contains_unsupported_asi(lhs: str) -> bool:
    unsupported_markers = [
        "$",
        ">=",
        "<=",
        "!=",
        "==",
        " TO ",
        "FROM ",
        "SELECT ",
        "EXCEPT ",
        "ATLEAST",
        "ATLEAST",
        "numberOf",
        "count",
    ]
    upper = lhs.upper()
    return any(marker.upper() in upper for marker in unsupported_markers)


def parse_lhs_to_expression(
    lhs: str,
    gene: str,
    reference_by_mutation: dict[tuple[str, int, str], str],
    reference_by_position: dict[tuple[str, int], str],
) -> ParsedExpression | None:
    if contains_unsupported_asi(lhs):
        return None

    # First collect mutation-looking tokens.
    matches = []
    occupied = [False] * len(lhs)

    for regex in (TERM_PAREN_RE, TERM_SIMPLE_RE):
        for match in regex.finditer(lhs):
            if any(occupied[i] for i in range(match.start(), match.end())):
                continue
            for i in range(match.start(), match.end()):
                occupied[i] = True
            matches.append(match)

    matches.sort(key=lambda m: m.start())
    if not matches:
        return None

    pieces = []
    members = []
    last = 0
    saw_boolean_or_group = False

    for match in matches:
        chunk = normalize_logic_chunk(lhs[last : match.start()])
        if chunk:
            # Only boolean syntax and parentheses are allowed in non-mutation chunks.
            cleaned = re.sub(r"\b(AND|OR|NOT|XOR)\b", "", chunk, flags=re.IGNORECASE)
            cleaned = cleaned.replace("(", "").replace(")", "").strip()
            if cleaned:
                return None
            pieces.append(chunk)
            if any(word in chunk.upper().split() for word in RESERVED_EXPR_WORDS):
                saw_boolean_or_group = True

        token_text = match.group(0)
        terms, negative = parse_term_text(token_text)
        if not terms:
            return None

        member_ids = []
        for term in terms:
            reference = choose_reference(
                gene=gene,
                term=term,
                by_mutation=reference_by_mutation,
                by_position=reference_by_position,
            )
            if not reference:
                return None

            member_id = make_shared_member_id(
                gene=gene,
                position=term.position,
                reference=reference,
                mutation=term.mutation,
            )

            members.append(
                {
                    "gene": gene,
                    "reference_identifier": REFERENCE_IDENTIFIER,
                    "position": str(term.position),
                    "reference": reference,
                    "mutation": term.mutation,
                    "member_id": member_id,
                }
            )
            member_ids.append(member_id)

        if len(member_ids) == 1:
            expr = member_ids[0]
        else:
            expr = "(" + " OR ".join(member_ids) + ")"
            saw_boolean_or_group = True

        if negative:
            expr = f"(NOT {expr})"
            saw_boolean_or_group = True

        pieces.append(expr)
        last = match.end()

    tail = normalize_logic_chunk(lhs[last:])
    if tail:
        cleaned = re.sub(r"\b(AND|OR|NOT|XOR)\b", "", tail, flags=re.IGNORECASE)
        cleaned = cleaned.replace("(", "").replace(")", "").strip()
        if cleaned:
            return None
        pieces.append(tail)
        if any(word in tail.upper().split() for word in RESERVED_EXPR_WORDS):
            saw_boolean_or_group = True

    expression = " ".join(piece for piece in pieces if piece)
    expression = re.sub(r"\s+", " ", expression).strip()

    # If the expression is just one member with no boolean/grouping, it can be atomic.
    is_formula = saw_boolean_or_group or len(members) != 1 or expression != members[0]["member_id"]

    return ParsedExpression(
        expression=expression,
        members=tuple(members),
        is_formula=is_formula,
    )


def tsv_from_rows(rows: list[dict], columns: list[str]) -> str:
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(norm(row.get(col, "")) for col in columns))
    return "\n".join(lines) + "\n"


def non_migrated_text(rows: list[dict]) -> str:
    lines = [
        "# Non-migrated Stanford HIVDB rules",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
        "\t".join(NON_MIGRATED_COLUMNS),
    ]

    for row in sorted(
        rows,
        key=lambda r: (
            norm(r.get("reason")),
            norm(r.get("gene")),
            norm(r.get("drug")),
            norm(r.get("raw_rule")),
        ),
    ):
        lines.append("\t".join(norm(row.get(col, "")) for col in NON_MIGRATED_COLUMNS))

    return "\n".join(lines) + "\n"


def add_non_migrated(
    bucket: list[dict],
    reason: str,
    drug: str,
    gene: str,
    score: str,
    raw_rule: str,
    details: str = "",
) -> None:
    bucket.append(
        {
            "reason": reason,
            "drug": drug,
            "gene": gene,
            "score": score,
            "raw_rule": raw_rule,
            "details": details,
        }
    )


def emit_atomic_or_formula(
    parsed: ParsedExpression,
    drug: str,
    score: str,
    raw_rule: str,
    rules_rows: list[dict],
    member_registry: dict[str, dict],
    formula_rows: list[dict],
    formula_comment: str = "",
) -> None:
    drug_norm = normalize_antiviral(drug)

    if not parsed.is_formula:
        member = parsed.members[0]
        rules_rows.append(
            {
                "gene": member["gene"],
                "reference_identifier": member["reference_identifier"],
                "position": member["position"],
                "reference": member["reference"],
                "mutation": member["mutation"],
                "antiviral": drug_norm,
                "group_id": "",
                "member_id": "",
                "phenotype": "",
                "score": score,
                "ic50": "",
                "publication": "",
                "source": SOURCE_LABEL,
                "comment": "HIVDB ASI atomic score rule",
            }
        )
        return

    group_id_match = re.match(r"(G[A-F0-9]{12})_", parsed.members[0]["member_id"])
    group_id = (
        group_id_match.group(1)
        if group_id_match
        else make_group_id(drug, parsed.members[0]["gene"], raw_rule, score)
    )

    for member in parsed.members:
        existing = member_registry.get(member["member_id"])
        if existing is None:
            member_registry[member["member_id"]] = {
                "gene": member["gene"],
                "reference_identifier": member["reference_identifier"],
                "position": member["position"],
                "reference": member["reference"],
                "mutation": member["mutation"],
                "member_id": member["member_id"],
                "group_ids": {group_id},
            }
        else:
            key_existing = (
                existing["gene"],
                existing["reference_identifier"],
                existing["position"],
                existing["reference"],
                existing["mutation"],
            )
            key_new = (
                member["gene"],
                member["reference_identifier"],
                member["position"],
                member["reference"],
                member["mutation"],
            )
            if key_existing != key_new:
                raise ValueError(
                    f"Shared member_id collision for {member['member_id']}: {key_existing} vs {key_new}"
                )
            existing["group_ids"].add(group_id)

    formula_rows.append(
        {
            "group_id": group_id,
            "antiviral": drug_norm,
            "expression": parsed.expression,
            "phenotype": "",
            "score": score,
            "ic50": "",
            "publication": "",
            "source": SOURCE_LABEL,
            "comment": formula_comment or "HIVDB ASI formula score rule",
        }
    )


def convert(
    source_info: SourceInfo,
    xml_bytes: bytes,
    mutation_types_text: str,
    output_dir: Path,
) -> tuple[int, int, int]:
    root = ET.fromstring(xml_bytes)

    alg_name = text_of(root, "ALGNAME")
    alg_version = text_of(root, "ALGVERSION")
    alg_date = text_of(root, "ALGDATE")

    if alg_name != "HIVDB":
        raise ValueError(f"Expected ALGNAME=HIVDB, found {alg_name!r}")

    if alg_version != source_info.source_version:
        raise ValueError(
            f"Filename version {source_info.source_version!r} does not match "
            f"XML ALGVERSION {alg_version!r}"
        )

    if alg_date != source_info.source_date:
        raise ValueError(
            f"versions.json date {source_info.source_date!r} does not match "
            f"XML ALGDATE {alg_date!r}"
        )

    drug_to_gene = build_drug_to_gene(root)
    reference_by_mutation = infer_reference_hints_from_comments(root)
    reference_by_position = infer_reference_by_position_from_mutation_types(
        mutation_types_text
    )

    rules_rows: list[dict] = []
    member_registry: dict[str, dict] = {}
    formula_rows: list[dict] = []
    non_migrated: list[dict] = []
    flattened_mixed_max = 0

    for drug_element in root.findall("DRUG"):
        drug = text_of(drug_element, "NAME")
        gene = drug_to_gene.get(drug, "")

        if not drug or not gene:
            add_non_migrated(
                non_migrated,
                reason="drug_without_gene_mapping",
                drug=drug,
                gene=gene,
                score="",
                raw_rule="",
                details="Could not map DRUG/NAME to gene via DEFINITIONS",
            )
            continue

        rule_node = drug_element.find("./RULE")
        condition = text_of(rule_node if rule_node is not None else drug_element, "CONDITION")
        items = extract_rule_items(condition)

        if not items:
            add_non_migrated(
                non_migrated,
                reason="no_score_items_extracted",
                drug=drug,
                gene=gene,
                score="",
                raw_rule=condition[:500],
            )
            continue

        # Process non-MAX rules directly.
        max_scopes: dict[str, list[RuleItem]] = defaultdict(list)

        for item in items:
            if item.max_scope:
                max_scopes[item.max_scope].append(item)
                continue

            parsed = parse_lhs_to_expression(
                item.lhs,
                gene=gene,
                reference_by_mutation=reference_by_mutation,
                reference_by_position=reference_by_position,
            )

            if parsed is None:
                add_non_migrated(
                    non_migrated,
                    reason="unsupported_or_unresolved_expression",
                    drug=drug,
                    gene=gene,
                    score=item.score,
                    raw_rule=item.raw,
                )
                continue

            emit_atomic_or_formula(
                parsed=parsed,
                drug=drug,
                score=item.score,
                raw_rule=item.raw,
                rules_rows=rules_rows,
                member_registry=member_registry,
                formula_rows=formula_rows,
            )

        # Process MAX groups. Equal-score MAX can be represented as OR.
        for scope, scoped_items in max_scopes.items():
            scores = {item.score for item in scoped_items}
            if len(scores) != 1:
                # Safe flattening case: all MAX branches are simple single-mutation
                # substitutions at exactly the same site, which are mutually exclusive.
                parsed_singletons = []
                common_site = None
                flattenable = True

                for item in scoped_items:
                    parsed = parse_lhs_to_expression(
                        item.lhs,
                        gene=gene,
                        reference_by_mutation=reference_by_mutation,
                        reference_by_position=reference_by_position,
                    )
                    if parsed is None or parsed.is_formula or len(parsed.members) != 1:
                        flattenable = False
                        break

                    member = parsed.members[0]
                    if member["mutation"] not in CANONICAL_AA:
                        flattenable = False
                        break

                    score_value = parse_finite_score(item.score)
                    if score_value is None:
                        flattenable = False
                        break

                    site = (member["gene"], member["position"], member["reference"])
                    if common_site is None:
                        common_site = site
                    elif common_site != site:
                        flattenable = False
                        break

                    parsed_singletons.append(
                        {
                            "member": member,
                            "score_str": item.score,
                            "score_value": score_value,
                        }
                    )

                if flattenable and parsed_singletons:
                    best_by_mutation: dict[str, dict] = {}
                    for candidate in parsed_singletons:
                        mutation = candidate["member"]["mutation"]
                        previous = best_by_mutation.get(mutation)
                        if previous is None or candidate["score_value"] > previous["score_value"]:
                            best_by_mutation[mutation] = candidate

                    for mutation in sorted(best_by_mutation):
                        candidate = best_by_mutation[mutation]
                        member = candidate["member"]
                        rules_rows.append(
                            {
                                "gene": member["gene"],
                                "reference_identifier": member["reference_identifier"],
                                "position": member["position"],
                                "reference": member["reference"],
                                "mutation": member["mutation"],
                                "antiviral": normalize_antiviral(drug),
                                "group_id": "",
                                "member_id": "",
                                "phenotype": "",
                                "score": candidate["score_str"],
                                "ic50": "",
                                "publication": "",
                                "source": SOURCE_LABEL,
                                "comment": "HIVDB ASI mixed-score MAX flattened (mutually exclusive substitutions)",
                            }
                        )

                    flattened_mixed_max += 1
                    continue

                add_non_migrated(
                    non_migrated,
                    reason="mixed_score_max_not_representable",
                    drug=drug,
                    gene=gene,
                    score=",".join(sorted(scores)),
                    raw_rule=" | ".join(item.raw for item in scoped_items),
                    details="MAX with mixed scores cannot be represented with plain additive boolean formulas",
                )
                continue

            score = next(iter(scores))
            parsed_parts = []
            member_rows = []
            ok = True

            for item in scoped_items:
                parsed = parse_lhs_to_expression(
                    item.lhs,
                    gene=gene,
                    reference_by_mutation=reference_by_mutation,
                    reference_by_position=reference_by_position,
                )
                if parsed is None:
                    ok = False
                    add_non_migrated(
                        non_migrated,
                        reason="unsupported_expression_inside_max",
                        drug=drug,
                        gene=gene,
                        score=item.score,
                        raw_rule=item.raw,
                    )
                    break
                parsed_parts.append(parsed.expression)
                member_rows.extend(parsed.members)

            if not ok:
                continue

            # Deduplicate members by member_id.
            unique_members = []
            seen_member_ids = set()
            for member in member_rows:
                if member["member_id"] in seen_member_ids:
                    continue
                seen_member_ids.add(member["member_id"])
                unique_members.append(member)

            # Deduplicate repeated parsed branches while preserving order.
            unique_parts = []
            seen_parts = set()
            for expr in parsed_parts:
                if expr in seen_parts:
                    continue
                seen_parts.add(expr)
                unique_parts.append(expr)
            parsed_parts = unique_parts

            max_expression = "(" + " OR ".join(f"({expr})" for expr in parsed_parts) + ")"
            expression_tokens = [
                token
                for token in re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", max_expression)
                if token.upper() not in RESERVED_EXPR_WORDS
            ]
            token_counts = Counter(expression_tokens)
            duplicate_tokens = [t for t, c in token_counts.items() if c > 1]
            if duplicate_tokens:
                add_non_migrated(
                    non_migrated,
                    reason="duplicate_member_in_max_expression",
                    drug=drug,
                    gene=gene,
                    score=score,
                    raw_rule="MAX:" + "|".join(item.raw for item in scoped_items),
                    details="duplicate member IDs in generated expression: "
                    + ",".join(sorted(duplicate_tokens)),
                )
                continue

            max_parsed = ParsedExpression(
                expression=max_expression,
                members=tuple(unique_members),
                is_formula=True,
            )

            emit_atomic_or_formula(
                parsed=max_parsed,
                drug=drug,
                score=score,
                raw_rule="MAX:" + "|".join(item.raw for item in scoped_items),
                rules_rows=rules_rows,
                member_registry=member_registry,
                formula_rows=formula_rows,
                formula_comment="HIVDB ASI equal-score MAX represented as OR",
            )

    # Materialize shared formula members, attaching all group memberships.
    for item in member_registry.values():
        rules_rows.append(
            {
                "gene": item["gene"],
                "reference_identifier": item["reference_identifier"],
                "position": item["position"],
                "reference": item["reference"],
                "mutation": item["mutation"],
                "antiviral": "",
                "group_id": ",".join(sorted(item["group_ids"])),
                "member_id": item["member_id"],
                "phenotype": "",
                "score": "",
                "ic50": "",
                "publication": "",
                "source": SOURCE_LABEL,
                "comment": "HIVDB ASI formula member",
            }
        )

    # Deduplicate rows exactly by emitted TSV columns.
    def dedupe(rows: list[dict], columns: list[str]) -> list[dict]:
        seen = set()
        out = []
        for row in rows:
            key = tuple(norm(row.get(col, "")) for col in columns)
            if key in seen:
                continue
            seen.add(key)
            out.append(row)
        return out

    rules_rows = dedupe(rules_rows, RULES_COLUMNS)
    formula_rows = dedupe(formula_rows, FORMULA_COLUMNS)

    rules_rows.sort(
        key=lambda r: (
            GENE_ORDER.get(norm(r["gene"]), 99),
            norm(r["gene"]),
            int(norm(r["position"]) or "0"),
            norm(r["reference"]),
            norm(r["mutation"]),
            norm(r["antiviral"]),
            norm(r["group_id"]),
            norm(r["member_id"]),
            norm(r["score"]),
        )
    )

    formula_rows.sort(
        key=lambda r: (
            norm(r["group_id"]),
            norm(r["antiviral"]),
            norm(r["expression"]),
            norm(r["score"]),
        )
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    rules_path = output_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")

    formula_path = output_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = tsv_from_rows(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
        eprint(f"Written {formula_path} ({len(formula_rows)} rows).")
    else:
        formula_path.unlink(missing_ok=True)

    metadata = {
        "maintainers": ["Stanford HIVDB Team"],
        "contact": "",
        "publication_pmid": "",
        "website": "https://hivdb.stanford.edu/",
        "description": "Stanford HIVDB ASI drug-resistance scoring rules converted to ResPro TSV artifacts.",
        "maintainer_update": source_info.source_date,
        "license": "Unlicense",
        "tsv_checksum": checksum_text(rules_content),
    }

    metadata_path = output_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    non_migrated_content = non_migrated_text(non_migrated)
    non_migrated_path = output_dir / "non-migrated-rules.txt"
    non_migrated_path.write_text(non_migrated_content, encoding="utf-8")

    validate_outputs(rules_path, formula_path if formula_rows else None)

    eprint(f"Written {rules_path} ({len(rules_rows)} rows).")
    eprint(f"Written {metadata_path}.")
    eprint(
        f"Written {non_migrated_path} ({len(non_migrated)} aggregated entries)."
    )
    eprint(f"Flattened mixed-score MAX rules: {flattened_mixed_max}.")

    return len(rules_rows), len(formula_rows), len(non_migrated)


def validate_outputs(rules_path: Path, formula_path: Path | None) -> None:
    with rules_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != RULES_COLUMNS:
            raise ValueError(f"rules.tsv header mismatch: {reader.fieldnames}")

        member_ids = set()
        group_ids = set()

        for line_no, row in enumerate(reader, start=2):
            if set(row) != set(RULES_COLUMNS):
                raise ValueError(f"rules.tsv line {line_no} has wrong columns")
            if not row["gene"]:
                raise ValueError(f"rules.tsv line {line_no} missing gene")
            if not row["position"].isdigit():
                raise ValueError(f"rules.tsv line {line_no} has non-integer position")

            member_id = row["member_id"]
            group_id = row["group_id"]

            if member_id:
                if member_id.upper() in RESERVED_EXPR_WORDS:
                    raise ValueError(
                        f"rules.tsv line {line_no} member_id uses reserved keyword"
                    )
                if not group_id:
                    raise ValueError(
                        f"rules.tsv line {line_no} member_id without group_id"
                    )
                member_ids.add(member_id)
                for gid in [p.strip() for p in group_id.split(",") if p.strip()]:
                    group_ids.add(gid)

    if formula_path is None:
        return

    with formula_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != FORMULA_COLUMNS:
            raise ValueError(f"formula-rules.tsv header mismatch: {reader.fieldnames}")

        formula_group_ids = set()

        for line_no, row in enumerate(reader, start=2):
            group_id = row["group_id"]
            expression = row["expression"]

            if not group_id:
                raise ValueError(f"formula-rules.tsv line {line_no} missing group_id")
            if not expression:
                raise ValueError(f"formula-rules.tsv line {line_no} missing expression")
            if group_id not in group_ids:
                raise ValueError(
                    f"formula-rules.tsv line {line_no} references unknown group_id"
                )
            if group_id in formula_group_ids:
                raise ValueError(
                    f"formula-rules.tsv line {line_no} duplicates group_id {group_id}"
                )
            formula_group_ids.add(group_id)

            tokens = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", expression)
            seen_formula_members = set()
            for token in tokens:
                if token.upper() in RESERVED_EXPR_WORDS:
                    continue
                if token not in member_ids:
                    raise ValueError(
                        f"formula-rules.tsv line {line_no} references unknown member_id {token}"
                    )
                if token in seen_formula_members:
                    raise ValueError(
                        f"formula-rules.tsv line {line_no} repeats member_id {token}"
                    )
                seen_formula_members.add(token)

        missing_formula_groups = sorted(group_ids - formula_group_ids)
        if missing_formula_groups:
            raise ValueError(
                "formula-rules.tsv missing group_id entries for grouped rules: "
                + ",".join(missing_formula_groups)
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert Stanford HIVDB ASI XML to ResPro TSV artifacts"
    )
    parser.add_argument("--source-url", default=DEFAULT_SOURCE_URL)
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parent.parent / "output"),
    )
    parser.add_argument("--mutation-types-url", default=DEFAULT_MUTATION_TYPES_URL)
    args = parser.parse_args()

    eprint(f"Downloading {args.source_url} …")
    source_info = resolve_source(args.source_url)

    if source_info.xml_url != args.source_url:
        eprint(f"Resolved latest upstream source: {source_info.xml_url}")

    xml_bytes = http_get_bytes(source_info.xml_url)

    eprint(f"Downloading {args.mutation_types_url} …")
    mutation_types_text = http_get_text(args.mutation_types_url)

    convert(
        source_info=source_info,
        xml_bytes=xml_bytes,
        mutation_types_text=mutation_types_text,
        output_dir=Path(args.output_dir),
    )

    eprint("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
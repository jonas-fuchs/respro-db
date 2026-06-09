#!/usr/bin/env python3
"""
Convert HBV-GLUE drug resistance tabular data into ResPro-compatible TSV artifacts.

Source:
  https://raw.githubusercontent.com/giffordlabcvr/HBV_DRUG_RESISTANCE/master/tabular/DrugResistanceAlessaReformatted.txt

Columns in source: Drug, Mutation, Subdomain, Genotypes, Publications

Architecture (two-pass):
  Pass 1 — Collect all atomic mutations by position:
    Parse every source row. For each drug/position, aggregate all alternative AAs,
    publications, and genotypes. Emit one deduplicated rule per (drug, position,
    ref, mut) with merged publications. Comment format: "REFpos ALT1/ALT2/..."
    (e.g. "L526 M/I" for L526M and L526I sharing a position).

  Pass 2 — Collect all formula (combination) rules:
    Parse combo rows (containing + or ±). Both + and ± are treated as simple
    AND (this matches GLUE's interpretation). For bare positions in combos
    (e.g. ±M250), look up all known alt AAs from Pass 1 to create OR groups.
    Each combo row produces one formula rule with hash-based group_id and
    member_ids referencing the atomic rules from Pass 1.

  Bare positions (e.g. T184, M250):
    No longer sent to non-migrated. Instead, they are expanded using already-
    known mutations at that position for the same drug (from Pass 1).

Operator semantics:
  +  -> AND (co-mutation in formula)
  ,  -> OR  (separate independent rules, each parsed and aggregated)
  /  -> alternative amino acids at the same position
  ±  -> AND (treated identically to +, per GLUE semantics)

Dirty tokens:
  rt prefix  -> stripped (e.g. rtL80IV -> L80IV)
  r prefix   -> stripped with comment (e.g. rL180M -> L180M)
  S/C256G    -> genotype polymorphism: C is the NC_003977 reference, S appears
                in some genotypes; emit only C256G, note genotype refs in comment
  B236T      -> B is invalid ref AA; use actual ref N from NC_003977 lookup

Outputs:
  - rules.tsv          (atomic mutation rules; member_id is shared across formulas)
  - formula-rules.tsv  (combination rules with AND/OR expressions, hash-based group_id)
  - metadata.json      (provenance, interpretation algorithms, checksum)
  - non-migrated-rules.txt (audit trail for unparseable rows)

Member ID scheme:
  Content-based hash: M + first 14 chars of SHA256(feature|position|ref|mut).upper()
  Same mutation always gets the same member_id regardless of which drug/formula.
  member_id is emitted only once — on the first (position, ref, mut) occurrence
  across all drugs. Subsequent rows for the same mutation across different drugs
  have an empty member_id, avoiding duplicate member_id definitions.
  Only rules referenced by at least one formula get a member_id at all.

Group ID scheme:
  Content-based hash: G + first 12 chars of SHA256(drug|raw_notation|expression).upper()
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SOURCE_URL = (
    "https://raw.githubusercontent.com/giffordlabcvr/HBV_DRUG_RESISTANCE/"
    "master/tabular/DrugResistanceAlessaReformatted.txt"
)

FEATURE = "P"
REFERENCE_IDENTIFIER = "NC_003977"
SOURCE_LABEL = "HBV-GLUE"

# RT numbering offset: P_protein_position = rt_position + RT_OFFSET
# Verified by aligning YMDD motif (rt204 = P539, M at P539) with NC_003977.1 P protein.
RT_OFFSET = 335

CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Publication alias mapping
PUBLICATION_ALIASES: dict[str, str] = {
    "Wang_2017": "31600164",
}

# Drug code -> canonical name
DRUG_CODE_MAP: dict[str, str] = {
    "LAM": "lamivudine",
    "ADV": "adefovir",
    "TDF": "tenofovir",
    "ETV": "entecavir",
    "LDT": "telbivudine",
}

# Reference amino acid lookup for NC_003977 P protein (832 aa).
# Generated from GenBank NC_003977.1 CDS gene=P, translation table 1,
# extracted on 2026-06-09. Previous version used an incorrect (843 aa)
# sequence from a different isolate; corrected to match the RefSeq annotation.
# Key position verifications (rt numbering -> P position -> ref AA):
#   rt80 -> P415 -> L    rt126 -> P461 -> H   rt153 -> P488 -> R
#   rt166 -> P501 -> F   rt169 -> P504 -> I   rt173 -> P508 -> V
#   rt180 -> P515 -> L   rt181 -> P516 -> A   rt184 -> P519 -> T
#   rt191 -> P526 -> V   rt194 -> P529 -> A   rt200 -> P535 -> A
#   rt202 -> P537 -> S   rt204 -> P539 -> M   rt207 -> P542 -> V
#   rt213 -> P548 -> S   rt214 -> P549 -> V   rt215 -> P550 -> Q
#   rt217 -> P552 -> L   rt218 -> P553 -> E   rt221 -> P556 -> F
#   rt233 -> P568 -> I   rt236 -> P571 -> N   rt237 -> P572 -> P
#   rt238 -> P573 -> N   rt239 -> P574 -> K   rt245 -> P580 -> Y
#   rt246 -> P581 -> S   rt248 -> P583 -> N   rt250 -> P585 -> M
# B-prefixed: B236 -> N571 (ref=N at rt236), B238 -> N573 (ref=N at rt238)
# S/C256: S591 (ref=C at rt256 on NC_003977; genotype polymorphism S/C)

_REF_AA_STR = (
    "MPLSYQHFRRLLLLDDEAGPLEEELPRLADEGLNRRVAEDLNLGNLNVSIPWTHKVGNFTGLYSSTVPVF"
    "NPHWKTPSFPNIHLHQDIIKKCEQFVGPLTVNEKRRLQLIMPARFYPKVTKYLPLDKGIKPYYPEHLVNH"
    "YFQTRHYLHTLWKAGILYKRETTHSASFCGSPYSWEQDLQHGAESFHQQSSGILSRPPVGSSLQSKHRK"
    "SRLGLQSQQGHLARRQQGRSWSIRAGFHPTARRPFGVEPSGSGHTTNFASKSASCLHQSPVRKAAYPAVS"
    "TFEKHSSSGHAVEFHNLPPNSARSQSERPVFPCWWLQFRNSKPCSDYCLSLIVNLLEDWGPCAEHGEHHI"
    "RIPRTPSRVTGGVFLVDKNPHNTAESRLVVDFSQFSRGNYRVSWPKFAVPNLQSLTNLLSSNLSWLSLDV"
    "SAAFYHLPLHPAAMPHLLVGSSGLSRYVARLSSNSRILNNQHGTMPDLHDYCSRNLYVSLLLLYQTFGRK"
    "LHLYSHPIILGFRKIPMGVGLSPFLLAQFTSAICSVVRRAFPHCLAFSYMDDVVLGAKSVQHLESLFTA"
    "VTNFLLSLGIHLNPNKTKRWGYSLNFMGYVIGCYGSLPQEHIIQKIKECFRKLPINRPIDWKVCQRIVGL"
    "LGFAAPFTQCGYPALMPLYACIQSKQAFTFSPTYKAFLCKQYLNLYPVARQRPGLCQVFADATPTGWGLV"
    "MGHQRMRGTFSAPLPIHTAELLAACFARSRSGANIIGTDNSVVLSRKYTSFPWLLGCAANWILRGTSFVY"
    "VPSALNPADDPSRGRLGLSRPLLRLPFRPTTGRTSLYADSPSVPSHLPDRVHFASPLHVAWRPP"
)

REF_AA: dict[int, str] = {i: aa for i, aa in enumerate(_REF_AA_STR, start=1)}

# Output TSV columns
RULES_COLUMNS = [
    "feature",
    "reference_identifier",
    "position",
    "reference",
    "mutation",
    "antiviral",
    "phenotype",
    "publication",
    "source",
    "member_id",
    "comment",
]

FORMULA_COLUMNS = [
    "group_id",
    "antiviral",
    "expression",
    "phenotype",
    "publication",
    "source",
    "comment",
]

NON_MIGRATED_COLUMNS = [
    "reason",
    "drug",
    "source_mutation",
    "subdomain",
    "genotypes",
    "publications",
    "details",
]

RESERVED_EXPR_WORDS: frozenset[str] = frozenset({"AND", "OR", "NOT", "XOR"})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_member_id(feature: str, position: int, reference: str, mutation: str) -> str:
    """Stable, hash-based member ID shared across drugs for the same mutation."""
    seed = f"{feature}|{position}|{reference}|{mutation}".encode("utf-8")
    mid = "M" + hashlib.sha256(seed).hexdigest()[:14].upper()
    return mid if mid.upper() not in RESERVED_EXPR_WORDS else "M_" + mid


def _make_group_id(drug: str, raw_mutation: str, expression: str) -> str:
    """Stable, hash-based group ID for a formula rule."""
    seed = f"{drug}|{raw_mutation}|{expression}".encode("utf-8")
    return "G" + hashlib.sha256(seed).hexdigest()[:12].upper()


def norm(v: object) -> str:
    """Normalize a value to a stripped string, empty for None."""
    if v is None:
        return ""
    return str(v).strip()


def rt_to_p_position(rt_pos: int) -> int:
    """Convert HBV RT numbering position to P protein position on NC_003977."""
    return rt_pos + RT_OFFSET


def p_to_rt_position(p_pos: int) -> int:
    """Convert P protein position back to HBV RT numbering."""
    return p_pos - RT_OFFSET


def get_ref_aa(p_pos: int) -> str:
    """Get reference amino acid at P protein position, or '*' if out of range."""
    return REF_AA.get(p_pos, "*")


def parse_publication_token(token: str) -> str:
    """Parse a single publication token."""
    cleaned = token.strip()
    if not cleaned:
        return ""
    if cleaned in PUBLICATION_ALIASES:
        return f"PMID:{PUBLICATION_ALIASES[cleaned]}"
    if re.fullmatch(r"\d+", cleaned):
        return f"PMID:{cleaned}"
    return cleaned


def parse_publications(raw: str) -> list[str]:
    """Parse the Publications column into a list of publication strings."""
    if not raw or not raw.strip():
        return []
    tokens = [t.strip() for t in raw.split(",") if t.strip()]
    return [parse_publication_token(t) for t in tokens]


def join_pub_ids(ids: list[str]) -> str:
    """Join publication IDs, deduplicating and sorting by PMID number."""
    seen: set[str] = set()
    out: list[str] = []
    for v in ids:
        cleaned = norm(v)
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        out.append(cleaned)
    return ",".join(sorted(
        out,
        key=lambda x: (int(x.replace("PMID:", "")) if x.startswith("PMID:") and x[5:].isdigit() else 0, x)
    ))


# ---------------------------------------------------------------------------
# Atomic mutation
# ---------------------------------------------------------------------------


class Mut:
    """A single resolved mutation at a specific position."""

    __slots__ = ("position", "ref_aa", "mut_aa")

    def __init__(self, position: int, ref_aa: str, mut_aa: str):
        self.position = position
        self.ref_aa = ref_aa
        self.mut_aa = mut_aa

    def __repr__(self) -> str:
        return f"{self.ref_aa}{self.position}{self.mut_aa}"


class ParseFailure:
    """A mutation token that could not be migrated."""

    __slots__ = ("reason", "source_mutation", "details")

    def __init__(self, reason: str, source_mutation: str, details: str = ""):
        self.reason = reason
        self.source_mutation = source_mutation
        self.details = details


# ---------------------------------------------------------------------------
# Token parsing — single mutation token (no +, ±, , operators)
# ---------------------------------------------------------------------------


def clean_token(token: str) -> tuple[str, str]:
    """Strip rt/r dirty prefixes from a mutation token.

    Returns (cleaned_token, comment_or_empty).
    """
    comment = ""
    cleaned = token.strip()
    cleaned = re.sub(r"\s+", "", cleaned)

    # Strip 'rt' prefix (e.g. rtL80IV -> L80IV)
    if cleaned.lower().startswith("rt") and len(cleaned) > 2 and cleaned[2].isupper():
        cleaned = cleaned[2:]

    # Strip 'r' prefix (e.g. rL180M -> L180M)
    if len(cleaned) > 1 and cleaned[0] == "r" and cleaned[1].isupper():
        comment = f"dirty_prefix: original={token}"
        cleaned = cleaned[1:]

    return cleaned, comment


def parse_mutation_token(token: str, original: str) -> tuple[list[Mut], list[ParseFailure], list[str]]:
    """Parse a single cleaned mutation token (no +/±/, operators).

    Bare positions (e.g. T184, M250) are returned as Mut with mut_aa="?BARE"
    to signal the caller to expand using known alts from Pass 1 aggregation.

    Returns:
      mutations: list of Mut objects (may include ?BARE markers)
      failures: list of ParseFailure for unparseable tokens
      comments: list of comment strings
    """
    comments: list[str] = []
    failures: list[ParseFailure] = []

    if not token:
        return [], [], []

    # Pattern: S/C256G (genotype-polymorphism notation)
    m = re.fullmatch(r"([A-Z]+)/([A-Z])(\d+)([A-Z/]+)", token)
    if m:
        _slash_alts, _genotype_ref, pos_str, mut_str = m.groups()
        rt_pos = int(pos_str)
        p_pos = rt_to_p_position(rt_pos)
        ref_aa = get_ref_aa(p_pos)

        if p_pos < 1:
            failures.append(ParseFailure("position_out_of_range", original, f"P position={p_pos}"))
            return [], failures, comments

        if "/" in mut_str:
            mut_alts = [a for a in mut_str.split("/") if a]
        else:
            mut_alts = list(mut_str)

        genotype_refs = f"{_slash_alts}/{_genotype_ref}"
        comments.append(f"genotype_polymorphism: genome_ref={ref_aa}, genotype_refs={genotype_refs}")

        mutations: list[Mut] = []
        for alt in mut_alts:
            if alt not in CANONICAL_AA and alt != "*":
                comments.append(f"non_canonical_alt:{alt}")
                continue
            if alt == ref_aa:
                continue
            mutations.append(Mut(position=p_pos, ref_aa=ref_aa, mut_aa=alt))
        return mutations, failures, comments

    # Pattern: standard mutation with alternatives A181T, M204V/I, M204VI
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z/]+)", token)
    if m:
        src_ref, pos_str, mut_str = m.groups()
        rt_pos = int(pos_str)
        p_pos = rt_to_p_position(rt_pos)
        ref_aa = get_ref_aa(p_pos)

        if src_ref not in CANONICAL_AA and src_ref != "*":
            if p_pos in REF_AA:
                comments.append(f"invalid_ref_aa_corrected: source={src_ref}, genome_ref={ref_aa}")
            else:
                failures.append(ParseFailure("invalid_reference_aa", original,
                                              f"ref AA {src_ref} not canonical and P pos {p_pos} not found"))
                return [], failures, comments
        elif src_ref in CANONICAL_AA and src_ref != ref_aa and p_pos in REF_AA:
            comments.append(f"ref_aa_corrected: source={src_ref}, genome_ref={ref_aa}")

        if p_pos < 1:
            failures.append(ParseFailure("position_out_of_range", original, f"P position={p_pos}"))
            return [], failures, comments

        if "/" in mut_str:
            alts = [a for a in mut_str.split("/") if a]
        else:
            alts = list(mut_str)

        mutations: list[Mut] = []
        for alt in alts:
            if alt not in CANONICAL_AA and alt != "*":
                comments.append(f"non_canonical_alt:{alt}")
                continue
            if alt == ref_aa:
                continue
            mutations.append(Mut(position=p_pos, ref_aa=ref_aa, mut_aa=alt))
        return mutations, failures, comments

    # Pattern: bare position (e.g. T184, M250, S202)
    # Returns a special Mut with mut_aa="?BARE" to signal bare position.
    # The caller expands using known alts from the PositionIndex.
    m = re.fullmatch(r"([A-Z])(\d+)", token)
    if m:
        src_ref, pos_str = m.groups()
        rt_pos = int(pos_str)
        p_pos = rt_to_p_position(rt_pos)
        ref_aa = get_ref_aa(p_pos)
        comments.append(f"bare_position: rt{pos_str}=P{p_pos}, ref={ref_aa}, source_ref={src_ref}")
        return [Mut(position=p_pos, ref_aa=ref_aa, mut_aa="?BARE")], failures, comments

    # Unparseable
    failures.append(ParseFailure("unparseable_mutation_token", original,
                                  f"Cannot parse: {token}"))
    return [], failures, comments


# ---------------------------------------------------------------------------
# Position index: aggregates all known mutations per (drug, position)
# ---------------------------------------------------------------------------

PositionKey = tuple[str, int]  # (drug_name, p_position)


class PositionIndex:
    """Aggregates all known mutations per (drug, position) from Pass 1.

    Used to:
      - Deduplicate/merge atomic rules
      - Provide alt AA lists for expanding bare positions in Pass 2
    """

    def __init__(self) -> None:
        self._data: dict[PositionKey, dict] = {}

    def add_mutation(
        self,
        drug: str,
        p_pos: int,
        ref_aa: str,
        mut_aa: str,
        publications: list[str],
        genotypes: str,
        source_notation: str,
        comments: list[str],
    ) -> None:
        """Register a single (drug, position, ref, mut) observation."""
        key = (drug, p_pos)
        if key not in self._data:
            self._data[key] = {
                "ref_aa": ref_aa,
                "alts": set(),
                "publications": [],
                "genotypes": set(),
                "source_notations": [],
                "comments": [],
            }
        entry = self._data[key]
        if entry["ref_aa"] != ref_aa:
            entry["comments"].append(f"ref_aa_conflict: expected={entry['ref_aa']}, got={ref_aa}")
        entry["alts"].add(mut_aa)
        entry["publications"].extend(publications)
        if genotypes:
            for g in genotypes.split(","):
                g = g.strip()
                if g:
                    entry["genotypes"].add(g)
        if source_notation:
            entry["source_notations"].append(source_notation)
        entry["comments"].extend(comments)

    def get_alts(self, drug: str, p_pos: int) -> list[str]:
        """Get all known alternative AAs for a (drug, position)."""
        key = (drug, p_pos)
        if key in self._data:
            return sorted(self._data[key]["alts"])
        return []

    def get_ref_aa(self, drug: str, p_pos: int) -> str:
        """Get the reference AA for a (drug, position)."""
        key = (drug, p_pos)
        if key in self._data:
            return self._data[key]["ref_aa"]
        return get_ref_aa(p_pos)

    def build_comment(self, drug: str, p_pos: int, extra: str = "") -> str:
        """Build the comment using RT notation: rtREFpos ALT1/ALT2/..."""
        key = (drug, p_pos)
        if key not in self._data:
            return extra
        entry = self._data[key]
        ref = entry["ref_aa"]
        alts = sorted(entry["alts"])
        rt_pos = p_to_rt_position(p_pos)
        parts = [f"rt{ref}{rt_pos} {'/'.join(alts)}"]
        if entry["genotypes"]:
            parts.append(f"genotypes:{','.join(sorted(entry['genotypes']))}")
        # Include all ref-correction and genotype-polymorphism comments
        special_comments = [c for c in entry["comments"]
                           if c.startswith("genotype_polymorphism")
                           or c.startswith("invalid_ref")
                           or c.startswith("ref_aa_corrected")]
        # Deduplicate special comments (same message can appear multiple times)
        seen_special: set[str] = set()
        for c in special_comments:
            if c not in seen_special:
                parts.append(c)
                seen_special.add(c)
        if extra:
            parts.append(extra)
        return "; ".join(parts)

    def get_publications(self, drug: str, p_pos: int) -> list[str]:
        """Get all publications for a (drug, position)."""
        key = (drug, p_pos)
        if key in self._data:
            return list(self._data[key]["publications"])
        return []

    def items(self) -> list[tuple[PositionKey, dict]]:
        """Return all (key, entry) pairs sorted by key."""
        return sorted(self._data.items())


# ---------------------------------------------------------------------------
# Combo tokenization
# ---------------------------------------------------------------------------


def _tokenize_combo(raw: str) -> list[tuple[str, str, str]]:
    """Tokenize a combo mutation string into (operator, cleaned_token, original) triples.

    The first token always has op='+' (mandatory).
    Both '+' and '±' are returned as-is for operator tracking, but both
    are treated as AND in Pass 2.
    """
    s = raw.strip()
    s = re.sub(r"\s*\+\s*", "+", s)
    s = re.sub(r"\s*±\s*", "±", s)

    parts = re.split(r"([+±])", s)
    parts = [p for p in parts if p.strip()]

    tokens: list[tuple[str, str, str]] = []
    i = 0
    while i < len(parts):
        part = parts[i]
        if part in ("+", "±"):
            if i + 1 < len(parts):
                next_part = parts[i + 1]
                cleaned, _ = clean_token(next_part)
                tokens.append((part, cleaned, next_part))
                i += 2
            else:
                i += 1
        else:
            cleaned, _ = clean_token(part)
            tokens.append(("+", cleaned, part))
            i += 1

    return tokens


# ---------------------------------------------------------------------------
# Pass 1: Collect atomic mutations from all source rows
# ---------------------------------------------------------------------------


def pass1_collect_atomic(source_rows: list[dict]) -> tuple[PositionIndex, list[ParseFailure]]:
    """Pass through all source rows, collecting every atomic mutation per position/drug.

    Returns:
      pos_index: PositionIndex with all known mutations aggregated
      failures: list of ParseFailure for unparseable tokens
    """
    pos_index = PositionIndex()
    all_failures: list[ParseFailure] = []

    for row in source_rows:
        drug_code = norm(row.get("Drug"))
        raw_mutation = norm(row.get("Mutation"))
        genotypes = norm(row.get("Genotypes"))
        raw_publications = norm(row.get("Publications"))

        if not drug_code or not raw_mutation:
            all_failures.append(ParseFailure(
                "missing_drug_or_mutation", f"{drug_code}:{raw_mutation}",
                f"drug={drug_code}, mutation={raw_mutation}",
            ))
            continue

        drug_name = DRUG_CODE_MAP.get(drug_code.upper(), "")
        if not drug_name:
            all_failures.append(ParseFailure(
                "unknown_drug_code", raw_mutation, f"code={drug_code}",
            ))
            continue

        publications = parse_publications(raw_publications)

        # Split on commas -> each alternative is independent
        alternatives = [a.strip() for a in raw_mutation.split(",") if a.strip()]

        for alt in alternatives:
            has_combo = "+" in alt or "±" in alt
            if has_combo:
                # Combos are handled in Pass 2; but we still collect their
                # individual member mutations to populate the position index
                _collect_combo_members(alt, drug_name, genotypes, publications, pos_index, all_failures)
            else:
                # Simple mutation (possibly with / alternatives)
                cleaned, _ = clean_token(alt)
                muts, failures, token_comments = parse_mutation_token(cleaned, alt)
                all_failures.extend(failures)
                for m in muts:
                    if m.mut_aa == "?BARE":
                        # Bare position — will be expanded later from position index
                        continue
                    pos_index.add_mutation(
                        drug=drug_name,
                        p_pos=m.position,
                        ref_aa=m.ref_aa,
                        mut_aa=m.mut_aa,
                        publications=publications,
                        genotypes=genotypes,
                        source_notation=alt,
                        comments=token_comments,
                    )

    return pos_index, all_failures


def _collect_combo_members(
    raw_combo: str,
    drug_name: str,
    genotypes: str,
    publications: list[str],
    pos_index: PositionIndex,
    all_failures: list[ParseFailure],
) -> None:
    """Collect individual mutations from a combo expression into the position index.

    This ensures that mutations appearing only in combos are still available
    for bare-position expansion in Pass 2.
    """
    tokens = _tokenize_combo(raw_combo)
    for _op, cleaned, original in tokens:
        muts, failures, token_comments = parse_mutation_token(cleaned, original)
        all_failures.extend(failures)
        for m in muts:
            if m.mut_aa == "?BARE":
                continue  # Will be expanded later
            pos_index.add_mutation(
                drug=drug_name,
                p_pos=m.position,
                ref_aa=m.ref_aa,
                mut_aa=m.mut_aa,
                publications=publications,
                genotypes=genotypes,
                source_notation=original,
                comments=token_comments,
            )


# ---------------------------------------------------------------------------
# Pass 2: Build rules.tsv and formula-rules.tsv
# ---------------------------------------------------------------------------


def _extract_member_ids_from_expression(expression: str) -> set[str]:
    """Extract all member_id tokens from a boolean expression string."""
    tokens = expression.replace("(", "").replace(")", "").split()
    reserved = {"AND", "OR", "NOT", "XOR"}
    return {t for t in tokens if t not in reserved}


def pass2_build_output(
    source_rows: list[dict],
    pos_index: PositionIndex,
) -> tuple[list[dict], list[dict], list[dict]]:
    """Build the final rules, formulas, and non-migrated rows.

    Uses the PositionIndex from Pass 1 to:
      - Emit deduplicated atomic rules
      - Expand bare positions using known alts
      - Build formula expressions with content-based member_ids
      - Only assign member_id to atomic rules referenced by formulas

    Architecture:
      Step 1: Build formula rules and collect all referenced member_ids.
      Step 2: Build atomic rules, assigning member_id only when needed by formulas.
      Step 3: Add placeholder rows for formula-referenced members missing from the index.
    """
    formula_rows: list[dict] = []
    non_migrated: list[dict] = []

    # Track member info for formula-referenced members (used for placeholders)
    formula_member_registry: dict[str, dict] = {}
    # The set of member_ids actually referenced by at least one formula
    formula_referenced_mids: set[str] = set()

    # =========================================================================
    # Step 1: Build formula rules
    # =========================================================================
    seen_formulas: set[str] = set()  # Deduplicate by (drug, expression)

    for row in source_rows:
        drug_code = norm(row.get("Drug"))
        raw_mutation = norm(row.get("Mutation"))
        genotypes = norm(row.get("Genotypes"))
        raw_publications = norm(row.get("Publications"))

        if not drug_code or not raw_mutation:
            continue

        drug_name = DRUG_CODE_MAP.get(drug_code.upper(), "")
        if not drug_name:
            continue

        publications = parse_publications(raw_publications)
        publication_str = join_pub_ids(publications)

        alternatives = [a.strip() for a in raw_mutation.split(",") if a.strip()]

        for alt in alternatives:
            has_combo = "+" in alt or "±" in alt
            if not has_combo:
                continue

            # --- Combo mutation ---
            tokens = _tokenize_combo(alt)
            if not tokens:
                continue

            member_entries: list[tuple[str, Mut]] = []
            expression_parts: list[str] = []
            skip = False

            for _op, cleaned, original in tokens:
                muts, failures, _ = parse_mutation_token(cleaned, original)

                # Emit parse failures as non-migrated
                for f in failures:
                    non_migrated.append({
                        "reason": f.reason,
                        "drug": drug_code,
                        "source_mutation": alt,
                        "subdomain": norm(row.get("Subdomain")),
                        "genotypes": genotypes,
                        "publications": raw_publications,
                        "details": f.details,
                    })

                # Separate bare and normal mutations
                bare_muts = [m for m in muts if m.mut_aa == "?BARE"]
                real_muts = [m for m in muts if m.mut_aa != "?BARE"]

                # Handle bare positions: expand using known alts from Pass 1
                if bare_muts:
                    bm = bare_muts[0]  # There should be exactly one per token
                    known_alts = pos_index.get_alts(drug_name, bm.position)
                    if not known_alts:
                        non_migrated.append({
                            "reason": "bare_position_no_known_alts",
                            "drug": drug_code,
                            "source_mutation": original,
                            "subdomain": norm(row.get("Subdomain")),
                            "genotypes": genotypes,
                            "publications": raw_publications,
                            "details": f"P{bm.position} has no known alts for {drug_name}",
                        })
                        skip = True
                        break
                    # Create OR group of all known alts at this position
                    or_mids: list[str] = []
                    for a in known_alts:
                        mid = _make_member_id(FEATURE, bm.position, bm.ref_aa, a)
                        member_entries.append((mid, Mut(bm.position, bm.ref_aa, a)))
                        or_mids.append(mid)
                        if mid not in formula_member_registry:
                            formula_member_registry[mid] = {
                                "feature": FEATURE,
                                "reference_identifier": REFERENCE_IDENTIFIER,
                                "position": bm.position,
                                "reference": bm.ref_aa,
                                "mutation": a,
                                "member_id": mid,
                            }
                    if len(or_mids) == 1:
                        expression_parts.append(or_mids[0])
                    else:
                        expression_parts.append("(" + " OR ".join(or_mids) + ")")
                    continue

                if not real_muts:
                    # No viable mutations from this token
                    skip = True
                    break

                if len(real_muts) == 1:
                    m = real_muts[0]
                    mid = _make_member_id(FEATURE, m.position, m.ref_aa, m.mut_aa)
                    member_entries.append((mid, m))
                    expression_parts.append(mid)
                else:
                    # Multiple alternative AAs at same position -> OR
                    or_mids: list[str] = []
                    for m in real_muts:
                        mid = _make_member_id(FEATURE, m.position, m.ref_aa, m.mut_aa)
                        member_entries.append((mid, m))
                        or_mids.append(mid)
                    if len(or_mids) == 1:
                        expression_parts.append(or_mids[0])
                    else:
                        expression_parts.append("(" + " OR ".join(or_mids) + ")")

            if skip or not member_entries:
                continue

            # Build AND expression (both + and ± treated the same)
            expression = " AND ".join(expression_parts)
            gid = _make_group_id(drug_name, alt, expression)

            # Deduplicate formulas by (drug, expression)
            formula_key = f"{drug_name}|{expression}"
            if formula_key in seen_formulas:
                continue
            seen_formulas.add(formula_key)

            # Collect referenced member_ids from this expression
            formula_referenced_mids.update(_extract_member_ids_from_expression(expression))

            # Register all member entries in the registry
            for mid, m in member_entries:
                if mid not in formula_member_registry:
                    formula_member_registry[mid] = {
                        "feature": FEATURE,
                        "reference_identifier": REFERENCE_IDENTIFIER,
                        "position": m.position,
                        "reference": m.ref_aa,
                        "mutation": m.mut_aa,
                        "member_id": mid,
                    }

            # Build comment
            comment_parts: list[str] = [alt]
            if genotypes:
                comment_parts.append(f"genotypes:{genotypes}")
            formula_comment = "; ".join(comment_parts)

            formula_rows.append({
                "group_id": gid,
                "antiviral": drug_name,
                "expression": expression,
                "phenotype": "resistant",
                "publication": publication_str,
                "source": SOURCE_LABEL,
                "comment": formula_comment,
            })

    # =========================================================================
    # Step 2: Build atomic rules from position index
    # Assign member_id only once per (position, ref, mut) — the first drug
    # that emits a rule for a given mutation gets the member_id; subsequent
    # drugs get an empty member_id to avoid duplicate definitions.
    # Only assign member_id when the mutation is referenced by a formula.
    # =========================================================================
    rules_rows: list[dict] = []
    emitted_member_ids: set[str] = set()  # Track which member_ids already have a row

    for (drug, p_pos), entry in pos_index.items():
        ref_aa = entry["ref_aa"]
        publications = join_pub_ids(entry["publications"])
        for mut_aa in sorted(entry["alts"]):
            mid = _make_member_id(FEATURE, p_pos, ref_aa, mut_aa)
            comment = pos_index.build_comment(drug, p_pos)
            # Only emit member_id when referenced by a formula AND not yet emitted
            if mid in formula_referenced_mids and mid not in emitted_member_ids:
                effective_mid = mid
                emitted_member_ids.add(mid)
            else:
                effective_mid = ""

            rules_rows.append({
                "feature": FEATURE,
                "reference_identifier": REFERENCE_IDENTIFIER,
                "position": p_pos,
                "reference": ref_aa,
                "mutation": mut_aa,
                "antiviral": drug,
                "phenotype": "resistant",
                "publication": publications,
                "source": SOURCE_LABEL,
                "member_id": effective_mid,
                "comment": comment,
            })

    # =========================================================================
    # Step 3: Add placeholder rows for formula-referenced members
    # not already present in the position index (e.g. bare-position expansions
    # for a different drug).
    # =========================================================================
    emitted_positions: set[tuple] = set()
    for r in rules_rows:
        emitted_positions.add((r["position"], r["reference"], r["mutation"]))

    for mid, info in formula_member_registry.items():
        pos_key = (info["position"], info["reference"], info["mutation"])
        if pos_key not in emitted_positions:
            rules_rows.append({
                "feature": info["feature"],
                "reference_identifier": info["reference_identifier"],
                "position": info["position"],
                "reference": info["reference"],
                "mutation": info["mutation"],
                "antiviral": "",
                "phenotype": "",
                "publication": "",
                "source": SOURCE_LABEL,
                "member_id": info["member_id"],
                "comment": "",
            })
            emitted_positions.add(pos_key)

    # --- Sort deterministically ---
    rules_rows.sort(key=lambda r: (
        r["antiviral"], r["position"], r["mutation"], r.get("member_id", ""),
    ))
    formula_rows.sort(key=lambda r: (r["antiviral"], r["group_id"]))

    return rules_rows, formula_rows, non_migrated


# ---------------------------------------------------------------------------
# Convert (top-level orchestrator)
# ---------------------------------------------------------------------------


def convert(source_rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    """Convert source rows into (rules_rows, formula_rows, non_migrated).

    Two-pass architecture:
      Pass 1: Collect all atomic mutations by position/drug
      Pass 2: Build deduplicated rules and formulas
    """
    # Pass 1: Collect all mutations per position
    pos_index, parse_failures = pass1_collect_atomic(source_rows)

    # Convert parse failures to non-migrated format
    non_migrated_from_pass1: list[dict] = []
    for f in parse_failures:
        non_migrated_from_pass1.append({
            "reason": f.reason,
            "drug": "",
            "source_mutation": f.source_mutation,
            "subdomain": "",
            "genotypes": "",
            "publications": "",
            "details": f.details,
        })

    # Pass 2: Build final output
    rules_rows, formula_rows, non_migrated_from_pass2 = pass2_build_output(source_rows, pos_index)

    # Merge non-migrated
    all_non_migrated = non_migrated_from_pass1 + non_migrated_from_pass2

    return rules_rows, formula_rows, all_non_migrated


# ---------------------------------------------------------------------------
# TSV output helpers
# ---------------------------------------------------------------------------


def tsv_from_rows(rows: list[dict], columns: list[str]) -> str:
    """Build a TSV string from a list of row dicts and column names."""
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(norm(row.get(col, "")) for col in columns))
    return "\n".join(lines) + "\n"


def checksum(content: str) -> str:
    """Compute sha256 checksum of TSV content."""
    return "sha256:" + hashlib.sha256(content.encode("utf-8")).hexdigest()


def non_migrated_text(rows: list[dict]) -> str:
    """Build the non-migrated-rules audit file text."""
    lines = [
        "# Non-migrated HBV-GLUE rows",
        "# Columns: " + "\t".join(NON_MIGRATED_COLUMNS),
        "",
        "\t".join(NON_MIGRATED_COLUMNS),
    ]
    for row in sorted(
        rows,
        key=lambda r: (
            norm(r.get("reason")),
            norm(r.get("drug")),
            norm(r.get("source_mutation")),
        ),
    ):
        lines.append("\t".join(norm(row.get(col, "")) for col in NON_MIGRATED_COLUMNS))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Source fetching and parsing
# ---------------------------------------------------------------------------


def fetch_source(url: str) -> str:
    """Fetch the source TSV content from a URL."""
    with urllib.request.urlopen(url) as resp:
        return resp.read().decode("utf-8")


def parse_source(content: str) -> list[dict]:
    """Parse the source TSV content into a list of row dicts."""
    lines = content.strip().splitlines()
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames is None:
        reader = csv.DictReader(lines, delimiter=",")
    return [dict(row) for row in reader]


def eprint(msg: str) -> None:
    """Print a message to stderr."""
    print(msg, file=sys.stderr)


def fetch_source_commit_date(source_url: str) -> str:
    """Query the GitHub API for the latest commit date of the source file.

    Returns the committer date as a YYYY-MM-DD string, or today's date as fallback.
    """
    # Parse owner/repo/path from a raw.githubusercontent.com URL
    # e.g. https://raw.githubusercontent.com/owner/repo/branch/path
    match = re.match(
        r"https://raw\.githubusercontent\.com/([^/]+)/([^/]+)/([^/]+)/(.+)",
        source_url,
    )
    if not match:
        eprint(
            f"WARNING: Cannot parse GitHub raw URL to fetch commit date: {source_url}"
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
        return commit_date
    except Exception as exc:
        eprint(
            f"WARNING: Could not fetch commit date from GitHub API ({exc}); using today."
        )
        return datetime.now(timezone.utc).strftime("%Y-%m-%d")


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------


def build_metadata(output_dir: Path, source_date: str) -> dict:
    """Build the metadata.json content."""
    rules_path = output_dir / "rules.tsv"
    formula_path = output_dir / "formula-rules.tsv"

    tsv_content = ""
    if rules_path.exists():
        tsv_content += rules_path.read_text(encoding="utf-8")
    if formula_path.exists():
        tsv_content += formula_path.read_text(encoding="utf-8")

    tsv_checksum = checksum(tsv_content)

    return {
        "maintainers": ["Robert J. Gifford"],
        "contact": "Robert.Gifford@glasgow.ac.uk",
        "publication_pmid": "30563445",
        "website": "https://hbv-glue.cvr.gla.ac.uk/#/home",
        "description": (
            "HBV polymerase drug resistance mutations curated from the "
            "HBV-GLUE drug resistance extension mapped to reference NC_003977."
        ),
        "maintainer_update": source_date,
        "license": "GNU General Public License v3.0",
        "tsv_checksum": tsv_checksum,
        "interpretation_algorithms": [
            {
                "name": "drug_groups",
                "groups": {
                    "Nucleoside Reverse Transcriptase Inhibitor": [
                        "lamivudine",
                        "adefovir",
                        "tenofovir",
                        "entecavir",
                        "telbivudine",
                    ],
                },
            },
            {
                "name": "drug_interpretation",
                "method": "by_phenotype",
                "thresholds": {
                    "resistant": 1,
                },
            },
            {
                "name": "drug_alias",
                "groups": {
                    "lamivudine": "LAM",
                    "adefovir": "ADV",
                    "tenofovir": "TDF",
                    "entecavir": "ETV",
                    "telbivudine": "LDT",
                },
            },
        ],
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-url",
        default=SOURCE_URL,
        help="URL to the source DrugResistanceAlessaReformatted.txt file",
    )
    parser.add_argument(
        "--source-date",
        default=None,
        help=(
            "Date of the source data (YYYY-MM-DD). "
            "Defaults to the last GitHub commit date of the source file."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: <script_dir>/../output)",
    )
    args = parser.parse_args()

    # Determine output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path(__file__).resolve().parent.parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Source date: use explicit --source-date, or fetch from GitHub API
    if args.source_date:
        source_date = args.source_date
    else:
        source_date = fetch_source_commit_date(args.source_url)
    eprint(f"Source date: {source_date}")

    # Fetch and parse source
    eprint(f"Downloading {args.source_url} …")
    source_content = fetch_source(args.source_url)
    source_rows = parse_source(source_content)
    eprint(f"Parsed {len(source_rows)} source rows")

    # Sanity check: RT_OFFSET verification
    assert rt_to_p_position(204) == 539, f"RT_OFFSET sanity check failed: rt204 -> {rt_to_p_position(204)}"
    assert get_ref_aa(539) == "M", f"RT_OFFSET sanity check failed: P539 ref = {get_ref_aa(539)}, expected M"

    # Convert
    rules_rows, formula_rows, non_migrated = convert(source_rows)

    # Write rules.tsv
    rules_content = tsv_from_rows(rules_rows, RULES_COLUMNS)
    rules_path = output_dir / "rules.tsv"
    rules_path.write_text(rules_content, encoding="utf-8")
    eprint(f"Written {rules_path} ({len(rules_rows)} rows).")

    # Write formula-rules.tsv
    formula_path = output_dir / "formula-rules.tsv"
    if formula_rows:
        formula_content = tsv_from_rows(formula_rows, FORMULA_COLUMNS)
        formula_path.write_text(formula_content, encoding="utf-8")
        eprint(f"Written {formula_path} ({len(formula_rows)} rows).")
    else:
        formula_path.unlink(missing_ok=True)

    # Build and write metadata.json
    metadata = build_metadata(output_dir, source_date)
    metadata_path = output_dir / "metadata.json"
    with metadata_path.open("w", encoding="utf-8") as fh:
        json.dump(metadata, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    eprint(f"Written {metadata_path}.")

    # Write non-migrated-rules.txt
    nm_content = non_migrated_text(non_migrated)
    nm_path = output_dir / "non-migrated-rules.txt"
    nm_path.write_text(nm_content, encoding="utf-8")
    eprint(f"Written {nm_path} ({len(non_migrated)} aggregated entries).")

    # Validate required columns
    for row in rules_rows:
        required_cols = ["feature", "reference_identifier", "position", "reference", "mutation"]
        for col in required_cols:
            if not row.get(col):
                eprint(f"WARNING: Empty required column '{col}' in row: {row}")

    # Validate formula consistency
    if formula_rows:
        member_ids_in_rules = set()
        for row in rules_rows:
            mid = row.get("member_id", "")
            if mid:
                member_ids_in_rules.add(mid)
        for row in formula_rows:
            gid = row.get("group_id", "")
            expr = row.get("expression", "")
            for token in expr.replace("(", "").replace(")", "").split():
                if token not in ("AND", "OR", "NOT", "XOR") and token not in member_ids_in_rules:
                    eprint(f"WARNING: Formula references unknown member_id '{token}' in group {gid}")

    eprint("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

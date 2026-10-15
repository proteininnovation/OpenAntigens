from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache


_MIN_BOUNDARY_MAPPING_COVERAGE = 0.70


@dataclass(slots=True)
class AlignmentResult:
    aligned_query: str
    aligned_subject: str
    matches: int
    aligned_positions: int
    query_covered: int

    @property
    def identity(self) -> float:
        if self.aligned_positions == 0:
            return 0.0
        return 100.0 * self.matches / self.aligned_positions

    @property
    def coverage(self) -> float:
        if self.query_covered == 0:
            return 0.0
        return 100.0 * self.aligned_positions / self.query_covered


@dataclass(slots=True)
class PositionSubsetIdentity:
    matches: int
    aligned_positions: int
    requested_positions: int

    @property
    def identity(self) -> float | None:
        if self.aligned_positions == 0:
            return None
        return 100.0 * self.matches / self.aligned_positions


# BLOSUM62 with affine gaps 11/1 and FREE END GAPS.
#
# F33: the previous scheme scored match +1, **mismatch 0**, gap -1 in all three
# backends. Mismatches being free lets a short query cherry-pick matches out of
# a long subject: 20 seeded random pairs of 83 vs 550 aa averaged 75.6%
# identity at 100% coverage, and 200 vs 1,866 aa averaged 86.9%. That is why a
# wrong-gene macaque row read 79.7%, and why the mouse family matrix put GPR63
# at 79.31% to LGR4 (87 aa vs 520 aa; honestly ~38% over the row).
#
# mismatch = -1 is not enough (83 vs 550 still scores 77%), and BLOSUM62 with
# ordinary end gaps still opens ~15 gaps and reaches 38-44%. Free end gaps are
# the part that removes the artefact: the same pairs then score 14% and 11%.
#
# One backend, not three. The parasail and pure-python backends disagreed with
# Biopython by up to 0.8 points on real rows, which makes any pinned number
# fragile, and a silent fallback between scoring schemes is exactly how this
# defect survived. Biopython is a hard dependency (pyproject.toml:16).
_GAP_OPEN = -11.0
_GAP_EXTEND = -1.0


@lru_cache(maxsize=1)
def _aligner():  # pragma: no cover - exercised through global_align
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = _GAP_OPEN
    aligner.extend_gap_score = _GAP_EXTEND
    # Free end gaps: a length difference at either terminus costs nothing, so
    # the aligner stops paying for alignment it cannot justify.
    aligner.end_insertion_score = 0.0
    aligner.end_deletion_score = 0.0
    return aligner


@lru_cache(maxsize=1)
def _scoring_alphabet() -> frozenset[str]:
    from Bio.Align import substitution_matrices

    return frozenset(substitution_matrices.load("BLOSUM62").alphabet)


def _scoring_sequence(sequence: str) -> str:
    """Map residues BLOSUM62 does not define onto its wildcard.

    BLOSUM62's alphabet is ARNDCQEGHILKMFPSTWYVBZX*; selenocysteine (U),
    pyrrolysine (O), the ambiguity code J and any placeholder letter in a test
    fixture are outside it. X is the matrix's own "unknown residue" and scores
    -1 against everything, which is the honest answer. The ORIGINAL characters
    are kept in the aligned strings and in the identity count, so a placeholder
    still matches itself.
    """
    alphabet = _scoring_alphabet()
    return "".join(residue if residue in alphabet else "X" for residue in sequence)


def global_align(query: str, subject: str) -> AlignmentResult:
    return _global_align_cached(query, subject)


@lru_cache(maxsize=1024)
def _global_align_cached(query: str, subject: str) -> AlignmentResult:
    if not query or not subject:
        return _alignment_from_strings(
            query + "-" * len(subject),
            "-" * len(query) + subject,
            query_length=len(query),
        )

    alignment = _aligner().align(_scoring_sequence(query), _scoring_sequence(subject))[0]
    coordinates = alignment.coordinates

    aligned_query_parts: list[str] = []
    aligned_subject_parts: list[str] = []
    for index in range(coordinates.shape[1] - 1):
        query_start = int(coordinates[0, index])
        query_end = int(coordinates[0, index + 1])
        subject_start = int(coordinates[1, index])
        subject_end = int(coordinates[1, index + 1])
        query_span = query_end - query_start
        subject_span = subject_end - subject_start
        if query_span and subject_span:
            if query_span != subject_span:
                raise ValueError("Biopython returned an unsupported unequal aligned block")
            aligned_query_parts.append(query[query_start:query_end])
            aligned_subject_parts.append(subject[subject_start:subject_end])
        elif query_span:
            aligned_query_parts.append(query[query_start:query_end])
            aligned_subject_parts.append("-" * query_span)
        elif subject_span:
            aligned_query_parts.append("-" * subject_span)
            aligned_subject_parts.append(subject[subject_start:subject_end])

    aligned_query = "".join(aligned_query_parts)
    aligned_subject = "".join(aligned_subject_parts)
    if len(aligned_query) != len(aligned_subject):
        raise ValueError("Biopython returned an invalid traceback alignment")
    return _alignment_from_strings(aligned_query, aligned_subject, query_length=len(query))


def _alignment_from_strings(query_text: str, subject_text: str, *, query_length: int) -> AlignmentResult:
    matches = sum(
        1
        for q, s in zip(query_text, subject_text, strict=True)
        if q != "-" and s != "-" and q == s
    )
    aligned_positions = sum(
        1 for q, s in zip(query_text, subject_text, strict=True) if q != "-" and s != "-"
    )
    return AlignmentResult(
        aligned_query=query_text,
        aligned_subject=subject_text,
        matches=matches,
        aligned_positions=aligned_positions,
        query_covered=query_length,
    )


# An identity over aligned positions is only a statement about the sequences if
# the aligned stretch is a real part of them.
#
# Free end gaps (PR-4) made this bite: the optimum is free to trim to a tiny
# high-scoring core, so matches/aligned_positions goes *up* for exactly the
# wrong-gene pairs the metric fix was about. Measured on the macaque ortholog
# table, 126 of 4,593 rows aligned fewer than 10 residues and OR8J3 against
# guanine deaminase would have published "100.00% identity" off one residue.
#
# The rule is proportional, with no absolute floor: a genuinely short protein
# that aligns end to end has an honest identity, and an absolute floor would
# blank it.
_MIN_REPORTABLE_ALIGNED_FRACTION = 0.30


def reportable_identity(alignment: AlignmentResult, subject_length: int) -> float | None:
    """Identity, or ``None`` where the alignment is too thin to mean anything.

    For consumers that publish the number on its own, with no coverage beside
    it to qualify it. A consumer that shows coverage next to identity can use
    ``alignment.identity`` directly.
    """
    shorter = min(alignment.query_covered, subject_length)
    if shorter <= 0 or alignment.aligned_positions == 0:
        return None
    if alignment.aligned_positions / shorter < _MIN_REPORTABLE_ALIGNED_FRACTION:
        return None
    return alignment.identity


def find_furin_sites(sequence: str, pattern: str) -> list[tuple[int, int, str]]:
    matches: list[tuple[int, int, str]] = []
    expression = re.compile(pattern)
    position = 0
    while (match := expression.search(sequence, position)) is not None:
        if match.end() == match.start():
            raise ValueError("Motif patterns must not match an empty sequence")
        position = match.start() + 1
        start = match.start() + 1
        end = match.end()
        matches.append((start, end, match.group(0)))
    return matches


def slice_sequence(sequence: str, start: int, end: int) -> str:
    return sequence[start - 1 : end]


def identity_for_query_positions(
    alignment: AlignmentResult,
    *,
    query_region_start: int,
    query_positions: set[int],
) -> PositionSubsetIdentity:
    query_position = query_region_start - 1
    matches = 0
    aligned_positions = 0
    requested = 0
    for query_residue, subject_residue in zip(
        alignment.aligned_query, alignment.aligned_subject, strict=True
    ):
        if query_residue != "-":
            query_position += 1
        if query_residue == "-" or query_position not in query_positions:
            continue
        requested += 1
        if subject_residue == "-":
            continue
        aligned_positions += 1
        if query_residue == subject_residue:
            matches += 1
    return PositionSubsetIdentity(
        matches=matches,
        aligned_positions=aligned_positions,
        requested_positions=requested,
    )


def map_query_region_to_subject(
    alignment: AlignmentResult,
    *,
    query_start: int,
    query_end: int,
    subject_sequence: str,
    requested_length: int | None = None,
) -> tuple[int | None, int | None, str | None, list[str]]:
    """Project a query range onto the subject through ``alignment``.

    ``requested_length`` is how many residues were *asked for*. Without it the
    coverage denominator is only the residues the alignment happens to contain,
    so projecting a full-length construct through an ectodomain-only alignment
    passed the 70% gate on the ectodomain's own length (OA-012: 422 human / 365
    mouse reports, 2,834 constructs with under 70% of their span inside the
    ectodomain).
    """
    query_position = 0
    subject_position = 0
    mapped_positions: list[int] = []
    gap_count = 0
    query_residue_count = 0

    for query_residue, subject_residue in zip(
        alignment.aligned_query, alignment.aligned_subject, strict=True
    ):
        if query_residue != "-":
            query_position += 1
        if subject_residue != "-":
            subject_position += 1
        if query_residue == "-" or not (query_start <= query_position <= query_end):
            continue
        query_residue_count += 1
        if subject_residue == "-":
            gap_count += 1
            continue
        mapped_positions.append(subject_position)

    if not mapped_positions:
        return None, None, None, ["No aligned residues mapped to the homolog ectodomain."]

    denominator = requested_length if requested_length else query_residue_count
    coverage = len(mapped_positions) / denominator if denominator else 0.0
    if coverage < _MIN_BOUNDARY_MAPPING_COVERAGE:
        return None, None, None, [
            "Boundary projection rejected: "
            f"{len(mapped_positions)}/{denominator} query residues aligned to the homolog "
            f"({coverage:.1%} < {_MIN_BOUNDARY_MAPPING_COVERAGE:.1%} minimum)."
        ]

    start = min(mapped_positions)
    end = max(mapped_positions)
    notes: list[str] = []
    if gap_count:
        notes.append(f"{gap_count} construct residues align to gaps in the homolog ectodomain.")
    return start, end, slice_sequence(subject_sequence, start, end), notes

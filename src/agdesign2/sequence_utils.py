from __future__ import annotations

import os
import re
from dataclasses import dataclass
from functools import lru_cache


_ALIGNMENT_BACKEND_ENV = "AGDESIGN2_ALIGNMENT_BACKEND"
_PARASAIL_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZXJUO*arn dcqeghilkmfpstwyvbzxjuo".replace(" ", "")


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


def global_align(query: str, subject: str) -> AlignmentResult:
    return _global_align_cached(query, subject, _alignment_backend())


def _alignment_backend() -> str:
    backend = os.environ.get(_ALIGNMENT_BACKEND_ENV, "auto").strip().lower()
    return backend if backend in {"auto", "parasail", "python"} else "auto"


@lru_cache(maxsize=1024)
def _global_align_cached(query: str, subject: str, backend: str) -> AlignmentResult:
    if backend in {"auto", "parasail"}:
        try:
            return _parasail_global_align(query, subject)
        except Exception:
            if backend == "parasail":
                raise
    if backend == "auto":
        try:
            return _biopython_global_align(query, subject)
        except Exception:
            pass
    return _python_global_align(query, subject)


def _parasail_global_align(query: str, subject: str) -> AlignmentResult:
    if not query or not subject:
        return _python_global_align(query, subject)
    if set(query) - set(_PARASAIL_ALPHABET) or set(subject) - set(_PARASAIL_ALPHABET):
        return _python_global_align(query, subject)

    import parasail  # type: ignore[import-not-found]

    matrix = parasail.matrix_create(_PARASAIL_ALPHABET, 1, 0)
    result = parasail.nw_trace_striped_16(query, subject, 1, 1, matrix)
    traceback = result.traceback
    aligned_query = str(traceback.query)
    aligned_subject = str(traceback.ref)
    if not aligned_query or not aligned_subject or len(aligned_query) != len(aligned_subject):
        raise ValueError("parasail returned an invalid traceback alignment")
    return _alignment_from_strings(aligned_query, aligned_subject, query_length=len(query))


def _biopython_global_align(query: str, subject: str) -> AlignmentResult:
    if not query or not subject:
        return _python_global_align(query, subject)

    from Bio.Align import PairwiseAligner  # type: ignore[import-not-found]

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    alignment = aligner.align(query, subject)[0]
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


def _python_global_align(query: str, subject: str) -> AlignmentResult:
    rows = len(query) + 1
    cols = len(subject) + 1
    score = [[0] * cols for _ in range(rows)]
    trace = [[""] * cols for _ in range(rows)]

    for i in range(1, rows):
        score[i][0] = -i
        trace[i][0] = "U"
    for j in range(1, cols):
        score[0][j] = -j
        trace[0][j] = "L"

    for i in range(1, rows):
        for j in range(1, cols):
            match_score = 1 if query[i - 1] == subject[j - 1] else 0
            diag = score[i - 1][j - 1] + match_score
            up = score[i - 1][j] - 1
            left = score[i][j - 1] - 1
            best = max(diag, up, left)
            score[i][j] = best
            trace[i][j] = "D" if best == diag else "U" if best == up else "L"

    aligned_query: list[str] = []
    aligned_subject: list[str] = []
    i = len(query)
    j = len(subject)
    while i > 0 or j > 0:
        direction = trace[i][j] if i >= 0 and j >= 0 else ""
        if direction == "D":
            aligned_query.append(query[i - 1])
            aligned_subject.append(subject[j - 1])
            i -= 1
            j -= 1
        elif direction == "U":
            aligned_query.append(query[i - 1])
            aligned_subject.append("-")
            i -= 1
        else:
            aligned_query.append("-")
            aligned_subject.append(subject[j - 1])
            j -= 1

    aligned_query.reverse()
    aligned_subject.reverse()
    return _alignment_from_strings("".join(aligned_query), "".join(aligned_subject), query_length=len(query))


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
) -> tuple[int | None, int | None, str | None, list[str]]:
    query_position = 0
    subject_position = 0
    mapped_positions: list[int] = []
    gap_count = 0

    for query_residue, subject_residue in zip(
        alignment.aligned_query, alignment.aligned_subject, strict=True
    ):
        if query_residue != "-":
            query_position += 1
        if subject_residue != "-":
            subject_position += 1
        if query_residue == "-" or not (query_start <= query_position <= query_end):
            continue
        if subject_residue == "-":
            gap_count += 1
            continue
        mapped_positions.append(subject_position)

    if not mapped_positions:
        return None, None, None, ["No aligned residues mapped to the homolog ectodomain."]

    start = min(mapped_positions)
    end = max(mapped_positions)
    notes: list[str] = []
    if gap_count:
        notes.append(f"{gap_count} construct residues align to gaps in the homolog ectodomain.")
    return start, end, slice_sequence(subject_sequence, start, end), notes

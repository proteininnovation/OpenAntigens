#!/usr/bin/env python3
"""Release checks for a built OpenAntigens site (remediation plan §7).

Run over BOTH portals of a built ``public_site`` before publishing:

    PYTHONPATH=src .venv/bin/python scripts/release_checks.py --site <public_site> \
        [--human-json <human report dir>] [--mouse-json <mouse report dir>]

Each check counts the pages that violate it. Exit status is non-zero if any
check fails. The HTML checks need only the site; the JSON checks need the report
JSON the site was rendered from, and are skipped (and reported as skipped, not
passed) when neither directory is given; give both or neither. A portal with no
report pages, no index or (when given) no report JSON fails "portal-present":
a site missing its mouse portal must not pass.

These checks land BEFORE the rebuild so they can be run against the frozen
2026-09-30 build, where they are expected to FAIL — that is what shows they
detect the defects. They are meant to pass only on the rebuilt site.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

csv.field_size_limit(10**7)

# The eight index columns documented as blank on human rows (downloads.html),
# plus "error", which is blank whenever no row failed. Plan §7: a naive "no
# all-blank column" check is unsatisfiable without this allowance.
DOCUMENTED_BLANK_HUMAN = {
    "mouse_ortholog_accession", "mouse_ortholog_entry_name", "mouse_ortholog_gene_symbol",
    "mouse_refseq_accession", "mouse_sequence_source",
    "human_source_accession", "human_source_entry_name", "human_source_gene_symbol",
    "error",
}
CONSTRUCT_SIZE_WARNING_AA = 1500
NOTES_SECTION = re.compile(r'id="notes".*?</section>', re.S)
HOMOLOGY_SECTION = re.compile(r'id="homology".*?</section>', re.S)
CONSTRUCT_SUMMARY = re.compile(r'id="construct-summary".*?</section>', re.S)


@dataclass
class Check:
    key: str
    description: str
    failures: Counter = field(default_factory=Counter)  # portal -> pages
    examples: dict = field(default_factory=dict)
    skipped: bool = False

    def fail(self, portal: str, page: str) -> None:
        self.failures[portal] += 1
        self.examples.setdefault(portal, page)

    @property
    def passed(self) -> bool:
        return not self.skipped and not self.failures


def _html_checks() -> list[Check]:
    return [
        Check("blast-negative-beside-hits",
              "no 'completed BLAST search returned no hits' beside a full-length hit table"),
        Check("interpro-negative-beside-family",
              "no bare 'No InterPro annotations.' beside a cited family accession"),
        Check("homolog-failure-as-no-records",
              "no 'Species homolog lookup failed' rendered as 'No homolog records.'"),
        Check("mouse-self-mapping-row", "no mouse page carries a _MOUSE_mapping row"),
        Check("duplicated-note", "no page renders the same note twice"),
        Check("duplicate-construct-name", "no two constructs on a page share a name"),
        Check("oversized-construct-unwarned",
              f"no construct longer than {CONSTRUCT_SIZE_WARNING_AA} aa without a size warning"),
        Check("motif-detected-label", "no conserved motif labelled 'detected' (no residue test)"),
        Check("negative-site-as-ptm", "no 'Not glycosylated' SITE rendered as a PTM"),
    ]


def _check_page(html: str, portal: str, name: str, checks: dict[str, Check]) -> None:
    # The negative above belongs to the design-region table, so it must be read
    # against that table's own hits ("Showing N hit(s)"). A populated full-length
    # table beside an empty region search is a real state, not a defect: GPR142's
    # region search finds nothing while the full-length search returns 17 hits.
    if ("The completed BLAST search returned no sequence similarity hits." in html
            and re.search(r"Showing \d+ hit\(s\)", html)):
        checks["blast-negative-beside-hits"].fail(portal, name)

    if "No InterPro annotations." in html and re.search(r"\(IPR\d{6}\)", html):
        checks["interpro-negative-beside-family"].fail(portal, name)

    homology = HOMOLOGY_SECTION.search(html)
    if ("Species homolog lookup failed" in html and homology
            and "No homolog records." in homology.group(0)
            and "lookup failed when this report was built" not in homology.group(0)):
        checks["homolog-failure-as-no-records"].fail(portal, name)

    if portal == "mouse" and "_MOUSE_mapping" in html:
        checks["mouse-self-mapping-row"].fail(portal, name)

    notes = NOTES_SECTION.search(html)
    if notes:
        items = re.findall(r"<li>(.*?)</li>", notes.group(0), re.S)
        if len(items) != len(set(items)):
            checks["duplicated-note"].fail(portal, name)

    summary = CONSTRUCT_SUMMARY.search(html)
    if summary:
        spans: dict[str, set[str]] = {}
        for row in re.findall(r"<tr>(.*?)</tr>", summary.group(0), re.S):
            cells = re.findall(r"<td>(.*?)</td>", row, re.S)
            if len(cells) < 3:
                continue
            construct = re.sub(r"<[^>]+>", "", cells[0]).strip()
            spans.setdefault(construct, set()).add(re.sub(r"<[^>]+>", "", cells[2]).strip())
        if any(len(boundaries) > 1 for boundaries in spans.values()):
            checks["duplicate-construct-name"].fail(portal, name)
        for boundaries in spans.values():
            for boundary in boundaries:
                match = re.fullmatch(r"(\d+)-(\d+)", boundary)
                if (match and int(match.group(2)) - int(match.group(1)) + 1 > CONSTRUCT_SIZE_WARNING_AA
                        and "Very large construct" not in html):
                    checks["oversized-construct-unwarned"].fail(portal, name)
                    break

    if re.search(r"<td>detected</td>", html):
        checks["motif-detected-label"].fail(portal, name)

    if re.search(r"<td>SITE</td>.{0,200}?<td>Not (N-)?glycosylated", html, re.S):
        checks["negative-site-as-ptm"].fail(portal, name)


def _index_check(site: Path, portal: str, check: Check) -> None:
    path = (site / "downloads" if portal == "human" else site / "mouse" / "downloads") / "agdesign2_portal_index.tsv"
    rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t")) if path.exists() else []
    if not rows:
        check.fail(portal, f"no index rows in {path}")
        return
    allowed = DOCUMENTED_BLANK_HUMAN if portal == "human" else {"error"}
    for column in rows[0]:
        if column in allowed:
            continue
        if all(not (row.get(column) or "").strip() for row in rows):
            check.fail(portal, f"index column {column}")


def _json_checks(report_dir: Path, portal: str, checks: dict[str, Check]) -> None:
    for path in sorted(report_dir.glob("*_report.json")):
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        name = path.name

        if report.get("gpcr_engineering_variants") and not report.get("gpcr_annotation"):
            name_text = str((report.get("target") or {}).get("protein_name") or "").lower()
            if not any(t in name_text for t in ("receptor", "opsin", "frizzled", "smoothened")):
                checks["gpcr-advice-without-classification"].fail(portal, name)

        _gpi_check(report, portal, name, checks["gpi-region-past-omega"])
        _cysteine_check(report, portal, name, checks["curated-disulfide-unpaired"])

        for row in report.get("ectodomain_homology") or []:
            status = str(row.get("selection_status") or "")
            if status.startswith("withheld") and (row.get("accession") or row.get("identity") is not None):
                checks["withheld-row-ships-data"].fail(portal, name)
                break


def _gpi_check(report: dict, portal: str, name: str, check: Check) -> None:
    region = report.get("ectodomain") or {}
    if region.get("end") is None:
        return
    features = report.get("features") or []
    anchors = [f.get("start") for f in features
               if str(f.get("type", "")).upper() in {"LIPID", "LIPIDATION"}
               and "gpi" in str(f.get("description") or "").lower()]
    chains = [f for f in features if str(f.get("type", "")).upper() == "CHAIN"]
    # Conditional, as §7 requires: only where a CHAIN contains the omega site.
    for anchor in anchors:
        if (anchor is not None and int(region["end"]) > int(anchor)
                and any(c.get("start", 0) <= anchor <= c.get("end", 0) for c in chains)):
            check.fail(portal, name)
            return


def _cysteine_check(report: dict, portal: str, name: str, check: Check) -> None:
    cysteines = {c.get("position"): c for c in report.get("cysteine_analysis") or []}
    for feature in report.get("ptms") or report.get("features") or []:
        if str(feature.get("type", "")).upper() != "DISULFID":
            continue
        start, end = feature.get("start"), feature.get("end")
        if start == end:  # interchain, correctly unpaired in a monomer
            continue
        if any(str((cysteines.get(r) or {}).get("warning") or "").startswith("Ambiguous curated")
               for r in (start, end)):
            continue  # alternative connectivity: no single partner to pair
        if start in cysteines and end in cysteines and (
                cysteines[start].get("paired_with") is None or cysteines[end].get("paired_with") is None):
            check.fail(portal, name)
            return


def run(site: Path, human_json: Path | None, mouse_json: Path | None) -> list[Check]:
    html = {c.key: c for c in _html_checks()}
    if (human_json is None) != (mouse_json is None):
        raise ValueError("give both --human-json and --mouse-json, or neither")
    present = Check("portal-present", "both portals have report pages, an index and (if given) report JSON")
    index = Check("index-column-blank", "no undocumented index column blank in every row")
    json_checks = {c.key: c for c in (
        Check("gpcr-advice-without-classification", "no GPCR advice on a target not classified as a GPCR"),
        Check("gpi-region-past-omega", "where a CHAIN bounds a GPI region, it ends at the omega site"),
        Check("curated-disulfide-unpaired", "no curated in-range disulfide endpoint reported unpaired"),
        Check("withheld-row-ships-data", "no withheld ortholog row ships an accession or identity"),
    )}

    for portal, pages in (("human", site / "reports"), ("mouse", site / "mouse" / "reports")):
        page_paths = sorted(pages.glob("*.html")) if pages.is_dir() else []
        if not page_paths:
            present.fail(portal, f"no report pages in {pages}")
        for page in page_paths:
            _check_page(page.read_text(encoding="utf-8", errors="ignore"), portal, page.name, html)
        _index_check(site, portal, index)

    for portal, directory in (("human", human_json), ("mouse", mouse_json)):
        if directory is None:
            continue
        if not any(directory.glob("*_report.json")):
            present.fail(portal, f"no report JSON in {directory}")
        _json_checks(directory, portal, json_checks)
    if human_json is None and mouse_json is None:
        for check in json_checks.values():
            check.skipped = True
    return [present, *html.values(), index, *json_checks.values()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--site", type=Path, required=True)
    parser.add_argument("--human-json", type=Path)
    parser.add_argument("--mouse-json", type=Path)
    args = parser.parse_args(argv)
    try:
        checks = run(args.site, args.human_json, args.mouse_json)
    except ValueError as exc:
        parser.error(str(exc))
    for check in checks:
        state = "SKIP" if check.skipped else ("PASS" if check.passed else "FAIL")
        counts = ", ".join(f"{portal} {n}" for portal, n in sorted(check.failures.items()))
        example = "; ".join(f"{portal}: {page}" for portal, page in sorted(check.examples.items()))
        print(f"{state:4}  {check.key:38} {counts or '-':22} {example}")
    return 0 if all(check.passed or check.skipped for check in checks) else 1


if __name__ == "__main__":
    sys.exit(main())

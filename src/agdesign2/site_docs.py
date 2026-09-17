"""Shared preprint citation and agent documentation for every portal backend."""
from __future__ import annotations

from html import escape
import json
from pathlib import Path

DATA = Path(__file__).with_name("data")
PAPER = json.loads((DATA / "openantigens_citation.json").read_text(encoding="utf-8"))
DOI_URL = "https://doi.org/" + PAPER["doi"]


def citation_text() -> str:
    authors = ", ".join(
        a["family"] + " " + "".join(part[0] for part in a["given"].split())
        for a in PAPER["authors"]
    )
    return f'{authors}. {PAPER["title"]}. {PAPER["journal"]} [preprint]. {PAPER["date"][:4]}. {DOI_URL}'


def citation_html() -> str:
    return escape(citation_text()).replace(DOI_URL, f'<a class="inline-link" href="{DOI_URL}">{DOI_URL}</a>')


def citation_short_html() -> str:
    label = f'{PAPER["authors"][0]["family"]} et al., {PAPER["journal"]} ({PAPER["date"][:4]})'
    return f'<a class="inline-link" href="{DOI_URL}">{escape(label)}</a>'


def citation_metadata() -> dict:
    return {**PAPER, "url": DOI_URL, "text": citation_text()}


def citation_downloads() -> dict[str, str]:
    authors = " and ".join(f'{a["family"]}, {a["given"]}' for a in PAPER["authors"])
    bib = (
        "@article{Teixeira2026OpenAntigens,\n"
        f"  author = {{{authors}}},\n"
        f'  title = {{{{{PAPER["title"]}}}}},\n'
        f'  journal = {{{PAPER["journal"]}}},\n'
        f'  year = {{{PAPER["date"][:4]}}},\n'
        f'  doi = {{{PAPER["doi"]}}},\n'
        f'  url = {{{DOI_URL}}},\n'
        f'  note = {{Preprint, version {PAPER["version"]}; posted {PAPER["date"]}}}\n'
        "}\n"
    )
    ris = "\n".join([
        "TY  - JOUR",
        *(f'AU  - {a["family"]}, {a["given"]}' for a in PAPER["authors"]),
        f'TI  - {PAPER["title"]}', f'JO  - {PAPER["journal"]}',
        f'PY  - {PAPER["date"][:4]}', f'DA  - {PAPER["date"].replace("-", "/")}',
        f'DO  - {PAPER["doi"]}', f'UR  - {DOI_URL}',
        f'N1  - Preprint, version {PAPER["version"]}; not peer reviewed.', "ER  -", "",
    ])
    return {"openantigens.bib": bib, "openantigens.ris": ris}


def citation_cff() -> str:
    def authors_yaml(authors: list[dict], indent: int) -> str:
        pad = " " * indent
        return "\n".join(
            f'{pad}- family-names: {json.dumps(a["family"])}\n'
            f'{pad}  given-names: {json.dumps(a["given"])}\n'
            f'{pad}  orcid: {json.dumps(a["orcid"])}'
            for a in authors
        )

    return f"""cff-version: 1.2.0
type: software
title: OpenAntigens
message: "If you use OpenAntigens in research, please cite the preprint below and record the release used."
authors:
{authors_yaml(PAPER['authors'][:1], 2)}
url: https://openantigens.org/
repository-code: https://github.com/proteininnovation/OpenAntigens
license: Apache-2.0
preferred-citation:
  type: article
  title: {json.dumps(PAPER['title'])}
  authors:
{authors_yaml(PAPER['authors'], 4)}
  journal: {PAPER['journal']}
  year: {PAPER['date'][:4]}
  date-released: {PAPER['date']}
  doi: {PAPER['doi']}
  url: {DOI_URL}
  notes: "Preprint, version {PAPER['version']}; not peer reviewed."
"""


def agent_guide_markdown() -> str:
    return (DATA / "agent_guide.md").read_text(encoding="utf-8").replace("{{citation}}", citation_text())


def agent_guide_html() -> str:
    import mistune

    # The page shell supplies its own title; retain the full guide below it.
    return mistune.create_markdown(escape=True)(agent_guide_markdown().removeprefix("# OpenAntigens\n"))


def site_document_files() -> dict[str, str]:
    return {
        "llms.txt": agent_guide_markdown(),
        "agent-guide.js": (DATA / "agent-guide.js").read_text(encoding="utf-8"),
        **{f"downloads/{name}": content for name, content in citation_downloads().items()},
    }


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    (root / "llms.txt").write_text(agent_guide_markdown(), encoding="utf-8")
    (root / "CITATION.cff").write_text(citation_cff(), encoding="utf-8")

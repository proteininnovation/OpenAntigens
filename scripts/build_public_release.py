#!/usr/bin/env python3
"""Build a public OpenAntigens static-site release bundle.

The release bundle is intentionally limited to the rendered portal directory
plus public license/manifest files. Local caches, BLAST databases, refresh
state, logs, and intermediate analysis outputs are not copied.
"""

from __future__ import annotations

import argparse
import hashlib
import gzip
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import json
import shutil
import sys
import tomllib
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from agdesign2.release_validation import validate_release_data, validate_structure_coverage


REQUIRED_PORTAL_FILES = (
    ".htaccess",
    "apple-touch-icon.png",
    "builder.html",
    "agent-guide.html",
    "agent-guide.js",
    "llms.txt",
    "downloads/openantigens.bib",
    "downloads/openantigens.ris",
    "calculator.html",
    "constructs.html",
    "downloads.html",
    "downloads/agdesign2_portal_index.json",
    "downloads/agdesign2_portal_index.tsv",
    "downloads/download_manifest.json",
    "favicon-32.png",
    "favicon.ico",
    "help.html",
    "icon-192.png",
    "icon-512.png",
    "index.html",
    "ipi-logo-dark-800.png",
    "ipi-logo-light-800.png",
    "methods.html",
    "portal-index-data.js",
    "portal-index.js",
    "portal.css",
    "portal_metadata.json",
    "privacy.html",
    "site.webmanifest",
    "terms.html",
)

REQUIRED_PORTAL_DIRECTORIES = (
    "downloads",
    "reports",
)

EXCLUDED_NAMES = {
    ".DS_Store",
    "Thumbs.db",
    "__pycache__",
}

PUBLIC_AF3_MARKERS = (
    "AlphaFold 3",
    "AlphaFold Server",
    "alphafold3",
)

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--portal-dir",
        default="outputs/surfy_batch/portal",
        help="Rendered static portal directory to publish.",
    )
    parser.add_argument(
        "--output-root",
        default="releases",
        help="Directory where versioned release folders are created.",
    )
    parser.add_argument(
        "--release-name",
        default=None,
        help="Release folder name. Defaults to openantigen-YYYY-MM-DD.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing release directory with the same name.",
    )
    parser.add_argument(
        "--checksums",
        action="store_true",
        help="Compute SHA-256 checksums for every published file. Slower for large releases.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and print the planned output path without copying files.",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    portal_dir = (repo_root / args.portal_dir).resolve()
    output_root = (repo_root / args.output_root).resolve()
    release_name = args.release_name or f"openantigen-{datetime.now(timezone.utc).date().isoformat()}"

    try:
        release_dir = resolve_release_dir(output_root, release_name)
        validate_portal(portal_dir)
        read_portal_metadata(portal_dir)
    except (ValueError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.dry_run:
        print(f"portal_dir={portal_dir}")
        print(f"release_dir={release_dir}")
        print("dry_run=true")
        return 0

    if release_dir.exists():
        if not args.force:
            print(
                f"error: release directory already exists: {release_dir}\n"
                "Use --force to replace it or choose --release-name.",
                file=sys.stderr,
            )
            return 2
        shutil.rmtree(release_dir)

    output_root.mkdir(parents=True, exist_ok=True)
    copy_public_portal(portal_dir, release_dir)
    copy_public_license_files(repo_root, release_dir)

    manifest = build_manifest(
        repo_root=repo_root,
        release_dir=release_dir,
        checksums=args.checksums,
    )
    (release_dir / "release_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (release_dir / "README_RELEASE.md").write_text(render_release_readme(manifest), encoding="utf-8")

    print(f"release_dir={release_dir}")
    print(f"files={manifest['file_count']}")
    print(f"size_bytes={manifest['total_size_bytes']}")
    print(f"size_human={manifest['total_size_human']}")
    print("serve_root=release_dir")
    return 0


def _read_portal_text(path: Path) -> str:
    data = path.read_bytes()
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    return data.decode("utf-8")


def validate_portal(portal_dir: Path) -> None:
    if not portal_dir.exists():
        raise ValueError(f"portal directory does not exist: {portal_dir}")
    if not portal_dir.is_dir():
        raise ValueError(f"portal path is not a directory: {portal_dir}")
    symlinks = [path.relative_to(portal_dir) for path in portal_dir.rglob("*") if path.is_symlink()]
    if symlinks:
        raise ValueError(f"portal directory contains symbolic links: {', '.join(map(str, symlinks))}")
    missing = [
        name
        for name in REQUIRED_PORTAL_FILES
        if not (portal_dir / name).is_file() or (portal_dir / name).stat().st_size == 0
    ]
    if missing:
        raise ValueError(f"portal directory has missing or empty required files: {', '.join(missing)}")
    missing_directories = [
        name for name in REQUIRED_PORTAL_DIRECTORIES if not (portal_dir / name).is_dir()
    ]
    if missing_directories:
        raise ValueError(f"portal directory is missing required directories: {', '.join(missing_directories)}")
    validate_structure_coverage(portal_dir)
    validate_release_data(portal_dir)
    index_data_paths = [portal_dir / "portal-index-data.js"]
    mouse_index_data = portal_dir / "mouse" / "portal-index-data.js"
    if mouse_index_data.is_file():
        index_data_paths.append(mouse_index_data)
    missing_pubtator = [
        path.relative_to(portal_dir)
        for path in index_data_paths
        if "pubtator3/docsum" not in _read_portal_text(path)
    ]
    if missing_pubtator:
        raise ValueError(
            "portal index data contains no populated PubTator links: "
            + ", ".join(map(str, missing_pubtator))
        )
    structure_images = [
        path.relative_to(portal_dir)
        for path in portal_dir.rglob("*_structure.png")
        if path.is_file()
    ]
    if structure_images:
        raise ValueError(
            "portal contains excluded static PyMOL structure images: "
            + ", ".join(map(str, structure_images[:10]))
        )
    for path in iter_files(portal_dir):
        if path.suffix.lower() not in {".css", ".html", ".js", ".json", ".md", ".tsv", ".txt"}:
            continue
        try:
            content = _read_portal_text(path)
        except UnicodeDecodeError:
            continue
        marker = next((value for value in PUBLIC_AF3_MARKERS if value in content), None)
        if marker is not None:
            raise ValueError(
                f"portal contains excluded AF3 content ({marker!r}): {path.relative_to(portal_dir)}"
            )

    validate_frontend_assets(portal_dir)

def resolve_release_dir(output_root: Path, release_name: str) -> Path:
    output_root = output_root.resolve()
    release_dir = (output_root / release_name).resolve()
    if release_dir.parent != output_root:
        raise ValueError("release name must be one directory name")
    return release_dir


def validate_frontend_assets(portal_dir: Path) -> None:
    class AssetParser(HTMLParser):
        def handle_starttag(self, tag, attrs):
            values = dict(attrs)
            url = values.get("src") if tag == "script" else values.get("href") if tag == "link" and values.get("rel") == "stylesheet" else None
            if not url:
                return
            parsed = urlsplit(url)
            if parsed.scheme or parsed.netloc:
                return
            asset = (portal_dir / unquote(parsed.path).lstrip("/")) if parsed.path.startswith("/") else page.parent / unquote(parsed.path)
            if not asset.resolve().is_relative_to(portal_dir.resolve()) or not asset.is_file():
                raise ValueError(f"Missing local frontend asset in {page.relative_to(portal_dir)}: {url}")

    for page in portal_dir.rglob("*.html"):
        AssetParser().feed(_read_portal_text(page))


def copy_public_portal(portal_dir: Path, release_dir: Path) -> None:
    validate_portal(portal_dir)
    shutil.copytree(
        portal_dir,
        release_dir,
        ignore=ignore_public_release_files,
        copy_function=shutil.copy2,
        symlinks=True,
    )


def ignore_public_release_files(directory: str, names: list[str]) -> set[str]:
    ignored = set()
    for name in names:
        if name in EXCLUDED_NAMES:
            ignored.add(name)
        elif name.startswith(".") and name not in {".htaccess", ".well-known"}:
            ignored.add(name)
        elif name.endswith((".log", ".tmp", ".bak")):
            ignored.add(name)
    return ignored


def copy_public_license_files(repo_root: Path, release_dir: Path) -> None:
    for name in ("LICENSE", "DATA_LICENSE.md"):
        source = repo_root / name
        if source.is_file():
            shutil.copy2(source, release_dir / name)


def build_manifest(
    *,
    repo_root: Path,
    release_dir: Path,
    checksums: bool,
) -> dict[str, object]:
    files = list(iter_files(release_dir))
    total_size = sum(path.stat().st_size for path in files)
    extensions = Counter(path.suffix.lower() or "[no extension]" for path in files)
    report_count = count_files(release_dir / "reports", "*.html")
    structure_count = count_files(release_dir / "structures", "*")
    download_count = count_files(release_dir / "downloads", "*")

    manifest: dict[str, object] = {
        "database": "OpenAntigens",
        "release_created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "software_version": read_project_version(repo_root),
        "portal_metadata": read_portal_metadata(release_dir),
        "serving_note": "Serve this directory as the web root; index.html is at the release root.",
        "file_count": len(files),
        "total_size_bytes": total_size,
        "total_size_human": human_size(total_size),
        "html_report_count": report_count,
        "structure_file_count": structure_count,
        "download_file_count": download_count,
        "file_extensions": dict(sorted(extensions.items())),
        "licenses": {
            "software": "Apache-2.0",
            "openantigen_generated_annotations": "CC BY 4.0",
            "third_party_source_data": "Retains original source database terms; see DATA_LICENSE.md.",
        },
        "excluded_from_release": [
            ".agdesign2 caches and BLAST databases",
            "outputs/surfy_batch intermediate JSON/Markdown/log/state files outside the rendered portal",
            "hidden desktop metadata such as .DS_Store",
        ],
    }
    if checksums:
        manifest["sha256"] = {
            str(path.relative_to(release_dir)): sha256_file(path) for path in files
        }
    return manifest


def iter_files(root: Path) -> Iterable[Path]:
    for path in sorted(root.rglob("*")):
        if path.is_file():
            yield path


def count_files(root: Path, pattern: str) -> int:
    if not root.exists():
        return 0
    return sum(1 for path in root.glob(pattern) if path.is_file())


def read_project_version(repo_root: Path) -> str:
    pyproject = repo_root / "pyproject.toml"
    if not pyproject.is_file():
        return "unknown"
    with pyproject.open("rb") as handle:
        data = tomllib.load(handle)
    return str(data.get("project", {}).get("version", "unknown"))


def read_portal_metadata(release_dir: Path) -> dict[str, object] | None:
    path = release_dir / "portal_metadata.json"
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"portal metadata must be a JSON object: {path}")
    return payload


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def human_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{size} B"


def render_release_readme(manifest: dict[str, object]) -> str:
    citation = (manifest.get("portal_metadata") or {}).get("citation")
    citation_section = (
        "## Recommended scholarly citation\n\n" + citation["text"] + "\n\n"
        "Record the portal build date, software version, artifact URL, and access date alongside the paper citation. "
        "Citation files: downloads/openantigens.bib and downloads/openantigens.ris. "
        "Third-party source terms and data attribution requirements still apply.\n\n"
        if citation else ""
    )
    return f"""# OpenAntigens Public Release

This directory is a static OpenAntigens portal release.

Serve this directory as the web root. The entry point is:

```text
index.html
```

## Contents

- HTML portal pages and per-target reports
- copied portal assets, structures, and downloadable index files
- `release_manifest.json` with release metadata and file counts
- `LICENSE` and `DATA_LICENSE.md`

## Release Summary

- Created UTC: `{manifest['release_created_utc']}`
- Software version: `{manifest['software_version']}`
- Files: `{manifest['file_count']}`
- Size: `{manifest['total_size_human']}`
- HTML reports: `{manifest['html_report_count']}`
- Structure files: `{manifest['structure_file_count']}`
- Download files: `{manifest['download_file_count']}`

{citation_section}## Not Included

This public bundle intentionally excludes local caches, BLAST databases,
refresh state, logs, and intermediate analysis files. Rebuild those from the
source repository if a full computational refresh is needed.
"""


if __name__ == "__main__":
    raise SystemExit(main())

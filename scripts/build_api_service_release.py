#!/usr/bin/env python3
"""Build a local, validated Node 22 API/MCP deployment ZIP; never deploy it."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from agdesign2.public_api_release import build_api_service_release


def main(argv: list[str] | None = None) -> int:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--service-dir", type=Path, default=repo / "services/public-api")
    parser.add_argument("--bundle-dir", type=Path, required=True, help="Validated frozen API data bundle.")
    parser.add_argument("--public-site-dir", type=Path, required=True, help="Final static artifact with all retained downloads.")
    parser.add_argument("--output", type=Path, required=True, help="New .zip path outside the input directories.")
    parser.add_argument("--repo-root", type=Path, default=repo, help="Canonical source checkout containing public licenses.")
    parser.add_argument("--dry-run", action="store_true", help="Validate and measure the real ZIP; preserve requested output paths.")
    parser.add_argument("--max-zip-bytes", type=int, default=100_000_000, help="Local compressed budget; hard ceiling is strictly below 100 MB.")
    parser.add_argument("--max-expanded-bytes", type=int, default=500_000_000, help="Local expanded budget; verify the actual host separately.")
    parser.add_argument("--reserve-fraction", type=float, default=0.20)
    parser.add_argument("--dependency-audit-report", type=Path, help="Record a local npm audit JSON report hash; this does not assert a passed gate.")
    args = parser.parse_args(argv)
    try:
        result = build_api_service_release(
            args.service_dir, args.bundle_dir, args.public_site_dir, args.output,
            repo_root=args.repo_root, dry_run=args.dry_run,
            max_zip_bytes=args.max_zip_bytes, max_expanded_bytes=args.max_expanded_bytes,
            reserve_fraction=args.reserve_fraction,
            dependency_audit_report=args.dependency_audit_report,
        )
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    summary = {key: value for key, value in result.items() if key != "files"}
    summary["artifact_record"] = None if args.dry_run else str(args.output.expanduser().resolve().with_suffix(".artifact.json"))
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

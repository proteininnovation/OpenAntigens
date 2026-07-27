from __future__ import annotations

import hashlib
import json
from importlib import resources
from pathlib import Path
from typing import Any

from .models import GPCREngineeringCassette


class GPCREngineeringRegistry:
    def __init__(self, cassettes: list[GPCREngineeringCassette], warnings: list[str] | None = None) -> None:
        self.cassettes = {cassette.cassette_id: cassette for cassette in cassettes}
        self.warnings = warnings or []

    def get(self, cassette_id: str) -> GPCREngineeringCassette | None:
        return self.cassettes.get(cassette_id)

    def validated_for_mode(self, cassette_id: str, mode: str) -> tuple[GPCREngineeringCassette | None, str | None]:
        cassette = self.get(cassette_id)
        if cassette is None:
            return None, f"GPCR engineering cassette `{cassette_id}` is not present in the registry."
        if mode not in cassette.allowed_use_modes:
            return None, f"GPCR engineering cassette `{cassette_id}` is not allowed for mode `{mode}`."
        if not cassette.sequence:
            return None, f"GPCR engineering cassette `{cassette_id}` has no approved sequence in the registry."
        expected = cassette.checksum_sha256
        observed = hashlib.sha256(cassette.sequence.encode("ascii")).hexdigest()
        if expected and observed != expected:
            return None, (
                f"GPCR engineering cassette `{cassette_id}` checksum mismatch "
                f"(expected {expected}, observed {observed})."
            )
        if not expected:
            return None, f"GPCR engineering cassette `{cassette_id}` has no checksum in the registry."
        return cassette, None


def load_gpcr_engineering_registry(path: str | Path | None = None) -> GPCREngineeringRegistry:
    payload = _load_registry_payload(path)
    cassettes: list[GPCREngineeringCassette] = []
    warnings: list[str] = []
    for item in payload.get("cassettes") or []:
        try:
            cassettes.append(
                GPCREngineeringCassette(
                    cassette_id=str(item.get("cassette_id") or "").strip(),
                    name=str(item.get("name") or "").strip(),
                    sequence=_clean_sequence(item.get("sequence")),
                    source_note=str(item.get("source_note") or "").strip(),
                    version=str(item.get("version") or payload.get("version") or "").strip(),
                    checksum_sha256=(str(item.get("checksum_sha256")).strip() if item.get("checksum_sha256") else None),
                    default_n_linker=_clean_sequence(item.get("default_n_linker")) or "",
                    default_c_linker=_clean_sequence(item.get("default_c_linker")) or "",
                    allowed_use_modes=[str(mode) for mode in item.get("allowed_use_modes") or []],
                )
            )
        except Exception as exc:
            warnings.append(f"Skipped invalid GPCR cassette registry entry: {exc}")
    return GPCREngineeringRegistry(cassettes=cassettes, warnings=warnings)


def _load_registry_payload(path: str | Path | None) -> dict[str, Any]:
    if path:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    resource = resources.files("agdesign2").joinpath("data/gpcr_engineering_cassettes.json")
    return json.loads(resource.read_text(encoding="utf-8"))


def _clean_sequence(value: Any) -> str | None:
    if value is None:
        return None
    sequence = "".join(str(value).split()).upper()
    return sequence or None

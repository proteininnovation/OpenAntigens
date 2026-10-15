from __future__ import annotations

from dataclasses import dataclass
import re

from .config import AnalysisConfig
from .models import AnalysisNote, Feature, Region, TopologyAnnotation


@dataclass(slots=True)
class TopologyResult:
    ectodomain: Region | None
    notes: list[AnalysisNote]
    topology: TopologyAnnotation | None = None


# Legacy spellings that reports written before the OA-018 alias fix carry.
_FEATURE_TYPE_ALIASES = {
    "PROPEPTIDE": "PROPEP",
    "LIPIDATION": "LIPID",
    "INTRAMEMBRANE": "INTRAMEM",
    "TRANSIT_PEPTIDE": "TRANSIT",
}


def feature_type(feature: Feature) -> str:
    """The canonical type, accepting the pre-OA-018 spellings on read."""
    raw = str(feature.type or "").upper()
    return _FEATURE_TYPE_ALIASES.get(raw, raw)


def is_beta_stranded_transmem(feature: Feature) -> bool:
    """A UniProt TRANSMEM that is a beta strand, not a helix.

    F10a: every TRANSMEM was counted as a helix and the description was never
    read. MACPF and gasdermin beta-hairpins and the VDAC1 beta-barrel are
    annotated "Beta stranded", so C6/C8A/C8B/C9/PRF1/MPEG1, GSDMA/B/D and VDAC1
    — 10 human / 9 mouse reports — became Track Multipass with a blank
    ectodomain, while C7, which has no TRANSMEM at all, is correctly Secreted
    23-843. These hairpins only insert into a target membrane on activation;
    the soluble protein is the antigen.
    """
    return "beta" in str(feature.description or "").lower()


def derive_ectodomain(
    features: list[Feature],
    sequence_length: int,
    config: AnalysisConfig,
) -> TopologyResult:
    """Derive the design region, then cap it at a GPI omega site.

    F16b residual: 7 reports carry a GPI anchor **and** a TRANSMEM, so they take
    the single-pass branch and never reach the cap inside
    ``_derive_gpi_anchored_region``. Capping once, after the branch has chosen,
    reaches all of them. Measured over both snapshots: 12 human / 116 mouse
    regions ran past the omega site, 6 / 6 survived the chain-branch cap alone,
    0 survive this.
    """
    result = _derive_ectodomain(features, sequence_length, config)
    return _cap_region_at_gpi_omega_site(result, features)


def _cap_region_at_gpi_omega_site(
    result: TopologyResult, features: list[Feature]
) -> TopologyResult:
    """Trim a design region that runs past a GPI omega site — conditionally.

    F16b residual: 6 human / 6 mouse reports carry a GPI anchor **and** a
    TRANSMEM, so they take the single-pass branch and never reach the cap
    inside ``_derive_gpi_anchored_region``.

    The condition matters, and it is the one the plan's §7 release check names:
    only cap **where a CHAIN bounds the GPI region**. An *untrimmed* region —
    no CHAIN at all, as on CD24 and in
    ``test_gpi_anchor_without_processed_chain_retains_untrimmed_construct`` —
    is left whole, because without a processed chain there is nothing saying
    the tail is a propeptide rather than part of the protein.

    Where a CHAIN does contain the omega site the biology is not in doubt: the
    transamidase cleaves there and replaces the C-terminal signal with the
    anchor, so those residues are not in the mature protein. That applies to
    TM+GPI proteins too (tetherin/BST2, ALPL, ULBP5), which take the
    single-pass branch and never reach the cap inside
    ``_derive_gpi_anchored_region``.
    """
    region = result.ectodomain
    if region is None:
        return result
    anchors = [
        feature.start
        for feature in features
        if feature_type(feature) == "LIPID" and "gpi" in (feature.description or "").lower()
    ]
    chains = [feature for feature in features if feature_type(feature) == "CHAIN"]
    if not anchors or not chains:
        return result
    inside = [
        anchor
        for anchor in anchors
        if region.start <= anchor < region.end
        and any(chain.start <= anchor <= chain.end for chain in chains)
    ]
    if not inside:
        return result
    omega = max(inside)
    result.notes.append(
        AnalysisNote(
            severity="info",
            message=(
                f"Design region trimmed to {omega}: the GPI transamidase removes everything "
                f"after the omega site, so residues {omega + 1}-{region.end} are not part of "
                "the mature protein."
            ),
            source="topology",
        )
    )
    result.ectodomain = replace_region(region, end=omega)
    return result


def _derive_ectodomain(
    features: list[Feature],
    sequence_length: int,
    config: AnalysisConfig,
) -> TopologyResult:
    notes: list[AnalysisNote] = []
    tm_regions = [
        feature
        for feature in features
        if feature.type.upper() == "TRANSMEM" and not is_beta_stranded_transmem(feature)
    ]
    intramembrane_regions = [feature for feature in features if feature.type.upper() == "INTRAMEMBRANE"]
    signal_regions = [feature for feature in features if feature.type.upper() == "SIGNAL"]
    topo_regions = [feature for feature in features if feature.type.upper() == "TOPO_DOM"]
    if len(tm_regions) == 0:
        gpi_anchors = [
            feature.start
            for feature in features
            # Both spellings: reports written before the OA-018 alias say
            # LIPIDATION, new ones say LIPID. Testing only one of them would
            # silently stop detecting GPI anchors.
            if feature_type(feature) == "LIPID" and "gpi" in (feature.description or "").lower()
        ]
        if gpi_anchors:
            gpi_region = _derive_gpi_anchored_region(features, sequence_length, gpi_anchors)
            used_processed_chain = gpi_region is not None
            if gpi_region is None:
                start = max((region.end for region in signal_regions), default=0) + 1
                gpi_region = Region(
                    start=start,
                    end=sequence_length,
                    label="GPI-anchored sequence (untrimmed)",
                    source="topology",
                    confidence=0.5,
                    metadata={"fallback_reason": "missing_processed_chain"},
                )
            topology = _build_topology_annotation(
                tm_regions=tm_regions,
                intramembrane_regions=intramembrane_regions,
                topo_regions=topo_regions,
                signal_regions=signal_regions,
                sequence_length=sequence_length,
                config=config,
                major_extracellular_region=gpi_region,
                topology_class="gpi_anchored" if used_processed_chain else "gpi_anchored_untrimmed",
            )
            topology.extracellular_regions = [gpi_region]
            notes.append(
                AnalysisNote(
                    severity="info" if used_processed_chain else "warning",
                    message=(
                        "Using the processed UniProt chain boundary for the GPI-anchored extracellular design region."
                        if used_processed_chain
                        else "A GPI anchor is annotated, but UniProt provides no processed chain containing the anchor; retaining the post-signal sequence without trimming the GPI propeptide."
                    ),
                    source="topology",
                )
            )
            return TopologyResult(ectodomain=gpi_region, notes=notes, topology=topology)
        secreted_region = _derive_secreted_region(signal_regions, topo_regions, sequence_length, features)
        topology = _build_topology_annotation(
            tm_regions=tm_regions,
            intramembrane_regions=intramembrane_regions,
            topo_regions=topo_regions,
            signal_regions=signal_regions,
            sequence_length=sequence_length,
            config=config,
            major_extracellular_region=secreted_region,
            topology_class="secreted" if secreted_region is not None else "non_membrane",
        )
        if secreted_region is not None:
            notes.append(
                AnalysisNote(
                    severity="info",
                    message="No transmembrane region detected; treating the mature secreted chain as the extracellular design region.",
                    source="topology",
                )
            )
            return TopologyResult(ectodomain=secreted_region, notes=notes, topology=topology)
        return TopologyResult(ectodomain=None, notes=notes, topology=topology)
    if len(tm_regions) != 1:
        rescued_region = _derive_from_dominant_extracellular_topology(
            tm_regions=tm_regions,
            signal_regions=signal_regions,
            topo_regions=topo_regions,
            sequence_length=sequence_length,
            config=config,
        )
        topology = _build_topology_annotation(
            tm_regions=tm_regions,
            intramembrane_regions=intramembrane_regions,
            topo_regions=topo_regions,
            signal_regions=signal_regions,
            sequence_length=sequence_length,
            config=config,
            major_extracellular_region=rescued_region,
        )
        if topology.major_extracellular_region is not None:
            notes.append(
                AnalysisNote(
                    severity="info",
                    message=(
                        f"Detected {len(tm_regions)} transmembrane annotations; using the major extracellular region "
                        "for soluble construct design and retaining the membrane core for future membrane-expression designs."
                    ),
                    source="topology",
                )
            )
        if rescued_region is not None:
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message=(
                        f"Detected {len(tm_regions)} transmembrane annotations; "
                        "using the dominant extracellular topological domain preceding the terminal helix."
                    ),
                    source="topology",
                )
            )
            return TopologyResult(ectodomain=rescued_region, notes=notes, topology=topology)
        if topology.topology_class == "multipass_mixed" and topology.major_extracellular_region is not None:
            return TopologyResult(
                ectodomain=topology.major_extracellular_region,
                notes=notes,
                topology=topology,
            )
        # F10b (1,531 human / 1,347 mouse pages): this described a failed
        # precondition of the single-pass derivation, not a problem with the
        # report. The class is already multipass_compact by this point, so the
        # page showed Track "Multipass" beside a WARNING saying the protein was
        # not what it had just been classified as. Say what happened instead.
        notes.append(
            AnalysisNote(
                severity="info",
                message=(
                    f"{len(tm_regions)} transmembrane regions are annotated; no single "
                    "extracellular domain could be assigned, so only full-length membrane "
                    "constructs are offered."
                ),
                source="topology",
            )
        )
        return TopologyResult(ectodomain=None, notes=notes, topology=topology)

    transmembrane = tm_regions[0]
    extracellular_n = _is_extracellular_before(transmembrane.start, topo_regions, config)
    extracellular_c = _is_extracellular_after(transmembrane.end, topo_regions, config)
    ectodomain: Region | None = None

    if extracellular_n and extracellular_c:
        notes.append(
            AnalysisNote(
                severity="warning",
                message="Both N- and C-terminal sides appear extracellular; ectodomain cannot be assigned confidently.",
                source="topology",
            )
        )
        topology = _build_topology_annotation(
            tm_regions=tm_regions,
            intramembrane_regions=intramembrane_regions,
            topo_regions=topo_regions,
            signal_regions=signal_regions,
            sequence_length=sequence_length,
            config=config,
        )
        return TopologyResult(ectodomain=None, notes=notes, topology=topology)

    if not extracellular_n and not extracellular_c:
        if signal_regions:
            extracellular_n = True
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message="Topology annotations are incomplete; inferred an N-terminal ectodomain from the signal peptide.",
                        source="topology",
                )
            )
        else:
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message="Unable to determine which side of the transmembrane helix is extracellular.",
                        source="topology",
                )
            )
            topology = _build_topology_annotation(
                tm_regions=tm_regions,
                intramembrane_regions=intramembrane_regions,
                topo_regions=topo_regions,
                signal_regions=signal_regions,
                sequence_length=sequence_length,
                config=config,
            )
            return TopologyResult(ectodomain=None, notes=notes, topology=topology)

    if extracellular_n:
        start = 1
        if signal_regions:
            start = max(region.end for region in signal_regions) + 1
        end = transmembrane.start - 1
        label = "N-terminal ectodomain"
    else:
        start = transmembrane.end + 1
        end = sequence_length
        label = "C-terminal ectodomain"

    if start >= end:
        notes.append(
            AnalysisNote(
                severity="warning",
                message="Computed ectodomain bounds are empty.",
                source="topology",
            )
        )
        topology = _build_topology_annotation(
            tm_regions=tm_regions,
            intramembrane_regions=intramembrane_regions,
            topo_regions=topo_regions,
            signal_regions=signal_regions,
            sequence_length=sequence_length,
            config=config,
        )
        return TopologyResult(ectodomain=None, notes=notes, topology=topology)

    ectodomain = Region(
        start=start,
        end=end,
        label=label,
        source="topology",
        confidence=0.8 if topo_regions else 0.5,
    )
    topology = _build_topology_annotation(
        tm_regions=tm_regions,
        intramembrane_regions=intramembrane_regions,
        topo_regions=topo_regions,
        signal_regions=signal_regions,
        sequence_length=sequence_length,
        config=config,
        major_extracellular_region=ectodomain,
        topology_class="single_pass",
    )
    return TopologyResult(
        ectodomain=ectodomain,
        notes=notes,
        topology=topology,
    )


def _derive_from_dominant_extracellular_topology(
    *,
    tm_regions: list[Feature],
    signal_regions: list[Feature],
    topo_regions: list[Feature],
    sequence_length: int,
    config: AnalysisConfig,
) -> Region | None:
    if len(tm_regions) < 2:
        return None
    extracellular_topologies = [
        feature
        for feature in topo_regions
        if _looks_extracellular(feature.description, config)
    ]
    if not extracellular_topologies:
        return None
    terminal_tm = max(tm_regions, key=lambda feature: feature.end)
    if terminal_tm.end < int(sequence_length * 0.75):
        return None
    signal_end = max((region.end for region in signal_regions), default=0)
    candidate_regions = [
        feature
        for feature in extracellular_topologies
        if feature.end >= terminal_tm.start - 5 and feature.start <= max(signal_end + 10, 25)
    ]
    if not candidate_regions:
        return None
    candidate = max(candidate_regions, key=lambda feature: (feature.end - feature.start, -feature.start))
    start = max(signal_end + 1, candidate.start)
    end = min(candidate.end, terminal_tm.start - 1)
    if end - start + 1 < max(config.multipass_mixed_ectodomain_min_length, config.min_construct_length):
        return None
    return Region(
        start=start,
        end=end,
        label="Dominant extracellular region",
        source="topology",
        confidence=0.55,
    )


def _build_topology_annotation(
    *,
    tm_regions: list[Feature],
    intramembrane_regions: list[Feature] | None = None,
    topo_regions: list[Feature],
    signal_regions: list[Feature],
    sequence_length: int,
    config: AnalysisConfig,
    major_extracellular_region: Region | None = None,
    topology_class: str | None = None,
) -> TopologyAnnotation:
    extracellular_regions = [
        _feature_region(feature, source="topology")
        for feature in topo_regions
        if _looks_extracellular(feature.description, config)
    ]
    cytoplasmic_regions = [
        _feature_region(feature, source="topology")
        for feature in topo_regions
        if _looks_cytoplasmic(feature.description, config)
    ]
    tm_region_annotations = [_feature_region(feature, source="UniProt") for feature in tm_regions]
    intramembrane_region_annotations = [
        _feature_region(feature, source="UniProt")
        for feature in (intramembrane_regions or [])
    ]
    if major_extracellular_region is None:
        major_extracellular_region = _select_major_extracellular_region(
            extracellular_regions=extracellular_regions,
            tm_regions=tm_region_annotations,
            signal_regions=signal_regions,
            sequence_length=sequence_length,
            config=config,
        )
    membrane_core_region = None
    if tm_region_annotations:
        membrane_core_region = Region(
            start=min(region.start for region in tm_region_annotations),
            end=max(region.end for region in tm_region_annotations),
            label="Membrane core",
            source="topology",
            confidence=0.8,
        )
    cytoplasmic_tail_region = _select_cytoplasmic_tail_region(
        cytoplasmic_regions=cytoplasmic_regions,
        tm_regions=tm_region_annotations,
    )
    if topology_class is None:
        if not tm_region_annotations:
            topology_class = "secreted" if major_extracellular_region is not None else "non_membrane"
        elif len(tm_region_annotations) == 1:
            topology_class = "single_pass"
        elif major_extracellular_region is not None:
            topology_class = "multipass_mixed"
        else:
            topology_class = "multipass_compact"
    return TopologyAnnotation(
        topology_class=topology_class,
        extracellular_regions=extracellular_regions,
        transmembrane_regions=tm_region_annotations,
        intramembrane_regions=intramembrane_region_annotations,
        cytoplasmic_regions=cytoplasmic_regions,
        major_extracellular_region=major_extracellular_region,
        membrane_core_region=membrane_core_region,
        cytoplasmic_tail_region=cytoplasmic_tail_region,
    )


def _select_major_extracellular_region(
    *,
    extracellular_regions: list[Region],
    tm_regions: list[Region],
    signal_regions: list[Feature],
    sequence_length: int,
    config: AnalysisConfig,
) -> Region | None:
    if not extracellular_regions:
        return None
    terminal_tm_start = max((region.start for region in tm_regions), default=sequence_length + 1)
    signal_end = max((region.end for region in signal_regions), default=0)
    ranked = sorted(
        extracellular_regions,
        key=lambda region: (
            region.end - region.start + 1,
            1 if region.start <= signal_end + 10 else 0,
            1 if region.end < terminal_tm_start else 0,
        ),
        reverse=True,
    )
    candidate = ranked[0]
    min_length = (
        max(config.multipass_mixed_ectodomain_min_length, config.min_construct_length)
        if len(tm_regions) > 1
        else max(config.min_construct_length, 50)
    )
    if (candidate.end - candidate.start + 1) < min_length:
        return None
    return replace_region(
        candidate,
        label="Major extracellular region" if len(tm_regions) > 1 else candidate.label,
    )


def _select_cytoplasmic_tail_region(
    *,
    cytoplasmic_regions: list[Region],
    tm_regions: list[Region],
) -> Region | None:
    if not cytoplasmic_regions or not tm_regions:
        return None
    terminal_tm_end = max(region.end for region in tm_regions)
    trailing = [region for region in cytoplasmic_regions if region.end > terminal_tm_end]
    if not trailing:
        return None
    candidate = max(trailing, key=lambda region: (region.end, region.end - region.start + 1))
    if candidate.start <= terminal_tm_end:
        candidate = replace_region(candidate, start=terminal_tm_end + 1)
    if candidate.start > candidate.end:
        return None
    return replace_region(candidate, label="Cytoplasmic tail")


def _feature_region(feature: Feature, *, source: str) -> Region:
    return Region(
        start=feature.start,
        end=feature.end,
        label=feature.description or feature.type.title(),
        source=source,
        confidence=0.8,
    )


def apply_surfy_topology_fallback(
    topology: TopologyAnnotation | None,
    surfy_topology: str | None,
    *,
    sequence_length: int,
) -> bool:
    """Fill missing extracellular/cytoplasmic topology regions from SURFY topology."""

    if topology is None or not surfy_topology:
        return False
    segments = _parse_surfy_topology_segments(surfy_topology, sequence_length=sequence_length)
    if not segments:
        return False
    changed = False
    if not topology.extracellular_regions:
        topology.extracellular_regions.extend(region for kind, region in segments if kind == "NC")
        changed = changed or bool(topology.extracellular_regions)
    if not topology.cytoplasmic_regions:
        topology.cytoplasmic_regions.extend(region for kind, region in segments if kind == "CY")
        changed = changed or bool(topology.cytoplasmic_regions)
    if not topology.transmembrane_regions:
        topology.transmembrane_regions.extend(region for kind, region in segments if kind == "TM")
        changed = changed or bool(topology.transmembrane_regions)
    if changed:
        # F02-B. These segments arrive AFTER _build_topology_annotation has
        # already computed the tail and the loops, so without this the tail
        # stays None (84 human / 227 mouse reports) and a C-terminal SURFY
        # segment is then offered as an internal cytoplasmic loop. SURFY
        # segments also straddle UniProt helices, which is how a fusion window
        # reached inside TM7 on or10ac1_mouse.
        _clip_regions_to_transmembrane(
            topology.cytoplasmic_regions, topology.transmembrane_regions
        )
        _clip_regions_to_transmembrane(
            topology.extracellular_regions, topology.transmembrane_regions
        )
        topology.cytoplasmic_tail_region = _select_cytoplasmic_tail_region(
            cytoplasmic_regions=topology.cytoplasmic_regions,
            tm_regions=topology.transmembrane_regions,
        )
    return changed


def _clip_regions_to_transmembrane(regions: list[Region], tm_regions: list[Region]) -> None:
    """Trim each region so it does not reach inside a transmembrane helix.

    A region that a helix splits in two is replaced by the longer half; one
    entirely inside a helix is dropped.
    """
    if not tm_regions:
        return
    clipped: list[Region] = []
    for region in regions:
        pieces = [(region.start, region.end)]
        for helix in tm_regions:
            next_pieces: list[tuple[int, int]] = []
            for start, end in pieces:
                if end < helix.start or start > helix.end:
                    next_pieces.append((start, end))
                    continue
                if start < helix.start:
                    next_pieces.append((start, helix.start - 1))
                if end > helix.end:
                    next_pieces.append((helix.end + 1, end))
            pieces = next_pieces
        if not pieces:
            continue
        start, end = max(pieces, key=lambda piece: piece[1] - piece[0])
        if (start, end) == (region.start, region.end):
            clipped.append(region)
        else:
            clipped.append(replace_region(region, start=start, end=end))
    regions[:] = clipped


def _parse_surfy_topology_segments(
    topology_text: str,
    *,
    sequence_length: int,
) -> list[tuple[str, Region]]:
    segments: list[tuple[str, Region]] = []
    labels = {
        "NC": "SURFY non-cytoplasmic region",
        "CY": "SURFY cytoplasmic region",
        "TM": "SURFY transmembrane region",
    }
    for token in str(topology_text or "").split(";"):
        match = re.fullmatch(r"\s*(NC|CY|TM)\s*:\s*(\d+)\s*-\s*(\d+)\s*", token)
        if not match:
            continue
        kind, start_text, end_text = match.groups()
        start = max(1, int(start_text))
        end = min(sequence_length, int(end_text))
        if start > end:
            continue
        segments.append(
            (
                kind,
                Region(
                    start=start,
                    end=end,
                    label=labels[kind],
                    source="SURFY",
                    confidence=0.65,
                    metadata={"surfy_segment_type": kind},
                ),
            )
        )
    return segments


def replace_region(region: Region, **changes: object) -> Region:
    data = region.to_dict()
    data.update(changes)
    return Region(**data)


def _is_extracellular_before(position: int, topo_regions: list[Feature], config: AnalysisConfig) -> bool:
    for feature in topo_regions:
        if feature.end < position and _looks_extracellular(feature.description, config):
            return True
    return False


def _is_extracellular_after(position: int, topo_regions: list[Feature], config: AnalysisConfig) -> bool:
    for feature in topo_regions:
        if feature.start > position and _looks_extracellular(feature.description, config):
            return True
    return False


def _looks_extracellular(description: str | None, config: AnalysisConfig) -> bool:
    text = (description or "").lower()
    return any(keyword in text for keyword in config.topology_keywords_extracellular)


def _looks_cytoplasmic(description: str | None, config: AnalysisConfig) -> bool:
    text = (description or "").lower()
    return any(keyword in text for keyword in config.topology_keywords_cytoplasmic)


def _derive_secreted_region(
    signal_regions: list[Feature],
    topo_regions: list[Feature],
    sequence_length: int,
    features: list[Feature] | None = None,
) -> Region | None:
    """The post-signal region.

    F16c: the label said "Secreted mature chain" on 306 human / 263 mouse
    regions that still contain a PROPEPTIDE — BMP8B 20-402 with a propeptide at
    20-263 and the mature chain at 264-402. The interval is deliberately kept
    (a pro-domain is often wanted for expression; trimming is a design decision,
    not a bug), but it is no longer called the mature chain.
    """
    start = 1
    confidence = 0.55
    label = "Secreted mature chain"
    if signal_regions:
        start = max(region.end for region in signal_regions) + 1
        confidence = 0.8
    extracellular_topology = [
        feature
        for feature in topo_regions
        if feature.start <= start and feature.end >= start
    ]
    if extracellular_topology:
        start = min(feature.start for feature in extracellular_topology)
        confidence = max(confidence, 0.85)
        label = "Secreted extracellular chain"
    if not signal_regions and not extracellular_topology:
        return None
    if start >= sequence_length:
        return None
    metadata: dict[str, object] = {}
    propeptides = [
        feature
        for feature in (features or [])
        if feature_type(feature) == "PROPEP" and feature.end >= start
    ]
    chains = [
        feature
        for feature in (features or [])
        if feature_type(feature) == "CHAIN" and feature.end >= start
    ]
    if propeptides:
        label = "Post-signal precursor region"
        metadata["propeptide_spans"] = [
            [feature.start, feature.end] for feature in propeptides
        ]
        if chains:
            metadata["processed_chain_spans"] = [
                [feature.start, feature.end] for feature in chains
            ]
    return Region(
        start=start,
        end=sequence_length,
        label=label,
        source="topology",
        confidence=confidence,
        metadata=metadata,
    )


def _derive_gpi_anchored_region(
    features: list[Feature], sequence_length: int, anchors: list[int]
) -> Region | None:
    chains = [
        feature
        for feature in features
        if feature.type.upper() == "CHAIN" and 1 <= feature.start <= feature.end <= sequence_length
    ]
    if not chains or not any(chain.start <= anchor <= chain.end for chain in chains for anchor in anchors):
        return None
    start = min(chain.start for chain in chains)
    end = max(chain.end for chain in chains)
    # F16b: UniProt gives some of these entries a CHAIN covering the whole
    # sequence, so the region ran past the omega site and kept the propeptide
    # that the GPI transamidase removes (CD14 20-367 with the site at 345).
    # Conditional, so an untrimmed region that is already correct is untouched.
    omega = max(anchor for anchor in anchors if start <= anchor <= end)
    if omega < end:
        end = omega
    if end <= start:
        return None
    return Region(
        start=start,
        end=end,
        label="GPI-anchored mature chain",
        source="UniProt",
        confidence=0.9,
    )

# OpenAntigens

> OpenAntigens, from the Institute for Protein Innovation, provides proposed antigen constructs and supporting sequence, topology, and structural evidence for human cell-surface and secreted proteins. A linked mouse portal covers resolved mouse orthologs. Use this guide to find records, interpret evidence, and preserve provenance.

Website: https://openantigens.org/

OpenAntigens supports research construct selection. A proposed construct is not evidence of successful expression, folding, antigenicity, antibody binding, or therapeutic efficacy. Available evidence varies by target and release.

## Start here

- [Browse human targets](https://openantigens.org/index.html): Search by gene symbol, UniProt entry name, protein name, alias, family, or topology. Use the separate disease filter for indexed Open Targets associations.
- [Help](https://openantigens.org/help.html): Report sections, search, evidence availability, and common workflows.
- [Construct builder guide](https://openantigens.org/builder.html): Select a proposed construct, inspect boundaries, and export the current selection.
- [Construct selection and methodology](https://openantigens.org/constructs.html): Construct classes, boundary generation, interpretation, and considerations after export.
- [Methods](https://openantigens.org/methods.html): Source databases, sequence matching, confidence metrics, ortholog mapping, and limitations.
- [Mouse portal](https://openantigens.org/mouse/index.html): Mouse targets derived from the human target set. Use mouse-specific reports and provenance when working with this portal.

## Find and verify a target

1. Search the gene or UniProt entry name, then open the report. Confirm the organism, accession, canonical sequence, and topology in the report before interpreting a construct. Similar names and paralogs are not interchangeable. Check the report's sequence source if a mapping uses RefSeq or another reference.
2. For programmatic lookup, use the JSON or TSV index below. Follow its `detail_page` path, resolved against the relevant portal root. Human index paths are relative to https://openantigens.org/, not to the downloads directory. Use links supplied by the site rather than inventing report filenames or API endpoints.
3. Inspect the report's evidence and warnings. `status` is a processing state; `ok` is not experimental validation. Other status values can describe reused or failed analyses. Check the linked report and any error field before interpreting them.
4. If the target or evidence is absent, report what was unavailable in the inspected release. Missing data, a blank field, an unavailable section, or zero indexed hits does not establish biological absence or specificity.

## Read boundaries and exports correctly

- Report and builder target boundaries use 1-based, inclusive coordinates in the canonical report sequence. For a single unmodified interval, length is end minus start plus one. Preserve accession, species, start, end, and sequence together. Do not silently renumber a mature chain from residue 1.
- PDB chain numbering and GPCR generic numbering are separate coordinate systems. Use the displayed mapping to the report sequence; do not substitute raw structure coordinates for UniProt boundaries.
- Confirm the selected species and current boundaries before copying TSV or FASTA. Browser edits can change the exported sequence. Retain any displayed mutation, sequence-source, or mapping annotations, and describe user-edited boundaries as edits.
- Use the report's export controls for exact sequences. The compact release index is not a construct-sequence catalogue. Do not reconstruct a sequence from the index's `ectodomain` text or infer a missing ortholog sequence.
- Species-equivalent intervals are alignment projections. Use each species' own accession, boundaries, and sequence; do not transfer human residue numbers directly. Sequence identity does not establish antibody cross-reactivity.
- Review topology, processing sites, domain boundaries, and the report's cysteine, cleavage, glycosylation, or assembly notes where present. Keep warnings about uncertain or untrimmed GPI boundaries with the selected record. Treat membrane-expression constructs separately from soluble constructs.
- Consult the after-export section of the construct guide for vector, secretion-leader, tag, and host considerations. OpenAntigens does not prescribe experimental expression or purification conditions.

## Interpret supporting evidence

- Full-region, PDB, strict, lenient, and annotated-domain constructs represent different sources or selection rules. Their display order is not a validated ranking of expression performance. Read the construct guide for the method used by the inspected release.
- AlphaFold structures are predictions. pLDDT describes local model confidence on a 0-100 scale; PAE describes uncertainty in relative residue placement in angstroms. Neither metric measures affinity, expression, or experimental structural accuracy. High local confidence alone does not establish interdomain orientation.
- A PDB link supplies experimental structure context. Mapped chain coverage does not establish that the proposed exported construct reproduces the deposited sequence, engineering, partners, or expression conditions. Inspect the original PDB record for those details.
- Surface-accessibility values come from the structural model and its calculation context. Missing glycans or partner chains can change biological exposure; model accessibility does not establish an accessible antibody epitope in cells.
- BLAST hits and family or paralog identity support similarity review. Interpret identity with coverage, alignment, species, and the region compared. `cross_reactivity_count` counts indexed BLAST hits, not measured cross-reactive antibodies. Whole-sequence, construct, and extracellular-accessible identity use different regions or denominators.
- Open Targets scores describe associations and their supporting evidence. Distinguish direct from indirect scores. They are not probabilities of clinical success or proof of a causal mechanism. PubTator counts describe literature retrieval, not evidence quality or target validation.

## Data access and release provenance

- [Downloads and field definitions](https://openantigens.org/downloads.html): Available release artifacts and index schema.
- [JSON target index](https://openantigens.org/downloads/agdesign2_portal_index.json): Compact target records with identifiers, processing states, evidence flags, counts, and report paths.
- [TSV target index](https://openantigens.org/downloads/agdesign2_portal_index.tsv): The same compact index for tabular analysis.
- [Download manifest](https://openantigens.org/downloads/download_manifest.json): Included data files, formats, columns, and source notes. Retrieve optional disease files only when listed or linked in the inspected release.
- [Portal metadata](https://openantigens.org/portal_metadata.json): Software version, portal build date, generation timestamp, counts, and license summary.

Prefer the release index for bulk lookup instead of scraping every report. The public site publishes static artifacts; do not assume development-server JSON APIs exist on the public host. If page text extraction omits a JavaScript-rendered table, use the downloadable index or a browser capable of rendering the page.

Record the portal build date and software version from metadata, the artifact URL or filename, and your retrieval date. Preserve original downloads when performing an analysis. Do not combine files from different releases without recording the difference. If counts, identifiers, dates, or sequences disagree between a page and a download, disclose the mismatch and verify the affected record before using it.

For a target-specific answer, include the gene or entry name, species, accession, report URL, and inspected release. When discussing a construct, also retain its name or class, boundaries, evidence source, and relevant warnings. Separate reported annotations from your interpretation and from recommendations for later experiments. Cite only records and underlying sources you actually inspected.

## How to cite OpenAntigens

When describing or using OpenAntigens in research, please cite:

{{citation}}

This is a preprint, version 1, posted August 4, 2026; it has not been peer reviewed. For target-specific claims, also cite the inspected report or artifact and record its release build date, software version, and your access date. Cite underlying databases or primary studies where their evidence is used. Do not use the database paper as evidence that an individual construct has been experimentally validated.

- [How to cite](https://openantigens.org/help.html#cite-openantigens): Citation and release-provenance guidance.
- [BibTeX citation](https://openantigens.org/downloads/openantigens.bib): Import into a bibliography manager.
- [RIS citation](https://openantigens.org/downloads/openantigens.ris): Import into a bibliography manager.

## Attribution and reuse

- [Terms](https://openantigens.org/terms.html): Software, annotation, and third-party data terms, attribution, limitations, and correction contact.
- [Privacy](https://openantigens.org/privacy.html): Website analytics and data handling.
- [Official source repository](https://github.com/proteininnovation/OpenAntigens): Code and project documentation.

Recommended attribution: OpenAntigens, Institute for Protein Innovation, created by Andre A. R. Teixeira. Include the build date, artifact name, software version when available, and relevant source databases.

The software is Apache-2.0; OpenAntigens-generated annotations are CC BY 4.0. Third-party sequences, structures, and other source-derived fields retain their original terms and citation requirements. Consult the current terms and the underlying providers before redistribution; the annotation license does not relicense every field in a download.

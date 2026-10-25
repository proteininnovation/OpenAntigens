import { mkdtemp, mkdir, writeFile, readFile, rm } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, dirname } from "node:path";
import { createHash } from "node:crypto";
import { gzipSync } from "node:zlib";
import { semanticHash } from "../src/provider.js";
import { SECTIONS } from "../src/contracts.js";
const sha = (b) => createHash("sha256").update(b).digest("hex");
const canonical = (v) =>
  JSON.stringify(
    v === null || typeof v !== "object"
      ? v
      : Array.isArray(v)
        ? v.map((x) => JSON.parse(canonical(x)))
        : Object.fromEntries(
            Object.keys(v)
              .sort()
              .map((k) => [k, JSON.parse(canonical(v[k]))]),
          ),
  );
export async function fixture(
  snapshot = "fixture",
  { legacyRootSchema = false, inputHashes = {} } = {},
) {
  const publicSchema = JSON.parse(
    await readFile(new URL("./schema.json", import.meta.url)),
  );
  if (legacyRootSchema)
    delete publicSchema.$defs.root_manifest.properties.retired_releases;
  const directory = await mkdtemp(join(tmpdir(), "oa-api-test-"));
  const core = {
    schema_version: 1,
    provenance: {
      software_commit: "fixture",
      snapshot_id: snapshot,
      input_sha256: inputHashes,
      notes: [],
      source_dates: {},
      source_versions: {},
    },
    payload_sha256: {},
    asset_sha256: {},
  };
  let release_id = sha(canonical(core));
  const files = {};
  const species = {};
  const envelopes = [];
  const construct_id = sha("construct");
  const put = async (
    path,
    value,
    media_type = "application/json",
    schema = null,
  ) => {
    const bytes = Buffer.isBuffer(value)
      ? value
      : Buffer.from(canonical(value) + "\n");
    await mkdir(dirname(join(directory, path)), { recursive: true });
    await writeFile(join(directory, path), bytes);
    files[path] = {
      bytes: bytes.length,
      sha256: sha(bytes),
      media_type,
      encoding: path.endsWith(".gz") ? "gzip" : null,
      schema,
    };
    return path;
  };
  const shard = async (s, path, data, kind) => {
    const payload = { release_id, schema_version: 1, species: s, data };
    await put(
      path,
      payload,
      "application/json",
      kind === "section"
        ? `section_${data.section}`
        : {
            index: "index",
            target: "target_shard",
            sequence: "sequence_shard",
            constructs: "construct_page",
          }[kind],
    );
    envelopes.push({ kind, path, payload });
    return path;
  };
  for (const [s, taxon_id, accession] of [
    ["human", 9606, "P12345"],
    ["mouse", 10090, "Q12345"],
  ]) {
    const prefix = `releases/${release_id}/${s}`;
    const target_id = `${taxon_id}:uniprot:${accession}`;
    const sequence = "ACDEFGHIKLMNPQRSTVWY";
    const sections = Object.fromEntries(
      SECTIONS.map((section) => [
        section,
        {
          state: "unknown",
          reason: "Source query history is unavailable.",
          total: 0,
          pages: [],
        },
      ]),
    );
    sections.topology = {
      state: "available",
      reason: null,
      total: 1,
      pages: [
        await shard(
          s,
          `${prefix}/sections/topology-0.json`,
          {
            target_id,
            section: "topology",
            state: "available",
            reason: null,
            items: [
              {
                topology_class: "extracellular",
                source: "fixture",
                extracellular_regions: [
                  {
                    label: "extracellular",
                    start: 1,
                    end: 20,
                    source: "fixture",
                    confidence: null,
                  },
                ],
                transmembrane_regions: [],
                intramembrane_regions: [],
                cytoplasmic_regions: [],
                major_extracellular_region: null,
                membrane_core_region: null,
                cytoplasmic_tail_region: null,
              },
            ],
          },
          "section",
        ),
      ],
    };
    const summary = {
      target_id,
      record_id: sha(target_id),
      species: s,
      taxon_id,
      identifier_source: "uniprot",
      identifier_accession: accession,
      entry_name: `EXAMPLE_${s}`,
      gene_symbol: "EXAMPLE",
      protein_name: "Example protein",
      aliases: ["Alias"],
      sequence_accession: accession,
      sequence_source: "uniprot",
      sequence_sha256: sha(sequence),
      sequence_length: 20,
      status: "ok",
      status_raw: "ok",
      track: "secreted",
      track_raw: "Secreted",
      analyzed_topology: "extracellular",
      analyzed_topology_source: "fixture",
      design_region: { start: 1, end: 20, source: "fixture", label: "region" },
      construct_count: 1,
      structure_state: "unavailable",
      construct_state: "available",
      section_availability: Object.fromEntries(
        Object.entries(sections).map(([k, v]) => [k, v.state]),
      ),
      disease_ids: s === "human" ? ["EFO:0001"] : [],
      family_names: ["Example"],
      live_report_url: `https://openantigens.org/${s}/example.html`,
      search_text: `EXAMPLE Alias ${accession} Example protein`,
    };
    const seq = {
      target_id,
      construct_id: null,
      species: s,
      taxon_id,
      identifier_accession: accession,
      identifier_source: "uniprot",
      sequence_accession: accession,
      sequence_source: "uniprot",
      sequence_basis: "reference",
      start: 1,
      end: 20,
      length: 20,
      sequence,
      sequence_sha256: sha(sequence),
      warnings: ["Warning with unicode α; spaces and = / ? %."],
      evidence: ["fixture:evidence?scope=1"],
    };
    const reference = await shard(
      s,
      `${prefix}/sequences/reference.json`,
      seq,
      "sequence",
    );
    const cseq = await shard(
      s,
      `${prefix}/sequences/construct.json`,
      { ...seq, construct_id, sequence_basis: "reference_slice" },
      "sequence",
    );
    const construct = {
      construct_id,
      target_id,
      name: "full",
      category: "full_region",
      classification: "soluble_ecd",
      start: 1,
      end: 20,
      length: 20,
      sequence_basis: "reference_slice",
      sequence_sha256: sha(sequence),
      sequence_accession: accession,
      sequence_source: "uniprot",
      export_state: "available",
      export_reason: null,
      score: null,
      rationale: null,
      warnings: [],
      evidence: [],
      pdb_id: null,
      structural_metrics: null,
    };
    const cpage = await shard(
      s,
      `${prefix}/constructs/0.json`,
      [construct],
      "constructs",
    );
    const target = await shard(
      s,
      `${prefix}/targets/0.json`,
      {
        summary,
        sections,
        constructs: {
          state: "available",
          reason: null,
          total: 1,
          pages: [cpage],
        },
        reference_sequence_path: reference,
        construct_sequence_paths: { [construct_id]: cseq },
      },
      "target",
    );
    const unresolved = {
      ...summary,
      target_id: null,
      record_id: sha(`${s}unresolved`),
      identifier_accession: null,
      identifier_source: null,
      entry_name: "unresolved",
      gene_symbol: "UNRESOLVED",
      search_text: "unresolved",
      construct_count: 0,
      construct_state: "unavailable",
    };
    await shard(
      s,
      `${prefix}/targets/unresolved.json`,
      {
        summary: unresolved,
        sections: Object.fromEntries(
          SECTIONS.map((section) => [
            section,
            {
              state: "unknown",
              reason: "Source query history is unavailable.",
              total: 0,
              pages: [],
            },
          ]),
        ),
        constructs: { state: "no_results", reason: null, total: 0, pages: [] },
        reference_sequence_path: null,
        construct_sequence_paths: {},
      },
      "target",
    );
    const index_path = await shard(
      s,
      `${prefix}/index.json`,
      [summary, unresolved],
      "index",
    );
    const downloadPath = `${prefix}/downloads/records.jsonl.gz`;
    await put(
      downloadPath,
      gzipSync(
        Buffer.from(
          envelopes
            .filter((v) => v.payload.species === s)
            .map((v) => canonical(v))
            .join("\n") + "\n",
        ),
        { mtime: 0 },
      ),
      "application/gzip",
      null,
    );
    const file = files[downloadPath];
    const download = {
      name: "records.jsonl.gz",
      url: `https://openantigens.org/data/releases/${release_id}/${s}/records.jsonl.gz`,
      path: downloadPath,
      ...file,
      row_count: 6,
      license: "fixture",
    };
    species[s] = {
      taxon_id,
      index_path,
      target_paths: { [target_id]: target },
      downloads: [download],
      counts: {
        targets: 2,
        resolved_targets: 1,
        constructs: 1,
        status: { ok: 2 },
      },
      capabilities: {
        sections: SECTIONS.filter((v) => v !== "disease" || s === "human"),
        disease: s === "human",
      },
      records_download: download,
    };
  }
  await put(
    `releases/${release_id}/schema.json`,
    publicSchema,
    "application/schema+json",
    null,
  );
  const manifest = {
    release_id,
    schema_version: 1,
    provenance: core.provenance,
    citation: {
      title: "Fixture",
      authors: [],
      date: null,
      doi: null,
      journal: null,
      type: null,
      version: null,
    },
    licenses: {
      software: "MIT",
      annotations: "CC BY 4.0",
      third_party: "fixture",
    },
    species,
    files,
    core,
  };
  const byPath = new Map(
    envelopes.map((entry) => [entry.path, entry.payload.data]),
  );
  for (const name of ["human", "mouse"]) {
    const summaries = byPath.get(species[name].index_path);
    const projected = summaries.map((summary) => {
      const target = envelopes.find(
        (entry) =>
          entry.kind === "target" &&
          entry.payload.species === name &&
          entry.payload.data.summary.record_id === summary.record_id,
      ).payload.data;
      return {
        summary,
        sections: Object.fromEntries(
          Object.entries(target.sections).map(([section, descriptor]) => [
            section,
            {
              state: descriptor.state,
              reason: descriptor.reason,
              items: descriptor.pages.flatMap((path) => byPath.get(path).items),
            },
          ]),
        ),
        constructs: target.constructs.pages.flatMap((path) => byPath.get(path)),
        reference_sequence: target.reference_sequence_path
          ? byPath.get(target.reference_sequence_path)
          : null,
        construct_sequences: Object.fromEntries(
          Object.entries(target.construct_sequence_paths).map(([id, path]) => [
            id,
            path ? byPath.get(path) : null,
          ]),
        ),
      };
    });
    core.payload_sha256[name] = semanticHash(projected);
  }
  core.payload_sha256.schema = semanticHash(publicSchema);
  core.payload_sha256.citation = semanticHash(manifest.citation);
  const previousId = release_id;
  release_id = sha(canonical(core));
  const replaceId = (value) =>
    JSON.parse(JSON.stringify(value).replaceAll(previousId, release_id));
  const finalManifest = replaceId(manifest);
  finalManifest.release_id = release_id;
  const finalEnvelopes = replaceId(envelopes);
  await rm(join(directory, "releases", previousId), {
    recursive: true,
    force: true,
  });
  for (const record of finalEnvelopes) {
    const originalPath = record.path.replace(release_id, previousId);
    const bytes = Buffer.from(canonical(record.payload) + "\n");
    await mkdir(dirname(join(directory, record.path)), { recursive: true });
    await writeFile(join(directory, record.path), bytes);
    finalManifest.files[record.path] = {
      ...files[originalPath],
      bytes: bytes.length,
      sha256: sha(bytes),
    };
  }
  const schemaPath = `releases/${release_id}/schema.json`;
  const schemaBytes = Buffer.from(canonical(publicSchema) + "\n");
  await writeFile(join(directory, schemaPath), schemaBytes);
  finalManifest.files[schemaPath] = {
    ...files[`releases/${previousId}/schema.json`],
    bytes: schemaBytes.length,
    sha256: sha(schemaBytes),
  };
  for (const name of ["human", "mouse"]) {
    const download = finalManifest.species[name].records_download;
    const entries = finalEnvelopes.filter(
      (entry) => entry.payload.species === name,
    );
    const bytes = gzipSync(
      Buffer.from(entries.map(canonical).join("\n") + "\n"),
      { mtime: 0 },
    );
    await mkdir(dirname(join(directory, download.path)), { recursive: true });
    await writeFile(join(directory, download.path), bytes);
    Object.assign(download, {
      bytes: bytes.length,
      sha256: sha(bytes),
      row_count: entries.length,
    });
    finalManifest.species[name].downloads[0] = download;
    finalManifest.files[download.path] = {
      ...finalManifest.files[download.path],
      bytes: bytes.length,
      sha256: sha(bytes),
    };
  }
  const path = `releases/${release_id}/manifest.json`;
  const bytes = Buffer.from(canonical(finalManifest) + "\n");
  await writeFile(join(directory, path), bytes);
  await writeFile(
    join(directory, "manifest.json"),
    canonical({
      schema_version: 1,
      current_release_id: release_id,
      releases: [{ release_id, manifest_path: path, sha256: sha(bytes) }],
      ...(legacyRootSchema ? {} : { retired_releases: [] }),
    }) + "\n",
  );
  return {
    directory,
    release_id,
    construct_id,
    target_id: "9606:uniprot:P12345",
    manifest: finalManifest,
  };
}

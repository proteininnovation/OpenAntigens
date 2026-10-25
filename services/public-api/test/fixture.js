import { mkdtemp, mkdir, writeFile, readFile, rm } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, dirname } from "node:path";
import { createHash } from "node:crypto";
import { gzipSync } from "node:zlib";
import { semanticHash, fasta } from "../src/provider.js";
import { SECTIONS, DOWNLOAD_ROLES } from "../src/contracts.js";
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
  {
    legacyRootSchema = false,
    extraUnresolved = 0,
    unresolvedFirst = false,
    legacyIndex = false,
    inputHashes = {},
    sequenceWarnings = ["Warning with unicode α; spaces and = / ? %."],
  } = {},
) {
  const publicSchema = JSON.parse(
    await readFile(new URL("./schema.json", import.meta.url)),
  );
  if (legacyRootSchema)
    delete publicSchema.$defs.root_manifest.properties.retired_releases;
  const rowFields = [
    "target_id",
    "record_id",
    "species",
    "taxon_id",
    "search_text",
    "track",
    "status",
    "structure_state",
    "construct_state",
  ];
  publicSchema.$defs.search_row = {
    type: "object",
    additionalProperties: false,
    required: [...rowFields, "summary_path"],
    properties: {
      ...Object.fromEntries(
        rowFields.map((key) => [
          key,
          publicSchema.$defs.summary.properties[key],
        ]),
      ),
      summary_path: { type: "string", maxLength: 1024 },
    },
  };
  if (!legacyIndex) {
    publicSchema["x-openantigens-index-version"] = 2;
    publicSchema.$defs.index.properties.data.items = {
      $ref: "#/$defs/search_row",
    };
  } else {
    delete publicSchema["x-openantigens-index-version"];
    publicSchema.$defs.index.properties.data.items = {
      $ref: "#/$defs/summary",
    };
  }
  publicSchema["x-openantigens-fasta-version"] = 2;
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
    const targetHash = sha(canonical(target_id));
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
          `${prefix}/sections/${targetHash}/topology-0.json`,
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
      search_text: `catalog EXAMPLE Alias ${accession} Example protein`,
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
      warnings: sequenceWarnings,
      evidence: ["fixture:evidence?scope=1"],
    };
    const reference = await shard(
      s,
      `${prefix}/sequences/${targetHash}-reference.json`,
      seq,
      "sequence",
    );
    const cseq = await shard(
      s,
      `${prefix}/sequences/${construct_id}.json`,
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
      `${prefix}/constructs/${targetHash}-0.json`,
      [construct],
      "constructs",
    );
    const target = await shard(
      s,
      `${prefix}/targets/${targetHash}.json`,
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
      search_text: "catalog unresolved",
      section_availability: Object.fromEntries(
        SECTIONS.map((section) => [section, "unknown"]),
      ),
      construct_count: 0,
      construct_state: "unavailable",
    };
    const unresolvedPath = await shard(
      s,
      `${prefix}/targets/${sha(canonical(unresolved.record_id))}.json`,
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
    const unresolvedRows = [[unresolved, unresolvedPath]];
    for (let index = 0; index < extraUnresolved; index++) {
      const record = {
        ...unresolved,
        record_id: sha(`${s}:unresolved-extra:${index}`),
        gene_symbol: `UNRESOLVED_${index}`,
      };
      const path = await shard(
        s,
        `${prefix}/targets/${sha(canonical(record.record_id))}.json`,
        {
          summary: record,
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
          constructs: {
            state: "no_results",
            reason: null,
            total: 0,
            pages: [],
          },
          reference_sequence_path: null,
          construct_sequence_paths: {},
        },
        "target",
      );
      unresolvedRows.push([record, path]);
    }
    const sourceRows = [[summary, target], ...unresolvedRows];
    if (unresolvedFirst) sourceRows.push(sourceRows.shift());
    const index_path = await shard(
      s,
      `${prefix}/index.json`,
      legacyIndex
        ? sourceRows.map(([record]) => record)
        : sourceRows.map(([record, path]) => ({
            ...Object.fromEntries(rowFields.map((key) => [key, record[key]])),
            summary_path: path,
          })),
      "index",
    );
    species[s] = {
      taxon_id,
      index_path,
      target_paths: { [target_id]: target },
      downloads: [],
      counts: {
        targets: 2 + extraUnresolved,
        resolved_targets: 1,
        constructs: 1,
        status: { ok: 2 + extraUnresolved },
      },
      capabilities: {
        sections: SECTIONS.filter((v) => v !== "disease" || s === "human"),
        disease: s === "human",
      },
      records_download: null,
    };
  }
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
    const projected = summaries.map((row) => {
      const target = envelopes.find(
        (entry) =>
          entry.kind === "target" &&
          entry.payload.species === name &&
          entry.payload.data.summary.record_id === row.record_id,
      ).payload.data;
      return {
        summary: target.summary,
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
  const jsonl = (records) =>
    Buffer.from(records.map(canonical).join("\n") + "\n");
  const tsv = (records, fields) =>
    Buffer.from(
      fields.join("\t") +
        "\n" +
        records
          .map((record) =>
            fields
              .map((key) => {
                const value = record[key];
                let text =
                  value === null
                    ? ""
                    : typeof value === "object"
                      ? canonical(value)
                      : String(value);
                if (/^[=+@-]/.test(text)) text = "'" + text;
                return /[\t\n\r"]/.test(text)
                  ? `"${text.replaceAll('"', '""')}"`
                  : text;
              })
              .join("\t"),
          )
          .join("\n") +
        "\n",
    );
  for (const name of ["human", "mouse"]) {
    const entries = finalEnvelopes.filter(
      (entry) => entry.payload.species === name,
    );
    const targets = entries
      .filter((entry) => entry.kind === "target")
      .map((entry) => entry.payload.data);
    const summaryByRecord = new Map(
      targets.map((target) => [target.summary.record_id, target.summary]),
    );
    const summaries = entries
      .find((entry) => entry.kind === "index")
      .payload.data.map((row) => summaryByRecord.get(row.record_id));
    const constructs = entries
      .filter((entry) => entry.kind === "constructs")
      .flatMap((entry) => entry.payload.data);
    const sequences = entries
      .filter((entry) => entry.kind === "sequence")
      .map((entry) => entry.payload.data);
    const provenance_url = `https://openantigens.org/data/releases/${release_id}/manifest.json`;
    const records_url = `https://openantigens.org/data/releases/${release_id}/${name}/records.jsonl.gz`;
    const contents = {
      "index.jsonl": jsonl(summaries),
      "index.tsv": tsv(
        summaries,
        Object.keys(publicSchema.$defs.summary.properties),
      ),
      "constructs.jsonl": jsonl(constructs),
      "constructs.tsv": tsv(
        constructs,
        Object.keys(publicSchema.$defs.construct.properties),
      ),
      "constructs.fasta": Buffer.from(
        sequences
          .filter((sequence) => sequence.construct_id)
          .map((sequence) =>
            fasta(sequence, { release_id, provenance_url, records_url }),
          )
          .join(""),
      ),
      "reference.fasta": Buffer.from(
        sequences
          .filter((sequence) => sequence.construct_id === null)
          .map((sequence) =>
            fasta(sequence, { release_id, provenance_url, records_url }),
          )
          .join(""),
      ),
      "schema.json": Buffer.from(canonical(publicSchema) + "\n"),
      "citation.json": Buffer.from(canonical(finalManifest.citation) + "\n"),
      "openantigens.bib": Buffer.from("@misc{fixture, title={Fixture}}\n"),
      "openantigens.ris": Buffer.from("TY  - GEN\nTI  - Fixture\nER  - \n"),
      "records.jsonl.gz": gzipSync(jsonl(entries), { mtime: 0 }),
    };
    finalManifest.species[name].downloads = [];
    for (const [filename, [media_type, schema, encoding]] of Object.entries(
      DOWNLOAD_ROLES,
    )) {
      const path = `releases/${release_id}/${name}/downloads/${filename}`;
      const bytes = contents[filename];
      await mkdir(dirname(join(directory, path)), { recursive: true });
      await writeFile(join(directory, path), bytes);
      const file = {
        bytes: bytes.length,
        sha256: sha(bytes),
        media_type,
        schema,
        encoding,
      };
      finalManifest.files[path] = file;
      const download = {
        name: filename,
        url: `https://openantigens.org/data/releases/${release_id}/${name}/${filename}`,
        path,
        ...file,
        row_count:
          filename === "records.jsonl.gz"
            ? entries.length
            : filename.startsWith("index.")
              ? summaries.length
              : filename.startsWith("constructs.")
                ? constructs.length
                : null,
        license: "fixture",
      };
      finalManifest.species[name].downloads.push(download);
      if (filename === "records.jsonl.gz")
        finalManifest.species[name].records_download = download;
    }
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

export async function packedFixture(sourceDirectory) {
  const directory = await mkdtemp(join(tmpdir(), "oa-packed-provider-"));
  const rootBytes = await readFile(join(sourceDirectory, "manifest.json"));
  const root = JSON.parse(rootBytes);
  const transport = {
    format: "openantigens-target-gzip-v1",
    root_manifest_sha256: sha(rootBytes),
    targets: {},
  };
  const write = async (path, bytes) => {
    await mkdir(dirname(join(directory, path)), { recursive: true });
    await writeFile(join(directory, path), bytes);
  };
  await write("manifest.json", rootBytes);
  await write(
    "manifest.sha256",
    Buffer.from(`${sha(rootBytes)}  manifest.json\n`),
  );
  for (const pointer of root.releases) {
    const bytes = await readFile(join(sourceDirectory, pointer.manifest_path));
    const manifest = JSON.parse(bytes);
    await write(pointer.manifest_path, bytes);
    for (const [path, file] of Object.entries(manifest.files)) {
      if (
        file.schema === "index" ||
        [
          "schema.json",
          "citation.json",
          "openantigens.bib",
          "openantigens.ris",
        ].some((name) => path.endsWith("/downloads/" + name))
      )
        await write(path, await readFile(join(sourceDirectory, path)));
      if (file.schema !== "target_shard") continue;
      const target = JSON.parse(
        await readFile(join(sourceDirectory, path)),
      ).data;
      const paths = [
        path,
        ...Object.values(target.sections).flatMap((section) => section.pages),
        ...target.constructs.pages,
        target.reference_sequence_path,
        ...Object.values(target.construct_sequence_paths),
      ]
        .filter(Boolean)
        .sort();
      const members = [];
      for (const member of paths)
        members.push([
          member,
          (await readFile(join(sourceDirectory, member))).toString("base64"),
        ]);
      const decoded = Buffer.from(
        canonical({ format: transport.format, members }) + "\n",
      );
      const compressed = gzipSync(decoded, { level: 9, mtime: 0 });
      const [, release_id, species, , leaf] = path.split("/");
      const pack = `packs/${release_id}/${species}/${leaf.slice(0, -5)}-0.json.gz`;
      await write(pack, compressed);
      transport.targets[path] = [
        {
          path: pack,
          sha256: sha(compressed),
          bytes: compressed.length,
          decoded_bytes: decoded.length,
          members: members.length,
          min_path: members[0][0],
          max_path: members.at(-1)[0],
        },
      ];
    }
  }
  const indexBytes = Buffer.from(canonical(transport) + "\n");
  await write("transport-index.json", indexBytes);
  const seal = {
    format: transport.format,
    root_manifest_sha256: sha(rootBytes),
    transport_index_sha256: sha(indexBytes),
  };
  return {
    directory,
    seal,
    pins: {
      expectedRootSha256: seal.root_manifest_sha256,
      expectedTransportSha256: seal.transport_index_sha256,
    },
  };
}

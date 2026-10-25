import { readFile, readdir, lstat } from "node:fs/promises";
import { join } from "node:path";
import { createHash } from "node:crypto";
import { isDeepStrictEqual } from "node:util";
import { z } from "zod";
import { parsePublicJson } from "./json.js";
import { validateRetirementMetadata } from "./public-boundary.js";
import {
  ApiError,
  DOWNLOAD_ROLES,
  LIMITS,
  SECTIONS,
  STATUSES,
} from "./contracts.js";
const FORMAT = "openantigens-own-static-query-v1";
const sha = (bytes) => createHash("sha256").update(bytes).digest("hex");
const fields = [
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
const fail = () => {
  throw new Error("Invalid local query catalog.");
};
const sourceRow = (summary, path) => ({
  ...Object.fromEntries(fields.map((key) => [key, summary[key]])),
  summary_path: path,
});
const targetPath = (rid, species, summary) =>
  `releases/${rid}/${species}/targets/${sha(Buffer.from(JSON.stringify(summary.target_id ?? summary.record_id)))}.json`;
export async function loadQueryCatalog(
  directory,
  seal,
  { readPhysical, readRange, semanticHash },
) {
  const schema = parsePublicJson(
    await readFile(new URL("./query-schema.json", import.meta.url)),
  );
  const validators = Object.fromEntries(
    Object.keys(schema.$defs).map((name) => [
      name,
      z.fromJSONSchema({ ...schema, $ref: `#/$defs/${name}` }),
    ]),
  );
  const check = (name, value) => {
    if (!validators[name].safeParse(value).success) fail();
  };
  check("seal", seal);
  const rootBytes = await readPhysical("query-root.json", LIMITS.shard);
  if (sha(rootBytes) !== seal.query_root_sha256) fail();
  const root = parsePublicJson(rootBytes);
  check("query_root", root);
  validateRetirementMetadata(root);
  if (
    root.format !== FORMAT ||
    root.static_origin !== seal.static_origin ||
    root.source_root_manifest_sha256 !== seal.source_root_manifest_sha256 ||
    root.transport_index.path !== "transport-index.json" ||
    root.transport_index.bytes > 16 * 1024 * 1024
  )
    fail();
  const expected = new Set(["query-root.json"]);
  const read = async (descriptor, max, collect = true) => {
    if (descriptor.bytes > max || expected.has(descriptor.path)) fail();
    expected.add(descriptor.path);
    return readPhysical(descriptor.path, max, descriptor, collect);
  };
  const transportIndex = parsePublicJson(
    await read(root.transport_index, 16 * 1024 * 1024),
  );
  check("transport_index", transportIndex);
  if (transportIndex.root_manifest_sha256 !== seal.source_root_manifest_sha256)
    fail();
  const infos = new Map();
  let indexBytes = 0;
  const allTargets = new Set();
  for (const descriptor of root.releases) {
    const rid = descriptor.release_id;
    if (
      infos.has(rid) ||
      descriptor.path !== `releases/${rid}/query-release.json`
    )
      fail();
    const meta = parsePublicJson(await read(descriptor, LIMITS.shard));
    check("query_release", meta);
    if (
      meta.release_id !== rid ||
      meta.source_manifest.url !==
        `https://openantigens.org/data/releases/${rid}/manifest.json` ||
      meta.provenance.input_hashes.manifest_url !== meta.source_manifest.url ||
      meta.provenance.input_hashes.manifest_sha256 !==
        meta.source_manifest.sha256
    )
      fail();
    if (
      !["human", "mouse"].some(
        (species) =>
          meta.public_schema.path ===
          `releases/${rid}/${species}/downloads/schema.json`,
      )
    )
      fail();
    const publicSchema = parsePublicJson(
      await read(meta.public_schema, LIMITS.shard),
    );
    if (
      semanticHash(publicSchema) !== meta.source_payload_sha256.schema ||
      semanticHash(meta.citation) !== meta.source_payload_sha256.citation
    )
      fail();
    const publicValidators = Object.fromEntries(
      Object.keys(publicSchema.$defs ?? {}).map((name) => [
        name,
        z.fromJSONSchema({ ...publicSchema, $ref: `#/$defs/${name}` }),
      ]),
    );
    if (!publicValidators.citation?.safeParse(meta.citation).success) fail();
    const manifest = {
      ...meta,
      manifest_url: meta.source_manifest.url,
      species: {},
    };
    const info = {
      query: true,
      manifest,
      pointer: { sha256: meta.source_manifest.sha256 },
      provenanceDTO: meta.provenance,
      schemas: publicSchema,
      validators: publicValidators,
      indexes: {},
      rowsByPath: new Map(),
      summaryLookups: {},
    };
    for (const species of ["human", "mouse"]) {
      const item = meta.species[species],
        prefix = `releases/${rid}/${species}`;
      if (
        item.taxon_id !== (species === "human" ? 9606 : 10090) ||
        item.index.path !== `${prefix}/index.json` ||
        item.summary_lines.path !== `${prefix}/downloads/index.jsonl` ||
        item.summary_lines.offsets.path !== `${prefix}/summary-offsets.json` ||
        new Set(item.capabilities.sections).size !==
          item.capabilities.sections.length ||
        item.capabilities.sections.some(
          (section) => !SECTIONS.includes(section),
        ) ||
        (species === "mouse" && item.capabilities.disease)
      )
        fail();
      const names = new Set();
      for (const download of item.downloads) {
        const role = DOWNLOAD_ROLES[download.name];
        if (
          !role ||
          names.has(download.name) ||
          download.path !== `${prefix}/downloads/${download.name}` ||
          download.url !==
            `https://openantigens.org/data/releases/${rid}/${species}/${download.name}` ||
          download.media_type !== role[0] ||
          download.schema !== role[1] ||
          download.encoding !== role[2] ||
          download.bytes >
            ([
              "schema.json",
              "citation.json",
              "openantigens.bib",
              "openantigens.ris",
            ].includes(download.name)
              ? LIMITS.shard
              : LIMITS.externalDownload)
        )
          fail();
        names.add(download.name);
      }
      const summaryDownload = item.downloads.find(
          (download) => download.name === "index.jsonl",
        ),
        schemaDownload = item.downloads.find(
          (download) => download.name === "schema.json",
        );
      if (
        summaryDownload.sha256 !== item.summary_lines.sha256 ||
        summaryDownload.bytes !== item.summary_lines.bytes ||
        schemaDownload.sha256 !== meta.public_schema.sha256 ||
        schemaDownload.bytes !== meta.public_schema.bytes
      )
        fail();
      const envelope = parsePublicJson(await read(item.index, LIMITS.index));
      indexBytes += item.index.bytes;
      if (
        indexBytes > LIMITS.indexes ||
        !publicValidators.index?.safeParse(envelope).success ||
        envelope.release_id !== rid ||
        envelope.species !== species ||
        envelope.schema_version !== 1
      )
        fail();
      const offsets = parsePublicJson(
        await read(item.summary_lines.offsets, LIMITS.index),
      );
      check("offsets", offsets);
      if (
        offsets.index_sha256 !== item.index.sha256 ||
        offsets.jsonl_sha256 !== item.summary_lines.sha256 ||
        offsets.rows.length !== envelope.data.length
      )
        fail();
      await read(item.summary_lines, LIMITS.externalDownload, false);
      let position = 0;
      for (const [offset, length] of offsets.rows.toSorted(
        (a, b) => a[0] - b[0],
      )) {
        if (
          offset !== position ||
          !Number.isSafeInteger(offset + length) ||
          offset + length > item.summary_lines.bytes
        )
          fail();
        position += length;
      }
      if (position !== item.summary_lines.bytes) fail();
      const seen = new Set(),
        rows = [],
        lookups = new Map(),
        paths = {};
      const statuses = {};
      for (const [index, source] of envelope.data.entries()) {
        const path = source.summary_path ?? targetPath(rid, species, source);
        const row = sourceRow(source, path);
        if (
          row.species !== species ||
          row.taxon_id !== item.taxon_id ||
          !STATUSES.includes(row.status) ||
          !/^[a-f0-9]{64}$/.test(row.record_id) ||
          path !== targetPath(rid, species, row) ||
          seen.has(row.record_id) ||
          allTargets.has(path) ||
          !transportIndex.targets[path]
        )
          fail();
        seen.add(row.record_id);
        allTargets.add(path);
        statuses[row.status] = (statuses[row.status] ?? 0) + 1;
        const parts = transportIndex.targets[path];
        let compressed = 0,
          decoded = 0,
          previous = "";
        const match =
          /^releases\/([a-f0-9]{64})\/(human|mouse)\/targets\/([a-f0-9]{64})\.json$/.exec(
            path,
          );
        for (const [number, part] of parts.entries()) {
          if (
            part.path !==
              `packs/${rid}/${species}/${match[3]}-${number}.json.gz` ||
            part.min_path > part.max_path ||
            (previous && previous >= part.min_path)
          )
            fail();
          previous = part.max_path;
          compressed += part.bytes;
          decoded += part.decoded_bytes;
        }
        if (compressed > 4 * 1024 * 1024 || decoded > 32 * 1024 * 1024) fail();
        if (row.target_id) {
          if (
            Object.hasOwn(paths, row.target_id) ||
            !row.target_id.startsWith(`${item.taxon_id}:`)
          )
            fail();
          paths[row.target_id] = path;
        }
        rows.push(row);
        lookups.set(row.record_id, {
          tuple: offsets.rows[index],
          file: item.summary_lines,
        });
        info.rowsByPath.set(path, row);
      }
      if (
        item.counts.targets !== rows.length ||
        item.counts.resolved_targets !== Object.keys(paths).length ||
        !isDeepStrictEqual(item.counts.status, statuses)
      )
        fail();
      info.indexes[species] = rows.toSorted((a, b) =>
        (a.target_id ?? `~${a.record_id}`) < (b.target_id ?? `~${b.record_id}`)
          ? -1
          : (a.target_id ?? `~${a.record_id}`) >
              (b.target_id ?? `~${b.record_id}`)
            ? 1
            : 0,
      );
      info.summaryLookups[species] = lookups;
      manifest.species[species] = {
        ...item,
        index_path: item.index.path,
        target_paths: paths,
        records_download: item.downloads.find(
          (download) => download.name === item.records_download_name,
        ),
      };
    }
    infos.set(rid, info);
  }
  if (
    !infos.has(root.current_release_id) ||
    allTargets.size !== Object.keys(transportIndex.targets).length ||
    Object.keys(transportIndex.targets).some((path) => !allTargets.has(path))
  )
    fail();
  const retired = new Map();
  for (const row of root.retired_releases) {
    if (infos.has(row.release_id) || retired.has(row.release_id)) fail();
    retired.set(row.release_id, row);
  }
  const directories = new Set();
  for (const path of expected) {
    const parts = path.split("/");
    for (let index = 1; index < parts.length; index++)
      directories.add(parts.slice(0, index).join("/"));
  }
  async function inventory(prefix = "") {
    for (const item of await readdir(join(directory, prefix), {
      withFileTypes: true,
    })) {
      const path = prefix ? `${prefix}/${item.name}` : item.name;
      const stat = await lstat(join(directory, path));
      if (stat.isSymbolicLink()) fail();
      if (stat.isDirectory()) {
        if (!directories.has(path)) fail();
        await inventory(path);
      } else if (!stat.isFile() || !expected.delete(path)) fail();
    }
  }
  await inventory();
  if (expected.size) fail();
  async function readSummary(info, row) {
    try {
      const lookup = info.summaryLookups[row.species]?.get(row.record_id);
      if (!lookup) fail();
      const [offset, length, expectedSha] = lookup.tuple;
      const bytes = await readRange(
        lookup.file.path,
        offset,
        length,
        lookup.file,
      );
      if (sha(bytes) !== expectedSha) fail();
      const summary = parsePublicJson(bytes);
      if (
        !info.validators.summary?.safeParse(summary).success ||
        !isDeepStrictEqual(sourceRow(summary, row.summary_path), row)
      )
        fail();
      return summary;
    } catch {
      const error = new ApiError(
        503,
        "unavailable",
        "The local released query data is unavailable.",
      );
      error.localQueryIntegrity = true;
      throw error;
    }
  }
  return { root, infos, retired, transportIndex, readSummary };
}

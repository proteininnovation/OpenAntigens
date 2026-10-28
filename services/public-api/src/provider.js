import { open, lstat, realpath } from "node:fs/promises";
import { constants, createReadStream } from "node:fs";
import { createGunzip } from "node:zlib";
import { resolve, relative, sep } from "node:path";
import { createHash } from "node:crypto";
import { z } from "zod";
import { createStorage } from "./storage.js";
import { loadQueryCatalog } from "./query-catalog.js";
import { createRemoteQueryStorage } from "./query-storage.js";
import { parsePublicJson } from "./json.js";
import {
  validateRetirementMetadata,
  isPublicArchiveUrl,
} from "./public-boundary.js";
import {
  outputSchema,
  provenanceOutputSchema,
  DOWNLOAD_ROLES,
  rootSchema,
  ApiError,
  bad,
  LIMITS,
  SECTIONS,
  parseInput,
} from "./contracts.js";
const parseJson = parsePublicJson;
export const semanticHash = (value) => digest(semantic(value));
const digest = (b) => createHash("sha256").update(b).digest("hex");
const hex = /^[a-f0-9]{64}$/;
const safePath = (p) =>
  typeof p === "string" &&
  p.length < 1024 &&
  !p.includes("\\") &&
  !p.startsWith("/") &&
  p.split("/").every((v) => v && v !== "." && v !== "..");
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
const semanticNumbers = (value) => {
  if (typeof value === "number") {
    if (!Number.isFinite(value)) throw new Error("Non-finite public number");
    const bytes = Buffer.alloc(8);
    bytes.writeDoubleBE(value);
    return { $number: bytes.toString("hex") };
  }
  if (Array.isArray(value)) return value.map(semanticNumbers);
  if (value && typeof value === "object")
    return Object.fromEntries(
      Object.entries(value).map(([key, item]) => [key, semanticNumbers(item)]),
    );
  return value;
};
const semantic = (value) => canonical(semanticNumbers(value));
const searchRow = (summary, path) => ({
  ...Object.fromEntries(
    [
      "target_id",
      "record_id",
      "species",
      "taxon_id",
      "search_text",
      "track",
      "status",
      "structure_state",
      "construct_state",
    ].map((key) => [key, summary[key]]),
  ),
  summary_path: path,
});
const notFound = () => {
  throw new ApiError(
    404,
    "not_found",
    "The requested released record does not exist.",
  );
};
export class Provider {
  constructor(directory, options = {}) {
    this.directory = resolve(directory);
    this.options = options;
    this.ready = false;
    this.releases = new Map();
    this.retired = new Map();
  }
  approvedFiles(manifest) {
    const roles = new Map(),
      cache = new Map();
    const role = (
      species,
      kind,
      schema,
      media_type = "application/json",
      encoding = null,
      bound = LIMITS.shard,
    ) => {
      const key = `${species}:${kind}:${schema}:${media_type}:${encoding}:${bound}`;
      if (!cache.has(key))
        cache.set(key, { species, kind, schema, media_type, encoding, bound });
      return cache.get(key);
    };
    const rid = manifest.release_id;
    for (const [path, file] of Object.entries(manifest.files)) {
      const parts = path.split("/");
      let expected;
      if (
        parts[0] === "releases" &&
        parts[1] === rid &&
        ["human", "mouse"].includes(parts[2])
      ) {
        const species = parts[2];
        if (parts.length === 4 && parts[3] === "index.json")
          expected = role(
            species,
            "storage",
            "index",
            "application/json",
            null,
            LIMITS.index,
          );
        else if (
          parts.length === 5 &&
          parts[3] === "downloads" &&
          Object.hasOwn(DOWNLOAD_ROLES, parts[4])
        ) {
          const [media, schema, encoding] = DOWNLOAD_ROLES[parts[4]];
          expected = role(
            species,
            parts[4] === "schema.json" ? "schema" : "download",
            schema,
            media,
            encoding,
            [
              "schema.json",
              "citation.json",
              "openantigens.bib",
              "openantigens.ris",
            ].includes(parts[4])
              ? LIMITS.shard
              : LIMITS.externalDownload,
          );
        } else if (
          parts.length === 5 &&
          parts[3] === "targets" &&
          /^[a-f0-9]{64}\.json$/.test(parts[4])
        )
          expected = role(
            species,
            "storage",
            "target_shard",
            "application/json",
            null,
            LIMITS.target,
          );
        else if (
          parts.length === 6 &&
          parts[3] === "sections" &&
          hex.test(parts[4])
        ) {
          const match = /^([a-z]+)-(0|[1-9][0-9]*)\.json$/.exec(parts[5]);
          if (match && SECTIONS.includes(match[1]))
            expected = role(species, "storage", `section_${match[1]}`);
        } else if (
          parts.length === 5 &&
          parts[3] === "constructs" &&
          /^[a-f0-9]{64}-(0|[1-9][0-9]*)\.json$/.test(parts[4])
        )
          expected = role(species, "storage", "construct_page");
        else if (
          parts.length === 5 &&
          parts[3] === "sequences" &&
          /^[a-f0-9]{64}(-reference)?\.json$/.test(parts[4])
        )
          expected = role(species, "storage", "sequence_shard");
      }
      if (
        !expected ||
        file.media_type !== expected.media_type ||
        file.schema !== expected.schema ||
        file.encoding !== expected.encoding
      )
        throw new Error("Unapproved artifact role");
      if (Number.isSafeInteger(file.bytes) && file.bytes > expected.bound) {
        const name =
          expected.kind === "storage" ? expected.schema : parts.at(-1);
        throw new Error(
          `Bundle size: role=${name} actual=${file.bytes} limit=${expected.bound}`,
        );
      }
      roles.set(path, expected);
    }
    for (const species of ["human", "mouse"]) {
      const downloads = manifest.species[species].downloads;
      if (
        downloads.length !== Object.keys(DOWNLOAD_ROLES).length ||
        new Set(downloads.map((item) => item.name)).size !== downloads.length
      )
        throw new Error("Download role coverage");
      for (const download of downloads) {
        const [media, schema, encoding] = DOWNLOAD_ROLES[download.name] ?? [];
        const path = `releases/${rid}/${species}/downloads/${download.name}`;
        const file = manifest.files[path];
        if (
          download.path !== path ||
          !file ||
          download.media_type !== media ||
          download.schema !== schema ||
          download.encoding !== encoding ||
          canonical(file) !==
            canonical(
              Object.fromEntries(
                ["bytes", "sha256", "media_type", "schema", "encoding"].map(
                  (key) => [key, download[key]],
                ),
              ),
            )
        )
          throw new Error("Download role disagreement");
      }
    }
    return roles;
  }
  readSafe(path, max, expected, collect = true, targetPath) {
    return this.storage
      ? this.storage.read(path, max, expected, collect, targetPath)
      : this.readPhysical(path, max, expected, collect);
  }
  async openPhysical(path, max, expected) {
    if (!safePath(path)) throw new Error("Invalid bundle path");
    let current = this.directory;
    for (const part of path.split("/")) {
      current = resolve(current, part);
      const stat = await lstat(current);
      if (stat.isSymbolicLink()) throw new Error("Bundle symlink");
    }
    const actual = await realpath(current);
    if (
      relative(this.directory, actual).startsWith(`..${sep}`) ||
      !actual.startsWith(`${this.directory}${sep}`)
    )
      throw new Error("Bundle containment");
    const handle = await open(
      actual,
      constants.O_RDONLY | constants.O_NOFOLLOW,
    );
    try {
      const stat = await handle.stat();
      if (
        !stat.isFile() ||
        stat.size > max ||
        (expected && stat.size !== expected.bytes)
      )
        throw new Error("Bundle size");
      return { handle, stat };
    } catch (error) {
      await handle.close();
      throw error;
    }
  }
  async readPhysical(path, max, expected, collect = true) {
    const { handle, stat } = await this.openPhysical(path, max, expected);
    try {
      const collected = collect ? Buffer.alloc(stat.size) : null;
      const checksum = createHash("sha256");
      let length = 0;
      for (;;) {
        const chunk =
          collected && length < stat.size
            ? collected.subarray(length, Math.min(length + 65536, stat.size))
            : Buffer.alloc(collect ? 1 : Math.min(65536, max - length + 1));
        const { bytesRead } = await handle.read(chunk);
        if (!bytesRead) break;
        length += bytesRead;
        if (
          length > stat.size ||
          length > max ||
          (expected && length > expected.bytes)
        )
          throw new Error("Bundle size");
        const bytes = chunk.subarray(0, bytesRead);
        checksum.update(bytes);
      }
      if (
        length !== stat.size ||
        (expected &&
          (length !== expected.bytes ||
            checksum.digest("hex") !== expected.sha256))
      )
        throw new Error("Bundle checksum");
      return collected;
    } finally {
      await handle.close();
    }
  }
  async readPhysicalRange(path, offset, length, expected) {
    if (
      !Number.isSafeInteger(offset) ||
      offset < 0 ||
      !Number.isSafeInteger(length) ||
      length < 1 ||
      length > LIMITS.shard ||
      !Number.isSafeInteger(offset + length) ||
      offset + length > expected.bytes
    )
      throw new Error("Invalid local summary range");
    const { handle, stat } = await this.openPhysical(
      path,
      LIMITS.externalDownload,
      expected,
    );
    try {
      const bytes = Buffer.alloc(length);
      let read = 0;
      while (read < length) {
        const result = await handle.read(
          bytes,
          read,
          length - read,
          offset + read,
        );
        if (!result.bytesRead) throw new Error("Invalid local summary range");
        read += result.bytesRead;
      }
      if ((await handle.stat()).size !== stat.size)
        throw new Error("Invalid local summary range");
      return bytes;
    } finally {
      await handle.close();
    }
  }
  async initializeQuery() {
    const catalog = await loadQueryCatalog(
      this.directory,
      this.options.querySeal,
      {
        readPhysical: this.readPhysical.bind(this),
        readRange: this.readPhysicalRange.bind(this),
        semanticHash,
      },
    );
    this.queryCatalog = catalog;
    this.releases = catalog.infos;
    this.retired = catalog.retired;
    this.current = catalog.root.current_release_id;
    this.storage = createRemoteQueryStorage({
      transportIndex: catalog.transportIndex,
      releases: this.releases,
      readSummary: catalog.readSummary,
      fetchImpl: this.options.fetchImpl,
    });
    const current = this.releases.get(this.current);
    const metadata = this.envelope(
      current,
      null,
      this.releaseMetadata(current),
      {
        discovery_state: "active",
        current_release_id: this.current,
        releases: [...this.releases.values()].map((info) =>
          this.releaseMetadata(info),
        ),
        retirement_history: {
          known_retired_count: this.retired.size,
          limit: 1000,
          lookup_endpoint: "/v1/releases/{release_id}",
        },
      },
    );
    if (Buffer.byteLength(JSON.stringify(metadata)) > LIMITS.json)
      throw new Error("Local release metadata exceeds response limit");
    this.ready = true;
    return this;
  }
  async initialize() {
    this.ready = false;
    this.directory = await realpath(this.directory);
    if (this.options.querySeal) return this.initializeQuery();
    this.storage = await createStorage(
      this.directory,
      this.options,
      this.readPhysical.bind(this),
    );
    return this.storage.withContext(() => this.initializeData());
  }
  async readReleaseManifest(pointer) {
    const bytes = await this.readSafe(pointer.manifest_path, LIMITS.manifest);
    if (digest(bytes) !== pointer.sha256) throw new Error("Release checksum");
    return parseJson(bytes);
  }
  async initializeData() {
    this.ready = false;
    this.releases.clear();
    this.retired.clear();
    this.directory = await realpath(this.directory);
    const root = parseJson(await this.readSafe("manifest.json", LIMITS.shard));
    validateRetirementMetadata(root);
    if (!rootSchema.safeParse(root).success) throw new Error("Bundle manifest");
    const retainedIds = new Set(
      root.releases.map((pointer) => pointer.release_id),
    );
    for (const tombstone of root.retired_releases ?? []) {
      if (
        retainedIds.has(tombstone.release_id) ||
        this.retired.has(tombstone.release_id)
      )
        throw new Error("Retired release identity");
      this.retired.set(tombstone.release_id, tombstone);
    }
    let totalIndexBytes = 0;
    for (const pointer of root.releases) {
      if (
        !hex.test(pointer.release_id) ||
        !hex.test(pointer.sha256) ||
        this.releases.has(pointer.release_id)
      )
        throw new Error("Release pointer");
      const manifest = await this.readReleaseManifest(pointer);
      if (
        manifest.schema_version !== 1 ||
        manifest.release_id !== pointer.release_id ||
        !manifest.files ||
        !manifest.species
      )
        throw new Error("Release manifest");
      if (
        manifest.core &&
        digest(canonical(manifest.core)) !== pointer.release_id
      )
        throw new Error("Release identity");
      if (
        !manifest.files ||
        typeof manifest.files !== "object" ||
        Array.isArray(manifest.files)
      )
        throw new Error("File descriptor map");
      const fileRoles = this.approvedFiles(manifest);
      this.storage.registerManifest(manifest);
      const schemaPaths = Object.keys(manifest.files).filter(
        (path) => manifest.files[path].media_type === "application/schema+json",
      );
      if (!schemaPaths.length) throw new Error("Missing release schema");
      const schemas = parseJson(
        await this.readSafe(
          schemaPaths[0],
          LIMITS.shard,
          manifest.files[schemaPaths[0]],
        ),
      );
      const validators = Object.fromEntries(
        Object.keys(schemas.$defs ?? {}).map((name) => [
          name,
          z.fromJSONSchema({ ...schemas, $ref: `#/$defs/${name}` }),
        ]),
      );
      if (
        !validators.release_manifest?.safeParse({ ...manifest, files: {} })
          .success
      )
        throw new Error("Manifest schema");
      const info = {
        manifest,
        pointer,
        indexes: {},
        schemas,
        validators,
        sourceIndexes: {},
        targetRecords: new Map(),
        envelopeDigests: new Map(),
        fileRoles,
        reachable: new Set(),
        legacyIndexDigests: {},
      };
      const processStorage = (path, content) => {
        const file = manifest.files[path];
        if (!file || fileRoles.get(path)?.kind !== "storage")
          throw new Error("Unapproved query member");
        const envelope = parseJson(content);
        if (!this.storage.packed)
          info.envelopeDigests.set(path, {
            sha256: semanticHash(envelope),
            species: envelope.species,
            kind:
              file.schema === "index"
                ? "index"
                : file.schema === "target_shard"
                  ? "target"
                  : file.schema === "construct_page"
                    ? "constructs"
                    : file.schema === "sequence_shard"
                      ? "sequence"
                      : "section",
            bytes: file.bytes,
          });
        if (file.schema === "target_shard") {
          const key = `${envelope.species}:${envelope.data.summary.record_id}`;
          if (info.targetRecords.has(key))
            throw new Error("Duplicate target record");
          info.targetRecords.set(key, path);
        }
        if (!validators[file.schema]?.safeParse(envelope).success)
          throw new Error("Shard schema");
        if (
          envelope.schema_version !== 1 ||
          envelope.release_id !== manifest.release_id ||
          !["human", "mouse"].includes(envelope.species) ||
          !Object.hasOwn(envelope, "data")
        )
          throw new Error("Data envelope");
      };
      for (const [path, file] of Object.entries(manifest.files)) {
        if (!validators.file?.safeParse(file).success)
          throw new Error("File descriptor");
        const role = fileRoles.get(path);
        if (
          file.media_type === "application/schema+json" &&
          file.sha256 !== manifest.files[schemaPaths[0]].sha256
        )
          throw new Error("Schema disagreement");
        if (this.storage.isExternal(path) || this.storage.isVirtual(path))
          continue;
        const jsonBoundary = [
          "application/json",
          "application/schema+json",
        ].includes(role.media_type);
        const content = await this.readSafe(
          path,
          role.bound,
          file,
          role.kind === "storage" || jsonBoundary,
        );
        if (role.kind === "storage") processStorage(path, content);
        else if (jsonBoundary) {
          const document = parsePublicJson(content);
          if (
            file.schema &&
            (!validators[file.schema]?.safeParse(document).success ||
              (file.schema === "citation" &&
                semanticHash(document) !== semanticHash(manifest.citation)))
          )
            throw new Error("Public JSON resource identity");
        }
      }
      if (this.storage.packed)
        for await (const { path, bytes } of this.storage.scanRelease(manifest))
          processStorage(path, bytes);
      const downloadOrigins = new Set(
        Object.values(manifest.species).flatMap((s) =>
          s.downloads.map((d) => new URL(d.url).origin),
        ),
      );
      if (downloadOrigins.size !== 1)
        throw new Error("Download origin disagreement");
      if (manifest.manifest_url) {
        const manifestUrl = new URL(manifest.manifest_url);
        if (
          manifestUrl.href !==
          `${[...downloadOrigins][0]}/data/releases/${manifest.release_id}/manifest.json`
        )
          throw new Error("Manifest URL identity");
      }
      for (const species of ["human", "mouse"]) {
        const entry = manifest.species[species];
        if (
          !entry ||
          entry.taxon_id !== (species === "human" ? 9606 : 10090) ||
          !entry.target_paths ||
          !entry.capabilities ||
          !Array.isArray(entry.capabilities.sections) ||
          entry.capabilities.sections.some((s) => !SECTIONS.includes(s))
        )
          throw new Error("Species manifest");
        const rawIndex = await this.json(info, entry.index_path, LIMITS.index);
        const version = schemas["x-openantigens-index-version"] ?? 1;
        if (
          ![1, 2].includes(version) ||
          entry.index_path !==
            `releases/${manifest.release_id}/${species}/index.json` ||
          !Array.isArray(rawIndex)
        )
          throw new Error("Index schema");
        totalIndexBytes += manifest.files[entry.index_path].bytes;
        const digests = new Map();
        const hadDiseaseIds =
          version === 2 &&
          Object.hasOwn(
            schemas.$defs.search_row?.properties ?? {},
            "disease_ids",
          );
        const index = rawIndex.map((row) => {
          const path =
            version === 2
              ? row.summary_path
              : info.targetRecords.get(`${species}:${row.record_id}`);
          if (
            !path ||
            info.targetRecords.get(`${species}:${row.record_id}`) !== path ||
            row.species !== species ||
            row.taxon_id !== entry.taxon_id ||
            !hex.test(row.record_id) ||
            (row.target_id && entry.target_paths[row.target_id] !== path)
          )
            throw new Error("Index identity");
          if (version === 1 || hadDiseaseIds)
            digests.set(row.record_id, semanticHash(row));
          return searchRow(row, path);
        });
        if (new Set(index.map((row) => row.record_id)).size !== index.length)
          throw new Error("Duplicate index record");
        if (
          Object.keys(entry.target_paths).length !==
          index.filter((row) => row.target_id).length
        )
          throw new Error("Resolved target catalogue disagreement");
        info.legacyIndexDigests[species] =
          version === 1 || hadDiseaseIds
            ? { kind: version === 1 ? "full" : "compact_disease", digests }
            : null;
        info.sourceIndexes[species] = index;
        info.indexes[species] = index.toSorted((a, b) =>
          (a.target_id ?? `~${a.record_id}`) <
          (b.target_id ?? `~${b.record_id}`)
            ? -1
            : (a.target_id ?? `~${a.record_id}`) >
                (b.target_id ?? `~${b.record_id}`)
              ? 1
              : 0,
        );
        if (
          !entry.records_download ||
          !entry.downloads?.some((d) => d.name === "records.jsonl.gz")
        )
          throw new Error("Missing complete fallback");
        for (const download of [...entry.downloads, entry.records_download]) {
          const file = manifest.files[download.path];
          const url = new URL(download.url);
          if (
            url.protocol !== "https:" ||
            url.username ||
            url.password ||
            download.url !==
              `${url.origin}/data/releases/${manifest.release_id}/${species}/${download.name}` ||
            !file ||
            file.sha256 !== download.sha256 ||
            file.bytes !== download.bytes
          )
            throw new Error("Download manifest");
        }
        info.reachable.add(entry.index_path);
        for (const download of entry.downloads)
          info.reachable.add(download.path);
        for (const row of index) {
          const key = digest(canonical(row.target_id ?? row.record_id));
          const prefix = `releases/${manifest.release_id}/${species}`;
          if (row.summary_path !== `${prefix}/targets/${key}.json`)
            throw new Error("Target path identity");
          const target = await this.json(info, row.summary_path, LIMITS.target);
          if (
            target.summary.target_id !== row.target_id ||
            target.summary.record_id !== row.record_id ||
            target.summary.species !== species
          )
            throw new Error("Target identity");
          const add = (path, schema) => {
            if (
              !path ||
              fileRoles.get(path)?.schema !== schema ||
              fileRoles.get(path)?.species !== species
            )
              throw new Error("Target shard role");
            info.reachable.add(path);
          };
          add(row.summary_path, "target_shard");
          for (const [section, descriptor] of Object.entries(target.sections))
            for (const [page, path] of descriptor.pages.entries()) {
              if (path !== `${prefix}/sections/${key}/${section}-${page}.json`)
                throw new Error("Section path identity");
              add(path, `section_${section}`);
            }
          for (const [page, path] of target.constructs.pages.entries()) {
            if (path !== `${prefix}/constructs/${key}-${page}.json`)
              throw new Error("Construct path identity");
            add(path, "construct_page");
          }
          if (target.reference_sequence_path) {
            if (
              target.reference_sequence_path !==
              `${prefix}/sequences/${key}-reference.json`
            )
              throw new Error("Reference path identity");
            add(target.reference_sequence_path, "sequence_shard");
          }
          for (const [cid, path] of Object.entries(
            target.construct_sequence_paths,
          ))
            if (path) {
              if (path !== `${prefix}/sequences/${cid}.json`)
                throw new Error("Sequence path identity");
              add(path, "sequence_shard");
            }
        }
      }
      if (
        info.reachable.size !== fileRoles.size ||
        [...fileRoles.keys()].some((path) => !info.reachable.has(path))
      )
        throw new Error("Unreachable artifact file");
      if (totalIndexBytes > LIMITS.indexes)
        throw new Error("Index memory budget");
      await this.verifySemanticIdentity(info);
      if (!this.storage.packed)
        for (const species of ["human", "mouse"])
          await this.verifyRecordsDownload(info, species);
      delete info.envelopeDigests;
      delete info.sourceIndexes;
      delete info.targetRecords;
      delete info.fileRoles;
      delete info.reachable;
      delete info.legacyIndexDigests;
      this.releases.set(manifest.release_id, info);
    }
    if (!this.releases.has(root.current_release_id))
      throw new Error("Current release");
    await this.storage.validateInventory(root);
    this.current = root.current_release_id;
    this.ready = true;
    return this;
  }
  async verifyRecordsDownload(info, species) {
    const sm = info.manifest.species[species];
    const download = sm.records_download;
    if (
      canonical(download) !==
      canonical(sm.downloads.find((d) => d.name === "records.jsonl.gz"))
    )
      throw new Error("Fallback download disagreement");
    const expected = new Map(
      [...info.envelopeDigests].filter(
        ([, entry]) => entry.species === species,
      ),
    );
    const maximumBytes = [...expected.values()].reduce(
      (sum, entry) => sum + entry.bytes + 2048,
      0,
    );
    const maximumLine =
      [...expected.values()].reduce(
        (max, entry) => Math.max(max, entry.bytes),
        0,
      ) + 2048;
    const source = createReadStream(resolve(this.directory, download.path), {
      flags: constants.O_RDONLY | constants.O_NOFOLLOW,
      highWaterMark: 65536,
    });
    const gunzip = createGunzip();
    source.on("error", (error) => gunzip.destroy(error));
    source.pipe(gunzip);
    const decoder = new TextDecoder("utf-8", { fatal: true });
    let pending = "",
      bytes = 0,
      count = 0;
    const seen = new Set();
    try {
      for await (const chunk of gunzip) {
        bytes += chunk.length;
        if (bytes > maximumBytes)
          throw new Error("Fallback decompression budget");
        pending += decoder.decode(chunk, { stream: true });
        let newline;
        while ((newline = pending.indexOf("\n")) !== -1) {
          const line = pending.slice(0, newline);
          pending = pending.slice(newline + 1);
          if (Buffer.byteLength(line) > maximumLine)
            throw new Error("Fallback record budget");
          const record = parsePublicJson(line);
          const entry = expected.get(record.path);
          if (
            !entry ||
            seen.has(record.path) ||
            Object.keys(record).sort().join(",") !== "kind,path,payload" ||
            record.kind !== entry.kind ||
            semanticHash(record.payload) !== entry.sha256
          )
            throw new Error("Fallback record disagreement");
          seen.add(record.path);
          count++;
        }
        if (Buffer.byteLength(pending) > maximumLine)
          throw new Error("Fallback record budget");
      }
      pending += decoder.decode();
      if (pending || count !== expected.size || count !== download.row_count)
        throw new Error("Fallback coverage disagreement");
    } finally {
      source.destroy();
      gunzip.destroy();
    }
  }
  async verifySemanticIdentity(info) {
    const { manifest } = info;
    const payload = manifest.core.payload_sha256;
    if (
      canonical(manifest.provenance) !== canonical(manifest.core.provenance) ||
      digest(semantic(info.schemas)) !== payload.schema ||
      digest(semantic(manifest.citation)) !== payload.citation
    )
      throw new Error("Release metadata identity");
    const assets = {};
    for (const species of ["human", "mouse"]) {
      const checksum = createHash("sha256");
      const write = (value) => checksum.update(value);
      const indexes = info.sourceIndexes[species];
      const sm = manifest.species[species];
      const status = {};
      let constructTotal = 0;
      write("[");
      for (let position = 0; position < indexes.length; position++) {
        const row = indexes[position];
        const path = row.summary_path;
        const target = await this.json(info, path, LIMITS.target);
        const summary = target.summary;
        if (
          canonical(searchRow(summary, path)) !== canonical(row) ||
          (info.legacyIndexDigests[species] &&
            semanticHash(
              info.legacyIndexDigests[species].kind === "full"
                ? summary
                : {
                    ...searchRow(summary, path),
                    disease_ids: summary.disease_ids,
                  },
            ) !== info.legacyIndexDigests[species].digests.get(row.record_id))
        )
          throw new Error("Index target disagreement");
        status[summary.status] = (status[summary.status] ?? 0) + 1;
        if (position) write(",");
        write('{"construct_sequences":{');
        let first = true;
        for (const cid of Object.keys(target.construct_sequence_paths).sort()) {
          if (!first) write(",");
          first = false;
          const sequencePath = target.construct_sequence_paths[cid];
          const sequence = sequencePath
            ? await this.json(info, sequencePath, LIMITS.shard, path)
            : null;
          if (
            sequence &&
            (sequence.construct_id !== cid ||
              sequence.target_id !== summary.target_id)
          )
            throw new Error("Construct sequence identity");
          write(`${canonical(cid)}:${semantic(sequence)}`);
        }
        write('},"constructs":[');
        let constructCount = 0;
        for (const page of target.constructs.pages) {
          for (const construct of await this.json(info, page)) {
            if (construct.target_id !== summary.target_id)
              throw new Error("Construct target identity");
            if (constructCount++) write(",");
            write(semantic(construct));
          }
        }
        if (
          constructCount !== target.constructs.total ||
          (summary.construct_count !== null &&
            constructCount !== summary.construct_count)
        )
          throw new Error("Construct count disagreement");
        constructTotal += constructCount;
        write('],"reference_sequence":');
        const reference = target.reference_sequence_path
          ? await this.json(info, target.reference_sequence_path)
          : null;
        if (
          reference &&
          (reference.target_id !== summary.target_id ||
            reference.construct_id !== null)
        )
          throw new Error("Reference sequence identity");
        write(semantic(reference));
        write(',"sections":{');
        first = true;
        for (const section of Object.keys(target.sections).sort()) {
          if (!first) write(",");
          first = false;
          const descriptor = target.sections[section];
          const coverage =
            section === "disease" ? summary.disease_coverage : undefined;
          write(
            `${canonical(section)}:{${coverage ? `"coverage":${semantic(coverage)},` : ""}"items":[`,
          );
          let count = 0;
          for (const page of descriptor.pages) {
            const evidence = await this.json(info, page);
            if (
              evidence.target_id !== summary.target_id ||
              evidence.section !== section ||
              evidence.state !== descriptor.state ||
              evidence.reason !== descriptor.reason
            )
              throw new Error("Evidence identity");
            if (
              coverage &&
              canonical(evidence.coverage) !== canonical(coverage)
            )
              throw new Error("Disease coverage disagreement");
            for (const item of evidence.items) {
              if (count++) write(",");
              write(semantic(item));
              if (
                section === "structures" &&
                item.sha256 &&
                item.live_asset_url
              )
                assets[new URL(item.live_asset_url).pathname.slice(1)] =
                  item.sha256;
            }
          }
          if (count !== descriptor.total)
            throw new Error("Evidence count disagreement");
          write(
            `],"reason":${semantic(descriptor.reason)},"state":${semantic(descriptor.state)}}`,
          );
        }
        write(`},"summary":${semantic(summary)}}`);
      }
      write("]");
      if (checksum.digest("hex") !== payload[species])
        throw new Error("Public payload identity");
      if (
        sm.counts.targets !== indexes.length ||
        sm.counts.resolved_targets !==
          indexes.filter((row) => row.target_id).length ||
        sm.counts.constructs !== constructTotal ||
        canonical(sm.counts.status) !== canonical(status)
      )
        throw new Error("Release count disagreement");
    }
    if (canonical(assets) !== canonical(manifest.core.asset_sha256))
      throw new Error("Asset payload identity");
  }
  async json(info, path, max = LIMITS.shard, targetPath) {
    if (info.query) {
      try {
        return await this.storage.readJson(info, path, targetPath);
      } catch (error) {
        if (error.localQueryIntegrity) this.ready = false;
        if (error instanceof ApiError) throw error;
        throw new ApiError(
          503,
          "unavailable",
          "The requested released data is temporarily unavailable.",
        );
      }
    }
    try {
      const file = info.manifest.files[path];
      if (!file) throw new Error("Unlisted shard");
      const result = parseJson(
        await this.readSafe(path, max, file, true, targetPath),
      );
      return result.data;
    } catch (error) {
      this.ready = false;
      throw error;
    }
  }
  resolve(id) {
    if (!this.ready)
      throw new ApiError(
        503,
        "unavailable",
        "The validated release bundle is unavailable.",
      );
    const info = this.releases.get(id);
    if (!info) {
      const retirement = this.retired.get(id);
      if (retirement && !isPublicArchiveUrl(retirement.archive_url)) {
        this.ready = false;
        throw new ApiError(
          503,
          "unavailable",
          "The released metadata is unavailable.",
        );
      }
      if (retirement)
        throw new ApiError(
          410,
          "retired_release",
          "This release is retired. Its recorded archive availability is included.",
          { release_id: id, retirement },
        );
      notFound();
    }
    return info;
  }
  manifestReference(info) {
    const m = info.manifest;
    return {
      url:
        m.manifest_url ??
        `${new URL(m.species.human.records_download.url).origin}/data/releases/${m.release_id}/manifest.json`,
      sha256: info.pointer.sha256,
    };
  }
  outputProvenance(info) {
    const {
      software_commit,
      snapshot_id,
      source_versions,
      source_dates,
      notes,
      input_sha256,
      projection_methods,
    } = info.provenanceDTO ?? info.manifest.provenance;
    const manifest = this.manifestReference(info);
    return {
      software_commit,
      snapshot_id,
      source_versions,
      source_dates,
      notes,
      projection_methods: projection_methods ?? null,
      input_hashes: {
        count: info.provenanceDTO
          ? info.provenanceDTO.input_hashes.count
          : Object.keys(input_sha256).length,
        manifest_url: manifest.url,
        manifest_sha256: manifest.sha256,
        json_pointer: "/provenance/input_sha256",
      },
    };
  }
  envelope(info, species, data, extra = {}) {
    const m = info.manifest;
    return {
      release_id: m.release_id,
      schema_version: 1,
      ...(species ? { species, taxon_id: m.species[species].taxon_id } : {}),
      data,
      manifest: this.manifestReference(info),
      provenance: this.outputProvenance(info),
      ...extra,
    };
  }
  releaseMetadata(info) {
    const m = info.manifest;
    return {
      release_id: m.release_id,
      schema_version: 1,
      provenance: this.outputProvenance(info),
      citation: m.citation,
      licenses: m.licenses,
      manifest: { sha256: info.pointer.sha256 },
      species: Object.fromEntries(
        Object.entries(m.species).map(([s, v]) => [
          s,
          {
            taxon_id: v.taxon_id,
            counts: v.counts,
            capabilities: v.capabilities,
            references: v.downloads.filter((download) =>
              [
                "schema.json",
                "citation.json",
                "openantigens.bib",
                "openantigens.ris",
              ].includes(download.name),
            ),
          },
        ]),
      ),
    };
  }
  call(name, raw, signal) {
    return this.storage
      ? this.storage.withContext(() => this.callData(name, raw, signal), {
          signal,
        })
      : this.callData(name, raw, signal);
  }
  async callData(name, raw, signal) {
    signal?.throwIfAborted();
    const args = parseInput(name, raw);
    if (
      name === "get_release" &&
      !args.release_id &&
      this.publicationState === "maintenance"
    )
      throw new ApiError(
        503,
        "publication_in_progress",
        "Current release discovery is temporarily unavailable during publication.",
      );
    const info = this.resolve(args.release_id ?? this.current);
    const species = args.species;
    if (name === "get_release")
      return this.envelope(
        info,
        null,
        this.releaseMetadata(info),
        this.publicationState === "maintenance"
          ? { discovery_state: "publication_in_progress" }
          : {
              discovery_state: "active",
              current_release_id: this.current,
              releases: [...this.releases.values()].map((i) =>
                this.releaseMetadata(i),
              ),
              retirement_history: {
                known_retired_count: this.retired.size,
                limit: 1000,
                lookup_endpoint: "/v1/releases/{release_id}",
              },
            },
      );
    const sm = info.manifest.species[species];
    if (name === "search_targets") {
      const q = (args.q ?? "").toLocaleLowerCase("en").trim();
      let rows = info.indexes[species].filter(
        (r) =>
          (!q || r.search_text?.toLocaleLowerCase("en").includes(q)) &&
          (!args.track || r.track === args.track) &&
          (!args.status || r.status === args.status) &&
          (!args.structure || r.structure_state === args.structure) &&
          (!args.constructs || r.construct_state === args.constructs),
      );
      const deadline = Date.now() + LIMITS.handlerMs;
      const pagination = this.pagination(info, species, args, "targets");
      if (pagination.position >= LIMITS.searchWindow)
        bad("refine_query", "Refine the query to inspect further matches.");
      const page = this.page(
        info,
        species,
        rows.slice(0, LIMITS.searchWindow),
        args,
        "targets",
        {
          search_window: LIMITS.searchWindow,
          refine_query: rows.length > LIMITS.searchWindow,
        },
      );
      page.total = rows.length;
      const selected = page.data;
      page.data = [];
      for (const row of selected) {
        signal?.throwIfAborted();
        if (Date.now() > deadline)
          throw new ApiError(
            503,
            "handler_timeout",
            "The request exceeded its time budget.",
          );
        let summary;
        if (info.query) {
          try {
            summary = await this.queryCatalog.readSummary(info, row);
          } catch (error) {
            this.ready = false;
            throw error;
          }
        } else
          summary = (await this.json(info, row.summary_path, LIMITS.target))
            .summary;
        if (
          summary.target_id !== row.target_id ||
          summary.record_id !== row.record_id
        )
          throw new Error("Target identity");
        page.data.push(summary);
      }
      return page;
    }
    if (!args.target_id.startsWith(`${sm.taxon_id}:`))
      bad(
        "species_mismatch",
        "The target ID does not match the requested species.",
      );
    const targetPath = sm.target_paths[args.target_id];
    if (!targetPath) notFound();
    const target = await this.json(info, targetPath, LIMITS.target);
    if (name === "get_target") {
      if (!args.section)
        return this.envelope(info, species, {
          ...target.summary,
          sections: Object.fromEntries(
            Object.entries(target.sections).map(([k, v]) => [
              k,
              { state: v.state, reason: v.reason, total: v.total },
            ]),
          ),
        });
      if (!sm.capabilities.sections.includes(args.section))
        throw new ApiError(
          400,
          "unsupported_capability",
          "This evidence section is not supported for this species and release.",
        );
      const descriptor = target.sections[args.section];
      if (!descriptor) throw new Error("Missing section descriptor");
      return this.shardPage(
        info,
        species,
        descriptor.pages,
        args,
        `section:${args.target_id}:${args.section}`,
        (data) => data.items,
        () => true,
        {
          state: descriptor.state,
          reason: descriptor.reason,
          ...(args.section === "disease" && target.summary.disease_coverage
            ? { coverage: target.summary.disease_coverage }
            : {}),
        },
        LIMITS.maxPage,
        signal,
      );
    }
    if (name === "list_constructs") {
      return this.shardPage(
        info,
        species,
        target.constructs.pages,
        args,
        `constructs:${args.target_id}`,
        (data) => data,
        (r) => !args.category || r.category === args.category,
        { state: target.constructs.state, reason: target.constructs.reason },
        LIMITS.maxPage,
        signal,
      );
    }
    if (name === "get_sequence") {
      const path = args.construct_id
        ? target.construct_sequence_paths[args.construct_id]
        : target.reference_sequence_path;
      if (!path) notFound();
      const sequence = await this.json(info, path, LIMITS.shard, targetPath);
      if (
        sequence.target_id !== args.target_id ||
        sequence.construct_id !== (args.construct_id ?? null) ||
        sequence.length !== sequence.sequence.length ||
        digest(sequence.sequence) !== sequence.sequence_sha256
      )
        throw new Error("Sequence identity");
      return this.envelope(info, species, sequence);
    }
    throw new Error("Unknown operation");
  }
  pagination(info, species, args, scope, maximum = LIMITS.maxPage) {
    const limit = args.limit ?? LIMITS.defaultPage;
    if (limit > maximum) bad();
    const normalized = { ...args };
    delete normalized.cursor;
    delete normalized.limit;
    if (scope === "targets")
      normalized.q = (args.q ?? "").toLocaleLowerCase("en").trim();
    const fingerprint = digest(canonical(normalized));
    let position = 0;
    if (args.cursor) {
      try {
        const c = parsePublicJson(Buffer.from(args.cursor, "base64url"));
        if (
          c.v !== 1 ||
          c.release_id !== info.manifest.release_id ||
          c.species !== species ||
          c.scope !== scope ||
          c.filter !== fingerprint ||
          c.sort !== "stable-id-v1" ||
          !Number.isSafeInteger(c.position) ||
          c.position < 1 ||
          Object.keys(c).length !== 7
        )
          bad("invalid_cursor");
        position = c.position;
      } catch {
        bad("invalid_cursor");
      }
    }
    return {
      limit,
      position,
      finish: (data, total, extra) => {
        if (position >= total && position !== 0) bad("invalid_cursor");
        const next = position + limit;
        const cursor =
          next < total
            ? Buffer.from(
                JSON.stringify({
                  v: 1,
                  release_id: info.manifest.release_id,
                  species,
                  scope,
                  filter: fingerprint,
                  sort: "stable-id-v1",
                  position: next,
                }),
              ).toString("base64url")
            : null;
        return this.envelope(info, species, data, {
          next_cursor: cursor,
          total,
          ...extra,
        });
      },
    };
  }
  page(info, species, rows, args, scope, extra = {}, maximum = LIMITS.maxPage) {
    const p = this.pagination(info, species, args, scope, maximum);
    return p.finish(
      rows.slice(p.position, p.position + p.limit),
      rows.length,
      extra,
    );
  }
  async shardPage(
    info,
    species,
    paths,
    args,
    scope,
    records,
    include,
    extra = {},
    maximum = LIMITS.maxPage,
    signal,
  ) {
    const p = this.pagination(info, species, args, scope, maximum);
    const selected = [];
    let total = 0;
    const deadline = Date.now() + LIMITS.handlerMs;
    for (const path of paths) {
      signal?.throwIfAborted();
      if (Date.now() > deadline)
        throw new ApiError(
          503,
          "handler_timeout",
          "The request exceeded its time budget.",
        );
      const rows = records(await this.json(info, path));
      if (!Array.isArray(rows) || rows.length > 25)
        throw new Error("Invalid collection shard");
      for (const row of rows)
        if (include(row)) {
          if (total >= p.position && selected.length < p.limit)
            selected.push(row);
          total++;
        }
    }
    return p.finish(selected, total, extra);
  }
  toolOutput(name) {
    const schemas = [...this.releases.values()].map((info) =>
      this.releaseToolOutput(name, info),
    );
    return schemas.length > 1 ? z.union(schemas) : (schemas[0] ?? outputSchema);
  }
  releaseToolOutput(name, info) {
    const defs = info.schemas.$defs;
    const ref = (name) => ({ $ref: `#/$defs/${name}` });
    let data;
    if (name === "search_targets")
      data = { type: "array", items: ref("summary") };
    else if (name === "list_constructs")
      data = { type: "array", items: ref("construct") };
    else if (name === "get_sequence") data = ref("sequence");
    else if (name === "get_target") {
      const summary = structuredClone(defs.summary);
      summary.properties.sections = {
        type: "object",
        additionalProperties: {
          type: "object",
          properties: {
            state: { type: "string" },
            reason: { type: ["string", "null"] },
            total: { type: "integer", minimum: 0 },
          },
          required: ["state", "reason", "total"],
          additionalProperties: false,
        },
      };
      summary.required.push("sections");
      data = {
        anyOf: [
          summary,
          ...SECTIONS.filter((s) => defs[`section_${s}`]).map((s) => ({
            type: "array",
            items: defs[`section_${s}`].properties.data.properties.items.items,
          })),
        ],
      };
    } else
      data = {
        type: "object",
        properties: {
          release_id: { type: "string" },
          schema_version: { const: 1 },
          provenance: z.toJSONSchema(provenanceOutputSchema),
          citation: ref("citation"),
          licenses: defs.release_manifest.properties.licenses,
          manifest: {
            type: "object",
            properties: { sha256: { type: "string" } },
            required: ["sha256"],
            additionalProperties: false,
          },
          species: {
            type: "object",
            properties: Object.fromEntries(
              ["human", "mouse"].map((s) => [
                s,
                {
                  type: "object",
                  properties: {
                    taxon_id: { type: "integer" },
                    counts: defs.species_manifest.properties.counts,
                    capabilities: defs.species_manifest.properties.capabilities,
                    references: {
                      type: "array",
                      minItems: 4,
                      maxItems: 4,
                      items: {
                        ...defs.download,
                        properties: {
                          ...defs.download.properties,
                          name: {
                            enum: [
                              "schema.json",
                              "citation.json",
                              "openantigens.bib",
                              "openantigens.ris",
                            ],
                          },
                        },
                      },
                    },
                  },
                  required: [
                    "taxon_id",
                    "counts",
                    "capabilities",
                    "references",
                  ],
                  additionalProperties: false,
                },
              ]),
            ),
            required: ["human", "mouse"],
            additionalProperties: false,
          },
        },
        required: [
          "release_id",
          "schema_version",
          "provenance",
          "citation",
          "licenses",
          "manifest",
          "species",
        ],
        additionalProperties: false,
      };
    const envelope = z.toJSONSchema(outputSchema);
    envelope.$defs = defs;
    envelope.properties.data = data;
    if (name === "get_release")
      envelope.properties.releases = { type: "array", items: data };
    return z.fromJSONSchema(envelope);
  }
  fallback(release_id, species) {
    return this.releases.get(release_id)?.manifest.species[species]
      ?.records_download;
  }
}
export function fasta(sequence, { release_id, provenance_url, records_url }) {
  const encode = (value) =>
    encodeURIComponent(value).replace(
      /[!'()*]/g,
      (c) => `%${c.charCodeAt(0).toString(16).toUpperCase()}`,
    );
  const fields = {
    release_id,
    target_id: sequence.target_id,
    construct_id: sequence.construct_id,
    species: sequence.species,
    taxon_id: sequence.taxon_id,
    identifier_source: sequence.identifier_source,
    identifier_accession: sequence.identifier_accession,
    sequence_accession: sequence.sequence_accession,
    sequence_source: sequence.sequence_source,
    sequence_basis: sequence.sequence_basis,
    start: sequence.start,
    end: sequence.end,
    length: sequence.length,
    sha256: sequence.sequence_sha256,
    warnings: canonical(sequence.warnings),
    evidence: canonical(sequence.evidence),
    provenance_url,
    records_url,
    record_selector: canonical({
      target_id: sequence.target_id,
      construct_id: sequence.construct_id,
    }),
  };
  const primary = encode(
    `${release_id}:${sequence.target_id}:${sequence.construct_id ?? "reference"}`,
  );
  return `>${primary} ${Object.entries(fields)
    .map(
      ([key, value]) =>
        `${key}=${encode(value === null ? "null" : String(value))}`,
    )
    .join(" ")}\n${sequence.sequence.match(/.{1,80}/g)?.join("\n") ?? ""}\n`;
}

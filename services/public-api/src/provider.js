import { open, lstat, realpath } from "node:fs/promises";
import { constants, createReadStream } from "node:fs";
import { createGunzip } from "node:zlib";
import { resolve, relative, sep } from "node:path";
import { createHash } from "node:crypto";
import { z } from "zod";
import {
  outputSchema,
  provenanceOutputSchema,
  rootSchema,
  ApiError,
  bad,
  LIMITS,
  SECTIONS,
  parseInput,
} from "./contracts.js";
const parseJson = (bytes) =>
  JSON.parse(new TextDecoder("utf-8", { fatal: true }).decode(bytes));
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
const notFound = () => {
  throw new ApiError(
    404,
    "not_found",
    "The requested released record does not exist.",
  );
};
export class Provider {
  constructor(directory) {
    this.directory = resolve(directory);
    this.ready = false;
    this.releases = new Map();
    this.retired = new Map();
  }
  async readSafe(path, max, expected, collect = true) {
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
      const chunks = [];
      const checksum = createHash("sha256");
      let length = 0;
      for (;;) {
        const chunk = Buffer.alloc(Math.min(65536, max - length + 1));
        const { bytesRead } = await handle.read(chunk);
        if (!bytesRead) break;
        length += bytesRead;
        if (length > max || (expected && length > expected.bytes))
          throw new Error("Bundle size");
        const bytes = chunk.subarray(0, bytesRead);
        checksum.update(bytes);
        if (collect) chunks.push(bytes);
      }
      if (
        expected &&
        (length !== expected.bytes ||
          checksum.digest("hex") !== expected.sha256)
      )
        throw new Error("Bundle checksum");
      return collect ? Buffer.concat(chunks, length) : null;
    } finally {
      await handle.close();
    }
  }
  async initialize() {
    this.ready = false;
    this.releases.clear();
    this.retired.clear();
    this.directory = await realpath(this.directory);
    const root = parseJson(await this.readSafe("manifest.json", LIMITS.shard));
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
      const bytes = await this.readSafe(pointer.manifest_path, LIMITS.manifest);
      if (digest(bytes) !== pointer.sha256) throw new Error("Release checksum");
      const manifest = parseJson(bytes);
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
      if (!validators.release_manifest?.safeParse(manifest).success)
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
      };
      for (const [path, file] of Object.entries(manifest.files)) {
        if (
          !hex.test(file.sha256) ||
          !Number.isSafeInteger(file.bytes) ||
          file.bytes < 0
        )
          throw new Error("File descriptor");
        const bound = path.endsWith("/index.json")
          ? LIMITS.index
          : path.includes("/targets/")
            ? LIMITS.target
            : file.media_type === "application/json" &&
                !path.includes("/schemas/")
              ? LIMITS.shard
              : 128 * 1024 * 1024;
        if (
          file.media_type === "application/schema+json" &&
          file.sha256 !== manifest.files[schemaPaths[0]].sha256
        )
          throw new Error("Schema disagreement");
        const parse =
          file.media_type === "application/json" &&
          file.schema &&
          !path.includes("/schemas/") &&
          !path.includes("/downloads/");
        const content = await this.readSafe(path, bound, file, Boolean(parse));
        if (parse) {
          const envelope = parseJson(content);
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
        }
      }
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
        const index = await this.json(info, entry.index_path, LIMITS.index);
        totalIndexBytes += manifest.files[entry.index_path].bytes;
        if (
          !Array.isArray(index) ||
          index.some(
            (row) =>
              row.species !== species ||
              row.taxon_id !== entry.taxon_id ||
              typeof row.record_id !== "string" ||
              (row.target_id && !entry.target_paths[row.target_id]),
          )
        )
          throw new Error("Index schema");
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
            !url.pathname.includes(`/data/releases/${manifest.release_id}/`) ||
            !file ||
            file.sha256 !== download.sha256 ||
            file.bytes !== download.bytes
          )
            throw new Error("Download manifest");
        }
        for (const [id, path] of Object.entries(entry.target_paths)) {
          if (!manifest.files[path]) throw new Error("Target mapping");
          const target = await this.json(info, path, LIMITS.target);
          if (
            target.summary.target_id !== id ||
            target.summary.species !== species ||
            !id.startsWith(`${entry.taxon_id}:`)
          )
            throw new Error("Target identity");
          const refs = [
            target.reference_sequence_path,
            ...Object.values(target.construct_sequence_paths),
            ...target.constructs.pages,
            ...Object.values(target.sections).flatMap((d) => d.pages),
          ].filter(Boolean);
          if (refs.some((p) => !manifest.files[p]))
            throw new Error("Target shard reference");
        }
      }
      if (totalIndexBytes > LIMITS.indexes)
        throw new Error("Index memory budget");
      await this.verifySemanticIdentity(info);
      for (const species of ["human", "mouse"])
        await this.verifyRecordsDownload(info, species);
      delete info.envelopeDigests;
      delete info.sourceIndexes;
      delete info.targetRecords;
      this.releases.set(manifest.release_id, info);
    }
    if (!this.releases.has(root.current_release_id))
      throw new Error("Current release");
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
          const record = JSON.parse(line);
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
        const summary = indexes[position];
        status[summary.status] = (status[summary.status] ?? 0) + 1;
        const path = info.targetRecords.get(`${species}:${summary.record_id}`);
        if (!path) throw new Error("Missing projected target");
        const target = await this.json(info, path, LIMITS.target);
        if (canonical(target.summary) !== canonical(summary))
          throw new Error("Index target disagreement");
        if (position) write(",");
        write('{"construct_sequences":{');
        let first = true;
        for (const cid of Object.keys(target.construct_sequence_paths).sort()) {
          if (!first) write(",");
          first = false;
          const sequencePath = target.construct_sequence_paths[cid];
          const sequence = sequencePath
            ? await this.json(info, sequencePath)
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
          write(`${canonical(section)}:{"items":[`);
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
  async json(info, path, max = LIMITS.shard) {
    try {
      const file = info.manifest.files[path];
      if (!file) throw new Error("Unlisted shard");
      const result = parseJson(await this.readSafe(path, max, file));
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
    } = info.manifest.provenance;
    const manifest = this.manifestReference(info);
    return {
      software_commit,
      snapshot_id,
      source_versions,
      source_dates,
      notes,
      input_hashes: {
        count: Object.keys(input_sha256).length,
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
            downloads: v.downloads,
          },
        ]),
      ),
    };
  }
  async call(name, raw, signal) {
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
              retired_releases: [...this.retired.values()],
            },
      );
    const sm = info.manifest.species[species];
    if (name === "get_downloads")
      return this.envelope(info, species, sm.downloads);
    if (name === "search_targets") {
      if (args.disease_id && !sm.capabilities.disease)
        throw new ApiError(
          400,
          "unsupported_capability",
          "Disease associations are not supported for this species and release.",
        );
      const q = (args.q ?? "").toLocaleLowerCase("en").trim();
      const rows = info.indexes[species].filter(
        (r) =>
          (!q || r.search_text?.toLocaleLowerCase("en").includes(q)) &&
          (!args.track || r.track === args.track) &&
          (!args.status || r.status === args.status) &&
          (!args.structure || r.structure_state === args.structure) &&
          (!args.constructs || r.construct_state === args.constructs) &&
          (!args.disease_id || r.disease_ids.includes(args.disease_id)),
      );
      return this.page(info, species, rows, args, "targets");
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
        { state: descriptor.state, reason: descriptor.reason },
        args.section === "similarity" ? 50 : 100,
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
        100,
        signal,
      );
    }
    if (name === "get_sequence") {
      const path = args.construct_id
        ? target.construct_sequence_paths[args.construct_id]
        : target.reference_sequence_path;
      if (!path) notFound();
      const sequence = await this.json(info, path);
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
  pagination(info, species, args, scope, maximum = 100) {
    const limit = args.limit ?? 25;
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
        const c = JSON.parse(Buffer.from(args.cursor, "base64url"));
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
  page(info, species, rows, args, scope, extra = {}, maximum = 100) {
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
    maximum = 100,
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
    else if (name === "get_downloads")
      data = { type: "array", items: ref("download") };
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
                    downloads: { type: "array", items: ref("download") },
                  },
                  required: ["taxon_id", "counts", "capabilities", "downloads"],
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

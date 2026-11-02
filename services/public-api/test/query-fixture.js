import { mkdtemp, mkdir, readFile, writeFile, rm } from "node:fs/promises";
import { join, dirname } from "node:path";
import { tmpdir } from "node:os";
import { createHash } from "node:crypto";
import { packedFixture } from "./fixture.js";
const sha = (bytes) => createHash("sha256").update(bytes).digest("hex");
export async function queryFixture(source) {
  const directory = await mkdtemp(join(tmpdir(), "oa-query-runtime-")),
    staticDirectory = await mkdtemp(join(tmpdir(), "oa-query-static-"));
  const packed = await packedFixture(source);
  const sourceRootBytes = await readFile(join(source, "manifest.json")),
    sourceRoot = JSON.parse(sourceRootBytes);
  const transportBytes = await readFile(
      join(packed.directory, "transport-index.json"),
    ),
    transport = JSON.parse(transportBytes);
  const write = async (path, bytes) => {
    await mkdir(dirname(join(directory, path)), { recursive: true });
    await writeFile(join(directory, path), bytes);
    return { path, bytes: bytes.length, sha256: sha(bytes) };
  };
  const root = {
    format: "openantigens-own-static-query-v1",
    source_root_manifest_sha256: sha(sourceRootBytes),
    current_release_id: sourceRoot.current_release_id,
    retired_releases: sourceRoot.retired_releases ?? [],
    transport_index: await write("transport-index.json", transportBytes),
    releases: [],
    static_origin: "https://openantigens.org",
  };
  for (const pointer of sourceRoot.releases) {
    const original = JSON.parse(
      await readFile(join(source, pointer.manifest_path)),
    );
    const rid = pointer.release_id;
    const url = `https://openantigens.org/data/releases/${rid}/manifest.json`;
    const provenance = {
      ...Object.fromEntries(
        [
          "software_commit",
          "snapshot_id",
          "source_versions",
          "source_dates",
          "notes",
        ].map((key) => [key, original.provenance[key]]),
      ),
      projection_methods: original.provenance.projection_methods ?? null,
      input_hashes: {
        count: Object.keys(original.provenance.input_sha256).length,
        manifest_url: url,
        manifest_sha256: pointer.sha256,
        json_pointer: "/provenance/input_sha256",
      },
    };
    const meta = {
      format: "openantigens-query-release-v1",
      release_id: rid,
      schema_version: 1,
      source_manifest: { url, sha256: pointer.sha256 },
      source_payload_sha256: original.core.payload_sha256,
      citation: original.citation,
      licenses: original.licenses,
      provenance,
      species: {},
      public_schema: null,
    };
    const schemaPath = Object.keys(original.files).find(
      (path) => original.files[path].media_type === "application/schema+json",
    );
    meta.public_schema = await write(
      schemaPath,
      await readFile(join(source, schemaPath)),
    );
    for (const species of ["human", "mouse"]) {
      const item = original.species[species],
        indexBytes = await readFile(join(source, item.index_path)),
        rows = JSON.parse(indexBytes).data;
      const summaryDownload = item.downloads.find(
        (download) => download.name === "index.jsonl",
      );
      const bytes = await readFile(join(source, summaryDownload.path));
      const lines = bytes.toString("utf8").trimEnd().split("\n");
      let position = 0;
      const byRecord = new Map();
      for (const text of lines) {
        const line = Buffer.from(text + "\n");
        byRecord.set(JSON.parse(line).record_id, [
          position,
          line.length,
          sha(line),
        ]);
        position += line.length;
      }
      const index = await write(item.index_path, indexBytes),
        summaryFile = await write(summaryDownload.path, bytes);
      const offsets = {
        format: "openantigens-summary-offsets-v1",
        index_sha256: index.sha256,
        jsonl_sha256: summaryFile.sha256,
        rows: rows.map((row) => byRecord.get(row.record_id)),
      };
      const offsetFile = await write(
        `releases/${rid}/${species}/summary-offsets.json`,
        Buffer.from(JSON.stringify(offsets) + "\n"),
      );
      meta.species[species] = {
        taxon_id: item.taxon_id,
        index,
        counts: item.counts,
        capabilities: item.capabilities,
        downloads: item.downloads,
        records_download_name: "records.jsonl.gz",
        summary_lines: { ...summaryFile, offsets: offsetFile },
      };
    }
    const descriptor = await write(
      `releases/${rid}/query-release.json`,
      Buffer.from(JSON.stringify(meta) + "\n"),
    );
    root.releases.push({ release_id: rid, ...descriptor });
  }
  for (const [target, parts] of Object.entries(transport.targets)) {
    const [, rid, species, , leaf] = target.split("/");
    for (const [number, part] of parts.entries()) {
      const path = `api-runtime/${rid}/${species}/${leaf.slice(0, -5)}-${number}-${part.sha256}.oap`;
      await mkdir(dirname(join(staticDirectory, path)), { recursive: true });
      await writeFile(
        join(staticDirectory, path),
        await readFile(join(packed.directory, part.path)),
      );
    }
  }
  const query = await write(
    "query-root.json",
    Buffer.from(JSON.stringify(root) + "\n"),
  );
  await rm(packed.directory, { recursive: true, force: true });
  return {
    directory,
    staticDirectory,
    seal: {
      format: root.format,
      query_root_sha256: query.sha256,
      source_root_manifest_sha256: root.source_root_manifest_sha256,
      static_origin: root.static_origin,
    },
  };
}

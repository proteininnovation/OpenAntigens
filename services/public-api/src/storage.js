import { AsyncLocalStorage } from "node:async_hooks";
import { createHash } from "node:crypto";
import { lstat, readdir } from "node:fs/promises";
import { join } from "node:path";
import { promisify } from "node:util";
import { gunzip } from "node:zlib";
import { parsePublicJson } from "./json.js";

const FORMAT = "openantigens-target-gzip-v1";
const DECODED_LIMIT = 16 * 1024 * 1024;
const COMPRESSED_LIMIT = 17 * 1024 * 1024;
const HEX = /^[a-f0-9]{64}$/;
const TARGET =
  /^releases\/([a-f0-9]{64})\/(human|mouse)\/targets\/([a-f0-9]{64})\.json$/;
const EXTERNAL =
  /^(?:index\.jsonl|index\.tsv|constructs\.jsonl|constructs\.tsv|constructs\.fasta|reference\.fasta|records\.jsonl\.gz)$/;
const inflate = promisify(gunzip);
const sha = (bytes) => createHash("sha256").update(bytes).digest("hex");
const fail = (message) => {
  throw new Error(`Packed storage: ${message}`);
};
const object = (v) => v !== null && typeof v === "object" && !Array.isArray(v);
const keys = (v, expected) =>
  object(v) && Object.keys(v).sort().join(",") === expected;
const positive = (v, max = Number.MAX_SAFE_INTEGER) =>
  Number.isSafeInteger(v) && v > 0 && v <= max;
const canonical = (v) =>
  JSON.stringify(v, (_, item) =>
    object(item)
      ? Object.fromEntries(
          Object.keys(item)
            .sort()
            .map((key) => [key, item[key]]),
        )
      : item,
  );
const parse = (bytes, strict = false) => {
  const value = parsePublicJson(bytes);
  if (strict && !Buffer.from(`${canonical(value)}\n`).equals(bytes))
    fail("noncanonical container");
  return value;
};

function virtual(path) {
  return (
    typeof path === "string" &&
    /^releases\/[a-f0-9]{64}\/(human|mouse)\/(targets|constructs|sequences|sections)\//.test(
      path,
    )
  );
}

function inferredTarget(path) {
  const target = TARGET.exec(path);
  if (target) return path;
  const match =
    /^releases\/([a-f0-9]{64})\/(human|mouse)\/(?:constructs\/([a-f0-9]{64})-(?:0|[1-9][0-9]*)\.json|sections\/([a-f0-9]{64})\/[a-z]+-(?:0|[1-9][0-9]*)\.json|sequences\/([a-f0-9]{64})-reference\.json)$/.exec(
      path,
    );
  return match
    ? `releases/${match[1]}/${match[2]}/targets/${match[3] ?? match[4] ?? match[5]}.json`
    : null;
}

export async function createStorage(directory, pins = {}, readPhysical) {
  const { expectedRootSha256, expectedTransportSha256 } = pins;
  let present = false;
  try {
    await lstat(join(directory, "transport-index.json"));
    present = true;
  } catch (error) {
    if (error.code !== "ENOENT") throw error;
  }
  const pinned =
    expectedRootSha256 !== undefined || expectedTransportSha256 !== undefined;
  if (!pinned) {
    if (present) fail("transport requires explicit pins");
    return {
      packed: false,
      read: readPhysical,
      registerManifest() {},
      isExternal: () => false,
      isVirtual: () => false,
      async *scanRelease() {},
      withContext: (fn) => fn(),
      async validateInventory() {},
    };
  }
  if (
    !HEX.test(expectedRootSha256 ?? "") ||
    !HEX.test(expectedTransportSha256 ?? "")
  )
    fail("both pins are required");
  const rootBytes = await readPhysical("manifest.json", 1024 * 1024);
  if (sha(rootBytes) !== expectedRootSha256) fail("root pin mismatch");
  const sidecar = await readPhysical("manifest.sha256", 256);
  if (sidecar.toString("utf8") !== `${expectedRootSha256}  manifest.json\n`)
    fail("root sidecar mismatch");
  const indexBytes = await readPhysical("transport-index.json", DECODED_LIMIT);
  if (sha(indexBytes) !== expectedTransportSha256)
    fail("transport pin mismatch");
  const index = parse(indexBytes, true);
  if (
    !keys(index, "format,root_manifest_sha256,targets") ||
    index.format !== FORMAT ||
    index.root_manifest_sha256 !== expectedRootSha256 ||
    !object(index.targets)
  )
    fail("transport index shape");
  for (const [target, parts] of Object.entries(index.targets)) {
    const match = TARGET.exec(target);
    if (!match || !Array.isArray(parts) || !parts.length)
      fail("target descriptor");
    let previous = "";
    for (const [i, part] of parts.entries()) {
      if (
        !keys(
          part,
          "bytes,decoded_bytes,max_path,members,min_path,path,sha256",
        ) ||
        part.path !==
          `packs/${match[1]}/${match[2]}/${match[3]}-${i}.json.gz` ||
        !HEX.test(part.sha256 ?? "") ||
        !positive(part.bytes, COMPRESSED_LIMIT) ||
        !positive(part.decoded_bytes, DECODED_LIMIT) ||
        !positive(part.members) ||
        typeof part.min_path !== "string" ||
        typeof part.max_path !== "string" ||
        !virtual(part.min_path) ||
        !virtual(part.max_path) ||
        part.min_path > part.max_path ||
        (previous && previous >= part.min_path)
      )
        fail("part descriptor or range");
      previous = part.max_path;
    }
  }
  const manifests = new Map();
  const contexts = new AsyncLocalStorage();
  const withContext = (fn) =>
    contexts.getStore() ? fn() : contexts.run(new Map(), fn);
  const sourceFile = (path) => manifests.get(path.split("/")[1])?.files[path];
  const isExternal = (path) => {
    const match =
      /^releases\/[a-f0-9]{64}\/(human|mouse)\/downloads\/([^/]+)$/.exec(path);
    return Boolean(match && EXTERNAL.test(match[2]));
  };

  async function loadPart(part) {
    const cache = contexts.getStore();
    if (cache.has(part.path)) {
      const cached = cache.get(part.path);
      cache.delete(part.path);
      cache.set(part.path, cached);
      return cached;
    }
    const compressed = await readPhysical(part.path, COMPRESSED_LIMIT, part);
    const decoded = await inflate(compressed, {
      maxOutputLength: DECODED_LIMIT,
    });
    if (decoded.length !== part.decoded_bytes) fail("decoded length");
    const container = parse(decoded, true);
    if (
      !keys(container, "format,members") ||
      container.format !== FORMAT ||
      !Array.isArray(container.members) ||
      container.members.length !== part.members
    )
      fail("member count or shape");
    let previous = "";
    for (const member of container.members) {
      if (!Array.isArray(member) || member.length !== 2) fail("member tuple");
      const [path, encoded] = member;
      const file = typeof path === "string" && sourceFile(path);
      if (
        !virtual(path) ||
        !file ||
        (previous && previous >= path) ||
        typeof encoded !== "string"
      )
        fail("unknown or duplicate member");
      if (
        !Number.isSafeInteger(file.bytes) ||
        file.bytes < 0 ||
        file.bytes > (TARGET.test(path) ? 8 : 1) * 1024 * 1024
      )
        fail("member size");
      const bytes = Buffer.from(encoded, "base64");
      if (
        bytes.toString("base64") !== encoded ||
        bytes.length !== file.bytes ||
        sha(bytes) !== file.sha256
      )
        fail("member bytes or base64");
      previous = path;
    }
    if (container.members[0][0] !== part.min_path || previous !== part.max_path)
      fail("member range");
    cache.set(part.path, container.members);
    while (cache.size > 2) cache.delete(cache.keys().next().value);
    return container.members;
  }

  async function read(path, max, expected, collect = true, targetPath) {
    if (isExternal(path)) fail("download bytes are external");
    if (path === "manifest.json")
      return readPhysical(
        path,
        max,
        { bytes: rootBytes.length, sha256: expectedRootSha256 },
        collect,
      );
    if (!virtual(path)) return readPhysical(path, max, expected, collect);
    return withContext(async () => {
      const file = sourceFile(path);
      if (
        !file ||
        file.bytes > max ||
        (expected &&
          (expected.bytes !== file.bytes || expected.sha256 !== file.sha256))
      )
        fail("unlisted member or size");
      const inferred = inferredTarget(path);
      if (targetPath && inferred && targetPath !== inferred)
        fail("target context mismatch");
      const target = targetPath ?? inferred;
      if (
        !target ||
        !TARGET.test(target) ||
        target.split("/").slice(0, 3).join("/") !==
          path.split("/").slice(0, 3).join("/")
      )
        fail("target context required");
      const parts = index.targets[target];
      if (!parts) fail("missing target group");
      const part = parts.find(
        (item) => item.min_path <= path && item.max_path >= path,
      );
      if (!part) fail("missing member range");
      const members = await loadPart(part);
      const member = members.find(([candidate]) => candidate === path);
      if (!member) fail("missing member");
      return collect ? Buffer.from(member[1], "base64") : null;
    });
  }

  function references(targetPath, target) {
    const paths = [targetPath];
    if (
      !object(target) ||
      !object(target.sections) ||
      !Array.isArray(target.constructs?.pages) ||
      !object(target.construct_sequence_paths)
    )
      fail("target references");
    paths.push(...target.constructs.pages);
    for (const section of Object.values(target.sections)) {
      if (!Array.isArray(section.pages)) fail("section references");
      paths.push(...section.pages);
    }
    if (target.reference_sequence_path)
      paths.push(target.reference_sequence_path);
    paths.push(
      ...Object.values(target.construct_sequence_paths).filter(
        (path) => path !== null,
      ),
    );
    const allowed = new Set();
    for (const path of paths) {
      if (
        !virtual(path) ||
        !sourceFile(path) ||
        allowed.has(path) ||
        path.split("/").slice(0, 3).join("/") !==
          targetPath.split("/").slice(0, 3).join("/")
      )
        fail("source references");
      const inferred = inferredTarget(path);
      if (inferred && inferred !== targetPath) fail("source group ownership");
      allowed.add(path);
    }
    return allowed;
  }

  async function* scanRelease(manifest) {
    const seen = new Set();
    for (const [targetPath, parts] of Object.entries(index.targets)) {
      if (targetPath.split("/")[1] !== manifest.release_id) continue;
      const targetBytes = await read(
        targetPath,
        8 * 1024 * 1024,
        manifest.files[targetPath],
      );
      const allowed = references(targetPath, parse(targetBytes).data);
      const groupSeen = new Set();
      for (const part of parts) {
        const members = await withContext(() => loadPart(part));
        for (const [path, encoded] of members) {
          if (!allowed.has(path) || seen.has(path) || groupSeen.has(path))
            fail("unreachable or duplicate member");
          seen.add(path);
          groupSeen.add(path);
          yield { path, bytes: Buffer.from(encoded, "base64") };
        }
      }
      if (groupSeen.size !== allowed.size) fail("missing source member");
    }
    for (const path of Object.keys(manifest.files)) {
      if (virtual(path) && !seen.has(path)) fail("missing manifest member");
    }
  }

  async function validateInventory(root) {
    const expected = new Set([
      "manifest.json",
      "manifest.sha256",
      "transport-index.json",
    ]);
    const releaseIds = new Set();
    for (const pointer of root.releases) {
      if (
        !manifests.has(pointer.release_id) ||
        pointer.manifest_path !== `releases/${pointer.release_id}/manifest.json`
      )
        fail("release inventory");
      releaseIds.add(pointer.release_id);
      expected.add(pointer.manifest_path);
    }
    if (releaseIds.size !== manifests.size) fail("release inventory");
    for (const manifest of manifests.values()) {
      for (const path of Object.keys(manifest.files)) {
        if (!virtual(path) && !isExternal(path)) expected.add(path);
      }
    }
    for (const [target, parts] of Object.entries(index.targets)) {
      if (!sourceFile(target) || !releaseIds.has(target.split("/")[1]))
        fail("unknown target inventory");
      for (const part of parts) expected.add(part.path);
    }
    const directories = new Set();
    for (const path of expected) {
      const parts = path.split("/");
      for (let i = 1; i < parts.length; i++)
        directories.add(parts.slice(0, i).join("/"));
    }
    async function visit(prefix = "") {
      for (const entry of await readdir(join(directory, prefix), {
        withFileTypes: true,
      })) {
        const path = prefix ? `${prefix}/${entry.name}` : entry.name;
        const stat = await lstat(join(directory, path));
        if (stat.isSymbolicLink()) fail("inventory symlink");
        if (stat.isDirectory()) {
          if (!directories.has(path)) fail("unknown directory");
          await visit(path);
        } else if (!stat.isFile() || !expected.delete(path))
          fail("unknown physical file");
      }
    }
    await visit();
    if (expected.size) fail("missing physical file");
  }
  return {
    packed: true,
    read,
    isExternal,
    isVirtual: virtual,
    withContext,
    scanRelease,
    validateInventory,
    registerManifest(manifest) {
      if (
        !HEX.test(manifest.release_id ?? "") ||
        !object(manifest.files) ||
        manifests.has(manifest.release_id)
      )
        fail("manifest registration");
      manifests.set(manifest.release_id, manifest);
    },
  };
}

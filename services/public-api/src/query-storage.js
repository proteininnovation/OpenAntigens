import { AsyncLocalStorage } from "node:async_hooks";
import { createHash } from "node:crypto";
import { isDeepStrictEqual, promisify } from "node:util";
import { gunzip } from "node:zlib";
import { ApiError, SECTIONS } from "./contracts.js";
import { parsePublicJson } from "./json.js";

const MiB = 1024 * 1024;
const TARGET =
  /^releases\/([a-f0-9]{64})\/(human|mouse)\/targets\/([a-f0-9]{64})\.json$/;
const HEX = /^[a-f0-9]{64}$/;
const inflate = promisify(gunzip);
const sha = (bytes) => createHash("sha256").update(bytes).digest("hex");
const unavailable = () =>
  new ApiError(
    503,
    "unavailable",
    "The requested released data is temporarily unavailable.",
  );
const requireValue = (condition) => {
  if (!condition) throw unavailable();
};
const keys = (value, names) =>
  value &&
  typeof value === "object" &&
  !Array.isArray(value) &&
  Object.keys(value).sort().join(",") === names;
const positive = (value, maximum) =>
  Number.isSafeInteger(value) && value > 0 && value <= maximum;

function inferTarget(path) {
  if (TARGET.test(path)) return path;
  const match =
    /^releases\/([a-f0-9]{64})\/(human|mouse)\/(?:constructs\/([a-f0-9]{64})-(?:0|[1-9][0-9]*)\.json|sections\/([a-f0-9]{64})\/[a-z]+-(?:0|[1-9][0-9]*)\.json|sequences\/([a-f0-9]{64})-reference\.json)$/.exec(
      path,
    );
  return match
    ? `releases/${match[1]}/${match[2]}/targets/${match[3] ?? match[4] ?? match[5]}.json`
    : null;
}

function abortable(promise, signal) {
  if (signal.aborted) return Promise.reject(unavailable());
  return new Promise((resolve, reject) => {
    const abort = () => reject(unavailable());
    signal.addEventListener("abort", abort, { once: true });
    Promise.resolve(promise)
      .then(resolve, reject)
      .finally(() => signal.removeEventListener("abort", abort));
  });
}

export function createRemoteQueryStorage({
  transportIndex,
  releases,
  readSummary,
  fetchImpl = globalThis.fetch,
  now = Date.now,
}) {
  const contexts = new AsyncLocalStorage();
  const cache = new Map();
  let cacheBytes = 0,
    activeTargets = 0,
    activeFetches = 0;

  const check = (context) =>
    requireValue(!context.signal.aborted && now() < context.deadline);
  async function withContext(fn, { signal } = {}) {
    if (contexts.getStore()) return fn();
    const controller = new AbortController();
    const abort = () => controller.abort();
    if (signal?.aborted) controller.abort();
    signal?.addEventListener("abort", abort, { once: true });
    const timer = setTimeout(abort, 3000);
    const context = {
      signal: controller.signal,
      deadline: now() + 3000,
      parts: new Map(),
      pending: new Map(),
      reserved: false,
    };
    try {
      return await contexts.run(context, fn);
    } finally {
      clearTimeout(timer);
      signal?.removeEventListener("abort", abort);
      controller.abort();
      if (context.reserved) activeTargets--;
    }
  }

  function store(key, bytes) {
    if (cache.has(key)) {
      cacheBytes -= cache.get(key).length;
      cache.delete(key);
    }
    cache.set(key, bytes);
    cacheBytes += bytes.length;
    while (cacheBytes > 8 * MiB || cache.size > 128) {
      const oldest = cache.keys().next().value;
      cacheBytes -= cache.get(oldest).length;
      cache.delete(oldest);
    }
  }

  async function compressed(context, descriptor) {
    check(context);
    const key = `${context.targetPath}:${descriptor.path}:${descriptor.sha256}`;
    if (cache.has(key)) {
      const bytes = cache.get(key);
      cache.delete(key);
      cache.set(key, bytes);
      return { key, bytes, cached: true };
    }
    requireValue(activeFetches < 2);
    activeFetches++;
    let reader, response;
    try {
      const match =
        /^packs\/([a-f0-9]{64})\/(human|mouse)\/([a-f0-9]{64})-(0|[1-9][0-9]*)\.json\.gz$/.exec(
          descriptor.path,
        );
      requireValue(match && HEX.test(descriptor.sha256));
      const url = `https://openantigens.org/api-runtime/${match[1]}/${match[2]}/${match[3]}-${match[4]}-${descriptor.sha256}.oap`;
      response = await abortable(
        fetchImpl(url, {
          method: "GET",
          credentials: "omit",
          redirect: "error",
          signal: context.signal,
          headers: {
            "accept-encoding": "identity",
            accept: "application/octet-stream",
            "user-agent": "OpenAntigens-Public-API/1",
          },
        }),
        context.signal,
      );
      requireValue(
        response.status === 200 &&
          !response.redirected &&
          (!response.url || response.url === url),
      );
      requireValue(
        response.headers.get("content-type")?.trim().toLowerCase() ===
          "application/octet-stream",
      );
      const encoding = response.headers.get("content-encoding");
      requireValue(
        encoding === null || encoding.trim().toLowerCase() === "identity",
      );
      const length = response.headers.get("content-length");
      requireValue(
        length === null ||
          (/^(0|[1-9][0-9]*)$/.test(length) &&
            Number(length) === descriptor.bytes),
      );
      requireValue(response.body?.getReader);
      reader = response.body.getReader();
      const chunks = [];
      let size = 0;
      for (;;) {
        check(context);
        const { done, value } = await abortable(reader.read(), context.signal);
        if (done) break;
        requireValue(value instanceof Uint8Array);
        size += value.byteLength;
        requireValue(size <= descriptor.bytes);
        chunks.push(Buffer.from(value));
      }
      const bytes = Buffer.concat(chunks, size);
      requireValue(
        size === descriptor.bytes && sha(bytes) === descriptor.sha256,
      );
      return { key, bytes, cached: false };
    } finally {
      activeFetches--;
      if (reader) {
        void reader.cancel().catch(() => {});
        reader.releaseLock();
      } else if (response?.body) void response.body.cancel().catch(() => {});
    }
  }

  function role(context, path) {
    if (path === context.targetPath) return "target_shard";
    const prefix = context.prefix;
    if (
      new RegExp(
        `^${prefix}/constructs/${context.key}-(0|[1-9][0-9]*)\\.json$`,
      ).test(path)
    )
      return "construct_page";
    if (
      new RegExp(
        `^${prefix}/sequences/(?:${context.key}-reference|[a-f0-9]{64})\\.json$`,
      ).test(path)
    )
      return "sequence_shard";
    const section = new RegExp(
      `^${prefix}/sections/${context.key}/([a-z]+)-(0|[1-9][0-9]*)\\.json$`,
    ).exec(path);
    return section && SECTIONS.includes(section[1])
      ? `section_${section[1]}`
      : null;
  }

  async function decode(context, descriptor, bytes) {
    const decoded = await inflate(bytes, { maxOutputLength: 16 * MiB });
    check(context);
    requireValue(decoded.length === descriptor.decoded_bytes);
    const container = parsePublicJson(decoded);
    requireValue(
      keys(container, "format,members") &&
        container.format === "openantigens-target-gzip-v1" &&
        Array.isArray(container.members) &&
        container.members.length === descriptor.members,
    );
    const members = new Map();
    let previous = "";
    for (const item of container.members) {
      requireValue(Array.isArray(item) && item.length === 2);
      const [path, encoded] = item,
        schema = typeof path === "string" && role(context, path);
      requireValue(
        schema && typeof encoded === "string" && (!previous || previous < path),
      );
      const raw = Buffer.from(encoded, "base64");
      requireValue(
        raw.toString("base64") === encoded &&
          raw.length <= (schema === "target_shard" ? 8 : 1) * MiB,
      );
      const envelope = parsePublicJson(raw);
      requireValue(
        envelope.release_id === context.info.manifest.release_id &&
          envelope.species === context.species &&
          envelope.schema_version === 1 &&
          context.info.validators[schema]?.safeParse(envelope).success,
      );
      members.set(path, envelope.data);
      previous = path;
    }
    requireValue(
      members.keys().next().value === descriptor.min_path &&
        previous === descriptor.max_path,
    );
    return members;
  }

  function graph(context) {
    const target = context.target,
      allowed = new Map([[context.targetPath, "target_shard"]]);
    const add = (path, schema) => {
      requireValue(
        typeof path === "string" &&
          !allowed.has(path) &&
          role(context, path) === schema,
      );
      allowed.set(path, schema);
    };
    const pages = (descriptor, pathFor, schema) => {
      requireValue(
        Array.isArray(descriptor.pages) &&
          Number.isSafeInteger(descriptor.total) &&
          descriptor.total >= 0 &&
          (descriptor.pages.length === 0) === (descriptor.total === 0),
      );
      requireValue(descriptor.total <= descriptor.pages.length * 25);
      descriptor.pages.forEach((path, i) => {
        requireValue(path === pathFor(i));
        add(path, schema);
      });
    };
    for (const [section, descriptor] of Object.entries(target.sections)) {
      requireValue(
        SECTIONS.includes(section) &&
          descriptor.state === target.summary.section_availability[section],
      );
      pages(
        descriptor,
        (i) => `${context.prefix}/sections/${context.key}/${section}-${i}.json`,
        `section_${section}`,
      );
    }
    pages(
      target.constructs,
      (i) => `${context.prefix}/constructs/${context.key}-${i}.json`,
      "construct_page",
    );
    requireValue(
      target.summary.construct_count === null ||
        target.summary.construct_count === target.constructs.total,
    );
    requireValue(
      Object.keys(target.construct_sequence_paths).length ===
        target.constructs.total,
    );
    if (target.reference_sequence_path !== null) {
      requireValue(
        target.reference_sequence_path ===
          `${context.prefix}/sequences/${context.key}-reference.json`,
      );
      add(target.reference_sequence_path, "sequence_shard");
    }
    for (const [cid, path] of Object.entries(target.construct_sequence_paths)) {
      requireValue(HEX.test(cid));
      if (path !== null) {
        requireValue(path === `${context.prefix}/sequences/${cid}.json`);
        add(path, "sequence_shard");
      }
    }
    context.expected = new Map(
      context.descriptors.map((part) => [part.path, []]),
    );
    for (const path of [...allowed.keys()].sort()) {
      const ranges = context.descriptors.filter(
        (part) => part.min_path <= path && path <= part.max_path,
      );
      requireValue(ranges.length === 1);
      context.expected.get(ranges[0].path).push(path);
    }
    for (const part of context.descriptors) {
      const paths = context.expected.get(part.path);
      requireValue(
        paths.length === part.members &&
          paths[0] === part.min_path &&
          paths.at(-1) === part.max_path,
      );
    }
    context.allowed = allowed;
  }

  function validateLoaded(context) {
    const target = context.target,
      summary = target.summary;
    const loaded = new Map();
    for (const [part, members] of context.parts) {
      requireValue(
        isDeepStrictEqual([...members.keys()], context.expected.get(part)),
      );
      for (const [path, data] of members) loaded.set(path, data);
    }
    const checkPages = (descriptor, inspect) => {
      let total = 0,
        complete = true;
      for (const path of descriptor.pages) {
        if (!loaded.has(path)) {
          complete = false;
          continue;
        }
        const count = inspect(loaded.get(path));
        requireValue(count <= 25);
        total += count;
      }
      requireValue(
        total <= descriptor.total && (!complete || total === descriptor.total),
      );
    };
    for (const [section, descriptor] of Object.entries(target.sections)) {
      checkPages(descriptor, (data) => {
        requireValue(
          data.target_id === summary.target_id &&
            data.section === section &&
            data.state === descriptor.state &&
            data.reason === descriptor.reason,
        );
        if (section === "disease" && summary.disease_coverage !== undefined)
          requireValue(
            isDeepStrictEqual(data.coverage, summary.disease_coverage),
          );
        return data.items.length;
      });
      if (section === "disease" && summary.disease_coverage !== undefined) {
        const coverage = summary.disease_coverage;
        requireValue(
          coverage.projected_count === descriptor.total &&
            descriptor.total <= 100 &&
            coverage.state === descriptor.state &&
            (coverage.source_count === null ||
              coverage.source_count >= descriptor.total),
        );
        if (descriptor.pages.every((path) => loaded.has(path))) {
          const ids = [
            ...new Set(
              descriptor.pages
                .flatMap((path) => loaded.get(path).items)
                .map((item) => item.disease_id)
                .filter(Boolean),
            ),
          ].sort();
          requireValue(isDeepStrictEqual(ids, summary.disease_ids));
        }
      }
    }
    const constructs = new Map();
    let previousCid = "";
    checkPages(target.constructs, (data) => {
      for (const item of data) {
        requireValue(
          item.target_id === summary.target_id &&
            Object.hasOwn(target.construct_sequence_paths, item.construct_id) &&
            !constructs.has(item.construct_id) &&
            (!previousCid || previousCid < item.construct_id),
        );
        constructs.set(item.construct_id, item);
        previousCid = item.construct_id;
      }
      return data.length;
    });
    for (const [path, data] of loaded) {
      if (context.allowed.get(path) !== "sequence_shard") continue;
      const cid =
        path === target.reference_sequence_path
          ? null
          : path.split("/").at(-1).slice(0, -5);
      requireValue(
        data.target_id === summary.target_id &&
          data.construct_id === cid &&
          data.species === context.species &&
          data.taxon_id === (context.species === "human" ? 9606 : 10090) &&
          data.length === data.sequence.length &&
          sha(data.sequence) === data.sequence_sha256,
      );
      if (cid === null)
        requireValue(
          data.sequence_sha256 === summary.sequence_sha256 &&
            data.length === summary.sequence_length &&
            data.start === 1 &&
            data.end === data.length,
        );
      const construct = constructs.get(cid);
      if (construct)
        for (const name of ["sequence_sha256", "length", "start", "end"])
          requireValue(construct[name] === data[name]);
    }
  }

  async function load(context, descriptor, first = false) {
    if (context.parts.has(descriptor.path))
      return context.parts.get(descriptor.path);
    if (context.pending.has(descriptor.path))
      return context.pending.get(descriptor.path);
    const task = (async () => {
      const { key, bytes, cached } = await compressed(context, descriptor);
      const members = await decode(context, descriptor, bytes);
      if (first) {
        context.target = members.get(context.targetPath);
        requireValue(context.target);
        const summary = await abortable(
          readSummary(context.info, context.row),
          context.signal,
        );
        requireValue(isDeepStrictEqual(context.target.summary, summary));
        graph(context);
      }
      context.parts.set(descriptor.path, members);
      try {
        validateLoaded(context);
        check(context);
      } catch (error) {
        context.parts.delete(descriptor.path);
        throw error;
      }
      if (!cached) store(key, bytes);
      return members;
    })();
    context.pending.set(descriptor.path, task);
    try {
      return await task;
    } finally {
      context.pending.delete(descriptor.path);
    }
  }

  async function initializeTarget(context, info, targetPath) {
    check(context);
    requireValue(
      activeTargets < 2 && releases.get(info.manifest.release_id) === info,
    );
    const match = TARGET.exec(targetPath),
      row = info.rowsByPath.get(targetPath),
      descriptors = transportIndex.targets[targetPath];
    requireValue(
      match &&
        match[1] === info.manifest.release_id &&
        row &&
        Array.isArray(descriptors) &&
        descriptors.length > 0 &&
        descriptors.length <= 2,
    );
    let compressedBytes = 0,
      decodedBytes = 0,
      previous = "";
    for (const [i, part] of descriptors.entries()) {
      requireValue(
        keys(
          part,
          "bytes,decoded_bytes,max_path,members,min_path,path,sha256",
        ) &&
          part.path ===
            `packs/${match[1]}/${match[2]}/${match[3]}-${i}.json.gz` &&
          HEX.test(part.sha256) &&
          positive(part.bytes, 4 * MiB) &&
          positive(part.decoded_bytes, 16 * MiB) &&
          positive(part.members, 1000000) &&
          typeof part.min_path === "string" &&
          typeof part.max_path === "string" &&
          part.min_path <= part.max_path &&
          (!previous || previous < part.min_path),
      );
      compressedBytes += part.bytes;
      decodedBytes += part.decoded_bytes;
      previous = part.max_path;
    }
    requireValue(compressedBytes <= 4 * MiB && decodedBytes <= 32 * MiB);
    const targetPart = descriptors.find(
      (part) => part.min_path <= targetPath && targetPath <= part.max_path,
    );
    requireValue(targetPart);
    activeTargets++;
    context.reserved = true;
    Object.assign(context, {
      info,
      row,
      descriptors,
      targetPath,
      species: match[2],
      key: match[3],
      prefix: `releases/${match[1]}/${match[2]}`,
    });
    await load(context, targetPart, true);
  }

  async function readJson(info, path, targetPath) {
    if (!contexts.getStore())
      return withContext(() => readJson(info, path, targetPath));
    const context = contexts.getStore();
    try {
      check(context);
      const inferred = typeof path === "string" && inferTarget(path);
      const target = targetPath ?? inferred ?? context.targetPath;
      requireValue(
        typeof target === "string" && (!inferred || inferred === target),
      );
      if (!context.initializing)
        context.initializing = initializeTarget(context, info, target);
      else requireValue(context.targetPath === target && context.info === info);
      await context.initializing;
      check(context);
      requireValue(context.allowed.has(path));
      const descriptor = context.descriptors.find(
        (part) => part.min_path <= path && path <= part.max_path,
      );
      const members = await load(context, descriptor);
      check(context);
      requireValue(members.has(path));
      return members.get(path);
    } catch (error) {
      if (error?.localQueryIntegrity) throw error;
      throw unavailable();
    }
  }

  return {
    withContext,
    readJson,
    snapshot: () => ({
      activeTargets,
      activeFetches,
      cacheBytes,
      cacheEntries: cache.size,
    }),
  };
}

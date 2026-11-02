import test from "node:test";
import assert from "node:assert/strict";
import { readFile, rm } from "node:fs/promises";
import { join } from "node:path";
import { createHash } from "node:crypto";
import { gzipSync, gunzipSync } from "node:zlib";
import { setImmediate as nextTurn } from "node:timers/promises";
import { z } from "zod";
import { fixture, packedFixture } from "./fixture.js";
import { createRemoteQueryStorage } from "../src/query-storage.js";
import { ApiError } from "../src/contracts.js";

const sha = (bytes) => createHash("sha256").update(bytes).digest("hex");
const MiB = 1024 * 1024;
const fail503 = (error) => { assert.equal(error.status, 503); assert.equal(error.code, "unavailable"); return true; };
const response = (bytes, headers = {}, status = 200) => new Response(bytes, { status, headers: { "content-type": "application/octet-stream", ...headers } });

async function setup(t, { extraUnresolved = 0 } = {}) {
  const source = await fixture("remote-storage", { extraUnresolved });
  const packed = await packedFixture(source.directory);
  t.after(async () => { await rm(source.directory, { recursive: true, force: true }); await rm(packed.directory, { recursive: true, force: true }); });
  const transportIndex = JSON.parse(await readFile(join(packed.directory, "transport-index.json")));
  const schemas = JSON.parse(await readFile(join(source.directory, `releases/${source.release_id}/human/downloads/schema.json`)));
  const validators = Object.fromEntries(Object.keys(schemas.$defs).map((name) => [name, z.fromJSONSchema({ ...schemas, $ref: `#/$defs/${name}` })]));
  const rowsByPath = new Map(), summaries = new Map(), assets = new Map();
  for (const species of ["human", "mouse"]) {
    const index = JSON.parse(await readFile(join(source.directory, source.manifest.species[species].index_path))).data;
    for (const row of index) {
      rowsByPath.set(row.summary_path, row);
      summaries.set(row.summary_path, JSON.parse(await readFile(join(source.directory, row.summary_path))).data.summary);
    }
  }
  const info = { manifest: source.manifest, schemas, validators, rowsByPath };
  const releases = new Map([[source.release_id, info]]);
  const urlFor = (part) => `https://openantigens.org/api-runtime/${part.path.slice(6, -8)}-${part.sha256}.oap`;
  for (const parts of Object.values(transportIndex.targets)) for (const part of parts) assets.set(urlFor(part), await readFile(join(packed.directory, part.path)));
  const calls = [];
  let behavior = (url) => response(assets.get(url));
  const fetchImpl = async (url, options) => {
    calls.push({ url, options });
    assert.ok(assets.has(url), "request URL must come from a sealed descriptor");
    assert.equal(options.method, "GET"); assert.equal(options.redirect, "error"); assert.equal(options.credentials, "omit");
    assert.deepEqual(options.headers, { "accept-encoding": "identity", accept: "application/octet-stream", "user-agent": "OpenAntigens-Public-API/1" });
    return behavior(url, options);
  };
  const target = source.manifest.species.human.target_paths[source.target_id];
  const mouseTarget = Object.values(source.manifest.species.mouse.target_paths)[0];
  const make = (options = {}) => createRemoteQueryStorage({ transportIndex, releases, readSummary: async (_, row) => summaries.get(row.summary_path), fetchImpl, ...options });
  const replace = (path, transform) => {
    const part = transportIndex.targets[path][0], oldUrl = urlFor(part);
    let decoded = gunzipSync(assets.get(oldUrl));
    const container = JSON.parse(decoded);
    const changed = transform(container);
    decoded = Buffer.isBuffer(changed) ? changed : Buffer.from(JSON.stringify(container) + "\n");
    const bytes = gzipSync(decoded);
    assets.delete(oldUrl);
    Object.assign(part, { bytes: bytes.length, decoded_bytes: decoded.length, sha256: sha(bytes), members: container.members.length, min_path: container.members[0][0], max_path: container.members.at(-1)[0] });
    assets.set(urlFor(part), bytes);
  };
  return { source, packed, info, target, mouseTarget, transportIndex, assets, calls, summaries, urlFor, make, replace, behavior: (fn) => { behavior = fn; } };
}

test("cold and warm target reads preserve public data and reuse only validated compressed bytes", async (t) => {
  const f = await setup(t), engine = f.make();
  await engine.withContext(async () => {
    const target = await engine.readJson(f.info, f.target);
    assert.deepEqual(target.summary, f.summaries.get(f.target));
    assert.equal((await engine.readJson(f.info, target.constructs.pages[0]))[0].construct_id, f.source.construct_id);
    const sequence = await engine.readJson(f.info, target.construct_sequence_paths[f.source.construct_id], f.target);
    assert.equal(sequence.sequence_sha256, sha(sequence.sequence));
  });
  assert.equal(f.calls.length, 1);
  f.behavior(() => { throw new Error("offline"); });
  assert.equal((await engine.readJson(f.info, f.target)).summary.target_id, f.source.target_id);
  assert.equal(f.calls.length, 1);
  assert.equal(engine.snapshot().cacheEntries, 1);
  assert.equal(engine.snapshot().activeTargets, 0);
});

test("redirects, HTML, content encoding, length, status and checksum failures are nonsticky", async (t) => {
  const f = await setup(t);
  for (const bad of [
    (bytes) => response(bytes, { location: "https://other.example" }, 302),
    (bytes) => response(bytes, { "content-type": "text/html" }),
    (bytes) => response(bytes, { "content-encoding": "gzip" }),
    (bytes) => response(bytes, { "content-length": "1" }),
    () => response(Buffer.from("private failure"), {}, 404),
    (bytes) => response(Buffer.concat([bytes, Buffer.from("extra")])),
    (bytes) => response(Buffer.from(bytes).fill(0, 0, 1)),
    (bytes) => ({ status: 200, redirected: true, headers: new Headers({ "content-type": "application/octet-stream" }), body: response(bytes).body }),
  ]) {
    const engine = f.make();
    f.behavior((url) => bad(f.assets.get(url)));
    await assert.rejects(engine.readJson(f.info, f.target), fail503);
    assert.equal(engine.snapshot().cacheEntries, 0);
    f.behavior((url) => response(f.assets.get(url)));
    await engine.readJson(f.info, f.target);
    assert.equal(engine.snapshot().cacheEntries, 1);
  }
});

test("stream overflow cancels the body before accepting or caching it", async (t) => {
  const f = await setup(t), engine = f.make(); let cancelled = false;
  f.behavior((url) => new Response(new ReadableStream({
    start(controller) { controller.enqueue(f.assets.get(url)); controller.enqueue(new Uint8Array([1])); },
    cancel() { cancelled = true; },
  }), { headers: { "content-type": "application/octet-stream" } }));
  await assert.rejects(engine.readJson(f.info, f.target), fail503);
  assert.equal(cancelled, true); assert.equal(engine.snapshot().cacheEntries, 0);
});

test("target summary must equal the locally sealed summary", async (t) => {
  const f = await setup(t);
  f.replace(f.target, (container) => {
    const member = container.members.find(([path]) => path === f.target), value = JSON.parse(Buffer.from(member[1], "base64"));
    value.data.summary.gene_symbol = "CHANGED"; member[1] = Buffer.from(JSON.stringify(value)).toString("base64");
  });
  await assert.rejects(f.make().readJson(f.info, f.target), fail503);
});

test("duplicate public JSON, noncanonical base64 and foreign ownership fail closed", async (t) => {
  for (const kind of ["duplicate", "base64", "foreign", "missing", "hash"]) {
    const f = await setup(t);
    f.replace(f.target, (container) => {
      if (kind === "duplicate") {
        const member = container.members.find(([path]) => path === f.target);
        const text = Buffer.from(member[1], "base64").toString("utf8");
        member[1] = Buffer.from(text.replace('{"data":', '{"data":"PRIVATE_CANARY","data":')).toString("base64");
      } else if (kind === "base64") container.members[0][1] += "\n";
      else if (kind === "foreign") container.members[0][0] = container.members[0][0].replace("/human/", "/mouse/");
      else if (kind === "missing") container.members.shift();
      else {
        const member = container.members.find(([path]) => path.includes("/sequences/")), value = JSON.parse(Buffer.from(member[1], "base64"));
        value.data.sequence = "A".repeat(value.data.sequence.length);
        member[1] = Buffer.from(JSON.stringify(value)).toString("base64");
      }
    });
    const engine = f.make();
    await assert.rejects(engine.readJson(f.info, f.target), fail503);
    assert.equal(engine.snapshot().cacheEntries, 0);
  }
});

test("one request context cannot fetch a second target", async (t) => {
  const f = await setup(t), engine = f.make();
  await assert.rejects(engine.withContext(async () => { await engine.readJson(f.info, f.target); await engine.readJson(f.info, f.mouseTarget); }), fail503);
  assert.equal(f.calls.length, 1);
});

test("two admitted targets have independent cancellation and a third has no queue", async (t) => {
  const f = await setup(t), engine = f.make(), pending = [];
  f.behavior((url, options) => new Promise((resolve, reject) => {
    pending.push(() => resolve(response(f.assets.get(url))));
    options.signal.addEventListener("abort", () => reject(new Error("aborted")), { once: true });
  }));
  const firstAbort = new AbortController();
  const first = engine.withContext(() => engine.readJson(f.info, f.target), { signal: firstAbort.signal });
  const firstRejected = assert.rejects(first, fail503);
  const second = engine.withContext(() => engine.readJson(f.info, f.mouseTarget));
  await nextTurn(); assert.equal(f.calls.length, 2); assert.equal(engine.snapshot().activeTargets, 2);
  await assert.rejects(engine.readJson(f.info, f.target), fail503);
  assert.equal(f.calls.length, 2);
  firstAbort.abort(); await firstRejected;
  assert.equal(engine.snapshot().activeTargets, 1);
  pending[1](); await second;
  assert.equal(engine.snapshot().activeTargets, 0); assert.equal(engine.snapshot().activeFetches, 0);
  assert.equal(engine.snapshot().cacheEntries, 1);
});

test("aborted and expired contexts cannot populate the cache", async (t) => {
  const f = await setup(t), controller = new AbortController(); controller.abort();
  let time = 0; const engine = f.make({ now: () => time });
  await assert.rejects(engine.withContext(() => engine.readJson(f.info, f.target), { signal: controller.signal }), fail503);
  assert.equal(f.calls.length, 0);
  f.behavior((url) => { time = 3001; return response(f.assets.get(url)); });
  await assert.rejects(engine.readJson(f.info, f.target), fail503);
  assert.equal(engine.snapshot().cacheEntries, 0);
});

test("cancelling a stalled response body releases both admission slots", async (t) => {
  const f = await setup(t), engine = f.make(), controller = new AbortController(); let cancelled = false;
  f.behavior(() => new Response(new ReadableStream({ cancel() { cancelled = true; } }), { headers: { "content-type": "application/octet-stream" } }));
  const request = engine.withContext(() => engine.readJson(f.info, f.target), { signal: controller.signal });
  const rejected = assert.rejects(request, fail503);
  await nextTurn(); controller.abort(); await rejected;
  assert.equal(cancelled, true);
  assert.equal(engine.snapshot().activeTargets, 0); assert.equal(engine.snapshot().activeFetches, 0);
});

test("caller errors and local summary integrity failures retain their identity", async (t) => {
  const f = await setup(t), callerError = new ApiError(404, "not_found", "Missing record");
  await assert.rejects(f.make().withContext(() => { throw callerError; }), (error) => error === callerError);
  const localError = Object.assign(new Error("Local integrity failure"), { localQueryIntegrity: true });
  const engine = f.make({ readSummary: async () => { throw localError; } });
  await assert.rejects(engine.readJson(f.info, f.target), (error) => error === localError);
  assert.equal(engine.snapshot().cacheEntries, 0);
});

test("warm decode work also uses the two target admission slots", async (t) => {
  const f = await setup(t), engine = f.make();
  await engine.readJson(f.info, f.target);
  const gates = [];
  const hold = () => engine.withContext(async () => { await engine.readJson(f.info, f.target); await new Promise((resolve) => gates.push(resolve)); });
  const first = hold(), second = hold();
  while (gates.length < 2) await nextTurn();
  await assert.rejects(engine.readJson(f.info, f.target), fail503);
  assert.equal(f.calls.length, 1);
  gates.forEach((resolve) => resolve()); await Promise.all([first, second]);
  assert.equal(engine.snapshot().activeTargets, 0);
});

test("gzip corruption and decoded overflow never enter the shared cache", async (t) => {
  for (const bytes of [Buffer.from("invalid gzip"), gzipSync(Buffer.alloc(16 * MiB + 1, 32))]) {
    const f = await setup(t), part = f.transportIndex.targets[f.target][0];
    f.assets.delete(f.urlFor(part)); part.bytes = bytes.length; part.sha256 = sha(bytes); f.assets.set(f.urlFor(part), bytes);
    const engine = f.make();
    await assert.rejects(engine.readJson(f.info, f.target), fail503);
    assert.equal(engine.snapshot().cacheEntries, 0);
  }
});

test("target aggregate and part limits reject before any network request", async (t) => {
  for (const mutate of [
    (parts) => { parts[0].bytes = 4 * MiB + 1; },
    (parts) => { parts[0].decoded_bytes = 16 * MiB + 1; },
    (parts) => { parts.push({ ...parts[0] }, { ...parts[0] }); },
  ]) {
    const f = await setup(t); mutate(f.transportIndex.targets[f.target]);
    await assert.rejects(f.make().readJson(f.info, f.target), fail503);
    assert.equal(f.calls.length, 0);
  }
});

test("split groups fetch the target part first and only needed additional parts", async (t) => {
  const f = await setup(t), original = f.transportIndex.targets[f.target][0];
  const members = JSON.parse(gunzipSync(f.assets.get(f.urlFor(original)))).members;
  f.assets.delete(f.urlFor(original));
  const split = [members.slice(0, -1), members.slice(-1)];
  f.transportIndex.targets[f.target] = split.map((entries, i) => {
    const decoded = Buffer.from(JSON.stringify({ format: "openantigens-target-gzip-v1", members: entries }) + "\n"), bytes = gzipSync(decoded);
    const part = { ...original, path: original.path.replace("-0.json.gz", `-${i}.json.gz`), bytes: bytes.length, sha256: sha(bytes), decoded_bytes: decoded.length, members: entries.length, min_path: entries[0][0], max_path: entries.at(-1)[0] };
    f.assets.set(f.urlFor(part), bytes); return part;
  });
  const engine = f.make();
  f.behavior((url) => response(url.includes("-0-") ? Buffer.from("broken") : f.assets.get(url)));
  await assert.rejects(engine.withContext(async () => {
    const target = await engine.readJson(f.info, f.target);
    assert.equal(f.calls.length, 1); assert.ok(f.calls[0].url.includes("-1-"));
    await engine.readJson(f.info, target.reference_sequence_path);
  }), fail503);
  assert.equal(engine.snapshot().cacheEntries, 1);
  f.behavior((url) => response(f.assets.get(url)));
  await engine.withContext(async () => {
    const target = await engine.readJson(f.info, f.target);
    assert.equal(f.calls.length, 2);
    await engine.readJson(f.info, target.reference_sequence_path);
    assert.equal(f.calls.length, 3);
    await engine.readJson(f.info, target.constructs.pages[0]);
    assert.equal(f.calls.length, 3);
  });
  assert.equal(engine.snapshot().cacheEntries, 2);
});

test("shared cache evicts by entry count and never retains target objects globally", async (t) => {
  const f = await setup(t, { extraUnresolved: 128 }), engine = f.make();
  const targets = [...f.info.rowsByPath.keys()].filter((path) => path.includes("/human/")).slice(0, 129);
  for (const path of targets) await engine.readJson(f.info, path);
  assert.equal(engine.snapshot().cacheEntries, 128);
  assert.ok(engine.snapshot().cacheBytes <= 8 * MiB);
  const before = f.calls.length; await engine.readJson(f.info, targets[0]);
  assert.equal(f.calls.length, before + 1);
});

test("shared cache also evicts by its eight MiB byte bound", async (t) => {
  const f = await setup(t), targets = [...f.info.rowsByPath.keys()].slice(0, 3);
  for (const path of targets) {
    const part = f.transportIndex.targets[path][0], oldUrl = f.urlFor(part), original = f.assets.get(oldUrl);
    const header = Buffer.from(original.subarray(0, 10)); header[3] |= 16;
    const bytes = Buffer.concat([header, Buffer.alloc(3 * MiB, 65), Buffer.from([0]), original.subarray(10)]);
    f.assets.delete(oldUrl); part.bytes = bytes.length; part.sha256 = sha(bytes); f.assets.set(f.urlFor(part), bytes);
  }
  const engine = f.make();
  for (const target of targets) await engine.readJson(f.info, target);
  assert.equal(engine.snapshot().cacheEntries, 2);
  assert.ok(engine.snapshot().cacheBytes <= 8 * MiB);
});

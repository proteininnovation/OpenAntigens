import test from "node:test";
import { createHash } from "node:crypto";
import { gzipSync, gunzipSync } from "node:zlib";
import { request } from "node:http";
import assert from "node:assert/strict";
import {
  readFile,
  writeFile,
  rm,
  symlink,
  cp,
  truncate,
} from "node:fs/promises";
import { join } from "node:path";
import {
  Client,
  StreamableHTTPClientTransport,
} from "@modelcontextprotocol/client";
import { Provider, semanticHash, fasta } from "../src/provider.js";
import { makeServer } from "../src/server.js";
import { LIMITS } from "../src/contracts.js";
import { fixture } from "./fixture.js";
import { z } from "zod";
const decodedFasta = (text) => {
  const [header, ...lines] = text.trimEnd().split("\n");
  return {
    id: decodeURIComponent(header.slice(1).split(" ")[0]),
    fields: Object.fromEntries(
      header
        .slice(1)
        .split(" ")
        .slice(1)
        .map((field) => {
          const equal = field.indexOf("=");
          return [
            field.slice(0, equal),
            decodeURIComponent(field.slice(equal + 1)),
          ];
        }),
    ),
    sequence: lines.join(""),
  };
};
const meta = {
  "io.modelcontextprotocol/protocolVersion": "2026-07-28",
  "io.modelcontextprotocol/clientCapabilities": {},
  "io.modelcontextprotocol/clientInfo": { name: "wire-test", version: "1" },
};
async function setup(t, config = {}) {
  const f = await fixture();
  const provider = await new Provider(f.directory).initialize();
  const server = makeServer(provider, {
    clientRate: 1000,
    globalRate: 10000,
    ...config,
  });
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  t.after(async () => {
    await new Promise((r) => {
      server.closeAllConnections();
      server.close(r);
    });
    await rm(f.directory, { recursive: true, force: true });
  });
  const base = `http://127.0.0.1:${server.address().port}`;
  const root = `${base}/v1/releases/${f.release_id}`;
  const target = `${root}/human/targets/${f.target_id}`;
  return { ...f, provider, server, base, root, target };
}
async function wire(base, method, params = {}, headers = {}, modern = true) {
  const body = {
    jsonrpc: "2.0",
    id: 1,
    method,
    params: { ...params, ...(modern ? { _meta: meta } : {}) },
  };
  return fetch(`${base}/mcp`, {
    method: "POST",
    headers: {
      "content-type": "application/json",
      accept: "application/json, text/event-stream",
      ...(modern
        ? {
            "mcp-protocol-version": "2026-07-28",
            "mcp-method": method,
            ...(params.name ? { "mcp-name": params.name } : {}),
          }
        : {}),
      ...headers,
    },
    body: JSON.stringify(body),
  });
}
async function protocolBody(response) {
  const text = await response.text();
  return text.startsWith("event:") || text.startsWith("data:")
    ? JSON.parse(
        text
          .split("\n")
          .find((v) => v.startsWith("data:"))
          .slice(5),
      )
    : JSON.parse(text);
}
test("HTTP routes preserve species, identity, filters, cursors, exact sequences and cache semantics", async (t) => {
  const f = await setup(t);
  let r = await fetch(`${f.base}/v1/releases/current`);
  assert.equal(r.status, 200);
  assert.equal(r.headers.get("cache-control"), "no-store");
  assert.equal((await r.json()).current_release_id, f.release_id);
  r = await fetch(`${f.root}/human/targets?limit=1`);
  const first = await r.json();
  assert.equal(first.total, 2);
  assert.equal(first.data[0].target_id, f.target_id);
  assert.ok(first.next_cursor);
  r = await fetch(
    `${f.root}/human/targets?limit=1&cursor=${first.next_cursor}`,
  );
  assert.equal((await r.json()).data[0].target_id, null);
  for (const args of [
    "q=ALIAS",
    "track=secreted&status=ok&structure=unavailable&constructs=available",
    "disease_id=EFO%3A0001",
  ]) {
    r = await fetch(`${f.root}/human/targets?${args}`);
    assert.equal(r.status, 200);
    assert.ok((await r.json()).total >= 1);
  }
  assert.equal((await fetch(f.root)).headers.get("cache-control"), "no-store");
  r = await fetch(f.target);
  assert.equal(r.status, 200);
  const etag = r.headers.get("etag");
  assert.ok(etag);
  assert.equal((await r.json()).data.sequence_accession, "P12345");
  assert.equal(
    (await fetch(f.target, { headers: { "if-none-match": etag } })).status,
    304,
  );
  r = await fetch(f.target, { method: "HEAD" });
  assert.equal(r.status, 200);
  assert.equal(await r.text(), "");
  assert.equal(r.headers.get("etag"), etag);
  r = await fetch(`${f.target}/sections/topology`);
  assert.equal((await r.json()).data[0].extracellular_regions[0].start, 1);
  r = await fetch(`${f.target}/constructs?category=full_region`);
  assert.equal((await r.json()).total, 1);
  r = await fetch(`${f.target}/constructs/${f.construct_id}/sequence`);
  const seq = (await r.json()).data;
  assert.equal(seq.sequence, "ACDEFGHIKLMNPQRSTVWY");
  r = await fetch(`${f.target}/sequence?format=fasta`, {
    headers: { accept: "text/x-fasta" },
  });
  const saved = decodedFasta(await r.text());
  assert.equal(saved.id, `${f.release_id}:${f.target_id}:reference`);
  assert.equal(saved.fields.release_id, f.release_id);
  assert.equal(saved.fields.target_id, f.target_id);
  assert.equal(saved.fields.construct_id, "null");
  assert.equal(saved.fields.sequence_source, seq.sequence_source);
  assert.equal(saved.fields.sequence_basis, "reference");
  assert.deepEqual(JSON.parse(saved.fields.warnings), seq.warnings);
  assert.deepEqual(JSON.parse(saved.fields.evidence), seq.evidence);
  assert.deepEqual(JSON.parse(saved.fields.record_selector), {
    target_id: f.target_id,
    construct_id: null,
  });
  assert.equal(saved.sequence, seq.sequence);
  assert.match(
    saved.fields.provenance_url,
    new RegExp(`/data/releases/${f.release_id}/manifest.json$`),
  );
  assert.match(r.headers.get("link"), /rel="describedby"/);
  assert.match(
    r.headers.get("content-disposition"),
    /^attachment; filename="openantigens-/,
  );
  r = await fetch(`${f.root}/mouse/targets`);
  assert.equal((await r.json()).data[0].taxon_id, 10090);
  r = await fetch(`${f.root}/human/downloads`);
  assert.equal((await r.json()).data[0].name, "records.jsonl.gz");
  for (const [path, status] of [
    [`${f.root}/human/targets?limit=101`, 400],
    [`${f.root}/human/targets?q=a&q=b`, 400],
    [`${f.target}?q=a`, 400],
    [`${f.root}/human/targets?cursor=${first.next_cursor}&q=changed`, 400],
    [`${f.root}/mouse/targets/${f.target_id}`, 400],
    [`${f.root}/mouse/targets?disease_id=EFO:1`, 400],
    [`${f.root}/mouse/targets/10090:uniprot:Q12345/sections/disease`, 400],
    [`${f.base}/v1/releases/current/human/targets`, 400],
    [`${f.target}/sequence?format=xml`, 400],
  ])
    assert.equal((await fetch(path)).status, status, path);
});
test("HTTP rejects untrusted origins, hosts, bodies, media and methods; exact CORS and request IDs", async (t) => {
  const f = await setup(t);
  await new Promise((resolve, reject) => {
    const req = request(
      f.target,
      { headers: { host: "attacker.invalid" } },
      (res) => {
        assert.equal(res.statusCode, 403);
        res.resume();
        res.on("end", resolve);
      },
    );
    req.on("error", reject);
    req.end();
  });
  for (const [headers, status] of [
    [{ origin: "https://attacker.invalid" }, 403],
    [{ origin: "null" }, 403],
    [{ accept: "application/xml" }, 406],
    [{ "content-encoding": "gzip" }, 415],
  ]) {
    const r = await fetch(f.target, { headers });
    assert.equal(r.status, status, JSON.stringify(headers));
    assert.ok((await r.json()).error.request_id);
    assert.equal(r.headers.get("cache-control"), "no-store");
  }
  let r = await fetch(f.target, { method: "POST", body: "x" });
  assert.equal(r.status, 405);
  r = await fetch(`${f.base}/mcp`, { method: "DELETE" });
  assert.equal(r.status, 405);
  r = await fetch(`${f.base}/mcp`);
  assert.equal(r.status, 405);
  r = await fetch(`${f.base}/mcp`, {
    method: "OPTIONS",
    headers: {
      origin: "https://openantigens.org",
      "access-control-request-method": "POST",
      "access-control-request-headers":
        "content-type, mcp-protocol-version, mcp-method, mcp-name",
    },
  });
  assert.equal(r.status, 204);
  assert.equal(
    r.headers.get("access-control-allow-origin"),
    "https://openantigens.org",
  );
  assert.ok(!r.headers.has("access-control-allow-credentials"));
  r = await fetch(`${f.base}/mcp`, {
    method: "OPTIONS",
    headers: {
      origin: "https://openantigens.org",
      "access-control-request-method": "POST",
      "access-control-request-headers": "authorization",
    },
  });
  assert.equal(r.status, 400);
  r = await fetch(`${f.base}/mcp`, {
    method: "POST",
    headers: { "content-type": "application/json" },
    body: "x".repeat(32769),
  });
  assert.equal(r.status, 413);
  assert.equal((await fetch(`${f.target}?${"a".repeat(4096)}`)).status, 413);
});
test("rate and concurrency admission limits return Retry-After and do not trust forwarded IPs", async (t) => {
  const f = await setup(t, { clientRate: 1 });
  assert.equal((await fetch(`${f.base}/healthz`)).status, 200);
  const r = await fetch(`${f.base}/healthz`, {
    headers: { "x-forwarded-for": "1.2.3.4" },
  });
  assert.equal(r.status, 429);
  assert.equal(r.headers.get("retry-after"), "60");
});
test("readiness fails closed; changed bytes and symlinks cannot become data", async (t) => {
  const f = await setup(t);
  const path = f.manifest.species.human.index_path;
  const original = await readFile(join(f.directory, path));
  await writeFile(join(f.directory, path), Buffer.from(original).fill(32));
  await assert.rejects(new Provider(f.directory).initialize());
  await writeFile(join(f.directory, path), original);
  await rm(join(f.directory, path));
  await symlink("/etc/passwd", join(f.directory, path));
  await assert.rejects(new Provider(f.directory).initialize());
  f.provider.ready = false;
  assert.equal((await fetch(`${f.base}/readyz`)).status, 503);
  assert.equal((await fetch(f.target)).status, 503);
  assert.equal((await fetch(`${f.base}/healthz`)).status, 200);
});
test("official MCP clients negotiate current and legacy, list six read-only tools and call provider", async (t) => {
  const f = await setup(t);
  for (const mode of ["auto", "legacy"]) {
    const client = new Client(
      { name: "test", version: "1" },
      { versionNegotiation: { mode }, capabilities: {} },
    );
    const transport = new StreamableHTTPClientTransport(
      new URL(`${f.base}/mcp`),
    );
    await client.connect(transport);
    try {
      const tools = await client.listTools();
      assert.equal(tools.tools.length, 6);
      for (const tool of tools.tools) {
        assert.equal(tool.annotations.readOnlyHint, true);
        assert.equal(tool.annotations.openWorldHint, false);
        assert.equal(tool.inputSchema.additionalProperties, false);
        assert.ok(tool.outputSchema);
      }
      const result = await client.callTool({
        name: "get_sequence",
        arguments: {
          release_id: f.release_id,
          species: "human",
          target_id: f.target_id,
        },
      });
      assert.equal(result.isError, undefined);
      assert.equal(result.structuredContent.data.sequence_accession, "P12345");
      assert.equal(transport.sessionId, undefined);
      const failure = await client.callTool({
        name: "get_target",
        arguments: {
          release_id: f.release_id,
          species: "mouse",
          target_id: f.target_id,
        },
      });
      assert.equal(failure.isError, true);
      assert.match(failure.content[0].text, /species_mismatch/);
    } finally {
      await client.close();
    }
  }
});
test("independent current and legacy MCP wire flow, protocol errors, notifications, batches and prohibited methods", async (t) => {
  const f = await setup(t);
  let r = await wire(f.base, "server/discover");
  assert.equal(r.status, 200);
  assert.match(r.headers.get("content-type"), /application\/json/);
  assert.equal(r.headers.has("mcp-session-id"), false);
  let result = await protocolBody(r);
  assert.ok(result.result.capabilities.tools);
  r = await wire(f.base, "tools/list");
  result = await protocolBody(r);
  assert.equal(result.result.tools.length, 6);
  r = await wire(f.base, "tools/call", {
    name: "get_downloads",
    arguments: { release_id: f.release_id, species: "mouse" },
  });
  result = await protocolBody(r);
  assert.equal(result.result.structuredContent.species, "mouse");
  r = await wire(
    f.base,
    "initialize",
    {
      protocolVersion: "2025-11-25",
      capabilities: {},
      clientInfo: { name: "wire-legacy", version: "1" },
    },
    {},
    false,
  );
  assert.equal(r.status, 200);
  assert.match(r.headers.get("content-type"), /text\/event-stream/);
  result = await protocolBody(r);
  assert.equal(result.result.protocolVersion, "2025-11-25");
  assert.equal(r.headers.has("mcp-session-id"), false);
  for (const [method, params, headers] of [
    ["tools/list", {}, { "mcp-protocol-version": "2099-01-01" }],
    ["tools/list", {}, { "mcp-method": "tools/call" }],
    [
      "tools/call",
      { name: "get_release", arguments: {} },
      { "mcp-name": "get_sequence" },
    ],
    ["resources/list", {}, {}],
    ["sampling/createMessage", {}, {}],
    ["subscriptions/listen", {}, {}],
  ]) {
    r = await wire(f.base, method, params, headers);
    result = await protocolBody(r);
    assert.ok(result.error, method);
  }
  r = await fetch(`${f.base}/mcp`, {
    method: "POST",
    headers: {
      "content-type": "application/json",
      accept: "application/json, text/event-stream",
    },
    body: "{",
  });
  assert.ok(r.status >= 400);
  r = await fetch(`${f.base}/mcp`, {
    method: "POST",
    headers: {
      "content-type": "application/json",
      accept: "application/json, text/event-stream",
    },
    body: JSON.stringify([{ jsonrpc: "2.0", id: 1, method: "tools/list" }]),
  });
  assert.ok((await protocolBody(r)).error);
  r = await fetch(`${f.base}/mcp`, {
    method: "POST",
    headers: {
      "content-type": "application/json",
      accept: "application/json, text/event-stream",
    },
    body: JSON.stringify({
      jsonrpc: "2.0",
      method: "notifications/initialized",
    }),
  });
  assert.equal(r.status, 202);
  assert.equal(await r.text(), "");
});

test("maintenance blocks discovery while literal release reads remain available", async (t) => {
  const f = await setup(t);
  f.provider.publicationState = "maintenance";
  for (const path of ["/v1/releases", "/v1/releases/current"]) {
    const response = await fetch(f.base + path);
    assert.equal(response.status, 503);
    assert.equal(response.headers.get("retry-after"), "60");
    assert.equal((await response.json()).error.code, "publication_in_progress");
  }
  const pinned = await fetch(f.root);
  assert.equal(pinned.status, 200);
  assert.equal(pinned.headers.get("cache-control"), "no-store");
  const metadata = await pinned.json();
  assert.equal(metadata.discovery_state, "publication_in_progress");
  for (const field of ["current_release_id", "releases", "retired_releases"])
    assert.equal(Object.hasOwn(metadata, field), false, field);
  const protocol = await protocolBody(
    await wire(f.base, "tools/call", {
      name: "get_release",
      arguments: { release_id: f.release_id },
    }),
  );
  assert.deepEqual(protocol.result.structuredContent, metadata);
  const omitted = await protocolBody(
    await wire(f.base, "tools/call", { name: "get_release", arguments: {} }),
  );
  assert.equal(omitted.result.isError, true);
  assert.equal(
    JSON.parse(omitted.result.content[0].text).code,
    "publication_in_progress",
  );
  assert.equal((await fetch(f.target)).status, 200);
});
test("global and concurrent quotas, deadlines and explicit overflow fallback", async (t) => {
  const global = await setup(t, { globalRate: 1 });
  assert.equal((await fetch(global.base + "/healthz")).status, 200);
  assert.equal((await fetch(global.base + "/healthz")).status, 429);
  const f = await setup(t, { concurrency: 1, handlerMs: 30 });
  const call = f.provider.call.bind(f.provider);
  let started;
  const signal = new Promise((r) => {
    started = r;
  });
  f.provider.call = async (...args) => {
    started();
    await new Promise((r) => setTimeout(r, 80));
    return call(...args);
  };
  const slow = fetch(f.target).catch(() => null);
  await signal;
  assert.equal((await fetch(f.target)).status, 429);
  const timeout = await slow;
  assert.equal(timeout.status, 503);
  assert.equal((await timeout.json()).error.code, "handler_timeout");
  await new Promise((r) => setTimeout(r, 90));
  f.provider.call = async () => ({
    ...(await call("get_target", {
      release_id: f.release_id,
      species: "human",
      target_id: f.target_id,
    })),
    data: { oversized: "x".repeat(262145) },
  });
  const overflow = await fetch(f.target);
  assert.equal(overflow.status, 413);
  const error = (await overflow.json()).error;
  assert.equal(error.code, "result_too_large");
  assert.equal(error.download.name, "records.jsonl.gz");
  assert.ok(!JSON.stringify(error).includes(f.directory));
});
test(
  "actual Python public bundle interoperates with both HTTP and official MCP clients",
  { skip: !process.env.OPENANTIGENS_TEST_DATA_DIR },
  async (t) => {
    const provider = await new Provider(
      process.env.OPENANTIGENS_TEST_DATA_DIR,
    ).initialize();
    const server = makeServer(provider, {
      clientRate: 1000,
      globalRate: 10000,
    });
    await new Promise((r) => server.listen(0, "127.0.0.1", r));
    t.after(
      () =>
        new Promise((r) => {
          server.closeAllConnections();
          server.close(r);
        }),
    );
    const base = `http://127.0.0.1:${server.address().port}`;
    const client = new Client(
      { name: "crosslanguage", version: "1" },
      { versionNegotiation: { mode: "auto" } },
    );
    await client.connect(
      new StreamableHTTPClientTransport(new URL(base + "/mcp")),
    );
    t.after(() => client.close());
    for (const species of ["human", "mouse"]) {
      const args = { release_id: provider.current, species };
      const http = await (
        await fetch(
          `${base}/v1/releases/${provider.current}/${species}/targets`,
        )
      ).json();
      const mcp = await client.callTool({
        name: "search_targets",
        arguments: args,
      });
      assert.deepEqual(mcp.structuredContent, http);
      const target = http.data.find((v) => v.target_id);
      assert.ok(target);
      const params = { ...args, target_id: target.target_id };
      for (const name of [
        "get_target",
        "list_constructs",
        "get_sequence",
        "get_downloads",
      ]) {
        const result = await client.callTool({
          name,
          arguments: name === "get_downloads" ? args : params,
        });
        assert.notEqual(result.isError, true, name);
        assert.deepEqual(
          result.structuredContent,
          await provider.call(name, name === "get_downloads" ? args : params),
        );
      }
      const downloads = (await provider.call("get_downloads", args)).data;
      const constructs = (await provider.call("list_constructs", params)).data;
      for (const name of ["reference.fasta", "constructs.fasta"]) {
        const download = downloads.find((d) => d.name === name);
        const bulk = await readFile(
          join(process.env.OPENANTIGENS_TEST_DATA_DIR, download.path),
          "utf8",
        );
        const identifiers = [...bulk.matchAll(/^>(\S+)/gm)].map(
          (match) => match[1],
        );
        assert.equal(
          new Set(identifiers).size,
          identifiers.length,
          `${species} ${name} unique primary identifiers`,
        );
      }
      for (const construct_id of [
        null,
        ...constructs
          .filter((c) => c.export_state === "available")
          .map((c) => c.construct_id),
      ]) {
        const sequence = (
          await provider.call("get_sequence", {
            ...params,
            ...(construct_id ? { construct_id } : {}),
          })
        ).data;
        const route = `${base}/v1/releases/${provider.current}/${species}/targets/${target.target_id}/${construct_id ? `constructs/${construct_id}/` : ""}sequence?format=fasta`;
        const response = await fetch(route);
        assert.equal(response.status, 200);
        const text = await response.text();
        const saved = decodedFasta(text);
        assert.equal(saved.fields.sequence_source, sequence.sequence_source);
        assert.equal(
          saved.fields.sequence_accession,
          String(sequence.sequence_accession),
        );
        assert.equal(saved.fields.sequence_basis, sequence.sequence_basis);
        assert.deepEqual(JSON.parse(saved.fields.warnings), sequence.warnings);
        assert.deepEqual(JSON.parse(saved.fields.evidence), sequence.evidence);
        assert.equal(saved.sequence, sequence.sequence);
        const download = downloads.find(
          (d) =>
            d.name === (construct_id ? "constructs.fasta" : "reference.fasta"),
        );
        const bulk = await readFile(
          join(process.env.OPENANTIGENS_TEST_DATA_DIR, download.path),
          "utf8",
        );
        assert.ok(
          bulk.includes(text),
          `Exact ${species} bulk/API FASTA parity`,
        );
      }
      const detail = await provider.call("get_target", params);
      for (const section of Object.keys(detail.data.sections)) {
        if (
          provider
            .resolve(provider.current)
            .manifest.species[species].capabilities.sections.includes(section)
        ) {
          const result = await client.callTool({
            name: "get_target",
            arguments: { ...params, section },
          });
          assert.notEqual(result.isError, true, section);
        }
      }
    }
    const r = await fetch(base + "/openapi.json");
    assert.equal(r.status, 200);
    const description = await r.json();
    assert.ok(description.components.schemas.get_sequenceOutput);
    assert.match(
      description.components.schemas.get_sequenceOutput.properties.data
        .properties.sequence_sha256.pattern,
      /a-f/,
    );
  },
);

test("portable numeric semantic hash matches the frozen Python vector", () => {
  assert.equal(
    semanticHash({
      integer: 50,
      decimal: 50.0,
      tiny: 1e-7,
      negative_zero: -0.0,
      boolean: false,
    }),
    "6541843bb6ae61c70aa65dfcb7aeb531d4a8736b3918b96aa24c930a4b1f9ede",
  );
});

test("release identity rejects rewritten shard checksums and corrupted rehashed fallback", async (t) => {
  const checksum = (bytes) => createHash("sha256").update(bytes).digest("hex");
  for (const fallback of [false, true]) {
    const f = await fixture();
    t.after(() => rm(f.directory, { recursive: true, force: true }));
    const rootPath = join(f.directory, "manifest.json");
    const root = JSON.parse(await readFile(rootPath));
    const manifestPath = join(f.directory, root.releases[0].manifest_path);
    const manifest = JSON.parse(await readFile(manifestPath));
    const path = fallback
      ? manifest.species.human.records_download.path
      : Object.keys(manifest.files).find((path) =>
          path.includes("/human/constructs/"),
        );
    let content;
    if (fallback) {
      const records = gunzipSync(await readFile(join(f.directory, path)))
        .toString()
        .trimEnd()
        .split("\n")
        .map((line) => JSON.parse(line));
      records[0].payload.data.items[0].source = "changed";
      content = gzipSync(
        Buffer.from(
          records.map((record) => JSON.stringify(record)).join("\n") + "\n",
        ),
        { mtime: 0 },
      );
    } else {
      const record = JSON.parse(await readFile(join(f.directory, path)));
      record.data[0].score = 50.0;
      content = Buffer.from(JSON.stringify(record) + "\n");
    }
    await writeFile(join(f.directory, path), content);
    Object.assign(manifest.files[path], {
      bytes: content.length,
      sha256: checksum(content),
    });
    if (fallback) {
      Object.assign(manifest.species.human.records_download, {
        bytes: content.length,
        sha256: checksum(content),
      });
      manifest.species.human.downloads[0] =
        manifest.species.human.records_download;
    }
    const manifestBytes = Buffer.from(JSON.stringify(manifest) + "\n");
    await writeFile(manifestPath, manifestBytes);
    root.releases[0].sha256 = checksum(manifestBytes);
    await writeFile(rootPath, JSON.stringify(root));
    await assert.rejects(
      new Provider(f.directory).initialize(),
      fallback ? /Fallback record disagreement/ : /Public payload identity/,
    );
  }
});

test("three release snapshots distinguish retired IDs, retained IDs and unknown IDs", async (t) => {
  const first = await fixture("first", { legacyRootSchema: true }),
    second = await fixture("second"),
    third = await fixture("third");
  for (const f of [first, second, third])
    t.after(() => rm(f.directory, { recursive: true, force: true }));
  const firstRoot = JSON.parse(
    await readFile(join(first.directory, "manifest.json")),
  );
  const secondRoot = JSON.parse(
    await readFile(join(second.directory, "manifest.json")),
  );
  await cp(
    join(first.directory, "releases", first.release_id),
    join(second.directory, "releases", first.release_id),
    { recursive: true },
  );
  secondRoot.releases.push(firstRoot.releases[0]);
  await writeFile(
    join(second.directory, "manifest.json"),
    JSON.stringify(secondRoot),
  );
  assert.equal(
    (await new Provider(second.directory).initialize()).releases.size,
    2,
  );
  const thirdRoot = JSON.parse(
    await readFile(join(third.directory, "manifest.json")),
  );
  await cp(
    join(second.directory, "releases", second.release_id),
    join(third.directory, "releases", second.release_id),
    { recursive: true },
  );
  thirdRoot.releases.push(secondRoot.releases[0]);
  thirdRoot.retired_releases = [
    {
      release_id: first.release_id,
      retired_by_release_id: third.release_id,
      archive_url: null,
    },
  ];
  await writeFile(
    join(third.directory, "manifest.json"),
    JSON.stringify(thirdRoot),
  );
  const provider = await new Provider(third.directory).initialize();
  const server = makeServer(provider, { clientRate: 1000, globalRate: 10000 });
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  t.after(
    () =>
      new Promise((r) => {
        server.closeAllConnections();
        server.close(r);
      }),
  );
  const base = `http://127.0.0.1:${server.address().port}`;
  for (const [id, status] of [
    [first.release_id, 410],
    [second.release_id, 200],
    [third.release_id, 200],
    ["f".repeat(64), 404],
  ]) {
    const response = await fetch(`${base}/v1/releases/${id}`);
    assert.equal(response.status, status);
    const protocol = await protocolBody(
      await wire(base, "tools/call", {
        name: "get_release",
        arguments: { release_id: id },
      }),
    );
    if (status === 410) {
      const error = (await response.json()).error;
      assert.equal(error.code, "retired_release");
      assert.deepEqual(error.retirement, thirdRoot.retired_releases[0]);
      assert.equal(
        JSON.parse(protocol.result.content[0].text).code,
        "retired_release",
      );
    } else if (status === 404)
      assert.equal(
        JSON.parse(protocol.result.content[0].text).code,
        "not_found",
      );
    else assert.notEqual(protocol.result.isError, true);
  }
  const current = await (await fetch(base + "/v1/releases/current")).json();
  assert.deepEqual(current.retired_releases, thirdRoot.retired_releases);
  delete thirdRoot.retired_releases;
  await writeFile(
    join(third.directory, "manifest.json"),
    JSON.stringify(thirdRoot),
  );
  assert.equal(
    (await new Provider(third.directory).initialize()).retired.size,
    0,
    "Legacy v1 roots remain readable",
  );
});
test("operational OpenAPI response schemas match actual bodies, status and cache headers", async (t) => {
  const f = await setup(t);
  const response = await fetch(f.base + "/openapi.json");
  const description = await response.json();
  assert.equal(
    description.externalDocs.url,
    "https://openantigens.org/agent-guide.html",
  );
  const validate = (path, status, body) => {
    const declaration =
      description.paths[path].get.responses[status].content["application/json"]
        .schema;
    const schema = { ...declaration, $defs: description.components.schemas };
    const json = JSON.parse(
      JSON.stringify(schema).replaceAll("#/components/schemas/", "#/$defs/"),
    );
    assert.equal(
      z.fromJSONSchema(json).safeParse(body).success,
      true,
      `${path} ${status}`,
    );
  };
  for (const path of ["/healthz", "/readyz", "/openapi.json"]) {
    const r = await fetch(f.base + path);
    assert.equal(r.status, 200);
    assert.equal(r.headers.get("cache-control"), "no-store");
    validate(path, 200, await r.json());
    assert.equal(
      Object.hasOwn(description.paths[path].get.responses, "304"),
      false,
    );
    assert.equal(
      description.paths[path].get.responses[200].headers.ETag,
      undefined,
    );
  }
  assert.equal(
    description.paths["/v1/releases/{release_id}"].get.responses[404].content[
      "application/json"
    ].schema.$ref,
    "#/components/schemas/Error",
  );
  assert.equal(
    description.paths["/healthz"].head.responses[200].content,
    undefined,
  );
  assert.equal(
    description.paths["/healthz"].head.responses[200].headers["Cache-Control"]
      .schema.const,
    "no-store",
  );
  f.provider.ready = false;
  const unavailable = await fetch(f.base + "/readyz");
  assert.equal(unavailable.status, 503);
  validate("/readyz", 503, await unavailable.json());
});

test("FASTA primary identifiers keep reference and multiple constructs distinct", async (t) => {
  const f = await fixture();
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const provider = await new Provider(f.directory).initialize();
  const envelope = await provider.call("get_sequence", {
    release_id: f.release_id,
    species: "human",
    target_id: f.target_id,
  });
  const options = {
    release_id: f.release_id,
    provenance_url: envelope.manifest.url,
    records_url: provider.fallback(f.release_id, "human").url,
  };
  const records = [
    envelope.data,
    ...["a", "b"].map((char) => ({
      ...envelope.data,
      construct_id: char.repeat(64),
      sequence_basis: "reference_slice",
    })),
  ];
  const texts = records.map((sequence) => fasta(sequence, options));
  assert.equal(new Set(texts.map((text) => text.split(/\s/)[0])).size, 3);
  for (const [position, text] of texts.entries())
    assert.equal(
      decodedFasta(text).id,
      `${f.release_id}:${f.target_id}:${records[position].construct_id ?? "reference"}`,
    );
});

test("measurement trial separates 64 MiB manifests and aggregate indexes from unchanged shard limits", async (t) => {
  const checksum = (bytes) => createHash("sha256").update(bytes).digest("hex");
  const f = await fixture("budget-trial");
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const rootPath = join(f.directory, "manifest.json");
  const root = JSON.parse(await readFile(rootPath));
  const manifestPath = join(f.directory, root.releases[0].manifest_path);
  const manifest = JSON.parse(await readFile(manifestPath));
  const writeManifest = async (padding) => {
    const json = Buffer.from(JSON.stringify(manifest) + "\n");
    const bytes = padding
      ? Buffer.concat([json, Buffer.alloc(padding - json.length, 32)])
      : json;
    await writeFile(manifestPath, bytes);
    root.releases[0].sha256 = checksum(bytes);
    await writeFile(rootPath, JSON.stringify(root));
  };
  await writeManifest(LIMITS.target + 1);
  assert.equal(
    (await new Provider(f.directory).initialize()).ready,
    true,
    "Release manifest may exceed target's 8 MiB cap",
  );
  await truncate(manifestPath, LIMITS.manifest + 1);
  await assert.rejects(
    new Provider(f.directory).initialize(),
    /Bundle size/,
    "Release manifest remains bounded at 64 MiB",
  );
  for (const species of ["human", "mouse"]) {
    const path = manifest.species[species].index_path;
    const original = await readFile(join(f.directory, path));
    const bytes = Buffer.concat([
      original,
      Buffer.alloc(LIMITS.index - original.length, 32),
    ]);
    await writeFile(join(f.directory, path), bytes);
    Object.assign(manifest.files[path], {
      bytes: bytes.length,
      sha256: checksum(bytes),
    });
  }
  await writeManifest();
  assert.equal(
    (await new Provider(f.directory).initialize()).ready,
    true,
    "Two 32 MiB indexes fit the aggregate 64 MiB trial cap",
  );
  const next = await fixture("budget-next");
  t.after(() => rm(next.directory, { recursive: true, force: true }));
  const nextRootPath = join(next.directory, "manifest.json");
  const nextRoot = JSON.parse(await readFile(nextRootPath));
  await cp(
    join(f.directory, "releases", f.release_id),
    join(next.directory, "releases", f.release_id),
    { recursive: true },
  );
  nextRoot.releases.push(root.releases[0]);
  await writeFile(nextRootPath, JSON.stringify(nextRoot));
  await assert.rejects(
    new Provider(next.directory).initialize(),
    /Index memory budget/,
    "The cap includes indexes from every retained release",
  );
  const targetPath = manifest.species.human.target_paths[f.target_id];
  await truncate(join(f.directory, targetPath), LIMITS.target + 1);
  manifest.files[targetPath].bytes = LIMITS.target + 1;
  await writeManifest();
  await assert.rejects(
    new Provider(f.directory).initialize(),
    /Bundle size/,
    "Target storage remains bounded at 8 MiB",
  );
  assert.equal(LIMITS.index, 32 * 1024 * 1024);
  assert.equal(LIMITS.shard, 1024 * 1024);
  assert.equal(LIMITS.manifest, 64 * 1024 * 1024);
  assert.equal(LIMITS.indexes, 64 * 1024 * 1024);
});

test("ten thousand input hashes stay in immutable manifests while HTTP and MCP responses stay bounded", async (t) => {
  const inputHashes = Object.fromEntries(
    Array.from({ length: 10001 }, (_, index) => [
      `input-record-${index}`,
      "a".repeat(64),
    ]),
  );
  const f = await fixture("large-input-provenance", { inputHashes });
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const provider = await new Provider(f.directory).initialize();
  const server = makeServer(provider, { clientRate: 1000, globalRate: 10000 });
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  t.after(
    () =>
      new Promise((r) => {
        server.closeAllConnections();
        server.close(r);
      }),
  );
  const base = `http://127.0.0.1:${server.address().port}`;
  const inspect = (result) => {
    const serialized = JSON.stringify(result);
    assert.ok(Buffer.byteLength(serialized) < LIMITS.json);
    assert.equal(serialized.includes('"input_sha256"'), false);
    assert.equal(result.provenance.input_hashes.count, 10001);
    assert.equal(
      result.provenance.input_hashes.json_pointer,
      "/provenance/input_sha256",
    );
    assert.equal(
      result.provenance.input_hashes.manifest_url,
      result.manifest.url,
    );
    assert.equal(
      result.provenance.input_hashes.manifest_sha256,
      result.manifest.sha256,
    );
  };
  for (const route of [
    "/v1/releases",
    "/v1/releases/current",
    `/v1/releases/${f.release_id}`,
    `/v1/releases/${f.release_id}/human/targets?limit=1`,
    `/v1/releases/${f.release_id}/human/targets/${f.target_id}`,
  ]) {
    const response = await fetch(base + route);
    assert.equal(response.status, 200, route);
    inspect(await response.json());
  }
  const client = new Client(
    { name: "large-provenance", version: "1" },
    { versionNegotiation: { mode: "auto" } },
  );
  await client.connect(
    new StreamableHTTPClientTransport(new URL(base + "/mcp")),
  );
  t.after(() => client.close());
  for (const [name, args] of [
    ["get_release", {}],
    [
      "search_targets",
      { release_id: f.release_id, species: "human", limit: 1 },
    ],
  ]) {
    const result = await client.callTool({ name, arguments: args });
    assert.notEqual(result.isError, true);
    inspect(result.structuredContent);
    assert.equal(
      provider.toolOutput(name).safeParse(result.structuredContent).success,
      true,
    );
  }
  const root = JSON.parse(await readFile(join(f.directory, "manifest.json")));
  const manifest = JSON.parse(
    await readFile(join(f.directory, root.releases[0].manifest_path)),
  );
  assert.equal(Object.keys(manifest.provenance.input_sha256).length, 10001);
  assert.deepEqual(manifest.provenance.input_sha256, inputHashes);
  assert.deepEqual(manifest.core.provenance.input_sha256, inputHashes);
  const metadata = await provider.call("get_release", {});
  for (const record of [metadata.data, ...metadata.releases]) {
    assert.equal(record.provenance.input_hashes.count, 10001);
    assert.equal(Object.hasOwn(record.provenance, "input_sha256"), false);
  }
});

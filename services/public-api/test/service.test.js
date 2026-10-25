import test from "node:test";
import { createHash } from "node:crypto";
import { gzipSync, gunzipSync } from "node:zlib";
import { request } from "node:http";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
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
import { makeServer, readBuildPins, publicationState } from "../src/server.js";
import { LIMITS } from "../src/contracts.js";
import { fixture, packedFixture } from "./fixture.js";
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
const testAdmission = {
  clientRate: 1000,
  clientBurst: 1000,
  globalRate: 10000,
  globalBurst: 10000,
  concurrency: 8,
  clientConcurrency: 8,
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
    ...testAdmission,
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
  r = await fetch(`${f.root}/human/targets?q=catalog&limit=1`);
  const first = await r.json();
  assert.equal(first.total, 2);
  assert.equal(first.data[0].target_id, f.target_id);
  assert.ok(first.next_cursor);
  r = await fetch(
    `${f.root}/human/targets?q=catalog&limit=1&cursor=${first.next_cursor}`,
  );
  assert.equal((await r.json()).data[0].target_id, null);
  for (const args of [
    "q=ALIAS",
    "q=catalog&track=secreted&status=ok&structure=unavailable&constructs=available",
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
  r = await fetch(`${f.root}/mouse/targets?q=catalog`);
  assert.equal((await r.json()).data[0].taxon_id, 10090);
  assert.equal((await fetch(`${f.root}/human/downloads`)).status, 404);
  for (const [path, status] of [
    [`${f.root}/human/targets?q=catalog&limit=101`, 400],
    [`${f.root}/human/targets?q=a&q=b`, 400],
    [`${f.target}?q=a`, 400],
    [`${f.root}/human/targets?cursor=${first.next_cursor}&q=changed`, 400],
    [`${f.root}/mouse/targets/${f.target_id}`, 400],
    [`${f.root}/mouse/targets?disease_id=EFO:1`, 400],
    [`${f.root}/human/targets?q=catalog&disease_id=EFO:0001`, 400],
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
  const f = await setup(t, { clientRate: 1, clientBurst: 1 });
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
test("official MCP clients negotiate current and legacy, list five read-only tools and call provider", async (t) => {
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
      assert.equal(tools.tools.length, 5);
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
  assert.equal(result.result.tools.length, 5);
  r = await wire(f.base, "tools/call", {
    name: "search_targets",
    arguments: { release_id: f.release_id, species: "mouse", q: "catalog" },
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
  const global = await setup(t, { globalRate: 1, globalBurst: 1 });
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
      ...testAdmission,
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
          `${base}/v1/releases/${provider.current}/${species}/targets?q=synthetic`,
        )
      ).json();
      const mcp = await client.callTool({
        name: "search_targets",
        arguments: { ...args, q: "synthetic" },
      });
      assert.deepEqual(mcp.structuredContent, http);
      const target = http.data.find((v) => v.target_id);
      assert.ok(target);
      const params = { ...args, target_id: target.target_id };
      for (const name of ["get_target", "list_constructs", "get_sequence"]) {
        const result = await client.callTool({
          name,
          arguments: params,
        });
        assert.notEqual(result.isError, true, name);
        assert.deepEqual(
          result.structuredContent,
          await provider.call(name, params),
        );
      }
      const downloads = provider.resolve(provider.current).manifest.species[
        species
      ].downloads;
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
    const releaseMetadata = await provider.call("get_release", {});
    assert.deepEqual(releaseMetadata.provenance.projection_methods, {
      disease: {
        selection: "portal_top_100",
        species: "human",
        limit: 100,
        order: "published_source_order",
      },
    });
    const firstDisease = await provider.call("get_target", {
      release_id: provider.current,
      species: "human",
      target_id: "9606:uniprot:P00001",
      section: "disease",
      limit: 10,
    });
    assert.equal(firstDisease.total, 100);
    assert.equal(firstDisease.coverage.source_count, 105);
    assert.equal(firstDisease.coverage.projected_count, 100);
    assert.equal(firstDisease.coverage.limit, 100);
    assert.equal(firstDisease.coverage.source_reference.state, "live");
    assert.ok(
      /^[a-f0-9]{64}$/.test(firstDisease.coverage.source_reference.sha256),
    );
    const allDisease = [];
    let diseasePage = firstDisease;
    for (;;) {
      allDisease.push(...diseasePage.data);
      if (!diseasePage.next_cursor) break;
      diseasePage = await provider.call("get_target", {
        release_id: provider.current,
        species: "human",
        target_id: "9606:uniprot:P00001",
        section: "disease",
        limit: 10,
        cursor: diseasePage.next_cursor,
      });
      assert.deepEqual(diseasePage.coverage, firstDisease.coverage);
    }
    assert.deepEqual(
      allDisease.map((row) => row.disease_id),
      Array.from(
        { length: 100 },
        (_, index) => `EFO_SYNTHETIC_${String(index).padStart(4, "0")}`,
      ),
    );
    const diseaseHttp = await (
      await fetch(
        `${base}/v1/releases/${provider.current}/human/targets/9606:uniprot:P00001/sections/disease?limit=10`,
      )
    ).json();
    const diseaseMcp = await client.callTool({
      name: "get_target",
      arguments: {
        release_id: provider.current,
        species: "human",
        target_id: "9606:uniprot:P00001",
        section: "disease",
        limit: 10,
      },
    });
    assert.deepEqual(diseaseMcp.structuredContent, diseaseHttp);
    assert.deepEqual(diseaseHttp.coverage, firstDisease.coverage);
    assert.equal(
      (
        await fetch(
          `${base}/v1/releases/${provider.current}/human/targets?q=SYNTHETIC&disease_id=EFO_SYNTHETIC_0000`,
        )
      ).status,
      400,
    );
    const mouseSummary = (
      await provider.call("get_target", {
        release_id: provider.current,
        species: "mouse",
        target_id: "10090:uniprot:P00002",
      })
    ).data;
    assert.equal(mouseSummary.disease_coverage.state, "not_applicable");
    for (const [species, q, target_id] of [
      ["mouse", "P00002", "10090:uniprot:P00002"],
      ["mouse", "NP_000001.1", "10090:uniprot:P00002"],
      ["human", "P00001-1", "9606:uniprot:P00001"],
    ]) {
      const response = await fetch(
        `${base}/v1/releases/${provider.current}/${species}/targets?q=${encodeURIComponent(q)}`,
      );
      assert.equal(response.status, 200);
      const http = await response.json();
      assert.equal(http.total, 1, q);
      assert.equal(http.data[0].target_id, target_id);
      const result = await client.callTool({
        name: "search_targets",
        arguments: { release_id: provider.current, species, q },
      });
      assert.notEqual(result.isError, true);
      assert.deepEqual(result.structuredContent, http);
      const sequence = (
        await provider.call("get_sequence", {
          release_id: provider.current,
          species,
          target_id,
        })
      ).data;
      assert.equal(sequence.sequence_sha256, http.data[0].sequence_sha256);
      assert.equal(
        sequence.sequence_accession,
        http.data[0].sequence_accession,
      );
      if (species === "mouse") {
        assert.equal(http.data[0].identifier_accession, "P00002");
        assert.equal(sequence.sequence_accession, "NP_000001.1");
      } else assert.equal(http.data[0].canonical_isoform_id, "P00001-1");
      assert.ok(
        provider
          .fallback(provider.current, species)
          .url.includes(provider.current),
      );
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
      const position = manifest.species.human.downloads.findIndex(
        (download) => download.name === "records.jsonl.gz",
      );
      manifest.species.human.downloads[position] =
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
  const first = await fixture("first", {
      legacyRootSchema: true,
      legacyIndex: true,
    }),
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
  const server = makeServer(provider, testAdmission);
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
  assert.equal(current.retirement_history.known_retired_count, 1);
  assert.equal(Object.hasOwn(current, "retired_releases"), false);
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
  for (const method of ["get", "head"]) {
    const parameters =
      description.paths["/v1/releases/{release_id}/{species}/targets"][method]
        .parameters;
    assert.equal(
      parameters.find((parameter) => parameter.name === "q").required,
      true,
    );
    assert.equal(
      parameters.find((parameter) => parameter.name === "q").schema.minLength,
      2,
    );
    assert.ok(
      parameters
        .filter(
          (parameter) => parameter.in === "query" && parameter.name !== "q",
        )
        .every((parameter) => parameter.required !== true),
    );
  }
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

test("measurement trial separates 96 MiB manifests and aggregate indexes from unchanged shard limits", async (t) => {
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
    "Release manifest remains bounded at 96 MiB",
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
  assert.equal(LIMITS.manifest, 96 * 1024 * 1024);
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
  const server = makeServer(provider, testAdmission);
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
    `/v1/releases/${f.release_id}/human/targets?q=catalog&limit=1`,
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
      { release_id: f.release_id, species: "human", q: "catalog", limit: 1 },
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

test("legal stored sequences with oversized encoded FASTA retain exact download recovery", async (t) => {
  const cases = [
    {
      name: "ASCII reserved",
      warnings: Array.from({ length: 40 }, () => "!".repeat(12000)),
    },
    {
      name: "Unicode",
      warnings: Array.from({ length: 50 }, () => "é".repeat(6000)),
    },
  ];
  if (process.env.OPENANTIGENS_TEST_OVERFLOW_DATA_DIR)
    cases.push({
      name: "actual Python producer",
      directory: process.env.OPENANTIGENS_TEST_OVERFLOW_DATA_DIR,
    });
  for (const item of cases) {
    const stored = item.directory
      ? { directory: item.directory }
      : await fixture(`stored-overflow-${item.name}`, {
          sequenceWarnings: item.warnings,
        });
    if (!item.directory)
      t.after(() => rm(stored.directory, { recursive: true, force: true }));
    const provider = await new Provider(stored.directory).initialize();
    const species = "human",
      release_id = provider.current;
    const target_id = (
      await provider.call("search_targets", {
        release_id,
        species,
        q: item.directory ? "synthetic" : "catalog",
      })
    ).data.find((row) => row.target_id).target_id;
    const target = await provider.json(
      provider.resolve(release_id),
      provider.resolve(release_id).manifest.species[species].target_paths[
        target_id
      ],
      LIMITS.target,
    );
    let envelope;
    for (const construct_id of [
      null,
      ...Object.keys(target.construct_sequence_paths),
    ]) {
      const candidate = await provider.call("get_sequence", {
        release_id,
        species,
        target_id,
        ...(construct_id ? { construct_id } : {}),
      });
      const encoded = fasta(candidate.data, {
        release_id,
        provenance_url: candidate.manifest.url,
        records_url: provider.fallback(release_id, species).url,
      });
      if (Buffer.byteLength(encoded) > LIMITS.sequence) {
        envelope = candidate;
        break;
      }
    }
    assert.ok(envelope, `${item.name} exercises encoded FASTA overflow`);
    assert.ok(
      Buffer.byteLength(JSON.stringify(envelope)) < LIMITS.sequence,
      "The actual JSON sequence result remains legal",
    );
    assert.equal(
      provider.toolOutput("get_sequence").safeParse(envelope).success,
      true,
    );
    const selector = {
      target_id: envelope.data.target_id,
      construct_id: envelope.data.construct_id,
    };
    const server = makeServer(provider, {
      ...testAdmission,
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
    const suffix = `targets/${target_id}/${selector.construct_id ? `constructs/${selector.construct_id}/` : ""}sequence`;
    const route = `${base}/v1/releases/${release_id}/${species}/${suffix}`;
    const normal = await fetch(route);
    assert.equal(normal.status, 200);
    assert.deepEqual(await normal.json(), envelope);
    const response = await fetch(route + "?format=fasta", {
      headers: { accept: "text/x-fasta" },
    });
    assert.equal(response.status, 413, item.name);
    assert.match(response.headers.get("content-type"), /^application\/json/);
    assert.equal(response.headers.get("cache-control"), "no-store");
    assert.equal(response.headers.has("content-disposition"), false);
    assert.equal(response.headers.has("link"), false);
    const error = (await response.json()).error;
    assert.equal(error.code, "result_too_large");
    assert.ok(error.request_id);
    assert.deepEqual(error.download, provider.fallback(release_id, species));
    assert.deepEqual(error.record_selector, selector);
    assert.match(
      error.download.url,
      new RegExp(`/data/releases/${release_id}/human/records\\.jsonl\\.gz$`),
    );
    const records = gunzipSync(
      await readFile(join(stored.directory, error.download.path)),
    )
      .toString("utf8")
      .trimEnd()
      .split("\n")
      .map((line) => JSON.parse(line));
    const recovered = records.filter(
      (record) =>
        record.kind === "sequence" &&
        record.payload.data.target_id === selector.target_id &&
        record.payload.data.construct_id === selector.construct_id,
    );
    assert.equal(recovered.length, 1);
    assert.deepEqual(recovered[0].payload.data, envelope.data);
    for (const [url, status] of [
      [route + "?format=xml", 400],
      [route + "?format=fasta&limit=1", 400],
      [route.replace("/human/", "/rat/") + "?format=fasta", 400],
      [route.replace("/human/", "/mouse/") + "?format=fasta", 400],
    ])
      assert.equal((await fetch(url)).status, status);
    const client = new Client(
      { name: "overflow-json-control", version: "1" },
      { versionNegotiation: { mode: "auto" } },
    );
    await client.connect(
      new StreamableHTTPClientTransport(new URL(base + "/mcp")),
    );
    t.after(() => client.close());
    const mcp = await client.callTool({
      name: "get_sequence",
      arguments: {
        release_id,
        species,
        target_id,
        ...(selector.construct_id
          ? { construct_id: selector.construct_id }
          : {}),
      },
    });
    assert.notEqual(mcp.isError, true);
    assert.deepEqual(mcp.structuredContent, envelope);
  }
  assert.equal(LIMITS.sequence, 1024 * 1024);
});

test("legal long retirement history keeps current and pinned HTTP/MCP discovery bounded", async (t) => {
  const f = await fixture("long-retirement-history");
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const rootPath = join(f.directory, "manifest.json");
  const root = JSON.parse(await readFile(rootPath));
  root.retired_releases = Array.from({ length: 300 }, (_, index) => ({
    release_id: createHash("sha256").update(`retired-${index}`).digest("hex"),
    retired_by_release_id: f.release_id,
    archive_url: `https://archive.example/${"x".repeat(740)}?record=${index}`,
  }));
  const bytes = Buffer.from(JSON.stringify(root));
  assert.ok(bytes.length > LIMITS.json && bytes.length < LIMITS.shard);
  await writeFile(rootPath, bytes);
  const provider = await new Provider(f.directory).initialize();
  const server = makeServer(provider, testAdmission);
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  t.after(
    () =>
      new Promise((r) => {
        server.closeAllConnections();
        server.close(r);
      }),
  );
  const base = `http://127.0.0.1:${server.address().port}`;
  for (const path of [
    "/v1/releases",
    "/v1/releases/current",
    `/v1/releases/${f.release_id}`,
  ]) {
    const response = await fetch(base + path);
    assert.equal(response.status, 200);
    const text = await response.text();
    assert.ok(Buffer.byteLength(text) < LIMITS.json);
    const result = JSON.parse(text);
    assert.deepEqual(result.retirement_history, {
      known_retired_count: 300,
      limit: 1000,
      lookup_endpoint: "/v1/releases/{release_id}",
    });
    assert.equal(Object.hasOwn(result, "retired_releases"), false);
  }
  for (const args of [{}, { release_id: f.release_id }]) {
    const wireResult = await protocolBody(
      await wire(base, "tools/call", { name: "get_release", arguments: args }),
    );
    assert.notEqual(wireResult.result.isError, true);
    assert.equal(
      wireResult.result.structuredContent.retirement_history
        .known_retired_count,
      300,
    );
    assert.ok(Buffer.byteLength(JSON.stringify(wireResult)) < LIMITS.json);
  }
  const retired = await fetch(
    `${base}/v1/releases/${root.retired_releases[0].release_id}`,
  );
  assert.equal(retired.status, 410);
  assert.deepEqual(
    (await retired.json()).error.retirement,
    root.retired_releases[0],
  );
});

test("compact indexes return full summaries through bounded interactive search", async (t) => {
  const f = await fixture("interactive-search-window", { extraUnresolved: 25 });
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const provider = await new Provider(f.directory).initialize();
  const index = provider.resolve(f.release_id).indexes.human;
  assert.equal(index.length, 27);
  for (const row of index) {
    assert.equal(Object.hasOwn(row, "protein_name"), false);
    assert.equal(Object.hasOwn(row, "aliases"), false);
    assert.equal(Object.hasOwn(row, "disease_ids"), false);
    assert.equal(Object.hasOwn(row, "design_region"), false);
    assert.ok(row.summary_path);
  }
  let reads = 0;
  const json = provider.json.bind(provider);
  provider.json = async (...args) => {
    if (args[1].includes("/targets/")) reads++;
    return json(...args);
  };
  const args = {
    release_id: f.release_id,
    species: "human",
    q: "catalog",
    limit: 10,
  };
  const first = await provider.call("search_targets", args);
  assert.equal(reads, 10);
  assert.equal(first.total, 27);
  assert.equal(first.data.length, 10);
  assert.equal(first.refine_query, true);
  assert.equal(first.search_window, 20);
  assert.ok(first.next_cursor);
  for (const summary of first.data) {
    assert.ok(Object.hasOwn(summary, "protein_name"));
    assert.ok(Object.hasOwn(summary, "disease_ids"));
    assert.equal(Object.hasOwn(summary, "summary_path"), false);
  }
  const second = await provider.call("search_targets", {
    ...args,
    cursor: first.next_cursor,
  });
  assert.equal(reads, 20);
  assert.equal(second.data.length, 10);
  assert.equal(second.next_cursor, null);
  assert.equal(
    new Set([...first.data, ...second.data].map((row) => row.record_id)).size,
    20,
  );
  const cursor = JSON.parse(Buffer.from(first.next_cursor, "base64url"));
  cursor.position = 20;
  await assert.rejects(
    provider.call("search_targets", {
      ...args,
      cursor: Buffer.from(JSON.stringify(cursor)).toString("base64url"),
    }),
    (error) => error.code === "refine_query",
  );
  reads = 0;
  await assert.rejects(
    provider.call("search_targets", { ...args, disease_id: "EFO:0001" }),
    (error) => error.code === "invalid_arguments",
  );
  assert.equal(
    reads,
    0,
    "Unsupported disease membership filters load no records",
  );
  const server = makeServer(provider, testAdmission);
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  t.after(
    () =>
      new Promise((r) => {
        server.closeAllConnections();
        server.close(r);
      }),
  );
  const base = `http://127.0.0.1:${server.address().port}/v1/releases/${f.release_id}/human`;
  for (const suffix of [
    "targets",
    "targets?q=a",
    "targets?q=%20%20",
    "targets?q=ca*&limit=1",
    "targets?q=catalog&limit=11",
    "targets?track=secreted",
  ]) {
    assert.equal((await fetch(`${base}/${suffix}`)).status, 400, suffix);
  }
  const normal = await (await fetch(`${base}/targets?q=catalog`)).json();
  assert.equal(normal.data.length, 5);
  assert.equal((await fetch(`${base}/downloads`)).status, 404);
  const metadata = await provider.call("get_release", {});
  for (const species of ["human", "mouse"]) {
    const info = metadata.data.species[species];
    assert.equal(Object.hasOwn(info, "downloads"), false);
    assert.deepEqual(
      info.references.map((reference) => reference.name).sort(),
      ["citation.json", "openantigens.bib", "openantigens.ris", "schema.json"],
    );
  }
});
test("correctly hashed private, unknown-download and unreachable typed artifacts fail readiness", async (t) => {
  const checksum = (value) => createHash("sha256").update(value).digest("hex");
  for (const kind of ["private", "unknown_download", "orphan_typed"]) {
    const f = await fixture(`artifact-rejection-${kind}`);
    t.after(() => rm(f.directory, { recursive: true, force: true }));
    const rootPath = join(f.directory, "manifest.json");
    const root = JSON.parse(await readFile(rootPath));
    const manifestPath = join(f.directory, root.releases[0].manifest_path);
    const manifest = JSON.parse(await readFile(manifestPath));
    let path, bytes, file;
    if (kind === "orphan_typed") {
      const original = Object.keys(manifest.files).find(
        (path) =>
          path.includes("/human/sections/") && path.endsWith("topology-0.json"),
      );
      path = original.replace("topology-0.json", "topology-1.json");
      bytes = await readFile(join(f.directory, original));
      file = { ...manifest.files[original] };
    } else {
      path = `releases/${f.release_id}/human/downloads/${kind === "private" ? "private.sqlite" : "extra.json"}`;
      bytes = Buffer.from("PRIVATE_CANARY synthetic raw confidential input\n");
      file = {
        media_type: "application/octet-stream",
        schema: null,
        encoding: null,
      };
    }
    await writeFile(join(f.directory, path), bytes);
    Object.assign(file, { bytes: bytes.length, sha256: checksum(bytes) });
    manifest.files[path] = file;
    if (kind === "unknown_download")
      manifest.species.human.downloads.push({
        name: "extra.json",
        url: `https://openantigens.org/data/releases/${f.release_id}/human/extra.json`,
        path,
        ...file,
        row_count: null,
        license: "fixture",
      });
    const manifestBytes = Buffer.from(JSON.stringify(manifest) + "\n");
    await writeFile(manifestPath, manifestBytes);
    root.releases[0].sha256 = checksum(manifestBytes);
    await writeFile(rootPath, JSON.stringify(root));
    await assert.rejects(
      new Provider(f.directory).initialize(),
      kind === "orphan_typed"
        ? /Unreachable artifact file/
        : /Unapproved artifact role|Download role coverage/,
    );
  }
});

test("sealed packed Provider preserves query results without shipping external bulk bytes", async (t) => {
  const source = await fixture("packed-provider-equivalence");
  t.after(() => rm(source.directory, { recursive: true, force: true }));
  const plain = await new Provider(source.directory).initialize();
  const packed = await packedFixture(source.directory);
  t.after(() => rm(packed.directory, { recursive: true, force: true }));
  await assert.rejects(
    new Provider(packed.directory).initialize(),
    /explicit pins/,
  );
  await assert.rejects(
    new Provider(packed.directory, {
      ...packed.pins,
      expectedRootSha256: "f".repeat(64),
    }).initialize(),
    /root pin mismatch/,
  );
  const transportPath = join(packed.directory, "transport-index.json");
  const transportBytes = await readFile(transportPath);
  await rm(transportPath);
  await assert.rejects(
    new Provider(packed.directory, packed.pins).initialize(),
    (error) => error.code === "ENOENT",
  );
  await writeFile(transportPath, transportBytes);
  const provider = await new Provider(
    packed.directory,
    packed.pins,
  ).initialize();
  assert.equal(provider.storage.packed, true);
  const args = {
    release_id: provider.current,
    species: "human",
    target_id: source.target_id,
  };
  for (const [name, parameters] of [
    ["get_release", {}],
    [
      "search_targets",
      { release_id: provider.current, species: "human", q: "catalog" },
    ],
    ["get_target", args],
    ["get_target", { ...args, section: "topology" }],
    ["list_constructs", args],
    ["get_sequence", args],
    ["get_sequence", { ...args, construct_id: source.construct_id }],
  ])
    assert.deepEqual(
      await provider.call(name, parameters),
      await plain.call(name, parameters),
      name,
    );
  const download = provider.fallback(provider.current, "human");
  await assert.rejects(
    readFile(join(packed.directory, download.path)),
    (error) => error.code === "ENOENT",
  );
  const server = makeServer(provider, testAdmission);
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  t.after(
    () =>
      new Promise((resolve) => {
        server.closeAllConnections();
        server.close(resolve);
      }),
  );
  const base = `http://127.0.0.1:${server.address().port}`;
  assert.equal((await fetch(base + "/readyz")).status, 200);
  const client = new Client(
    { name: "packed-client", version: "1" },
    { versionNegotiation: { mode: "auto" } },
  );
  await client.connect(
    new StreamableHTTPClientTransport(new URL(base + "/mcp")),
  );
  t.after(() => client.close());
  const result = await client.callTool({
    name: "get_sequence",
    arguments: { ...args, construct_id: source.construct_id },
  });
  assert.notEqual(result.isError, true);
  assert.deepEqual(
    result.structuredContent,
    await plain.call("get_sequence", {
      ...args,
      construct_id: source.construct_id,
    }),
  );
});
test("trusted app-code build seal yields independent pins and rejects malformed or linked seals", async (t) => {
  const f = await fixture("build-seal-reader");
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const path = join(f.directory, "service-build.json");
  const data_transport = {
    format: "openantigens-target-gzip-v1",
    root_manifest_sha256: "a".repeat(64),
    transport_index_sha256: "b".repeat(64),
  };
  assert.deepEqual(await readBuildPins(path), {});
  await writeFile(
    path,
    JSON.stringify({ service_build_id: "fixture", data_transport }),
  );
  assert.deepEqual(await readBuildPins(path), {
    expectedRootSha256: data_transport.root_manifest_sha256,
    expectedTransportSha256: data_transport.transport_index_sha256,
  });
  await writeFile(
    path,
    JSON.stringify({
      data_transport: { ...data_transport, format: "unknown" },
    }),
  );
  await assert.rejects(readBuildPins(path), /Invalid service build seal/);
  await rm(path);
  await symlink(join(f.directory, "manifest.json"), path);
  await assert.rejects(readBuildPins(path), /Invalid service build seal/);
});

test("approved external download roles require their exact immutable public URL", async (t) => {
  const f = await fixture("external-download-url");
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const rootPath = join(f.directory, "manifest.json");
  const root = JSON.parse(await readFile(rootPath));
  const path = join(f.directory, root.releases[0].manifest_path);
  const manifest = JSON.parse(await readFile(path));
  const download = manifest.species.human.downloads.find(
    (item) => item.name === "records.jsonl.gz",
  );
  download.url = `https://openantigens.org/data/releases/${f.release_id}/human/private.sqlite`;
  manifest.species.human.records_download = download;
  const bytes = Buffer.from(JSON.stringify(manifest) + "\n");
  await writeFile(path, bytes);
  root.releases[0].sha256 = createHash("sha256").update(bytes).digest("hex");
  await writeFile(rootPath, JSON.stringify(root));
  await assert.rejects(
    new Provider(f.directory).initialize(),
    /Download manifest/,
  );
});

test("publication state configuration accepts only the documented enum", () => {
  for (const value of [undefined, "", "active"])
    assert.equal(publicationState(value), "active");
  assert.equal(publicationState("maintenance"), "maintenance");
  for (const value of [
    "maintenence",
    " maintenance ",
    "ACTIVE",
    null,
    false,
    0,
  ])
    assert.throws(() => publicationState(value), {
      message: "Invalid publication state configuration.",
    });
});
test("a publication-state typo exits before data access or any listener without logging raw configuration", () => {
  const result = spawnSync(
    process.execPath,
    [fileURLToPath(new URL("../src/server.js", import.meta.url))],
    {
      encoding: "utf8",
      timeout: 2000,
      maxBuffer: 16384,
      env: {
        ...process.env,
        OPENANTIGENS_PUBLICATION_STATE: "PRIVATE_BAD_PUBLICATION_STATE",
        OPENANTIGENS_DATA_DIR: "/private/tmp/DO_NOT_READ_PRIVATE_DATA_ROOT",
        OPENANTIGENS_ANALYTICS_PLAUSIBLE: "0",
        PORT: "43129",
      },
    },
  );
  assert.equal(result.status, 1);
  assert.equal(result.signal, null);
  assert.equal(result.stdout, "");
  assert.equal(result.stderr, "public_api startup_failed\n");
  for (const value of [
    "PRIVATE_BAD_PUBLICATION_STATE",
    "DO_NOT_READ_PRIVATE_DATA_ROOT",
    "listening",
    "bundle_validation_failed",
  ])
    assert.equal((result.stdout + result.stderr).includes(value), false, value);
});

test("only the seven external bulk roles accept metadata above 128 MiB and unpacked bodies stream", async (t) => {
  const f = await fixture("external-role-budget");
  t.after(() => rm(f.directory, { recursive: true, force: true }));
  const provider = new Provider(f.directory);
  const external = [
    "index.jsonl",
    "index.tsv",
    "constructs.jsonl",
    "constructs.tsv",
    "constructs.fasta",
    "reference.fasta",
    "records.jsonl.gz",
  ];
  const mutate = (name, size) => {
    const manifest = structuredClone(f.manifest);
    const download = manifest.species.human.downloads.find(
      (item) => item.name === name,
    );
    download.bytes = size;
    manifest.files[download.path].bytes = size;
    if (name === "records.jsonl.gz")
      manifest.species.human.records_download.bytes = size;
    return { manifest, download };
  };
  for (const name of external) {
    const { manifest, download } = mutate(name, 129 * 1024 * 1024);
    assert.equal(
      provider.approvedFiles(manifest).get(download.path).bound,
      1024 * 1024 * 1024,
      name,
    );
    const tooLarge = mutate(name, LIMITS.externalDownload + 1);
    assert.throws(
      () => provider.approvedFiles(tooLarge.manifest),
      (error) =>
        error.message ===
        `Bundle size: role=${name} actual=${LIMITS.externalDownload + 1} limit=${LIMITS.externalDownload}`,
    );
  }
  for (const name of [
    "schema.json",
    "citation.json",
    "openantigens.bib",
    "openantigens.ris",
  ]) {
    const { manifest } = mutate(name, LIMITS.shard + 1);
    assert.throws(() => provider.approvedFiles(manifest), /Bundle size/);
  }
  for (const [schema, limit] of [
    ["target_shard", LIMITS.target],
    ["sequence_shard", LIMITS.shard],
    ["section_topology", LIMITS.shard],
    ["index", LIMITS.index],
  ]) {
    const manifest = structuredClone(f.manifest);
    const path = Object.keys(manifest.files).find(
      (path) => manifest.files[path].schema === schema,
    );
    manifest.files[path].bytes = limit + 1;
    assert.throws(() => provider.approvedFiles(manifest), /Bundle size/);
  }
  const original = provider.readPhysical.bind(provider);
  let streams = 0;
  provider.readPhysical = async (path, max, expected, collect = true) => {
    if (external.some((name) => path.endsWith("/downloads/" + name))) {
      assert.equal(max, LIMITS.externalDownload);
      assert.equal(collect, false);
      streams++;
    }
    return original(path, max, expected, collect);
  };
  await provider.initialize();
  assert.equal(streams, 14);
  assert.equal(LIMITS.json, 256 * 1024);
  assert.equal(LIMITS.sequence, 1024 * 1024);
  assert.equal(LIMITS.manifest, 96 * 1024 * 1024);
  assert.equal(LIMITS.index, 32 * 1024 * 1024);
  assert.equal(LIMITS.indexes, 64 * 1024 * 1024);
});

import test from "node:test";
import assert from "node:assert/strict";
import { request } from "node:http";
import { rm } from "node:fs/promises";
import { Provider } from "../src/provider.js";
import { makeServer } from "../src/server.js";
import { makeAnalytics } from "../src/analytics.js";
import { fixture } from "./fixture.js";
const meta = {
  "io.modelcontextprotocol/protocolVersion": "2026-07-28",
  "io.modelcontextprotocol/clientCapabilities": {},
};
async function setup(t, options = {}) {
  const f = await fixture("analytics-integration");
  const provider = await new Provider(f.directory).initialize();
  const lines = [],
    sent = [];
  const analytics = makeAnalytics({
    clientSecret: "unit-test-secret-at-least-thirty-two-characters",
    write: (line) => {
      lines.push(JSON.parse(line));
      return true;
    },
    plausible: { enabled: options.collector === true, identityVerified: true },
    fetch: async (url, init) => {
      sent.push({ url, init });
      return new Response(null, { status: 202 });
    },
  });
  const server = makeServer(provider, {
    clientRate: 1000,
    clientBurst: 1000,
    globalRate: 10000,
    globalBurst: 10000,
    analytics,
    identityVerified: true,
    ...options,
  });
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  t.after(async () => {
    server.closeAllConnections();
    await new Promise((resolve) => server.close(resolve));
    await analytics.settled();
    await rm(f.directory, { recursive: true, force: true });
  });
  const base = `http://127.0.0.1:${server.address().port}`;
  return {
    ...f,
    provider,
    analytics,
    lines,
    sent,
    base,
    requests: () => lines.filter((line) => line.type === "api_request"),
  };
}
async function tool(base, name, args) {
  const response = await fetch(base + "/mcp", {
    method: "POST",
    headers: {
      "content-type": "application/json",
      accept: "application/json, text/event-stream",
      "mcp-protocol-version": "2026-07-28",
      "mcp-method": "tools/call",
      "mcp-name": name,
      "user-agent": "secret-unit-UA",
    },
    body: JSON.stringify({
      jsonrpc: "2.0",
      id: 1,
      method: "tools/call",
      params: { name, arguments: args, _meta: meta },
    }),
  });
  return { response, body: await response.json() };
}
test("HTTP and MCP completion hooks record sanitized successes, logical errors and rejected requests once", async (t) => {
  const f = await setup(t);
  const target = `${f.base}/v1/releases/${f.release_id}/human/targets/${f.target_id}`;
  let response = await fetch(target);
  assert.equal(response.status, 200);
  const text = await response.text();
  response = await fetch(target, { method: "HEAD" });
  assert.equal(response.status, 200);
  await response.text();
  response = await fetch(
    `${f.base}/v1/releases/${f.release_id}/human/targets?q=secret-probe`,
    { headers: { "user-agent": "secret-unit-UA" } },
  );
  assert.equal(response.status, 200);
  await response.text();
  const result = await tool(f.base, "get_target", {
    release_id: f.release_id,
    species: "human",
    target_id: "9606:uniprot:NOT_FOUND",
  });
  assert.equal(result.response.status, 200);
  assert.equal(result.body.result.isError, true);
  await new Promise((resolve, reject) => {
    const req = request(
      target,
      { headers: { host: "invalid.example" } },
      (response) => {
        assert.equal(response.statusCode, 403);
        response.resume();
        response.on("end", resolve);
      },
    );
    req.on("error", reject);
    req.end();
  });
  const records = f.requests();
  assert.equal(records.length, 5);
  assert.equal(new Set(records.map((record) => record.request_id)).size, 5);
  assert.equal(records[0].bytes, Buffer.byteLength(text));
  assert.equal(records[1].bytes, 0);
  assert.equal(records[3].transport, "mcp");
  assert.equal(records[3].operation, "get_target");
  assert.equal(records[3].error_code, "not_found");
  assert.equal(records[3].status, 200);
  assert.equal(records[4].admitted, false);
  assert.equal(records[4].error_code, "invalid_host");
  assert.equal(f.sent.length, 0);
  const ledger = JSON.stringify(f.lines);
  for (const secret of [
    "secret-probe",
    "secret-unit-UA",
    f.target_id,
    "127.0.0.1",
    f.directory,
  ])
    assert.equal(ledger.includes(secret), false, secret);
});
test("rate-limit and aborted completion records are retained while optional collector excludes rejected traffic", async (t) => {
  const f = await setup(t, { collector: true, clientBurst: 1 });
  let response = await fetch(`${f.base}/v1/releases/current`);
  assert.equal(response.status, 200);
  await response.text();
  response = await fetch(`${f.base}/v1/releases/current`);
  assert.equal(response.status, 429);
  await response.text();
  await f.analytics.settled();
  assert.equal(f.requests().length, 2);
  assert.equal(f.requests()[1].admitted, false);
  assert.equal(f.requests()[1].error_code, "rate_limited");
  assert.equal(f.sent.length, 1);
  assert.ok(
    f.lines.some(
      (line) =>
        line.type === "analytics_delivery" && line.outcome === "policy_skipped",
    ),
  );
  const aborted = await setup(t);
  const original = aborted.provider.call.bind(aborted.provider);
  let release, started;
  const waiting = new Promise((resolve) => {
    release = resolve;
  });
  const began = new Promise((resolve) => {
    started = resolve;
  });
  aborted.provider.call = async (...args) => {
    started();
    await waiting;
    return original(...args);
  };
  const req = request(
    `${aborted.base}/v1/releases/${aborted.release_id}/human/targets/${aborted.target_id}`,
  );
  req.on("error", () => {});
  req.end();
  await began;
  req.destroy();
  await new Promise((resolve) => setTimeout(resolve, 20));
  release();
  await new Promise((resolve) => setTimeout(resolve, 20));
  assert.equal(aborted.requests().length, 1);
  assert.equal(aborted.requests()[0].outcome, "aborted");
  assert.equal(aborted.requests()[0].error_code, "request_aborted");
});
test("opt-in collector records valid MCP business failures and protocol handshakes only in their proper channels", async (t) => {
  const f = await setup(t, { collector: true });
  await tool(f.base, "get_target", {
    release_id: f.release_id,
    species: "human",
    target_id: "9606:uniprot:NOT_FOUND",
  });
  const response = await fetch(f.base + "/mcp", {
    method: "POST",
    headers: {
      "content-type": "application/json",
      accept: "application/json",
      "mcp-protocol-version": "2026-07-28",
      "mcp-method": "server/discover",
    },
    body: JSON.stringify({
      jsonrpc: "2.0",
      id: 2,
      method: "server/discover",
      params: { _meta: meta },
    }),
  });
  assert.equal(response.status, 200);
  await response.text();
  await f.analytics.settled();
  assert.equal(f.requests().length, 2);
  assert.equal(f.requests()[1].operation, "discover");
  assert.equal(f.sent.length, 1);
  const outgoing = JSON.parse(f.sent[0].init.body);
  assert.equal(outgoing.url, "https://openantigens.org/_mcp/human/get_target");
  assert.equal(JSON.stringify(outgoing).includes("NOT_FOUND"), false);
});

test("verified forwarding shares exact client identity between admission and transient analytics", async (t) => {
  const f = await setup(t, {
    collector: true,
    clientBurst: 1,
    trustedProxyPeers: ["127.0.0.1"],
    proxyHopCount: 2,
  });
  const read = async (chain) => {
    const response = await fetch(`${f.base}/v1/releases/current`, {
      headers: {
        "x-forwarded-for": chain,
        "user-agent": "authentic-unit-agent",
      },
    });
    await response.text();
    return response;
  };
  assert.equal((await read("198.51.100.7, 10.0.0.1")).status, 200);
  assert.equal(
    (await read("203.0.113.123, 198.51.100.7, 10.0.0.1")).status,
    429,
    "Spoofed left prefix cannot change the verified client bucket",
  );
  assert.equal(
    (await read("198.51.100.8, 10.0.0.1")).status,
    200,
    "A distinct verified client has its own bucket",
  );
  assert.equal(
    (await read("unknown, 10.0.0.1")).status,
    200,
    "Malformed forwarding uses only the raw peer rate bucket",
  );
  await f.analytics.settled();
  assert.equal(f.sent.length, 2);
  assert.deepEqual(
    f.sent.map((call) => call.init.headers["x-forwarded-for"]),
    ["198.51.100.7", "198.51.100.8"],
  );
  assert.ok(
    f.sent.every(
      (call) => call.init.headers["user-agent"] === "authentic-unit-agent",
    ),
  );
  const records = f.requests();
  assert.equal(records.length, 4);
  assert.equal(records[0].client_id, records[1].client_id);
  assert.notEqual(records[0].client_id, records[2].client_id);
  assert.equal(records[3].client_id, null);
  const ledger = JSON.stringify(f.lines);
  for (const value of [
    "198.51.100.7",
    "198.51.100.8",
    "203.0.113.123",
    "10.0.0.1",
    "127.0.0.1",
    "authentic-unit-agent",
  ])
    assert.equal(ledger.includes(value), false, value);
});

test("SDK-handled current and legacy tool/protocol errors are counted once without diagnostic leakage", async (t) => {
  const f = await setup(t);
  const responses = [];
  const decode = async (response) => {
    const text = await response.text();
    if (response.headers.get("content-type").startsWith("application/json"))
      return JSON.parse(text);
    return JSON.parse(
      text
        .split(/\r?\n/)
        .find((line) => line.startsWith("data:"))
        .slice(5),
    );
  };
  for (const modern of [true, false])
    for (const [name, args] of [
      ["search_targets", { release_id: f.release_id, species: "human" }],
      [
        "search_targets",
        {
          release_id: f.release_id,
          species: "human",
          q: "secret-probe",
          unknown_private: "DO_NOT_LOG",
        },
      ],
      ["not_a_tool_PRIVATE", {}],
    ]) {
      const response = await fetch(f.base + "/mcp", {
        method: "POST",
        headers: {
          "content-type": "application/json",
          accept: "application/json, text/event-stream",
          ...(modern
            ? {
                "mcp-protocol-version": "2026-07-28",
                "mcp-method": "tools/call",
                "mcp-name": name,
              }
            : {}),
        },
        body: JSON.stringify({
          jsonrpc: "2.0",
          id: 1,
          method: "tools/call",
          params: { name, arguments: args, ...(modern ? { _meta: meta } : {}) },
        }),
      });
      assert.equal(response.status, 200);
      const body = await decode(response);
      assert.ok(body.error || body.result.isError === true);
      responses.push(body);
    }
  const records = f.requests();
  assert.equal(
    records.length,
    6,
    "Each handled request has exactly one ledger completion",
  );
  assert.ok(
    records.every(
      (record) =>
        record.status === 200 && record.error_code === "invalid_arguments",
    ),
  );
  assert.equal(
    records.filter((record) => record.error_code !== "none").length,
    6,
  );
  const valid = await tool(f.base, "get_release", {});
  assert.equal(valid.response.status, 200);
  assert.equal(valid.body.result.isError, undefined);
  assert.equal(f.requests().length, 7);
  assert.equal(f.requests()[6].error_code, "none");
  assert.equal(f.requests()[6].operation, "get_release");
  const ledger = JSON.stringify(f.lines);
  for (const secret of [
    "secret-probe",
    "DO_NOT_LOG",
    "not_a_tool_PRIVATE",
    "Input validation",
    "Invalid arguments",
  ])
    assert.equal(ledger.includes(secret), false, secret);
});

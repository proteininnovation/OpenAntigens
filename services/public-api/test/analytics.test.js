import test from "node:test";
import assert from "node:assert/strict";
import { createServer } from "node:http";
import { setTimeout as delay } from "node:timers/promises";
import {
  makeAnalytics,
  ANALYTICS_OPERATIONS,
  ANALYTICS_ERRORS,
  PLAUSIBLE_DAILY_LIMIT,
} from "../src/analytics.js";

const id = "12345678-1234-1234-1234-123456789abc";
const client = {
  ip: "192.0.2.10",
  userAgent: "ScientificAgent/1.0",
  trusted: true,
};
const metadata = {
  request_id: id,
  transport: "api",
  operation: "get_target",
  species: "human",
  method: "GET",
  status: 200,
  duration_ms: 12,
  bytes: 42,
  admitted: true,
};
function fixture(options = {}) {
  const records = [];
  const analytics = makeAnalytics({
    write: (line) => {
      records.push(JSON.parse(line));
      return true;
    },
    now: () => new Date("2026-10-25T12:00:00Z"),
    ...options,
  });
  return { analytics, records };
}
test("default ledger strips secrets and normalizes hostile dimensions without making network calls", async () => {
  let called = false;
  const { analytics, records } = fixture({
    fetch: () => {
      called = true;
      throw new Error("must not send");
    },
  });
  analytics.record(
    {
      ...metadata,
      request_id: "secret-token",
      operation: "/target/private-sequence?q=secret",
      species: "secret",
      method: "secret",
      error_code: "private-body",
      status: Infinity,
      bytes: -1,
      duration_ms: NaN,
      query: "secret",
      body: "sequence",
      url: "secret",
      headers: "secret",
    },
    { ...client, userAgent: "secret-user-agent" },
  );
  await analytics.settled();
  assert.equal(called, false);
  assert.deepEqual(records[0], {
    type: "api_request",
    version: 1,
    timestamp: "2026-10-25T12:00:00.000Z",
    request_id: null,
    transport: "api",
    operation: "unknown",
    species: "none",
    method: "OTHER",
    status: 0,
    error_code: "unknown",
    outcome: "finished",
    duration_ms: 0,
    bytes: 0,
    client_id: null,
    admitted: true,
  });
  assert.equal(records[1].outcome, "disabled");
  assert.doesNotMatch(JSON.stringify(records), /secret|sequence|192\.0\.2/);
});
test("daily anonymous client IDs require a secret and verified identity, and change each UTC day", () => {
  let date = "2026-10-25T12:00:00Z";
  const { analytics, records } = fixture({
    clientSecret: "x".repeat(32),
    now: () => new Date(date),
  });
  const first = analytics.record(metadata, client).client_id;
  assert.match(first, /^[a-f0-9]{32}$/);
  assert.equal(analytics.record(metadata, client).client_id, first);
  assert.notEqual(
    analytics.record(metadata, { ...client, userAgent: "python-httpx/0.28.1" })
      .client_id,
    first,
  );
  assert.notEqual(
    analytics.record(metadata, { ...client, ip: "192.0.2.11" }).client_id,
    first,
  );
  assert.equal(
    analytics.record(metadata, { ...client, trusted: false }).client_id,
    null,
  );
  for (const userAgent of [
    undefined,
    "",
    "ResearchSDK token=private-value",
    "contact@example.org",
  ]) {
    assert.equal(
      analytics.record(metadata, { ...client, userAgent }).client_id,
      null,
    );
  }
  date = "2026-10-26T00:00:00Z";
  assert.notEqual(analytics.record(metadata, client).client_id, first);
  assert.doesNotMatch(
    JSON.stringify(records),
    /192\.0\.2|ScientificAgent|python-httpx|private-value|example\.org/,
  );
});
test("opt-in collector uses authentic transient identity and fixed redacted URLs for both species/transports", async () => {
  const calls = [];
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: async (url, options) => {
      calls.push({ url, options });
      return new Response("{}", { status: 202 });
    },
  });
  analytics.record(
    {
      ...metadata,
      url: "https://attacker.invalid/secret?q=sequence",
      target_id: "private",
    },
    client,
  );
  analytics.record(
    {
      ...metadata,
      transport: "mcp",
      operation: "get_sequence",
      species: "mouse",
      method: "POST",
      status: 200,
      error_code: "not_found",
      admitted: true,
    },
    client,
  );
  await analytics.settled();
  assert.equal(calls.length, 2);
  for (const call of calls) {
    assert.equal(call.url, "https://plausible.io/api/event");
    assert.equal(call.options.headers["user-agent"], client.userAgent);
    assert.equal(call.options.headers["x-forwarded-for"], client.ip);
    const payload = JSON.parse(call.options.body);
    assert.equal(payload.domain, "openantigens.org");
    assert.equal(payload.interactive, false);
    assert.equal(Object.hasOwn(payload, "props"), false);
    assert.doesNotMatch(call.options.body, /private|attacker|sequence\?/);
  }
  assert.deepEqual(
    calls.map((v) => JSON.parse(v.options.body)),
    [
      {
        domain: "openantigens.org",
        name: "API Request",
        url: "https://openantigens.org/_api/human/get_target",
        interactive: false,
      },
      {
        domain: "openantigens.org",
        name: "MCP Request",
        url: "https://openantigens.org/_mcp/mouse/get_sequence",
        interactive: false,
      },
    ],
  );
  assert.equal(records.filter((v) => v.type === "api_request").length, 2);
  assert.equal(
    records.filter(
      (v) => v.type === "analytics_delivery" && v.outcome === "accepted",
    ).length,
    2,
  );
  assert.doesNotMatch(JSON.stringify(records), /192\.0\.2|ScientificAgent/);
});
test("trusted identity and explicit enablement are mandatory", async () => {
  let calls = 0;
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: async () => {
      calls++;
      return new Response(null, { status: 202 });
    },
  });
  for (const identity of [
    { ...client, trusted: false },
    { ...client, ip: "spoofed,192.0.2.10" },
    { ...client, userAgent: "" },
    { ...client, userAgent: "x\r\ninjected" },
  ])
    analytics.record(metadata, identity);
  const other = fixture({
    plausible: { enabled: true },
    fetch: async () => {
      calls++;
    },
  });
  other.analytics.record(metadata, client);
  await analytics.settled();
  assert.equal(calls, 0);
  assert.equal(
    records.filter((v) => v.outcome === "identity_unavailable").length,
    4,
  );
  assert.equal(other.records[1].outcome, "disabled");
});
test("collector refuses credential-bearing, contact-email and control-character user agents without leaking them", async () => {
  let calls = 0;
  let sentUserAgent;
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: async (_url, options) => {
      calls++;
      sentUserAgent = options.headers["user-agent"];
      return new Response(null, { status: 202 });
    },
  });
  const unsafe = [
    "ResearchSDK/1.0 Bearer private-value",
    "ResearchSDK/1.0 token=private-value",
    "ResearchSDK/1.0 access_token: private-value",
    "ResearchSDK/1.0 api_key=private-value",
    "ResearchSDK/1.0 API-Key private-value",
    "ResearchSDK/1.0 password=private-value",
    "ResearchSDK/1.0 Authorization: Basic cHJpdmF0ZTprZXk=",
    "ResearchSDK/1.0 sk-proj-private-value",
    "ResearchSDK/1.0 sk-privatevalue1234",
    "ResearchSDK/1.0 (contact: private@example.org)",
    "ResearchSDK/1.0\tprivate-value",
    "ResearchSDK/1.0\x7fprivate-value",
    "ResearchSDK/1.0\x01private-value",
    "ResearchSDK/1.0\x85private-value",
    "x".repeat(1025),
    "é".repeat(513),
  ];
  for (const userAgent of unsafe)
    analytics.record(metadata, { ...client, userAgent });
  await analytics.settled();
  assert.equal(calls, 0);
  assert.equal(analytics.snapshot().plausible_budget_used, 0);
  assert.equal(
    records.filter((r) => r.type === "api_request").length,
    unsafe.length,
  );
  assert.equal(
    records.filter(
      (r) =>
        r.type === "analytics_delivery" && r.outcome === "identity_unavailable",
    ).length,
    unsafe.length,
  );
  assert.doesNotMatch(
    JSON.stringify(records),
    /private|example\.org|Bearer|api_key|Authorization|sk-proj/,
  );
  analytics.record(metadata, { ...client, userAgent: "python-httpx/0.28.1" });
  await analytics.settled();
  assert.equal(calls, 1);
  assert.equal(sentUserAgent, "python-httpx/0.28.1");
});
test("rejected, aborted, health and handshake traffic never spends Plausible quota", async () => {
  let calls = 0;
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: async () => {
      calls++;
      return new Response(null, { status: 202 });
    },
  });
  const excluded = [
    [
      { admitted: false, status: 429, error_code: "rate_limited" },
      "not_admitted",
    ],
    [
      { admitted: false, status: 403, error_code: "invalid_host" },
      "not_admitted",
    ],
    [
      { outcome: "aborted", status: 0, error_code: "request_aborted" },
      "not_finished",
    ],
    ...[
      "health",
      "ready",
      "openapi",
      "unknown",
      "initialize",
      "initialized",
      "ping",
      "list_tools",
      "call_tool",
      "get_downloads",
      "discover",
    ].map((operation) => [{ operation }, "operation_not_exported"]),
    [{ method: "OPTIONS" }, "operation_not_exported"],
    ...[
      "invalid_arguments",
      "invalid_cursor",
      "species_mismatch",
      "request_too_large",
      "method_not_allowed",
      "unsupported_media_type",
    ].map((error_code) => [{ status: 400, error_code }, "request_failed"]),
    [{ status: 429, error_code: "unavailable" }, "request_failed"],
  ];
  for (const [fields] of excluded)
    analytics.record({ ...metadata, ...fields }, client);
  await analytics.settled();
  assert.equal(calls, 0);
  assert.equal(
    records.filter((r) => r.type === "api_request").length,
    excluded.length,
  );
  const deliveries = records.filter((r) => r.type === "analytics_delivery");
  assert.deepEqual(
    deliveries.map((r) => [r.outcome, r.reason]),
    excluded.map(([, reason]) => ["policy_skipped", reason]),
  );
});
test("admitted valid lookups and known business failures remain eligible", async () => {
  let calls = 0;
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: async () => {
      calls++;
      return new Response(null, { status: 202 });
    },
  });
  const valid = [
    { operation: "search_targets" },
    { operation: "get_release" },
    { operation: "get_section", method: "HEAD", status: 304 },
    { error_code: "not_found", status: 404 },
    { error_code: "refine_query", status: 400 },
    { error_code: "unsupported_capability", status: 400 },
    { error_code: "retired_release", status: 410 },
    { error_code: "publication_in_progress", status: 503 },
    { error_code: "unavailable", status: 503 },
    { error_code: "result_too_large", status: 413 },
    {
      transport: "mcp",
      method: "POST",
      operation: "list_constructs",
      error_code: "result_too_large",
      status: 200,
    },
  ];
  for (const fields of valid) {
    analytics.record({ ...metadata, ...fields }, client);
    await analytics.settled();
  }
  assert.equal(calls, valid.length);
  assert.ok(
    records
      .filter((r) => r.type === "analytics_delivery")
      .every((r) => r.outcome === "accepted"),
  );
});
test("provider filtering, rejection, failures and timeouts are visible without retries", async () => {
  for (const [fetch, expected] of [
    [
      async () =>
        new Response(null, {
          status: 202,
          headers: { "x-plausible-dropped": "1" },
        }),
      "provider_dropped",
    ],
    [async () => new Response(null, { status: 500 }), "failed"],
    [
      async () => {
        throw new Error("private failure body");
      },
      "failed",
    ],
    [async () => new Promise(() => {}), "timeout"],
  ]) {
    let calls = 0;
    const { analytics, records } = fixture({
      plausible: { enabled: true, identityVerified: true },
      timeoutMs: 5,
      fetch: (...args) => {
        calls++;
        return fetch(...args);
      },
    });
    analytics.record(metadata, client);
    await Promise.all([analytics.settled(), delay(10)]);
    assert.equal(calls, 1);
    assert.equal(records[1].outcome, expected);
    assert.doesNotMatch(JSON.stringify(records), /private failure/);
  }
});
test("actual fetch refuses collector redirects without forwarding identity to another origin", async (t) => {
  let destinationRequests = 0;
  let sourceRequests = 0;
  const destination = createServer((_req, res) => {
    destinationRequests++;
    res.writeHead(202).end();
  });
  await new Promise((resolve) => destination.listen(0, "127.0.0.1", resolve));
  const source = createServer((req, res) => {
    sourceRequests++;
    assert.equal(req.headers["user-agent"], client.userAgent);
    assert.equal(req.headers["x-forwarded-for"], client.ip);
    res
      .writeHead(307, {
        location: `http://127.0.0.1:${destination.address().port}/other-collector`,
      })
      .end();
  });
  await new Promise((resolve) => source.listen(0, "127.0.0.1", resolve));
  t.after(async () => {
    for (const server of [source, destination]) {
      server.closeAllConnections();
      await new Promise((resolve) => server.close(resolve));
    }
  });
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: (_url, options) => {
      assert.equal(options.redirect, "error");
      return globalThis.fetch(
        `http://127.0.0.1:${source.address().port}/collector`,
        options,
      );
    },
  });
  analytics.record(metadata, client);
  await analytics.settled();
  assert.equal(sourceRequests, 1);
  assert.equal(destinationRequests, 0);
  assert.equal(records.filter((r) => r.type === "api_request").length, 1);
  assert.equal(
    records.find((r) => r.type === "analytics_delivery").outcome,
    "failed",
  );
  assert.doesNotMatch(
    JSON.stringify(records),
    /192\.0\.2|ScientificAgent|other-collector/,
  );
});
test("slow collectors have a fixed in-flight bound and explicit skipped deliveries", async () => {
  let release;
  const gate = new Promise((resolve) => {
    release = resolve;
  });
  const { analytics, records } = fixture({
    plausible: { enabled: true, identityVerified: true },
    fetch: async () => {
      await gate;
      return new Response(null, { status: 202 });
    },
  });
  for (let i = 0; i < 20; i++)
    analytics.record(
      {
        ...metadata,
        outcome: "finished",
      },
      client,
    );
  assert.equal(analytics.snapshot().in_flight, 4);
  assert.equal(records.filter((v) => v.type === "api_request").length, 20);
  assert.equal(
    records.filter((v) => v.outcome === "capacity_skipped").length,
    16,
  );
  release();
  await analytics.settled();
  assert.equal(analytics.snapshot().in_flight, 0);
});
test("collector attempts have an atomic UTC daily budget that clock rollback cannot reopen", async () => {
  let date = "2026-10-25T23:59:59Z";
  let calls = 0;
  const { analytics, records } = fixture({
    now: () => new Date(date),
    plausible: { enabled: true, identityVerified: true },
    fetch: async () => {
      calls++;
      return new Response(null, { status: 202 });
    },
  });
  assert.equal(PLAUSIBLE_DAILY_LIMIT, 100);
  for (let i = 0; i < PLAUSIBLE_DAILY_LIMIT - 2; i++) {
    analytics.record(metadata, client);
    await analytics.settled();
  }
  for (let i = 0; i < 20; i++) analytics.record(metadata, client);
  await analytics.settled();
  assert.equal(calls, 100);
  assert.equal(analytics.snapshot().plausible_budget_used, 100);
  assert.equal(
    records.filter(
      (r) => r.outcome === "capacity_skipped" && r.reason === "daily_budget",
    ).length,
    18,
  );
  assert.equal(records.filter((r) => r.type === "api_request").length, 118);
  date = "2026-10-26T00:00:00Z";
  for (let i = 0; i < 100; i++) {
    analytics.record(metadata, client);
    await analytics.settled();
  }
  analytics.record(metadata, client);
  await analytics.settled();
  assert.equal(calls, 200);
  assert.equal(analytics.snapshot().plausible_budget_day, "2026-10-26");
  date = "2026-10-25T23:59:59Z";
  analytics.record(metadata, client);
  await analytics.settled();
  assert.equal(calls, 200);
  assert.equal(records.at(-1).reason, "daily_budget");
});
test("stdout backpressure retains bounded aggregates and exact all-request totals", () => {
  const records = [];
  let drain;
  let blocked = true;
  const analytics = makeAnalytics({
    write: (line) => {
      records.push(JSON.parse(line));
      return !blocked;
    },
    onDrain: (callback) => {
      drain = callback;
    },
  });
  for (let i = 0; i < 2000; i++)
    analytics.record({
      ...metadata,
      operation: ANALYTICS_OPERATIONS[i % ANALYTICS_OPERATIONS.length],
      error_code:
        ANALYTICS_ERRORS[Math.floor(i / 16) % ANALYTICS_ERRORS.length],
      status: 400 + (i % 100),
      admitted: i % 2 === 0,
    });
  assert.ok(analytics.snapshot().buffered_aggregate_groups <= 537);
  blocked = false;
  drain();
  assert.equal(analytics.snapshot().buffered_requests, 0);
  const requests = records.filter(
    (v) => v.type === "api_request" || v.type === "api_request_aggregate",
  );
  assert.equal(
    requests.reduce((n, v) => n + (v.count ?? 1), 0),
    2000,
  );
  assert.equal(
    requests.reduce((n, v) => n + (v.bytes_total ?? v.bytes), 0),
    2000 * 42,
  );
  assert.ok(requests.some((v) => v.dimensions_collapsed === true));
  assert.ok(
    requests
      .filter((v) => v.type === "api_request_aggregate")
      .every((v) => v.client_ids_available === false),
  );
  assert.equal(
    records
      .filter((v) => v.type.startsWith("analytics_delivery"))
      .reduce((n, v) => n + (v.count ?? 1), 0),
    2000,
  );
});
test("failed log writes retain counts for recovery and expose the failure", () => {
  let failing = true;
  const records = [];
  const { analytics } = fixture({
    write: (line) => {
      if (failing) throw new Error("stdout unavailable");
      records.push(JSON.parse(line));
      return true;
    },
  });
  analytics.record(metadata, client);
  assert.equal(analytics.snapshot().log_write_errors, 1);
  assert.equal(analytics.snapshot().buffered_requests, 1);
  failing = false;
  analytics.flush();
  assert.equal(
    records.find((v) => v.type === "api_request_aggregate").count,
    1,
  );
  assert.equal(analytics.snapshot().buffered_requests, 0);
});

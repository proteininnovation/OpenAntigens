import { createHmac } from "node:crypto";
import { isIP } from "node:net";

export const ANALYTICS_OPERATIONS = Object.freeze([
  "unknown",
  "get_release",
  "search_targets",
  "get_target",
  "get_section",
  "list_constructs",
  "get_sequence",
  "get_downloads",
  "health",
  "ready",
  "openapi",
  "initialize",
  "initialized",
  "ping",
  "list_tools",
  "call_tool",
  "discover",
]);
export const PLAUSIBLE_DAILY_LIMIT = 100;
export const ANALYTICS_ERRORS = Object.freeze([
  "none",
  "unknown",
  "invalid_arguments",
  "invalid_cursor",
  "body_not_allowed",
  "species_mismatch",
  "refine_query",
  "handler_timeout",
  "invalid_host",
  "invalid_origin",
  "method_not_allowed",
  "not_acceptable",
  "not_found",
  "rate_limited",
  "request_too_large",
  "result_too_large",
  "unavailable",
  "unsupported_media_type",
  "publication_in_progress",
  "retired_release",
  "unsupported_capability",
  "request_aborted",
]);
const pick = (value, allowed, fallback) =>
  allowed.includes(value) ? value : fallback;
const integer = (value) =>
  Number.isSafeInteger(value) && value >= 0 ? value : 0;
const uuid = (value) =>
  typeof value === "string" &&
  /^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$/i.test(value)
    ? value
    : null;
const endpoint = "https://plausible.io/api/event";
const exportedOperations = [
  "get_release",
  "search_targets",
  "get_target",
  "get_section",
  "list_constructs",
  "get_sequence",
];
const exportedErrors = [
  "not_found",
  "refine_query",
  "unsupported_capability",
  "retired_release",
  "publication_in_progress",
  "unavailable",
  "result_too_large",
];
function safeUserAgent(value) {
  return (
    typeof value === "string" &&
    value.length > 0 &&
    Buffer.byteLength(value) <= 1024 &&
    !/[\x00-\x1f\x7f-\x9f]/.test(value) &&
    !/(?:\b(?:bearer|authorization|password|api[_-]?key)\b|token\s*[:=]|sk-proj-|\bsk-[a-z0-9_-]{12,}|[a-z0-9._%+-]+@[a-z0-9.-]+\.[a-z]{2,})/i.test(
      value,
    )
  );
}

function exportExclusion(record) {
  if (!record.admitted) return "not_admitted";
  if (record.outcome !== "finished") return "not_finished";
  if (
    !exportedOperations.includes(record.operation) ||
    !(
      (record.transport === "api" && ["GET", "HEAD"].includes(record.method)) ||
      (record.transport === "mcp" && record.method === "POST")
    )
  )
    return "operation_not_exported";
  if (
    record.status === 429 ||
    record.status === 0 ||
    !(
      (record.status >= 200 &&
        record.status < 400 &&
        record.error_code === "none") ||
      exportedErrors.includes(record.error_code)
    )
  )
    return "request_failed";
  return null;
}

// Client identity comes from the server's verified resolver, shared with admission.
export function makeAnalytics(options = {}) {
  const now = options.now ?? (() => new Date());
  const write = options.write ?? ((line) => process.stdout.write(line));
  const onDrain =
    options.onDrain ?? ((callback) => process.stdout.once("drain", callback));
  const fetchEvent = options.fetch ?? globalThis.fetch;
  const secret =
    typeof options.clientSecret === "string" &&
    options.clientSecret.length >= 32
      ? options.clientSecret
      : null;
  const plausible =
    options.plausible?.enabled === true &&
    options.plausible?.identityVerified === true;
  const pending = new Set();
  const aggregates = new Map();
  const totals = { requests: 0, deliveries: 0, log_write_errors: 0 };
  let budgetDay = "";
  let budgetUsed = 0;
  let blocked = false;
  let retry;
  let waitingDrain = false;

  function schedule() {
    if (retry) return;
    retry = setTimeout(() => {
      retry = undefined;
      flush();
    }, 1000);
    retry.unref();
  }
  function output(record) {
    try {
      if (write(`${JSON.stringify(record)}\n`) === false) {
        blocked = true;
        if (!waitingDrain) {
          waitingDrain = true;
          onDrain(() => {
            waitingDrain = false;
            blocked = false;
            flush();
          });
        }
      }
      return true;
    } catch {
      totals.log_write_errors++;
      blocked = true;
      schedule();
      return false;
    }
  }
  function retain(record) {
    const request = record.type === "api_request";
    let dimensions = request
      ? {
          transport: record.transport,
          operation: record.operation,
          species: record.species,
          method: record.method,
          status: record.status,
          error_code: record.error_code,
          outcome: record.outcome,
          admitted: record.admitted,
          dimensions_collapsed: false,
        }
      : {
          transport: record.transport,
          outcome: record.outcome,
          ...(record.reason ? { reason: record.reason } : {}),
        };
    let key = JSON.stringify([record.type, dimensions]);
    if (
      request &&
      !aggregates.has(key) &&
      [...aggregates.values()].filter((v) => v.type === "api_request_aggregate")
        .length >= 512
    ) {
      dimensions = {
        transport: "unknown",
        operation: "unknown",
        species: "none",
        method: "OTHER",
        status: 0,
        error_code: "unknown",
        outcome: "finished",
        admitted: false,
        dimensions_collapsed: true,
      };
      key = "request_overflow";
    }
    let aggregate = aggregates.get(key);
    if (!aggregate) {
      aggregate = {
        type: request
          ? "api_request_aggregate"
          : "analytics_delivery_aggregate",
        version: 1,
        timestamp: record.timestamp,
        first_timestamp: record.timestamp,
        last_timestamp: record.timestamp,
        ...dimensions,
        count: 0,
        ...(request
          ? {
              duration_ms_total: 0,
              bytes_total: 0,
              client_ids_available: false,
            }
          : {}),
      };
      aggregates.set(key, aggregate);
    }
    aggregate.count++;
    aggregate.last_timestamp = record.timestamp;
    if (request) {
      aggregate.duration_ms_total += record.duration_ms;
      aggregate.bytes_total += record.bytes;
    }
    schedule();
  }
  function log(record) {
    if (blocked || !output(record)) retain(record);
  }
  function flush() {
    if (waitingDrain) return;
    blocked = false;
    for (const [key, aggregate] of aggregates) {
      aggregate.timestamp = now().toISOString();
      if (!output(aggregate)) break;
      aggregates.delete(key);
      if (blocked) break;
    }
    if (aggregates.size) schedule();
  }
  function delivery(record, outcome, reason) {
    totals.deliveries++;
    log({
      type: "analytics_delivery",
      version: 1,
      timestamp: now().toISOString(),
      request_id: record.request_id,
      transport: record.transport,
      outcome,
      ...(reason ? { reason } : {}),
    });
  }
  async function send(record, client) {
    const controller = new AbortController();
    let timer;
    try {
      const response = await Promise.race([
        Promise.resolve().then(async () => {
          const response = await fetchEvent(endpoint, {
            method: "POST",
            redirect: "error",
            signal: controller.signal,
            headers: {
              "content-type": "application/json",
              "user-agent": client.userAgent,
              "x-forwarded-for": client.ip,
            },
            body: JSON.stringify({
              domain: "openantigens.org",
              name: record.transport === "mcp" ? "MCP Request" : "API Request",
              url: `https://openantigens.org/_${record.transport}/${record.species}/${record.operation}`,
              interactive: false,
            }),
          });
          await response.body?.cancel();
          return response;
        }),
        new Promise((_, reject) => {
          timer = setTimeout(() => {
            controller.abort();
            reject(new Error("timeout"));
          }, options.timeoutMs ?? 1500);
          timer.unref();
        }),
      ]);
      // Plausible returns 202 even when its bot filter rejects an event.
      delivery(
        record,
        response.headers.get("x-plausible-dropped") === "1"
          ? "provider_dropped"
          : response.ok
            ? "accepted"
            : "failed",
      );
    } catch {
      delivery(record, controller.signal.aborted ? "timeout" : "failed");
    } finally {
      clearTimeout(timer);
    }
  }
  function record(metadata = {}, client = {}) {
    const timestamp = now().toISOString();
    const ip =
      typeof client.ip === "string" && isIP(client.ip) ? client.ip : null;
    const userAgent = safeUserAgent(client.userAgent) ? client.userAgent : null;
    const value = {
      type: "api_request",
      version: 1,
      timestamp,
      request_id: uuid(metadata.request_id),
      transport: pick(metadata.transport, ["api", "mcp", "unknown"], "unknown"),
      operation: pick(metadata.operation, ANALYTICS_OPERATIONS, "unknown"),
      species: pick(metadata.species, ["human", "mouse", "none"], "none"),
      method: pick(
        metadata.method,
        ["GET", "HEAD", "POST", "OPTIONS"],
        "OTHER",
      ),
      status:
        Number.isInteger(metadata.status) &&
        metadata.status >= 100 &&
        metadata.status <= 599
          ? metadata.status
          : 0,
      error_code: pick(
        metadata.error_code ?? "none",
        ANALYTICS_ERRORS,
        "unknown",
      ),
      outcome: pick(metadata.outcome, ["finished", "aborted"], "finished"),
      duration_ms: integer(metadata.duration_ms),
      bytes: integer(metadata.bytes),
      client_id:
        secret && ip && userAgent && client.trusted === true
          ? createHmac("sha256", secret)
              .update(`${timestamp.slice(0, 10)}\0${ip}\0${userAgent}`)
              .digest("hex")
              .slice(0, 32)
          : null,
      admitted: metadata.admitted === true,
    };
    totals.requests++;
    log(value);
    const exclusion = exportExclusion(value);
    if (!plausible) delivery(value, "disabled");
    else if (exclusion) delivery(value, "policy_skipped", exclusion);
    else if (
      !ip ||
      client.trusted !== true ||
      !userAgent ||
      value.transport === "unknown"
    )
      delivery(value, "identity_unavailable");
    else if (pending.size >= 4) delivery(value, "capacity_skipped");
    else {
      const day = timestamp.slice(0, 10);
      if (day > budgetDay) {
        budgetDay = day;
        budgetUsed = 0;
      }
      if (day < budgetDay || budgetUsed >= PLAUSIBLE_DAILY_LIMIT)
        delivery(value, "capacity_skipped", "daily_budget");
      else {
        // This process-local budget resets on restart; deployment must verify account quota.
        budgetUsed++;
        const task = send(value, { ip, userAgent: client.userAgent });
        pending.add(task);
        void task.finally(() => pending.delete(task));
      }
    }
    return value;
  }
  return {
    record,
    flush,
    settled: async () => {
      await Promise.allSettled([...pending]);
      flush();
    },
    snapshot: () => ({
      ...totals,
      in_flight: pending.size,
      plausible_budget_day: budgetDay || null,
      plausible_budget_used: budgetUsed,
      plausible_daily_limit: PLAUSIBLE_DAILY_LIMIT,
      buffered_aggregate_groups: aggregates.size,
      buffered_requests: [...aggregates.values()]
        .filter((v) => v.type === "api_request_aggregate")
        .reduce((n, v) => n + v.count, 0),
      logging_blocked: blocked,
    }),
  };
}

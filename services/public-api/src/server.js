import { createServer } from "node:http";
import { createHash, randomUUID } from "node:crypto";
import { pathToFileURL } from "node:url";
import { resolve } from "node:path";
import { Provider, fasta } from "./provider.js";
import { makeMcpHandler } from "./mcp.js";
import { ApiError, bad, LIMITS } from "./contracts.js";
import { openapiFor } from "./openapi.js";
const hash = (b) => `"${createHash("sha256").update(b).digest("hex")}"`;
const defaults = {
  hosts: ["api.openantigens.org", "localhost", "127.0.0.1", "[::1]"],
  origins: ["https://openantigens.org", "https://www.openantigens.org"],
  clientRate: 60,
  globalRate: 300,
  concurrency: 8,
  handlerMs: 5000,
};
const allowHeaders = [
  "content-type",
  "accept",
  "mcp-protocol-version",
  "mcp-method",
  "mcp-name",
];
function argumentsFrom(url, allowed) {
  const args = {};
  for (const [key, value] of url.searchParams) {
    if (!allowed.includes(key) || Object.hasOwn(args, key)) bad();
    args[key] = value;
  }
  if (args.limit !== undefined) {
    if (!/^[1-9]\d{0,2}$/.test(args.limit)) bad();
    args.limit = Number(args.limit);
  }
  return args;
}
function accepts(req, type) {
  return (
    !req.headers.accept ||
    req.headers.accept.split(",").some((v) => {
      const [media, ...params] = v.trim().split(";");
      return (
        (media === type ||
          media === "*/*" ||
          media === `${type.split("/")[0]}/*`) &&
        !params.some((p) => /^\s*q=0(?:\.0*)?\s*$/.test(p))
      );
    })
  );
}
export function makeServer(provider, options = {}) {
  const config = { ...defaults, ...options };
  const mcp = makeMcpHandler(provider);
  let active = 0;
  let global = { start: Date.now(), count: 0 };
  const clients = new Map();
  const server = createServer(
    { maxHeaderSize: 8192, requestTimeout: 6000, headersTimeout: 6000 },
    async (req, res) => {
      const id = randomUUID();
      let timer;
      let admitted = false;
      const abort = new AbortController();
      res.once("close", () => {
        if (!res.writableEnded) abort.abort();
      });
      res.setHeader("X-Request-ID", id);
      res.setHeader("X-Content-Type-Options", "nosniff");
      res.setHeader("Cache-Control", "no-store");
      const send = (
        status,
        data,
        {
          immutable = false,
          type = "application/json; charset=utf-8",
          budget = LIMITS.json,
        } = {},
      ) => {
        if (res.destroyed || res.writableEnded) return;
        const body = Buffer.from(
          typeof data === "string" ? data : JSON.stringify(data),
        );
        if (body.length > budget)
          throw new ApiError(
            413,
            "result_too_large",
            "Use the complete released records download.",
            { download: provider.fallback(data?.release_id, data?.species) },
          );
        res.statusCode = status;
        res.setHeader("Content-Type", type);
        res.setHeader("Content-Length", body.length);
        if (immutable) {
          res.setHeader("Cache-Control", "public, max-age=31536000, immutable");
          res.setHeader("ETag", hash(body));
          if (
            req.headers["if-none-match"]
              ?.split(",")
              .map((s) => s.trim())
              .some(
                (v) => v === hash(body) || v === `W/${hash(body)}` || v === "*",
              )
          ) {
            res.statusCode = 304;
            res.removeHeader("Content-Length");
            res.end();
            return;
          }
        }
        res.end(req.method === "HEAD" ? undefined : body);
      };
      const fail = (error) => {
        if (res.headersSent) {
          res.destroy();
          return;
        }
        res.removeHeader("Content-Disposition");
        res.removeHeader("Link");
        const known = error instanceof ApiError;
        const status = known ? error.status : 503;
        if (
          status === 429 ||
          (known && error.code === "publication_in_progress")
        )
          res.setHeader("Retry-After", "60");
        res.setHeader("Cache-Control", "no-store");
        send(status, {
          error: {
            code: known ? error.code : "unavailable",
            message: known
              ? error.message
              : "The released data is unavailable.",
            request_id: id,
            ...(known && error.details ? error.details : {}),
          },
        });
      };
      try {
        if (Buffer.byteLength(req.url ?? "") > LIMITS.url)
          throw new ApiError(
            413,
            "request_too_large",
            "The URL exceeds the request limit.",
          );
        if (
          !req.url.startsWith("/") ||
          req.url.startsWith("//") ||
          !/^(\[[A-Fa-f0-9:.]+\]|[A-Za-z0-9.-]+)(:[0-9]{1,5})?$/.test(
            req.headers.host ?? "",
          ) ||
          !req.headers.host ||
          req.headers.host.includes("@") ||
          req.rawHeaders.filter(
            (v, i) => i % 2 === 0 && v.toLowerCase() === "host",
          ).length !== 1
        )
          throw new ApiError(
            403,
            "invalid_host",
            "The request host is not allowed.",
          );
        let host;
        try {
          host = new URL(`http://${req.headers.host}`).hostname;
        } catch {
          throw new ApiError(
            403,
            "invalid_host",
            "The request host is not allowed.",
          );
        }
        if (!config.hosts.includes(host))
          throw new ApiError(
            403,
            "invalid_host",
            "The request host is not allowed.",
          );
        const origin = req.headers.origin;
        if (origin && !config.origins.includes(origin))
          throw new ApiError(
            403,
            "invalid_origin",
            "The request origin is not allowed.",
          );
        if (origin) {
          res.setHeader("Access-Control-Allow-Origin", origin);
          res.setHeader("Vary", "Origin");
          res.setHeader(
            "Access-Control-Expose-Headers",
            "ETag, X-Request-ID, Retry-After, Link, Content-Disposition",
          );
        }
        const url = new URL(req.url, `http://${req.headers.host}`);
        if (
          req.headers["content-encoding"] &&
          req.headers["content-encoding"] !== "identity"
        )
          throw new ApiError(
            415,
            "unsupported_media_type",
            "Compressed request bodies are not supported.",
          );
        if (req.method === "OPTIONS") {
          if (!origin) bad();
          const method = req.headers["access-control-request-method"];
          const allowed = url.pathname === "/mcp" ? ["POST"] : ["GET", "HEAD"];
          const headers = (req.headers["access-control-request-headers"] ?? "")
            .split(",")
            .map((s) => s.trim().toLowerCase())
            .filter(Boolean);
          if (
            !allowed.includes(method) ||
            headers.some((h) => !allowHeaders.includes(h))
          )
            bad();
          res.setHeader("Access-Control-Allow-Methods", allowed.join(", "));
          res.setHeader(
            "Access-Control-Allow-Headers",
            allowHeaders.join(", "),
          );
          res.statusCode = 204;
          res.end();
          return;
        }
        const now = Date.now();
        if (now - global.start >= 60000) {
          global = { start: now, count: 0 };
          for (const [key, b] of clients)
            if (now - b.start >= 60000) clients.delete(key);
        }
        const key = req.socket.remoteAddress ?? "unknown";
        let bucket = clients.get(key);
        if (!bucket || now - bucket.start >= 60000) {
          if (clients.size >= 10000 && !clients.has(key))
            throw new ApiError(
              429,
              "rate_limited",
              "The request rate limit is reached.",
            );
          bucket = { start: now, count: 0 };
          clients.set(key, bucket);
        }
        if (
          ++bucket.count > config.clientRate ||
          ++global.count > config.globalRate ||
          active >= config.concurrency
        )
          throw new ApiError(
            429,
            "rate_limited",
            "The request admission limit is reached.",
          );
        active++;
        admitted = true;
        timer = setTimeout(() => {
          fail(
            new ApiError(
              503,
              "handler_timeout",
              "The request exceeded its time budget.",
            ),
          );
          abort.abort();
          req.destroy();
        }, config.handlerMs);
        timer.unref();
        if (url.pathname === "/mcp") {
          argumentsFrom(url, []);
          if (req.method !== "POST") {
            res.setHeader("Allow", "POST, OPTIONS");
            throw new ApiError(
              405,
              "method_not_allowed",
              "Use POST for this stateless endpoint.",
            );
          }
          if (Number(req.headers["content-length"] ?? 0) > LIMITS.body)
            throw new ApiError(
              413,
              "request_too_large",
              "The request body exceeds the limit.",
            );
          if (!provider.ready)
            throw new ApiError(
              503,
              "unavailable",
              "The validated release bundle is unavailable.",
            );
          await mcp(req, res);
          return;
        }
        if (!["GET", "HEAD"].includes(req.method)) {
          res.setHeader("Allow", "GET, HEAD, OPTIONS");
          throw new ApiError(
            405,
            "method_not_allowed",
            "Use GET or HEAD for public data.",
          );
        }
        if (
          req.headers["transfer-encoding"] ||
          Number(req.headers["content-length"] ?? 0) !== 0
        )
          bad(
            "body_not_allowed",
            "Public data reads do not accept a request body.",
          );
        if (url.pathname === "/healthz") {
          argumentsFrom(url, []);
          send(200, { status: "ok" });
          return;
        }
        if (url.pathname === "/readyz") {
          argumentsFrom(url, []);
          send(provider.ready ? 200 : 503, {
            status: provider.ready ? "ready" : "unavailable",
          });
          return;
        }
        if (url.pathname === "/openapi.json") {
          argumentsFrom(url, []);
          send(200, openapiFor(provider));
          return;
        }
        const parts = url.pathname
          .split("/")
          .slice(1)
          .map((v) => {
            try {
              return decodeURIComponent(v);
            } catch {
              bad();
            }
          });
        if (parts[0] !== "v1" || parts[1] !== "releases")
          throw new ApiError(
            404,
            "not_found",
            "The requested route does not exist.",
          );
        const read = (name, args) => provider.call(name, args, abort.signal);
        let value;
        let immutable = false;
        let isSequence = false;
        let format = "json";
        if (parts.length === 2) {
          argumentsFrom(url, []);
          value = await read("get_release", {});
        } else if (parts.length === 3) {
          argumentsFrom(url, []);
          value = await read(
            "get_release",
            parts[2] === "current" ? {} : { release_id: parts[2] },
          );
          // Discovery metadata includes the current marker and retained release list.
          immutable = false;
        } else {
          const [, , release_id, species, resource, target_id, ...rest] = parts;
          const base = { release_id, species };
          immutable = true;
          if (resource === "downloads" && parts.length === 5) {
            argumentsFrom(url, []);
            value = await read("get_downloads", base);
          } else if (resource === "targets" && parts.length === 5)
            value = await read("search_targets", {
              ...base,
              ...argumentsFrom(url, [
                "q",
                "track",
                "status",
                "structure",
                "constructs",
                "disease_id",
                "limit",
                "cursor",
              ]),
            });
          else if (resource === "targets" && parts.length === 6) {
            argumentsFrom(url, []);
            value = await read("get_target", { ...base, target_id });
          } else if (
            resource === "targets" &&
            rest[0] === "sections" &&
            rest.length === 2
          )
            value = await read("get_target", {
              ...base,
              target_id,
              section: rest[1],
              ...argumentsFrom(url, ["limit", "cursor"]),
            });
          else if (
            resource === "targets" &&
            rest[0] === "constructs" &&
            rest.length === 1
          )
            value = await read("list_constructs", {
              ...base,
              target_id,
              ...argumentsFrom(url, ["category", "limit", "cursor"]),
            });
          else if (
            resource === "targets" &&
            ((rest.length === 1 && rest[0] === "sequence") ||
              (rest.length === 3 &&
                rest[0] === "constructs" &&
                rest[2] === "sequence"))
          ) {
            const args = argumentsFrom(url, ["format"]);
            format = args.format ?? "json";
            if (!["json", "fasta"].includes(format)) bad();
            value = await read("get_sequence", {
              ...base,
              target_id,
              ...(rest.length === 3 ? { construct_id: rest[1] } : {}),
            });
            isSequence = true;
          } else
            throw new ApiError(
              404,
              "not_found",
              "The requested route does not exist.",
            );
        }
        const type = format === "fasta" ? "text/x-fasta" : "application/json";
        if (!accepts(req, type))
          throw new ApiError(
            406,
            "not_acceptable",
            "The requested response media type is not supported.",
          );
        let representation = value;
        if (format === "fasta") {
          const download = provider.fallback(value.release_id, value.species);
          res.setHeader(
            "Link",
            `<${value.manifest.url}>; rel="describedby"; type="application/json", <${download.url}>; rel="alternate"; type="application/gzip"`,
          );
          const identity = (value.data.construct_id ?? "reference").replace(
            /[^A-Za-z0-9_.-]/g,
            "_",
          );
          res.setHeader(
            "Content-Disposition",
            `attachment; filename="openantigens-${value.release_id}-${value.species}-${identity}.fasta"`,
          );
          representation = fasta(value.data, {
            release_id: value.release_id,
            provenance_url: value.manifest.url,
            records_url: download.url,
          });
        }
        send(200, representation, {
          immutable,
          type: `${type}; charset=utf-8`,
          budget: isSequence ? LIMITS.sequence : LIMITS.json,
        });
      } catch (error) {
        fail(error);
      } finally {
        clearTimeout(timer);
        if (admitted) active--;
      }
    },
  );
  server.maxRequestsPerSocket = 1000;
  server.timeout = config.handlerMs + 1000;
  return server;
}
export async function start() {
  const provider = new Provider(
    process.env.OPENANTIGENS_DATA_DIR ?? resolve("data"),
  );
  provider.publicationState = process.env.OPENANTIGENS_PUBLICATION_STATE;
  try {
    await provider.initialize();
  } catch {
    process.stderr.write("public_api bundle_validation_failed\n");
  }
  const csv = (name, fallback) =>
    process.env[name]
      ? process.env[name]
          .split(",")
          .map((v) => v.trim())
          .filter(Boolean)
      : fallback;
  const server = makeServer(provider, {
    hosts: csv("OPENANTIGENS_ALLOWED_HOSTS", defaults.hosts),
    origins: csv("OPENANTIGENS_ALLOWED_ORIGINS", defaults.origins),
  });
  const port = Number(process.env.PORT ?? 3000);
  if (!Number.isInteger(port) || port < 1 || port > 65535)
    throw new Error("Invalid service port");
  server.listen(port, process.env.OPENANTIGENS_HOST ?? "0.0.0.0", () =>
    process.stdout.write(
      `public_api listening port=${port} ready=${provider.ready}\n`,
    ),
  );
  const stop = () => {
    server.close();
    setTimeout(() => server.closeAllConnections(), 5000).unref();
  };
  process.once("SIGTERM", stop);
  process.once("SIGINT", stop);
  return server;
}
if (
  process.argv[1] &&
  import.meta.url === pathToFileURL(resolve(process.argv[1])).href
)
  start().catch(() => {
    process.stderr.write("public_api startup_failed\n");
    process.exitCode = 1;
  });

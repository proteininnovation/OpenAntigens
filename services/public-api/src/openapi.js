import { z } from "zod";
import {
  inputs,
  provenanceOutputSchema,
  GUIDANCE,
  LIMITS,
  SECTIONS,
  CATEGORIES,
  TRACKS,
  STATUSES,
} from "./contracts.js";
const release = {
  name: "release_id",
  in: "path",
  required: true,
  schema: { type: "string", pattern: "^[a-f0-9]{64}$" },
};
const species = {
  name: "species",
  in: "path",
  required: true,
  schema: { type: "string", enum: ["human", "mouse"] },
};
const target = {
  name: "target_id",
  in: "path",
  required: true,
  schema: {
    type: "string",
    pattern: "^(9606|10090):(uniprot|refseq):[A-Za-z0-9_.-]+$",
  },
};
const query = (name, schema) => ({ name, in: "query", schema });
const pages = [
  query("limit", { type: "integer", minimum: 1, maximum: 100, default: 25 }),
  query("cursor", { type: "string", maxLength: 2048 }),
];
const response = {
  description:
    "Released public data. Collection envelopes contain total and next_cursor; all envelopes include release_id, schema_version, manifest and provenance.",
  content: {
    "application/json": { schema: { $ref: "#/components/schemas/Envelope" } },
  },
  headers: {
    ETag: { schema: { type: "string" } },
    "X-Request-ID": { schema: { type: "string" } },
  },
};
const errors = Object.fromEntries(
  [400, 403, 404, 405, 406, 410, 413, 415, 429, 503].map((code) => [
    code,
    {
      description: {
        400: "invalid_arguments, invalid_cursor, species_mismatch or unsupported_capability",
        403: "invalid_host or invalid_origin",
        404: "not_found",
        405: "method_not_allowed",
        406: "not_acceptable",
        410: "retired_release (only with a published tombstone)",
        413: "request_too_large or result_too_large; data overflow supplies a complete records download",
        415: "unsupported_media_type",
        429: "rate_limited, Retry-After: 60",
        503: "unavailable or handler_timeout",
      }[code],
      content: {
        "application/json": { schema: { $ref: "#/components/schemas/Error" } },
      },
    },
  ]),
);
const operation = (operationId, summary, parameters = []) => ({
  get: {
    operationId,
    summary,
    description: GUIDANCE,
    parameters,
    responses: {
      200: structuredClone(response),
      304: { description: "ETag matches; empty body." },
      ...errors,
    },
  },
  head: {
    summary: `Headers for ${summary}`,
    parameters,
    responses: {
      200: { description: "Same representation headers as GET; empty body." },
      304: { description: "ETag matches; empty body." },
      ...errors,
    },
  },
});
const base = "/v1/releases/{release_id}/{species}";
const detail = `${base}/targets/{target_id}`;
export const openapi = {
  openapi: "3.1.0",
  info: {
    title: "OpenAntigens public read-only API",
    version: "1.0.0",
    description: GUIDANCE,
  },
  servers: [{ url: "https://api.openantigens.org" }],
  externalDocs: {
    url: "https://openantigens.org/agent-guide.html",
    description: "Scientific interpretation and citation guidance",
  },
  "x-agent-guides": {
    human: "https://openantigens.org/agent-guide.html",
    mouse: "https://openantigens.org/mouse/agent-guide.html",
  },
  "x-limits": LIMITS,
  "x-query-policy":
    "Unknown or duplicate query arguments are rejected. Literal release IDs required on data routes. Species and target ID must agree. Data is frozen and processing state is historical.",
  "x-mcp": {
    endpoint: "/mcp",
    protocol: "2026-07-28",
    legacy: "stateless 2025-11-25, one-shot SSE from official SDK",
    tools: Object.keys(inputs),
    sessions: false,
    annotations: {
      readOnlyHint: true,
      destructiveHint: false,
      idempotentHint: true,
      openWorldHint: false,
    },
  },
  paths: {
    "/healthz": operation("health", "Read process liveness"),
    "/readyz": operation("ready", "Read validated bundle readiness"),
    "/openapi.json": operation("openapi", "Read API description"),
    "/v1/releases": operation("releases", "Discover retained releases"),
    "/v1/releases/current": operation(
      "currentRelease",
      "Resolve current immutable release; Cache-Control: no-store",
    ),
    "/v1/releases/{release_id}": operation(
      "getRelease",
      "Read immutable release metadata",
      [release],
    ),
    [`${base}/targets`]: operation(
      "searchTargets",
      "Search or list released catalogue",
      [
        release,
        species,
        query("q", {
          type: "string",
          maxLength: 200,
          description:
            "Case-insensitive substring over public aliases, accessions, names, family and topology.",
        }),
        query("track", { type: "string", enum: TRACKS }),
        query("status", { type: "string", enum: STATUSES }),
        query("structure", {
          type: "string",
          enum: ["available", "unavailable", "unknown"],
        }),
        query("constructs", {
          type: "string",
          enum: ["available", "unavailable", "unknown"],
        }),
        query("disease_id", {
          type: "string",
          minLength: 1,
          maxLength: 128,
          description:
            "Exact recorded association source ID; human capability only.",
        }),
        ...pages,
      ],
    ),
    [detail]: operation("getTarget", "Read one bounded target summary", [
      release,
      species,
      target,
    ]),
    [`${detail}/sections/{section}`]: operation(
      "getSection",
      "Read one evidence section page; similarity limit maximum 50",
      [
        release,
        species,
        target,
        {
          name: "section",
          in: "path",
          required: true,
          schema: { type: "string", enum: SECTIONS },
        },
        ...pages,
      ],
    ),
    [`${detail}/constructs`]: operation(
      "listConstructs",
      "Read released constructs",
      [
        release,
        species,
        target,
        query("category", { type: "string", enum: CATEGORIES }),
        ...pages,
      ],
    ),
    [`${detail}/sequence`]: operation(
      "getReferenceSequence",
      "Read exact released reference sequence",
      [
        release,
        species,
        target,
        query("format", {
          type: "string",
          enum: ["json", "fasta"],
          default: "json",
        }),
      ],
    ),
    [`${detail}/constructs/{construct_id}/sequence`]: operation(
      "getConstructSequence",
      "Read exact released construct sequence",
      [
        release,
        species,
        target,
        {
          name: "construct_id",
          in: "path",
          required: true,
          schema: { type: "string", pattern: "^[a-f0-9]{64}$" },
        },
        query("format", {
          type: "string",
          enum: ["json", "fasta"],
          default: "json",
        }),
      ],
    ),
    [`${base}/downloads`]: operation(
      "getDownloads",
      "Discover complete immutable bulk records, index, constructs, sequences, schemas and citations",
      [release, species],
    ),
    "/mcp": {
      post: {
        summary: "Stateless tools-only MCP endpoint",
        description:
          "Official SDK validates current per-request negotiation and legacy initialization. Current responses use application/json; stateless legacy responses can use finite text/event-stream. No GET streams or session store. Host, Origin, size and admission guards apply.",
        requestBody: {
          required: true,
          content: { "application/json": { schema: { type: "object" } } },
        },
        responses: {
          200: {
            description: "Official MCP JSON or finite legacy SSE response",
          },
          202: { description: "Notification accepted" },
          ...errors,
        },
      },
    },
  },
  components: {
    schemas: {
      Envelope: {
        type: "object",
        required: [
          "release_id",
          "schema_version",
          "data",
          "manifest",
          "provenance",
        ],
        properties: {
          release_id: { type: "string", pattern: "^[a-f0-9]{64}$" },
          schema_version: { const: 1 },
          species: { enum: ["human", "mouse"] },
          taxon_id: { enum: [9606, 10090] },
          data: {
            description:
              "Exact versioned release schema, supplied in the immutable schema downloads.",
          },
          manifest: {
            type: "object",
            required: ["url", "sha256"],
            additionalProperties: false,
            properties: {
              url: { type: "string", format: "uri" },
              sha256: { type: "string", pattern: "^[a-f0-9]{64}$" },
            },
          },
          provenance: { $ref: "#/components/schemas/ProvenanceOutput" },
          next_cursor: { type: ["string", "null"] },
          total: { type: "integer", minimum: 0 },
          state: {
            enum: [
              "available",
              "no_results",
              "unavailable",
              "not_applicable",
              "error",
              "unknown",
            ],
          },
          reason: { type: ["string", "null"] },
          current_release_id: { type: "string" },
          releases: { type: "array", items: { type: "object" } },
          retired_releases: {
            type: "array",
            maxItems: 1000,
            items: { $ref: "#/components/schemas/Retirement" },
          },
          discovery_state: { enum: ["active", "publication_in_progress"] },
        },
        additionalProperties: false,
      },
      ProvenanceOutput: z.toJSONSchema(provenanceOutputSchema),
      Retirement: {
        type: "object",
        required: ["release_id", "retired_by_release_id", "archive_url"],
        additionalProperties: false,
        properties: {
          release_id: { type: "string", pattern: "^[a-f0-9]{64}$" },
          retired_by_release_id: { type: "string", pattern: "^[a-f0-9]{64}$" },
          archive_url: { type: ["string", "null"], format: "uri" },
        },
      },
      OpenApiDocument: {
        type: "object",
        required: ["openapi", "info", "paths"],
        properties: {
          openapi: { const: "3.1.0" },
          info: {
            type: "object",
            required: ["title", "version"],
            properties: {
              title: { type: "string" },
              version: { type: "string" },
            },
          },
          paths: { type: "object" },
          components: { type: "object" },
        },
        additionalProperties: true,
      },
      Error: {
        type: "object",
        required: ["error"],
        additionalProperties: false,
        properties: {
          error: {
            type: "object",
            required: ["code", "message", "request_id"],
            additionalProperties: false,
            properties: {
              code: { type: "string" },
              message: { type: "string" },
              request_id: { type: "string", format: "uuid" },
              download: { type: "object" },
              release_id: { type: "string", pattern: "^[a-f0-9]{64}$" },
              retirement: { $ref: "#/components/schemas/Retirement" },
            },
          },
        },
      },
      ...Object.fromEntries(
        Object.entries(inputs).map(([name, schema]) => [
          `${name}Input`,
          z.toJSONSchema(schema, { unrepresentable: "any" }),
        ]),
      ),
    },
  },
  "x-examples": {
    human_target_id: "9606:uniprot:P12345",
    mouse_target_id: "10090:refseq:NP_001234.1",
    coordinate_system:
      "1-based inclusive in the explicitly named sequence basis",
    missing_state: {
      state: "unknown",
      reason: "The source history does not distinguish empty from unqueried.",
    },
  },
};
for (const key of [
  `${detail}/sequence`,
  `${detail}/constructs/{construct_id}/sequence`,
])
  openapi.paths[key].get.responses["200"].content["text/x-fasta"] = {
    schema: { type: "string" },
  };

const statusResponse = (status) => ({
  description: `Process status: ${status}.`,
  headers: {
    "Cache-Control": { schema: { const: "no-store" } },
    "X-Request-ID": { schema: { type: "string", format: "uuid" } },
  },
  content: {
    "application/json": {
      schema: {
        type: "object",
        required: ["status"],
        additionalProperties: false,
        properties: { status: { const: status } },
      },
    },
  },
});
for (const path of ["/healthz", "/readyz", "/openapi.json"]) {
  const item = openapi.paths[path];
  delete item.get.responses[304];
  delete item.head.responses[304];
  item.get.responses[200] =
    path === "/openapi.json"
      ? {
          description: "This OpenAPI document.",
          headers: statusResponse("ok").headers,
          content: {
            "application/json": {
              schema: { $ref: "#/components/schemas/OpenApiDocument" },
            },
          },
        }
      : statusResponse(path === "/healthz" ? "ok" : "ready");
  if (path === "/readyz")
    item.get.responses[503] = {
      description:
        "The validated bundle is unavailable, or a request-level failure occurred.",
      headers: statusResponse("unavailable").headers,
      content: {
        "application/json": {
          schema: {
            oneOf: [
              statusResponse("unavailable").content["application/json"].schema,
              { $ref: "#/components/schemas/Error" },
            ],
          },
        },
      },
    };
}
for (const path of [
  "/v1/releases",
  "/v1/releases/current",
  "/v1/releases/{release_id}",
]) {
  const item = openapi.paths[path];
  delete item.get.responses[304];
  delete item.head.responses[304];
  item.get.responses[200].headers = {
    "Cache-Control": { schema: { const: "no-store" } },
    "X-Request-ID": { schema: { type: "string", format: "uuid" } },
  };
}
for (const path of [
  `${detail}/sequence`,
  `${detail}/constructs/{construct_id}/sequence`,
])
  Object.assign(openapi.paths[path].get.responses[200].headers, {
    Link: {
      description: "Immutable manifest metadata and complete records download.",
      schema: { type: "string" },
    },
    "Content-Disposition": {
      description: "FASTA attachment with a sanitized filename.",
      schema: { type: "string" },
    },
  });

export function openapiFor(provider) {
  const result = structuredClone(openapi);
  if (!provider.ready) return result;
  const tools = Object.fromEntries(
    Object.keys(inputs).map((name) => [
      name,
      z.toJSONSchema(provider.toolOutput(name)),
    ]),
  );
  for (const [name, schema] of Object.entries(tools))
    result.components.schemas[`${name}Output`] = schema;
  const operations = {
    releases: "get_release",
    currentRelease: "get_release",
    getRelease: "get_release",
    searchTargets: "search_targets",
    getTarget: "get_target",
    getSection: "get_target",
    listConstructs: "list_constructs",
    getReferenceSequence: "get_sequence",
    getConstructSequence: "get_sequence",
    getDownloads: "get_downloads",
  };
  for (const item of Object.values(result.paths)) {
    const name = operations[item.get?.operationId];
    if (name)
      item.get.responses["200"].content["application/json"].schema = {
        $ref: `#/components/schemas/${name}Output`,
      };
  }
  return result;
}

for (const item of Object.values(openapi.paths)) {
  if (item.head)
    for (const [status, response] of Object.entries(item.head.responses)) {
      const headersOnly = structuredClone(
        item.get.responses[status] ?? response,
      );
      delete headersOnly.content;
      item.head.responses[status] = headersOnly;
    }
}

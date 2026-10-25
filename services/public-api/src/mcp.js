import { McpServer, createMcpHandler } from "@modelcontextprotocol/server";
import { toNodeHandler } from "@modelcontextprotocol/node";
import { randomUUID } from "node:crypto";
import { AsyncLocalStorage } from "node:async_hooks";
import { inputs, GUIDANCE, LIMITS, ApiError } from "./contracts.js";
const descriptions = {
  get_release:
    "Read immutable release metadata, retained releases, capabilities, sources and citations.",
  search_targets:
    "Find up to 20 matching target candidates in one explicit species and release using at least two literal query characters. Refine the query when more matches exist. Aliases can yield several candidates; unresolved records retain null target IDs.",
  get_target:
    "Read one target summary or one explicit paged evidence section. Family evidence projects the stored report family_context; portal evidence views provide full-length family matrices and construct homolog exports. Availability states preserve missingness; human disease evidence requires the declared capability.",
  list_constructs:
    "Read released construct classes, boundaries, exact sequence identities, recorded scores, evidence and warnings. Classification and category remain distinct.",
  get_sequence:
    "Read one exact released reference or construct sequence with the actual sequence accession/source, coordinates, checksum and warnings. Coordinates are 1-based inclusive in the declared sequence basis.",
};
// Inspect only bounded terminal messages; forward every SDK response byte unchanged.
function observeSdkResponse(response, exchange, observe) {
  const type = response.headers.get("content-type")?.split(";")[0];
  if (
    !response.body ||
    !["application/json", "text/event-stream"].includes(type)
  )
    return response;
  let bytes = 0,
    complete = true;
  const chunks = [];
  const inspect = (message) => {
    if (
      (message?.error || message?.result?.isError === true) &&
      !exchange.error
    ) {
      exchange.error = exchange.entered ? "unavailable" : "invalid_arguments";
      observe({ error_code: exchange.error });
    }
  };
  const observer = new TransformStream({
    transform(chunk, controller) {
      bytes += chunk.byteLength;
      if (complete && bytes <= LIMITS.sequence + LIMITS.json)
        chunks.push(Buffer.from(chunk));
      else {
        complete = false;
        chunks.length = 0;
      }
      controller.enqueue(chunk);
    },
    flush() {
      if (!complete) return;
      try {
        const text = new TextDecoder("utf-8", { fatal: true }).decode(
          Buffer.concat(chunks),
        );
        if (!text.trim()) return;
        if (type === "application/json") inspect(JSON.parse(text));
        else
          for (const event of text.split(/\r?\n\r?\n/)) {
            const data = event
              .split(/\r?\n/)
              .filter((line) => line.startsWith("data:"))
              .map((line) => line.slice(5).replace(/^ /, ""))
              .join("\n");
            if (data) inspect(JSON.parse(data));
          }
      } catch {
        if (!exchange.error) observe({ error_code: "unavailable" });
      }
    },
  });
  return new Response(response.body.pipeThrough(observer), {
    status: response.status,
    statusText: response.statusText,
    headers: response.headers,
  });
}
export function makeMcpHandler(provider, observe = () => {}) {
  const exchanges = new AsyncLocalStorage();
  const handler = createMcpHandler(
    () => {
      const server = new McpServer(
        { name: "openantigens", version: "1.0.0" },
        { instructions: GUIDANCE },
      );
      for (const [name, inputSchema] of Object.entries(inputs))
        server.registerTool(
          name,
          {
            description: `${descriptions[name]} ${GUIDANCE}`,
            inputSchema,
            outputSchema: provider.toolOutput(name),
            annotations: {
              readOnlyHint: true,
              destructiveHint: false,
              idempotentHint: true,
              openWorldHint: false,
            },
          },
          async (args, context) => {
            const exchange = exchanges.getStore();
            if (exchange) exchange.entered = true;
            observe({
              operation:
                name === "get_target" && args.section ? "get_section" : name,
              species: args.species ?? "none",
            });
            try {
              const value = await provider.call(name, args, context.signal);
              const budget =
                name === "get_sequence" ? LIMITS.sequence : LIMITS.json;
              if (Buffer.byteLength(JSON.stringify(value)) > budget)
                throw new ApiError(
                  413,
                  "result_too_large",
                  "Use the complete released records download.",
                  {
                    download: provider.fallback(args.release_id, args.species),
                  },
                );
              return {
                structuredContent: value,
                content: [
                  {
                    type: "text",
                    text: `Release ${value.release_id}${value.species ? `; species ${value.species}` : ""}. Sources and provenance are in structuredContent; availability states preserve missingness. ${GUIDANCE}`,
                  },
                ],
              };
            } catch (error) {
              const known = error instanceof ApiError;
              const error_code = known ? error.code : "unavailable";
              if (exchange) exchange.error = error_code;
              observe({ error_code });
              const failure = {
                code: known ? error.code : "unavailable",
                message: known
                  ? error.message
                  : "The released data could not be read.",
                request_id: randomUUID(),
                ...(known && error.details ? error.details : {}),
              };
              return {
                isError: true,
                content: [{ type: "text", text: JSON.stringify(failure) }],
              };
            }
          },
        );
      return server;
    },
    {
      responseMode: "json",
      legacy: "stateless",
      maxRequestBodySize: LIMITS.body,
      maxSubscriptions: 0,
      keepAliveMs: 0,
      onerror: () => observe({ error_code: "unknown" }),
    },
  );
  return toNodeHandler(
    {
      fetch: async (request, options) => {
        let parsed;
        try {
          parsed = await request.clone().json();
        } catch {}
        const protocol = {
          "server/discover": "discover",
          initialize: "initialize",
          "notifications/initialized": "initialized",
          ping: "ping",
          "tools/list": "list_tools",
          "tools/call": "call_tool",
        };
        observe({
          operation:
            typeof parsed?.method === "string" &&
            Object.hasOwn(protocol, parsed.method)
              ? protocol[parsed.method]
              : "unknown",
        });
        if (Array.isArray(parsed)) {
          observe({ error_code: "invalid_arguments" });

          return Response.json(
            {
              error: {
                code: "invalid_arguments",
                message: "JSON-RPC batches are not supported.",
                request_id: randomUUID(),
              },
            },
            { status: 400, headers: { "cache-control": "no-store" } },
          );
        }
        const exchange = { entered: false, error: null };
        const response = await exchanges.run(exchange, () =>
          handler.fetch(request, options),
        );
        return observeSdkResponse(response, exchange, observe);
      },
    },
    { maxRequestBodySize: LIMITS.body, onerror: () => {} },
  );
}

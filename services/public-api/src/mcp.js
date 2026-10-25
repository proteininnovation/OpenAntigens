import { McpServer, createMcpHandler } from "@modelcontextprotocol/server";
import { toNodeHandler } from "@modelcontextprotocol/node";
import { randomUUID } from "node:crypto";
import { inputs, GUIDANCE, LIMITS, ApiError } from "./contracts.js";
const descriptions = {
  get_release:
    "Read immutable release metadata, retained releases, capabilities, sources and citations.",
  search_targets:
    "Find released target candidates in one explicit species and release. Aliases can yield several candidates; unresolved records retain null target IDs.",
  get_target:
    "Read one target summary or one explicit paged evidence section. Family evidence projects the stored report family_context; portal evidence views provide full-length family matrices and construct homolog exports. Availability states preserve missingness; human disease evidence requires the declared capability.",
  list_constructs:
    "Read released construct classes, boundaries, exact sequence identities, recorded scores, evidence and warnings. Classification and category remain distinct.",
  get_sequence:
    "Read one exact released reference or construct sequence with the actual sequence accession/source, coordinates, checksum and warnings. Coordinates are 1-based inclusive in the declared sequence basis.",
  get_downloads:
    "Discover complete immutable bulk downloads and their sizes, SHA-256 checksums, sources and attribution.",
};
export function makeMcpHandler(provider) {
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
      onerror: () => {},
    },
  );
  return toNodeHandler(
    {
      fetch: async (request, options) => {
        let parsed;
        try {
          parsed = await request.clone().json();
        } catch {}
        if (Array.isArray(parsed))
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
        return handler.fetch(request, options);
      },
    },
    { maxRequestBodySize: LIMITS.body, onerror: () => {} },
  );
}

import { z } from "zod";
export const SECTIONS = [
  "topology",
  "domains",
  "ptms",
  "orthologs",
  "paralogs",
  "family",
  "assembly",
  "structures",
  "similarity",
  "literature",
  "disease",
];
export const CATEGORIES = [
  "full_region",
  "pdb",
  "domain_repeat",
  "strict",
  "lenient",
  "multipass_full_length",
  "multipass_trimmed",
];
export const TRACKS = [
  "secreted",
  "gpi",
  "single_pass",
  "multipass",
  "other",
  "unknown",
];
export const STATUSES = [
  "ok",
  "error",
  "pending",
  "running",
  "started",
  "skipped_existing",
  "unknown",
];
export const STATES = ["available", "unavailable", "unknown"];
export const GUIDANCE =
  "Processing success does not establish experimental validation. Model scores and predictions retain their recorded meaning. Sequence similarity does not establish specificity, epitope identity, or cross-reactivity. Empty results do not establish biological absence. Cite the release, source accession and supplied attribution.";
export const LIMITS = Object.freeze({
  url: 4096,
  body: 32768,
  json: 262144,
  sequence: 1048576,
  index: 33554432,
  indexes: 67108864,
  manifest: 67108864,
  target: 8388608,
  shard: 1048576,
  defaultPage: 25,
  maxPage: 100,
  similarity: 50,
  concurrency: 8,
  handlerMs: 5000,
});
const release = z.string().regex(/^[a-f0-9]{64}$/);
const species = z.enum(["human", "mouse"]);
export const tombstoneSchema = z
  .object({
    release_id: release,
    retired_by_release_id: release,
    archive_url: z.string().max(1024).url().nullable(),
  })
  .strict();
export const rootSchema = z
  .object({
    schema_version: z.literal(1),
    current_release_id: release,
    releases: z
      .array(
        z
          .object({
            release_id: release,
            manifest_path: z.string().min(1).max(1024),
            sha256: release,
          })
          .strict(),
      )
      .min(1)
      .max(2),
    retired_releases: z.array(tombstoneSchema).max(1000).optional(),
  })
  .strict();

const target_id = z
  .string()
  .regex(/^(9606|10090):(uniprot|refseq):[A-Za-z0-9_.-]+$/)
  .max(200);
const page = {
  limit: z.number().int().min(1).max(100).optional(),
  cursor: z.string().max(2048).optional(),
};
const filters = {
  q: z
    .string()
    .refine((v) => [...v].length <= 200)
    .optional(),
  track: z.enum(TRACKS).optional(),
  status: z.enum(STATUSES).optional(),
  structure: z.enum(STATES).optional(),
  constructs: z.enum(STATES).optional(),
  disease_id: z.string().min(1).max(128).optional(),
};
const base = { release_id: release, species };
export const inputs = {
  get_release: z.object({ release_id: release.optional() }).strict(),
  search_targets: z.object({ ...base, ...filters, ...page }).strict(),
  get_target: z
    .object({
      ...base,
      target_id,
      section: z.enum(SECTIONS).optional(),
      ...page,
    })
    .strict()
    .refine(
      (v) =>
        v.section !== undefined ||
        (v.limit === undefined && v.cursor === undefined),
      { message: "Page arguments require section" },
    ),
  list_constructs: z
    .object({
      ...base,
      target_id,
      category: z.enum(CATEGORIES).optional(),
      ...page,
    })
    .strict(),
  get_sequence: z
    .object({ ...base, target_id, construct_id: release.optional() })
    .strict(),
  get_downloads: z.object(base).strict(),
};
export const provenanceOutputSchema = z
  .object({
    software_commit: z.string().max(256),
    snapshot_id: z.string().max(256),
    source_versions: z.record(z.string(), z.string().max(16384).nullable()),
    source_dates: z.record(z.string(), z.string().max(16384).nullable()),
    notes: z.array(z.string().max(16384)),
    input_hashes: z
      .object({
        count: z.number().int().nonnegative(),
        manifest_url: z.string().max(1024).url(),
        manifest_sha256: release,
        json_pointer: z.literal("/provenance/input_sha256"),
      })
      .strict(),
  })
  .strict();
// Output references the complete manifest provenance without repeating its input map.
export const outputSchema = z
  .object({
    release_id: release,
    schema_version: z.literal(1),
    species: species.optional(),
    taxon_id: z.number().int().optional(),
    data: z.unknown(),
    manifest: z.object({ url: z.string(), sha256: release }),
    provenance: provenanceOutputSchema,
    next_cursor: z.string().nullable().optional(),
    total: z.number().int().nonnegative().optional(),
    state: z.string().optional(),
    reason: z.string().nullable().optional(),
    current_release_id: release.optional(),
    releases: z.array(z.unknown()).optional(),
    retired_releases: z.array(tombstoneSchema).max(1000).optional(),
    discovery_state: z.enum(["active", "publication_in_progress"]).optional(),
  })
  .strict();
export class ApiError extends Error {
  constructor(status, code, message, details) {
    super(message);
    this.status = status;
    this.code = code;
    this.details = details;
  }
}
export function bad(
  code = "invalid_arguments",
  message = "The request arguments are invalid.",
) {
  throw new ApiError(400, code, message);
}
export function parseInput(name, args) {
  const result = inputs[name].safeParse(args);
  if (!result.success) bad();
  return result.data;
}

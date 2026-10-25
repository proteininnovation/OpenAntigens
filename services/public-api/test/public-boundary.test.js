import test from "node:test";
import assert from "node:assert/strict";
import { readFile, writeFile, rm } from "node:fs/promises";
import { join } from "node:path";
import { createHash } from "node:crypto";
import {
  isPublicArchiveUrl,
  validateRetirementMetadata,
} from "../src/public-boundary.js";
import { Provider } from "../src/provider.js";
import { fixture, packedFixture } from "./fixture.js";
import { queryFixture } from "./query-fixture.js";
const hash = (value) => createHash("sha256").update(value).digest("hex");
const bad = [
  "https://archive.org/PRIVATE_CANARY_RETIREMENT",
  "https://archive.org/%2550RIVATE_CANARY_RETIREMENT",
  "https://archive.org/%2FUsers%2Foperator%2Fdata",
  "http://archive.org/details/record",
  "HTTPS://archive.org/details/record",
  "https://user:password@archive.org/details/record",
  "https://@archive.org/details/record",
  "https://:@archive.org/details/record",
  "https://127.0.0.1/details/record",
  "https://[::1]/record",
  "https://localhost/record",
  "https://archive.localdomain/record",
  "https://records.home.arpa/record",
  "https://archive.example/record",
  "https://archive.org/details/record?token=credential",
  "https://archive.org/details/record?",
  "https://archive.org/details/record#",
  "https://archive.org/details/%253Ftoken%253Dcredential",
  "https://archive.org/a/../record",
  "https://archive.org/a/%252e%252e/record",
  "https://archive.org/a%255Cb",
  "https://archive.org/%ZZ",
  "https://archive.org/%C0%AF",
  "https://archive.org/%C2%A0",
  "https://archive.org/\u0085",
  "https://archive.org:444/record",
  "https://archive.org/" + "a".repeat(1024),
];
test("one shared archive guard allows public permalinks and rejects private or ambiguous links", () => {
  for (const value of [
    null,
    "https://archive.org/details/record",
    "https://zenodo.org:443/records/123",
    "https://doi.org/10.64898/2026.07.30.741735",
    "https://Archive.ORG/details/secreted-protein-study",
    "https://archive.org/encoded%25percent",
  ])
    assert.equal(isPublicArchiveUrl(value), true, String(value));
  for (const value of bad)
    assert.equal(
      isPublicArchiveUrl(value),
      false,
      "Unsafe archive link must be rejected",
    );
  assert.doesNotThrow(() =>
    validateRetirementMetadata({
      retired_releases: [],
      note: "A secreted protein study.",
    }),
  );
  assert.throws(
    () =>
      validateRetirementMetadata({
        retired_releases: [],
        note: "PRIVATE_CANARY hidden metadata",
      }),
    { message: "Invalid public retirement metadata." },
  );
});
async function sourceWithArchive(t, archive_url) {
  const source = await fixture("retirement-boundary");
  t.after(() => rm(source.directory, { recursive: true, force: true }));
  const path = join(source.directory, "manifest.json");
  const root = JSON.parse(await readFile(path));
  const retired = {
    release_id: hash("retired-source"),
    retired_by_release_id: source.release_id,
    archive_url,
  };
  root.retired_releases = [retired];
  await writeFile(path, JSON.stringify(root));
  return { source, retired };
}
test("resealed unsafe retirement links fail old, packed and own-static startup without raw metadata errors", async (t) => {
  for (const archive_url of [bad[0], bad[1], bad[5], bad[14]]) {
    const { source } = await sourceWithArchive(t, archive_url);
    const packed = await packedFixture(source.directory),
      query = await queryFixture(source.directory);
    t.after(async () => {
      await rm(packed.directory, { recursive: true, force: true });
      await rm(query.directory, { recursive: true, force: true });
      await rm(query.staticDirectory, { recursive: true, force: true });
    });
    let fetches = 0;
    const providers = [
      new Provider(source.directory),
      new Provider(packed.directory, packed.pins),
      new Provider(query.directory, {
        querySeal: query.seal,
        fetchImpl: () => {
          fetches++;
          throw Error("Forbidden origin call");
        },
      }),
    ];
    for (const provider of providers) {
      await assert.rejects(
        provider.initialize(),
        (error) =>
          !String(error).includes(archive_url) &&
          !String(error).includes("PRIVATE_CANARY"),
      );
      assert.equal(provider.ready, false);
    }
    assert.equal(fetches, 0);
  }
});
test("null and safe archived links preserve the literal 410 metadata in all runtime modes", async (t) => {
  for (const archive_url of [
    null,
    "https://archive.org/details/public-release",
    "https://zenodo.org:443/records/123",
  ]) {
    const { source, retired } = await sourceWithArchive(t, archive_url);
    const packed = await packedFixture(source.directory),
      query = await queryFixture(source.directory);
    t.after(async () => {
      await rm(packed.directory, { recursive: true, force: true });
      await rm(query.directory, { recursive: true, force: true });
      await rm(query.staticDirectory, { recursive: true, force: true });
    });
    for (const provider of [
      new Provider(source.directory),
      new Provider(packed.directory, packed.pins),
      new Provider(query.directory, {
        querySeal: query.seal,
        fetchImpl: () => {
          throw Error("No origin");
        },
      }),
    ]) {
      await provider.initialize();
      await assert.rejects(
        provider.call("get_release", { release_id: retired.release_id }),
        (error) =>
          error.status === 410 &&
          isPublicArchiveUrl(error.details.retirement.archive_url) &&
          error.details.retirement.archive_url === archive_url,
      );
      assert.equal(provider.ready, true);
    }
  }
});

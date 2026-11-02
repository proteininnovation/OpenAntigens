import test from "node:test";
import assert from "node:assert/strict";
import {
  mkdtemp,
  mkdir,
  readFile,
  writeFile,
  rm,
  symlink,
} from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, dirname } from "node:path";
import { createHash } from "node:crypto";
import { gzipSync } from "node:zlib";
import { createStorage } from "../src/storage.js";

const hash = (bytes) => createHash("sha256").update(bytes).digest("hex");
const canonical = (value) =>
  Buffer.from(
    JSON.stringify(value, (_, item) =>
      item && typeof item === "object" && !Array.isArray(item)
        ? Object.fromEntries(
            Object.keys(item)
              .sort()
              .map((key) => [key, item[key]]),
          )
        : item,
    ) + "\n",
  );
const FORMAT = "openantigens-target-gzip-v1";

async function fixture(t) {
  const directory = await mkdtemp(join(tmpdir(), "oa-storage-"));
  t.after(() => rm(directory, { recursive: true, force: true }));
  const rid = "a".repeat(64),
    key = "b".repeat(64),
    cid = "c".repeat(64);
  const prefix = `releases/${rid}/human`;
  const target = `${prefix}/targets/${key}.json`;
  const construct = `${prefix}/constructs/${key}-0.json`;
  const sequence = `${prefix}/sequences/${cid}.json`;
  const reference = `${prefix}/sequences/${key}-reference.json`;
  const members = new Map([
    [construct, Buffer.from('{"data":[{"value":1.0}]}\n')],
    [sequence, Buffer.from('{"data":{"value":1e-06}}\n')],
    [reference, canonical({ data: { sequence: "AAA" } })],
    [
      target,
      canonical({
        data: {
          sections: {},
          constructs: { pages: [construct] },
          construct_sequence_paths: { [cid]: sequence },
          reference_sequence_path: reference,
        },
      }),
    ],
  ]);
  const manifest = {
    release_id: rid,
    files: Object.fromEntries(
      [...members].map(([path, bytes]) => [
        path,
        { bytes: bytes.length, sha256: hash(bytes) },
      ]),
    ),
  };
  const root = {
    releases: [
      { release_id: rid, manifest_path: `releases/${rid}/manifest.json` },
    ],
  };
  const rootBytes = canonical(root);
  const index = {
    format: FORMAT,
    root_manifest_sha256: hash(rootBytes),
    targets: { [target]: [] },
  };
  const counts = new Map();
  const write = async (path, bytes) => {
    await mkdir(dirname(join(directory, path)), { recursive: true });
    await writeFile(join(directory, path), bytes);
  };
  const setPart = async (i, entries, encodedBytes) => {
    const decoded =
      encodedBytes ?? canonical({ format: FORMAT, members: entries });
    const bytes = gzipSync(decoded);
    const path = `packs/${rid}/human/${key}-${i}.json.gz`;
    await write(path, bytes);
    index.targets[target][i] = {
      path,
      sha256: hash(bytes),
      bytes: bytes.length,
      decoded_bytes: decoded.length,
      members: entries.length,
      min_path: entries[0][0],
      max_path: entries.at(-1)[0],
    };
  };
  let i = 0;
  for (const [path, bytes] of [...members].sort(([a], [b]) => (a < b ? -1 : 1)))
    await setPart(i++, [[path, bytes.toString("base64")]]);
  await write("manifest.json", rootBytes);
  await write(
    "manifest.sha256",
    Buffer.from(`${hash(rootBytes)}  manifest.json\n`),
  );
  await write(root.releases[0].manifest_path, canonical(manifest));
  const readPhysical = async (path, max, expected, collect = true) => {
    counts.set(path, (counts.get(path) ?? 0) + 1);
    const bytes = await readFile(join(directory, path));
    assert.ok(bytes.length <= max, "physical size");
    if (expected) {
      assert.equal(bytes.length, expected.bytes);
      assert.equal(hash(bytes), expected.sha256);
    }
    return collect ? bytes : null;
  };
  const seal = async () => {
    const bytes = canonical(index);
    await write("transport-index.json", bytes);
    return {
      expectedRootSha256: hash(rootBytes),
      expectedTransportSha256: hash(bytes),
    };
  };
  const open = async () => {
    const reader = await createStorage(directory, await seal(), readPhysical);
    reader.registerManifest(manifest);
    return reader;
  };
  return {
    directory,
    manifest,
    root,
    index,
    target,
    construct,
    sequence,
    reference,
    members,
    counts,
    write,
    setPart,
    readPhysical,
    seal,
    open,
  };
}

test("packed source bytes and split groups preserve lexical numbers", async (t) => {
  const f = await fixture(t),
    reader = await f.open();
  await reader.withContext(async () => {
    const scanned = new Map();
    for await (const item of reader.scanRelease(f.manifest))
      scanned.set(item.path, item.bytes);
    assert.equal(scanned.size, f.members.size);
    for (const [path, expected] of f.members)
      assert.deepEqual(scanned.get(path), expected);
  });
  await reader.validateInventory(f.root);
  assert.deepEqual(
    await reader.read(
      f.sequence,
      1024,
      f.manifest.files[f.sequence],
      true,
      f.target,
    ),
    f.members.get(f.sequence),
  );
  await assert.rejects(reader.read(f.sequence, 1024), /context required/);
  await assert.rejects(reader.read(f.target, 1), /size/);
});

test("packed mode needs both independent pins and matching root sidecar", async (t) => {
  const f = await fixture(t),
    pins = await f.seal();
  await assert.rejects(
    createStorage(f.directory, {}, f.readPhysical),
    /explicit pins/,
  );
  await assert.rejects(
    createStorage(
      f.directory,
      { expectedRootSha256: pins.expectedRootSha256 },
      f.readPhysical,
    ),
    /both pins/,
  );
  await assert.rejects(
    createStorage(
      f.directory,
      { ...pins, expectedRootSha256: "d".repeat(64) },
      f.readPhysical,
    ),
    /root pin/,
  );
  await assert.rejects(
    createStorage(
      f.directory,
      { ...pins, expectedTransportSha256: "d".repeat(64) },
      f.readPhysical,
    ),
    /transport pin/,
  );
  await f.write("manifest.sha256", Buffer.from("incorrect\n"));
  await assert.rejects(
    createStorage(f.directory, pins, f.readPhysical),
    /sidecar/,
  );
});

test("unpacked mode delegates reads without a transport index", async (t) => {
  const f = await fixture(t);
  const reader = await createStorage(f.directory, {}, f.readPhysical);
  assert.equal(reader.packed, false);
  assert.deepEqual(
    await reader.read("manifest.json", 1024),
    await f.readPhysical("manifest.json", 1024),
  );
});

test("packed root reads remain pinned after storage creation", async (t) => {
  const f = await fixture(t),
    reader = await f.open();
  await f.write("manifest.json", Buffer.from("{}\n"));
  await assert.rejects(reader.read("manifest.json", 1024));
});

test("closed physical inventory rejects extras, directories, and symlinks", async (t) => {
  for (const kind of ["file", "directory", "symlink"]) {
    const f = await fixture(t),
      reader = await f.open();
    if (kind === "file")
      await f.write("private.sqlite", Buffer.from("private"));
    if (kind === "directory") await mkdir(join(f.directory, "unknown"));
    if (kind === "symlink")
      await symlink(
        join(f.directory, "manifest.json"),
        join(f.directory, "alias.json"),
      );
    await assert.rejects(reader.validateInventory(f.root), /unknown|symlink/);
  }
});

test("transport descriptor bounds and numbered contiguous ranges fail closed", async (t) => {
  for (const mutate of [
    (parts) => {
      parts[0].path = parts[0].path.replace("-0.json.gz", "-3.json.gz");
    },
    (parts) => {
      parts[0].decoded_bytes = 16 * 1024 * 1024 + 1;
    },
    (parts) => {
      parts[0].bytes = 17 * 1024 * 1024 + 1;
    },
    (parts) => {
      parts[1].min_path = parts[0].max_path;
    },
    (parts) => {
      parts[0].extra = true;
    },
  ]) {
    const f = await fixture(t);
    mutate(f.index.targets[f.target]);
    await assert.rejects(f.open(), /descriptor|range/);
  }
});

test("corrupt gzip and noncanonical base64 are rejected", async (t) => {
  const f = await fixture(t);
  const part = f.index.targets[f.target][0];
  const bad = Buffer.from("invalid gzip");
  await f.write(part.path, bad);
  part.bytes = bad.length;
  part.sha256 = hash(bad);
  let reader = await f.open();
  await assert.rejects(reader.read(f.construct, 1024));
  await f.setPart(0, [
    [f.construct, `${f.members.get(f.construct).toString("base64")}\n`],
  ]);
  reader = await f.open();
  await assert.rejects(reader.read(f.construct, 1024), /base64/);
});

test("source-unreachable, missing, and changed members are rejected", async (t) => {
  const f = await fixture(t);
  await f.setPart(0, [[f.construct, Buffer.from("wrong").toString("base64")]]);
  let reader = await f.open();
  await assert.rejects(reader.read(f.construct, 1024), /member bytes/);
  await f.setPart(0, [
    [f.construct, f.members.get(f.construct).toString("base64")],
  ]);
  f.manifest.files[
    `releases/${f.manifest.release_id}/human/sequences/${"d".repeat(64)}.json`
  ] = { bytes: 0, sha256: hash("") };
  reader = await f.open();
  await assert.rejects(async () => {
    for await (const ignored of reader.scanRelease(f.manifest)) void ignored;
  }, /missing manifest member/);
});

test("request cache holds at most two containers and is not global", async (t) => {
  const f = await fixture(t),
    reader = await f.open();
  const part = f.index.targets[f.target][0].path;
  f.counts.clear();
  await reader.withContext(async () => {
    await reader.read(f.construct, 1024);
    await reader.read(f.construct, 1024);
    assert.equal(f.counts.get(part), 1);
    await reader.read(f.reference, 1024);
    await reader.read(f.target, 1024);
    await reader.read(f.construct, 1024);
    assert.equal(f.counts.get(part), 2);
  });
  await reader.read(f.construct, 1024);
  await reader.read(f.construct, 1024);
  assert.equal(f.counts.get(part), 4);
});

test("source references control membership independently of signed pack metadata", async (t) => {
  const f = await fixture(t);
  const target = JSON.parse(f.members.get(f.target));
  target.data.constructs.pages = [];
  const bytes = canonical(target);
  f.manifest.files[f.target] = { bytes: bytes.length, sha256: hash(bytes) };
  await f.setPart(3, [[f.target, bytes.toString("base64")]]);
  const reader = await f.open();
  await assert.rejects(async () => {
    for await (const ignored of reader.scanRelease(f.manifest)) void ignored;
  }, /unreachable/);
});

test("decoded part count and range must agree with signed descriptors", async (t) => {
  for (const change of [
    (part) => {
      part.members = 2;
    },
    (part) => {
      part.max_path = part.max_path.replace("-0.json", "-1.json");
    },
  ]) {
    const f = await fixture(t);
    change(f.index.targets[f.target][0]);
    const reader = await f.open();
    await assert.rejects(reader.read(f.construct, 1024), /count|range/);
  }
});

test("gzip decoded output is bounded before JSON parsing", async (t) => {
  const f = await fixture(t),
    part = f.index.targets[f.target][0];
  const bytes = gzipSync(Buffer.alloc(16 * 1024 * 1024 + 1, 32));
  await f.write(part.path, bytes);
  part.bytes = bytes.length;
  part.sha256 = hash(bytes);
  const reader = await f.open();
  await assert.rejects(
    reader.read(f.construct, 1024),
    /larger than|buffer|size/i,
  );
});

test(
  "optional real Python transport fixture validates complete inventory",
  { skip: !process.env.OPENANTIGENS_PACKED_FIXTURE },
  async () => {
    const directory = process.env.OPENANTIGENS_PACKED_FIXTURE;
    const readPhysical = async (path, max, expected, collect = true) => {
      const bytes = await readFile(join(directory, path));
      assert.ok(bytes.length <= max);
      if (expected) {
        assert.equal(bytes.length, expected.bytes);
        assert.equal(hash(bytes), expected.sha256);
      }
      return collect ? bytes : null;
    };
    const root = JSON.parse(
      await readPhysical("manifest.json", 16 * 1024 * 1024),
    );
    const reader = await createStorage(
      directory,
      {
        expectedRootSha256: hash(
          await readPhysical("manifest.json", 16 * 1024 * 1024),
        ),
        expectedTransportSha256: hash(
          await readPhysical("transport-index.json", 16 * 1024 * 1024),
        ),
      },
      readPhysical,
    );
    for (const pointer of root.releases) {
      const manifest = JSON.parse(
        await readPhysical(pointer.manifest_path, 96 * 1024 * 1024),
      );
      reader.registerManifest(manifest);
      await reader.withContext(async () => {
        for await (const ignored of reader.scanRelease(manifest)) void ignored;
      });
    }
    await reader.validateInventory(root);
  },
);

import test from "node:test";
import assert from "node:assert/strict";
import { readFile, writeFile, rm, cp, mkdtemp } from "node:fs/promises";
import { join } from "node:path";
import { tmpdir } from "node:os";
import {
  Client,
  StreamableHTTPClientTransport,
} from "@modelcontextprotocol/client";
import { Provider } from "../src/provider.js";
import { makeServer, readBuildPins } from "../src/server.js";
import { fixture } from "./fixture.js";
import { queryFixture } from "./query-fixture.js";
const admission = {
  clientRate: 1000,
  clientBurst: 1000,
  globalRate: 10000,
  globalBurst: 10000,
  concurrency: 8,
  clientConcurrency: 8,
};
function origin(directory) {
  const calls = [];
  const control = { bad: false };
  const fetchImpl = async (url, options) => {
    calls.push(url);
    const parsed = new URL(url);
    assert.equal(parsed.origin, "https://openantigens.org");
    assert.match(
      parsed.pathname,
      /^\/api-runtime\/[a-f0-9]{64}\/(human|mouse)\/[a-f0-9]{64}-[01]-[a-f0-9]{64}\.oap$/,
    );
    assert.equal(options.redirect, "error");
    assert.equal(options.credentials, "omit");
    assert.equal(options.headers["accept-encoding"], "identity");
    if (control.bad)
      return new Response("challenge", {
        headers: { "content-type": "text/html" },
      });
    const bytes = await readFile(join(directory, parsed.pathname));
    return new Response(bytes, {
      headers: {
        "content-type": "application/octet-stream",
        "content-length": String(bytes.length),
      },
    });
  };
  return { calls, control, fetchImpl };
}
async function own(t, options = {}) {
  const source = await fixture("query-integration", options);
  const query = await queryFixture(source.directory);
  t.after(async () => {
    await rm(source.directory, { recursive: true, force: true });
    await rm(query.directory, { recursive: true, force: true });
    await rm(query.staticDirectory, { recursive: true, force: true });
  });
  return { source, query };
}
async function serve(t, provider) {
  const server = makeServer(provider, admission);
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  t.after(
    () =>
      new Promise((resolve) => {
        server.closeAllConnections();
        server.close(resolve);
      }),
  );
  return `http://127.0.0.1:${server.address().port}`;
}
test("query mode boots and searches with zero origin work and preserves unresolved source-offset association", async (t) => {
  const { source, query } = await own(t, {
    unresolvedFirst: true,
    extraUnresolved: 2,
  });
  const remote = origin(query.staticDirectory);
  const provider = await new Provider(query.directory, {
    querySeal: query.seal,
    fetchImpl: remote.fetchImpl,
  }).initialize();
  assert.equal(remote.calls.length, 0);
  const args = { release_id: provider.current, species: "human", q: "catalog" };
  const page = await provider.call("search_targets", args);
  assert.equal(page.total, 4);
  assert.equal(page.data[0].target_id, source.target_id);
  assert.equal(page.data.filter((row) => row.target_id === null).length, 3);
  assert.equal(remote.calls.length, 0);
  assert.ok(
    page.data.every(
      (row) =>
        Object.hasOwn(row, "protein_name") &&
        !Object.hasOwn(row, "summary_path"),
    ),
  );
  await provider.call("get_release", {});
  assert.equal(remote.calls.length, 0);
  assert.deepEqual(provider.storage.snapshot(), {
    activeTargets: 0,
    activeFetches: 0,
    cacheBytes: 0,
    cacheEntries: 0,
  });
});
test("remote target failure is request-scoped, repairs recover, and known-good compressed cache survives", async (t) => {
  const { source, query } = await own(t);
  const remote = origin(query.staticDirectory);
  remote.control.bad = true;
  const provider = await new Provider(query.directory, {
    querySeal: query.seal,
    fetchImpl: remote.fetchImpl,
  }).initialize();
  const args = {
    release_id: provider.current,
    species: "human",
    target_id: source.target_id,
  };
  await assert.rejects(
    provider.call("get_target", args),
    (error) => error.status === 503 && error.code === "unavailable",
  );
  assert.equal(provider.ready, true);
  assert.equal(provider.storage.snapshot().cacheEntries, 0);
  await provider.call("search_targets", {
    release_id: provider.current,
    species: "human",
    q: "catalog",
  });
  await provider.call("get_release", {});
  assert.equal(remote.calls.length, 1);
  remote.control.bad = false;
  const repaired = await provider.call("get_target", args);
  assert.equal(repaired.data.target_id, source.target_id);
  assert.equal(provider.storage.snapshot().cacheEntries, 1);
  remote.control.bad = true;
  assert.deepEqual(await provider.call("get_target", args), repaired);
  assert.equal(remote.calls.length, 2);
  assert.equal(provider.ready, true);
  const base = await serve(t, provider);
  assert.equal((await fetch(base + "/readyz")).status, 200);
});
test("local catalog, summary bytes and code seals are independently pinned without remote fallback", async (t) => {
  const { query } = await own(t);
  const remote = origin(query.staticDirectory);
  await assert.rejects(
    new Provider(query.directory, {
      querySeal: { ...query.seal, query_root_sha256: "f".repeat(64) },
      fetchImpl: remote.fetchImpl,
    }).initialize(),
  );
  assert.equal(remote.calls.length, 0);
  await assert.rejects(new Provider(query.directory).initialize());
  const provider = await new Provider(query.directory, {
    querySeal: query.seal,
    fetchImpl: remote.fetchImpl,
  }).initialize();
  const row = provider.resolve(provider.current).indexes.human[0];
  const file = provider.resolve(provider.current).manifest.species.human
    .summary_lines;
  const bytes = await readFile(join(query.directory, file.path));
  bytes[10] = bytes[10] === 32 ? 33 : 32;
  await writeFile(join(query.directory, file.path), bytes);
  await assert.rejects(
    provider.call("search_targets", {
      release_id: provider.current,
      species: "human",
      q: "catalog",
    }),
    (error) => error.status === 503 && error.localQueryIntegrity === true,
  );
  assert.equal(provider.ready, false);
  assert.equal(remote.calls.length, 0);
  const sealPath = join(query.directory, "app-seal.json");
  await writeFile(sealPath, JSON.stringify({ data_transport: query.seal }));
  assert.deepEqual(await readBuildPins(sealPath), { querySeal: query.seal });
  assert.ok(row.record_id);
});
test(
  "actual Python query runtime has cold/warm HTTP and official MCP parity for all five tools",
  { skip: !process.env.OPENANTIGENS_QUERY_FIXTURE },
  async (t) => {
    const root = process.env.OPENANTIGENS_QUERY_FIXTURE;
    const seal = JSON.parse(await readFile(join(root, "remote-seal.json")));
    const remote = origin(join(root, "remote-static"));
    const provider = await new Provider(join(root, "remote-runtime"), {
      querySeal: seal,
      fetchImpl: remote.fetchImpl,
    }).initialize();
    const plain = await new Provider(join(root, "bundle")).initialize();
    assert.equal(remote.calls.length, 0);
    const base = await serve(t, provider);
    const client = new Client(
      { name: "own-static-fixture", version: "1" },
      { versionNegotiation: { mode: "auto" } },
    );
    await client.connect(
      new StreamableHTTPClientTransport(new URL(base + "/mcp")),
    );
    t.after(() => client.close());
    for (const name of ["get_release", "search_targets"]) {
      const args =
        name === "get_release"
          ? {}
          : {
              release_id: provider.current,
              species: "mouse",
              q: "NP_000001.1",
            };
      const result = await client.callTool({ name, arguments: args });
      assert.notEqual(result.isError, true);
      const expected = await plain.call(name, args);
      assert.deepEqual(result.structuredContent, expected);
      assert.equal(
        JSON.stringify(result.structuredContent),
        JSON.stringify(expected),
      );
    }
    assert.equal(remote.calls.length, 0);
    for (const species of ["human", "mouse"]) {
      const target_id =
        species === "human" ? "9606:uniprot:P00001" : "10090:uniprot:P00002";
      const args = { release_id: provider.current, species, target_id };
      for (const name of ["get_target", "list_constructs", "get_sequence"]) {
        const mcp = await client.callTool({ name, arguments: args });
        assert.notEqual(mcp.isError, true);
        const expected = await plain.call(name, args);
        assert.deepEqual(mcp.structuredContent, expected);
        assert.equal(
          JSON.stringify(mcp.structuredContent),
          JSON.stringify(expected),
        );
      }
      const before = remote.calls.length;
      const again = await client.callTool({
        name: "get_target",
        arguments: args,
      });
      assert.notEqual(again.isError, true);
      assert.equal(
        remote.calls.length,
        before,
        "Warm reads use verified compressed cache",
      );
    }
    assert.ok(remote.calls.length <= 4);
    assert.equal(provider.ready, true);
  },
);
test(
  "actual Python warning-heavy query runtime retains static overflow recovery",
  { skip: !process.env.OPENANTIGENS_QUERY_OVERFLOW_FIXTURE },
  async (t) => {
    const root = process.env.OPENANTIGENS_QUERY_OVERFLOW_FIXTURE;
    const seal = JSON.parse(await readFile(join(root, "remote-seal.json")));
    const remote = origin(join(root, "remote-static"));
    const provider = await new Provider(join(root, "remote-runtime"), {
      querySeal: seal,
      fetchImpl: remote.fetchImpl,
    }).initialize();
    const base = await serve(t, provider);
    const target_id = "9606:uniprot:P00001";
    const args = { release_id: provider.current, species: "human", target_id };
    const constructs = await provider.call("list_constructs", args);
    const construct = constructs.data.find((row) =>
      row.warnings.some((warning) => warning.length > 1000),
    );
    assert.ok(construct);
    const response = await fetch(
      `${base}/v1/releases/${provider.current}/human/targets/${target_id}/constructs/${construct.construct_id}/sequence?format=fasta`,
    );
    assert.equal(response.status, 413);
    const error = (await response.json()).error;
    assert.deepEqual(
      error.download,
      provider.fallback(provider.current, "human"),
    );
    assert.deepEqual(error.record_selector, {
      target_id,
      construct_id: construct.construct_id,
    });
    assert.equal(response.headers.get("content-disposition"), null);
    assert.equal(provider.ready, true);
  },
);

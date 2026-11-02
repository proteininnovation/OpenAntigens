import test from "node:test";
import assert from "node:assert/strict";
import { rm } from "node:fs/promises";
import { Provider } from "../src/provider.js";
import { makeServer } from "../src/server.js";
import { fixture } from "./fixture.js";

const origins = [undefined, "https://openantigens.org", "https://www.openantigens.org"];
const varyTokens = (response) => response.headers.get("vary").split(",").map((value) => value.trim().toLowerCase());

async function setup(t, priorVary) {
  const f = await fixture();
  const provider = await new Provider(f.directory).initialize();
  const server = makeServer(provider, {
    clientRate: 1000,
    clientBurst: 1000,
    globalRate: 1000,
    globalBurst: 1000,
  });
  if (priorVary !== undefined)
    server.prependListener("request", (_req, res) => res.setHeader("Vary", priorVary));
  await new Promise((resolve) => server.listen(0, "127.0.0.1", resolve));
  t.after(async () => {
    server.closeAllConnections();
    await new Promise((resolve) => server.close(resolve));
    await rm(f.directory, { recursive: true, force: true });
  });
  const base = `http://127.0.0.1:${server.address().port}`;
  const target = (species) => `${base}/v1/releases/${f.release_id}/${species}/targets/${species === "human" ? "9606:uniprot:P12345" : "10090:uniprot:Q12345"}`;
  return { target };
}

function assertOrigin(response, origin) {
  assert.equal(response.headers.get("access-control-allow-origin"), origin ?? null);
  assert.ok(varyTokens(response).includes("origin"));
  assert.equal(response.headers.get("cache-control"), "public, max-age=31536000, immutable");
}

test("both species vary immutable GET, HEAD and cached 304 responses by Origin", async (t) => {
  const f = await setup(t);
  for (const species of ["human", "mouse"]) {
    let body;
    let etag;
    for (const origin of origins) {
      const headers = origin ? { origin } : {};
      const response = await fetch(f.target(species), { headers });
      assert.equal(response.status, 200);
      assertOrigin(response, origin);
      const representation = await response.text();
      body ??= representation;
      etag ??= response.headers.get("etag");
      assert.ok(etag);
      assert.equal(representation, body);
      assert.equal(response.headers.get("etag"), etag);
      const head = await fetch(f.target(species), { method: "HEAD", headers });
      assert.equal(head.status, 200);
      assertOrigin(head, origin);
      assert.equal(head.headers.get("etag"), etag);
      assert.equal(await head.text(), "");
      const cached = await fetch(f.target(species), {
        headers: { ...headers, "if-none-match": etag },
      });
      assert.equal(cached.status, 304);
      assertOrigin(cached, origin);
      assert.equal(cached.headers.get("etag"), etag);
      assert.equal(await cached.text(), "");
    }
  }
});

test("CORS preserves existing Vary tokens, arrays and wildcard semantics", async (t) => {
  for (const [prior, expected] of [
    ["Accept-Encoding", ["accept-encoding", "origin"]],
    [["Accept-Encoding", "Accept-Language"], ["accept-encoding", "accept-language", "origin"]],
    ["Accept-Encoding, oRiGiN", ["accept-encoding", "origin"]],
    ["*", ["*"]],
    [["Accept-Encoding", "*"], ["accept-encoding", "*"]],
  ]) {
    const f = await setup(t, prior);
    for (const origin of origins) {
      const response = await fetch(f.target("human"), { headers: origin ? { origin } : {} });
      assert.equal(response.status, 200);
      assert.equal(response.headers.get("access-control-allow-origin"), origin ?? null);
      assert.deepEqual(varyTokens(response), expected);
      await response.text();
    }
  }
});

test("denied origins retain no-store errors and do not receive allow-origin", async (t) => {
  const f = await setup(t, "Accept-Encoding");
  const response = await fetch(f.target("human"), { headers: { origin: "https://untrusted.example" } });
  assert.equal(response.status, 403);
  assert.equal(response.headers.get("access-control-allow-origin"), null);
  assert.equal(response.headers.get("cache-control"), "no-store");
  assert.deepEqual(varyTokens(response), ["accept-encoding", "origin"]);
  assert.equal((await response.json()).error.code, "invalid_origin");
});

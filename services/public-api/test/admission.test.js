import test from "node:test";
import assert from "node:assert/strict";
import { request } from "node:http";
import { rm } from "node:fs/promises";
import { Provider } from "../src/provider.js";
import { makeServer } from "../src/server.js";
import { fixture } from "./fixture.js";

async function setup(t, options = {}) {
  const f = await fixture();
  const provider = await new Provider(f.directory).initialize();
  let now = 0;
  const server = makeServer(provider, { now: () => now, ...options });
  await new Promise((resolve) => server.listen(0, "::", resolve));
  t.after(async () => {
    server.closeAllConnections();
    await new Promise((resolve) => server.close(resolve));
    await rm(f.directory, { recursive: true, force: true });
  });
  const read = (path = "/healthz", headers = {}, hostname = "127.0.0.1") =>
    new Promise((resolve, reject) => {
      const req = request(
        {
          hostname,
          port: server.address().port,
          path,
          headers: { host: "localhost", ...headers },
          agent: false,
        },
        (res) => {
          let body = "";
          res.setEncoding("utf8");
          res.on("data", (chunk) => {
            body += chunk;
          });
          res.on("end", () =>
            resolve({ status: res.statusCode, headers: res.headers, body }),
          );
        },
      );
      req.on("error", reject);
      req.end();
    });
  return {
    ...f,
    provider,
    read,
    advance: (milliseconds) => {
      now += milliseconds;
    },
    target: `/v1/releases/${f.release_id}/human/targets/${f.target_id}`,
  };
}

function holdCalls(provider) {
  const call = provider.call.bind(provider);
  let release;
  const held = new Promise((resolve) => {
    release = resolve;
  });
  let started;
  const starts = [];
  provider.call = async (...args) => {
    starts.push(args);
    started?.();
    await held;
    return call(...args);
  };
  return {
    release,
    started: (count) =>
      starts.length >= count
        ? Promise.resolve()
        : new Promise((resolve) => {
            started = () => {
              if (starts.length >= count) resolve();
            };
          }),
  };
}

test("default client bucket allows eight requests and refills at twenty per minute", async (t) => {
  const f = await setup(t);
  for (let n = 0; n < 8; n++) assert.equal((await f.read()).status, 200);
  let response = await f.read();
  assert.equal(response.status, 429);
  assert.equal(response.headers["retry-after"], "3");
  assert.equal(response.headers["cache-control"], "no-store");
  f.advance(2000);
  response = await f.read();
  assert.equal(response.status, 429);
  assert.equal(response.headers["retry-after"], "1");
  f.advance(1000);
  assert.equal((await f.read()).status, 200);
  assert.equal((await f.read()).status, 429);
  f.advance(60000);
  for (let n = 0; n < 8; n++) assert.equal((await f.read()).status, 200);
  assert.equal((await f.read()).status, 429);
});

test("forwarded addresses cannot evade the socket client bucket", async (t) => {
  const f = await setup(t, { clientBurst: 1 });
  assert.equal(
    (await f.read("/healthz", { "x-forwarded-for": "192.0.2.1" })).status,
    200,
  );
  const response = await f.read("/healthz", {
    "x-forwarded-for": "192.0.2.2",
    forwarded: "for=192.0.2.3",
    "x-real-ip": "192.0.2.4",
  });
  assert.equal(response.status, 429);
  assert.equal(response.headers["retry-after"], "3");
  assert.equal((await f.read("/healthz", {}, "::1")).status, 200);
});

test("default global bucket caps the burst at ten and refills at ninety per minute", async (t) => {
  const f = await setup(t, { clientBurst: 100 });
  for (let n = 0; n < 10; n++) assert.equal((await f.read()).status, 200);
  const response = await f.read("/healthz", {}, "::1");
  assert.equal(response.status, 429);
  assert.equal(response.headers["retry-after"], "1");
  f.advance(667);
  assert.equal((await f.read("/healthz", {}, "::1")).status, 200);
  assert.equal((await f.read()).status, 429);
});

test("one active request per socket client is enforced with no waiting queue", async (t) => {
  const f = await setup(t);
  const held = holdCalls(f.provider);
  t.after(held.release);
  const first = f.read(f.target);
  await held.started(1);
  const rejected = await f.read(f.target, { "x-forwarded-for": "192.0.2.1" });
  assert.equal(rejected.status, 429);
  assert.equal(rejected.headers["retry-after"], "1");
  const other = f.read(f.target, {}, "::1");
  await held.started(2);
  held.release();
  assert.equal((await first).status, 200);
  assert.equal((await other).status, 200);
  assert.equal((await f.read(f.target)).status, 200);
});

test("two active requests globally are enforced and busy rejections preserve tokens", async (t) => {
  const f = await setup(t, {
    clientBurst: 2,
    globalBurst: 3,
    clientConcurrency: 3,
  });
  const held = holdCalls(f.provider);
  t.after(held.release);
  const first = f.read(f.target);
  await held.started(1);
  const second = f.read(f.target, {}, "::1");
  await held.started(2);
  for (let n = 0; n < 4; n++) {
    const response = await f.read(f.target);
    assert.equal(response.status, 429);
    assert.equal(response.headers["retry-after"], "1");
  }
  held.release();
  assert.equal((await first).status, 200);
  assert.equal((await second).status, 200);
  assert.equal((await f.read(f.target)).status, 200);
  assert.equal((await f.read(f.target, {}, "::1")).status, 429);
});

test("interactive API has no species downloads route", async (t) => {
  const f = await setup(t);
  for (const species of ["human", "mouse"]) {
    const response = await f.read(
      `/v1/releases/${f.release_id}/${species}/downloads`,
    );
    assert.equal(response.status, 404);
    assert.equal(JSON.parse(response.body).error.code, "not_found");
  }
});

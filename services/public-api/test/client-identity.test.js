import test from "node:test";
import assert from "node:assert/strict";
import { makeClientResolver } from "../src/client-identity.js";

const req = (peer = "127.0.0.1", value, rawHeaders) => ({
  socket: { remoteAddress: peer },
  headers: value === undefined ? {} : { "x-forwarded-for": value },
  rawHeaders:
    rawHeaders ?? (value === undefined ? [] : ["X-Forwarded-For", value]),
});
const proxy = (options = {}) =>
  makeClientResolver({
    trustedProxyPeers: ["127.0.0.1"],
    proxyHopCount: 1,
    identityVerified: true,
    ...options,
  });
const peerResult = (ip, trusted = false) => ({
  rateKey: ip,
  ip,
  trusted,
  source: "peer",
});

test("direct mode ignores forged forwarded addresses and requires explicit identity verification", () => {
  for (const identityVerified of [false, true]) {
    const resolve = makeClientResolver({ identityVerified });
    assert.deepEqual(
      resolve(req("127.0.0.1", "192.0.2.1")),
      peerResult("127.0.0.1", identityVerified),
    );
  }
  assert.deepEqual(makeClientResolver()(req()), peerResult("127.0.0.1"));
});

test("configured forwarding requires both a verified identity policy and an exact proxy peer", () => {
  assert.deepEqual(
    proxy({ identityVerified: false })(req("127.0.0.1", "192.0.2.1")),
    peerResult("127.0.0.1"),
  );
  assert.deepEqual(
    proxy()(req("127.0.0.2", "192.0.2.1")),
    peerResult("127.0.0.2"),
  );
  assert.deepEqual(proxy()(req("127.0.0.1")), peerResult("127.0.0.1"));
});

test("verified forwarding selects from the right at the configured hop count", () => {
  assert.deepEqual(proxy()(req("127.0.0.1", "192.0.2.9, 198.51.100.7")), {
    rateKey: "198.51.100.7",
    ip: "198.51.100.7",
    trusted: true,
    source: "forwarded",
  });
  const resolve = proxy({ proxyHopCount: 2 });
  assert.deepEqual(
    resolve(req("127.0.0.1", "192.0.2.9, 198.51.100.7, 203.0.113.8")),
    {
      rateKey: "198.51.100.7",
      ip: "198.51.100.7",
      trusted: true,
      source: "forwarded",
    },
  );
  assert.deepEqual(
    resolve(req("127.0.0.1", "198.51.100.7")),
    peerResult("127.0.0.1"),
  );
});

test("duplicate or malformed forwarded headers fall back to an untrusted socket peer", () => {
  const resolve = proxy();
  const duplicates = [
    "X-Forwarded-For",
    "192.0.2.1",
    "x-forwarded-for",
    "198.51.100.1",
  ];
  assert.deepEqual(
    resolve(req("127.0.0.1", "192.0.2.1, 198.51.100.1", duplicates)),
    peerResult("127.0.0.1"),
  );
  for (const value of [
    "",
    "unknown",
    "192.0.2.1,",
    ",192.0.2.1",
    "spoofed, 192.0.2.1",
    "192.0.2.1:8080",
    "[2001:db8::1]",
    "[2001:db8::1]:8080",
    "fe80::1%en0",
    "192.168.001.1",
    "192.0.2.1\n",
    "\u00a0192.0.2.1",
  ])
    assert.deepEqual(
      resolve(req("127.0.0.1", value)),
      peerResult("127.0.0.1"),
      value,
    );
  assert.deepEqual(
    resolve(req("127.0.0.1", "192.0.2.1", [])),
    peerResult("127.0.0.1"),
  );
});

test("IPv6 and IPv4-mapped addresses share canonical peer and forwarded identities", () => {
  const direct = makeClientResolver({ identityVerified: true });
  for (const [input, output] of [
    ["2001:0DB8:0000:0000:0000:0000:0000:0001", "2001:db8::1"],
    ["::FFFF:192.0.2.1", "192.0.2.1"],
    ["0000:0000:0000:0000:0000:FFFF:C000:0201", "192.0.2.1"],
    ["::ffff:0:1", "0.0.0.1"],
    ["::1", "::1"],
  ])
    assert.deepEqual(direct(req(input)), peerResult(output, true));
  const resolve = proxy({ trustedProxyPeers: ["::ffff:7f00:1"] });
  assert.deepEqual(
    resolve(req("::ffff:127.0.0.1", "192.0.2.1,\t2001:0DB8::A ")),
    {
      rateKey: "2001:db8::a",
      ip: "2001:db8::a",
      trusted: true,
      source: "forwarded",
    },
  );
  assert.deepEqual(resolve(req("127.0.0.1", "::ffff:c000:201")), {
    rateKey: "192.0.2.1",
    ip: "192.0.2.1",
    trusted: true,
    source: "forwarded",
  });
});

test("missing or invalid socket identities cannot establish a forwarded client", () => {
  for (const peer of [undefined, "unknown", "fe80::1%en0", "[::1]"]) {
    const request = {
      socket: { remoteAddress: peer },
      rawHeaders: ["X-Forwarded-For", "192.0.2.1"],
    };
    assert.deepEqual(proxy()(request), {
      rateKey: "unavailable",
      ip: null,
      trusted: false,
      source: "unavailable",
    });
  }
});

test("invalid proxy configurations use one fixed error without input details", () => {
  for (const options of [
    null,
    [],
    "192.0.2.1",
    { trustedProxyPeers: "127.0.0.1" },
    { trustedProxyPeers: Array(1), proxyHopCount: 1 },
    { trustedProxyPeers: ["127.0.0.1"], proxyHopCount: 0 },
    { trustedProxyPeers: [], proxyHopCount: 1 },
    { trustedProxyPeers: ["bad-secret-value"], proxyHopCount: 1 },
    { trustedProxyPeers: ["127.0.0.1/8"], proxyHopCount: 1 },
    { trustedProxyPeers: ["fe80::1%en0"], proxyHopCount: 1 },
    { trustedProxyPeers: ["127.0.0.1"], proxyHopCount: 9 },
    { trustedProxyPeers: ["127.0.0.1"], proxyHopCount: 1.5 },
    { proxyHopCount: -1 },
    { identityVerified: "true" },
  ])
    assert.throws(() => makeClientResolver(options), {
      message: "Invalid client identity configuration.",
    });
  assert.doesNotThrow(() => proxy({ proxyHopCount: 8 }));
});

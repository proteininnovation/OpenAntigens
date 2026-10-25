import { isIP } from "node:net";

function normalizeIP(value) {
  if (typeof value !== "string" || value.includes("%")) return null;
  const version = isIP(value);
  if (version === 4) return value;
  if (version !== 6) return null;
  const normalized = new URL(`http://[${value}]/`).hostname.slice(1, -1);
  const mapped = /^::ffff:([a-f0-9]{1,4}):([a-f0-9]{1,4})$/.exec(normalized);
  if (!mapped) return normalized;
  const high = parseInt(mapped[1], 16);
  const low = parseInt(mapped[2], 16);
  return `${high >> 8}.${high & 255}.${low >> 8}.${low & 255}`;
}

export function makeClientResolver(options = {}) {
  const invalid = () => {
    throw new Error("Invalid client identity configuration.");
  };
  if (!options || typeof options !== "object" || Array.isArray(options))
    invalid();
  const {
    trustedProxyPeers = [],
    proxyHopCount = 0,
    identityVerified = false,
  } = options;
  if (
    !Array.isArray(trustedProxyPeers) ||
    typeof identityVerified !== "boolean" ||
    !Number.isInteger(proxyHopCount) ||
    proxyHopCount < 0 ||
    proxyHopCount > 8 ||
    (trustedProxyPeers.length === 0) !== (proxyHopCount === 0)
  )
    invalid();
  const peers = Array.from(trustedProxyPeers, normalizeIP);
  if (peers.some((peer) => peer === null)) invalid();
  const allowedPeers = new Set(peers);

  return (req) => {
    const peer = normalizeIP(req?.socket?.remoteAddress);
    const fallback = {
      rateKey: peer ?? "unavailable",
      ip: peer,
      trusted: Boolean(peer && identityVerified && allowedPeers.size === 0),
      source: peer ? "peer" : "unavailable",
    };
    if (!peer || !identityVerified || !allowedPeers.has(peer)) return fallback;
    const forwarded = [];
    const headers = req.rawHeaders;
    if (!Array.isArray(headers)) return fallback;
    for (let i = 0; i < headers.length; i += 2) {
      if (
        typeof headers[i] === "string" &&
        headers[i].toLowerCase() === "x-forwarded-for"
      )
        forwarded.push(headers[i + 1]);
    }
    if (forwarded.length !== 1 || typeof forwarded[0] !== "string")
      return fallback;
    const chain = forwarded[0]
      .split(",")
      .map((value) => normalizeIP(value.replace(/^[ \t]+|[ \t]+$/g, "")));
    if (chain.length < proxyHopCount || chain.some((ip) => ip === null))
      return fallback;
    const ip = chain[chain.length - proxyHopCount];
    return { rateKey: ip, ip, trusted: true, source: "forwarded" };
  };
}

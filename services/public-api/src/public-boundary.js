import { readFileSync } from "node:fs";
import { parsePublicJson } from "./json.js";
const policy = parsePublicJson(
  readFileSync(new URL("./query-schema.json", import.meta.url)),
)["x-public-text-policy"];
if (
  !policy ||
  policy.max_length !== 1024 ||
  policy.percent_decode_passes !== 4 ||
  policy.scheme !== "https" ||
  policy.reject_unicode_whitespace !== true ||
  policy.allow_query ||
  policy.allow_fragment ||
  !Array.isArray(policy.private_markers) ||
  !Array.isArray(policy.blocked_host_suffixes)
)
  throw new Error("Invalid public text policy.");
const dns = new RegExp(`^(?:${policy.dns_pattern})$`);
const markers = policy.private_markers.map((marker) => marker.toLowerCase());
const reject = () => {
  throw new Error("Invalid public retirement metadata.");
};
function scanPublic(value, key = false) {
  if (typeof value === "string") {
    if (
      value.length > (key ? 16384 : 1000000) ||
      policy.private_markers.some((marker) => value.includes(marker)) ||
      /[\u0000-\u0008\u000b\u000c\u000e-\u001f]/u.test(value)
    )
      reject();
  } else if (Array.isArray(value)) for (const item of value) scanPublic(item);
  else if (value && typeof value === "object")
    for (const [name, item] of Object.entries(value)) {
      scanPublic(name, true);
      scanPublic(item);
    }
}
export function isPublicArchiveUrl(value) {
  if (value === null) return true;
  if (
    typeof value !== "string" ||
    !value ||
    value.length > policy.max_length ||
    /%(?![a-f0-9]{2})/iu.test(value)
  )
    return false;
  let text = value;
  try {
    for (let layer = 0; layer <= policy.percent_decode_passes; layer++) {
      scanPublic(text);
      if (
        !text.startsWith(`${policy.scheme}://`) ||
        /[\p{White_Space}\u0000-\u0020\u007f]/u.test(text) ||
        text.includes("\\") ||
        text.includes("?") ||
        text.includes("#") ||
        markers.some((marker) => text.toLowerCase().includes(marker))
      )
        return false;
      const authority = text.slice(8).split("/")[0];
      if (
        !authority ||
        /[^\u0000-\u007f]/u.test(authority) ||
        authority.includes("@")
      )
        return false;
      const rawHost = authority.split(":")[0].toLowerCase();
      const parsed = new URL(text);
      const host = parsed.hostname.toLowerCase();
      if (
        !dns.test(rawHost) ||
        rawHost !== host ||
        parsed.username ||
        parsed.password ||
        (parsed.port && !policy.allowed_ports.includes(Number(parsed.port))) ||
        policy.blocked_host_suffixes.some(
          (suffix) => host === suffix || host.endsWith("." + suffix),
        )
      )
        return false;
      const rawPath = text.slice(8 + authority.length);
      if (
        rawPath
          .split("/")
          .some((segment) => segment === "." || segment === "..")
      )
        return false;
      // urllib-compatible decoding preserves a literal percent introduced by a prior layer.
      const decoded = decodeURIComponent(
        text.replace(/%(?![a-f0-9]{2})/giu, "%25"),
      );
      if (decoded === text) return true;
      if (layer === policy.percent_decode_passes) return false;
      text = decoded;
    }
  } catch {
    return false;
  }
  return false;
}
export function validateRetirementMetadata(root) {
  if (
    !root ||
    typeof root !== "object" ||
    Array.isArray(root) ||
    !Array.isArray(root.retired_releases ?? [])
  )
    reject();
  scanPublic(root);
  for (const row of root.retired_releases ?? [])
    if (
      !row ||
      typeof row !== "object" ||
      Array.isArray(row) ||
      !isPublicArchiveUrl(row.archive_url)
    )
      reject();
}

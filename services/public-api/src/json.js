// Public JSON keeps its original numeric spelling; only object-key ambiguity is rejected.
const invalid = () => {
  throw new Error("Invalid public JSON.");
};
const whitespace = (code) =>
  code === 32 || code === 9 || code === 10 || code === 13;
function stringEnd(text, start) {
  for (let position = start + 1; position < text.length; position++) {
    if (text[position] === "\\") position++;
    else if (text[position] === '"') return position + 1;
  }
  invalid();
}
function keyAt(text, start, end) {
  const key = text.slice(start + 1, end - 1);
  return key.includes("\\") ? JSON.parse(text.slice(start, end)) : key;
}
function uniqueMembers(text) {
  const objects = [];
  let depth = 0;
  for (let position = 0; position < text.length; ) {
    const character = text[position];
    if (character === '"') {
      const start = position,
        end = stringEnd(text, start);
      position = end;
      let next = end;
      while (whitespace(text.charCodeAt(next))) next++;
      if (text[next] !== ":" || !objects.length) continue;
      const key = keyAt(text, start, end),
        frame = objects.at(-1);
      if (!frame.seen && frame.last !== undefined && key <= frame.last) {
        if (key === frame.last) invalid();
        frame.seen = new Set(
          frame.offsets.map((offset) =>
            keyAt(text, offset, stringEnd(text, offset)),
          ),
        );
        frame.offsets = null;
      }
      if (frame.seen) {
        if (frame.seen.has(key)) invalid();
        frame.seen.add(key);
      } else frame.offsets.push(start);
      frame.last = key;
      continue;
    }
    if (character === "{" || character === "[") {
      if (++depth > 64) invalid();
      if (character === "{")
        objects.push({ last: undefined, offsets: [], seen: null });
    } else if (character === "}" || character === "]") {
      if (--depth < 0) invalid();
      if (character === "}") objects.pop();
    }
    position++;
  }
}
export function parsePublicJson(input) {
  try {
    const text =
      typeof input === "string"
        ? input
        : new TextDecoder("utf-8", { fatal: true }).decode(input);
    uniqueMembers(text);
    return JSON.parse(text);
  } catch {
    invalid();
  }
}

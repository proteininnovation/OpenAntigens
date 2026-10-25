import test from "node:test";
import assert from "node:assert/strict";
import { parsePublicJson } from "../src/json.js";

const CANARY = "PRIVATE_CANARY";
const reject = (source) =>
  assert.throws(
    () => parsePublicJson(source),
    (error) => {
      assert.ok(error instanceof Error);
      assert.ok(
        !error.message.includes(CANARY),
        "parser errors must not expose input values",
      );
      return true;
    },
  );

test("duplicate sequence data cannot hide an earlier private payload", () => {
  reject(
    Buffer.from(`{"data":{"sequence":"${CANARY}"},"data":{"sequence":"AAA"}}`),
  );
  reject(`{"data":{"sequence":"AAA"},"data":{"sequence":"${CANARY}"}}`);
});

test("nested evidence object keys are unique at every level", () => {
  reject(`{"data":{"evidence":[{"source":"${CANARY}","source":"public"}]}}`);
  reject(
    `{"data":{"evidence":{"source":{"url":"${CANARY}","url":"https://example.org"}}}}`,
  );
  reject(
    `{"data":{"evidence":{"source":"public"},"evidence":{"source":"${CANARY}"}}}`,
  );
});

test("escaped, unsorted, and nonadjacent duplicate keys are rejected", () => {
  for (const source of [
    String.raw`{"data":"PRIVATE_CANARY","\u0064ata":"public"}`,
    String.raw`{"\u0064ata":"PRIVATE_CANARY","data":"public"}`,
    String.raw`{"z":"PRIVATE_CANARY","a":{"z":"nested"},"m":[1,2],"z":"public"}`,
    String.raw`{"a/b":"PRIVATE_CANARY","a\/b":"public"}`,
    String.raw`{"a\n":"PRIVATE_CANARY","a\u000a":"public"}`,
    String.raw`{"😀":"PRIVATE_CANARY","\ud83d\ude00":"public"}`,
    String.raw`{"__proto__":"PRIVATE_CANARY","__proto__":"public"}`,
    String.raw`{"constructor":"PRIVATE_CANARY","constructor":"public"}`,
  ])
    reject(source);
});

test("duplicates use fixed errors independent of the input key and value", () => {
  const messages = [];
  for (const source of [
    '{"a":0,"a":1}',
    `{"${CANARY}":"${CANARY}","${CANARY}":0}`,
    '{"x":{"z":1,"z":2}}',
  ]) {
    assert.throws(
      () => parsePublicJson(source),
      (error) => {
        assert.ok(!error.message.includes(CANARY));
        messages.push(error.message);
        return true;
      },
    );
  }
  assert.equal(new Set(messages).size, 1);
});

test("the same key in different objects is valid", () => {
  const source =
    '{"data":{"value":1},"children":[{"data":{"value":2}},{"data":{"value":3}}]}';
  assert.deepEqual(parsePublicJson(source), JSON.parse(source));
  const special = '{"__proto__":{"value":1},"constructor":2,"toString":3}';
  const parsed = parsePublicJson(special);
  assert.deepEqual(parsed, JSON.parse(special));
  assert.equal(Object.getPrototypeOf(parsed), Object.prototype);
  assert.ok(Object.hasOwn(parsed, "__proto__"));
});

test("quoted braces, escaped quotes and backslashes do not change object scope", () => {
  const values = [
    '{"data":"inside a string","data":"also inside a string"}',
    '},],{[,:"',
    'quote: "; slash: \\; escaped quote: \\"; end',
    "\\",
    "\\\\",
    '"',
    '\\"',
    '\\\\"',
    "\b\f\n\r\t",
    "line\u2028separator\u2029end",
  ];
  const source = JSON.stringify({ values, nested: { values }, final: true });
  assert.deepEqual(parsePublicJson(source), JSON.parse(source));
  assert.deepEqual(
    parsePublicJson(' \r\n\t [ { "a" : 1 }, { "a" : 2 } ] \t\n'),
    [{ a: 1 }, { a: 2 }],
  );
});

test("Unicode keys use decoded exact code points without normalization", () => {
  const source = String.raw`{"\u0061":1,"é":2,"\ud83d\ude00":3,"Å":4,"Å":5,"a\\b":6,"a\u0008":7}`;
  assert.deepEqual(parsePublicJson(source), JSON.parse(source));
  reject(String.raw`{"é":"PRIVATE_CANARY","\u00e9":"public"}`);
});

test("numbers retain native JSON semantics and source bytes remain unchanged", () => {
  const bytes = Buffer.from(
    '{"whole":1.0,"small":1e-06,"negativeZero":-0,"fraction":0.10000000000000001,"large":1e+20}\n',
  );
  const original = Buffer.from(bytes);
  const value = parsePublicJson(bytes);
  assert.deepEqual(value, JSON.parse(bytes));
  assert.ok(Object.is(value.negativeZero, -0));
  assert.deepEqual(bytes, original);
  for (const source of [
    "null",
    "true",
    "false",
    "0",
    "-0",
    "1.0",
    "1e-06",
    '"text"',
    "[]",
    "{}",
  ]) {
    assert.deepEqual(parsePublicJson(source), JSON.parse(source));
  }
});

test("UTF-8 Buffer and Uint8Array views use exactly their supplied bytes", () => {
  const encoded = Buffer.from('{"é":"😀"}');
  assert.deepEqual(parsePublicJson(encoded), { é: "😀" });
  const padded = Buffer.concat([
    Buffer.from("PRIVATE_CANARY"),
    encoded,
    Buffer.from("PRIVATE_CANARY"),
  ]);
  const view = new Uint8Array(
    padded.buffer,
    padded.byteOffset + CANARY.length,
    encoded.length,
  );
  assert.deepEqual(parsePublicJson(view), { é: "😀" });
  for (const bad of [
    [0xc0, 0xaf],
    [0xed, 0xa0, 0x80],
    [0xf4, 0x90, 0x80, 0x80],
    [0x80],
    [0xe2, 0x82],
  ]) {
    reject(Buffer.from([0x22, ...bad, 0x22]));
    reject(new Uint8Array([0x22, ...bad, 0x22]));
  }
});

test("nesting depth counts both objects and arrays with a 64-container limit", () => {
  const arrays = (depth) => "[".repeat(depth) + "0" + "]".repeat(depth);
  const objects = (depth) => '{"a":'.repeat(depth) + "0" + "}".repeat(depth);
  const mixed = (depth) => {
    let source = "0";
    for (let i = 0; i < depth; i++)
      source = i % 2 ? `[${source}]` : `{"a":${source}}`;
    return source;
  };
  for (const build of [arrays, objects, mixed]) {
    assert.deepEqual(parsePublicJson(build(64)), JSON.parse(build(64)));
    reject(build(65));
  }
});

test("large sorted objects parse with all keys intact", () => {
  const original = Object.fromEntries(
    Array.from({ length: 10000 }, (_, i) => [
      `key_${String(i).padStart(5, "0")}`,
      i,
    ]),
  );
  const source = Buffer.from(JSON.stringify(original));
  assert.deepEqual(parsePublicJson(source), original);
  const duplicate = `${source.toString("utf8").slice(0, -1)},"key_00000":"PRIVATE_CANARY"}`;
  reject(duplicate);
});

test("malformed JSON errors do not disclose input", () => {
  for (const source of [
    `PRIVATE_CANARY`,
    `{"value":"${CANARY}",}`,
    `{"value":"${CANARY}"} trailing`,
    `{"value":"${CANARY}`,
    `{"value":"${CANARY}\\x20"}`,
    `{"value":"${CANARY}\n"}`,
    `{"value":"${CANARY}","n":01}`,
    `{"value":"${CANARY}","n":+1}`,
    `{"value":"${CANARY}","n":1.}`,
    `{"value":"${CANARY}","n":1e}`,
    `{"value":"${CANARY}","n":undefined}`,
    `{"value":"${CANARY}","items":[1,]}`,
  ])
    reject(source);
});

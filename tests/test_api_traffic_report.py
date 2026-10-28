from __future__ import annotations

import csv
import contextlib
import importlib.util
import io
import json
import os
import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("api_traffic_report", ROOT / "scripts/api_traffic_report.py")
reporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reporter)


def request(**changes):
    row = dict(type="api_request", version=1, timestamp="2026-10-25T12:00:00Z",
               request_id="e7e5afae-7b12-4ab4-8efc-8842892d852a", transport="api",
               operation="get_target", species="human", method="GET", status=200,
               error_code="none", outcome="finished", admitted=True,
               duration_ms=4, bytes=100, client_id="a" * 32)
    row.update(changes)
    return row


def aggregate(**changes):
    row = request(type="api_request_aggregate", timestamp="2026-10-25T12:05:00Z",
                  first_timestamp="2026-10-25T12:00:00Z", last_timestamp="2026-10-25T12:04:00Z",
                  count=3, duration_ms_total=12, bytes_total=300,
                  client_ids_available=False, dimensions_collapsed=False)
    row.pop("client_id")
    row.pop("request_id")
    row.pop("duration_ms")
    row.pop("bytes")
    row.update(changes)
    return row


def delivery(**changes):
    row = dict(type="analytics_delivery", version=1, timestamp="2026-10-25T12:00:01Z",
               transport="api", outcome="accepted", request_id=None)
    row.update(changes)
    return row


class ApiTrafficReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=Path(tempfile.gettempdir()).resolve())
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "export.log"

    def summarize(self, records, **options):
        self.path.write_text("\n".join(json.dumps(row) for row in records) + "\n", encoding="utf-8")
        return reporter.summarize([self.path], **options)

    def test_requests_aggregates_refusals_and_delivery_have_independent_totals(self):
        report = self.summarize([
            request(), request(transport="mcp", operation="search_targets", species="mouse", client_id="b" * 32),
            request(status=429, error_code="rate_limited", admitted=False, client_id=None),
            request(status=0, error_code="request_aborted", outcome="aborted", client_id=None),
            aggregate(transport="mcp", operation="get_sequence", species="mouse"),
            delivery(), delivery(outcome="provider_dropped"), delivery(outcome="failed"),
            delivery(outcome="policy_skipped", reason="not_admitted"),
            delivery(type="analytics_delivery_aggregate", count=4, outcome="capacity_skipped",
                     first_timestamp="2026-10-25T12:00:00Z", last_timestamp="2026-10-25T12:00:01Z"),
        ])
        self.assertEqual(report["totals"]["requests"], 7)
        self.assertEqual(report["totals"]["errors"], 2)
        self.assertEqual(report["totals"]["rate_limited"], 1)
        self.assertEqual(report["totals"]["rejected"], 1)
        self.assertEqual(report["totals"]["aborted"], 1)
        self.assertEqual(report["totals"]["bytes"], 700)
        self.assertEqual(report["totals"]["duration_ms"], 28)
        self.assertEqual(report["coverage"]["client_identity_percent"], 28.57)
        self.assertEqual(report["coverage"]["daily_client_sum"], 2)
        self.assertEqual(report["delivery_totals"]["accepted"], 1)
        self.assertEqual(report["delivery_totals"]["provider_dropped"], 1)
        self.assertEqual(report["delivery_totals"]["capacity_skipped"], 4)
        self.assertEqual(report["delivery_totals"]["policy_skipped"], 1)
        self.assertTrue(any(row["reason"] == "not_admitted" for row in report["deliveries"]))
        self.assertEqual(sum(row["requests"] for row in report["requests"]), 7)

    def test_client_identity_is_daily_utc_even_with_display_timezone(self):
        report = self.summarize([
            request(timestamp="2026-10-25T00:01:00Z"),
            request(timestamp="2026-10-25T00:02:00Z"),
            request(timestamp="2026-10-26T00:01:00Z"),
            request(timestamp="2026-10-26T00:02:00Z", client_id=None),
        ], tz="America/New_York")
        self.assertEqual(report["daily_clients_utc"], [
            {"day": "2026-10-25", "approximate_clients": 1},
            {"day": "2026-10-26", "approximate_clients": 1},
        ])
        self.assertEqual(report["coverage"]["client_identity_percent"], 75)
        self.assertTrue(report["requests"][0]["hour"].startswith("2026-10-24T20:00:00-04:00"))

    def test_mcp_semantic_error_in_http_200_counts_as_error_across_formats(self):
        # Same controlled record emitted by the producer for an MCP lookup failure.
        report = self.summarize([request(transport="mcp", operation="get_target", method="POST",
                                         status=200, error_code="not_found", species="mouse")])
        self.assertEqual(report["totals"]["requests"], 1)
        self.assertEqual(report["totals"]["errors"], 1)
        self.assertEqual(report["requests"][0]["status"], 200)
        self.assertEqual(report["requests"][0]["transport"], "mcp")
        rows = list(reporter.csv_rows(report))
        self.assertIn(("totals", "{}", "errors", 1), rows)
        markup = reporter.render_html(report)
        self.assertIn("sum('errors')", markup)

    def test_sdk_protocol_error_is_an_error_even_with_http_200(self):
        report = self.summarize([request(transport="mcp", operation="call_tool", method="POST",
                                         status=200, error_code="protocol_error", client_id=None)])
        self.assertEqual(report["totals"]["requests"], 1)
        self.assertEqual(report["totals"]["errors"], 1)
        self.assertEqual(report["requests"][0]["errors"], 1)
        self.assertEqual(report["input_quality"]["malformed_lines"], 0)

    def test_discovery_and_daily_budget_keep_complete_counts_and_finite_reasons(self):
        report = self.summarize([
            request(operation="discover", transport="mcp"),
            delivery(outcome="policy_skipped", reason="operation_not_exported"),
            delivery(outcome="capacity_skipped", reason="daily_budget"),
            delivery(type="analytics_delivery_aggregate", count=3, outcome="capacity_skipped",
                     reason="daily_budget", first_timestamp="2026-10-25T12:00:00Z",
                     last_timestamp="2026-10-25T12:00:01Z"),
            delivery(outcome="capacity_skipped"),
            delivery(outcome="capacity_skipped", reason="PRIVATE_REASON"),
        ])
        self.assertEqual(report["totals"]["requests"], 1)
        self.assertEqual(report["delivery_totals"]["capacity_skipped"], 5)
        self.assertEqual(sum(row["count"] for row in report["deliveries"] if row["reason"] == "daily_budget"), 4)
        self.assertEqual(report["input_quality"]["malformed_lines"], 1)
        self.assertNotIn("PRIVATE_REASON", json.dumps(report))

    def test_strict_whitelist_and_all_formats_never_export_private_source_fields(self):
        secret = '=HYPERLINK("https://private.invalid?q=PRIVATE_SEQUENCE")</script><img src=x onerror=alert(1)>'
        report = self.summarize([
            request(raw_ip="198.51.100.4", ua=secret, search=secret, query=secret,
                    body=secret, token=secret, target_id=secret, country=secret, city=secret),
            request(operation=secret), request(error_code=secret), request(species=secret),
            request(client_id=secret), request(status=True), request(admitted="false"),
            request(duration_ms=-1), request(version=True), request(type=[secret]),
            {"message": secret},
        ])
        self.assertEqual(report["totals"]["requests"], 1)
        self.assertEqual(report["input_quality"]["malformed_lines"], 8)
        reporter.write_reports(report, Path(self.tmp.name) / "out")
        for name in ("report.json", "report.csv", "report.html"):
            data = (Path(self.tmp.name) / "out" / name).read_text()
            for value in ("PRIVATE_SEQUENCE", "198.51.100.4", "private.invalid", "a" * 32,
                          "e7e5afae-7b12-4ab4-8efc-8842892d852a"):
                self.assertNotIn(value, data)

    def test_prefixed_host_lines_malformed_and_oversized_are_bounded(self):
        self.path.write_bytes(b"host startup\n" + b"2026-10-25 INFO " + json.dumps(request()).encode() + b"\n"
                             + b'{"type":"api_request",bad}\n' + b'{"type":"api_request","version":1}\n'
                             + b"x" * (reporter.MAX_LINE_BYTES * 3) + b"\n"
                             + b"\xff\n" + json.dumps(request()).encode() + b"\n")
        report = reporter.summarize([self.path])
        self.assertEqual(report["input_quality"]["lines"], 7)
        self.assertEqual(report["totals"]["requests"], 2)
        self.assertEqual(report["input_quality"]["oversized_lines"], 1)
        self.assertEqual(report["input_quality"]["malformed_lines"], 3)
        self.assertEqual(report["input_quality"]["ignored_lines"], 1)

    def test_filter_boundaries_do_not_guess_aggregate_distribution(self):
        report = self.summarize([
            request(timestamp="2026-10-25T11:59:59Z"), request(),
            request(timestamp="2026-10-25T12:03:00Z"), aggregate(),
            delivery(type="analytics_delivery_aggregate", count=2, timestamp="2026-10-25T12:05:00Z",
                     first_timestamp="2026-10-25T12:00:00Z", last_timestamp="2026-10-25T12:04:00Z"),
        ], since=reporter.instant("2026-10-25T12:00:00Z"), until=reporter.instant("2026-10-25T12:03:00Z"))
        self.assertEqual(report["totals"]["requests"], 1)
        self.assertEqual(report["input_quality"]["filtered_records"], 2)
        self.assertEqual(report["input_quality"]["boundary_requests"], 3)
        self.assertEqual(report["input_quality"]["boundary_deliveries"], 2)

    def test_collapsed_server_dimensions_do_not_manufacture_refused_requests(self):
        report = self.summarize([request(), aggregate(dimensions_collapsed=True,
            transport="unknown", operation="unknown", species="none", method="OTHER",
            status=0, error_code="unknown", admitted=False)])
        self.assertEqual(report["totals"]["requests"], 4)
        self.assertEqual(report["totals"]["rejected"], 0)
        self.assertEqual(report["coverage"]["request_dimensions_percent"], 25)

    def test_cardinality_caps_preserve_request_and_delivery_counts(self):
        with patch.object(reporter, "MAX_BUCKETS", 1), patch.object(reporter, "MAX_CLIENTS", 1):
            report = self.summarize([
                request(), request(species="mouse", client_id="b" * 32),
                request(transport="mcp", client_id="c" * 32),
                delivery(), delivery(transport="mcp"), delivery(outcome="failed"),
            ])
        self.assertLessEqual(len(report["requests"]), 2)
        self.assertEqual(sum(row["requests"] for row in report["requests"]), 3)
        self.assertEqual(report["totals"]["requests"], 3)
        self.assertEqual(report["totals"]["bucket_overflow_requests"], 2)
        self.assertTrue(report["coverage"]["client_counts_are_lower_bounds"])
        self.assertEqual(report["coverage"]["daily_client_sum"], 1)
        self.assertEqual(sum(report["delivery_totals"].values()), 3)

    def test_json_csv_and_html_embed_identical_report_and_no_external_assets(self):
        report = self.summarize([request(), aggregate(), delivery()])
        reporter.write_reports(report, Path(self.tmp.name) / "out")
        output = Path(self.tmp.name) / "out"
        self.assertEqual(json.loads((output / "report.json").read_text()), report)
        with (output / "report.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        totals = {row["metric"]: int(row["value"]) for row in rows if row["section"] == "totals"}
        self.assertEqual(totals, report["totals"])
        self.assertEqual(sum(int(row["value"]) for row in rows if row["section"] == "requests" and row["metric"] == "requests"), 4)
        markup = (output / "report.html").read_text()
        embedded = markup.split('<script type="application/json" id="report-data">')[1].split("</script>")[0]
        self.assertEqual(json.loads(embedded), report)
        self.assertIn("connect-src 'none'", markup)
        for external in ("<script src=", "<link ", "fetch(", "XMLHttpRequest", "https://cdn", "<img", "<iframe"):
            self.assertNotIn(external, markup)
        for control in ('id="transport"', 'id="species"', 'id="operation"', 'id="status"', 'id="from"', 'id="to"'):
            self.assertIn(control, markup)

    def test_html_json_and_csv_defend_against_injected_text(self):
        report = self.summarize([])
        attack = '</script><script>alert("private")</script>&'
        report["notes"] = [attack]
        markup = reporter.render_html(report)
        self.assertNotIn(attack, markup)
        self.assertIn("\\u003c/script\\u003e", markup)
        self.assertIn("&lt;/script&gt;", markup)
        for value in ("=cmd()", " +cmd()", "-cmd()", "@cmd()", "\t=cmd()", "\r=cmd()"):
            self.assertTrue(reporter.csv_safe(value).startswith("'"))

    def test_bad_intervals_and_duplicate_input_files_fail(self):
        for row in (request(timestamp="2026-10-25T12:00:00"), aggregate(count=True),
                    aggregate(count=0), aggregate(client_ids_available=True),
                    aggregate(first_timestamp="2026-10-26T12:00:00Z"),
                    delivery(outcome="made_up"), request(status=99),
                    delivery(outcome="policy_skipped", reason="made_up")):
            with self.assertRaises(ValueError):
                reporter.validate(row)
        self.path.write_text("")
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            reporter.main([str(self.path), str(self.path), "--output-dir", str(Path(self.tmp.name) / "out")])

    def test_output_cannot_replace_source_and_timezone_overflow_is_malformed(self):
        source = Path(self.tmp.name) / "report.json"
        original = json.dumps(request()) + "\n"
        source.write_text(original)
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            reporter.main([str(source), "--output-dir", self.tmp.name])
        self.assertEqual(source.read_text(), original)
        report = self.summarize([request(timestamp="9999-12-31T23:59:59Z")], tz="Pacific/Kiritimati")
        self.assertEqual(report["input_quality"]["malformed_lines"], 1)

    def test_writer_preflights_every_destination_and_preserves_existing_link_targets(self):
        report = self.summarize([request()])
        source_bytes = self.path.read_bytes()
        for kind in ("regular", "symlink", "hardlink"):
            for name in ("report.json", "report.csv", "report.html"):
                with self.subTest(kind=kind, name=name):
                    directory = Path(self.tmp.name) / (kind + name)
                    directory.mkdir()
                    outside = directory.parent / (kind + name + ".outside")
                    outside.write_bytes(b"PRIVATE_USER_MATERIAL_PRESERVE")
                    destination = directory / name
                    if kind == "symlink":
                        destination.symlink_to(outside)
                    elif kind == "hardlink":
                        os.link(outside, destination)
                    else:
                        destination.write_bytes(b"PREVIOUS_REPORT_PRESERVE")
                    before = destination.read_bytes()
                    with self.assertRaises(FileExistsError):
                        reporter.write_reports(report, directory)
                    self.assertEqual(outside.read_bytes(), b"PRIVATE_USER_MATERIAL_PRESERVE")
                    self.assertEqual(destination.read_bytes(), before)
                    self.assertEqual(self.path.read_bytes(), source_bytes)
                    self.assertEqual([entry.name for entry in directory.iterdir()], [name])

    def test_writer_refuses_directory_symlinks_and_symlinked_parents(self):
        report = self.summarize([request()])
        target = Path(self.tmp.name) / "outside"
        target.mkdir()
        (target / "canary").write_bytes(b"PRESERVE")
        link = Path(self.tmp.name) / "linked-directory"
        link.symlink_to(target, target_is_directory=True)
        for output in (link, link / "nested" / "out"):
            with self.subTest(output=output), self.assertRaises(OSError):
                reporter.write_reports(report, output)
        self.assertEqual([entry.name for entry in target.iterdir()], ["canary"])
        self.assertEqual((target / "canary").read_bytes(), b"PRESERVE")

    def test_writer_creates_private_files_and_rolls_back_failed_publication(self):
        report = self.summarize([request()])
        success = Path(self.tmp.name) / "private-output"
        reporter.write_reports(report, success)
        for path in success.iterdir():
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600)
            self.assertEqual(path.stat().st_nlink, 1)
        failure = Path(self.tmp.name) / "failed-output"
        real_link = os.link
        def fail_third(source, destination, **kwargs):
            if destination == "report.html":
                raise OSError("synthetic publication failure")
            return real_link(source, destination, **kwargs)
        with patch.object(reporter.os, "link", side_effect=fail_third), self.assertRaises(OSError):
            reporter.write_reports(report, failure)
        self.assertEqual(list(failure.iterdir()), [])
        partial = Path(self.tmp.name) / "failed-stage"
        with patch.object(reporter, "csv_rows", side_effect=ValueError("synthetic stage failure")), self.assertRaises(ValueError):
            reporter.write_reports(report, partial)
        self.assertEqual(list(partial.iterdir()), [])

    def test_writer_refuses_a_late_destination_without_overwriting_or_leaving_partial_reports(self):
        report = self.summarize([request()])
        output = Path(self.tmp.name) / "raced-output"
        outside = Path(self.tmp.name) / "outside-raced"
        outside.write_bytes(b"USER_MATERIAL")
        real_link = os.link
        def late_hardlink(source, destination, **kwargs):
            if destination == "report.html":
                real_link(outside, destination, dst_dir_fd=kwargs["dst_dir_fd"])
            return real_link(source, destination, **kwargs)
        with patch.object(reporter.os, "link", side_effect=late_hardlink), self.assertRaises(FileExistsError):
            reporter.write_reports(report, output)
        self.assertEqual(outside.read_bytes(), b"USER_MATERIAL")
        self.assertEqual([entry.name for entry in output.iterdir()], ["report.html"])
        self.assertEqual((output / "report.html").read_bytes(), b"USER_MATERIAL")

    @unittest.skipUnless(shutil.which("node"), "Node is needed for the offline DOM regression")
    def test_overflow_dom_keeps_known_outcomes_and_marks_partial_filters(self):
        with patch.object(reporter, "MAX_BUCKETS", 1), patch.object(reporter, "MAX_CLIENTS", 1):
            report = self.summarize([
                request(), request(species="mouse", status=429, error_code="rate_limited", admitted=False, client_id="b" * 32),
                aggregate(error_code="protocol_error", status=200, transport="mcp", method="POST"),
            ], until=reporter.instant("2026-10-25T12:03:00Z"))
        self.assertEqual(report["totals"]["errors"], 1)
        self.assertEqual(report["input_quality"]["boundary_requests"], 3)
        self.assertEqual(report["coverage"]["request_dimensions_percent"], 100)
        self.assertEqual(report["coverage"]["report_grouping_percent"], 50)
        self.assertTrue(report["coverage"]["client_counts_are_lower_bounds"])
        for metric in ("requests", "errors", "rate_limited", "rejected", "aborted"):
            self.assertEqual(sum(row[metric] for row in report["requests"]), report["totals"][metric])
            self.assertEqual(sum(row[3] for row in reporter.csv_rows(report) if row[0] == "requests" and row[2] == metric), report["totals"][metric])
        markup = reporter.render_html(report)
        source = "const html=" + json.dumps(markup) + ";\n" + r'''
const vm=require('node:vm'), assert=require('node:assert/strict');
class El {constructor(){this.children=[];this.listeners={};this.textContent='';this.value='';this.style={};}append(...items){this.children.push(...items);}replaceChildren(...items){this.children=items;}addEventListener(key,value){this.listeners[key]=value;}get lastChild(){return this.children.at(-1);}}
const elements=new Map(), byId=id=>{if(!elements.has(id))elements.set(id,new El());return elements.get(id);};
for(const id of ['trend','operations','clients','delivery','policy']){byId(id).tBodies=[new El()];byId(id).tHead={rows:[new El()]};}
const embedded=html.split('<script type="application/json" id="report-data">')[1].split('</script>')[0];byId('report-data').textContent=embedded;
const document={getElementById:byId,createElement:()=>new El(),createTextNode:textContent=>({textContent})};
vm.runInNewContext(html.split('<script>')[1].split('</script>')[0],{document},{timeout:1000});
const counts=()=>byId('cards').children.map(el=>el.children[0].textContent);
assert.deepEqual(counts(),['2','1','1','1','0','0']);
assert.ok(byId('coverage').textContent.includes('Source outcome coverage: 100%'));
assert.ok(byId('coverage').textContent.includes('Retained report grouping coverage: 50%'));
byId('species').value='mouse';byId('species').listeners.change();
assert.deepEqual(counts(),['0','0','0','0','0','0']);
assert.ok(byId('filter-coverage').textContent.includes('Partial filter coverage: 1'));
byId('species').value='human';byId('species').listeners.change();
assert.deepEqual(counts(),['1','0','0','0','0','0']);
byId('reset').listeners.click();assert.deepEqual(counts(),['2','1','1','1','0','0']);
'''
        result = subprocess.run([shutil.which("node"), "-"], input=source, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()

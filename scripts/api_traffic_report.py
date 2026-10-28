#!/usr/bin/env python3
"""Summarize exported OpenAntigens application logs locally, with no network calls."""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import secrets
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

OPERATIONS = frozenset("unknown get_release search_targets get_target get_section list_constructs get_sequence get_downloads health ready openapi initialize initialized ping list_tools call_tool discover".split())
ERRORS = frozenset("none unknown invalid_arguments invalid_cursor body_not_allowed species_mismatch refine_query handler_timeout invalid_host invalid_origin method_not_allowed not_acceptable not_found rate_limited request_too_large result_too_large unavailable unsupported_media_type publication_in_progress retired_release unsupported_capability request_aborted protocol_error".split())
DELIVERIES = frozenset("accepted provider_dropped failed timeout capacity_skipped identity_unavailable disabled policy_skipped".split())
DELIVERY_REASONS = frozenset("not_admitted operation_not_exported request_failed not_finished".split())
TRANSPORTS = frozenset(("api", "mcp", "unknown"))
TYPES = frozenset(("api_request", "api_request_aggregate", "analytics_delivery", "analytics_delivery_aggregate"))
MAX_LINE_BYTES = 65536
MAX_BUCKETS = 20000
MAX_CLIENTS = 100000
MAX_COUNT = 10**15
DIMENSIONS = ("hour", "transport", "operation", "species", "method", "status", "error_code", "outcome", "admitted", "dimensions_collapsed")
REQUEST_METRICS = ("requests", "errors", "rate_limited", "rejected", "aborted", "duration_ms", "bytes", "client_covered_requests", "outcome_covered_requests", "grouping_covered_requests")
CLIENT_RE = re.compile(r"[a-f0-9]{32}\Z")


def instant(value: object) -> datetime:
    if not isinstance(value, str) or len(value) > 40:
        raise ValueError("invalid timestamp")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("timestamp needs a UTC offset")
    return parsed.astimezone(timezone.utc)


def utc_text(value: datetime | None) -> str | None:
    return value.isoformat().replace("+00:00", "Z") if value else None


def integer(value: object, minimum: int = 0, maximum: int = MAX_COUNT) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError("invalid integer")
    return value


def choice(value: object, allowed: frozenset) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ValueError("invalid category")
    return value


def boolean(value: object) -> bool:
    if type(value) is not bool:
        raise ValueError("invalid boolean")
    return value


def validate(record: dict) -> dict:
    """Copy only controlled fields; untrusted extras never reach a report."""
    kind = choice(record.get("type"), TYPES)
    if type(record.get("version")) is not int or record["version"] != 1:
        raise ValueError("unsupported record version")
    stamp = instant(record.get("timestamp"))
    aggregate = kind.endswith("_aggregate")
    first = instant(record.get("first_timestamp")) if aggregate else stamp
    last = instant(record.get("last_timestamp")) if aggregate else stamp
    if first > last or last > stamp:
        raise ValueError("invalid interval")
    clean = {"type": kind, "first": first, "last": last,
             "transport": choice(record.get("transport"), TRANSPORTS),
             "count": integer(record.get("count"), 1) if aggregate else 1}
    if kind.startswith("analytics_delivery"):
        clean["outcome"] = choice(record.get("outcome"), DELIVERIES)
        clean["reason"] = choice(record.get("reason"), DELIVERY_REASONS) if clean["outcome"] == "policy_skipped" else "none"
        if clean["outcome"] == "capacity_skipped" and record.get("reason") is not None:
            clean["reason"] = choice(record["reason"], frozenset(("daily_budget",)))
        return clean
    clean.update({
        "operation": choice(record.get("operation"), OPERATIONS),
        "species": choice(record.get("species"), frozenset(("human", "mouse", "none"))),
        "method": choice(record.get("method"), frozenset(("GET", "HEAD", "POST", "OPTIONS", "OTHER"))),
        "error_code": choice(record.get("error_code"), ERRORS),
        "outcome": choice(record.get("outcome"), frozenset(("finished", "aborted"))),
        "admitted": boolean(record.get("admitted")),
        "dimensions_collapsed": boolean(record.get("dimensions_collapsed", False)),
        "duration_ms": integer(record.get("duration_ms_total" if aggregate else "duration_ms")),
        "bytes": integer(record.get("bytes_total" if aggregate else "bytes")),
    })
    status = integer(record.get("status"), 0, 599)
    if status and status < 100:
        raise ValueError("invalid status")
    clean["status"] = status
    client = None if aggregate else record.get("client_id")
    if client is not None and (not isinstance(client, str) or not CLIENT_RE.fullmatch(client)):
        raise ValueError("invalid client ID")
    if aggregate and record.get("client_ids_available") is not False:
        raise ValueError("aggregate cannot carry clients")
    clean["client_id"] = client
    return clean


def bounded_lines(source):
    """Drain oversized lines in bounded chunks without allocating the full line."""
    while line := source.readline(MAX_LINE_BYTES + 1):
        if len(line) > MAX_LINE_BYTES:
            while not line.endswith(b"\n"):
                line = source.readline(MAX_LINE_BYTES + 1)
                if not line:
                    break
            yield None
        else:
            yield line


def summarize(paths: list[Path], *, since: datetime | None = None,
              until: datetime | None = None, tz: str = "UTC") -> dict:
    zone = ZoneInfo(tz)
    totals = Counter(requests=0, errors=0, rate_limited=0, rejected=0, aborted=0,
                     client_covered_requests=0, dimension_covered_requests=0,
                     report_dimension_covered_requests=0,
                     duration_ms=0, bytes=0, aggregate_requests=0,
                     bucket_overflow_requests=0, client_capacity_requests=0)
    quality = Counter(lines=0, ignored_lines=0, malformed_lines=0, oversized_lines=0,
                      filtered_records=0, boundary_records=0, boundary_requests=0,
                      boundary_deliveries=0)
    buckets: dict[tuple, Counter] = {}
    delivery_buckets: Counter = Counter()
    clients: set[tuple[str, str]] = set()
    daily_clients: Counter = Counter()
    first_seen = last_seen = None
    for path in paths:
        with path.open("rb") as source:
            for line in bounded_lines(source):
                quality["lines"] += 1
                if line is None:
                    quality["oversized_lines"] += 1
                    continue
                try:
                    text = line.decode("utf-8")
                    # GoDaddy's text view may prefix each JSON message with a timestamp.
                    start = text.find("{")
                    if start < 0:
                        quality["ignored_lines"] += 1
                        continue
                    record = json.loads(text[start:])
                except (ValueError, UnicodeError, RecursionError):
                    quality["malformed_lines"] += 1
                    continue
                if not isinstance(record, dict) or not isinstance(record.get("type"), str) or record["type"] not in TYPES:
                    quality["ignored_lines"] += 1
                    continue
                try:
                    event = validate(record)
                    hour = event["last"].astimezone(zone).replace(minute=0, second=0, microsecond=0).isoformat()
                except (ValueError, TypeError, OverflowError, RecursionError):
                    quality["malformed_lines"] += 1
                    continue
                first, last, count = event["first"], event["last"], event["count"]
                if (since and last < since) or (until and first >= until):
                    quality["filtered_records"] += 1
                    continue
                if (since and first < since) or (until and last >= until):
                    quality["boundary_records"] += 1
                    quality["boundary_deliveries" if event["type"].startswith("analytics") else "boundary_requests"] += count
                    continue
                if event["type"].startswith("analytics"):
                    key = (hour, event["transport"], event["outcome"], event["reason"])
                    if key not in delivery_buckets and len(delivery_buckets) >= MAX_BUCKETS:
                        key = ("overflow", "unknown", event["outcome"], event["reason"])
                    delivery_buckets[key] += count
                    continue
                first_seen = min(first_seen, first) if first_seen else first
                last_seen = max(last_seen, last) if last_seen else last
                totals["requests"] += count
                totals["duration_ms"] += event["duration_ms"]
                totals["bytes"] += event["bytes"]
                if event["type"].endswith("_aggregate"):
                    totals["aggregate_requests"] += count
                if not event["dimensions_collapsed"]:
                    totals["dimension_covered_requests"] += count
                    totals["errors"] += count * (event["status"] >= 400 or event["error_code"] != "none")
                    totals["rate_limited"] += count * (event["status"] == 429)
                    totals["rejected"] += count * (not event["admitted"])
                    totals["aborted"] += count * (event["outcome"] == "aborted")
                client = event["client_id"]
                if client is not None:
                    totals["client_covered_requests"] += count
                    day = first.date().isoformat()  # HMAC identities reset at UTC midnight.
                    client_key = (day, client)
                    if client_key not in clients:
                        if len(clients) < MAX_CLIENTS:
                            clients.add(client_key)
                            daily_clients[day] += 1
                        else:
                            totals["client_capacity_requests"] += count
                key = (hour,) + tuple(event[k] for k in DIMENSIONS[1:])
                grouping_known = not event["dimensions_collapsed"]
                if key not in buckets and len(buckets) >= MAX_BUCKETS:
                    key = ("overflow", "unknown", "unknown", "none", "OTHER", 0, "unknown", "finished", False, True)
                    totals["bucket_overflow_requests"] += count
                    grouping_known = False
                totals["report_dimension_covered_requests"] += count * grouping_known
                bucket = buckets.setdefault(key, Counter({metric: 0 for metric in REQUEST_METRICS}))
                bucket["requests"] += count
                bucket["duration_ms"] += event["duration_ms"]
                bucket["bytes"] += event["bytes"]
                bucket["client_covered_requests"] += count * (client is not None)
                bucket["outcome_covered_requests"] += count * (not event["dimensions_collapsed"])
                bucket["grouping_covered_requests"] += count * grouping_known
                if not event["dimensions_collapsed"]:
                    bucket["errors"] += count * (event["status"] >= 400 or event["error_code"] != "none")
                    bucket["rate_limited"] += count * (event["status"] == 429)
                    bucket["rejected"] += count * (not event["admitted"])
                    bucket["aborted"] += count * (event["outcome"] == "aborted")
    requests = totals["requests"]
    delivery_totals = Counter({key: 0 for key in sorted(DELIVERIES)})
    for (_, _, outcome, _), count in delivery_buckets.items():
        delivery_totals[outcome] += count
    return {
        "schema_version": 1,
        "parameters": {"timezone": tz, "since_inclusive": utc_text(since), "until_exclusive": utc_text(until),
                       "max_line_bytes": MAX_LINE_BYTES, "max_buckets": MAX_BUCKETS, "max_daily_client_entries": MAX_CLIENTS},
        "observed": {"first_request_utc": utc_text(first_seen), "last_request_utc": utc_text(last_seen)},
        "totals": dict(totals), "input_quality": dict(quality),
        "coverage": {"client_identity_percent": round(100 * totals["client_covered_requests"] / requests, 2) if requests else 0,
                     "request_dimensions_percent": round(100 * totals["dimension_covered_requests"] / requests, 2) if requests else 0,
                     "report_grouping_percent": round(100 * totals["report_dimension_covered_requests"] / requests, 2) if requests else 0,
                     "daily_client_sum": sum(daily_clients.values()),
                     "client_counts_are_lower_bounds": bool(totals["client_capacity_requests"])},
        "daily_clients_utc": [{"day": day, "approximate_clients": count} for day, count in sorted(daily_clients.items())],
        "requests": [dict(zip(DIMENSIONS, key), **dict(value)) for key, value in sorted(buckets.items())],
        "delivery_totals": dict(delivery_totals),
        "deliveries": [{"hour": hour, "transport": transport, "outcome": outcome, "reason": reason, "count": count}
                       for (hour, transport, outcome, reason), count in sorted(delivery_buckets.items())],
        "notes": ["Counts cover validated exported application records, including bots and refused requests. Host omissions, restarts and missing exports limit completeness.",
                  "Anonymous daily clients are approximate and do not establish people. Shared identities merge clients; changing identities split them. Counts reset at UTC midnight and do not add to unique visitors.",
                  "Plausible delivery acceptance does not prove dashboard inclusion. Plausible can filter bots. Country, city and eligible visitors remain in the existing Plausible dashboard.",
                  "The application ledger is on by default; Plausible export is off by default. Enabled export is partial: eligible operations share a maximum of 100 delivery attempts per process per UTC day. Restarts and replicas can increase the account total; this is not an account quota.",
                  "Aggregates count each contained request once and carry no client identities. Aggregates crossing a filter boundary are excluded and counted separately. Aggregates spanning hours are assigned to their last request hour.",
                  "Input files must contain non-overlapping exports from one environment. Repeated or overlapping exports double-count traffic. No source paths, request IDs or client IDs are included."],
    }


def csv_safe(value: object) -> str:
    text = str(value)
    return "'" + text if text.lstrip().startswith(("=", "+", "-", "@", "\t", "\r")) else text


def csv_rows(report: dict):
    yield ("section", "dimensions_json", "metric", "value")
    for section in ("parameters", "observed", "totals", "input_quality", "coverage", "delivery_totals"):
        for metric, value in report[section].items():
            yield (section, "{}", metric, value)
    for row in report["requests"]:
        dimensions = {key: row[key] for key in DIMENSIONS}
        for metric in REQUEST_METRICS:
            yield ("requests", json.dumps(dimensions, sort_keys=True), metric, row[metric])
    for section, metrics in (("deliveries", ("count",)), ("daily_clients_utc", ("approximate_clients",))):
        for row in report[section]:
            dimensions = {key: value for key, value in row.items() if key not in metrics}
            for metric in metrics:
                yield (section, json.dumps(dimensions, sort_keys=True), metric, row[metric])


def render_html(report: dict) -> str:
    # Escape HTML delimiters in JSON as well as using textContent for all rendered data.
    data = json.dumps(report, ensure_ascii=True, separators=(",", ":")).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    notes = "".join(f"<li>{html.escape(note)}</li>" for note in report["notes"])
    return HTML.replace("__REPORT_DATA__", data).replace("__NOTES__", notes)


HTML = '''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; connect-src 'none'; img-src 'none'; base-uri 'none'; form-action 'none'">
<title>OpenAntigens API and MCP traffic</title>
<style>body{font:16px system-ui,sans-serif;max-width:1250px;margin:2rem auto;padding:0 1rem;color:#142f40;background:#fafcfd}h1{font-size:1.8rem}section{margin:1.5rem 0}label{display:inline-block;margin:.4rem .8rem .4rem 0}select,input,button{font:inherit;padding:.35rem}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th,td{text-align:left;padding:.5rem;border-bottom:1px solid #cad7df}th button{border:0;background:transparent;font-weight:bold;cursor:pointer}p,li{line-height:1.5}.cards{display:flex;flex-wrap:wrap;gap:1rem}.card{padding:1rem;background:#e6f1f6;border-radius:.4rem;min-width:145px}.card strong{display:block;font-size:1.6rem}.scroll{overflow:auto;max-height:32rem}#trend td:last-child{width:65%}.bar{display:block;background:#237b9a;min-width:2px;height:1rem}small{color:#365365}.warn{background:#fff1ce;padding:1rem}a{color:#126184}</style>
<h1>OpenAntigens API and MCP traffic</h1>
<p>Application request counts include bots. Anonymous daily clients are approximate and do not establish human visitors.</p>
<p><a href="https://plausible.io/openantigens.org" target="_blank" rel="noreferrer noopener">Open existing Plausible dashboard</a> for eligible visitors, country, city and pages; select the <strong>API Request</strong> or <strong>MCP Request</strong> goal. This local report makes no network requests.</p>
<section aria-label="Report details"><p id="details"></p><p id="quality" class="warn"></p></section>
<section aria-label="Request filters"><h2>Explore requests</h2>
<label>Transport <select id="transport"></select></label><label>Species <select id="species"></select></label>
<label>Operation <select id="operation"></select></label><label>Status <select id="status"></select></label>
<label>First hour <select id="from"></select></label><label>Last hour <select id="to"></select></label>
<button id="reset">Reset filters</button><p><small>Filters apply to request cards, the time trend and operation rows. Coverage, daily clients and delivery totals describe the full imported interval.</small></p>
<p id="filter-coverage" aria-live="polite"></p><div class="cards" id="cards" aria-live="polite"></div></section>
<section><h2>Requests by hour</h2><div class="scroll"><table id="trend"><thead><tr><th>Hour</th><th>Requests</th><th>Relative count</th></tr></thead><tbody></tbody></table></div></section>
<section><h2>Request operations and outcomes</h2><p><small>Click a column heading to sort. Collapsed dimensions have unavailable status and admission information.</small></p><div class="scroll"><table id="operations"><thead><tr></tr></thead><tbody></tbody></table></div></section>
<section><h2>Anonymous daily clients and coverage</h2><p id="coverage"></p><div class="scroll"><table id="clients"><thead><tr><th>Day (UTC)</th><th>Approximate clients</th></tr></thead><tbody></tbody></table></div></section>
<section><h2>Plausible event delivery</h2><p>Accepted means the ingestion endpoint accepted the event. Plausible may still exclude it from visitor and location reports. Export covers eligible requests up to 100 attempts per process per UTC day; restarts and replicas can increase account usage.</p><table id="delivery"><thead><tr><th>Delivery outcome</th><th>Events</th></tr></thead><tbody></tbody></table><h3>Delivery exclusions and limits</h3><table id="policy"><thead><tr><th>Reason</th><th>Requests</th></tr></thead><tbody></tbody></table></section>
<section><h2>Reading the report</h2><ul>__NOTES__</ul></section>
<script type="application/json" id="report-data">__REPORT_DATA__</script>
<script>
'use strict';
const report=JSON.parse(document.getElementById('report-data').textContent), rows=report.requests;
const byId=id=>document.getElementById(id), fmt=n=>Number(n).toLocaleString('en-US');
const text=(id,value)=>{byId(id).textContent=value;};
text('details',`Time zone: ${report.parameters.timezone}. Filter: ${report.parameters.since_inclusive||'start of export'} (inclusive) to ${report.parameters.until_exclusive||'end of export'} (exclusive). Observed requests: ${report.observed.first_request_utc||'none'} to ${report.observed.last_request_utc||'none'}.`);
const q=report.input_quality;
text('quality',`Input quality: ${fmt(q.lines)} lines; ${fmt(q.ignored_lines)} unrelated, ${fmt(q.malformed_lines)} malformed, ${fmt(q.oversized_lines)} oversized; ${fmt(q.filtered_records)} records outside the interval; ${fmt(q.boundary_requests)} requests and ${fmt(q.boundary_deliveries)} deliveries excluded in boundary-spanning aggregates. Count and grouping limits: ${fmt(report.totals.bucket_overflow_requests)} requests in report overflow; ${fmt(report.totals.aggregate_requests)} requests in server aggregates.`);
const appendRow=(body,values)=>{const tr=document.createElement('tr');for(const value of values){const td=document.createElement('td');td.textContent=value;tr.append(td);}body.append(tr);return tr;};
for(const key of ['transport','species','operation','status','from','to']){const select=byId(key),field=['from','to'].includes(key)?'hour':key;const all=document.createElement('option');all.value='';all.textContent='All';select.append(all);for(const v of [...new Set(rows.filter(row=>!row.dimensions_collapsed).map(row=>String(row[field])))].sort()){const option=document.createElement('option');option.value=v;option.textContent=v;select.append(option);}select.addEventListener('change',update);}
const columns=['transport','operation','species','status','error_code','outcome','admitted','requests'];let sortKey='requests',direction=-1;
for(const column of columns){const th=document.createElement('th'),button=document.createElement('button');button.textContent=column.replaceAll('_',' ');button.addEventListener('click',()=>{direction=sortKey===column?-direction:-1;sortKey=column;update();});th.append(button);byId('operations').tHead.rows[0].append(th);}
function update(){const filtered=['transport','species','operation','status','from','to'].some(key=>byId(key).value);const selected=rows.filter(row=>(!filtered||!row.dimensions_collapsed)&&['transport','species','operation','status'].every(key=>!byId(key).value||String(row[key])===byId(key).value)&&(!byId('from').value||row.hour>=byId('from').value)&&(!byId('to').value||row.hour<=byId('to').value));const sum=key=>selected.reduce((n,row)=>n+row[key],0);const unavailable=report.totals.requests-report.totals.report_dimension_covered_requests;text('filter-coverage',filtered&&unavailable?`Partial filter coverage: ${fmt(unavailable)} imported requests have unavailable filter dimensions and are excluded from filtered cards. Clear filters to include their known outcome counts.`:`Known outcome counts include all selected records. ${fmt(unavailable)} imported requests have unavailable filter dimensions.`);const cards=[['Requests',sum('requests')],['HTTP errors / application errors',sum('errors')],['429 responses',sum('rate_limited')],['Rejected before admission',sum('rejected')],['Aborted',sum('aborted')],['Requests with unknown outcome counts',sum('requests')-sum('outcome_covered_requests')]];byId('cards').replaceChildren();for(const [label,value]of cards){const card=document.createElement('div'),strong=document.createElement('strong');card.className='card';strong.textContent=fmt(value);card.append(strong,document.createTextNode(label));byId('cards').append(card);}const hours=new Map(),groups=new Map();for(const row of selected){hours.set(row.hour,(hours.get(row.hour)||0)+row.requests);const values=columns.slice(0,-1).map(key=>row.dimensions_collapsed?'unavailable':row[key]),key=JSON.stringify(values);if(!groups.has(key))groups.set(key,{values,requests:0});groups.get(key).requests+=row.requests;}const body=byId('operations').tBodies[0];body.replaceChildren();const index=columns.indexOf(sortKey);for(const row of [...groups.values()].sort((a,b)=>{const av=index===7?a.requests:a.values[index],bv=index===7?b.requests:b.values[index];return direction*(typeof av==='number'&&typeof bv==='number'?av-bv:String(av).localeCompare(String(bv)));})){appendRow(body,[...row.values,fmt(row.requests)]);}const trend=byId('trend').tBodies[0];trend.replaceChildren();const max=Math.max(1,...hours.values());for(const [hour,n]of [...hours].sort((a,b)=>a[0].localeCompare(b[0]))){const tr=appendRow(trend,[hour,fmt(n),'']),bar=document.createElement('span');bar.className='bar';bar.style.width=`${100*n/max}%`;tr.lastChild.append(bar);}}
byId('reset').addEventListener('click',()=>{for(const id of ['transport','species','operation','status','from','to'])byId(id).value='';update();});
text('coverage',`${fmt(report.totals.client_covered_requests)} of ${fmt(report.totals.requests)} requests carry an anonymous daily client identity (${report.coverage.client_identity_percent}%). Source outcome coverage: ${report.coverage.request_dimensions_percent}%. Retained report grouping coverage: ${report.coverage.report_grouping_percent}%. Known outcome totals remain available when report grouping overflows. Sum of daily clients: ${fmt(report.coverage.daily_client_sum)}${report.coverage.client_counts_are_lower_bounds?' (lower bound: client capacity reached)':''}. Client counts cannot be filtered by request operation or added into unique people.`);
for(const row of report.daily_clients_utc)appendRow(byId('clients').tBodies[0],[row.day,fmt(row.approximate_clients)]);
for(const [outcome,n]of Object.entries(report.delivery_totals))appendRow(byId('delivery').tBodies[0],[outcome.replaceAll('_',' '),fmt(n)]);
const policies=new Map();for(const row of report.deliveries)if(row.reason!=='none')policies.set(row.reason,(policies.get(row.reason)||0)+row.count);for(const [reason,n]of [...policies].sort())appendRow(byId('policy').tBodies[0],[reason.replaceAll('_',' '),fmt(n)]);
update();
</script></html>'''


def write_reports(report: dict, output: Path) -> None:
    """Publish three new private files, preserving existing files and linked targets."""
    names = ("report.json", "report.csv", "report.html")
    # Walk via directory descriptors so a symlink in any parent cannot redirect writes.
    directory = os.open("/", os.O_RDONLY | os.O_DIRECTORY)
    staged: dict[str, tuple[str, tuple[int, int]]] = {}
    published: list[str] = []
    try:
        for part in Path(os.path.abspath(output)).parts[1:]:
            try:
                child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory)
            except FileNotFoundError:
                os.mkdir(part, mode=0o700, dir_fd=directory)
                child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory)
            os.close(directory)
            directory = child
        for name in names:
            try:
                os.stat(name, dir_fd=directory, follow_symlinks=False)
            except FileNotFoundError:
                continue
            raise FileExistsError("report destination exists; use a new output directory")
        for name in names:
            temporary = f".traffic-{secrets.token_hex(16)}"
            descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                                 0o600, dir_fd=directory)
            info = os.fstat(descriptor)
            staged[name] = (temporary, (info.st_dev, info.st_ino))
            with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
                os.fchmod(stream.fileno(), 0o600)
                if name == "report.json":
                    json.dump(report, stream, indent=2, sort_keys=True)
                    stream.write("\n")
                elif name == "report.csv":
                    csv.writer(stream).writerows(tuple(csv_safe(value) for value in row) for row in csv_rows(report))
                else:
                    stream.write(render_html(report))
                stream.flush()
                os.fsync(stream.fileno())
        for name, (temporary, _) in staged.items():
            # A hard link publishes the completed file atomically and refuses any existing name.
            os.link(temporary, name, src_dir_fd=directory, dst_dir_fd=directory, follow_symlinks=False)
            published.append(name)
    except BaseException:
        for name in published:
            try:
                info = os.stat(name, dir_fd=directory, follow_symlinks=False)
                if (info.st_dev, info.st_ino) == staged[name][1]:
                    os.unlink(name, dir_fd=directory)
            except FileNotFoundError:
                pass
        raise
    finally:
        for temporary, _ in staged.values():
            try:
                os.unlink(temporary, dir_fd=directory)
            except FileNotFoundError:
                pass
        os.close(directory)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", type=Path, nargs="+", help="Non-overlapping exported application JSONL/text logs")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--since", help="Inclusive ISO timestamp with UTC offset")
    parser.add_argument("--until", help="Exclusive ISO timestamp with UTC offset")
    parser.add_argument("--timezone", default="UTC", help="IANA time zone for hourly trends (default UTC); client days stay UTC")
    args = parser.parse_args(argv)
    try:
        since = instant(args.since) if args.since else None
        until = instant(args.until) if args.until else None
        if since and until and since >= until:
            raise ValueError("since must precede until")
        if len({path.resolve() for path in args.logs}) != len(args.logs):
            raise ValueError("each input file must be supplied once")
        sources = {path.resolve() for path in args.logs}
        if any((args.output_dir / name).resolve() in sources for name in ("report.json", "report.csv", "report.html")):
            raise ValueError("report output would replace an input")
        report = summarize(args.logs, since=since, until=until, tz=args.timezone)
        write_reports(report, args.output_dir)
    except (OSError, ValueError, ZoneInfoNotFoundError):
        # Source file paths and input contents can themselves contain private data.
        parser.error("cannot generate report: check readable inputs, unique files, output directory, timestamp offsets and IANA time zone")
    print("Wrote report.json, report.csv and report.html locally.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

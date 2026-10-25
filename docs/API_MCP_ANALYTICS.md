# API and MCP traffic reports

Use the existing [Plausible dashboard](https://plausible.io/openantigens.org) for eligible visitors, countries, cities and pages. Use the local report below for application request counts, including bots, errors and requests refused before admission. Both human and mouse traffic use the same records and report controls.

Plausible filters some automated traffic. An accepted event can still be absent from its visitor or location reports. For API/MCP calls, locations describe the verified network client, which may be an agent or gateway. The local report keeps those two measures separate. It reports approximate anonymous daily clients when an identity is available; it does not infer people, country or city from those identities. The application ledger is on by default; Plausible delivery is off by default.

## Export production logs

1. Sign in to the existing GoDaddy account and open **Node.js Hosting**, then **My Apps**. Select the OpenAntigens service and open **Runtime Logs**. These controls are described in [GoDaddy's Node.js Hosting help](https://www.godaddy.com/help/deploy-my-loveable-app-with-godaddy-nodejs-hosting-42981).
2. Select production (**PUBLISH**) rather than the preview app. Select standard output (**STDOUT**) and the complete required time range, with all severity levels. Keep preview and production exports separate. GoDaddy documents a default lookback of 24 hours and a maximum lookback of seven days in [How to read logs](https://developer.godaddy.com/en/docs/api-users/hosting/read-logs). Archive exports locally before that window expires; this report cannot recover missing logs.
3. Save the complete log messages as a UTF-8 text file. If the current panel has no download control, copy its complete log messages into a plain text file. Preserve one JSON message per line. Timestamp or severity prefixes before each JSON object are accepted. Paginate or load the full selected interval before copying; a visible screen of recent messages may be only part of the export.
4. Use non-overlapping files from the same environment. Repeated or overlapping exports count the same traffic again. Keep original exports private: host logs can contain fields beyond the application's sanitized records.

The signed-in log panel's current export controls and completeness must be checked by the operator. The public documentation establishes the log window and the Runtime Logs entry point, but does not establish a particular download button. This workflow needs no new account, API key, hosting change or report server.

## Generate and open a report

Python 3.10 or later is sufficient. Run from the repository root:

```sh
python3 scripts/api_traffic_report.py /private/path/production-stdout.log \
  --output-dir /private/path/traffic-report
```

Open `/private/path/traffic-report/report.html` in a browser. The page is self-contained and reads its embedded aggregate data. It makes no network requests. The Plausible link opens the existing dashboard only when selected. Filter by transport, species, operation, status or hour; click operation-table headings to sort. Filtered cards show requests, errors, 429 responses, refusals and aborted requests. Daily client counts and Plausible delivery totals describe the full imported interval.

The same run writes `report.json` and `report.csv`. All three formats contain the same counts. CSV uses one metric per row, with a section and JSON dimensions, so it can be grouped in a spreadsheet. CSV cells are protected against formula interpretation.

To select an interval, give timestamps with explicit UTC offsets. The start is inclusive and the end is exclusive. Hourly trends use UTC by default. Use the existing dashboard's time zone when comparing periods:

```sh
python3 scripts/api_traffic_report.py /private/path/production-stdout.log \
  --since 2026-10-24T00:00:00-04:00 \
  --until 2026-10-25T00:00:00-04:00 \
  --timezone America/New_York \
  --output-dir /private/path/traffic-report-2026-10-24
```

Several non-overlapping files can precede `--output-dir`. Use a new output directory for each report. The writer checks all three destinations before staging files and refuses existing files, symlinks and hardlinks. It also rejects symlinks in the output directory and its parents; use physical directory paths, such as `/private/tmp` rather than `/tmp` on macOS. Report files have private `0600` permissions and appear only after their contents have been written. A staging or publication failure removes files created by that run while preserving existing material. Reports retain normalized interval parameters, bounds and input-quality counters; they omit input filenames and paths.

## Read counts and coverage

Total requests include each validated individual request once and each aggregate's `count` once. They include requests that finish, abort or fail admission. API and MCP operation names are controlled categories; requested targets, search text and tool arguments are omitted. Status, error code and outcome distinguish HTTP errors, application errors, 429 responses and aborted requests. A refusal means `admitted: false`; it includes rejection before the request enters the handler.

Anonymous daily clients use a secret-keyed hash of the UTC date, verified IP address and safe user agent. Clients with the same address and user agent share an identity; a different address or user agent creates another identity. The service emits no identity when the required secret or verified client information is unavailable. Daily totals always use UTC, even when hourly trends use another time zone. Their sum describes daily clients, not unique visitors across days. Client coverage is the fraction of requests with an available identity. The report omits the hashes themselves.

During output backpressure, the service merges requests with the same dimensions into bounded aggregate records. These aggregates retain request counts, duration and bytes but carry no client identities. If server dimension capacity is exhausted, a collapsed aggregate retains counts while its dimensions are unavailable. The report exposes dimension coverage and does not interpret its placeholder status or admission value as an observed outcome. A time filter excludes any aggregate that crosses its boundary and reports the excluded request and delivery counts. An aggregate that spans several hours appears in its last request hour.

Input is streamed with a 64 KiB line limit. Oversized and malformed lines are counted and excluded; unrelated host messages are ignored. Request and delivery tables each retain at most 20,000 ordinary groups plus fixed overflow groups. Known outcome counts remain in overflow when that limit is reached: unfiltered request, error, 429, refusal and abort cards retain the same counts as JSON and CSV. Source outcome coverage and retained report grouping coverage are shown separately. Filters exclude groups whose dimensions are unavailable, with an explicit partial-coverage notice. The report retains at most 100,000 distinct daily client identities. It labels client totals as lower bounds if that capacity is reached. Source fields outside the strict contract never enter a report.

Complete counts depend on complete exports and successful host log capture. A process termination can lose buffered aggregates. Missing exports, malformed records, host truncation and boundary exclusions reduce coverage. Check the report's input-quality notice before comparing intervals. No report can establish that an absent request never reached an upstream proxy or hosting edge.

Plausible delivery outcomes are independent of server request counts:

| Outcome | Meaning |
| --- | --- |
| `accepted` | The ingestion endpoint accepted the event; dashboard inclusion is still subject to Plausible's processing and bot rules. |
| `provider_dropped` | Plausible explicitly indicated that it dropped the event. |
| `failed` | Delivery returned an error or failed. |
| `timeout` | The bounded delivery timeout expired. |
| `capacity_skipped` | Concurrent delivery capacity was full, or the process exhausted its delivery budget. The `daily_budget` reason identifies the latter. |
| `identity_unavailable` | The verified network identity, user agent or transport required for delivery was unavailable. |
| `disabled` | Server-side Plausible delivery was off. The local request record still counts. |
| `policy_skipped` | The request was outside the permitted Plausible export policy. The reason identifies refusal, an excluded operation, request failure or an unfinished request. The local request record still counts. |

## Prepare optional Plausible goals

The existing dashboard is the destination. The **API Request** and **MCP Request** custom-event goals require operator setup; their presence and active delivery must be verified. After approval to change analytics settings:

1. Check the existing Plausible plan's event allowance and expected request volume. Server request events add usage. Delivery is restricted to admitted, finished data lookup operations: successful requests and specified lookup failures, such as a missing target or unavailable release. Health checks, protocol handshakes, refusals, 429 responses and aborted requests remain in the local report. With delivery enabled, excluded requests have a `policy_skipped` outcome and a controlled reason. With delivery off, every request has a `disabled` outcome. Enabled export has a fixed maximum of 100 attempts per process per UTC day, so Plausible covers only part of eligible server traffic once that budget is exhausted. Restarts reset the process budget; replicas each have their own budget. This limit does not establish an account-wide daily or monthly cap. Verify restart behavior, replica count and the actual plan allowance before enabling delivery.
2. In the existing site's settings, create custom-event goals named exactly **API Request** and **MCP Request**. See [Plausible's custom-event documentation](https://plausible.io/docs/custom-event-goals). Do not create a second analytics site.
3. Complete the trusted-client and proxy checks in [release operations](API_MCP_DEPLOYMENT.md). Enable server-side delivery only in the approved deployment configuration. The service's local client secret and any host credentials stay outside reports and browser code. This report requires no Plausible API key.
4. After the approved deployment, compare the local delivery outcomes with the dashboard goals. Check both human and mouse events. Use eligible dashboard traffic for country, city and visitor views, and the exported server records for total API/MCP counts. Avoid repeated or bulk requests to manufacture analytics traffic.

The application ledger contract and bounded delivery behavior are implemented in `services/public-api/src/analytics.js`. The local report validates that contract independently. Run its regression checks with:

```sh
python3 -m unittest discover -s tests -p test_api_traffic_report.py -v
```

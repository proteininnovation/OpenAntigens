# API/MCP host verification — 2026-10-25

## Public checks

Read-only checks confirmed that public DNS resolves for the existing portal, the public edge presents valid TLS, and both human and mouse portal metadata endpoints return HTTP 200. DNS for the planned API subdomain did not resolve in the October 25 check.

The authenticated hosting dashboard confirmed that Node.js Apps are included with the account's plan as of October 25. No app was created. GoDaddy documentation describes Node 22 and ZIP uploads under 100 MB; its Node.js Hosting FAQ labels managed Node hosting Beta and states that preview apps are removed after 30 days.

## Previous all-data package measurement

The previous all-data package failed the managed Node upload gate. After full source, static-copy and transport validation, the measured ZIP was **245,426,866 bytes**, above both the local 80,000,000-byte budget and GoDaddy's documented hard limit below 100,000,000 bytes. The packager refused to publish the ZIP or its artifact record. No upload was attempted.

The separate local runtime diagnostic contains all 9,949 targets and 78,125 constructs. Its physical data occupies 316,899,235 bytes in 9,963 files, including metadata and compressed target containers. Node 22.23.4 reached readiness in 59.901 seconds in maintenance mode and 60.929 seconds in active mode. Maximum observed peak RSS was 1,213,923,328 bytes. Both modes passed all five tools through HTTP and MCP for selected human and mouse queries, with exact sequence hashes and matching results. These are local macOS measurements, not managed-app allowances; the cPanel memory display does not establish the Node app's budget.

## Fixed-GoDaddy revision

The owner confirmed GoDaddy as the provider. The revised plan keeps sealed metadata, indexes and exact full-summary lookup files in the Node app, with immutable target containers on the existing GoDaddy static site. Search makes zero target-origin requests. Target reads use fixed-origin, hash-qualified paths with strict byte, concurrency, cache and deadline bounds. Scientific IDs, bytes, schemas and the approved top-100 disease scope remain unchanged.

The actual current-plus-previous artifact measured **12,909,859 bytes ZIP / 82,198,385 bytes expanded**, below the unchanged 80 MB reserved upload budget and documented 100 MB limit. Each release contains 9,949 targets and 78,125 constructs. Their 19,898 static target containers total 447,159,606 bytes. The final code repack preserves all 18 data files and seals from the successful full-validation baseline; current policy/static checks and inherited scientific validation are labelled separately in its receipt.

Local Node 22.23.4 profiles passed all five HTTP/MCP tools for both species and releases in active and maintenance modes. Reported app startup was 436.5 ms and 434.9 ms. Peak RSS was 257,327,104 bytes across the core runs and 348,897,280 bytes for the largest-target diagnostic. Startup, metadata and search fetched zero target assets; cold target reads fetched one, immediate warm reads none. All 87 paced API attempts completed with no retry or failed check. The fixture served nine asset reads, with no external requests or live Plausible delivery.

The local static tree after the recorded header refresh contains 70,054 regular files and 17,539,445,070 bytes, excluding other hosting-account data. These are local capacity measurements, not account-wide free-space evidence. The old all-data startup peak does not describe this reader. GoDaddy remains fixed; no different provider or paid change is authorized.

Profiles used a strict loopback fixture for the sealed static URLs, so they do not verify GoDaddy TLS/WAF, actual wire headers, outbound routing or production latency. Genuine two-client concurrency could not be measured because the OS would not bind a second source address. External `ps` sampling was unavailable; Node RSS samples and its OS peak counter were available. These limits remain explicit before launch.

## Launch gates

The live managed Node runtime, app-root behavior, app-specific resource limits, application isolation, subscription and publishing terms, and any related costs remain unverified. Confirm API subdomain routing, proxy and MCP transport behavior, and outbound access to the sealed static assets. Verify TLS from the edge to the origin and actual `.oap` MIME, absent HTTP content encoding, WAF and cache behavior. Measure app and static storage separately, including retained/rollback asset unions and other account data. Verify additive publication and rollback before launch. cPanel account displays do not establish managed Node app budgets. These are checks within the fixed GoDaddy decision.

These checks were read-only. No service app was created or deployed. The full account-specific verification is retained in the ignored local evidence folder `outputs/host-verification-2026-10-25/` and is not part of the public documentation.

Source: [GoDaddy Node.js Hosting FAQ](https://www.godaddy.com/help/godaddy-nodejs-hosting-faq-42915); provider runtime and upload documentation summarized in [API_MCP_PLAN.md](API_MCP_PLAN.md) §10.

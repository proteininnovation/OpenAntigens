# API/MCP host verification — 2026-10-25

## Public checks

Read-only checks confirmed that public DNS resolves for the existing portal, the public edge presents valid TLS, and both human and mouse portal metadata endpoints return HTTP 200. DNS for the planned API subdomain did not resolve in the October 25 check.

The authenticated hosting dashboard confirmed that Node.js Apps are included with the account's plan as of October 25. No app was created. GoDaddy documentation describes Node 22 and ZIP uploads under 100 MB; its Node.js Hosting FAQ labels managed Node hosting Beta and states that preview apps are removed after 30 days.

## Full-catalogue capacity result

The current self-contained catalogue fails the managed Node upload gate. After full source, static-copy and transport validation, the measured ZIP was **245,426,866 bytes**, above both the local 80,000,000-byte budget and GoDaddy's documented hard limit below 100,000,000 bytes. The packager refused to publish the ZIP or its artifact record. No upload was attempted.

The separate local runtime diagnostic contains all 9,949 targets and 78,125 constructs. Its physical data occupies 316,899,235 bytes in 9,963 files, including metadata and compressed target containers. Node 22.23.4 reached readiness in 59.901 seconds in maintenance mode and 60.929 seconds in active mode. Maximum observed peak RSS was 1,213,923,328 bytes. Both modes passed all five tools through HTTP and MCP for selected human and mouse queries, with exact sequence hashes and matching results. These are local macOS measurements, not managed-app allowances; the cPanel memory display does not establish the Node app's budget.

The included app slot therefore does not establish that this complete release can run there. Launch remains blocked. The recommended next deployment decision is a Node 22 host that accepts the complete artifact and provides startup time and memory above the measured peaks with headroom. A new provider/VPS, spending or a different data-hosting design requires separate authorization under the plan.

## Launch gates

The live managed Node runtime, app-root behavior, app-specific resource limits, application isolation, subscription and publishing terms, and any related costs remain unverified. Confirm API subdomain routing, proxy and MCP transport behavior. Verify TLS from the edge to the origin. Measure the complete release package against current host limits and verify rollback before publication. cPanel account displays do not establish managed Node app budgets.

These checks were read-only. No service app was created or deployed. The full account-specific verification is retained in the ignored local evidence folder `outputs/host-verification-2026-10-25/` and is not part of the public documentation.

Source: [GoDaddy Node.js Hosting FAQ](https://www.godaddy.com/help/godaddy-nodejs-hosting-faq-42915); provider runtime and upload documentation summarized in [API_MCP_PLAN.md](API_MCP_PLAN.md) §10.

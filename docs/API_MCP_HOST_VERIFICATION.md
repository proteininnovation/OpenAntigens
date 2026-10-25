# API/MCP host verification — 2026-10-25

## Public checks

Read-only checks confirmed that public DNS resolves for the existing portal, the public edge presents valid TLS, and both human and mouse portal metadata endpoints return HTTP 200. DNS for the planned API subdomain did not resolve in the October 25 check.

The authenticated hosting dashboard confirmed that Node.js Apps are included with the account's plan as of October 25. No app was created. GoDaddy documentation describes Node 22 and ZIP uploads under 100 MB; its Node.js Hosting FAQ labels managed Node hosting Beta and states that preview apps are removed after 30 days.

## Launch gates

The live managed Node runtime, app-root behavior, app-specific resource limits, application isolation, subscription and publishing terms, and any related costs remain unverified. Confirm API subdomain routing, proxy and MCP transport behavior. Verify TLS from the edge to the origin. Measure the complete release package against current host limits and verify rollback before publication. cPanel account displays do not establish managed Node app budgets.

These checks were read-only. No service app was created or deployed. The full account-specific verification is retained in the ignored local evidence folder `outputs/host-verification-2026-10-25/` and is not part of the public documentation.

Source: [GoDaddy Node.js Hosting FAQ](https://www.godaddy.com/help/godaddy-nodejs-hosting-faq-42915); provider runtime and upload documentation summarized in [API_MCP_PLAN.md](API_MCP_PLAN.md) §10.

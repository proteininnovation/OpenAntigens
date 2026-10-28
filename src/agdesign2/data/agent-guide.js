const agentCopy = document.getElementById('agent-copy');
agentCopy.addEventListener('click', async () => {
  const text = document.getElementById('agent-prompt-text');
  const status = document.getElementById('agent-copy-status');
  try {
    await navigator.clipboard.writeText(text.value);
    status.textContent = 'Prompt copied. Paste it into your assistant and add your question.';
  } catch {
    document.getElementById('agent-prompt').open = true;
    text.focus();
    text.select();
    status.textContent = 'Clipboard access is unavailable. Copy the selected prompt manually.';
  }
});

const connection = document.getElementById('mcp-connection');
if (connection) {
  const status = document.getElementById('mcp-connection-status');
  const address = document.getElementById('mcp-server-url');
  const copyStatus = document.getElementById('mcp-copy-status');
  document.getElementById('mcp-copy').addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(address.value);
      copyStatus.textContent = 'Server address copied. Paste it into your MCP client.';
    } catch {
      address.focus();
      address.select();
      copyStatus.textContent = 'Clipboard access is unavailable. Copy the selected address manually.';
    }
  });
  (async () => {
    try {
      const response = await fetch('portal_metadata.json', {
        cache: 'no-store', credentials: 'omit', redirect: 'error', signal: AbortSignal.timeout(5000),
      });
      if (!response.ok) throw new Error('Portal metadata unavailable');
      const { public_api: api } = await response.json();
      if (api?.state !== 'active' || typeof api.mcp_url !== 'string') {
        status.textContent = 'MCP is not published for this portal yet. Use Copy AI prompt to work with the website and downloads.';
        return;
      }
      const url = new URL(api.mcp_url);
      if (url.protocol !== 'https:' || url.username || url.password || url.search || url.hash ||
          url.pathname !== '/mcp') throw new Error('Invalid MCP address');
      address.value = url.href;
      document.getElementById('mcp-connection-controls').hidden = false;
      status.textContent = 'Published MCP server. Use this address to connect your AI client.';
    } catch {
      status.textContent = 'The MCP connection address could not be loaded. Use Copy AI prompt to work with the website and downloads, or try this page again.';
    }
  })();
}

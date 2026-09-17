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

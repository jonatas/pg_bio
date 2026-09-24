# IDE and Assistant Integration

The true power of `pg_bio` extends beyond raw Python scripts. We have bundled `scripts/mcp_server.py`, an open-standard **Model Context Protocol (MCP)** server. Because it strictly adheres to the MCP specification, any modern assistant can mount this server and gain the ability to fold proteins natively on your local hardware.

You are not locked into any specific ecosystem. You can plug this bio-engine into Claude, Cursor, Zed, or any other MCP-compliant client.

## 1. Connecting to Claude Desktop (Anthropic)
To give Claude the ability to interact with your local database and run predictions, edit your Claude Desktop configuration file (usually located at `~/Library/Application Support/Claude/claude_desktop_config.json` on macOS):

```json
{
  "mcpServers": {
    "pg_bio_server": {
      "command": "uv",
      "args": [
        "run",
        "/absolute/path/to/pg_bio/scripts/mcp_server.py"
      ]
    }
  }
}
```
Restart Claude Desktop. You can now ask Claude: *"Please generate a novel enzyme sequence and fold it for me."* Claude will invisibly trigger your local database via the MCP tool.

## 2. Connecting to Cursor IDE
1. Open Cursor Settings > Features > MCP.
2. Click **+ Add New MCP Server**.
3. Type: `stdio`
4. Command: `uv run /absolute/path/to/pg_bio/scripts/mcp_server.py`
5. Cursor's internal assistant will now have the `fold_sequence` tool available for autonomous biological engineering.

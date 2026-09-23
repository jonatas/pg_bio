# pg_bio: The PostgreSQL Bioinformatics OS

`pg_bio` is a high-performance PostgreSQL extension and local AI ecosystem designed for synthetic biology. It natively handles 3D Z-Order spatial indexing, high-dimensional vector embeddings, and sparse neural network attention maps.

## Model Context Protocol (MCP) Integration
The true power of `pg_bio` is its **Vendor-Agnostic AI Integration**. 

We have bundled `mcp_server.py`, an open-standard Model Context Protocol (MCP) server. Because it strictly adheres to the MCP specification, **any modern AI assistant can mount this server and gain the ability to fold proteins natively on your local GPU.**

You are not locked into Antigravity. You can plug this bio-engine into Claude, Cursor, Zed, or any other MCP-compliant client.

### 1. Connecting to Claude Desktop (Anthropic)
To give Claude the ability to fold proteins locally, edit your Claude Desktop configuration file (usually located at `~/Library/Application Support/Claude/claude_desktop_config.json` on macOS):

```json
{
  "mcpServers": {
    "pg_bio_esmfold": {
      "command": "uv",
      "args": [
        "run",
        "/absolute/path/to/pg_bio/mcp_server.py"
      ]
    }
  }
}
```
Restart Claude Desktop. You can now ask Claude: *"Please generate a novel enzyme sequence and fold it for me."* Claude will invisibly trigger your local Mac GPU to fold the protein via the MCP tool.

### 2. Connecting to Cursor IDE
1. Open Cursor Settings > Features > MCP.
2. Click **+ Add New MCP Server**.
3. Type: `stdio`
4. Command: `uv run /absolute/path/to/pg_bio/mcp_server.py`
5. Cursor's internal AI (whether you use Claude 3.5 Sonnet, GPT-4o, or Codex) will now have the `fold_sequence` tool available for autonomous biological engineering.

### 3. Connecting to Antigravity
Antigravity automatically discovers this server via the `.agents/mcp_config.json` file we already configured in this repository. No extra setup is required.

---

## Extension Architecture
* **`src/lib.rs`**: Core Rust implementation of the `pg_bio` Postgres extension (Z-Order, Vector Distances, Sparse Attention Maps).
* **`experiment_*.sql`**: SQL scripts demonstrating end-to-end biological computations natively in the database.
* **`benchmark.py`**: Python benchmarking proving that Postgres Z-Order spatial queries outperform standard Python/Numpy math by 4.5x.

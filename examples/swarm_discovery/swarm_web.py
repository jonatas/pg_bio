"""Forwarding shim: Unified Pipeline System Web UI & Telemetry moved to jonatas/bio-pipelines.

Canonical repository & entrypoint:
    https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_web.py
"""
import os
import sys

BIO_PIPELINES_REPO = "https://github.com/jonatas/bio-pipelines"
BIO_PIPELINES_ROOT = os.getenv(
    "BIO_PIPELINES_ROOT",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "bio-pipelines")),
)

if not os.path.isdir(BIO_PIPELINES_ROOT):
    raise RuntimeError(
        f"Unified Pipeline System Web UI has moved to {BIO_PIPELINES_REPO}. "
        f"Clone it alongside pg_bio or set BIO_PIPELINES_ROOT to your local checkout."
    )

if BIO_PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, BIO_PIPELINES_ROOT)
    sys.path.insert(0, os.path.join(BIO_PIPELINES_ROOT, "base"))

from base.swarm_web import app  # noqa: E402,F401

if __name__ == "__main__":
    port = int(os.getenv("SWARM_WEB_PORT", "8080"))
    app.run(host="127.0.0.1", port=port, debug=False)

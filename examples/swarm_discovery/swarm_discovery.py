"""Forwarding shim: Swarm Discovery & Autonomous Intelligence moved to jonatas/bio-pipelines.

Canonical repository & entrypoint:
    https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_discovery.py
"""
import os
import sys
import runpy

BIO_PIPELINES_REPO = "https://github.com/jonatas/bio-pipelines"
BIO_PIPELINES_ROOT = os.getenv(
    "BIO_PIPELINES_ROOT",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "bio-pipelines")),
)
TARGET_SCRIPT = os.path.join(BIO_PIPELINES_ROOT, "base", "swarm_discovery.py")

if not os.path.isdir(BIO_PIPELINES_ROOT):
    raise RuntimeError(
        f"Swarm Discovery has moved to {BIO_PIPELINES_REPO}. "
        f"Clone it alongside pg_bio or set BIO_PIPELINES_ROOT to your local checkout."
    )

if BIO_PIPELINES_ROOT not in sys.path:
    sys.path.insert(0, BIO_PIPELINES_ROOT)
    sys.path.insert(0, os.path.join(BIO_PIPELINES_ROOT, "base"))

if __name__ == "__main__":
    os.chdir(BIO_PIPELINES_ROOT)
    runpy.run_path(TARGET_SCRIPT, run_name="__main__")
else:
    from base.swarm_discovery import *  # noqa: F401,F403

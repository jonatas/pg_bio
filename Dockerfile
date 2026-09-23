FROM postgres:16-bookworm

# Install build dependencies
RUN apt-get update && apt-get install -y \
    curl \
    build-essential \
    libclang-dev \
    postgresql-server-dev-16 \
    pkg-config \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install cargo-pgrx (matching the version in Cargo.toml)
RUN cargo install --locked cargo-pgrx --version 0.17.0

# Initialize pgrx pointing to the system Postgres 16
RUN cargo pgrx init --pg16 /usr/lib/postgresql/16/bin/pg_config

WORKDIR /usr/src/pg_bio
COPY . .

# Build and install the extension into the system postgres
RUN cargo pgrx install --pg16 --release

# Ensure the extension is loaded on startup (optional but helpful)
# We can just let users CREATE EXTENSION pg_bio;

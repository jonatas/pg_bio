FROM postgres:18-bookworm

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    build-essential \
    libclang-dev \
    clang \
    libssl-dev \
    postgresql-server-dev-18 \
    postgresql-18-pgvector \
    pkg-config \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"

# Install cargo-pgrx (matching the version in Cargo.toml)
RUN cargo install --locked cargo-pgrx --version 0.17.0

# Initialize pgrx pointing to the system Postgres 18
RUN cargo pgrx init --pg18 /usr/lib/postgresql/18/bin/pg_config

WORKDIR /usr/src/pg_bio
COPY . .

# Build and install the extension into the system Postgres 18
RUN cargo pgrx install --pg18 --release --no-default-features --features pg18

FROM postgres:18-bookworm

# Install build dependencies (including meson, ninja-build, and libopenblas-dev for pg_vectorsearch)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    build-essential \
    libclang-dev \
    clang \
    libssl-dev \
    libopenblas-dev \
    meson \
    ninja-build \
    postgresql-server-dev-18 \
    postgresql-18-pgvector \
    pkg-config \
    git \
    && rm -rf /var/lib/apt/lists/*

# Build and install Tiger Data's pg_vectorsearch (PRISM + RaBitQ) for PostgreSQL 18
RUN git clone --depth 1 https://github.com/timescale/pg_vectorsearch.git /usr/src/pg_vectorsearch \
    && make -C /usr/src/pg_vectorsearch install PG_CONFIG=/usr/lib/postgresql/18/bin/pg_config \
    && rm -rf /usr/src/pg_vectorsearch

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

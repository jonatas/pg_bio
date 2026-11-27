{
  lib,
  buildPgrxExtension,
  cargo-pgrx_0_17_0,
  postgresql,
}:

(buildPgrxExtension {
  pname = "pg_bio";
  inherit ((lib.importTOML ../Cargo.toml).package) version;

  # Only what the build reads: no videos/, data/, notebooks or docs.
  src = lib.fileset.toSource {
    root = ./..;
    fileset = lib.fileset.unions [
      ../Cargo.toml
      ../Cargo.lock
      ../pg_bio.control
      ../.cargo
      ../src
    ];
  };

  cargoLock.lockFile = ../Cargo.lock;

  inherit postgresql;
  # Must equal the pgrx version Cargo.toml pins (=0.17.0): bump both
  # together, or the build stops on a cargo-pgrx version mismatch.
  cargo-pgrx = cargo-pgrx_0_17_0;

  # The default feature is pg13; build for the postgresql given, only.
  buildNoDefaultFeatures = true;
  buildFeatures = [ "pg${lib.versions.major postgresql.version}" ];

  # `cargo pgrx test pg18` needs a running server inside the build; the
  # flake's pg-bio-load check loads the extension instead.
  doCheck = false;

  meta = {
    description = "Bioinformatics primitives for PostgreSQL";
    homepage = "https://github.com/jonatas/pg_bio";
    platforms = postgresql.meta.platforms;
  };
}).overrideAttrs
  (old: {
    # buildPgrxExtension marks every pgrx extension broken on macOS, because
    # PostgreSQL started inside the macOS sandbox can leak SysV shared memory
    # (NixOS/nix#12548). The build itself starts no server here (doCheck is
    # off); the flake's pg-bio-load check does, and stops it on every exit.
    meta = old.meta // {
      broken = false;
    };
  })

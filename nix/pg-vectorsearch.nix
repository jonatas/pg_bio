{
  lib,
  stdenv,
  fetchFromGitHub,
  meson,
  ninja,
  pkg-config,
  openblas,
  postgresql,
  postgresqlBuildExtension,
}:

# Tiger Data's pg_vectorsearch, which pg_bio.control requires. Remove this
# file once nixpkgs ships postgresqlPackages.pg_vectorsearch.
postgresqlBuildExtension (finalAttrs: {
  pname = "pg_vectorsearch";
  version = "0.2.0";

  src = fetchFromGitHub {
    owner = "timescale";
    repo = "pg_vectorsearch";
    tag = "v${finalAttrs.version}";
    hash = "sha256-l2j1oYmSEvmEggj5JqfgzfJ5743Ok/ivhn1qX1X8ABE=";
  };

  # macOS: meson keeps -isysroot from pg_config but drops its SDK path, and
  # nothing compiles ('math.h' file not found). Fix offered upstream.
  patches = [ ./pg-vectorsearch-darwin-isysroot.patch ];

  nativeBuildInputs = [
    meson
    ninja
    pkg-config
  ];

  # BLAS speeds up RaBitQ encoding: meson finds Accelerate on macOS.
  buildInputs = lib.optionals (!stdenv.hostPlatform.isDarwin) [ openblas ];

  mesonBuildType = "release";
  mesonFlags = [
    "-Dpostgresql=enabled"
    "-Dpg_config=${lib.getExe' postgresql.pg_config "pg_config"}"
    "-Dtools=disabled"
    "-Dtap_tests=disabled"
    "-Dhdf5=disabled"
    "-Dfaiss=false"
    "-Dnative=false"
  ];

  # meson installs into pg_config's directories, inside the postgresql
  # store path; DESTDIR lands them under $out, and postgresqlBuildExtension
  # moves them into place.
  preInstall = ''
    export DESTDIR=$out
  '';

  meta = {
    description = "Vector search for PostgreSQL with PRISM and RaBitQ";
    homepage = "https://github.com/timescale/pg_vectorsearch";
    license = lib.licenses.postgresql;
    platforms = postgresql.meta.platforms;
  };
})

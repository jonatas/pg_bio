{
  description = "pg_bio: bioinformatics primitives for PostgreSQL 18, as a nix flake";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-26.05";

  outputs =
    { self, nixpkgs }:
    let
      systems = [
        "aarch64-darwin"
        "x86_64-darwin"
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAll = f: nixpkgs.lib.genAttrs systems (system: f (pkgsFor system));
      pkgsFor =
        system:
        import nixpkgs {
          inherit system;
          overlays = [ self.overlays.default ];
        };
    in
    {
      # Both build in postgresql_18's own scope, which supplies postgresql,
      # its stdenv (clang on Linux, which pg_vectorsearch's meson requires)
      # and postgresqlBuildExtension. The PostgreSQL major is chosen here
      # only.
      overlays.default =
        final: _prev:
        let
          pg = final.postgresql_18;
        in
        {
          pg_vectorsearch = pg.pkgs.callPackage ./nix/pg-vectorsearch.nix { };
          pg_bio = pg.pkgs.callPackage ./nix/pg-bio.nix { };
          # PostgreSQL 18 with everything pg_bio.control requires.
          postgresql-pg_bio = pg.withPackages (p: [
            p.pgvector
            final.pg_vectorsearch
            final.pg_bio
          ]);
        };

      packages = forAll (pkgs: {
        inherit (pkgs) pg_bio pg_vectorsearch postgresql-pg_bio;
        default = pkgs.postgresql-pg_bio;
      });

      # CREATE EXTENSION pg_bio CASCADE in a throwaway cluster, offline.
      checks = forAll (pkgs: {
        pg-bio-load = pkgs.runCommand "pg-bio-load" {
          nativeBuildInputs = [ pkgs.postgresql-pg_bio ];
        } (builtins.readFile ./nix/pg-bio-load.sh);
      });
    };
}

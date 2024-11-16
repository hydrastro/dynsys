{
  description = "dynsys";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    ds.url = "github:hydrastro/ds";
    lizard.url = "github:hydrastro/lizard";
  };

  outputs = { self, nixpkgs, ds, lizard }: {
    packages = nixpkgs.lib.genAttrs [ "x86_64-linux" ] (system:
      let
        pkgs = import nixpkgs { inherit system; };
      in rec {
        dynsys = pkgs.stdenv.mkDerivation {
          pname = "dynsys";
          version = "0.0.0";

          src = ./.;

          # (pkgs.callPackage ds { })
          buildInputs = [ pkgs.gmp pkgs.stdenv.cc ds.defaultPackage.x86_64-linux lizard.defaultPackage.x86_64-linux pkgs.xorg.libX11 ];

          buildPhase = ''
            mkdir -p $out
            gcc -c hm.c -lX11 -o dynsys -lgmp -lm
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp dynsys $out/bin
          '';
          postInstall = ''
            cp $out/bin/dynsys ./dynsys
          '';

          meta = with pkgs.lib; {
            description = "hm";
            #license = licenses.mit;
            #maintainers = [ maintainers.yourname ];
            platforms = platforms.unix;
          };
        };
      });

    defaultPackage = { x86_64-linux = self.packages.x86_64-linux.dynsys; };

  };
}

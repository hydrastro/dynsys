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
          buildInputs = [ pkgs.gmp pkgs.stdenv.cc ds.defaultPackage.x86_64-linux lizard.defaultPackage.x86_64-linux pkgs.virtualgl pkgs.freeglut pkgs.glew pkgs.glibc pkgs.cglm pkgs.glfw3 pkgs.freetype];

          buildPhase = ''
            mkdir -p $out
            gcc -c dynsys.c -o dynsys -lgmp -lm -lGL -lGLU -lglut -lds -llizard
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp dynsys $out/bin
          '';

          postInstall = ''
            cp $out/bin/dynsys ./dynsys
          '';

          meta = with pkgs.lib; {
            description = "dynsys";
            #license = licenses.mit;
            #maintainers = [ maintainers.yourname ];
            platforms = platforms.unix;
          };
        };
      });

    defaultPackage = { x86_64-linux = self.packages.x86_64-linux.dynsys; };

  };
}

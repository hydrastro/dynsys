{
  description = "Dynamical-system visualizer using TPCAS infix equations";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs systems;
    in {
      packages = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
          fontPath = "${pkgs.dejavu_fonts}/share/fonts/truetype/DejaVuSansMono.ttf";
        in {
          default = pkgs.stdenv.mkDerivation {
            pname = "dynsys";
            version = "0.1.0";
            src = ./.;

            nativeBuildInputs = [ pkgs.gnumake pkgs.pkg-config ];
            buildInputs = [
              pkgs.cglm
              pkgs.freetype
              pkgs.glew
              pkgs.glfw3
              pkgs.libGL
            ];

            buildPhase = ''
              make MODE=release FONT_PATH=${fontPath}
            '';

            installPhase = ''
              make install PREFIX=$out
            '';
          };
        });

      apps = forAllSystems (system: {
        default = {
          type = "app";
          program = "${self.packages.${system}.default}/bin/dynsys";
        };
      });

      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
          fontPath = "${pkgs.dejavu_fonts}/share/fonts/truetype/DejaVuSansMono.ttf";
        in {
          default = pkgs.mkShell {
            nativeBuildInputs = [ pkgs.gnumake pkgs.pkg-config pkgs.clang-tools ];
            buildInputs = [ pkgs.cglm pkgs.freetype pkgs.glew pkgs.glfw3 pkgs.libGL ];
            DYNSYS_FONT_PATH = fontPath;
          };
        });
    };
}

{
  description = "Development shell for the dynsys TPCAS + Dear ImGui visualizer";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    imgui = {
      url = "github:ocornut/imgui/v1.92.8";
      flake = false;
    };
  };

  outputs = { nixpkgs, imgui, ... }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs systems;
    in {
      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
        in {
          default = pkgs.mkShell {
            nativeBuildInputs = with pkgs; [
              clang
              clang-tools
              gdb
              gnumake
              pkg-config
            ];

            buildInputs = with pkgs; [
              cglm
              glew
              glfw3
              libGL
              xorg.libX11
              xorg.libXcursor
              xorg.libXi
              xorg.libXinerama
              xorg.libXrandr
            ];

            CC = "clang";
            CXX = "clang++";
            IMGUI_DIR = "${imgui}";

            shellHook = ''
              export CC=clang
              export CXX=clang++
              echo "dynsys Dear ImGui dev shell"
              echo "  build: make"
              echo "  run:   make run"
              echo "  clean: make clean"
              echo "  IMGUI_DIR=$IMGUI_DIR"
            '';
          };
        });

      formatter = forAllSystems (system:
        let pkgs = import nixpkgs { inherit system; };
        in pkgs.nixfmt-rfc-style);
    };
}

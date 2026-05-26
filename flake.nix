{
  description = "Development shells for the dynsys TPCAS + Dear ImGui visualizer";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    imgui = {
      url = "github:ocornut/imgui/v1.92.8";
      flake = false;
    };
    # Use the official release archive, not the GitHub source checkout.
    # The GitHub checkout does not contain the generated src/glew.c file
    # needed by this project's simple Windows static GLEW build.
    glew-src = {
      url = "https://downloads.sourceforge.net/project/glew/glew/2.2.0/glew-2.2.0.tgz";
      flake = false;
    };
  };

  outputs = { nixpkgs, imgui, glew-src, ... }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs systems;
    in {
      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
          # Some useful cross packages are blocked only by nixpkgs metadata,
          # not by the actual MinGW build. Keep the normal native shell strict,
          # but allow unsupported host platforms for the Windows cross package set.
          pkgsAllowUnsupported = import nixpkgs {
            inherit system;
            config.allowUnsupportedSystem = true;
          };
          lib = pkgs.lib;
          winPkgs = pkgsAllowUnsupported.pkgsCross.mingwW64;
          winTarget = winPkgs.stdenv.hostPlatform.config;
          winTargetPkgs = [
            winPkgs.glfw3
          ];
          winTargetPkgOutputs = winTargetPkgs ++ (map lib.getDev winTargetPkgs);
          winPkgConfigPath = lib.makeSearchPath "lib/pkgconfig" winTargetPkgOutputs;
          winPkgConfigSharePath = lib.makeSearchPath "share/pkgconfig" winTargetPkgOutputs;
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
              echo "dynsys Dear ImGui native dev shell"
              echo "  build: make"
              echo "  run:   make run"
              echo "  clean: make clean"
              echo "  IMGUI_DIR=$IMGUI_DIR"
              echo "  GLEW_DIR=$GLEW_DIR"
            '';
          };

          # Native NixOS shell containing a MinGW-w64 cross compiler and
          # Windows-target GLFW pkg-config files. GLEW is built from source by
          # the Makefile for Windows because nixpkgs marks winPkgs.glew
          # unsupported for MinGW. cglm is also used as headers only.
          windows = pkgs.mkShell {
            nativeBuildInputs = [
              pkgs.gnumake
              pkgs.pkg-config
              winPkgs.stdenv.cc
            ];

            IMGUI_DIR = "${imgui}";
            GLEW_DIR = "${glew-src}";
            CGLM_INCLUDE_DIR = "${lib.getDev pkgs.cglm}/include";

            PKG_CONFIG_ALLOW_CROSS = "1";
            PKG_CONFIG_LIBDIR = lib.concatStringsSep ":" [ winPkgConfigPath winPkgConfigSharePath ];
            PKG_CONFIG_SYSROOT_DIR = "";

            CC = "${winTarget}-gcc";
            CXX = "${winTarget}-g++";
            AR = "${winTarget}-ar";
            WIN_TRIPLE = winTarget;
            WIN_CC = "${winTarget}-gcc";
            WIN_CXX = "${winTarget}-g++";
            WIN_AR = "${winTarget}-ar";
            WIN_PKG_CONFIG = "pkg-config";

            shellHook = ''
              export CC=${winTarget}-gcc
              export CXX=${winTarget}-g++
              export AR=${winTarget}-ar
              export WIN_TRIPLE=${winTarget}
              export WIN_CC=${winTarget}-gcc
              export WIN_CXX=${winTarget}-g++
              export WIN_AR=${winTarget}-ar
              export WIN_PKG_CONFIG=pkg-config
              export PKG_CONFIG_ALLOW_CROSS=1
              export PKG_CONFIG_LIBDIR=${lib.escapeShellArg (lib.concatStringsSep ":" [ winPkgConfigPath winPkgConfigSharePath ])}
              export PKG_CONFIG_SYSROOT_DIR=
              export GLEW_DIR=${lib.escapeShellArg "${glew-src}"}
              export CGLM_INCLUDE_DIR=${lib.escapeShellArg "${lib.getDev pkgs.cglm}/include"}
              echo "dynsys Windows cross-build shell"
              echo "  target: ${winTarget}"
              echo "  build:  make windows"
              echo "  output: build/windows/dynsys.exe"
              echo "  IMGUI_DIR=$IMGUI_DIR"
              echo "  GLEW_DIR=$GLEW_DIR"
            '';
          };
        });

      formatter = forAllSystems (system:
        let pkgs = import nixpkgs { inherit system; };
        in pkgs.nixfmt-rfc-style);
    };
}

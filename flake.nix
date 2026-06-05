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
    # Phase E (exact symbolic analysis): the Sangaku proof-carrying CAS and the
    # Lizard interpreter it runs on. The dynsys cas_bridge shells out to the
    # `lizard` binary with Sangaku on the module path for exact eigenvalues,
    # equilibria, and certified derivatives. Lizard's own flake pulls its C
    # dependencies (ds + gmp), so we don't re-declare them here.
    lizard = {
      url = "github:hydrastro/lizard";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    sangaku = {
      url = "github:hydrastro/sangaku";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.lizard.follows = "lizard";
    };
  };

  outputs = { nixpkgs, imgui, glew-src, lizard, sangaku, ... }:
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
          # Phase E CAS bridge: the Lizard interpreter binary and the Sangaku
          # library source. lizardBin/bin/lizard runs Sangaku Lisp; sangakuSrc
          # is the library root (so SANGAKU_ROOT points the launcher at it).
          lizardBin = lizard.packages.${system}.default or lizard.defaultPackage.${system};
          sangakuSrc = sangaku;
          winPkgs = pkgsAllowUnsupported.pkgsCross.mingwW64;
          winTarget = winPkgs.stdenv.hostPlatform.config;
          # Force GLFW to be built as a static MinGW library so the final
          # executable does not need a glfw3.dll next to it.
          winGlfwStatic = winPkgs.glfw3.overrideAttrs (old: {
            cmakeFlags = (old.cmakeFlags or [ ]) ++ [
              "-DBUILD_SHARED_LIBS=OFF"
              "-DGLFW_BUILD_DOCS=OFF"
              "-DGLFW_BUILD_EXAMPLES=OFF"
              "-DGLFW_BUILD_TESTS=OFF"
            ];
          });
          winTargetPkgs = [
            winGlfwStatic
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
              lizardBin
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
            # Phase E: where the cas_bridge finds the CAS. LIZARD is the
            # interpreter binary; SANGAKU_ROOT is the library the launcher
            # feeds (prelude + script). The bridge degrades gracefully to the
            # numeric path if these are unset / the binary is absent.
            LIZARD = "${lizardBin}/bin/lizard";
            SANGAKU_ROOT = "${sangakuSrc}";

            shellHook = ''
              export CC=clang
              export CXX=clang++
              export LIZARD=${lib.escapeShellArg "${lizardBin}/bin/lizard"}
              export SANGAKU_ROOT=${lib.escapeShellArg "${sangakuSrc}"}
              echo "dynsys Dear ImGui native dev shell"
              echo "  build: make"
              echo "  run:   make run"
              echo "  clean: make clean"
              echo "  IMGUI_DIR=$IMGUI_DIR"
              echo "  GLEW_DIR=$GLEW_DIR"
              echo "  LIZARD=$LIZARD  (Phase E CAS bridge)"
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
            PKG_CONFIG_ALLOW_SYSTEM_LIBS = "1";
            PKG_CONFIG_ALLOW_SYSTEM_CFLAGS = "1";

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
              export PKG_CONFIG_ALLOW_SYSTEM_LIBS=1
              export PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1
              export GLEW_DIR=${lib.escapeShellArg "${glew-src}"}
              export CGLM_INCLUDE_DIR=${lib.escapeShellArg "${lib.getDev pkgs.cglm}/include"}
              echo "dynsys Windows cross-build shell"
              echo "  target: ${winTarget}"
              echo "  build:  make windows"
              echo "  link:   static MinGW runtime + static GLFW/GLEW; Windows system DLLs remain system-provided"
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

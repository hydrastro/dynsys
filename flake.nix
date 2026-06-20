{
  description = "dynsys — dynamical-system visualizer and analysis toolkit";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    imgui = {
      url = "github:ocornut/imgui/v1.92.8";
      flake = false;
    };

    # The official release archive contains the generated src/glew.c needed by
    # the MinGW static build. The GitHub source checkout does not.
    glew-src = {
      url = "https://downloads.sourceforge.net/project/glew/glew/2.2.0/glew-2.2.0.tgz";
      flake = false;
    };

    # dynsys currently compiles the small TPCAS frontend directly from source.
    tpcas = {
      url = "git+https://github.com/hydrastro/tpcas?rev=4ede7411a5f82f91c114b9b09784a4d06d154b25&submodules=1";
      flake = false;
    };

    # TPCAS no longer vendors ds. Keep the source revision aligned with the ds
    # revision in the locked TPCAS flake and pass it to the existing Makefile as
    # DS_ROOT.
    ds-src = {
      url = "github:hydrastro/ds/9e3224b3aef9b6e271d02b76ff671b1db4301601";
      flake = false;
    };

    # Optional exact-symbolic runtime. It is deliberately kept out of the
    # default shell/package so a Lizard failure cannot block the core build.
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

  outputs =
    {
      self,
      nixpkgs,
      imgui,
      glew-src,
      tpcas,
      ds-src,
      lizard,
      sangaku,
      ...
    }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAllSystems = nixpkgs.lib.genAttrs systems;

      mkPkgs = system: import nixpkgs { inherit system; };

      commonBuildInputs = pkgs: [
        pkgs.cglm
        pkgs.glew
        pkgs.glfw3
        pkgs.libGL
        pkgs.xorg.libX11
        pkgs.xorg.libXcursor
        pkgs.xorg.libXi
        pkgs.xorg.libXinerama
        pkgs.xorg.libXrandr
      ];

      mkDynsys =
        system:
        let
          pkgs = mkPkgs system;
        in
        pkgs.stdenv.mkDerivation {
          pname = "dynsys";
          version = "0.1.0";
          src = ./.;

          strictDeps = true;
          enableParallelBuilding = true;

          nativeBuildInputs = [
            pkgs.gnumake
            pkgs.pkg-config
          ];
          buildInputs = commonBuildInputs pkgs;

          IMGUI_DIR = "${imgui}";
          TPCAS_DIR = "${tpcas}";
          DS_ROOT = "${ds-src}";

          makeFlags = [ "MODE=release" ];

          installPhase = ''
            runHook preInstall
            make install MODE=release PREFIX="$out"
            runHook postInstall
          '';

          meta = {
            description = "Interactive visualizer and analyzer for dynamical systems";
            homepage = "https://github.com/hydrastro/dynsys";
            mainProgram = "dynsys";
            platforms = systems;
          };
        };
    in
    {
      packages = forAllSystems (system: rec {
        dynsys = mkDynsys system;
        default = dynsys;
      });

      apps = forAllSystems (system: {
        default = {
          type = "app";
          program = "${self.packages.${system}.default}/bin/dynsys";
        };
      });

      checks = forAllSystems (system: {
        package = self.packages.${system}.default;
      });

      devShells = forAllSystems (
        system:
        let
          pkgs = mkPkgs system;
          lib = pkgs.lib;
          hasLizard = builtins.hasAttr system lizard.packages;
          lizardBin =
            if hasLizard then
              lizard.packages.${system}.default or lizard.defaultPackage.${system}
            else
              null;

          nativeTools = [
            pkgs.clang
            pkgs.clang-tools
            pkgs.gdb
            pkgs.gnumake
            pkgs.pkg-config
          ];

          mkNativeShell =
            extraPackages: extraAttrs:
            pkgs.mkShell (
              {
                inputsFrom = [ self.packages.${system}.default ];
                packages = nativeTools ++ extraPackages;

                CC = "clang";
                CXX = "clang++";
                IMGUI_DIR = "${imgui}";
                TPCAS_DIR = "${tpcas}";
                DS_ROOT = "${ds-src}";
              }
              // extraAttrs
            );

          pkgsAllowUnsupported = import nixpkgs {
            inherit system;
            config.allowUnsupportedSystem = true;
          };
          winPkgs = pkgsAllowUnsupported.pkgsCross.mingwW64;
          winTarget = winPkgs.stdenv.hostPlatform.config;
          winGlfwStatic = winPkgs.glfw3.overrideAttrs (old: {
            cmakeFlags = (old.cmakeFlags or [ ]) ++ [
              "-DBUILD_SHARED_LIBS=OFF"
              "-DGLFW_BUILD_DOCS=OFF"
              "-DGLFW_BUILD_EXAMPLES=OFF"
              "-DGLFW_BUILD_TESTS=OFF"
            ];
          });
          winTargetPkgOutputs = [
            winGlfwStatic
            (lib.getDev winGlfwStatic)
          ];
          winPkgConfigPath = lib.makeSearchPath "lib/pkgconfig" winTargetPkgOutputs;
          winPkgConfigSharePath = lib.makeSearchPath "share/pkgconfig" winTargetPkgOutputs;
        in
        {
          default = mkNativeShell [ ] {
            shellHook = ''
              echo "dynsys native development shell"
              echo "  build package: nix build"
              echo "  development:   make"
              echo "  run:           make run"
              echo "  IMGUI_DIR=$IMGUI_DIR"
              echo "  TPCAS_DIR=$TPCAS_DIR"
              echo "  DS_ROOT=$DS_ROOT"
            '';
          };

          windows = pkgs.mkShell {
            nativeBuildInputs = [
              pkgs.gnumake
              pkgs.pkg-config
              winPkgs.stdenv.cc
            ];

            IMGUI_DIR = "${imgui}";
            GLEW_DIR = "${glew-src}";
            CGLM_INCLUDE_DIR = "${lib.getDev pkgs.cglm}/include";
            TPCAS_DIR = "${tpcas}";
            DS_ROOT = "${ds-src}";

            PKG_CONFIG_ALLOW_CROSS = "1";
            PKG_CONFIG_LIBDIR = lib.concatStringsSep ":" [
              winPkgConfigPath
              winPkgConfigSharePath
            ];
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
              echo "dynsys Windows cross-build shell"
              echo "  target: ${winTarget}"
              echo "  build:  make windows"
              echo "  output: build/windows/dynsys.exe"
              echo "  IMGUI_DIR=$IMGUI_DIR"
              echo "  GLEW_DIR=$GLEW_DIR"
              echo "  TPCAS_DIR=$TPCAS_DIR"
              echo "  DS_ROOT=$DS_ROOT"
            '';
          };
        }
        // lib.optionalAttrs hasLizard {
          cas = mkNativeShell [ lizardBin ] {
            LIZARD = "${lizardBin}/bin/lizard";
            SANGAKU_ROOT = "${sangaku}";

            shellHook = ''
              echo "dynsys development shell with exact-symbolic CAS support"
              echo "  build:        make"
              echo "  run:          make run"
              echo "  LIZARD=$LIZARD"
              echo "  SANGAKU_ROOT=$SANGAKU_ROOT"
              echo "  TPCAS_DIR=$TPCAS_DIR"
              echo "  DS_ROOT=$DS_ROOT"
            '';
          };
        }
      );

      formatter = forAllSystems (system: (mkPkgs system).nixfmt-rfc-style);
    };
}

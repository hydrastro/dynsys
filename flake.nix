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

    # Requires the TPCAS DSA migration. Refresh this input after publishing it.
    tpcas = {
      url = "github:hydrastro/tpcas";
      flake = false;
    };

    # Build DSA separately; after the rename change the owner/repo to
    # hydrastro/dsa and refresh this input's lock.
    dsa-src = {
      url = "github:hydrastro/ds/909b4a462410b190db2488cb0fac7fca088af81c";
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
      dsa-src,
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

      # Use DSA's build system, keeping its Unicode generator on the build
      # machine when targeting MinGW. Windows omits the two POSIX I/O adapters.
      mkDsa = pkgs: pkgs.stdenv.mkDerivation {
        pname = "dsa";
        version = "909b4a4";
        src = dsa-src;
        nativeBuildInputs = [ pkgs.buildPackages.gnumake ];
        dontConfigure = true;
        buildPhase = ''
          runHook preBuild
          make build/tools/gen_ucd CC=${pkgs.buildPackages.stdenv.cc}/bin/cc
          ${pkgs.lib.optionalString pkgs.stdenv.hostPlatform.isWindows ''
            sources=$(find src/core src/structures src/algorithms -name '*.c' ! -name posix_store.c ! -name str_io_posix.c | sort | tr '\n' ' ')
          ''}
          make build/libdsa.a CC="$CC" AR="$AR" \
            ${pkgs.lib.optionalString pkgs.stdenv.hostPlatform.isWindows ''LIB_SRC="$sources"''}
          runHook postBuild
        '';
        installPhase = ''
          runHook preInstall
          mkdir -p "$out/include" "$out/lib"
          cp -R include/. "$out/include/"
          cp build/libdsa.a "$out/lib/"
          runHook postInstall
        '';
      };

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
            pkgs.patch
            pkgs.pkg-config
          ];
          buildInputs = commonBuildInputs pkgs;

          IMGUI_DIR = "${imgui}";
          TPCAS_DIR = "${tpcas}";
          DSA_ROOT = "${mkDsa pkgs}";

          makeFlags = [ "MODE=release" ];
          doCheck = true;
          checkTarget = "test";

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
            pkgs.patch
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
                DSA_ROOT = "${mkDsa pkgs}";
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
              echo "  DSA_ROOT=$DSA_ROOT"
            '';
          };

          windows = pkgs.mkShell {
            nativeBuildInputs = [
              pkgs.gnumake
              pkgs.patch
              pkgs.pkg-config
              winPkgs.stdenv.cc
            ];

            IMGUI_DIR = "${imgui}";
            GLEW_DIR = "${glew-src}";
            CGLM_INCLUDE_DIR = "${lib.getDev pkgs.cglm}/include";
            TPCAS_DIR = "${tpcas}";
            DSA_ROOT = "${mkDsa winPkgs}";

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
              echo "  DSA_ROOT=$DSA_ROOT"
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
              echo "  DSA_ROOT=$DSA_ROOT"
            '';
          };
        }
      );

      formatter = forAllSystems (system: (mkPkgs system).nixfmt-rfc-style);
    };
}

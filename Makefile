# dynsys — GLFW/OpenGL dynamical-system visualizer with Dear ImGui + TPCAS
#
# Native Linux workflow:
#   nix develop
#   make
#   make run
#
# Windows cross-build workflow from NixOS:
#   nix develop .#windows
#   make windows
#
# The Windows build writes build/windows/dynsys.exe so it does not reuse
# native Linux object files.

CC         ?= cc
CXX        ?= c++
AR         ?= ar
RM         ?= rm -f
RMDIR      ?= rm -rf
MKDIR_P    ?= mkdir -p
INSTALL    ?= install
PKG_CONFIG ?= pkg-config

# Set TARGET_OS=windows directly, or use `make windows` / `make build-windows`.
TARGET_OS ?= native
HOST_TRIPLE := $(shell $(CC) -dumpmachine 2>/dev/null || echo unknown)
WINDOWS_BUILD := 0
ifeq ($(TARGET_OS),windows)
  WINDOWS_BUILD := 1
else ifneq (,$(findstring mingw,$(HOST_TRIPLE)))
  WINDOWS_BUILD := 1
else ifneq (,$(findstring windows,$(HOST_TRIPLE)))
  WINDOWS_BUILD := 1
endif

BUILD_ROOT ?= build
ifeq ($(WINDOWS_BUILD),1)
  EXEEXT ?= .exe
  BUILD_DIR := $(BUILD_ROOT)/windows
else
  EXEEXT ?=
  BUILD_DIR := $(BUILD_ROOT)
endif

SRC_DIR    := src
TPCAS_DIR  := vendor/tpcas
DS_DIR     := $(TPCAS_DIR)/vendor/ds/lib
OBJ_DIR    := $(BUILD_DIR)/obj
C_OBJ_DIR  := $(BUILD_DIR)/obj-c
CXX_OBJ_DIR:= $(BUILD_DIR)/obj-cxx
DEP_DIR    := $(BUILD_DIR)/dep
C_DEP_DIR  := $(BUILD_DIR)/dep-c
CXX_DEP_DIR:= $(BUILD_DIR)/dep-cxx
TARGET     := $(BUILD_DIR)/dynsys$(EXEEXT)
IR_TEST_TARGET := $(BUILD_DIR)/ir_smoke$(EXEEXT)
ANALYSIS_TEST_TARGET := $(BUILD_DIR)/analysis_smoke$(EXEEXT)
AD_TEST_TARGET := $(BUILD_DIR)/ad_smoke$(EXEEXT)
NULLCLINE_TEST_TARGET := $(BUILD_DIR)/nullcline_smoke$(EXEEXT)
DIM_TEST_TARGET := $(BUILD_DIR)/dim_detect_smoke$(EXEEXT)
FP_TEST_TARGET := $(BUILD_DIR)/fixedpoints_smoke$(EXEEXT)
LYAP_TEST_TARGET := $(BUILD_DIR)/lyapunov_smoke$(EXEEXT)
FRACTAL_TEST_TARGET := $(BUILD_DIR)/fractal_smoke$(EXEEXT)
BRIDGE_TEST_TARGET := $(BUILD_DIR)/bridge_smoke$(EXEEXT)
BASIN_TEST_TARGET := $(BUILD_DIR)/basin_smoke$(EXEEXT)
SOLVER_TEST_TARGET := $(BUILD_DIR)/solver_smoke$(EXEEXT)
SCAN_TEST_TARGET := $(BUILD_DIR)/scan_smoke$(EXEEXT)
ODEBIF_TEST_TARGET := $(BUILD_DIR)/odebif_smoke$(EXEEXT)
PROG_TEST_TARGET := $(BUILD_DIR)/progressive_smoke$(EXEEXT)
BASINCHAOS_TEST_TARGET := $(BUILD_DIR)/basinchaos_smoke$(EXEEXT)
CONT_TEST_TARGET := $(BUILD_DIR)/continuation_smoke$(EXEEXT)
PNG_TEST_TARGET := $(BUILD_DIR)/png_smoke$(EXEEXT)
PARAMSYNC_TEST_TARGET := $(BUILD_DIR)/paramsync_smoke$(EXEEXT)
PERIOD_TEST_TARGET := $(BUILD_DIR)/period_smoke$(EXEEXT)
TEST_OBJ_DIR := $(BUILD_DIR)/obj-test

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

MODE ?= debug
PKGS_NATIVE ?= glfw3 glew cglm
PKGS_WINDOWS ?= glfw3
ifeq ($(WINDOWS_BUILD),1)
  PKGS ?= $(PKGS_WINDOWS)
else
  PKGS ?= $(PKGS_NATIVE)
endif
IMGUI_DIR ?=
CGLM_INCLUDE_DIR ?=
GLEW_DIR ?=

# Used by the recursive `make windows` target. Override these if your MinGW
# toolchain uses different names.
WIN_TRIPLE ?= x86_64-w64-mingw32
WIN_CC ?= $(WIN_TRIPLE)-gcc
WIN_CXX ?= $(WIN_TRIPLE)-g++
WIN_AR ?= $(WIN_TRIPLE)-ar
WIN_PKG_CONFIG ?= pkg-config
WIN_BUILD_DIR ?= $(BUILD_ROOT)/windows
WIN_TARGET ?= $(WIN_BUILD_DIR)/dynsys.exe

WARNINGS_C   := -Wall -Wextra -Wpedantic
WARNINGS_CXX := -Wall -Wextra
CSTD         := -std=c11
CXXSTD       := -std=c++17

ifeq ($(WINDOWS_BUILD),1)
  # For the Windows build, ask pkg-config for static link flags so GLFW's
  # transitive system libraries are included. This still links against normal
  # Windows system DLLs such as kernel32/user32/gdi32/opengl32, which are
  # supplied by Windows.
  PKG_CONFIG_LIBS_ARGS ?= --libs --static
else
  PKG_CONFIG_LIBS_ARGS ?= --libs
endif

PKG_CFLAGS := $(shell $(PKG_CONFIG) --cflags $(PKGS) 2>/dev/null)
PKG_LIBS   := $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_ARGS) $(PKGS) 2>/dev/null)

CPPFLAGS ?=
CPPFLAGS += -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds $(PKG_CFLAGS)
ifneq ($(strip $(IMGUI_DIR)),)
  CPPFLAGS += -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends
endif
ifeq ($(WINDOWS_BUILD),1)
  CPPFLAGS += -D_WIN32_WINNT=0x0601 -DGLEW_STATIC
  ifneq ($(strip $(CGLM_INCLUDE_DIR)),)
    CPPFLAGS += -I$(CGLM_INCLUDE_DIR)
  endif
  ifneq ($(strip $(GLEW_DIR)),)
    CPPFLAGS += -I$(GLEW_DIR)/include
  endif
endif

CFLAGS ?=
CFLAGS += $(CSTD) $(WARNINGS_C)
CXXFLAGS ?=
CXXFLAGS += $(CXXSTD) $(WARNINGS_CXX) -pthread
LDFLAGS ?=
LDLIBS ?=

ifeq ($(WINDOWS_BUILD),1)
  # MinGW/Windows OpenGL + GLFW backend system libraries. GLEW is compiled
  # from source into this project for Windows, so do not link -lglew32 here.
  LDLIBS += $(PKG_LIBS) -lopengl32 -lgdi32 -limm32 -lole32 -luuid -lwinmm -lm
  # Produce a self-contained MinGW executable for third-party/runtime
  # libraries: libstdc++, libgcc, winpthread, and static GLFW/GLEW where
  # available. Windows system DLLs remain dynamically loaded by design.
  WIN_STATIC_RUNTIME ?= 1
  ifeq ($(WIN_STATIC_RUNTIME),1)
    LDFLAGS += -static -static-libgcc -static-libstdc++
  endif
else
  LDLIBS += $(PKG_LIBS) -lGL -ldl -lm -pthread
endif

ifeq ($(MODE),release)
  CFLAGS += -O2 -DNDEBUG
  CXXFLAGS += -O2 -DNDEBUG
else ifeq ($(MODE),asan)
  CFLAGS += -O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
  CXXFLAGS += -O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
  LDFLAGS += -fsanitize=address,undefined
else
  # -Og keeps full debuggability while providing the minimal optimization
  # that glibc's _FORTIFY_SOURCE (injected by the nix stdenv) requires;
  # plain -O0 triggers a "_FORTIFY_SOURCE requires compiling with
  # optimization" warning on every translation unit.
  CFLAGS += -Og -g3
  CXXFLAGS += -Og -g3
endif

DYNSYS_CPP_SRCS := $(SRC_DIR)/dynsys.cpp $(SRC_DIR)/expr_ir.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/expr_ir_ad.cpp
DYNSYS_OBJS := $(patsubst %.cpp,$(CXX_OBJ_DIR)/%.o,$(DYNSYS_CPP_SRCS))
DYNSYS_DEPS := $(patsubst %.cpp,$(CXX_DEP_DIR)/%.d,$(DYNSYS_CPP_SRCS))

IR_TEST_CPP_SRCS := $(SRC_DIR)/expr_ir.cpp test/ir_smoke.cpp
IR_TEST_CPP_OBJS := $(patsubst %.cpp,$(TEST_OBJ_DIR)/%.o,$(IR_TEST_CPP_SRCS))

TPCAS_SRCS := \
  $(TPCAS_DIR)/src/arena.c \
  $(TPCAS_DIR)/src/ast.c \
  $(TPCAS_DIR)/src/lex.c \
  $(TPCAS_DIR)/src/pratt.c
DS_SRCS := \
  $(DS_DIR)/common.c \
  $(DS_DIR)/status.c \
  $(DS_DIR)/error.c \
  $(DS_DIR)/diagnostic.c \
  $(DS_DIR)/context.c \
  $(DS_DIR)/allocators.c
C_SRCS := $(TPCAS_SRCS) $(DS_SRCS)
C_OBJS := $(patsubst %.c,$(C_OBJ_DIR)/%.o,$(C_SRCS))
C_DEPS := $(patsubst %.c,$(C_DEP_DIR)/%.d,$(C_SRCS))
ifeq ($(WINDOWS_BUILD),1)
  GLEW_OBJ := $(C_OBJ_DIR)/glew/glew.o
  GLEW_DEP := $(C_DEP_DIR)/glew/glew.d
else
  GLEW_OBJ :=
  GLEW_DEP :=
endif
IR_TEST_C_OBJS := $(patsubst %.c,$(TEST_OBJ_DIR)/%.o,$(C_SRCS))
IR_TEST_OBJS := $(IR_TEST_CPP_OBJS) $(IR_TEST_C_OBJS)
PARAMSYNC_OBJS := $(TEST_OBJ_DIR)/test/paramsync_smoke.o $(IR_TEST_C_OBJS) $(TEST_OBJ_DIR)/src/expr_ir.o
ESHADOW_OBJS := $(TEST_OBJ_DIR)/test/e_shadow_smoke.o $(IR_TEST_C_OBJS) $(TEST_OBJ_DIR)/src/expr_ir.o
IFSLIT_OBJS := $(TEST_OBJ_DIR)/test/ifslit_smoke.o $(IR_TEST_C_OBJS) $(TEST_OBJ_DIR)/src/expr_ir.o
BOXDIM_TEST_TARGET := $(BUILD_DIR)/boxdim_smoke$(EXEEXT)
IFS_TEST_TARGET := $(BUILD_DIR)/ifs_smoke$(EXEEXT)
LC_TEST_TARGET := $(BUILD_DIR)/limitcycle_smoke$(EXEEXT)
LCSWEEP_TEST_TARGET := $(BUILD_DIR)/lcsweep_smoke$(EXEEXT)
IFSMODEL_TEST_TARGET := $(BUILD_DIR)/ifsmodel_smoke$(EXEEXT)
IFSPARAM_TEST_TARGET := $(BUILD_DIR)/ifsparam_smoke$(EXEEXT)
IFSLIT_TEST_TARGET := $(BUILD_DIR)/ifslit_smoke$(EXEEXT)
CAS_TEST_TARGET := $(BUILD_DIR)/cas_bridge_smoke$(EXEEXT)
HOPFL1_TEST_TARGET := $(BUILD_DIR)/hopf_l1_smoke$(EXEEXT)

IMGUI_CORE_SRCS := imgui.cpp imgui_draw.cpp imgui_tables.cpp imgui_widgets.cpp
IMGUI_BACKEND_SRCS := imgui_impl_glfw.cpp imgui_impl_opengl3.cpp
IMGUI_CORE_OBJS := $(addprefix $(CXX_OBJ_DIR)/imgui/,$(IMGUI_CORE_SRCS:.cpp=.o))
IMGUI_BACKEND_OBJS := $(addprefix $(CXX_OBJ_DIR)/imgui/backends/,$(IMGUI_BACKEND_SRCS:.cpp=.o))
IMGUI_OBJS := $(IMGUI_CORE_OBJS) $(IMGUI_BACKEND_OBJS)
IMGUI_DEPS := $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$(IMGUI_OBJS))

OBJS := $(DYNSYS_OBJS) $(C_OBJS) $(GLEW_OBJ) $(IMGUI_OBJS)
DEPS := $(DYNSYS_DEPS) $(C_DEPS) $(GLEW_DEP) $(IMGUI_DEPS)

.PHONY: all build check-deps check-legacy prune-legacy run headless headless-ast headless-smoke bench test ir-smoke test-analysis test-ad test-nullcline test-dim test-fp test-lyap test-fractal test-bridge test-bridgefamily test-lpcarc test-branchsw test-btcodim2 test-basinsmt test-homoclinic test-homocont test-validation test-codim2coef test-homoseed test-zhhh test-lcseed test-pdns test-pdcurve test-hetero test-bpc test-codim2cyc test-lindiag test-findhomo debug release asan windows build-windows clean distclean install uninstall format print-vars help

all: check-deps check-legacy $(TARGET)

build: all

windows build-windows:
	$(MAKE) TARGET_OS=windows \
	  BUILD_DIR="$(WIN_BUILD_DIR)" \
	  TARGET="$(WIN_TARGET)" \
	  CC="$(WIN_CC)" \
	  CXX="$(WIN_CXX)" \
	  AR="$(WIN_AR)" \
	  PKG_CONFIG="$(WIN_PKG_CONFIG)" \
	  PKGS="$(PKGS_WINDOWS)"

check-deps:
	@if [ "$(WINDOWS_BUILD)" = "1" ] && ! command -v "$(CXX)" >/dev/null 2>&1; then \
	  echo "error: Windows cross compiler not found: $(CXX)" >&2; \
	  echo "hint: enter the shell first: nix develop .#windows" >&2; \
	  exit 1; \
	fi
	@$(PKG_CONFIG) --exists $(PKGS) || { \
	  echo "error: missing pkg-config dependencies: $(PKGS)" >&2; \
	  if [ "$(WINDOWS_BUILD)" = "1" ]; then \
	    echo "hint: run 'nix develop .#windows' first, then run 'make windows' again." >&2; \
	  else \
	    echo "hint: run 'nix develop' first, then run 'make' again." >&2; \
	  fi; \
	  exit 1; \
	}
	@if [ "$(WINDOWS_BUILD)" = "1" ]; then \
	  printf '#include <cglm/cglm.h>\n' | $(CC) $(CPPFLAGS) -x c -E - >/dev/null 2>&1 || { \
	    echo "error: cglm headers not found for the Windows build" >&2; \
	    echo "hint: enter the updated shell with: nix develop .#windows" >&2; \
	    echo "      or set CGLM_INCLUDE_DIR=/path/to/cglm/include" >&2; \
	    exit 1; \
	  }; \
	  test -f "$(GLEW_DIR)/src/glew.c" || { \
	    echo "error: GLEW_DIR does not contain generated GLEW source: $(GLEW_DIR)/src/glew.c" >&2; \
	    echo "hint: use the updated flake that points glew-src at the official glew-2.2.0.tgz release archive," >&2; \
	    echo "      then refresh the lock file with: nix flake lock --update-input glew-src" >&2; \
	    echo "      or set GLEW_DIR=/path/to/unpacked/glew-2.2.0-release" >&2; \
	    exit 1; \
	  }; \
	fi
	@if [ -z "$(strip $(IMGUI_DIR))" ]; then \
	  echo "error: IMGUI_DIR is not set" >&2; \
	  echo "hint: run 'nix develop' or 'nix develop .#windows' first; the flake exports IMGUI_DIR." >&2; \
	  echo "      or run: make IMGUI_DIR=/path/to/imgui" >&2; \
	  exit 1; \
	fi
	@test -f "$(IMGUI_DIR)/imgui.cpp" || { \
	  echo "error: IMGUI_DIR does not point at a Dear ImGui source tree: $(IMGUI_DIR)" >&2; \
	  exit 1; \
	}

check-legacy:
	@if [ -f "$(SRC_DIR)/dynsys.c" ]; then \
	  echo "warning: legacy $(SRC_DIR)/dynsys.c exists but is intentionally ignored; run 'make prune-legacy CLEAN_APPLY=1' to remove old FreeType code." >&2; \
	fi

prune-legacy:
	@if [ "$(CLEAN_APPLY)" = "1" ]; then \
	  scripts/prune-legacy.sh --apply; \
	else \
	  scripts/prune-legacy.sh; \
	fi

$(TARGET): $(OBJS)
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(C_OBJ_DIR)/%.o: %.c
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(C_OBJ_DIR)/%.o,$(C_DEP_DIR)/%.d,$@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -MF $(patsubst $(C_OBJ_DIR)/%.o,$(C_DEP_DIR)/%.d,$@) -c $< -o $@

$(CXX_OBJ_DIR)/%.o: %.cpp
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@) -c $< -o $@

$(GLEW_OBJ): $(GLEW_DIR)/src/glew.c
	@$(MKDIR_P) $(dir $@) $(dir $(GLEW_DEP))
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -MF $(GLEW_DEP) -c $< -o $@

$(CXX_OBJ_DIR)/imgui/%.o: $(IMGUI_DIR)/%.cpp
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@) -c $< -o $@

$(CXX_OBJ_DIR)/imgui/backends/%.o: $(IMGUI_DIR)/backends/%.cpp
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@) -c $< -o $@

test: test-analysis test-ad test-nullcline test-dim test-fp test-lyap test-fractal test-bridge test-bridgefamily test-basin test-solver test-scan test-odebif test-progressive test-basinchaos test-continuation test-period test-png test-paramsync test-boxdim test-ifs test-limitcycle test-lcsweep test-ifsmodel test-ifsparam test-ifslit test-cas test-hopfl1 test-foldnf test-codim2 test-twoparam test-lccolloc test-tpc2 test-lpc test-lpcarc test-branchsw test-btcodim2 test-basinsmt test-homoclinic test-homocont test-validation test-codim2coef test-homoseed test-zhhh test-lcseed test-pdns test-pdcurve test-hetero test-bpc test-codim2cyc test-lindiag test-findhomo test-lpccurve test-eshadow test-bridgealign test-projsolid

test-analysis: $(ANALYSIS_TEST_TARGET)
	./$(ANALYSIS_TEST_TARGET)

$(ANALYSIS_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/analysis_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/analysis_smoke.cpp -o $@ -lm

test-nullcline: $(NULLCLINE_TEST_TARGET)
	./$(NULLCLINE_TEST_TARGET)

$(NULLCLINE_TEST_TARGET): test/nullcline_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/nullcline_smoke.cpp -o $@ -lm

test-dim: $(DIM_TEST_TARGET)
	./$(DIM_TEST_TARGET)

$(DIM_TEST_TARGET): test/dim_detect_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/dim_detect_smoke.cpp -o $@ -lm

test-fp: $(FP_TEST_TARGET)
	./$(FP_TEST_TARGET)

$(FP_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/fixedpoints_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/fixedpoints_smoke.cpp -o $@ -lm

test-lyap: $(LYAP_TEST_TARGET)
	./$(LYAP_TEST_TARGET)

$(LYAP_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/lyapunov_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/lyapunov_smoke.cpp -o $@ -lm

test-fractal: $(FRACTAL_TEST_TARGET)
	./$(FRACTAL_TEST_TARGET)

$(FRACTAL_TEST_TARGET): test/fractal_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/fractal_smoke.cpp -o $@ -lm

test-bridge: $(BRIDGE_TEST_TARGET)
	./$(BRIDGE_TEST_TARGET)

test-bridgefamily: $(BUILD_DIR)/bridge_family_smoke$(EXEEXT)
	./$(BUILD_DIR)/bridge_family_smoke$(EXEEXT)

$(BUILD_DIR)/bridge_family_smoke$(EXEEXT): test/bridge_family_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/bridge_family_smoke.cpp -o $@ -lm

$(BRIDGE_TEST_TARGET): test/bridge_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/bridge_smoke.cpp -o $@ -lm

test-basin: $(BASIN_TEST_TARGET)
	./$(BASIN_TEST_TARGET)

$(BASIN_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/basin_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/basin_smoke.cpp -o $@ -lm

test-solver: $(SOLVER_TEST_TARGET)
	./$(SOLVER_TEST_TARGET)

$(SOLVER_TEST_TARGET): test/solver_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/solver_smoke.cpp -o $@ -lm

test-scan: $(SCAN_TEST_TARGET)
	./$(SCAN_TEST_TARGET)

$(SCAN_TEST_TARGET): test/scan_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/scan_smoke.cpp -o $@ -lm

test-odebif: $(ODEBIF_TEST_TARGET)
	./$(ODEBIF_TEST_TARGET)

$(ODEBIF_TEST_TARGET): test/odebif_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/odebif_smoke.cpp -o $@ -lm

test-progressive: $(PROG_TEST_TARGET)
	./$(PROG_TEST_TARGET)

$(PROG_TEST_TARGET): test/progressive_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/progressive_smoke.cpp -o $@ -lm

test-basinchaos: $(BASINCHAOS_TEST_TARGET)
	./$(BASINCHAOS_TEST_TARGET)

$(BASINCHAOS_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/basinchaos_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/basinchaos_smoke.cpp -o $@ -lm

test-continuation: $(CONT_TEST_TARGET)
	./$(CONT_TEST_TARGET)

$(CONT_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/continuation_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/continuation_smoke.cpp -o $@ -lm

test-png: $(PNG_TEST_TARGET)
	./$(PNG_TEST_TARGET)

$(PNG_TEST_TARGET): test/png_smoke.cpp $(SRC_DIR)/png_writer.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/png_smoke.cpp -o $@ -lm

test-period: $(PERIOD_TEST_TARGET)
	./$(PERIOD_TEST_TARGET)

$(PERIOD_TEST_TARGET): test/period_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/period_smoke.cpp -o $@ -lm

test-ad: $(AD_TEST_TARGET)
	./$(AD_TEST_TARGET)

AD_TPCAS_C := $(TPCAS_DIR)/src/ast.c $(TPCAS_DIR)/src/arena.c \
  $(DS_DIR)/common.c $(DS_DIR)/status.c $(DS_DIR)/error.c \
  $(DS_DIR)/diagnostic.c $(DS_DIR)/context.c $(DS_DIR)/allocators.c
AD_INCLUDES := -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds
$(AD_TEST_TARGET): $(SRC_DIR)/expr_ir.cpp $(SRC_DIR)/expr_ir_ad.cpp test/ad_smoke.cpp $(AD_TPCAS_C)
	@$(MKDIR_P) $(dir $@) $(BUILD_DIR)/ad-cobj
	@for c in $(AD_TPCAS_C); do \
	  $(CC) $(CSTD) -O2 $(AD_INCLUDES) -c $$c -o $(BUILD_DIR)/ad-cobj/`basename $${c%.c}`.o || exit 1; \
	done
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 $(AD_INCLUDES) $(SRC_DIR)/expr_ir.cpp $(SRC_DIR)/expr_ir_ad.cpp test/ad_smoke.cpp $(BUILD_DIR)/ad-cobj/*.o -o $@ -lm

ir-smoke: $(IR_TEST_TARGET)
	./$(IR_TEST_TARGET)

$(IR_TEST_TARGET): $(IR_TEST_OBJS)
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(IR_TEST_OBJS) -lm

$(TEST_OBJ_DIR)/%.o: %.c
	@$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TEST_OBJ_DIR)/%.o: %.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

run: all
	./$(TARGET)

headless: all
	./$(TARGET) --headless $(ARGS)

headless-ast: all
	./$(TARGET) --headless $(ARGS) --use-ast

headless-smoke: all
	./$(TARGET) --headless examples/lorenz.dyn --steps 10000

bench:
	$(MAKE) MODE=release
	./$(TARGET) --headless examples/lorenz.dyn --steps 200000
	./$(TARGET) --headless examples/lorenz.dyn --steps 200000 --use-ast

debug:
	$(MAKE) MODE=debug

release:
	$(MAKE) MODE=release

asan:
	$(MAKE) MODE=asan

install: all
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 0755 $(TARGET) $(DESTDIR)$(BINDIR)/dynsys$(EXEEXT)

uninstall:
	$(RM) $(DESTDIR)$(BINDIR)/dynsys $(DESTDIR)$(BINDIR)/dynsys.exe

clean:
	$(RMDIR) $(BUILD_ROOT)

distclean: clean
	@# nix build symlink
	$(RM) result
	@# windows build root (if built out-of-tree separately)
	$(RMDIR) $(WIN_BUILD_DIR)
	@# runtime artifacts dynsys writes into the working tree:
	@#   - timestamped screenshots dynsys_YYYYMMDD_HHMMSS.png (Save PNG)
	@#   - the default trajectory/Poincare export (dynsys_export.csv)
	$(RM) dynsys_export.csv
	$(RM) dynsys_[0-9]*.png

format:
	clang-format -i $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h test/*.cpp

print-vars:
	@echo "CC=$(CC)"
	@echo "CXX=$(CXX)"
	@echo "AR=$(AR)"
	@echo "PKG_CONFIG=$(PKG_CONFIG)"
	@echo "HOST_TRIPLE=$(HOST_TRIPLE)"
	@echo "TARGET_OS=$(TARGET_OS)"
	@echo "WINDOWS_BUILD=$(WINDOWS_BUILD)"
	@echo "MODE=$(MODE)"
	@echo "BUILD_DIR=$(BUILD_DIR)"
	@echo "TARGET=$(TARGET)"
	@echo "C_OBJ_DIR=$(C_OBJ_DIR)"
	@echo "CXX_OBJ_DIR=$(CXX_OBJ_DIR)"
	@echo "PKGS=$(PKGS)"
	@echo "IMGUI_DIR=$(IMGUI_DIR)"
	@echo "PKG_CFLAGS=$(PKG_CFLAGS)"
	@echo "PKG_CONFIG_LIBS_ARGS=$(PKG_CONFIG_LIBS_ARGS)"
	@echo "PKG_LIBS=$(PKG_LIBS)"
	@echo "LDFLAGS=$(LDFLAGS)"
	@echo "LDLIBS=$(LDLIBS)"

help:
	@echo "Targets:"
	@echo "  make                 build native $(TARGET)"
	@echo "  make run             build and run native executable"
	@echo "  make windows         cross-build build/windows/dynsys.exe"
	@echo "  make build-windows   same as make windows"
	@echo "  make test            build and run standalone IR smoke test"
	@echo "  make headless        run headless; pass ARGS='examples/lorenz.dyn --steps 10000'"
	@echo "  make bench           release IR/AST headless comparison on Lorenz"
	@echo "  make release         optimized native build"

-include $(DEPS)

test-paramsync: $(PARAMSYNC_TEST_TARGET)
	./$(PARAMSYNC_TEST_TARGET)

$(TEST_OBJ_DIR)/test/paramsync_smoke.o: test/paramsync_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) -O2 -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds -c test/paramsync_smoke.cpp -o $@

$(PARAMSYNC_TEST_TARGET): $(PARAMSYNC_OBJS)
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) -o $@ $(PARAMSYNC_OBJS) -lm

ESHADOW_TEST_TARGET := $(BUILD_DIR)/e_shadow_smoke$(EXEEXT)
test-eshadow: $(ESHADOW_TEST_TARGET)
	./$(ESHADOW_TEST_TARGET)

$(TEST_OBJ_DIR)/test/e_shadow_smoke.o: test/e_shadow_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) -O2 -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds -c test/e_shadow_smoke.cpp -o $@

$(ESHADOW_TEST_TARGET): $(ESHADOW_OBJS)
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) -o $@ $(ESHADOW_OBJS) -lm

test-boxdim: $(BOXDIM_TEST_TARGET)
	./$(BOXDIM_TEST_TARGET)

$(BOXDIM_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/boxdim_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/boxdim_smoke.cpp -o $@ -lm

test-ifs: $(IFS_TEST_TARGET)
	./$(IFS_TEST_TARGET)

$(IFS_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/ifs_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/ifs_smoke.cpp -o $@ -lm

test-limitcycle: $(LC_TEST_TARGET)
	./$(LC_TEST_TARGET)

$(LC_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/limitcycle_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/limitcycle_smoke.cpp -o $@ -lm

test-lcsweep: $(LCSWEEP_TEST_TARGET)
	./$(LCSWEEP_TEST_TARGET)

$(LCSWEEP_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/lcsweep_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/lcsweep_smoke.cpp -o $@ -lm

test-ifsmodel: $(IFSMODEL_TEST_TARGET)
	./$(IFSMODEL_TEST_TARGET)

$(IFSMODEL_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/ifsmodel_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/ifsmodel_smoke.cpp -o $@ -lm

test-ifsparam: $(IFSPARAM_TEST_TARGET)
	./$(IFSPARAM_TEST_TARGET)

$(IFSPARAM_TEST_TARGET): $(SRC_DIR)/analysis.cpp test/ifsparam_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) $(SRC_DIR)/analysis.cpp test/ifsparam_smoke.cpp -o $@ -lm

test-ifslit: $(IFSLIT_TEST_TARGET)
	./$(IFSLIT_TEST_TARGET)

$(TEST_OBJ_DIR)/test/ifslit_smoke.o: test/ifslit_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) -O2 -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds -c test/ifslit_smoke.cpp -o $@

$(IFSLIT_TEST_TARGET): $(IFSLIT_OBJS)
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) -o $@ $(IFSLIT_OBJS) -lm

test-cas: $(CAS_TEST_TARGET)
	./$(CAS_TEST_TARGET)

$(CAS_TEST_TARGET): test/cas_bridge_smoke.cpp $(SRC_DIR)/cas_bridge.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -I$(SRC_DIR) test/cas_bridge_smoke.cpp -o $@

test-hopfl1: $(HOPFL1_TEST_TARGET)
	./$(HOPFL1_TEST_TARGET)

$(HOPFL1_TEST_TARGET): test/hopf_l1_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/hopf_l1_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

FOLDNF_TEST_TARGET := $(BUILD_DIR)/fold_nf_smoke$(EXEEXT)
test-foldnf: $(FOLDNF_TEST_TARGET)
	./$(FOLDNF_TEST_TARGET)

$(FOLDNF_TEST_TARGET): test/fold_nf_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/fold_nf_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

CODIM2_TEST_TARGET := $(BUILD_DIR)/codim2_smoke$(EXEEXT)
test-codim2: $(CODIM2_TEST_TARGET)
	./$(CODIM2_TEST_TARGET)

$(CODIM2_TEST_TARGET): test/codim2_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/codim2_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

TWOPARAM_TEST_TARGET := $(BUILD_DIR)/twoparam_smoke$(EXEEXT)
test-twoparam: $(TWOPARAM_TEST_TARGET)
	./$(TWOPARAM_TEST_TARGET)

$(TWOPARAM_TEST_TARGET): test/twoparam_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/twoparam_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

LCCOLLOC_TEST_TARGET := $(BUILD_DIR)/lc_colloc_smoke$(EXEEXT)
test-lccolloc: $(LCCOLLOC_TEST_TARGET)
	./$(LCCOLLOC_TEST_TARGET)

$(LCCOLLOC_TEST_TARGET): test/lc_colloc_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/lc_colloc_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

TPC2_TEST_TARGET := $(BUILD_DIR)/twoparam_codim2_smoke$(EXEEXT)
test-tpc2: $(TPC2_TEST_TARGET)
	./$(TPC2_TEST_TARGET)

$(TPC2_TEST_TARGET): test/twoparam_codim2_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/twoparam_codim2_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

LPC_TEST_TARGET := $(BUILD_DIR)/lpc_smoke$(EXEEXT)
test-lpc: $(LPC_TEST_TARGET)
	./$(LPC_TEST_TARGET)

$(LPC_TEST_TARGET): test/lpc_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/lpc_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

LPCARC_TEST_TARGET := $(BUILD_DIR)/lpc_arclength_smoke$(EXEEXT)
test-lpcarc: $(LPCARC_TEST_TARGET)
	./$(LPCARC_TEST_TARGET)

$(LPCARC_TEST_TARGET): test/lpc_arclength_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/lpc_arclength_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

BRANCHSW_TEST_TARGET := $(BUILD_DIR)/branch_switch_smoke$(EXEEXT)
test-branchsw: $(BRANCHSW_TEST_TARGET)
	./$(BRANCHSW_TEST_TARGET)

$(BRANCHSW_TEST_TARGET): test/branch_switch_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/branch_switch_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

BTCODIM2_TEST_TARGET := $(BUILD_DIR)/bt_codim2_smoke$(EXEEXT)
test-btcodim2: $(BTCODIM2_TEST_TARGET)
	./$(BTCODIM2_TEST_TARGET)

$(BTCODIM2_TEST_TARGET): test/bt_codim2_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/bt_codim2_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

BASINSMT_TEST_TARGET := $(BUILD_DIR)/basins_mt_smoke$(EXEEXT)
test-basinsmt: $(BASINSMT_TEST_TARGET)
	./$(BASINSMT_TEST_TARGET)

$(BASINSMT_TEST_TARGET): test/basins_mt_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/basins_mt_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

HOMOCLINIC_TEST_TARGET := $(BUILD_DIR)/homoclinic_smoke$(EXEEXT)
test-homoclinic: $(HOMOCLINIC_TEST_TARGET)
	./$(HOMOCLINIC_TEST_TARGET)

$(HOMOCLINIC_TEST_TARGET): test/homoclinic_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/homoclinic_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

HOMOCONT_TEST_TARGET := $(BUILD_DIR)/homoclinic_cont_smoke$(EXEEXT)
test-homocont: $(HOMOCONT_TEST_TARGET)
	./$(HOMOCONT_TEST_TARGET)

$(HOMOCONT_TEST_TARGET): test/homoclinic_cont_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/homoclinic_cont_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

VALIDATION_TEST_TARGET := $(BUILD_DIR)/validation_smoke$(EXEEXT)
test-validation: $(VALIDATION_TEST_TARGET)
	./$(VALIDATION_TEST_TARGET)

$(VALIDATION_TEST_TARGET): test/validation_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/validation_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

CODIM2COEF_TEST_TARGET := $(BUILD_DIR)/codim2_coeffs_smoke$(EXEEXT)
test-codim2coef: $(CODIM2COEF_TEST_TARGET)
	./$(CODIM2COEF_TEST_TARGET)

$(CODIM2COEF_TEST_TARGET): test/codim2_coeffs_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/codim2_coeffs_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

HOMOSEED_TEST_TARGET := $(BUILD_DIR)/homoclinic_seed_smoke$(EXEEXT)
test-homoseed: $(HOMOSEED_TEST_TARGET)
	./$(HOMOSEED_TEST_TARGET)

$(HOMOSEED_TEST_TARGET): test/homoclinic_seed_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/homoclinic_seed_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

ZHHH_TEST_TARGET := $(BUILD_DIR)/zhhh_smoke$(EXEEXT)
test-zhhh: $(ZHHH_TEST_TARGET)
	./$(ZHHH_TEST_TARGET)

$(ZHHH_TEST_TARGET): test/zhhh_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/zhhh_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

LCSEED_TEST_TARGET := $(BUILD_DIR)/lc_selfseed_smoke$(EXEEXT)
test-lcseed: $(LCSEED_TEST_TARGET)
	./$(LCSEED_TEST_TARGET)

$(LCSEED_TEST_TARGET): test/lc_selfseed_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/lc_selfseed_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

PDNS_TEST_TARGET := $(BUILD_DIR)/cycle_pdns_smoke$(EXEEXT)
test-pdns: $(PDNS_TEST_TARGET)
	./$(PDNS_TEST_TARGET)

$(PDNS_TEST_TARGET): test/cycle_pdns_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/cycle_pdns_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

PDCURVE_TEST_TARGET := $(BUILD_DIR)/pd_curve_smoke$(EXEEXT)
test-pdcurve: $(PDCURVE_TEST_TARGET)
	./$(PDCURVE_TEST_TARGET)

$(PDCURVE_TEST_TARGET): test/pd_curve_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/pd_curve_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

HETERO_TEST_TARGET := $(BUILD_DIR)/heteroclinic_smoke$(EXEEXT)
test-hetero: $(HETERO_TEST_TARGET)
	./$(HETERO_TEST_TARGET)

$(HETERO_TEST_TARGET): test/heteroclinic_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/heteroclinic_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

BPC_TEST_TARGET := $(BUILD_DIR)/bpc_smoke$(EXEEXT)
test-bpc: $(BPC_TEST_TARGET)
	./$(BPC_TEST_TARGET)

$(BPC_TEST_TARGET): test/bpc_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/bpc_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

CODIM2CYC_TEST_TARGET := $(BUILD_DIR)/codim2_cycle_smoke$(EXEEXT)
test-codim2cyc: $(CODIM2CYC_TEST_TARGET)
	./$(CODIM2CYC_TEST_TARGET)

$(CODIM2CYC_TEST_TARGET): test/codim2_cycle_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/codim2_cycle_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

LINDIAG_TEST_TARGET := $(BUILD_DIR)/lin_diag_smoke$(EXEEXT)
test-lindiag: $(LINDIAG_TEST_TARGET)
	./$(LINDIAG_TEST_TARGET)

$(LINDIAG_TEST_TARGET): test/lin_diag_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/lin_diag_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

FINDHOMO_TEST_TARGET := $(BUILD_DIR)/find_homoclinic_smoke$(EXEEXT)
test-findhomo: $(FINDHOMO_TEST_TARGET)
	./$(FINDHOMO_TEST_TARGET)

$(FINDHOMO_TEST_TARGET): test/find_homoclinic_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/find_homoclinic_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm


LPCCURVE_TEST_TARGET := $(BUILD_DIR)/lpc_curve_smoke$(EXEEXT)
test-lpccurve: $(LPCCURVE_TEST_TARGET)
	./$(LPCCURVE_TEST_TARGET)

$(LPCCURVE_TEST_TARGET): test/lpc_curve_smoke.cpp $(SRC_DIR)/analysis.cpp $(SRC_DIR)/analysis.h
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 -pthread -I$(SRC_DIR) test/lpc_curve_smoke.cpp $(SRC_DIR)/analysis.cpp -o $@ -lm

BRIDGEALIGN_TEST_TARGET := $(BUILD_DIR)/bridge_align_smoke$(EXEEXT)
test-bridgealign: $(BRIDGEALIGN_TEST_TARGET)
	./$(BRIDGEALIGN_TEST_TARGET)

PROJSOLID_TEST_TARGET := $(BUILD_DIR)/projsolid_smoke$(EXEEXT)
test-projsolid: $(PROJSOLID_TEST_TARGET)
	./$(PROJSOLID_TEST_TARGET)

$(PROJSOLID_TEST_TARGET): test/projsolid_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/projsolid_smoke.cpp -o $@ -lm

$(BRIDGEALIGN_TEST_TARGET): test/bridge_align_smoke.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXSTD) $(WARNINGS_CXX) -O2 test/bridge_align_smoke.cpp -o $@ -lm


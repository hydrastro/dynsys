# dynsys — GLFW/OpenGL dynamical-system visualizer with Dear ImGui + TPCAS
#
# Normal workflow:
#   nix develop
#   make
#   make run

CC         ?= cc
CXX        ?= c++
AR         ?= ar
RM         ?= rm -f
RMDIR      ?= rm -rf
MKDIR_P    ?= mkdir -p
INSTALL    ?= install
PKG_CONFIG ?= pkg-config

SRC_DIR    := src
TPCAS_DIR  := vendor/tpcas
DS_DIR     := $(TPCAS_DIR)/vendor/ds/lib
BUILD_DIR  := build
OBJ_DIR    := $(BUILD_DIR)/obj
C_OBJ_DIR  := $(BUILD_DIR)/obj-c
CXX_OBJ_DIR:= $(BUILD_DIR)/obj-cxx
DEP_DIR    := $(BUILD_DIR)/dep
C_DEP_DIR  := $(BUILD_DIR)/dep-c
CXX_DEP_DIR:= $(BUILD_DIR)/dep-cxx
TARGET     := $(BUILD_DIR)/dynsys

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

MODE ?= debug
PKGS ?= glfw3 glew cglm
IMGUI_DIR ?=

WARNINGS_C   := -Wall -Wextra -Wpedantic
WARNINGS_CXX := -Wall -Wextra
CSTD         := -std=c11
CXXSTD       := -std=c++17

PKG_CFLAGS := $(shell $(PKG_CONFIG) --cflags $(PKGS) 2>/dev/null)
PKG_LIBS   := $(shell $(PKG_CONFIG) --libs $(PKGS) 2>/dev/null)

CPPFLAGS ?=
CPPFLAGS += -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds $(PKG_CFLAGS)
ifneq ($(strip $(IMGUI_DIR)),)
  CPPFLAGS += -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends
endif

CFLAGS ?=
CFLAGS += $(CSTD) $(WARNINGS_C)
CXXFLAGS ?=
CXXFLAGS += $(CXXSTD) $(WARNINGS_CXX)
LDFLAGS ?=
LDLIBS ?=
LDLIBS += $(PKG_LIBS) -lGL -ldl -lm

ifeq ($(MODE),release)
  CFLAGS += -O2 -DNDEBUG
  CXXFLAGS += -O2 -DNDEBUG
else ifeq ($(MODE),asan)
  CFLAGS += -O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
  CXXFLAGS += -O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
  LDFLAGS += -fsanitize=address,undefined
else
  CFLAGS += -O0 -g3
  CXXFLAGS += -O0 -g3
endif

DYNSYS_CPP_SRCS := $(SRC_DIR)/dynsys.cpp
DYNSYS_OBJS := $(patsubst %.cpp,$(CXX_OBJ_DIR)/%.o,$(DYNSYS_CPP_SRCS))
DYNSYS_DEPS := $(patsubst %.cpp,$(CXX_DEP_DIR)/%.d,$(DYNSYS_CPP_SRCS))

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

IMGUI_CORE_SRCS := imgui.cpp imgui_draw.cpp imgui_tables.cpp imgui_widgets.cpp
IMGUI_BACKEND_SRCS := imgui_impl_glfw.cpp imgui_impl_opengl3.cpp
IMGUI_CORE_OBJS := $(addprefix $(CXX_OBJ_DIR)/imgui/,$(IMGUI_CORE_SRCS:.cpp=.o))
IMGUI_BACKEND_OBJS := $(addprefix $(CXX_OBJ_DIR)/imgui/backends/,$(IMGUI_BACKEND_SRCS:.cpp=.o))
IMGUI_OBJS := $(IMGUI_CORE_OBJS) $(IMGUI_BACKEND_OBJS)
IMGUI_DEPS := $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$(IMGUI_OBJS))

OBJS := $(DYNSYS_OBJS) $(C_OBJS) $(IMGUI_OBJS)
DEPS := $(DYNSYS_DEPS) $(C_DEPS) $(IMGUI_DEPS)

.PHONY: all check-deps check-legacy prune-legacy run debug release asan clean distclean install uninstall format print-vars help

all: check-deps check-legacy $(TARGET)

check-deps:
	@$(PKG_CONFIG) --exists $(PKGS) || { \
	  echo "error: missing pkg-config dependencies: $(PKGS)" >&2; \
	  echo "hint: run 'nix develop' first, then run 'make' again." >&2; \
	  exit 1; \
	}
	@if [ -z "$(strip $(IMGUI_DIR))" ]; then \
	  echo "error: IMGUI_DIR is not set" >&2; \
	  echo "hint: run 'nix develop' first; the flake exports IMGUI_DIR." >&2; \
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

$(CXX_OBJ_DIR)/imgui/%.o: $(IMGUI_DIR)/%.cpp
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@) -c $< -o $@

$(CXX_OBJ_DIR)/imgui/backends/%.o: $(IMGUI_DIR)/backends/%.cpp
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -MF $(patsubst $(CXX_OBJ_DIR)/%.o,$(CXX_DEP_DIR)/%.d,$@) -c $< -o $@

run: all
	./$(TARGET)

debug:
	$(MAKE) MODE=debug

release:
	$(MAKE) MODE=release

asan:
	$(MAKE) MODE=asan

install: all
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 0755 $(TARGET) $(DESTDIR)$(BINDIR)/dynsys

uninstall:
	$(RM) $(DESTDIR)$(BINDIR)/dynsys

clean:
	$(RMDIR) $(BUILD_DIR)

distclean: clean
	$(RM) result

format:
	clang-format -i $(SRC_DIR)/*.cpp

print-vars:
	@echo "CC=$(CC)"
	@echo "CXX=$(CXX)"
	@echo "MODE=$(MODE)"
	@echo "TARGET=$(TARGET)"
	@echo "C_OBJ_DIR=$(C_OBJ_DIR)"
	@echo "CXX_OBJ_DIR=$(CXX_OBJ_DIR)"
	@echo "PKGS=$(PKGS)"
	@echo "IMGUI_DIR=$(IMGUI_DIR)"
	@echo "PKG_CFLAGS=$(PKG_CFLAGS)"
	@echo "PKG_LIBS=$(PKG_LIBS)"

help:
	@echo "Targets:"
	@echo "  make             build $(TARGET)"
	@echo "  make run         build and run"
	@echo "  make release     optimized build"
	@echo "  make asan        Address/UBSan build"
	@echo "  make clean       remove build artifacts"
	@echo "  make prune-legacy dry-run removal of old C/FreeType tree artifacts"
	@echo "  make prune-legacy CLEAN_APPLY=1 apply legacy cleanup"
	@echo "  make format      clang-format src/*.cpp"

-include $(DEPS)

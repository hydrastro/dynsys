# dynsys — GLFW/OpenGL dynamical-system visualizer using TPCAS-mode equations

CC       ?= cc
AR       ?= ar
RM       ?= rm -f
RMDIR    ?= rm -rf
MKDIR_P  ?= mkdir -p
INSTALL  ?= install
PKG_CONFIG ?= pkg-config

SRC_DIR    := src
TPCAS_DIR  := vendor/tpcas
DS_DIR     := $(TPCAS_DIR)/vendor/ds/lib
BUILD_DIR  := build
OBJ_DIR    := $(BUILD_DIR)/obj
DEP_DIR    := $(BUILD_DIR)/dep
TARGET     := $(BUILD_DIR)/dynsys

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

MODE ?= debug
PKGS ?= glfw3 glew cglm freetype2
FONT_PATH ?= /usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf

WARNINGS := -Wall -Wextra -Wpedantic
STD      := -std=c11

PKG_CFLAGS := $(shell $(PKG_CONFIG) --cflags $(PKGS) 2>/dev/null)
PKG_LIBS   := $(shell $(PKG_CONFIG) --libs $(PKGS) 2>/dev/null)

CPPFLAGS ?=
CPPFLAGS += -I$(SRC_DIR) -I$(TPCAS_DIR)/src -I$(TPCAS_DIR)/vendor/ds $(PKG_CFLAGS)
CPPFLAGS += -DDYNSYS_DEFAULT_FONT_PATH='"$(FONT_PATH)"'
CFLAGS ?=
CFLAGS += $(STD) $(WARNINGS)
LDFLAGS ?=
LDLIBS ?=
LDLIBS += $(PKG_LIBS) -lGL -lm

ifeq ($(MODE),release)
  CFLAGS += -O2 -DNDEBUG
else ifeq ($(MODE),asan)
  CFLAGS += -O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined
  LDFLAGS += -fsanitize=address,undefined
else
  CFLAGS += -O0 -g3
endif

DYNSYS_SRCS := $(SRC_DIR)/dynsys.c
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

SRCS := $(DYNSYS_SRCS) $(TPCAS_SRCS) $(DS_SRCS)
OBJS := $(patsubst %.c,$(OBJ_DIR)/%.o,$(SRCS))
DEPS := $(patsubst %.c,$(DEP_DIR)/%.d,$(SRCS))

.PHONY: all run clean distclean install uninstall format print-vars

all: $(TARGET)

$(TARGET): $(OBJS)
	@$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(OBJ_DIR)/%.o: %.c
	@$(MKDIR_P) $(dir $@) $(dir $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@))
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

install: $(TARGET)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 0755 $(TARGET) $(DESTDIR)$(BINDIR)/dynsys

uninstall:
	$(RM) $(DESTDIR)$(BINDIR)/dynsys

clean:
	$(RMDIR) $(BUILD_DIR)

distclean: clean
	$(RM) result

format:
	clang-format -i $(SRC_DIR)/*.c

print-vars:
	@echo "CC=$(CC)"
	@echo "MODE=$(MODE)"
	@echo "TARGET=$(TARGET)"
	@echo "PKGS=$(PKGS)"
	@echo "FONT_PATH=$(FONT_PATH)"
	@echo "PKG_CFLAGS=$(PKG_CFLAGS)"
	@echo "PKG_LIBS=$(PKG_LIBS)"

-include $(DEPS)

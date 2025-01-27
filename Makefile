CC = gcc
AR = ar
CFLAGS = -lm -lgmp -lcglm -lglut -lGL -lGLEW -lglfw  -llizard -lds -lfreetype
CFLAGS += -O2 -march=native -mtune=native
#CFLAGS += -Wall -Wextra -Werror -pedantic -pedantic-errors
CFLAGS += -fno-common -Wl,--gc-sections -Wredundant-decls -Wno-unused-parameter
CFLAGS += -fstack-protector-strong -fPIE -fpie
CFLAGS += -D_FORTIFY_SOURCE=2 -Wformat -Wformat-security
CFLAGS += -fstack-clash-protection -z noexecstack -z relro -z now
#CFLAGS += -Wl,-z,relro,-z,now -Wl,-pie #-fpie
#CFLAGS += -Waggregate-return -Wbad-function-cast -Wcast-align -Wcast-qual -Wdeclaration-after-statement
#CLFAGS += -Wfloat-equal -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wmissing-prototypes -Wnested-externs
#CFLAGS += -Wpointer-arith -Wredundant-decls -Wsequence-point -Wstrict-prototypes
#CFLAGS += -Wswitch -Wundef -Wunreachable-code -Wwrite-strings -Wconversion

CFLAGS_DEBUG = -g -fno-omit-frame-pointer -fsanitize=address,undefined -fsanitize=leak

SRC = dynsys.c
OUT = dynsys

all: $(OUT)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(OUT)

rebuild: clean all

run: $(OUT)
	./$(OUT)

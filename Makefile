# -*- Makefile -*-

-include make.inc
CC ?= gcc
CFLAGS ?= -std=c11 -ggdb -O3 -fopenmp
LDFLAGS ?= -fopenmp
LDLIBS ?= -lgraphblas

OBJS = GrB-mxm-timer.o cmdline.o generator.o prng.o globals.o
ifndef TARGET_MWX
OBJS += hooks.o
endif

CPPFLAGS += -Irandom123/include
LDLIBS += -lm

ifdef TARGET_MWX
TARGET_EXECUTABLE = GrB-mxm-timer.mwx
else
TARGET_EXECUTABLE = GrB-mxm-timer
endif

all: $(TARGET_EXECUTABLE)

GrB-mxm-timer: $(OBJS)

%.mwx: %
	cp $< $@

cmdline.c cmdline.h &: cmdline.ggo
	gengetopt < $^

GrB-mxm-timer.o: GrB-mxm-timer.c globals.h generator.h prng.h
cmdline.o: cmdline.c
generator.o: generator.c globals.h prng.h compat.h
prng.o: prng.c prng.h globals.h
globals.o: globals.c globals.h
hooks.o: hooks.c hooks.h

# Ugly hack.
GrB-mxm-timer.mwx: GrB-mxm-timer

.PHONY: clean
clean:
	rm -f GrB-mxm-timer.mwx GrB-mxm-timer $(OBJS)

# -*- Makefile -*-

-include make.inc
CC ?= gcc
CFLAGS ?= -std=c11 -ggdb -O3 -fopenmp
CPPFLAGS ?= -Irandom123/include
LDFLAGS ?= -fopenmp
LDLIBS ?= -lgraphblas

OBJS = GrB-mxm-timer.o cmdline.o generator.o prng.o globals.o
ifndef TARGET_MWX
OBJS += hooks.o
endif

LDLIBS += -lm -lrt

GrB-mxm-timer:	$(OBJS)

cmdline.c cmdline.h &: cmdline.ggo
	gengetopt < $^

GrB-mxm-timer.o: GrB-mxm-timer.c globals.h generator.h prng.h
cmdline.o: cmdline.c
generator.o: generator.c globals.h prng.h compat.h
prng.o: prng.c prng.h globals.h
globals.o: globals.c globals.h
hooks.o: hooks.c hooks.h

.PHONY: clean
clean:
	rm -f GrB-mxm-timer GrB-mxm-timer.o cmdline.o generator.o prng.o globals.o

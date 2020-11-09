CC = gcc
CFLAGS = -O3 -Wall -DDEBUG_LEVEL=$(DEBUG_LEVEL)
LFLAGS = -pthread -lm
DSFMTFLAGS = -DHAVE_SSE2 -DDSFMT_MEXP=521
INCDIR = $(ROOTINC)/include/ $(ROOTINC)/external/
INCDIRPARAMS = $(INCDIR:%=-I%)
OBJ = debug.o linked_list.o particle.o simulate.o count_cores.o output.o simulate_rdme.o simulate_threads.o model.o pthread_barrier.o
#OBJ = linked_list.o particle.o simulate.o count_cores.o output.o simulate_rdme.o binheap.o simulate_threads.o model.o pthread_barrier.o


.PHONY: all

all: ssa_sdpd

main.o:
	$(CC) -c -o main.o $(MODEL) $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

dSFMT.o:
	$(CC) -c -o dSFMT.o $(ROOTINC)/external/dSFMT/dSFMT.c $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

%.o: $(ROOT)/src/%.c
	$(CC) -c -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

ssa_sdpd: main.o dSFMT.o $(OBJ)
	$(CC)  -o ssa_sdpd $(OBJ)  main.o dSFMT.o $(LFLAGS)

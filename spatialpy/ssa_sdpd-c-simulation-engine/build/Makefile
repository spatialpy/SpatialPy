CC = g++
CFLAGS = -c -std=c++14 -Wall -O3 -g
SIMFLAGS = -l. -std=c++14 -Wall -O3
LFLAGS = -pthread -lm
INCDIR = $(ROOTINC)/include/ $(ROOTINC)/external/ $(ROOTINC)/external/ANN/include/
INCDIRPARAMS = $(INCDIR:%=-I%)
OBJ = NRMConstant_v5.o particle.o simulate.o count_cores.o output.o simulate_rdme.o simulate_threads.o model.o pthread_barrier.o
ANNOBJECTS = ANN.o brute.o kd_tree.o kd_util.o kd_split.o \
	kd_dump.o kd_search.o kd_pr_search.o kd_fix_rad_search.o \
	bd_tree.o bd_search.o bd_pr_search.o bd_fix_rad_search.o \
	perf.o

.PHONY: all

all: ssa_sdpd

ANN:
	make -f $(ROOTINC)/external/ANN/src/Makefile.spatialpy ROOTINC=$(ROOTINC)

main.o:
	$(CC) -c $(GPROFFLAG) $(GDB_FLAG) -o main.o $(MODEL) $(INCDIRPARAMS) $(CFLAGS)

%.o: $(ROOT)/src/%.cpp
	$(CC) -c $(GPROFFLAG) $(GDB_FLAG) -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS)

%.o: $(ROOT)/src/%.c
	$(CC) -c $(GPROFFLAG) $(GDB_FLAG) -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS)

ssa_sdpd: ANN main.o $(OBJ)
	$(CC) $(GPROFFLAG) $(GDB_FLAG) -o ssa_sdpd $(OBJ)  main.o $(ANNOBJECTS) $(LFLAGS)


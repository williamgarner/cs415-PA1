## On Bridges we will check versus your performance versus Intel MKL library's BLAS.
#
#CC = cc
#OPT = -O3
#CFLAGS = -Wall -std=gnu99 $(OPT)
##MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
##LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
#LDLIBS = -lrt  -I$(MKLROOT)/include -Wl,-L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
#
#targets = benchmark-naive benchmark-blocked benchmark-blas
#objects = benchmark.o dgemm-naive.o dgemm-blocked.o dgemm-blas.o
#
#.PHONY : default
#default : all
#
#.PHONY : all
#all : clean $(targets)
#
#benchmark-naive : benchmark.o dgemm-naive.o
#	$(CC) -o $@ $^ $(LDLIBS)
#benchmark-blocked : benchmark.o dgemm-blocked.o
#	$(CC) -o $@ $^ $(LDLIBS)
#benchmark-blas : benchmark.o dgemm-blas.o
#	$(CC) -o $@ $^ $(LDLIBS)
#
#%.o : %.c
#	$(CC) -c $(CFLAGS) $<
#
#.PHONY : clean
#clean:
#	rm -f $(targets) $(objects) *.stdout


# on Hopper, we will benchmark you against Cray LibSci, the default vendor-tuned BLAS. The Cray compiler wrappers handle all the linking. If you wish to compare with other BLAS implementations, check the NERSC documentation.
# This makefile is intended for the GNU C compiler. On Hopper, the Portland compilers are default, so you must instruct the Cray compiler wrappers to switch to GNU: type "module swap PrgEnv-pgi PrgEnv-gnu"

CC = cc
OPT = -Ofast -ffast-math -funroll-loops
CFLAGS = -Wall -std=gnu99 -march=opteron -msse -msse2 -msse3 $(OPT)
LDFLAGS = -Wall
# librt is needed for clock_gettime
LDLIBS = -lrt

targets = benchmark-naive benchmark-blocked benchmark-blas
objects = benchmark.o dgemm-naive.o dgemm-blocked.o dgemm-blas.o

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-naive : benchmark.o dgemm-naive.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects)
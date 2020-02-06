# On Bridges we will check versus your performance versus Intel MKL library's BLAS.

CC = icc
OPT = -Ofast -guide
CFLAGS = -Wall -funroll-all-loops -global-hoist -march=haswell -qopt-prefetch=5 -unroll-aggressive -vec -no-multibyte-chars -std=gnu99 $(OPT)
#MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
#LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
LDLIBS = -lrt  -I$(MKLROOT)/include -Wl,-L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

targets = benchmark-naive
objects = benchmark.o dgemm-naive.o

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-naive : benchmark.o dgemm-naive.o
	$(CC) -o $@ $^ $(LDLIBS)
#benchmark-blocked : benchmark.o dgemm-blocked.o
#	$(CC) -o $@ $^ $(LDLIBS)
#benchmark-blas : benchmark.o dgemm-blas.o
#	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects) *.stdout

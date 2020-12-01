CFLAGS = -g -Wall -Wstrict-prototypes 
PROGS = make-matrix print-matrix mm-serial mm-parallel
LDFLAGS = -lm 
CC = gcc
MCC = mpicc

all: $(PROGS)

make-matrix: make-matrix.o 

print-matrix: print-matrix.o 

mm-serial: mm-serial.o 

mm-parallel: mm-parallel.o matrix_checkerboard_io.o MyMPI.o
	$(MCC) $(LDFLAGS) -o $@ $+

clean:
	rm -f $(PROGS) *.o core*

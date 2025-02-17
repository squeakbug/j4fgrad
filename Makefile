CC = mpicc

all: bench

bench: libmatrix-bare
	$(CC) $(CFLAGS) -I include -o allgather.elf bin/bench.c libmatrix.o

libmatrix-bare: src/matrix.c src/bare/*.c
	$(CC) $(CFLAGS) -c @^ -I include libmatrix.o

libmatrix-mpi:
	$(CC) $(CFLAGS) -I include main.c -o allgather.elf

libmatrix-omp:
	$(CC) $(CFLAGS) -I include main.c -o allgather.elf

run:
	mpirun -n 4 allgather.elf

.PHONY: clean
	clean: rm *.elf

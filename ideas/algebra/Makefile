# ideas/algebra: KP-3 (de Bruijn adjacency), KP-14 (spaced seeds), root assignment.
CC = gcc
CXX = g++
OPT = -O3 -march=native
CFLAGS = $(OPT) -Wall -Wextra
HASH = ../../hash
NTHASH ?= /home/rcb112/project/claude/bench/data/kmer_poly_work/ntHash

all: test_spaced bench_spaced

kp7.o: $(HASH)/kp7.c $(HASH)/kp7.h
	$(CC) $(CFLAGS) -c $(HASH)/kp7.c -o $@

spaced.o: spaced.c spaced.h field61.h
	$(CC) $(CFLAGS) -c spaced.c -o $@

test_spaced: test_spaced.c spaced.o kp7.o
	$(CC) $(CFLAGS) test_spaced.c spaced.o kp7.o -o $@

test: test_spaced
	./test_spaced

mutation: test_spaced.c spaced.c spaced.h kp7.o
	$(CC) $(CFLAGS) -DSP_MUTATE test_spaced.c spaced.c kp7.o -o test_spaced_mutant
	@if ./test_spaced_mutant > mutant.out; then echo "SPACED MUTANT SURVIVED"; exit 1; else echo "spaced mutant killed: $$(tail -1 mutant.out)"; fi

bench_spaced: bench_spaced.cpp spaced.o kp7.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include bench_spaced.cpp spaced.o kp7.o $(NTHASH)/build_o3/libnthash.a -o $@

clean:
	rm -f *.o test_spaced test_spaced_mutant mutant.out bench_spaced dbg eisen

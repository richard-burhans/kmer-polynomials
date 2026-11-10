# Same flags for kp7 and the ntHash baseline (see bench.cpp); recorded in results/.
CC = gcc
CXX = g++
OPT = -O3 -march=native
CFLAGS = $(OPT) -Wall -Wextra
NTHASH ?= /home/rcb112/project/claude/bench/data/kmer_poly_work/ntHash

all: test_kp7 bench

kp7.o: kp7.c kp7.h
	$(CC) $(CFLAGS) -c kp7.c -o $@

test_kp7: test_kp7.c kp7.o
	$(CC) $(CFLAGS) test_kp7.c kp7.o -o $@

test: test_kp7
	./test_kp7

# break one rolling term; the tests must fail
mutation: test_kp7.c kp7.c kp7.h
	$(CC) $(CFLAGS) -DKP7_MUTATE test_kp7.c kp7.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "MUTANT SURVIVED"; exit 1; else echo "mutant killed: $$(tail -1 mutant.out)"; fi
	$(CC) $(CFLAGS) -DKP7_MUTATE_K test_kp7.c kp7.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "KARATSUBA MUTANT SURVIVED"; exit 1; else echo "karatsuba mutant killed: $$(tail -1 mutant.out)"; fi
	$(CC) $(CFLAGS) -DKP7_MUTATE_OPT test_kp7.c kp7.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "OPT MUTANT SURVIVED"; exit 1; else echo "opt mutant killed: $$(tail -1 mutant.out)"; fi
	$(CC) $(CFLAGS) -DKP7_MUTATE_SEQ test_kp7.c kp7.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "SEQ MUTANT SURVIVED"; exit 1; else echo "seq mutant killed: $$(tail -1 mutant.out)"; fi
	$(CC) $(CFLAGS) -DKP7_MUTATE_GATE test_kp7.c kp7.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "GATE MUTANT SURVIVED"; exit 1; else echo "gate mutant killed: $$(tail -1 mutant.out)"; fi

bench: bench.cpp kp7.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include bench.cpp kp7.o $(NTHASH)/build_o3/libnthash.a -o $@

clean:
	rm -f *.o test_kp7 test_kp7_mutant mutant.out bench

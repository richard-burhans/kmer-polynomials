# ideas/eisenstein: the Eisenstein revisit of KP-3 (de Bruijn) and KP-14 (spaced seeds). See results/ideas_eisenstein_*.md.
CC = gcc
CXX = g++
OPT = -O3 -march=native
CFLAGS = $(OPT) -Wall -Wextra
HASH = ../../hash
NTHASH ?= /home/rcb112/project/claude/bench/data/kmer_poly_work/ntHash

all: dbg_e31 dbg_e47 dbg_e63 test_spaced_e bench_spaced_e

kp7e.o: $(HASH)/kp7e.c $(HASH)/kp7e.h
	$(CC) $(CFLAGS) -c $(HASH)/kp7e.c -o $@

dbg_e31 dbg_e47 dbg_e63: dbg_e.c kp7e.o
	$(CC) $(CFLAGS) -DDBG_K=$(subst dbg_e,,$@) dbg_e.c kp7e.o -o $@

spaced_e.o: spaced_e.c spaced_e.h
	$(CC) $(CFLAGS) -c spaced_e.c -o $@

kp7.o: $(HASH)/kp7.c $(HASH)/kp7.h
	$(CC) $(CFLAGS) -c $(HASH)/kp7.c -o $@

test_spaced_e: test_spaced_e.c spaced_e.o kp7e.o
	$(CC) $(CFLAGS) test_spaced_e.c spaced_e.o kp7e.o -o $@

test: test_spaced_e
	./test_spaced_e

mutation: test_spaced_e.c spaced_e.c spaced_e.h kp7e.o
	$(CC) $(CFLAGS) -DSPE_MUTATE test_spaced_e.c spaced_e.c kp7e.o -o test_spaced_e_mutant
	@if ./test_spaced_e_mutant > mutant.out; then echo "SPACED_E MUTANT SURVIVED"; exit 1; else echo "spaced_e mutant killed: $$(tail -1 mutant.out)"; fi

# Eisenstein spaced seeds against ntHash2's SeedNtHash and the Gaussian spaced roller (ideas/algebra)
bench_spaced_e: bench_spaced_e.cpp spaced_e.o kp7e.o ../algebra/spaced.c ../algebra/spaced.h kp7.o
	$(CC) $(CFLAGS) -c ../algebra/spaced.c -o spaced_g.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include bench_spaced_e.cpp spaced_e.o spaced_g.o kp7e.o kp7.o $(NTHASH)/build_o3/libnthash.a -o $@

clean:
	rm -f *.o dbg_e31 dbg_e47 dbg_e63 test_spaced_e test_spaced_e_mutant bench_spaced_e mutant.out

# direct h(w) == h(rc w) test of SeedNtHash and the Eisenstein S_m, palindromic and non-palindromic masks (round 1)
seedrc: seedrc.cpp spaced_e.o kp7e.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include seedrc.cpp spaced_e.o kp7e.o $(NTHASH)/build_o3/libnthash.a -o $@

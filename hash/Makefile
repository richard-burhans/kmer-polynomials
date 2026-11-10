# Same flags for kp7 and the ntHash baseline (see bench.cpp); recorded in results/.
CC = gcc
CXX = g++
OPT = -O3 -march=native
CFLAGS = $(OPT) -Wall -Wextra
NTHASH ?= /home/rcb112/project/claude/bench/data/kmer_poly_work/ntHash

all: test_kp7 test_kp7m31 bench

kp7.o: kp7.c kp7.h
	$(CC) $(CFLAGS) -c kp7.c -o $@

test_kp7: test_kp7.c kp7.o
	$(CC) $(CFLAGS) test_kp7.c kp7.o -o $@

# CM31 (p = 2^31 - 1) scalar / AVX2 / AVX-512; the AVX-512 path is built only when -march targets it
kp7m31.o: kp7m31.c kp7m31.h
	$(CC) $(CFLAGS) -c kp7m31.c -o $@

test_kp7m31: test_kp7m31.c kp7m31.o
	$(CC) $(CFLAGS) test_kp7m31.c kp7m31.o -o $@

test: test_kp7 test_kp7m31
	./test_kp7
	./test_kp7m31

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
	$(CC) $(CFLAGS) -DKP7M31_MUTATE_ODD test_kp7m31.c kp7m31.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "M31 ODD-LANE MUTANT SURVIVED"; exit 1; else echo "m31 odd-lane mutant killed: $$(tail -1 mutant.out)"; fi
	$(CC) $(CFLAGS) -DKP7M31_MUTATE_BLEND test_kp7m31.c kp7m31.c -o test_kp7_mutant
	@if ./test_kp7_mutant > mutant.out; then echo "M31 BLEND MUTANT SURVIVED"; exit 1; else echo "m31 blend mutant killed: $$(tail -1 mutant.out)"; fi

bench: bench.cpp kp7.o kp7m31.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include bench.cpp kp7.o kp7m31.o $(NTHASH)/build_o3/libnthash.a -o $@

clean:
	rm -f *.o test_kp7 test_kp7m31 test_kp7_mutant mutant.out bench

# ---- hybrid v2 (kp7pack: packed-code gate, provable KP-7 hash for kept k-mers) ----
kp7pack.o: kp7pack.c kp7pack.h kp7.h kp7m31.h
	$(CC) $(CFLAGS) -c kp7pack.c -o $@

test_kp7pack: test_kp7pack.c kp7pack.o kp7.o kp7m31.o
	$(CC) $(CFLAGS) test_kp7pack.c kp7pack.o kp7.o kp7m31.o -o $@

test_pack: test_kp7pack
	./test_kp7pack

mutation_pack: test_kp7pack.c kp7pack.c kp7pack.h kp7.c kp7.h kp7m31.o
	$(CC) $(CFLAGS) -DKP7PACK_MUTATE_GATE test_kp7pack.c kp7pack.c kp7.c kp7m31.o -o test_kp7pack_mutant
	@if ./test_kp7pack_mutant > mutant_pack.out; then echo "PACK GATE MUTANT SURVIVED"; exit 1; else echo "pack gate mutant killed: $$(tail -1 mutant_pack.out)"; fi
	$(CC) $(CFLAGS) -DKP7PACK_MUTATE_RC test_kp7pack.c kp7pack.c kp7.c kp7m31.o -o test_kp7pack_mutant
	@if ./test_kp7pack_mutant > mutant_pack.out; then echo "PACK RC MUTANT SURVIVED"; exit 1; else echo "pack rc mutant killed: $$(tail -1 mutant_pack.out)"; fi
	$(CC) $(CFLAGS) -DKP7PACK_MUTATE_CHUNK test_kp7pack.c kp7pack.c kp7.c kp7m31.o -o test_kp7pack_mutant
	@if ./test_kp7pack_mutant > mutant_pack.out; then echo "PACK CHUNK MUTANT SURVIVED"; exit 1; else echo "pack chunk mutant killed: $$(tail -1 mutant_pack.out)"; fi
	$(CC) $(CFLAGS) -DKP7PACK_MUTATE_LAM test_kp7pack.c kp7pack.c kp7.c kp7m31.o -o test_kp7pack_mutant
	@if ./test_kp7pack_mutant > mutant_pack.out; then echo "PACK LAM MUTANT SURVIVED"; exit 1; else echo "pack lam mutant killed: $$(tail -1 mutant_pack.out)"; fi

bench_pack: bench_pack.cpp kp7pack.o kp7.o kp7m31.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include bench_pack.cpp kp7pack.o kp7.o kp7m31.o $(NTHASH)/build_o3/libnthash.a -o $@

# multi-seed sampling statistics (KP-7 next steps, item 3)
stats_seeds: stats_seeds.cpp kp7pack.o kp7.o kp7m31.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include stats_seeds.cpp kp7pack.o kp7.o kp7m31.o $(NTHASH)/build_o3/libnthash.a -o $@

# ---- Eisenstein KP-7 (kp7e: A=1, C=w, G=-w, T=-1 with w in F_q, q = 2^61 - 1; one field element per hash) ----
kp7e.o: kp7e.c kp7e.h
	$(CC) $(CFLAGS) -c kp7e.c -o $@

test_kp7e: test_kp7e.c kp7e.o kp7pack.o kp7.o kp7m31.o
	$(CC) $(CFLAGS) test_kp7e.c kp7e.o kp7pack.o kp7.o kp7m31.o -lm -o $@

test_e: test_kp7e
	./test_kp7e

mutation_e: test_kp7e.c kp7e.c kp7e.h kp7pack.o kp7.o kp7m31.o
	@for m in ROLL LANE TAB; do \
	  $(CC) $(CFLAGS) -DKP7E_MUTATE_$$m test_kp7e.c kp7e.c kp7pack.o kp7.o kp7m31.o -lm -o test_kp7e_mutant || exit 1; \
	  if ./test_kp7e_mutant > mutant_e.out; then echo "E $$m MUTANT SURVIVED"; exit 1; else echo "e $$m mutant killed: $$(tail -1 mutant_e.out)"; fi; \
	done

bench_e: bench_e.cpp kp7e.o kp7pack.o kp7.o kp7m31.o
	$(CXX) $(OPT) -std=c++17 -I$(NTHASH)/include bench_e.cpp kp7e.o kp7pack.o kp7.o kp7m31.o $(NTHASH)/build_o3/libnthash.a -o $@

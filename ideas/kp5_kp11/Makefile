# ideas/kp5_kp11: KP-5 (RC palindromes / inverted repeats as a polynomial identity) and KP-11 (generating functions;
# the useful form is mismatch-neighbourhood counting by convolution over (Z_4)^k). See results/ideas_kp5_kp11_*.md.
CC = gcc
OPT = -O3 -march=native
CFLAGS = $(OPT) -Wall -Wextra
HASH = ../../hash

all: ir mismatch

kp7e.o: $(HASH)/kp7e.c $(HASH)/kp7e.h
	$(CC) $(CFLAGS) -c $(HASH)/kp7e.c -o $@

ir: ir.c fasta.h kp7e.o
	$(CC) $(CFLAGS) ir.c kp7e.o -o $@

mismatch: mismatch.c fasta.h
	$(CC) $(CFLAGS) mismatch.c -o $@

test: ir mismatch
	./test.sh

# each mutant must make test.sh fail
mutation: ir.c mismatch.c fasta.h kp7e.o
	$(CC) $(CFLAGS) -DKP5_MUTATE ir.c kp7e.o -o ir_mut && $(CC) $(CFLAGS) mismatch.c -o mismatch_ok
	@cp ir ir.keep; cp ir_mut ir; if ./test.sh >/dev/null 2>&1; then cp ir.keep ir; echo "KP5 MUTANT SURVIVED"; exit 1; fi; cp ir.keep ir; echo "kp5 mutant killed"
	$(CC) $(CFLAGS) -DKP11_MUTATE mismatch.c -o mismatch_mut
	@cp mismatch mismatch.keep; cp mismatch_mut mismatch; if ./test.sh >/dev/null 2>&1; then cp mismatch.keep mismatch; echo "KP11 MUTANT SURVIVED"; exit 1; fi; cp mismatch.keep mismatch; echo "kp11 mutant killed"
	@rm -f ir_mut mismatch_mut mismatch_ok ir.keep mismatch.keep

clean:
	rm -f *.o ir mismatch

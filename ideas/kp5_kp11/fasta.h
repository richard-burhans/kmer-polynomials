/* Tiny FASTA loader shared by ir.c and mismatch.c: one record (by name, or the first), upper-cased, in one buffer. */
#ifndef KP511_FASTA_H
#define KP511_FASTA_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

static char* load_record(const char* path, const char* want, size_t* len_out)
{
  FILE* f = fopen(path, "r");
  if (!f) { perror(path); exit(2); }
  size_t cap = 1 << 26, n = 0;
  char* s = malloc(cap);
  char line[1 << 16];
  int on = 0, done = 0;
  while (!done && fgets(line, sizeof line, f)) {
    if (line[0] == '>') {
      if (on) { done = 1; break; }
      char name[4096];
      sscanf(line + 1, "%4095s", name);
      on = !want || !strcmp(name, want);
      continue;
    }
    if (!on) continue;
    for (char* c = line; *c && *c != '\n'; c++) {
      if (n + 1 >= cap) { cap *= 2; s = realloc(s, cap); }
      s[n++] = (char)toupper((unsigned char)*c);
    }
  }
  fclose(f);
  if (!n) { fprintf(stderr, "record %s not found in %s\n", want ? want : "(first)", path); exit(2); }
  s[n] = 0;
  *len_out = n;
  return s;
}

static double now_s(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}
#endif

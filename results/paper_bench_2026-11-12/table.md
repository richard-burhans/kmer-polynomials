Load baselines (median of each binary's __load_only__ runs): bench 7.69e9 cycles (n = 6, spread 2.7%); bench_e 6.98e9 cycles (n = 6, spread 1.2%); bench_kept 6.99e9 cycles (n = 6, spread 0.8%)

| method | binary | cycles/k-mer (median) | spread | min | instructions/k-mer | branch-miss rate | n | vs batched ntHash (same binary) |
|---|---|---|---|---|---|---|---|---|
| ntHash batched (ntbatch) | bench_e | 11.37 | 0.4% | 11.34 | 47.0 | 0.34% | 5 | 1.00× |
| ntHash library roll | bench_e | 23.92 | 1.5% | 23.90 | 110.1 | 0.23% | 5 | 0.48× |
| e sym_seq L1 | bench_e | 27.85 | 1.6% | 27.74 | 72.6 | 0.32% | 5 | 0.41× |
| e sym_seq L2 | bench_e | 24.67 | 2.8% | 24.44 | 78.5 | 0.33% | 5 | 0.46× |
| gauss61 sym_seq L2 | bench_e | 47.19 | 0.6% | 47.12 | 179.3 | 0.33% | 5 | 0.24× |
| e branchy 1/1000 | bench_e | 9.43 | 1.5% | 9.38 | 43.0 | 0.35% | 5 | 1.21× |
| e branchy 1/100 | bench_e | 9.94 | 0.7% | 9.92 | 43.8 | 0.42% | 5 | 1.14× |
| e batched 1/16 | bench_e | 11.74 | 3.8% | 11.65 | 51.9 | 0.33% | 5 | 0.97× |
| e batched 1/4 | bench_e | 17.12 | 1.5% | 17.04 | 70.1 | 0.33% | 5 | 0.66× |
| m31 avx2 8-lane | bench | 26.44 | 0.6% | 26.35 | 61.7 | 0.32% | 5 | 0.43× |
| m31 avx512 16-lane | bench | 25.94 | 1.3% | 25.77 | 60.2 | 0.34% | 5 | 0.44× |
| ntHash batched | bench_kept | 11.43 | 1.8% | 11.28 | 48.0 | 0.34% | 5 | 1.00× |
| uh all | bench_kept | 10.24 | 3.7% | 10.14 | 46.9 | 0.25% | 5 | 1.12× |
| lib-e 8 br 1/1000 | bench_kept | 9.39 | 4.3% | 9.35 | 43.0 | 0.35% | 5 | 1.22× |
| lib-e 8 br 1/100 | bench_kept | 9.90 | 43.0% | 9.88 | 43.8 | 0.42% | 5 | 1.15× |
| lib-e 8 bat 1/1000 | bench_kept | 9.93 | 0.4% | 9.92 | 46.0 | 0.34% | 5 | 1.15× |
| lib-e 8 bat 1/100 | bench_kept | 10.19 | 2.2% | 10.09 | 46.8 | 0.34% | 5 | 1.12× |
| lib-e 8 bat 1/16 | bench_kept | 11.21 | 1.1% | 11.21 | 51.6 | 0.34% | 5 | 1.02× |
| lib-e 8 bat 1/4 | bench_kept | 15.43 | 2.0% | 15.28 | 68.9 | 0.33% | 5 | 0.74× |
| ec 13 br 1/1000 | bench_kept | 9.15 | 0.9% | 9.09 | 40.0 | 0.35% | 5 | 1.25× |
| ec 13 br 1/100 | bench_kept | 9.65 | 0.3% | 9.64 | 40.7 | 0.42% | 5 | 1.18× |
| nt 8 br 1/1000 | bench_kept | 9.72 | 0.5% | 9.70 | 42.0 | 0.35% | 5 | 1.18× |
| nt 8 br 1/100 | bench_kept | 10.23 | 0.5% | 10.20 | 42.8 | 0.42% | 5 | 1.12× |
| nt 8 bat 1/16 | bench_kept | 11.08 | 1.8% | 10.96 | 50.4 | 0.34% | 5 | 1.03× |
| nt 8 bat 1/4 | bench_kept | 14.70 | 2.4% | 14.68 | 67.1 | 0.33% | 5 | 0.78× |
| uh 0 br 1/1000 | bench_kept | 9.15 | 1.2% | 9.11 | 39.9 | 0.35% | 5 | 1.25× |
| uh 0 br 1/100 | bench_kept | 9.52 | 1.8% | 9.48 | 40.1 | 0.42% | 5 | 1.20× |
| uh 0 bat 1/16 | bench_kept | 10.25 | 1.8% | 10.15 | 48.3 | 0.34% | 5 | 1.12× |
| uh 0 bat 1/4 | bench_kept | 11.36 | 1.2% | 11.31 | 52.7 | 0.33% | 5 | 1.01× |

- bench_e `e batched`: 28.7 cycles per kept k-mer, intercept 9.95; beats batched ntHash (11.37) below 1/20.2
- bench_kept `lib-e 8 bat`: 22.5 cycles per kept k-mer, intercept 9.80; beats batched ntHash (11.43) below 1/13.8
- bench_kept `nt 8 bat`: 19.3 cycles per kept k-mer, intercept 9.87; beats batched ntHash (11.43) below 1/12.3
- bench_kept `uh 0 bat`: 5.9 cycles per kept k-mer, intercept 9.88; beats batched ntHash (11.43) below 1/3.8

Spaced seeds (run_165, its own wall-clock medians):
```
bases 60000000, rounds 5
masks 1: Gaussian 42088079 hashes 3.933 s (10.7 M/s)  Eisenstein 42088079 1.933 s (21.8 M/s)  SeedNtHash 42667269 1.461 s (29.2 M/s)  Gaussian/nt 0.37x  Eisenstein/nt 0.76x  (sinks 4b17e88121b9f673 1fab40382d7699e8 31e4636cdf9dcd4b)
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
masks 2: Gaussian 84176158 hashes 8.576 s (9.8 M/s)  Eisenstein 84176158 4.145 s (20.3 M/s)  SeedNtHash 85334538 3.016 s (28.3 M/s)  Gaussian/nt 0.35x  Eisenstein/nt 0.73x  (sinks 40ade31507deae24 13448dd61b276347 a24e9c470ef16a29)
SNP rate 0.001 mask 0: windows 42088079, unchanged: Gaussian 0.97632  Eisenstein 0.97632  SeedNtHash 0.97632  exact(care positions of both strands) 0.97632
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
SNP rate 0.001 mask 1: windows 42088079, unchanged: Gaussian 0.97146  Eisenstein 0.97146  SeedNtHash 0.98023  exact(care positions of both strands) 0.97146
SNP rate 0.010 mask 0: windows 42088079, unchanged: Gaussian 0.78573  Eisenstein 0.78573  SeedNtHash 0.78573  exact(care positions of both strands) 0.78573
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
SNP rate 0.010 mask 1: windows 42088079, unchanged: Gaussian 0.74722  Eisenstein 0.74722  SeedNtHash 0.81794  exact(care positions of both strands) 0.74722
SNP rate 0.050 mask 0: windows 42088079, unchanged: Gaussian 0.29184  Eisenstein 0.29184  SeedNtHash 0.29184  exact(care positions of both strands) 0.29184
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
SNP rate 0.050 mask 1: windows 42088079, unchanged: Gaussian 0.22567  Eisenstein 0.22567  SeedNtHash 0.35839  exact(care positions of both strands) 0.22567
```

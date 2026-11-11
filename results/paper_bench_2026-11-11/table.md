Load baselines (mean of each binary's __load_only__ runs): bench 7.65e9 cycles, 20.75e9 instructions (n = 3, spread 0.1%); bench_e 6.98e9 cycles, 20.75e9 instructions (n = 3, spread 0.7%)

| method | binary | cycles/k-mer | spread | instructions/k-mer | branch-miss rate | vs batched ntHash |
|---|---|---|---|---|---|---|
| ntHash batched (ntbatch) | bench_e | 11.37 | 0.7% | 47.0 | 0.34% | 1.00× |
| ntHash library roll | bench_e | 23.53 | 0.0% | 110.1 | 0.23% | 0.48× |
| e scalar roller (reference) | bench_e | 77.30 | 4.0% | 172.0 | 3.83% | 0.15× |
| e sym_seq L1 | bench_e | 27.82 | 0.3% | 72.6 | 0.32% | 0.41× |
| e sym_seq L2 | bench_e | 24.69 | 0.4% | 78.5 | 0.33% | 0.46× |
| gauss61 sym_seq L1 | bench_e | 50.00 | 0.0% | 160.2 | 0.31% | 0.23× |
| gauss61 sym_seq L2 | bench_e | 47.37 | 0.2% | 179.3 | 0.33% | 0.24× |
| m31 simd (avx512 else avx2) | bench_e | 25.99 | 0.6% | 60.2 | 0.37% | 0.44× |
| e branchy 1/1000 | bench_e | 9.40 | 0.2% | 43.0 | 0.35% | 1.21× |
| e branchy 1/100 | bench_e | 10.01 | 0.4% | 43.8 | 0.42% | 1.14× |
| e branchy 1/16 | bench_e | 13.32 | 0.7% | 48.8 | 0.86% | 0.85× |
| e branchy 1/4 | bench_e | 25.39 | 0.6% | 66.6 | 2.35% | 0.45× |
| e batched 1/1000 | bench_e | 9.95 | 0.5% | 46.0 | 0.34% | 1.14× |
| e batched 1/100 | bench_e | 10.22 | 0.3% | 46.9 | 0.34% | 1.11× |
| e batched 1/16 | bench_e | 11.74 | 0.8% | 51.9 | 0.34% | 0.97× |
| e batched 1/4 | bench_e | 17.08 | 0.3% | 70.1 | 0.33% | 0.67× |
| m31 branchy 1/1000 | bench_e | 9.52 | 0.5% | 43.1 | 0.35% | 1.19× |
| m31 branchy 1/100 | bench_e | 10.44 | 1.3% | 44.8 | 0.42% | 1.09× |
| m31 batched 1/16 | bench_e | 13.67 | 1.0% | 58.0 | 0.33% | 0.83× |
| m31 batched 1/4 | bench_e | 25.72 | 0.6% | 94.3 | 0.34% | 0.44× |
| m31 avx2 8-lane | bench | 26.43 | 0.3% | 61.7 | 0.32% | 0.43× |
| m31 avx512 16-lane | bench | 25.86 | 0.6% | 60.2 | 0.34% | 0.44× |

e batched: 28.5 cycles per kept k-mer, base 9.96; beats batched ntHash (11.37) below 1/20.3; beats e sym_seq L2 (24.69) below 1/1.9

m31 batched: 64.3 cycles per kept k-mer, base 9.66; beats batched ntHash (11.37) below 1/37.6; beats e sym_seq L2 (24.69) below 1/4.3

Spaced seeds (run_51, its own wall-clock rates):
```
bases 60000000, rounds 3
masks 1: Gaussian 42088079 hashes 3.792 s (11.1 M/s)  Eisenstein 42088079 1.953 s (21.5 M/s)  SeedNtHash 42667269 1.482 s (28.8 M/s)  Gaussian/nt 0.39x  Eisenstein/nt 0.76x  (sinks 98df426c5d3a5923 13d43d33e37c7ac9 31e4636cdf9dcd4b)
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
[ntHash::SeedNtHash] [33mWARNING: [0mSeed 1101100111010111011100101110101 is not symmetric, reverse-complement hashing will be inconsistent
masks 2: Gaussian 84176158 hashes 8.280 s (10.2 M/s)  Eisenstein 84176158 4.184 s (20.1 M/s)  SeedNtHash 85334538 2.991 s (28.5 M/s)  Gaussian/nt 0.36x  Eisenstein/nt 0.71x  (sinks 7a6c4c7d53d8845b 0264b0f0e2737252 a24e9c470ef16a29)
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

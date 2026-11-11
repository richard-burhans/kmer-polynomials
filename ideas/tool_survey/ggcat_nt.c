// Exhaustive check of GGCAT's canonical ntHash (crates/hashes/src/nthash_base.rs:49-61, cn_nthash.rs:29-31,98)
// over all m-mers: count distinct canonical m-mers whose (min(fw,rc)<<1) values collide.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static const uint64_t MUL=0x397f178c6ae330f9ULL;
static inline uint64_t rol(uint64_t x,unsigned r){r&=63;return r?(x<<r)|(x>>(64-r)):x;}
static inline uint64_t h(unsigned c){return (uint64_t)((c<<1)+1)*MUL;}        // compressed code c: A0 C1 T2 G3
static inline uint64_t rcv(unsigned c){return (uint64_t)(((c<<1)^4)+1)*MUL;}
static int cmp(const void*a,const void*b){uint64_t x=*(uint64_t*)a,y=*(uint64_t*)b;return x<y?-1:x>y;}
int main(int argc,char**argv){
  int m=atoi(argv[1]); uint64_t n=1ULL<<(2*m); uint64_t *v=malloc(n*8); uint64_t cnt=0;
  for(uint64_t x=0;x<n;x++){
    unsigned b[64]; for(int i=0;i<m;i++) b[i]=(x>>(2*i))&3;
    // canonical representative only (skip if revcomp code smaller), comp in A0C1T2G3 = c^2
    uint64_t r=0; for(int i=0;i<m;i++) r|=(uint64_t)(b[m-1-i]^2)<<(2*i);
    if(r<x) continue;
    uint64_t fh=0,bw=0;
    for(int i=0;i<m;i++){ fh^=rol(h(b[i]),m-1-i); bw^=rol(rcv(b[i]),i); }
    uint64_t c=(fh<bw?fh:bw)<<1; v[cnt++]=c;
  }
  qsort(v,cnt,8,cmp); uint64_t coll=0; for(uint64_t i=1;i<cnt;i++) coll+=(v[i]==v[i-1]);
  printf("m=%d canonical m-mers=%llu colliding pairs=%llu (random 63-bit expectation %.2e)\n",m,(unsigned long long)cnt,(unsigned long long)coll,(double)cnt*cnt/2/9.22e18);
}

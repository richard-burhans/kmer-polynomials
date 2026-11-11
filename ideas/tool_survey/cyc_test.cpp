// Tests: (1) metagraph KmerBloomFilter hash = min(CyclicHash(fwd), CyclicHash(rc)), wordsize 64,
//            seeds (0,1) as in metagraph common/hashers/rolling_hasher.hpp:18-19.
//            Claim: for k>64, swapping bases at positions i and i+64 gives identical fwd AND rc hash.
//        (2) Bifrost RepHash (MAX_KMER_SIZE<=64): canonical = wyhash(sorted{h,ht}).
//            Check the KmerStream defect (h^ht) on Bifrost's own state, and that Bifrost's hash() avoids it.
#include <cstdio>
#include <cstring>
#include <string>
#include <random>
#include <algorithm>
#include "cyclichash.h"
#define MAX_KMER_SIZE 32
#include "RepHash.hpp"
using namespace std;
static string rc(const string&s){string r(s.rbegin(),s.rend());for(char&c:r)c=(c=='A'?'T':c=='C'?'G':c=='G'?'C':'A');return r;}
static uint64_t cyc(const string&s,int k){CyclicHash<uint64_t,unsigned char> h(k,0,1,64);for(int i=0;i<k;i++)h.eat((unsigned char)s[i]);return h.hashvalue;}
static uint64_t mg_canon(const string&s,int k){return min(cyc(s,k),cyc(rc(s),k));}
int main(){
  mt19937_64 g(42); const char*B="ACGT";
  // (1)
  for(int k: {65, 81, 127}){
    int coll=0,trials=1000;
    for(int t=0;t<trials;t++){
      string x(k,'A'); for(auto&c:x)c=B[g()&3];
      int i=g()%(k-64); if(x[i]==x[i+64]) x[i+64]=B[(strchr(B,x[i])-B+1)&3];
      string y=x; swap(y[i],y[i+64]);
      if(x!=y && mg_canon(x,k)==mg_canon(y,k)) coll++;
    }
    printf("metagraph CyclicHash canonical, k=%d: swap(i,i+64) collisions %d/%d\n",k,coll,trials);
  }
  { // control k=63: swap(i,i+? ) has no alias
    int k=63,coll=0; for(int t=0;t<1000;t++){string x(k,'A');for(auto&c:x)c=B[g()&3];int i=g()%(k-1); if(x[i]==x[i+1]) x[i+1]=B[(strchr(B,x[i])-B+1)&3]; string y=x; swap(y[i],y[i+1]); if(mg_canon(x,k)==mg_canon(y,k)) coll++;}
    printf("metagraph CyclicHash canonical, k=63 control adjacent swap: %d/1000\n",coll);
  }
  // (2)
  for(int k: {21,31,23}){
    int xorc=0,bif=0,trials=10000; RepHash hf(k);
    for(int t=0;t<trials;t++){
      string x(k,'A'); for(auto&c:x)c=B[g()&3];
      string y=x; int m=k/2; y[m]=(y[m]=='A'?'T':y[m]=='C'?'G':y[m]=='G'?'C':'A');
      hf.init(x.c_str()); uint64_t hx=hf.hash(); uint64_t xx=hf.h^hf.ht;
      hf.init(y.c_str()); uint64_t hy=hf.hash(); uint64_t xy=hf.h^hf.ht;
      xorc += (xx==xy); bif += (hx==hy);
    }
    printf("Bifrost RepHash k=%d middle-base-complement: XOR(h,ht) collisions %d/%d; Bifrost hash() collisions %d/%d\n",k,xorc,trials,bif,trials);
  }
}

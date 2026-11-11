#include "RepHash.hpp"
#include <cstdio>
#include <string>
#include <random>
static char comp(char c){return c=='A'?'T':c=='T'?'A':c=='C'?'G':'C';}
static std::string rc(const std::string&s){std::string r(s.rbegin(),s.rend());for(auto&c:r)c=comp(c);return r;}
int main(){
  std::mt19937 g(1); const char*B="ACGT"; int coll=0,trials=10000,k=31;
  RepHash h(k); h.seed(12345);
  for(int t=0;t<trials;t++){
    std::string s(k,'A'); for(auto&c:s)c=B[g()%4];
    std::string s2=s; s2[k/2]=comp(s[k/2]);
    if(s2==rc(s)) continue;
    h.init(s.c_str()); uint64_t a=h.hash(); h.init(s2.c_str()); uint64_t b=h.hash();
    coll += (a==b);
  }
  // also check rolling == init and rc symmetry
  std::string s(k,'A'); for(auto&c:s)c=B[g()%4];
  h.init(s.c_str()); uint64_t a=h.hash(); h.init(rc(s).c_str());
  printf("k=%d middle-complement collisions: %d/%d ; rc-symmetric: %d\n",k,coll,trials,a==h.hash());
}

M=(1<<64)-1; C=0xbea225f9eb34556d; Ci=pow(C,-1,1<<64)
def mix(x):
    x^=x>>32; x=(x*C)&M; x^=x>>29; x=(x*C)&M; x^=x>>32; x=(x*C)&M; return x^(x>>29)
def unxs(y,s):
    x=y
    for _ in range(64//s+1): x=y^(x>>s)
    return x
def inv(y):
    x=unxs(y,29); x=(x*Ci)&M; x=unxs(x,32); x=(x*Ci)&M; x=unxs(x,29); x=(x*Ci)&M; return unxs(x,32)
p=inv(M); assert mix(p)==M
print(hex(p), "bits", p.bit_length())
for k in (15,21,25,31):
    if p < 4**k: print("valid packed k-mer for k =",k)

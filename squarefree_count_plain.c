#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;
typedef int64_t i64;
static u64 isqrt128(u128 n){
    u64 r=(u64)sqrtl((long double)n);
    while((u128)r*r>n) r--;
    while((u128)(r+1)*(r+1)<=n) r++;
    return r;
}
static u128 parse(const char*s){
    u128 v=0;
    for(;*s;s++) v=v*10+(u128)(*s-'0');
    return v;
}
static void print128(i128 v){
    char b[64];
    int n=0;
    int neg=v<0;
    u128 u=neg?(u128)(-v):(u128)v;
    if(u==0) b[n++]='0';
    while(u){ b[n++]='0'+(int)(u%10); u/=10; }
    if(neg) putchar('-');
    while(n) putchar(b[--n]);
}
static u64 *primes_upto(u64 n,size_t*cnt){
    char*c=calloc(n+1,1);
    u64*p=malloc(sizeof(u64)*(n/2+16));
    size_t k=0;
    for(u64 i=2;i<=n;i++){
        if(!c[i]){ p[k++]=i; for(u64 j=i*i;j<=n;j+=i) c[j]=1; }
    }
    free(c);
    *cnt=k;
    return p;
}
i128 squarefree_count(u128 x,u64 D,int32_t**Mout){
    u64 r=isqrt128(x);
    if(D>r) D=r;
    if(D<1) D=1;
    int32_t*M=malloc(sizeof(int32_t)*(D+1));
    M[0]=0;
    size_t np;
    u64 *pr=primes_upto(isqrt128(D)+1,&np);
    const u64 SEG=1u<<16;
    int8_t*mu=malloc(SEG);
    uint32_t*prod=malloc(sizeof(uint32_t)*SEG);
    i128 S1=0;
    int32_t run=0;
    for(u64 lo=1;lo<=D;lo+=SEG){
        u64 hi=lo+SEG-1;
        if(hi>D) hi=D;
        u64 len=hi-lo+1;
        for(u64 i=0;i<len;i++){ mu[i]=1; prod[i]=1; }
        for(size_t q=0;q<np;q++){
            u64 p=pr[q];
            if(p*p>hi) break;
            u64 st=((lo+p-1)/p)*p;
            for(u64 j=st;j<=hi;j+=p){ mu[j-lo]=-mu[j-lo]; prod[j-lo]*=(uint32_t)p; }
            u64 p2=p*p;
            u64 st2=((lo+p2-1)/p2)*p2;
            for(u64 j=st2;j<=hi;j+=p2) mu[j-lo]=0;
        }
        for(u64 i=0;i<len;i++){
            u64 n=lo+i;
            int m=mu[i];
            if(m&&prod[i]!=(uint32_t)n) m=-m;
            run+=m;
            M[n]=run;
            if(m){
                u128 d2=(u128)n*n;
                u128 q;
                if(x<((u128)1<<64)) q=(u64)x/(u64)d2;
                else if((u64)(x>>64)<(u64)d2){ u64 qq,rr; u64 lo64=(u64)x,hi64=(u64)(x>>64); __asm__("divq %4":"=a"(qq),"=d"(rr):"a"(lo64),"d"(hi64),"r"((u64)d2)); q=qq; }
                else q=x/d2;
                S1+=m>0?(i128)q:-(i128)q;
            }
        }
    }
    free(mu);free(prod);free(pr);
    u128 Dp=(u128)(D+1)*(D+1);
    u64 T=(u64)(x/Dp);
    i64*Mb=malloc(sizeof(i64)*(T+2));
    for(u64 t=T;t>=1;t--){
        u64 y=isqrt128(x/t);
        u64 s=(u64)sqrt((double)y);
        while(s*s>y) s--;
        while((s+1)*(s+1)<=y) s++;
        u64 kb=(u64)sqrt((double)(T/t));
        while(kb*kb>T/t) kb--;
        while((kb+1)*(kb+1)<=T/t) kb++;
        u64 sum=0;
        for(u64 k=2;k<=kb;k++) sum+=(u64)Mb[t*k*k];
        u64 sk=s>kb?s:kb;
        u64 vmax=y/(sk+1);
        if(y<((u64)1<<50)){
            double yd=(double)y;
            for(u64 k=kb+1;k<=s;k++){
                u64 q=(u64)(yd/(double)k);
                sum+=(u64)(i64)M[q];
            }
            u64 prev=y;
            for(u64 v=1;v<=vmax;v++){
                u64 nx=(u64)(yd/(double)(v+1));
                u64 lo2=nx<sk?sk:nx;
                sum+=(u64)(i64)M[v]*(prev-lo2);
                prev=nx;
            }
        } else {
            if(kb<s) for(u64 k=kb+1;k<=s;k++) sum+=(u64)(i64)M[y/k];
            u64 prev=y;
            for(u64 v=1;v<=vmax;v++){
                u64 nx=y/(v+1);
                u64 lo2=nx<sk?sk:nx;
                sum+=(u64)(i64)M[v]*(prev-lo2);
                prev=nx;
            }
        }
        Mb[t]=1-(i64)sum;
    }
    i128 S2=0;
    for(u64 t=1;t<=T;t++) S2+=(i128)Mb[t]-(i128)M[D];
    free(Mb);
    if(Mout) *Mout=M; else free(M);
    return S1+S2;
}
int main(int argc,char**argv){
    if(argc<2) return 1;
    u128 x=parse(argv[1]);
    double c=argc>2?atof(argv[2]):0.5;
    u64 D=(u64)(c*pow((double)x,0.4));
    if(D<2) D=2;
    clock_t t0=clock();
    i128 q=squarefree_count(x,D,NULL);
    double dt=(double)(clock()-t0)/CLOCKS_PER_SEC;
    print128(q);
    printf(" %.3f\n",dt);
    return 0;
}

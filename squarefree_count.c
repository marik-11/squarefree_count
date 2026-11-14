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
typedef uint32_t u32;
#ifdef W11
enum { WH = 2310, PHI = 480, NE = 32, NPR = 5, SEGA = 1024 };
static const int wpr[NPR] = {2, 3, 5, 7, 11};
#else
enum { WH = 210, PHI = 48, NE = 16, NPR = 4, SEGA = 2048 };
static const int wpr[NPR] = {2, 3, 5, 7};
#endif
static int res[PHI], Rr[WH], ne, eg[NE], emu[NE];
static u128 X;
static double XD;
static u64 DMIN;
static inline u64 cnt(u64 m) { return (m / WH) * PHI + Rr[m % WH]; }
static u64 isqrt64(u64 n) {
    u64 r = (u64)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}
static u64 isqrt128(u128 n) {
    u64 r = (u64)sqrtl((long double)n);
    while ((u128)r * r > n) r--;
    while ((u128)(r + 1) * (r + 1) <= n) r++;
    return r;
}
static u128 parse(const char *s) {
    u128 v = 0;
    for (; *s; s++) v = v * 10 + (u128)(*s - '0');
    return v;
}
static void print128(i128 v) {
    char b[64];
    int n = 0;
    int neg = v < 0;
    u128 u = neg ? (u128)(-v) : (u128)v;
    if (u == 0) b[n++] = '0';
    while (u) { b[n++] = '0' + (int)(u % 10); u /= 10; }
    if (neg) putchar('-');
    while (n) putchar(b[--n]);
}
static void wheel_init(void) {
    int c = 0;
    for (int r = 1; r < WH; r++) {
        int ok = 1;
        for (int b = 0; b < NPR; b++) if (r % wpr[b] == 0) ok = 0;
        if (ok) res[c++] = r;
    }
    for (int r = 0; r < WH; r++) {
        int k = 0;
        for (int j = 0; j < PHI; j++) if (res[j] <= r) k++;
        Rr[r] = k;
    }
    int tmpe[NE], tmpm[NE];
    for (int m = 0; m < NE; m++) {
        int e = 1, sg = 1;
        for (int b = 0; b < NPR; b++) if (m >> b & 1) { e *= wpr[b]; sg = -sg; }
        tmpe[m] = e; tmpm[m] = sg;
    }
    ne = NE;
    for (int a = 0; a < NE; a++) {
        int best = -1;
        for (int b = 0; b < NE; b++) if (tmpe[b] > 0 && (best < 0 || tmpe[b] < tmpe[best])) best = b;
        eg[a] = tmpe[best]; emu[a] = tmpm[best]; tmpe[best] = -1;
    }
}
static inline u128 qfun(u64 d) {
    u64 d2 = d * d;
    if (d >= DMIN) {
        u64 est = (u64)(XD / (double)d2);
        u128 p = (u128)est * d2;
        if (p > X) est--;
        else if (X - p >= d2) est++;
        return est;
    }
    if ((X >> 64) == 0) return (u64)X / d2;
    if ((u64)(X >> 64) < d2) {
        u64 qq, rr, lo = (u64)X, hi = (u64)(X >> 64);
        __asm__("divq %4" : "=a"(qq), "=d"(rr) : "a"(lo), "d"(hi), "r"(d2));
        return qq;
    }
    return X / d2;
}
static i64 powmod(i64 b, i64 e, i64 m) {
    i64 r = 1; b %= m;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    return r;
}
static i64 invmod(i64 a, i64 m) {
    i64 g = m, x0 = 0, x1 = 1, aa = a % m;
    i64 r0 = m, r1 = aa;
    x0 = 0; x1 = 1;
    while (r1) {
        i64 q = r0 / r1;
        i64 t = r0 - q * r1; r0 = r1; r1 = t;
        t = x0 - q * x1; x0 = x1; x1 = t;
    }
    (void)g;
    x0 %= m; if (x0 < 0) x0 += m;
    return x0;
}
i128 squarefree_count(u128 x, u64 D, double *tsieve) {
    wheel_init();
    X = x; XD = (double)x;
    u64 r = isqrt128(x);
    u64 q4 = (u64)sqrt((double)r) + 3;
    if (D < q4) D = q4;
    if (D > r) D = r;
    if (D < 2) D = 2;
    DMIN = isqrt128(x >> 40) + 1;
    u64 ncop = cnt(D);
    int32_t *Mc = malloc(sizeof(int32_t) * (ncop + 2));
    Mc[0] = 0;
    u64 sq = isqrt64(D) + 1;
    char *cs = calloc(sq + 1, 1);
    u64 *pl = malloc(sizeof(u64) * (sq + 8));
    size_t np = 0;
    for (u64 i = 2; i <= sq; i++) {
        if (!cs[i]) { if (i > (u64)wpr[NPR - 1]) pl[np++] = i; for (u64 j = i * i; j <= sq; j += i) cs[j] = 1; }
    }
    free(cs);
    u64 *nx1 = malloc(sizeof(u64) * np * PHI + 8), *nx2 = malloc(sizeof(u64) * np * PHI + 8);
    for (size_t k = 0; k < np; k++) {
        i64 p = (i64)pl[k], p2 = p * p;
        i64 wi = powmod(WH % p, p - 2, p);
        i64 wi2 = invmod(WH, p2);
        for (int i = 0; i < PHI; i++) {
            nx1[k * PHI + i] = (u64)((p - (res[i] % p) * wi % p) % p);
            nx2[k * PHI + i] = (u64)((p2 - (i64)res[i] * wi2 % p2) % p2);
        }
    }
    enum { S = SEGA };
    int8_t *mu = malloc((size_t)PHI * S);
    u32 *prod = malloc(sizeof(u32) * (size_t)PHI * S);
    u64 NV = (u64)sqrt(sqrt((double)x)) + 4;
    u64 *nl = malloc(sizeof(u64) * (NV + 8));
    int *ml = malloc(sizeof(int) * (NV + 8));
    size_t nn = 0;
    i128 S1 = 0;
    int32_t run = 0;
    u64 amax = D / WH;
    for (u64 A0 = 0; A0 <= amax; A0 += S) {
        u64 A1 = A0 + S;
        memset(mu, 1, (size_t)PHI * S);
        for (size_t i = 0; i < (size_t)PHI * S; i++) prod[i] = 1;
        for (size_t k = 0; k < np; k++) {
            u64 p = pl[k], p2 = p * p;
            for (int i = 0; i < PHI; i++) {
                u64 a = nx1[k * PHI + i];
                int8_t *m = mu + (size_t)i * S;
                u32 *pr = prod + (size_t)i * S;
                for (; a < A1; a += p) { m[a - A0] = -m[a - A0]; pr[a - A0] *= (u32)p; }
                nx1[k * PHI + i] = a;
                a = nx2[k * PHI + i];
                for (; a < A1; a += p2) m[a - A0] = 0;
                nx2[k * PHI + i] = a;
            }
        }
        for (u64 aa = 0; aa < S; aa++) {
            u64 a = A0 + aa;
            for (int i = 0; i < PHI; i++) {
                u64 n = a * WH + res[i];
                if (n > D) goto done;
                int m = mu[(size_t)i * S + aa];
                if (m && prod[(size_t)i * S + aa] != (u32)n) m = -m;
                run += m;
                Mc[a * PHI + i + 1] = run;
                if (m) {
                    if (n <= NV) { nl[nn] = n; ml[nn] = m; nn++; }
                    for (int e = 0; e < ne; e++) {
                        u64 d = n * (u64)eg[e];
                        if (d > D) break;
                        i128 qv = (i128)qfun(d);
                        S1 += (m * emu[e] > 0) ? qv : -qv;
                    }
                }
            }
        }
    }
done:
    free(mu); free(prod); free(nx1); free(nx2); free(pl);
    if (tsieve) *tsieve = (double)clock() / CLOCKS_PER_SEC;
    u128 Dp = (u128)(D + 1) * (D + 1);
    u64 T = (u64)(x / Dp);
    i64 *Mb = malloc(sizeof(i64) * (T + 2));
    u64 *Yb = malloc(sizeof(u64) * (T + 2));
    for (u64 t = T; t >= 1; t--) {
        u64 y = isqrt128(x / t);
        Yb[t] = y;
        u64 s = isqrt64(y), kb = isqrt64(T / t);
        u64 sum = 0;
        u64 jA = cnt(kb), cs_ = cnt(s);
        {
            u64 j = 2, a = 0; int i = 1;
            for (; j <= jA; j++) {
                u64 k = a * WH + res[i];
                sum += (u64)Mb[t * k * k];
                if (++i == PHI) { i = 0; a++; }
            }
        }
        u64 V = y / (s + 1);
        if (y < ((u64)1 << 50)) {
            double yd = (double)y;
            u64 j = jA + 1, a = jA / PHI; int i = (int)(jA % PHI);
            u64 ib[128];
            while (j <= cs_) {
                u64 m = cs_ - j + 1;
                if (m > 128) m = 128;
                for (u64 w = 0; w < m; w++) {
                    u64 k = a * WH + res[i];
                    u64 q = (u64)(yd / (double)k);
                    u64 id = cnt(q);
                    ib[w] = id;
                    __builtin_prefetch(Mc + id);
                    if (++i == PHI) { i = 0; a++; }
                }
                for (u64 w = 0; w < m; w++) sum += (u64)(i64)Mc[ib[w]];
                j += m;
            }
            for (size_t z = 0; z < nn && nl[z] <= V; z++) {
                u64 q = (u64)(yd / (double)nl[z]);
                sum += (u64)((i64)ml[z] * (i64)(cnt(q) - cs_));
            }
        } else {
            u64 j = jA + 1, a = jA / PHI; int i = (int)(jA % PHI);
            for (; j <= cs_; j++) {
                u64 k = a * WH + res[i];
                sum += (u64)(i64)Mc[cnt(y / k)];
                if (++i == PHI) { i = 0; a++; }
            }
            for (size_t z = 0; z < nn && nl[z] <= V; z++)
                sum += (u64)((i64)ml[z] * (i64)(cnt(y / nl[z]) - cs_));
        }
        Mb[t] = 1 - (i64)sum;
    }
    i128 MD = 0;
    for (int e = 0; e < ne; e++) MD += (i128)emu[e] * Mc[cnt(D / eg[e])];
    i128 S2 = 0;
    for (u64 t = 1; t <= T; t++) {
        i128 Mv = 0;
        for (int e = 0; e < ne; e++) {
            u64 te = t * (u64)eg[e] * (u64)eg[e];
            i64 v = te <= T ? Mb[te] : (i64)Mc[cnt(Yb[t] / eg[e])];
            Mv += (i128)emu[e] * v;
        }
        S2 += Mv - MD;
    }
    free(Mb); free(Yb); free(Mc); free(nl); free(ml);
    return S1 + S2;
}
int main(int argc, char **argv) {
    if (argc < 2) return 1;
    u128 x = parse(argv[1]);
    double c = argc > 2 ? atof(argv[2]) : 0.5;
    u64 D = (u64)(c * pow((double)x, 0.4));
    double ts = 0;
    clock_t t0 = clock();
    i128 q = squarefree_count(x, D, &ts);
    double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
    print128(q);
    printf(" %.3f %.3f\n", dt, ts);
    return 0;
}

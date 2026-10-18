/* Exhaustive check over all labeled digraphs on n vertices (n<=6):
   for strongly connected, non-symmetric digraphs, test the back-arc condition (Thm 10.7 of ALC)
   and the clearing criterion (with normality propagation). Report counts and any digraph that
   passes the back-arc test but is not fully clearable. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n, full;
static unsigned out[8], inn[8];
static unsigned char can[1 << 8];
static int backv[1 << 8];
static int path[8], baok;
static void dfs(int len, int mask, int back) {
    if (!baok) return;
    if (backv[mask] < 0) backv[mask] = back; else if (backv[mask] != back) { baok = 0; return; }
    int v = path[len - 1];
    unsigned o = out[v] & ~mask;
    while (o) {
        int w = __builtin_ctz(o); o &= o - 1;
        int b = __builtin_popcount(out[w] & mask);   /* arcs from w back to earlier path vertices */
        path[len] = w;
        dfs(len + 1, mask | (1 << w), back + b);
        if (!baok) return;
    }
}
static int strong(void) {
    unsigned seen = 1, frontier = 1;
    while (frontier) { int v = __builtin_ctz(frontier); frontier &= frontier - 1; unsigned nw = out[v] & ~seen; seen |= nw; frontier |= nw; }
    if (seen != (unsigned)full) return 0;
    seen = 1; frontier = 1;
    while (frontier) { int v = __builtin_ctz(frontier); frontier &= frontier - 1; unsigned nw = inn[v] & ~seen; seen |= nw; frontier |= nw; }
    return seen == (unsigned)full;
}
static int clearable(int refine) {
    memset(can, 0, sizeof(can));
    for (int v = 0; v < n; v++) can[1 << v] = 1 << v;
    for (int S = 1; S <= full; S++) {
        unsigned e = can[S];
        while (e) { int v = __builtin_ctz(e); e &= e - 1; unsigned o = out[v] & ~S; while (o) { int w = __builtin_ctz(o); o &= o - 1; can[S | (1 << w)] |= 1 << w; } }
    }
    unsigned normal = 0; int dk = can[full] != 0, changed = 1;
    while (changed) {
        changed = 0;
        if (refine && !dk) for (int S = 1; S <= full; S++) if (can[S] && ((full & ~S) & ~normal) == 0) { dk = 1; break; }
        for (int x = 0; x < n; x++) {
            if (normal >> x & 1) continue;
            for (int S = 1; S <= full; S++) {
                if (!(can[S] >> x & 1)) continue;
                unsigned rest = full & ~S, r = rest & ~out[x]; int ok = 1;
                while (r) {
                    int w = __builtin_ctz(r); r &= r - 1;
                    if (can[S | (1 << w)]) continue;
                    if (refine && dk && (((rest & ~(1u << w)) & ~normal) == 0)) continue;
                    ok = 0; break;
                }
                if (ok) { normal |= 1u << x; changed = 1; break; }
            }
        }
    }
    return normal == (unsigned)full;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); full = (1 << n) - 1;
    int np = n * (n - 1); int pi[64], pj[64], k = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { pi[k] = i; pj[k] = j; k++; }
    long long nstrong = 0, nonsym = 0, ba = 0, ba_bad = 0, notclear = 0, notclear_basic = 0;
    for (long long mask = 0; mask < (1LL << np); mask++) {
        if ((mask & 0x3FFFFFF) == 0) { fprintf(stderr, "progress %lld / %lld\n", mask >> 26, (1LL << np) >> 26); fflush(stderr); }
        memset(out, 0, sizeof(out)); memset(inn, 0, sizeof(inn));
        for (int e = 0; e < np; e++) if (mask >> e & 1) { out[pi[e]] |= 1u << pj[e]; inn[pj[e]] |= 1u << pi[e]; }
        if (!strong()) continue;
        nstrong++;
        int sym = 1; for (int v = 0; v < n; v++) if (out[v] != inn[v]) { sym = 0; break; }
        if (sym) continue;
        nonsym++;
        if (clearable(1)) continue;
        notclear++;
        for (int m = 0; m <= full; m++) backv[m] = -1;
        baok = 1;
        for (int s = 0; s < n && baok; s++) { path[0] = s; dfs(1, 1 << s, 0); }
        if (!baok) continue;
        ba_bad++;
        if (ba_bad <= 20) { printf("BA-pass, not clearable:"); for (int e = 0; e < np; e++) if (mask >> e & 1) printf(" %d>%d", pi[e], pj[e]); printf("\n"); }
    }
    printf("n=%d labeled: strongly connected %lld, non-symmetric %lld, not clearable (refined) %lld, of which pass the back-arc test %lld\n", n, nstrong, nonsym, notclear, ba_bad);
    return 0;
}

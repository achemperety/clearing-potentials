/* Read digraphs in directg -T format ("n ne a b a b ...") and, for the strongly connected ones that are not
   symmetric, apply the clearing criterion (with normality propagation). For the non-clearable ones apply the
   back-arc test; print those that pass it (the only ones that can carry a set potential). Counts on stderr. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n, full;
static unsigned out[16], inn[16];
static unsigned short can[1 << 12];
static int backv[1 << 12];
static int path[16], baok;
static void dfs(int len, int mask, int back) {
    if (!baok) return;
    if (backv[mask] < 0) backv[mask] = back; else if (backv[mask] != back) { baok = 0; return; }
    int v = path[len - 1];
    unsigned o = out[v] & ~mask;
    while (o) {
        int w = __builtin_ctz(o); o &= o - 1;
        int b = __builtin_popcount(out[w] & mask);
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
static unsigned reachm[16];
static void reach_all(void) {
    for (int v = 0; v < n; v++) {
        unsigned seen = 1u << v, frontier = 1u << v;
        while (frontier) { int x = __builtin_ctz(frontier); frontier &= frontier - 1; unsigned nw = out[x] & ~seen; seen |= nw; frontier |= nw; }
        reachm[v] = seen;
    }
}
static int unilateral(void) {
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (!((reachm[u] >> v & 1) || (reachm[v] >> u & 1))) return 0;
    return 1;
}
static int scc_symmetric(void) {
    for (int u = 0; u < n; u++) { unsigned o = out[u]; while (o) { int v = __builtin_ctz(o); o &= o - 1; if ((reachm[v] >> u & 1) && !(out[v] >> u & 1)) return 0; } }
    return 1;
}
static int clearable(void) {
    memset(can, 0, sizeof(unsigned short) * (full + 1));
    for (int v = 0; v < n; v++) can[1 << v] = 1 << v;
    for (int S = 1; S <= full; S++) {
        unsigned e = can[S];
        while (e) { int v = __builtin_ctz(e); e &= e - 1; unsigned o = out[v] & ~S; while (o) { int w = __builtin_ctz(o); o &= o - 1; can[S | (1 << w)] |= 1 << w; } }
    }
    unsigned normal = 0; int dk = can[full] != 0, changed = 1;
    while (changed) {
        changed = 0;
        if (!dk) for (int S = 1; S <= full; S++) if (can[S] && ((full & ~S) & ~normal) == 0) { dk = 1; break; }
        for (int x = 0; x < n; x++) {
            if (normal >> x & 1) continue;
            for (int S = 1; S <= full; S++) {
                if (!(can[S] >> x & 1)) continue;
                unsigned rest = full & ~S, r = rest & ~out[x]; int ok = 1;
                while (r) {
                    int w = __builtin_ctz(r); r &= r - 1;
                    if (can[S | (1 << w)]) continue;
                    if (dk && (((rest & ~(1u << w)) & ~normal) == 0)) continue;
                    ok = 0; break;
                }
                if (ok) { normal |= 1u << x; changed = 1; break; }
            }
        }
    }
    return normal == (unsigned)full;
}
int main(void) {
    int ne; static int ed[2 * 256];
    long long total = 0, nstrong = 0, nonsym = 0, notclear = 0, ba_bad = 0;
    while (scanf("%d %d", &n, &ne) == 2) {
        full = (1 << n) - 1; total++;
        memset(out, 0, sizeof(out)); memset(inn, 0, sizeof(inn));
        for (int i = 0; i < ne; i++) { int a, b; if (scanf("%d %d", &a, &b) != 2) return 1; ed[2*i] = a; ed[2*i+1] = b; out[a] |= 1u << b; inn[b] |= 1u << a; }
        reach_all();
        if (!unilateral()) continue;
        nstrong++;
        if (scc_symmetric()) continue;
        nonsym++;
        if (clearable()) continue;
        notclear++;
        for (int m = 0; m <= full; m++) backv[m] = -1;
        baok = 1;
        for (int s = 0; s < n && baok; s++) { path[0] = s; dfs(1, 1 << s, 0); }
        if (!baok) continue;
        ba_bad++;
        printf("%d %d", n, ne); for (int i = 0; i < ne; i++) printf(" %d %d", ed[2*i], ed[2*i+1]); printf("\n");
        fflush(stdout);
    }
    fprintf(stderr, "digraphs %lld, unilateral %lld, not SCC-symmetric %lld, not clearable %lld, of which pass the back-arc test %lld\n",
            total, nstrong, nonsym, notclear, ba_bad);
    return 0;
}

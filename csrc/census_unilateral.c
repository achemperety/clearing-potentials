/*
 * census_unilateral.c -- unilateral variant of census.c for
 *   "Clearing positions and the rigidity of complex matching potentials on digraphs".
 *
 * Reads digraphs in nauty's simple text format (one per line: "n m a1 b1 ... am bm"),
 * as produced by   geng -c -q n | directg -T -q
 * and, for every unilateral digraph that is not SCC-symmetric, decides
 *   (1) whether it is fully clearable (Definition 3.2 of the paper), and, if not,
 *   (2) whether it passes the back-arc test (Lemma 2.3: any two directed paths with
 *       the same vertex set have the same number of back arcs).
 * Digraphs that are not fully clearable and pass the back-arc test are echoed to stdout;
 * counts go to stderr.  n <= 8.
 *
 * The program is written from the definitions in the paper and shares no code with the
 * earlier programs in csrc/earlier; agreement of the counts is the point.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 8

static int n;
static unsigned full;
static unsigned outn[MAXN], inn[MAXN];
static unsigned char ends_[1u << MAXN];   /* ends_[S]: vertices v such that a directed path with vertex set S ends at v */

static unsigned closure(unsigned start, const unsigned *adj) {
    unsigned seen = start, todo = start;
    while (todo) {
        int v = __builtin_ctz(todo);
        todo &= todo - 1;
        unsigned nw = adj[v] & ~seen;
        seen |= nw;
        todo |= nw;
    }
    return seen;
}

static void path_table(void) {
    memset(ends_, 0, (size_t)(full + 1));
    for (int v = 0; v < n; v++) ends_[1u << v] = (unsigned char)(1u << v);
    /* subsets in increasing numerical order: every proper subset of S precedes S */
    for (unsigned S = 1; S <= full; S++) {
        unsigned e = ends_[S];
        while (e) {
            int v = __builtin_ctz(e);
            e &= e - 1;
            unsigned o = outn[v] & ~S;
            while (o) {
                int w = __builtin_ctz(o);
                o &= o - 1;
                ends_[S | (1u << w)] |= (unsigned char)(1u << w);
            }
        }
    }
}

/* Normality propagation N_0 = {} subset N_1 subset ... ; returns N_infinity. */
static unsigned propagate(void) {
    unsigned N = 0;
    for (;;) {
        /* global condition of clause (b): some path set X with V \ X contained in N */
        int haveX = 0;
        for (unsigned X = 1; X <= full && !haveX; X++)
            if (ends_[X] && ((full & ~X) & ~N) == 0) haveX = 1;
        unsigned newN = N;
        for (int x = 0; x < n; x++) {
            if (N >> x & 1) continue;
            int found = 0;
            for (unsigned S = 1; S <= full && !found; S++) {
                if (!(ends_[S] >> x & 1)) continue;          /* a path with vertex set S ending at x */
                unsigned rest = full & ~S;
                unsigned need = rest & ~outn[x];              /* w outside S with x -/-> w */
                int ok = 1;
                while (need && ok) {
                    int w = __builtin_ctz(need);
                    need &= need - 1;
                    if (ends_[S | (1u << w)]) continue;                       /* clause (a) */
                    if (haveX && (((rest & ~(1u << w)) & ~N) == 0)) continue; /* clause (b) */
                    ok = 0;
                }
                if (ok) found = 1;
            }
            if (found) newN |= 1u << x;   /* uses N_k only: a Jacobi-style round */
        }
        if (newN == N) return N;
        N = newN;
    }
}

/* back-arc test by depth-first search over all directed paths */
static int bk[1u << MAXN];
static int ok_back;
static int pth[MAXN];
static void dfs(int len, unsigned mask, int back) {
    if (!ok_back) return;
    if (bk[mask] == -1) bk[mask] = back;
    else if (bk[mask] != back) { ok_back = 0; return; }
    int v = pth[len - 1];
    unsigned o = outn[v] & ~mask;
    while (o && ok_back) {
        int w = __builtin_ctz(o);
        o &= o - 1;
        int add = __builtin_popcount(outn[w] & mask);   /* arcs from the new vertex back into the path */
        pth[len] = w;
        dfs(len + 1, mask | (1u << w), back + add);
    }
}
static int back_arc_test(void) {
    for (unsigned S = 0; S <= full; S++) bk[S] = -1;
    ok_back = 1;
    for (int s = 0; s < n && ok_back; s++) { pth[0] = s; dfs(1, 1u << s, 0); }
    return ok_back;
}

int main(void) {
    static char line[4096];
    long long total = 0, strong = 0, nonsym = 0, hampath = 0, notclear = 0, survivors = 0;
    while (fgets(line, sizeof line, stdin)) {
        char *p = line, *q;
        long nv = strtol(p, &q, 10);
        if (q == p) continue;
        p = q;
        long m = strtol(p, &q, 10);
        p = q;
        if (nv < 1 || nv > MAXN) { fprintf(stderr, "bad n\n"); return 1; }
        n = (int)nv;
        full = (1u << n) - 1;
        memset(outn, 0, sizeof outn);
        memset(inn, 0, sizeof inn);
        for (long i = 0; i < m; i++) {
            long a = strtol(p, &q, 10); p = q;
            long b = strtol(p, &q, 10); p = q;
            outn[a] |= 1u << b;
            inn[b] |= 1u << a;
        }
        total++;
        /* unilateral: any two vertices are comparable */
        unsigned reach[MAXN];
        for (int v = 0; v < n; v++) reach[v] = closure(1u << v, outn);
        int uni = 1;
        for (int u = 0; u < n && uni; u++)
            for (int v = 0; v < n; v++)
                if (!((reach[u] >> v) & 1) && !((reach[v] >> u) & 1)) { uni = 0; break; }
        if (!uni) continue;
        strong++;                                   /* counts unilateral digraphs in this variant */
        /* SCC-symmetric: every arc a->b with b reaching a has its reverse */
        int sccsym = 1;
        for (int a = 0; a < n && sccsym; a++) {
            unsigned o = outn[a];
            while (o) {
                int b = __builtin_ctz(o);
                o &= o - 1;
                if (((reach[b] >> a) & 1) && !((outn[b] >> a) & 1)) { sccsym = 0; break; }
            }
        }
        if (sccsym) continue;
        nonsym++;                                   /* counts non-SCC-symmetric unilateral digraphs */
        path_table();
        if (ends_[full]) { hampath++; continue; }      /* Theorem 4.4: Hamiltonian path => fully clearable */
        if (propagate() == full) continue;
        notclear++;
        if (!back_arc_test()) continue;
        survivors++;
        fputs(line, stdout);
    }
    fprintf(stderr,
            "digraphs %lld | unilateral %lld | not SCC-symmetric %lld | with Hamiltonian path %lld | "
            "not fully clearable %lld | of these passing the back-arc test %lld\n",
            total, strong, nonsym, hampath, notclear, survivors);
    return 0;
}

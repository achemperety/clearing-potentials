# clearing-potentials

Companion code for

> R. Baiandin, *Clearing positions and the rigidity of complex matching potentials on digraphs*
> (Graduate School of Economics and Management, Ural Federal University).

The paper continues R. Baiandin, *Alternation, invariants, and the limits of collapse* (ALC), and uses its notation:

- `R(W,v)`: root resolvents of directed vertex geography;
- `(U,a)`: local weighted matching potentials;
- `M(S)`: weighted matching polynomials;
- `Gamma(P)`: products of resolvents along plays;
- `N_inf`: normality propagation.

Everything is exact.

## Requirements

- Python 3.9 or later. The package `potentials` uses only the standard library.
- A C compiler and nauty (`geng`, `directg`), for the exhaustive censuses.
- Optional: `networkx`, for Remark 5.4 in `verify_all.py`.
- Optional: `passagemath-singular`, for the Groebner computations of Section 6 in `groebner/`.
- `pytest` for the tests.

## Reproducing the paper

| Command | Time | Statements |
|---|---|---|
| `python3 -m pytest -q tests` | 1 s | unit tests |
| `python3 scripts/verify_all.py [--quick]` | 1 min | Lemmas 2.1-2.5, Thm 2.7, Lemma 3.1, Lemma 3.5, Prop. 3.8, Example 3.7, Thms 4.3, 4.5, Cor. 4.2, Remark 4.4, Example 4.6, Table 1, the census data of Thm 5.1 and Prop. 5.2 re-read and certified, Remark 5.4, Table 2 |
| `cc -O2 -o census csrc/census.c`, then `geng -c -q n \| directg -T -q \| ./census` for n = 3..6 | seconds | Thm 5.1(1); counts on stderr, see `logs/census_n*.err` |
| the same for n = 7, split as `directg -T -q -s0/2` and `-s1/2` | 1 hour on 2 cores | Thm 5.1(2): 705,421,509 classes, 1,041 not fully clearable, 72 back-arc survivors (`data/n7_survivors.txt`) |
| `python3 scripts/certify_survivors.py data/n7_survivors.txt --out data/n7_gamma_certificates.json` | 1 min | the 72 Gamma-conflicts at z = 7/3 |
| `csrc/census_unilateral.c` (same pipeline), then `certify_survivors.py data/unilateral6_survivors.txt --unilateral` | seconds | Prop. 5.2: 145 survivors = 132 Gamma-conflicts + 13 set potentials |
| `csrc/census_backarc.c` for n = 6 | seconds | the 79 six-vertex back-arc survivors (`data/strong6_backarc_survivors.txt`) |
| `cd groebner; python3 scf_groebner.py branches-pendants 4 1,1,1,1` (and the other rows of Table 3) | 0 s - 9 min | Thm 6.9 and Table 3; outputs in `logs/reruns/groebner_*.txt` |
| `cd groebner; python3 scf_groebner.py triangle-pendants 2,2,2 800` | 2.5 min | the triangle with two pendant vertices at every vertex, with the comparison on 800 graphs |
| `cd groebner; python3 scf_c4pendants.py` | minutes | Prop. 6.5, Cor. 6.6 (modular cross-check) |
| `cd groebner; python3 unilateral_open.py OUT 3 4` | 20 min | Prop. 5.2: digraphs #03 and #04 have no potential over any field (`logs/reruns/unilateral6_groebner_03_04.txt`) |
| `cd groebner; python3 unilateral_open.py OUT 5 12` | hours | attempts the eight open cases of Prop. 5.2 |

## Two independent censuses

`csrc/census.c` and its two variants were written from the definitions in the paper. The programs in `csrc/earlier/` were written independently, earlier:

- `clearall.c` runs over all labelled digraphs with n <= 6;
- `clearfilter.c` and `clearfilter_uni.c` read nauty output.

The two sets of programs share no code, and their counts agree:

- `logs/earlier/c_directg_*.txt`;
- `logs/reruns/labelled_n3_6_clearall.txt`: 14, 1,568, 564,352 and 734,748,072 labelled non-symmetric strongly connected digraphs, all fully clearable;
- the 72 seven-vertex survivors of both programs are isomorphic.

## Layout

```
potentials/   qz.py (exact Q[z], Q(z)), dvg.py (digraphs, path sets, positions, exact resolvents, Gamma,
              back-arc and Gamma-conflict tests, set potentials), clearing.py (N_inf with certificates),
              potential.py (matching polynomials, potentials, directed collapse), families.py (blow-ups,
              tournaments, cycles with pendants and legs, random SCC-symmetric digraphs)
scripts/      verify_all.py, certify_survivors.py
csrc/         census.c, census_unilateral.c, census_backarc.c; earlier/ independent programs
groebner/     exact Groebner computations over Q(z) (passagemath-singular): scf_core.py, scf_groebner.py,
              scf_c4pendants.py, unilateral_open.py
data/         survivors, certificates, set-potential lists
logs/         census counts; earlier/ outputs of the earlier programs; reruns/ re-runs made for this version
tests/        pytest suite
```

## License

MIT (see `LICENSE`). `MANIFEST.md` lists SHA-256 checksums of all files.

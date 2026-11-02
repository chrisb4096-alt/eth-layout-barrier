"""Exact checks for the far-pair counterexample (Section 6 of the paper).

1. The base families (q = 3, 4, 5): 2q+1 pairwise anticommuting Pauli labels (w, v)
   with distinct X-parts w, and nonzero determinant of the h_w^2 coefficient system
   at the times of directions 0..q-1. Translation invariance of that determinant is
   checked at a second placement.
2. The recursive families for q = 6..10: family axioms and nonzero determinant.
3. The Clifford cluster for q = 3, 4, 5 with 2q+2 prescribed outside roots: M(0) = I,
   deg M <= q+1, every entry measured at an inside time vanishes, N^2 = kappa I,
   g equals the prescribed norm, and its roots are simple (so rank M = 2^(q-1) there).
All arithmetic is exact (rationals and radicals). Run from the repository root:
    uv run --with sympy python verify/far_pair.py
"""
import sympy as S

T = S.Symbol('T')
pc = lambda x: bin(x).count('1')

BASE = {
    3: [(1, 0), (2, 1), (3, 5), (4, 3), (5, 1), (6, 3), (7, 5)],
    4: [(0, 3), (1, 0), (2, 1), (5, 1), (6, 11), (9, 5), (10, 3), (13, 5), (14, 11)],
    5: [(5, 26), (6, 18), (7, 23), (10, 30), (14, 17), (16, 0), (17, 31), (19, 19),
        (22, 27), (24, 22), (28, 17)],
}


def check_family(fam, q):
    assert len(fam) == 2*q + 1 and len({w for w, _ in fam}) == 2*q + 1
    for i in range(len(fam)):
        assert fam[i] != (0, 0)
        for k in range(i):
            (w, v), (w2, v2) = fam[i], fam[k]
            assert (pc(w & v2) + pc(w2 & v)) % 2 == 1


def times(P):
    return [2*p + 2 for p in P], [2*p + 1 for p in P]  # alpha (antisymmetric), sigma (symmetric)


def h(w, P):
    al, si = times(P)
    return S.expand(S.prod(T - (si[i] if w >> i & 1 else al[i]) for i in range(len(P))))


def det_system(fam, P):
    q = len(P)
    rows = [S.Poly(h(w, P)**2, T) for w, _ in fam]
    return S.Matrix([[r.nth(e) for r in rows] for e in range(2*q + 1)]).det(method='bareiss')


def combine(f1, q1, f2):
    g = f1[0]
    return f1[1:] + [(g[0] | (w << q1), g[1] | (v << q1)) for w, v in f2]


def family(q):
    if q in BASE: return BASE[q]
    fours = q % 3
    parts = [3]*((q - 4*fours)//3) + [4]*fours
    fam, q1 = BASE[parts[0]], parts[0]
    for q2 in parts[1:]:
        fam = combine(fam, q1, BASE[q2]); q1 += q2
    return fam


def cluster(q, fam, roots):
    """Build M for directions P = 0..q-1 and verify it exactly. roots: outside directions."""
    P = list(range(q)); n = 2**q
    al, si = times(P)
    F = S.expand(S.prod(1 - T/S.Integer(a) for a in al))
    H = [S.expand(T*h(w, P)) for w, _ in fam]
    g = S.Poly(S.prod(1 - T/S.Integer(2*l + 2) for l in roots), T)
    b = (g.nth(1) - 2*S.diff(F, T).subs(T, 0))/2
    rhs = S.Poly(S.expand(F**2*(1 + b*T)**2 - g.as_expr()), T)
    V = S.Matrix([[S.Poly(x*x, T).nth(e) for x in H] for e in range(2, 2*q + 3)])
    k = list(V.LUsolve(S.Matrix([rhs.nth(e) for e in range(2, 2*q + 3)])))
    kappa = S.expand(sum(ki*x*x for ki, x in zip(k, H)))
    assert S.Poly(S.expand(F**2*(1 + b*T)**2 - kappa), T) == g
    eps = [(-1)**pc(w & v) for w, v in fam]
    c = [S.sqrt(ki*e) for ki, e in zip(k, eps)]
    d = S.expand(F*(1 + b*T))
    sector = {w: i for i, (w, _) in enumerate(fam)}

    def N(x, y):  # entry of N = sum c_w H_w G_{w,v}; G has (-1)^{v.y} at (y^w, y)
        i = sector.get(x ^ y)
        return 0 if i is None else c[i]*H[i]*(-1)**pc(fam[i][1] & y)

    def M(x, y):
        return N(x, y) + (d if x == y else 0)

    assert S.expand(d.subs(T, 0)) == 1 and all(x.subs(T, 0) == 0 for x in H)
    assert S.degree(d, T) <= q + 1 and all(S.degree(x, T) <= q + 1 for x in H)
    # measured inside entries: at alpha_p entries with x_p == y_p, at sigma_p with x_p != y_p
    for i in range(q):
        for x in range(n):
            for y in range(n):
                t = al[i] if (x >> i & 1) == (y >> i & 1) else si[i]
                e = M(x, y)
                assert e == 0 or S.expand(e.subs(T, t)) == 0, (i, x, y, t)
    # N^2 = kappa I
    for x in range(n):
        for z in range(n):
            s = sum(N(x, y)*N(y, z) for y in range(n))
            assert S.expand(s - (kappa if x == z else 0)) == 0
    # simple roots of g at the prescribed outside times
    for l in roots:
        t = 2*l + 2
        assert g.eval(t) == 0 and S.diff(g.as_expr(), T).subs(T, t) != 0
    return n


if __name__ == '__main__':
    for q, fam in BASE.items():
        check_family(fam, q)
        D0, D1 = det_system(fam, list(range(q))), det_system(fam, list(range(3, 3 + q)))
        assert D0 != 0 and D0 == D1
        print(f'q={q}: base family ok, determinant {D0} (same at a shifted placement)', flush=True)
    for q in range(6, 11):
        fam = family(q); check_family(fam, q)
        assert det_system(fam, list(range(q))) != 0
        print(f'q={q}: recursive family ok, determinant nonzero', flush=True)
    for q in (3, 4, 5):
        n = cluster(q, BASE[q], list(range(q, 3*q + 2)))
        print(f'q={q}: {n}x{n} cluster verified: rank {n//2} at {2*q+2} outside times, '
              f'zero measured inside entries', flush=True)
    print('DONE: all far-pair checks passed')

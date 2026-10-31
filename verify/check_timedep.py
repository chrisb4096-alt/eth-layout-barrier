#!/usr/bin/env python3
"""Exact rational layout/determinant certificates. No numerical rank tests.
Run: uv run --with sympy python verify/check_timedep.py
"""
from __future__ import annotations

import itertools
import json
from datetime import datetime
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
from zoneinfo import ZoneInfo

import sympy as s

T = s.Symbol('T')
HERE = Path(__file__).resolve().parent


def expansion(n, scopes):
    ratios = []
    for mask in range(1, 1 << n):
        w = {v for v in range(n) if mask >> v & 1}
        if n <= 4 * len(w) and 2 * len(w) <= n:
            crossed = sum(bool(w.intersection(e)) and not set(e) <= w for e in scopes)
            ratios.append(Fraction(crossed, len(w)))
    return str(min(ratios))


def layout(name, scopes, q, r, sigma, measured, roots=(1,)):
    n = max(max(e) for e in scopes) + 1
    parts = [[[v] if v < n else [] for v in range(r)] for _ in range(q)]
    c = [[s.eye(1) if g in measured[t] else s.zeros(0, len(parts[t][g]))
          for g in range(r)] for t in range(q)]
    p = {(f, v): s.prod(1 - T / s.Integer(t) for t in roots) * s.eye(n)[v, :]
         for f, e in enumerate(scopes) for v in e}
    return dict(name=name, n=n, q=q, r=r, R=len(roots), D=1, scopes=scopes,
                parts=parts, sigma=sigma, C=c, P=p,
                tau=Fraction(max(map(len, measured)), 2*n))


def verify(z):
    n, q, r, degree = (z[k] for k in ('n', 'q', 'r', 'R'))
    assert 1 <= degree < q and 0 < z['tau'] <= Fraction(1, 2)
    den = factorial(degree-1) * factorial(q)
    max_coeff = 0
    for t in range(q):
        assert sorted(v for group in z['parts'][t] for v in group) == list(range(n))
        assert all(len(group) * r <= 2*n for group in z['parts'][t])
        assert sum(c.rows for c in z['C'][t]) <= 2*z['tau']*n
        for g, group in enumerate(z['parts'][t]):
            c = z['C'][t][g]
            assert c.cols == len(group)
            assert all(x.is_Integer for x in c)
            assert all(sum(abs(x) for x in c.row(j)) <= z['D'] for j in range(c.rows))
    for (f, v), p in z['P'].items():
        assert p.subs(T, 0) == s.eye(n)[v, :]
        for x in p:
            poly = s.Poly(x, T)
            assert poly.degree() <= degree
            for a in poly.all_coeffs():
                assert (a * den).is_Integer
                max_coeff = max(max_coeff, abs(int(a * den)))
        for t in range(q):
            for g, group in enumerate(z['parts'][t]):
                if g != z['sigma'][f][t]:
                    row = p.subs(T, t+1).extract([0], group)
                    for null in z['C'][t][g].nullspace():
                        assert (row*null)[0] == 0
    return dict(label='checked (exact)', name=z['name'], n=n, Q=q, r=r, R=degree,
                D=z['D'], tau=str(z['tau']), balanced_cut_expansion=expansion(n, z['scopes']),
                maximum_incidence=max(sum(v in e for e in z['scopes']) for v in range(n)),
                denominator=den, maximum_scaled_coefficient=max_coeff,
                measurement_rows=[sum(c.rows for c in cs) for cs in z['C']])


def det_certificate(z, signs, L):
    q, r = z['q'], z['r']
    a = [[next(g for g in range(r) if v in z['parts'][t][g]) for t in range(q)]
         for v in range(z['n'])]
    score = lambda w: sum(signs[t][w[t]] for t in range(q))
    chosen = []
    for v in range(z['n']):
        if score(a[v]) <= -L:
            for f, edge in enumerate(z['scopes']):
                if v in edge and score(z['sigma'][f]) >= L:
                    chosen.append((f, v))
                    break
    m = len(chosen)
    A = s.Matrix(m, m, lambda i,j: z['P'][chosen[i]][chosen[j][1]])
    assert A.subs(T, 0) == s.eye(m)
    det = s.Poly(A.det(), T)
    deficits, deltas, ct = [], [], []
    for t in range(q):
        u = [sum(z['sigma'][f][t] == g for f, v in chosen) for g in range(r)]
        vv = [sum(a[v][t] == g for f, v in chosen) for g in range(r)]
        delta = sum(abs(x-y) for x,y in zip(u,vv)) // 2
        point = A.subs(T, t+1)
        M = point.copy()
        for i, (f,v) in enumerate(chosen):
            for j, (_,w) in enumerate(chosen):
                if z['sigma'][f][t] == a[w][t]:
                    M[i,j] = 0
        c = sum(mat.rows for mat in z['C'][t])
        assert M.rank() <= c
        assert (point-M).rank() <= sum(min(x,y) for x,y in zip(u,vv))
        deficit = m - point.rank()
        assert deficit >= delta-c
        _, rem = s.div(det, s.Poly((T-(t+1))**deficit, T))
        assert rem.is_zero
        deficits.append(deficit)
        deltas.append(delta)
        ct.append(c)
    assert sum(deltas) >= m*L
    assert sum(deficits) <= det.degree() <= m*z['R']
    assert m*(L-z['R']) <= sum(ct)
    return dict(label='checked (exact)', L=L, selected=chosen, m=m,
                determinant=str(det.as_expr()), point_rank_deficits=deficits,
                deltas=deltas, c_t=ct, lhs=m*(L-z['R']), rhs=sum(ct),
                equality=m*(L-z['R']) == sum(ct))


def infeasible_cycle():
    # Singleton groups, four times, no measurements, alternate endpoint assignments.
    # Solve every polynomial coefficient, not just diagonal ansatz.
    edges = [[0,1],[1,2],[2,3],[3,0]]
    answer = []
    for degree in (1,2):
        certs = []
        for f,e in enumerate(edges):
            sigma = [e[t%2] for t in range(4)]
            for v in e:
                rows, rhs = [], []
                for t in range(1,5):
                    for w in range(4):
                        if w == sigma[t-1]:
                            continue
                        rows.append([s.Integer(t)**k if y == w else 0
                                     for k in range(1,degree+1) for y in range(4)])
                        rhs.append(-int(w == v))
                A, b = s.Matrix(rows), s.Matrix(rhs)
                certs.append(dict(f=f,v=v,rank=A.rank(),augmented=A.row_join(b).rank()))
        ok = all(c['rank']==c['augmented'] for c in certs)
        assert ok == (degree == 2)
        answer.append(dict(label='checked (exact)', degree=degree, feasible=ok, ranks=certs))
    return answer


def swap_cycle():
    # Actual C4 incidence embedding in the Q=2 binary full swap array.
    q, r, radius = 2, 2, 1
    words = list(itertools.product(range(r), repeat=q))
    labels = list(itertools.product(words, words))
    index = {label:i for i,label in enumerate(labels)}
    dim = len(labels)
    I = s.eye(dim)
    swaps = []
    for t in range(q):
        H = s.zeros(dim)
        for i,(a,d) in enumerate(labels):
            aa,dd = list(a),list(d)
            aa[t],dd[t] = dd[t],aa[t]
            H[i,index[tuple(aa),tuple(dd)]] = 1
        swaps.append(H)
    W = sum(((I-H)/2 for H in swaps), s.zeros(dim))
    M = swaps[0]*swaps[1]*(I-2*W)  # interpolates (-1)^j on j=0,1
    B = I-M
    assert B.rank() == 1
    template = B.copy()
    for i,(a,d) in enumerate(labels):
        for j,(aa,dd) in enumerate(labels):
            if M[i,j]:
                roots = [t+1 for t in range(q) if d[t] != aa[t]]
                assert len(roots) <= radius
                template[i,j] += M[i,j]*s.prod(1-T/s.Integer(t) for t in roots)
    edges = [[0,2],[2,1],[1,3],[3,0]]
    aa = [(0,0),(0,0),(1,1),(1,1)]
    sigma = words
    embedding = s.zeros(dim,4)
    occupied = {}
    for f,e in enumerate(edges):
        for v in e:
            i = index[aa[v],sigma[f]]
            assert i not in occupied
            occupied[i] = v
            embedding[i,v] = 1
    p = {(f,v): template[index[aa[v],sigma[f]],:]*embedding
         for f,e in enumerate(edges) for v in e}
    basis = B.rowspace()[0]
    den = s.ilcm(*(x.q for x in basis))
    measurements = den*basis*embedding
    parts = [[[0,1],[2,3]] for _ in range(q)]
    c = [[measurements.extract([0], group) for group in part] for part in parts]
    z = dict(name='C4, embedded Q=2 swap construction', n=4,q=2,r=2,R=1,
             D=max(int(sum(abs(x) for x in mat)) for cs in c for mat in cs),
             scopes=edges,parts=parts,sigma=sigma,C=c,P=p,tau=Fraction(1,4))
    return z


def all_sign_scan(z,L):
    # Exhaustive sign scan at fixed L. Checks the determinant bound and score/variation step.
    n,q,r = z['n'],z['q'],z['r']
    used = sorted({g for t in range(q) for g in range(r)
                   if z['parts'][t][g] or any(w[t]==g for w in z['sigma'])})
    a = [[next(g for g in used if v in z['parts'][t][g]) for t in range(q)] for v in range(n)]
    rhs = sum(c.rows for cs in z['C'] for c in cs)
    max_m, grids, nonempty = 0,0,0
    best = None
    for flat in itertools.product((-1,1), repeat=q*len(used)):
        signs = [[1]*r for _ in range(q)]
        for t in range(q):
            for j,g in enumerate(used):
                signs[t][g] = flat[t*len(used)+j]
        scores = [sum(signs[t][w[t]] for t in range(q)) for w in z['sigma']]
        selected = []
        for v in range(n):
            if sum(signs[t][a[v][t]] for t in range(q)) <= -L:
                for f,e in enumerate(z['scopes']):
                    if v in e and scores[f] >= L:
                        selected.append((f,v))
                        break
        m = len(selected)
        assert m*(L-z['R']) <= rhs
        delta = sum(sum(abs(sum(z['sigma'][f][t]==g for f,v in selected)
                               -sum(a[v][t]==g for f,v in selected)) for g in used)//2
                    for t in range(q))
        assert delta >= m*L
        grids += 1
        nonempty += bool(m)
        if m > max_m:
            max_m,best = m,signs
    return dict(label='checked (exact)', L=L, exhaustive_sign_grids=grids,
                nonempty_grids=nonempty, largest_m=max_m,
                best_determinant=det_certificate(z,best,L) if best else None)


def probability(q,h,L):
    answer = Fraction(0)
    for k in range(h+1):
        for cp in range(q-k+1):
            C = 2*cp-(q-k)
            for zp in range(k+1):
                Z = 2*zp-k
                if C-Z <= -L and C+Z >= L:
                    answer += Fraction(comb(h,k)*comb(q-k,cp)*comb(k,zp),2**(h+q))
    return answer


def main():
    cycle = [[0,1],[1,2],[2,3],[3,0]]
    fixed = layout('C4, saturating determinant',cycle,4,8,
                   [[1]*4,[1]*4,[3]*4,[3]*4],[[],[0,2],[0,2],[0,2]])
    k4 = list(map(list,itertools.combinations(range(4),2)))
    measured = [[],[0,1,2],[1,2,3],[0,2,3]]
    sigma = [[next((v for v in e if v not in measured[t]),e[0]) for t in range(4)] for e in k4]
    moving = layout('K4, rotating measurement spaces',k4,4,8,sigma,measured)
    swap = swap_cycle()
    results = []
    for z in [fixed,moving,swap]:
        record = verify(z)
        record['P'] = {f'{f},{v}':[str(x) for x in p] for (f,v),p in z['P'].items()}
        record['partitions'] = z['parts']
        record['factor_schedules'] = z['sigma']
        record['measurements'] = [[c.tolist() for c in cs] for cs in z['C']]
        if z is fixed:
            record['equality_certificate'] = det_certificate(z,[[-1,1,-1,1,1,1,1,1]]*4,4)
        record['sign_scan'] = all_sign_scan(z, 2)
        results.append(record)
        print(record['label'], record['name'], 'D=',z['D'],'tau=',z['tau'], flush=True)
    probs = []
    for q in range(1,5):
        for h in range(q+1):
            for L in range(1,q+1):
                count = sum(sum(signs[2*t] for t in range(q)) <= -L
                            and sum(signs[2*t+int(t<h)] for t in range(q)) >= L
                            for signs in itertools.product((-1,1),repeat=2*q))
                p = probability(q,h,L)
                assert p == Fraction(count,2**(2*q))
                probs.append(dict(Q=q,distance=h,L=L,probability=str(p)))
    output = dict(generated_central=datetime.now(ZoneInfo('America/Chicago')).isoformat(),
                  label='checked (exact)', layouts=results,
                  fixed_schedule_feasibility=infeasible_cycle(),probability_checks=probs)
    # Serialize exact SymPy integers without turning them into floats.
    (HERE/'certificates.json').write_text(json.dumps(output,ensure_ascii=False,indent=2,
                                                    default=lambda x:int(x) if x.is_Integer else str(x))+'\n')
    print('checked (exact):',len(probs),'probability comparisons; degree 1 infeasible / degree 2 feasible')
    print('Wrote certificates.json')


if __name__ == '__main__':
    main()

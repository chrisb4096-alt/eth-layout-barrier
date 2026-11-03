# Exact check of the interpolation identity behind Proposition 10 (paired j-cube, j=4).
# A = I + diagonal part of degree R + off-diagonal part of degree d, random integer coefficients.
# For each direction l and d-1 helper directions H, P^-_l A(alpha_l) equals the sum over sign
# vectors eps of Lagrange interpolation of Q_eps P^-_l A(T)/T at sigma_l and the times t_h,
# where the data at t_h are the measured parts P^{eps_h}_h A(t_h).
import random, itertools
import sympy as sp
j=4; m=2**j
xs=list(itertools.product([0,1],repeat=j))
Z=[sp.diag(*[(-1)**x[l] for x in xs]) for l in range(j)]
def P(l,s,M): return (M+s*Z[l]*M*Z[l])/2
al=[sp.Integer(2*l+2) for l in range(j)]; sg=[sp.Integer(2*l+1) for l in range(j)]
random.seed(1); ok=True
for d,R in ((1,3),(2,3),(3,4)):
    Dc=[sp.diag(*[random.randint(-3,3) for _ in range(m)]) for _ in range(R)]
    Oc=[sp.Matrix(m,m,lambda a,b: 0 if a==b else random.randint(-3,3)) for _ in range(d)]
    def A(t): return sp.eye(m)+sum((t**(k+1)*Dc[k] for k in range(R)),sp.zeros(m))+sum((t**(k+1)*Oc[k] for k in range(d)),sp.zeros(m))
    for l in range(j):
        H=[h for h in range(j) if h!=l][:d-1]
        rhs=sp.zeros(m)
        for eps in itertools.product([1,-1],repeat=d-1):
            t=[al[h] if e==1 else sg[h] for h,e in zip(H,eps)]
            nodes=[sg[l]]+t
            data=[P(l,-1,A(sg[l]))]+[P(h,e,A(th)) for h,e,th in zip(H,eps,t)]
            for i,ti in enumerate(nodes):
                L=sp.prod([(al[l]-tk)/(ti-tk) for k,tk in enumerate(nodes) if k!=i])
                M=P(l,-1,data[i])
                for h,e in zip(H,eps): M=P(h,e,M)
                rhs+=al[l]/ti*L*M
        ok &= rhs==P(l,-1,A(al[l]))
print("ok" if ok else "FAIL")

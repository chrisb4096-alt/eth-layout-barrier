# Exact check of P^-_{y+y'}(X) = P^-_y(X) + Z_y P^-_{y'}(X) Z_y and the block-rank identity (k=3).
import random, itertools
from fractions import Fraction
import sympy as sp
k=3; m=2**k
xs=list(itertools.product([0,1],repeat=k))
dot=lambda a,b: sum(i*j for i,j in zip(a,b))%2
def Z(y): return sp.diag(*[(-1)**dot(x,y) for x in xs])
def Pm(M,y): return (M - Z(y)*M*Z(y))/2
def Pp(M,y): return (M + Z(y)*M*Z(y))/2
random.seed(1); ok=True
for _ in range(20):
    M=sp.Matrix(m,m,lambda i,j: random.randint(-3,3))
    y=random.choice(xs[1:]); y2=random.choice(xs[1:]); s=tuple((a+b)%2 for a,b in zip(y,y2))
    ok &= Pm(M,s)==Pm(M,y)+Z(y)*Pm(M,y2)*Z(y)
    # cost via blocks: antisym (e=1): rows c_t==g, cols c_t==g ; sym: rows c_t!=g
    for e in (0,1):
        cost=0
        for g in (0,1):
            R_=[i for i,x in enumerate(xs) if (dot(x,y)+e)%2!=g]; C_=[j for j,x in enumerate(xs) if dot(x,y)==g]
            cost+=M.extract(R_,C_).rank()
        ok &= cost==(Pp(M,y) if e else Pm(M,y)).rank()
print("ok" if ok else "FAIL")

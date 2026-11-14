import sys,time
from math import isqrt
from sympy import mobius
def A071172(n): return sum(mobius(k)*(10**n//k**2) for k in range(1, isqrt(10**n)+1))
for n in map(int,sys.argv[1:]):
    t=time.time(); v=A071172(n); print(n,v,round(time.time()-t,3))

"""Independent numerical identity checks; no physical measurements or device simulation.
Run with Python 3. Uses only the standard library. All scales are dimensionless:
c=epsilon0=mu0=A=k=1. Restoring SI units is described in the manuscript.
"""
import cmath
import math
import random

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def norm2(a): return sum(abs(v)**2 for v in a)
def inner(a,b): return sum(x.conjugate()*y for x,y in zip(a,b))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def scale(a,z): return tuple(z*x for x in a)
def close(a,b): assert abs(a-b)<2e-11,(a,b)

rng=random.Random(17012026)
for _ in range(500):
    th=rng.uniform(.001,math.pi/2-.001)
    s,h=math.sin(th),math.cos(th)
    x,z,phi=(rng.uniform(-6,6) for _ in range(3))
    q=s*x+phi/2
    ep=(0j,cmath.exp(1j*(h*z+q)),0j)
    em=(0j,cmath.exp(1j*(h*z-q)),0j)
    bp=cross((s,0,h),ep);bm=cross((-s,0,h),em)
    e,b=add(ep,em),add(bp,bm)
    close(norm2(e),4*math.cos(q)**2)
    close(norm2(b),4*(h*h*math.cos(q)**2+s*s*math.sin(q)**2))
    close((norm2(e)+norm2(b))/4,1+h*h*math.cos(2*q))
    flux=tuple(v.real/2 for v in cross(e,tuple(v.conjugate() for v in b)))
    close(flux[0],0);close(flux[1],0);close(flux[2],2*h*math.cos(q)**2)
    # Maxwell transversality and curl identities checked separately for each mode.
    for khat,ef,bf in [((s,0,h),ep,bp),((-s,0,h),em,bm)]:
        close(sum(kh*v for kh,v in zip(khat,ef)),0)
        close(sum(kh*v for kh,v in zip(khat,bf)),0)
        for lhs,rhs in zip(cross(khat,bf),scale(ef,-1)):close(lhs,rhs)
    # Quadratic mismatch bound for a randomly weighted positive diagonal detector.
    f=tuple(complex(rng.gauss(0,1),rng.gauss(0,1)) for _ in range(6))
    df=tuple(.01*complex(rng.gauss(0,1),rng.gauss(0,1)) for _ in range(6))
    weights=tuple(rng.uniform(0,3) for _ in range(6))
    response=lambda v:sum(w*abs(a)**2 for w,a in zip(weights,v))
    eps=math.sqrt(norm2(df))
    bound=max(weights)*(2*math.sqrt(norm2(f))*eps+eps**2)
    assert abs(response(add(f,df))-response(f))<=bound+1e-12
print('PASS: 500 randomized vector-field, Maxwell, energy/flux, and mismatch-bound checks.')
print('These checks verify algebra and conventions only, not experimental feasibility.')

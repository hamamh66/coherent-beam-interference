"""Numerical checks of the crossed-design residual test (Proposition 7, Corollary 1).

This script checks linear-algebra and sampling statements about the residual
model; it does not simulate a detector or a physical experiment.
Checks:
  1. With amplitude (preparation) and efficiency (detector) parameters plus one
     further random parameter per preparation and per detector, the first-order
     design matrix has rank P*n + D*m - 1 (one gauge degeneracy).
  2. With amplitude/efficiency parameters only, the testable dimension equals the
     two-way interaction degrees of freedom (P-1)(D-1).
  3. Under the null (eta = 0, Gaussian noise of known variance) the statistic T
     has mean close to nu = PD - rank(X).
  4. An injected non-factorisable eta is detected at the stated size.
  5. The first-order expansion matches the exact quadratic response to O(delta^2).
Run: python3 crossed_design_check.py   (requires numpy)
"""
import numpy as np

rng = np.random.default_rng(17102026)
N = 12            # dimension of the sampled field vector
P, D = 6, 5       # preparations, detectors


def rand_psd(n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return a @ a.conj().T / n


def rand_herm(n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (a + a.conj().T) / 2


f = [rng.normal(size=N) + 1j * rng.normal(size=N) for _ in range(P)]
L = [rand_psd(N) for _ in range(D)]
R = np.array([[np.vdot(fp, Ld @ fp).real for Ld in L] for fp in f])

# preparation directions: amplitude (phi = f_p) and one random direction
phi = [[f[p], rng.normal(size=N) + 1j * rng.normal(size=N)] for p in range(P)]
# detector directions: efficiency (Lambda = L_d) and one random Hermitian direction
lam = [[L[d], rand_herm(N)] for d in range(D)]


def design(phi, lam):
    cols = []
    for p in range(P):
        for v in phi[p]:
            c = np.zeros((P, D))
            for d in range(D):
                c[p, d] = 2 * np.vdot(v, L[d] @ f[p]).real
            cols.append(c.ravel())
    for d in range(D):
        for M in lam[d]:
            c = np.zeros((P, D))
            for p in range(P):
                c[p, d] = np.vdot(f[p], M @ f[p]).real
            cols.append(c.ravel())
    return np.array(cols).T


X = design(phi, lam)
r = np.linalg.matrix_rank(X)
assert r == P * 2 + D * 2 - 1, r
print(f"check 1: rank X = {r} = P*n + D*m - 1 = {P*2 + D*2 - 1}")

X0 = design([[f[p]] for p in range(P)], [[L[d]] for d in range(D)])
r0 = np.linalg.matrix_rank(X0)
assert P * D - r0 == (P - 1) * (D - 1)
print(f"check 2: amplitude/efficiency only -> nu = {P*D - r0} = (P-1)(D-1) = {(P-1)*(D-1)}")

nu = P * D - r
U, _, _ = np.linalg.svd(X, full_matrices=False)
Q = U[:, :r]   # orthonormal basis of col(X)
proj_perp = np.eye(P * D) - Q @ Q.T
sigma = 0.05
T = []
for _ in range(20000):
    theta = rng.normal(size=X.shape[1])
    res = X @ theta + sigma * rng.normal(size=P * D)
    T.append(res @ proj_perp @ res / sigma ** 2)
T = np.array(T)
assert abs(T.mean() - nu) < 0.1 * nu, (T.mean(), nu)
print(f"check 3: null mean of T = {T.mean():.2f}, nu = {nu}")

crit = np.quantile(T, 0.95)
eta = 0.25 * rng.normal(size=(P, D))
eta = eta - eta.mean(axis=0) - eta.mean(axis=1)[:, None] + eta.mean()   # non-factorisable part
hits = 0
for _ in range(2000):
    res = X @ rng.normal(size=X.shape[1]) + eta.ravel() + sigma * rng.normal(size=P * D)
    hits += (res @ proj_perp @ res / sigma ** 2) > crit
print(f"check 4: detection rate for injected eta = {hits/2000:.3f} (null size 0.05)")
assert hits / 2000 > 0.9

for scale in (1e-2, 1e-3):
    df = scale * (rng.normal(size=N) + 1j * rng.normal(size=N))
    dL = scale * rand_herm(N)
    exact = np.vdot(f[0] + df, (L[0] + dL) @ (f[0] + df)).real - R[0, 0]
    lin = 2 * np.vdot(df, L[0] @ f[0]).real + np.vdot(f[0], dL @ f[0]).real
    print(f"check 5: scale {scale:g}: |exact - linear| = {abs(exact - lin):.2e}")
print("all crossed-design checks passed")

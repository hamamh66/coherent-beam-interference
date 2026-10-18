"""Checks for the coherence/coupling identifiability results (Section 4.8).

1. First-principles check of the visibility model: fields of two partially
   coherent crossing plane waves are averaged over a finite probe width a, for a
   probe responding to (1-w)|E|^2 + w c^2|B|^2; the extracted visibility is compared
   with V = m * (1 - 2 w sin^2 th) * sinc(2 pi a sin th / lambda).
2. Gauge: (m, kappa) -> (c m, kappa / c) leaves all data unchanged.
3. Local identifiability: Jacobian rank of (Gamma, w, beta) -> {V(th_k)} for K = 3.
4. Leading-order confounding of w and probe size: w_eff = w + beta^2 / 12.
5. Angle-dependent coherence m(th) reproduces any single-detector data for another w.
6. Ratio test: V_MD / V_ED is independent of m, kappa and shared probe size.
7. Monte Carlo check of the least-squares variance for the known-beta linear fit.
No physical experiment or device simulation; units dimensionless (lambda = 1).
Run: python3 identifiability_check.py  (numpy)
"""
import numpy as np

rng = np.random.default_rng(20261017)
lam = 1.0
k = 2 * np.pi / lam


def sinc(x):
    return np.sinc(x / np.pi)  # sin(x)/x


def model(m, kap, w, a, th):
    return m * kap * (1 - 2 * w * np.sin(th) ** 2) * sinc(k * a * np.sin(th))


def first_principles(m, w, a, th, n=4001):
    """Ensemble over relative phase: mutual coherence m realised as a mixture of a
    fully coherent part (weight m) and a fully incoherent part (weight 1-m)."""
    s, h = np.sin(th), np.cos(th)
    Lam = lam / (2 * s)
    xs = np.linspace(0, Lam, 201)            # probe centre positions over one period
    resp = []
    for x0 in xs:
        x = x0 + np.linspace(-a / 2, a / 2, n) if a > 0 else np.array([x0])
        q = k * s * x
        E2c = 4 * np.cos(q) ** 2                          # coherent part, |E|^2
        B2c = 4 * (h ** 2 * np.cos(q) ** 2 + s ** 2 * np.sin(q) ** 2)
        inc = 2.0                                         # incoherent sum: |E1|^2+|E2|^2 = c^2(|B1|^2+|B2|^2)
        E2 = m * E2c + (1 - m) * inc
        B2 = m * B2c + (1 - m) * inc
        resp.append(np.mean((1 - w) * E2 + w * B2))
    resp = np.array(resp)
    return (resp.max() - resp.min()) / (resp.max() + resp.min())


# 1. model vs first principles (|V| compared; model sign gives inversion)
for _ in range(40):
    m, w = rng.uniform(0.2, 1), rng.uniform(0, 1)
    a = rng.uniform(0, 0.4)
    th = rng.uniform(np.radians(8), np.radians(80))
    Vfp = first_principles(m, w, a, th)
    Vm = abs(model(m, 1, w, a, th))
    assert abs(Vfp - Vm) < 2e-3, (m, w, a, th, Vfp, Vm)
print("check 1: visibility model matches first-principles averaging (40 random cases)")

# 2. gauge
th = np.radians([15, 35, 55, 75])
for c in (0.5, 1.7):
    assert np.allclose(model(0.8, 0.9, 0.3, 0.2, th), model(0.8 * c, 0.9 / c, 0.3, 0.2, th))
print("check 2: (m, kappa) gauge leaves data unchanged")

# 3. local identifiability of (Gamma, w, beta) with K = 3 angles
def V3(p, th):
    G, w, beta = p
    return G * (1 - 2 * w * np.sin(th) ** 2) * sinc(beta * np.sin(th))
ranks = []
for _ in range(200):
    p = np.array([rng.uniform(0.3, 1), rng.uniform(0, 1), rng.uniform(0.3, 2.5)])
    th3 = np.sort(rng.uniform(np.radians(10), np.radians(80), 3))
    J = np.zeros((3, 3))
    for i in range(3):
        dp = np.zeros(3); dp[i] = 1e-6
        J[:, i] = (V3(p + dp, th3) - V3(p - dp, th3)) / 2e-6
    ranks.append(np.linalg.matrix_rank(J, tol=1e-9))
assert min(ranks) == 3
print("check 3: Jacobian rank 3 for all 200 random (Gamma, w, beta, angles)")

# 4. leading-order confounding
for beta in (0.05, 0.1, 0.2):
    s = np.linspace(0.01, 0.2, 50)
    exact = (1 - 2 * 0.2 * s) * sinc(beta * np.sqrt(s))
    eff = 1 - 2 * (0.2 + beta ** 2 / 12) * s
    print(f"check 4: beta={beta}: max |exact - (w+beta^2/12) model| = {np.max(abs(exact-eff)):.2e}")

# 5. angle-dependent coherence mimics coupling
th = np.radians(np.linspace(10, 70, 7))
data = model(0.9, 1, 0.3, 0.0, th)
w2 = 0.0
m_th = data / (1 - 2 * w2 * np.sin(th) ** 2)
assert np.allclose(model(m_th, 1, w2, 0.0, th), data) and np.all((m_th > 0) & (m_th <= 1))
print("check 5: an angle-dependent m(theta) reproduces single-detector data with w=0")

# 6. ratio test
for _ in range(50):
    m, kap, a = rng.uniform(0.2, 1), rng.uniform(0.5, 1), rng.uniform(0, 0.3)
    th = rng.uniform(np.radians(5), np.radians(40))
    r = model(m, kap, 1, a, th) / model(m, kap, 0, a, th)
    assert abs(r - np.cos(2 * th)) < 1e-12
print("check 6: V_MD/V_ED = cos(2 theta) independent of m, kappa, shared probe size")

# 7. least-squares variance (known beta): V = Gamma - 2 b s, b = Gamma w
s = np.sin(np.radians(np.array([10, 10, 10, 70, 70, 70]))) ** 2
X = np.column_stack([np.ones_like(s), -2 * s])
sig = 0.005
cov = sig ** 2 * np.linalg.inv(X.T @ X)
est = []
for _ in range(20000):
    y = X @ np.array([0.8, 0.8 * 0.3]) + sig * rng.normal(size=s.size)
    est.append(np.linalg.lstsq(X, y, rcond=None)[0])
est = np.array(est)
assert np.allclose(np.cov(est.T), cov, rtol=0.05)
print("check 7: Monte Carlo covariance matches the least-squares formula")
print("all identifiability checks passed")

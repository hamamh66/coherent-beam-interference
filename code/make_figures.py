"""Analytical figures for the crossing-plane-wave benchmark (no measured data).

Dimensionless units: c = epsilon0 = mu0 = A = 1. Quantities follow the closed-form
expressions of the Results section:
  electric response      |E|^2      = 4 cos^2 q
  magnetic response      c^2|B|^2   = 4 (h^2 cos^2 q + s^2 sin^2 q)
  energy density         u          = 1 + h^2 cos 2q
  longitudinal flux      S_z        = 2 h cos^2 q      (mean transverse flux is zero)
with s = sin(theta), h = cos(theta), q = k s x (phi = 0), fringe period Lambda = lambda/(2 s).
Run: python3 make_figures.py  (writes ../Figures/*.pdf). Requires numpy and matplotlib.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.linewidth": 0.6,
                     "lines.linewidth": 1.2, "legend.frameon": False})
C = {"E": "#1f4e9c", "B": "#c0392b", "u": "#222222", "S": "#2e8b57"}


def fields(theta, xi):
    """xi = x / Lambda, so q = pi * xi."""
    s, h = np.sin(theta), np.cos(theta)
    q = np.pi * xi
    E2 = 4 * np.cos(q) ** 2
    B2 = 4 * (h ** 2 * np.cos(q) ** 2 + s ** 2 * np.sin(q) ** 2)
    u = 1 + h ** 2 * np.cos(2 * q)
    Sz = 2 * h * np.cos(q) ** 2
    return E2, B2, u, Sz


# Figure 1: profiles across two fringe periods for three half-angles.
xi = np.linspace(-1, 1, 801)
fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=True)
for ax, deg in zip(axes, (30, 45, 70)):
    th = np.radians(deg)
    E2, B2, u, Sz = fields(th, xi)
    # normalise each response by its own spatial mean to compare fringe contrast
    ax.plot(xi, E2 / E2.mean(), color=C["E"], label=r"$|\mathcal{E}|^2$ (electric)")
    ax.plot(xi, B2 / B2.mean(), color=C["B"], ls="--", label=r"$c^2|\mathcal{B}|^2$ (magnetic)")
    ax.plot(xi, u / u.mean(), color=C["u"], ls=":", label=r"$\bar u$ (energy density)")
    ax.set_title(rf"$\theta={deg}^\circ$", fontsize=8)
    ax.set_xlabel(r"$x/\Lambda$")
    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 2.15)
    for n in (-0.5, 0.5):
        ax.axvline(n, color="0.8", lw=0.5, zorder=0)
axes[0].set_ylabel("response / spatial mean")
h_, l_ = axes[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="upper center", ncol=3, fontsize=7)
fig.tight_layout(rect=(0, 0, 1, 0.88))
fig.savefig(os.path.join(OUT, "fig_benchmark_profiles.pdf"))

# Figure 2: fringe visibility of each functional versus half-angle.
th = np.linspace(0.0, np.pi / 2, 400)
deg = np.degrees(th)
VE = np.ones_like(th)
VB = np.cos(2 * th)          # signed: negative means contrast inverted relative to electric
Vu = np.cos(th) ** 2
fig, ax = plt.subplots(figsize=(3.4, 2.3))
ax.plot(deg, VE, color=C["E"], label=r"electric $V_E$ (also $\bar S_z$)")
ax.plot(deg, VB, color=C["B"], ls="--", label=r"magnetic, signed: $\cos 2\theta$")
ax.plot(deg, Vu, color=C["u"], ls=":", label=r"energy density $V_u=\cos^2\theta$")
ax.axhline(0, color="0.6", lw=0.5)
ax.axvline(45, color="0.8", lw=0.5)
ax.set_xlabel(r"half crossing angle $\theta$ (deg)")
ax.set_ylabel("signed fringe visibility")
ax.set_xlim(0, 90)
ax.set_ylim(-1.05, 1.1)
ax.legend(loc="lower left", fontsize=6.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_visibility_angle.pdf"))
print("figures written to", os.path.abspath(OUT))

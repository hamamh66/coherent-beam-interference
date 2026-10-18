"""Figure for the coherence/coupling identifiability results (analytical model; no data).

(a) Two situations indistinguishable at one crossing angle with one detector:
    A: m = 0.7, electric probe (w = 0);  B: m = 1, mixed probe with w chosen so that
    V_A = V_B at theta0 = 40 deg. An angle scan separates them.
(b) Two-channel ratio test on a shared probe (electric- and magnetic-dipole channels):
    a loss of coherence leaves V_MD / V_ED = cos 2 theta unchanged; an admixture of
    electric coupling into the magnetic channel (w: 1 -> 0.85) changes it.
Run: python3 make_identifiability_figure.py  (numpy, matplotlib)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures")
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                     "font.serif": ["STIXGeneral", "DejaVu Serif"], "font.size": 8,
                     "axes.linewidth": 0.6, "lines.linewidth": 1.3, "legend.frameon": False})
th0 = np.radians(40)
wB = 0.3 / (2 * np.sin(th0) ** 2)
th = np.radians(np.linspace(1, 75, 300))
s = np.sin(th) ** 2
VA = 0.7 * np.ones_like(th)
VB = 1 - 2 * wB * s

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.6))
ax1.plot(s, VA, color="#1f4e9c", label=r"A: $m=0.7$, electric probe ($w=0$)")
ax1.plot(s, VB, color="#c0392b", ls="--", label=rf"B: $m=1$, mixed probe ($w={wB:.2f}$)")
ax1.axvline(np.sin(th0) ** 2, color="0.75", lw=0.6)
ax1.annotate(r"$\theta_0=40^\circ$", (np.sin(th0) ** 2, 0.98), xytext=(4, 0),
             textcoords="offset points", fontsize=7, color="0.35")
ax1.set_xlabel(r"$\sin^2\theta$")
ax1.set_ylabel("fringe visibility $V$")
ax1.set_ylim(0.3, 1.05)
ax1.set_title("(a) one detector, angle scan", fontsize=8)
ax1.legend(loc="lower left", fontsize=6.5)

th2 = np.radians(np.linspace(1, 44, 300))
r0 = np.cos(2 * th2)
r_drift = (1 - 2 * 0.85 * np.sin(th2) ** 2) / 1.0
ax2.plot(np.degrees(th2), r0, color="#1f4e9c", label=r"reference, or any loss of coherence $m$")
ax2.plot(np.degrees(th2), r_drift, color="#c0392b", ls="--", label=r"coupling change, $w_{\rm MD}: 1\to0.85$")
ax2.set_xlabel(r"half crossing angle $\theta$ (deg)")
ax2.set_ylabel(r"$V_{\rm MD}/V_{\rm ED}$")
ax2.set_title("(b) two channels of one probe, ratio test", fontsize=8)
ax2.legend(loc="lower left", fontsize=6.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Figure_Identifiability.pdf"))
print("wB =", round(wB, 4))

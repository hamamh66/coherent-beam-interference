"""Amplitude-addition figure for a selected detection channel (schematic, exact algebra).

Panels show the signed single-excitation detection amplitude for a channel
b = kappa1 a1 + kappa2 a2 sampled across a detector plane,
    A_j(x, y) = kappa c_j u(x, y),   A = A_1 + A_2,
with the same mode profile u for both inputs (balanced channel, kappa1 = kappa2).
Row 1: c2 = c1  -> A = 2 A_1 (rate |A|^2 = 4 |A_1|^2).
Row 2: c2 = -c1 -> A = 0 everywhere: the state is dark to this channel, while the
excitation is carried by the orthogonal combination of the modes.
Run: python3 make_amplitude_figure.py  (numpy, matplotlib)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures")
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                     "font.serif": ["STIXGeneral", "DejaVu Serif"]})

x = np.linspace(-3, 3, 121)
X, Y = np.meshgrid(x, x)
u = (1 - X ** 2) * np.exp(-(X ** 2 + Y ** 2) / 1.6)       # mode profile: central lobe, side lobes of opposite sign
c1 = 1.0  # amplitudes shown in units of kappa*c1 (normalisation omitted for display)
rows = [("Constructive: $c_2=+c_1$", [c1 * u, c1 * u]),
        ("Destructive: $c_2=-c_1$", [c1 * u, -c1 * u])]
norm = colors.TwoSlopeNorm(vmin=-1.0, vcenter=0, vmax=1.0)
cmap = cm.RdBu  # red negative, blue positive

fig = plt.figure(figsize=(7.2, 4.4))
titles = [r"$A_1=\kappa c_1 u(x,y)$", r"$A_2=\kappa c_2 u(x,y)$", r"$A=A_1+A_2$"]
for r, (rt, (a1, a2)) in enumerate(rows):
    total = a1 + a2
    for k, Z in enumerate((a1, a2, total)):
        ax = fig.add_axes([-0.01 + 0.33 * k, 0.47 - 0.45 * r, 0.36, 0.50], projection="3d")
        ax.plot_surface(X, Y, Z, facecolors=cmap(norm(Z)), rstride=2, cstride=2,
                        linewidth=0.2, edgecolor=(1, 1, 1, 0.45), antialiased=True, shade=False)
        ax.set_zlim(-2.0, 2.0)
        ax.view_init(elev=20, azim=-62)
        ax.set_axis_off()
        ax.set_box_aspect((1, 1, 0.75), zoom=1.25)
        fig.text(0.17 + 0.33 * k, 0.895 - 0.45 * r, titles[k], ha="center", fontsize=9)
        if k < 2:
            fig.text(0.345 + 0.33 * k, 0.70 - 0.45 * r, "+" if k == 0 else "=",
                     ha="center", va="center", fontsize=16, color="0.35")
    fig.text(0.02, 0.955 - 0.45 * r, rt, fontsize=10.5, weight="bold", color="#193b53")
fig.text(0.83, 0.075, "rate $\\propto |A|^2=0$:\ndark to this channel;\nexcitation remains in the\northogonal mode combination",
         ha="center", fontsize=7.5, color="#a86220")
fig.text(0.83, 0.535, "rate $\\propto |A|^2=4|A_1|^2$", ha="center", fontsize=7.5, color="#193b53")
fig.text(0.5, 0.012, "Signed detection amplitudes for one excitation (blue $>0$, red $<0$); not field energy and not photon energy",
         ha="center", fontsize=7.5, color="0.35")
fig.savefig(os.path.join(OUT, "Figure_Amplitude_Addition.pdf"))
fig.savefig(os.path.join(OUT, "Figure_Amplitude_Addition.png"), dpi=200)
assert np.allclose(rows[0][1][0] + rows[0][1][1], 2 * rows[0][1][0])
assert np.allclose(rows[1][1][0] + rows[1][1][1], 0)
print("written", os.path.abspath(OUT))

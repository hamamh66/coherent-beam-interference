"""Framework overview figure (schematic; no data).

Regenerates Figures/Figure_1_Framework_Overview.pdf and .svg with matplotlib only
(serif text, mathtext for equations), so it runs on any OS without external fonts.
Run: python3 make_framework_overview.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Figures")
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                     "font.serif": ["STIXGeneral", "DejaVu Serif"], "svg.fonttype": "path"})
NAVY, GREY, ORANGE, FILL, EDGE = "#193b53", "#52616d", "#a86220", "#f5f8fa", "#b6c6d0"

fig = plt.figure(figsize=(7.2, 4.3))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 720)
ax.set_ylim(0, 430)
ax.axis("off")


def box(x, y, w, h):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=8",
                                fc=FILL, ec=EDGE, lw=0.8))


def arrow(p, q, color=GREY, conn="arc3"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=10,
                                 color=color, lw=1.0, connectionstyle=conn))


def t(x, y, s, size=9, color=NAVY, weight="normal"):
    ax.text(x, y, s, fontsize=size, color=color, weight=weight, ha="left", va="center")


t(20, 410, "(a)  Preparation, interaction, and recorded response", 11, weight="bold")
for x in (20, 260, 500):
    box(x, 215, 200, 170)

t(32, 365, "PREPARED LIGHT", 9.5, "#246b8c", "bold")
t(32, 340, r"Vector fields or quantum state $\rho$", 8.5)
t(32, 312, r"$\mathbf{E}=\mathbf{E}_1+\mathbf{E}_2,\quad \mathbf{B}=\mathbf{B}_1+\mathbf{B}_2$", 10)
t(32, 278, "Phase, polarization, coherence;", 8)
t(32, 262, "mode structure and correlations", 8)

t(272, 365, "CHOSEN INTERACTION", 9.5, "#277b70", "bold")
t(272, 340, "Calibrated probe and its boundaries", 8.5)
t(272, 312, "Electric / magnetic sensitivity;", 8)
t(272, 296, "spatial and spectral selection", 8)
t(272, 262, "Absorption, scattering, or dispersion", 8)

t(512, 365, "MEASUREMENT RECORD", 9.5, "#246b8c", "bold")
t(512, 340, "Counts, absorbed power, or force", 8.5)
t(512, 312, "Fringe phase and visibility depend", 8)
t(512, 296, "on preparation and response", 8)
t(512, 262, "A dark channel need not mean", 8)
t(512, 246, "zero electromagnetic energy", 8)

arrow((221, 300), (259, 300))
arrow((461, 300), (499, 300))
# Back-action loops from the probe onto the incident field at the probe
# (the arrow entering the interaction box), not onto the source preparation.
arrow((360, 215), (240, 296), color=ORANGE, conn="arc3,rad=-0.6")
t(395, 195, "Back-action on the field in the", 8.5, ORANGE)
t(395, 180, "interaction region: reflection,", 8.5, ORANGE)
t(395, 165, "rescattering, absorption, mode changes", 8.5, ORANGE)
t(395, 150, "(source preparation unchanged)", 8, ORANGE)

t(20, 122, "(b)  Quantum description of a selected detection channel", 11, weight="bold")
box(20, 14, 680, 90)
t(36, 80, r"$\hat b=\kappa_1\hat a_1+\kappa_2\hat a_2,\qquad R_b=\mathcal{K}\,\langle \hat b^{\dagger}\hat b\rangle$", 12)
t(36, 42, r"$R_b/\mathcal{K}=|\kappa_1|^2\langle\hat n_1\rangle+|\kappa_2|^2\langle\hat n_2\rangle"
          r"+2\,\mathrm{Re}\left[\kappa_1^{*}\kappa_2\langle\hat a_1^{\dagger}\hat a_2\rangle\right]$", 12)

for ext in ("pdf", "svg"):
    fig.savefig(os.path.join(OUT, f"Figure_1_Framework_Overview.{ext}"))
print("written", os.path.abspath(OUT))

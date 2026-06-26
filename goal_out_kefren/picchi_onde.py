#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
picchi_onde.py — disegna SOLO i punti di picco delle forme d'onda della
tomografia verticale (anti-FT delle 6 armoniche, base 268x276 m, prof. 0..1000 m).
Per ogni pixel cerca i massimi locali di |w(z)| (gli anti-nodi dell'onda
stazionaria) e li mostra come nuvola di punti 3D, colorati per ampiezza.
"""
import sys, os
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = "rettangolo"
ZTOP, ZBOT = 210.0, -200.0  # finestra quota reale s.l.m. [m]
CONF = dict(scrollZoom=True, displayModeBar=True)   # zoom abilitato
d = np.load(os.path.join(OUT, "tomografia_verticale.npz"), allow_pickle=True)
w = d["w"]; z = d["z"]; x = d["x"]; y = d["y"]
Lx = float(d["Lx"]); Ly = float(d["Ly"]); ZD = float(d["zdepth"])
z0 = d["z0"]                        # quota reale del terreno (DEM) [m s.l.m.]
ny, nx, nz = w.shape

px, py, pz, pv = [], [], [], []
for i in range(ny):
    for j in range(nx):
        a = np.abs(w[i, j])
        idx, _ = find_peaks(a)                  # massimi locali = anti-nodi
        for k in idx:
            px.append(x[j]); py.append(y[i]); pz.append(z0[i, j] - z[k]); pv.append(w[i, j, k])
px, py, pz, pv = map(np.array, (px, py, pz, pv))
m = (pz >= ZBOT) & (pz <= ZTOP)                       # tieni la finestra +200/-200 m
px, py, pz, pv = px[m], py[m], pz[m], pv[m]
print(f"picchi trovati: {len(px)}  (nella finestra +200/-200 m)")

# ---- nuvola 3D dei soli punti di picco -------------------------------------
Xg, Yg = np.meshgrid(x, y)
ground = go.Scatter3d(x=Xg.ravel(), y=Yg.ravel(), z=z0.ravel(), mode="markers",
                      marker=dict(size=2.5, color="red"), name="altezza suolo (DEM)")
fig = go.Figure([go.Scatter3d(
    x=px, y=py, z=pz, mode="markers", name="picchi",
    marker=dict(size=2.5, color=pv, colorscale="Turbo", opacity=0.85,
                colorbar=dict(title="ampiezza al picco"))), ground])
fig.update_layout(
    title=f"Punti di picco delle forme d'onda (anti-nodi) — base {Lx:.0f}×{Ly:.0f} m, "
          f"profondità 0..{ZD:.0f} m",
    scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
               zaxis_title="quota reale s.l.m. [m]", zaxis=dict(range=[ZBOT, ZTOP]),
               aspectmode="manual",
               aspectratio=dict(x=1.0, y=Ly / Lx, z=(ZTOP - ZBOT) / Lx)))
fig.write_html(os.path.join(OUT, "picchi_onde.html"), include_plotlyjs="cdn",
               config=CONF)

# anteprima statica
figm = plt.figure(figsize=(8, 9)); ax = figm.add_subplot(111, projection="3d")
ax.scatter(Xg.ravel(), Yg.ravel(), z0.ravel(), c="red", s=5, depthshade=False)
s = ax.scatter(px, py, pz, c=pv, cmap="turbo", s=4)
ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZTOP - ZBOT))
ax.set_title("Solo punti di picco delle forme d'onda (anti-nodi)")
figm.colorbar(s, ax=ax, shrink=0.5, label="ampiezza")
figm.savefig(os.path.join(OUT, "picchi_onde.png"), dpi=130, bbox_inches="tight")
plt.close(figm)

np.savez_compressed(os.path.join(OUT, "picchi_onde.npz"),
                    px=px, py=py, pz=pz, pv=pv)
print("salvati in", OUT + "/: picchi_onde.html, picchi_onde.png, picchi_onde.npz")

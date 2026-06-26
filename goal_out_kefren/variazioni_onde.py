#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variazioni_onde.py — come picchi_onde.py ma sui punti di VARIAZIONE DI FREQUENZA
delle forme d'onda della tomografia verticale (anti-FT delle 6 armoniche,
base 268x276 m, profondita' 0..1000 m).
Per ogni pixel: frequenza istantanea lungo z (Hilbert) -> |d(freq)/dz|;
i punti con variazione marcata = dove cambia la densita' del mezzo (Biondi-Malanga).
"""
import sys, os
import numpy as np
from scipy.signal import hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = "rettangolo"
DMIN, DMAX = 0.03, 0.06     # banda densità: tieni 0.03 < |Δfreq| <= 0.06
ZTOP, ZBOT = 210.0, -200.0  # finestra quota reale s.l.m. [m]
CONF = dict(scrollZoom=True, displayModeBar=True)   # zoom abilitato
d = np.load(os.path.join(OUT, "tomografia_verticale.npz"), allow_pickle=True)
w = d["w"]; z = d["z"]; x = d["x"]; y = d["y"]
Lx = float(d["Lx"]); Ly = float(d["Ly"]); ZD = float(d["zdepth"])
z0 = d["z0"]                        # quota reale del terreno (DEM) [m s.l.m.]
ny, nx, nz = w.shape
dz = z[1] - z[0]

# frequenza istantanea lungo la profondita' (z)
an = hilbert(w, axis=2)
inst_f = np.gradient(np.unwrap(np.angle(an), axis=2), dz, axis=2) / (2 * np.pi)
dfreq = np.abs(np.gradient(inst_f, dz, axis=2))     # variazione di frequenza

ii, jj, kk = np.where((dfreq > DMIN) & (dfreq <= DMAX))   # banda densità 0.03..0.06
px = x[jj]; py = y[ii]; pz = z0[ii, jj] - z[kk]; pv = dfreq[ii, jj, kk]   # quota reale
m = (pz >= ZBOT) & (pz <= ZTOP)                       # tieni la finestra +210/-200 m
px, py, pz, pv = px[m], py[m], pz[m], pv[m]
print(f"banda densità {DMIN}<|Δfreq|<={DMAX}; punti tenuti: {len(px)} su {dfreq.size} "
      f"(max densità={dfreq.max():.4f})")

# ---- nuvola 3D dei punti di variazione di frequenza ------------------------
Xg, Yg = np.meshgrid(x, y)
ground = go.Scatter3d(x=Xg.ravel(), y=Yg.ravel(), z=z0.ravel(), mode="markers",
                      marker=dict(size=2.5, color="red"), name="altezza suolo (DEM)")
fig = go.Figure([go.Scatter3d(
    x=px, y=py, z=pz, mode="markers", name="variazioni freq.",
    marker=dict(size=2.5, color=pv, colorscale="Inferno", opacity=0.85,
                colorbar=dict(title="|Δfreq| (densità)"))), ground])
fig.update_layout(
    title=f"Variazioni di frequenza (=> densità) delle forme d'onda — "
          f"base {Lx:.0f}×{Ly:.0f} m, profondità 0..{ZD:.0f} m",
    scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
               zaxis_title="quota reale s.l.m. [m]", zaxis=dict(range=[ZBOT, ZTOP]),
               aspectmode="manual",
               aspectratio=dict(x=1.0, y=Ly / Lx, z=(ZTOP - ZBOT) / Lx)))
fig.write_html(os.path.join(OUT, "variazioni_onde.html"), include_plotlyjs="cdn",
               config=CONF)

# anteprima statica
figm = plt.figure(figsize=(8, 9)); ax = figm.add_subplot(111, projection="3d")
ax.scatter(Xg.ravel(), Yg.ravel(), z0.ravel(), c="red", s=5, depthshade=False)
s = ax.scatter(px, py, pz, c=pv, cmap="inferno", s=4)
ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZTOP - ZBOT))
ax.set_title("Punti di variazione di frequenza (=> densità) delle forme d'onda")
figm.colorbar(s, ax=ax, shrink=0.5, label="|Δfreq|")
figm.savefig(os.path.join(OUT, "variazioni_onde.png"), dpi=130, bbox_inches="tight")
plt.close(figm)

np.savez_compressed(os.path.join(OUT, "variazioni_onde.npz"),
                    px=px, py=py, pz=pz, pv=pv, dmin=DMIN, dmax=DMAX)
print("salvati in", OUT + "/: variazioni_onde.html, variazioni_onde.png, variazioni_onde.npz")

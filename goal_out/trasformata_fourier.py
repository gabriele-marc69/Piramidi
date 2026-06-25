#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
trasformata_fourier.py — trasformata di Fourier (spaziale) della forma d'onda W
costruita con le 6 armoniche nel rettangolo 268x276 m.
W(x,y) = somma_{k=1..6} D_k(x,y) sin(2*pi*k*x/Lx)  ->  FFT lungo x.
Lo spettro |F(kx, y)| deve mostrare 6 creste alle frequenze k/Lx (k=1..6).

Grafici (cartella rettangolo/):
 - trasformata_3d.html/.png : superficie 3D dello spettro |F| su (freq kx, y N-S)
 - trasformata_spettro.png  : spettro medio 1D con i 6 picchi armonici marcati
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = "rettangolo"
d = np.load(os.path.join(OUT, "forma_donda.npz"), allow_pickle=True)
W = d["W"]                      # (ny, NXX)
x = d["x"]; y = d["y"]; Lx = float(d["Lx"]); Ly = float(d["Ly"])
ny, nxx = W.shape
dx = x[1] - x[0]

# FFT spaziale lungo x con zero-padding -> spettro interpolato (forma d'onda liscia)
NPAD = 8192
F = np.fft.rfft(W, n=NPAD, axis=1)              # (ny, nf) complesso
mag = np.abs(F) * (2.0 / nxx)                   # ampiezza fisica [mm]
freq = np.fft.rfftfreq(NPAD, d=dx)              # cicli/m
harm = np.arange(1, 7) / Lx                     # le 6 frequenze armoniche attese

# limita la vista alle prime ~8 armoniche
fmax = 8.0 / Lx
sel = freq <= fmax
freqs = freq[sel]; M = mag[:, sel]
print(f"spettro: {M.shape}  freq 0..{freqs[-1]*1000:.2f} cicli/km ; "
      f"6 armoniche a k/Lx = {np.round(harm*1000,3)} cicli/km")

# ---- 1) superficie 3D dello spettro: X=freq, Y=N-S, Z=|F| ------------------
Fg, Yg = np.meshgrid(freqs, y)
fig = go.Figure(go.Surface(x=Fg * 1000, y=Yg, z=M, colorscale="Viridis",
                           colorbar=dict(title="|F| [mm]")))
# linee verticali sulle 6 armoniche attese
zmax = float(M.max())
for k, hk in enumerate(harm, 1):
    fig.add_trace(go.Scatter3d(x=[hk * 1000, hk * 1000], y=[y[0], y[-1]],
                               z=[zmax, zmax], mode="lines",
                               line=dict(color="crimson", width=4),
                               name=f"armonica {k}", showlegend=(k == 1)))
fig.update_layout(
    title="Trasformata di Fourier della forma d'onda (6 armoniche) — rettangolo 268×276 m",
    scene=dict(xaxis_title="frequenza spaziale kx [cicli/km]",
               yaxis_title="y N-S [m]", zaxis_title="|F(kx,y)| [mm]",
               aspectratio=dict(x=1.2, y=0.8, z=0.5)))
fig.write_html(os.path.join(OUT, "trasformata_3d.html"), include_plotlyjs="cdn")

# statico png della superficie
figm = plt.figure(figsize=(10, 7)); ax = figm.add_subplot(111, projection="3d")
ax.plot_surface(Fg * 1000, Yg, M, cmap="viridis", linewidth=0, antialiased=True)
for hk in harm:
    ax.plot([hk * 1000, hk * 1000], [y[0], y[-1]], [zmax, zmax], color="crimson", lw=2)
ax.set_xlabel("kx [cicli/km]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("|F| [mm]")
ax.set_title("Trasformata di Fourier della forma d'onda (creste = 6 armoniche)")
figm.savefig(os.path.join(OUT, "trasformata_3d.png"), dpi=130, bbox_inches="tight")
plt.close(figm)

# ---- 2) spettro medio 1D con i 6 picchi armonici ---------------------------
mean_spec = M.mean(0)
fig2, ax = plt.subplots(figsize=(9, 5))
ax.plot(freqs * 1000, mean_spec, color="navy", lw=1.8)
for k, hk in enumerate(harm, 1):
    ax.axvline(hk * 1000, color="crimson", ls="--", lw=1)
    ax.text(hk * 1000, mean_spec.max() * (0.9 - 0.06 * (k % 2)), f"k={k}",
            color="crimson", ha="center", fontsize=8)
ax.set_xlabel("frequenza spaziale kx [cicli/km]")
ax.set_ylabel("ampiezza media |F| [mm]")
ax.set_title("Spettro di Fourier della forma d'onda — 6 armoniche nel rettangolo")
fig2.tight_layout(); fig2.savefig(os.path.join(OUT, "trasformata_spettro.png"), dpi=130)
plt.close(fig2)

print("salvati in", OUT + "/: trasformata_3d.html/.png, trasformata_spettro.png")

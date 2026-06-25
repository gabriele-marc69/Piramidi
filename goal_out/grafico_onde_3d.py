#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grafico_onde_3d.py — TUTTE le funzioni d'onda di Fourier dello step 4
sovrapposte in un unico grafico 3D: una linea per ogni pixel.
Assi: X = s (range stirato 0..1000 m), Y = indice pixel (raster az,rng),
Z = ampiezza f(s) [mm]. Colore = ampiezza.
"""
import sys, argparse
import numpy as np
import plotly.graph_objects as go
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ap = argparse.ArgumentParser()
ap.add_argument("--npz", default="s1d/array3d_fourier.npz")
ap.add_argument("--out", default="s1d/step4_onde_3d.html")
ap.add_argument("--sub", type=int, default=2, help="sottocampiona s (1=tutti i 400)")
a = ap.parse_args()

d = np.load(a.npz, allow_pickle=True)
F = d["fourier"]                 # (ny, nx, Ns)
s = d["s"]
ny, nx, ns = F.shape
ss = s[::a.sub]
flat = F.reshape(ny * nx, ns)[:, ::a.sub]    # (Npix, len(ss))
npix, nsub = flat.shape

# una sola traccia: ogni pixel = un segmento, separato dal successivo con NaN
nan = np.full(npix, np.nan)
X = np.concatenate([np.tile(ss, (npix, 1)),
                    nan[:, None]], axis=1).ravel()           # s ripetuto + gap
Y = np.concatenate([np.repeat(np.arange(npix), nsub).reshape(npix, nsub),
                    nan[:, None]], axis=1).ravel()           # indice pixel
Z = np.concatenate([flat, nan[:, None]], axis=1).ravel()     # ampiezza
col = Z.copy()

fig = go.Figure(go.Scatter3d(
    x=X, y=Y, z=Z, mode="lines",
    line=dict(width=2, color=col, colorscale="Turbo",
              colorbar=dict(title="f(s) [mm]")),
    connectgaps=False, opacity=0.7))
fig.update_layout(
    title=f"Step 4 — {npix} funzioni d'onda di Fourier sovrapposte (una per pixel)",
    scene=dict(xaxis_title="s — range stirato [m] (0..1000)",
               yaxis_title="indice pixel (raster az, rng)",
               zaxis_title="ampiezza f(s) [mm]",
               aspectratio=dict(x=1, y=1.4, z=0.5)))
fig.write_html(a.out, include_plotlyjs="cdn")
print(f"salvato {a.out}  ({npix} onde, {nsub} campioni s)")

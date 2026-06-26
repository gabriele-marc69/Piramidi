#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grafico_onde_pixel.py — disegna le funzioni d'onda di Fourier dello step 4
UNA PER UNA, una per ogni pixel del box. Animazione: ogni frame mostra
la funzione d'onda f(s), s in 0..1000 m, di un singolo pixel.
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
ap.add_argument("--out", default="s1d/step4_onde_pixel.html")
ap.add_argument("--sub", type=int, default=2, help="sottocampiona s (1=tutti i 400)")
a = ap.parse_args()

d = np.load(a.npz, allow_pickle=True)
F = d["fourier"]                 # (ny, nx, Ns)
s = d["s"]; x = d["x"]; y = d["y"]
ny, nx, ns = F.shape
ss = s[::a.sub]
zmin, zmax = float(F.min()), float(F.max())

# un frame per pixel, in ordine raster (i azimuth, j range)
frames = []
for i in range(ny):
    for j in range(nx):
        frames.append(go.Frame(
            data=[go.Scatter(x=ss, y=F[i, j, ::a.sub], mode="lines",
                             line=dict(color="crimson", width=2))],
            name=f"{i}_{j}",
            layout=go.Layout(title=f"Funzione d'onda — pixel (az {i}, rng {j})  "
                                   f"x={x[j]:.0f} m  y={y[i]:.0f} m")))

f0 = frames[0].data[0]
fig = go.Figure(data=[f0], frames=frames)
fig.update_layout(
    title="Step 4 — funzioni d'onda di Fourier, una per pixel (az 0, rng 0)",
    xaxis_title="s — range stirato [m] (0..1000)",
    yaxis_title="ampiezza f(s) [mm]",
    yaxis=dict(range=[zmin, zmax]),
    updatemenus=[dict(type="buttons", showactive=False, x=0.05, y=1.12,
                      buttons=[dict(label="▶ Play", method="animate",
                                    args=[None, dict(frame=dict(duration=60, redraw=True),
                                                     fromcurrent=True)]),
                               dict(label="⏸ Pause", method="animate",
                                    args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                       mode="immediate")])])],
    sliders=[dict(currentvalue=dict(prefix="pixel "),
                  steps=[dict(method="animate", label=fr.name,
                              args=[[fr.name], dict(mode="immediate",
                                                    frame=dict(duration=0, redraw=True))])
                         for fr in frames])])
fig.write_html(a.out, include_plotlyjs="cdn")
print(f"salvato {a.out}  ({len(frames)} pixel = frame, {len(ss)} campioni s)")

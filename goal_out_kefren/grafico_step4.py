#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
grafico_step4.py — grafico 3D dinamico dell'array dello STEP 4 (array3d_fourier).
L'array e' f(x, y, s): superficie spaziale (x range, y azimuth) che evolve
mentre si scorre lungo l'asse s (range stirato 0..1000 m).
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
ap.add_argument("--out", default="s1d/step4_array3d_dinamico.html")
ap.add_argument("--nframes", type=int, default=80)
a = ap.parse_args()

d = np.load(a.npz, allow_pickle=True)
F = d["fourier"]          # (ny, nx, Ns)
s = d["s"]; x = d["x"]; y = d["y"]
ny, nx, ns = F.shape
zmin, zmax = float(F.min()), float(F.max())

idx = np.linspace(0, ns - 1, min(a.nframes, ns)).astype(int)
Xg, Yg = np.meshgrid(x, y)     # (ny, nx)

def surf(t):
    return go.Surface(x=Xg, y=Yg, z=F[:, :, t], cmin=zmin, cmax=zmax,
                      colorscale="Turbo", colorbar=dict(title="f [mm]"))

frames = [go.Frame(data=[surf(t)], name=str(t),
                   layout=go.Layout(title=f"Step 4 — array 3D di Fourier  |  s = {s[t]:.0f} m"))
          for t in idx]
fig = go.Figure(data=[surf(idx[0])], frames=frames)
fig.update_layout(
    title=f"Step 4 — array 3D di Fourier (dinamico su s 0..1000 m)  |  s = {s[idx[0]]:.0f} m",
    scene=dict(xaxis_title="x range [m]", yaxis_title="y azimuth [m]",
               zaxis_title="f(x,y,s) [mm]", zaxis=dict(range=[zmin, zmax]),
               aspectratio=dict(x=1, y=1.1, z=0.7)),
    updatemenus=[dict(type="buttons", showactive=False, x=0.05, y=0.05,
                      buttons=[dict(label="▶ Play", method="animate",
                                    args=[None, dict(frame=dict(duration=80, redraw=True),
                                                     fromcurrent=True)]),
                               dict(label="⏸ Pause", method="animate",
                                    args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                       mode="immediate")])])],
    sliders=[dict(currentvalue=dict(prefix="s = ", suffix=" m"),
                  steps=[dict(method="animate", label=f"{s[t]:.0f}",
                              args=[[str(t)], dict(mode="immediate",
                                                   frame=dict(duration=0, redraw=True))])
                         for t in idx])])
fig.write_html(a.out, include_plotlyjs="cdn")
print(f"salvato {a.out}  ({len(idx)} frame, surface {ny}x{nx})")

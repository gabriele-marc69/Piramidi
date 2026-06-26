#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
forma_donda_rettangolo.py
Tutti i grafici hanno il footprint REALE del rettangolo geografico (m, dai corner).
Le 6 armoniche (una per ognuno dei 6 strati .tiff) vengono SOVRAPPOSTE per
generare la forma d'onda contenuta nel rettangolo:

    W(x,y) = somma_{k=1..6}  D_k(x,y) * sin(2*pi*k*x / Lx)

D_k = strato di spostamento LOS della scena k (riferito alla media temporale,
cosi' tutti e 6 gli strati sono non nulli => 6 armoniche reali).
Output (cartella rettangolo/): superficie 3D interattiva, animazione della
sovrapposizione armonica, mappa di densita' (variazione di frequenza), .npz.
"""
import sys, os, re
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

C = 299_792_458.0
OUT = "rettangolo"
os.makedirs(OUT, exist_ok=True)
NXX = 400          # campioni fini lungo x (E-W) per disegnare le onde

# ---- dati e rettangolo geografico ------------------------------------------
d = np.load("box.npz", allow_pickle=True)
vv = d["vv"]                                   # (nS, ny, nx) complex
nomi = [str(x) for x in d["vv_scene"]]
inc = float(d["incid_mid"]); f_em = float(d["f_em"])
lam = C / f_em; K = lam / (4 * np.pi) * 1000.0     # mm/rad
c = d["corners_lonlat"]; lonW, latN = c[0]; lonE, latS = c[2]
latm = np.radians((latN + latS) / 2)
Lx = (lonE - lonW) * 111320 * np.cos(latm)     # E-W [m]
Ly = (latN - latS) * 110540                    # N-S [m]
nS, ny, nx = vv.shape

# ordina per data
date = [re.search(r"-vv-(\d{8})t", n).group(1) for n in nomi]
order = np.argsort(date); vv = vv[order]; date = [date[i] for i in order]

# 6 strati di spostamento riferiti alla MEDIA temporale (tutti non nulli)
mean = vv.mean(0)
phi = np.angle(vv * np.conj(mean)[None])       # (nS, ny, nx)
D = (-K * phi).astype(np.float32)              # mm, 6 strati = 6 ampiezze armoniche
print(f"Rettangolo Lx={Lx:.1f} m (E-W) x Ly={Ly:.1f} m (N-S); {nS} strati {date}")

# assi reali del rettangolo
x_px = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
xf = np.linspace(0, Lx, NXX)                    # x fine per le onde

# interpola i 6 strati sul x fine, poi sovrapponi le 6 armoniche
Df = np.empty((nS, ny, NXX), np.float32)
for k in range(nS):
    for i in range(ny):
        Df[k, i] = np.interp(xf, x_px, D[k, i])
comp = np.stack([Df[k] * np.sin(2 * np.pi * (k + 1) * xf / Lx)[None, :]
                 for k in range(nS)], 0)        # (nS, ny, NXX) armoniche
W = comp.sum(0).astype(np.float32)             # (ny, NXX) forma d'onda nel rettangolo
cum = np.cumsum(comp, 0)                        # sovrapposizione progressiva
print(f"W {W.shape}  ampiezza [{W.min():.1f}, {W.max():.1f}] mm")

np.savez_compressed(os.path.join(OUT, "forma_donda.npz"),
                    W=W, comp=comp, x=xf, y=y, Lx=Lx, Ly=Ly,
                    D=D, x_px=x_px, dates=np.array(date))

Xg, Yg = np.meshgrid(xf, y)                     # (ny, NXX)
zlim = float(np.abs(W).max())
asp = dict(x=1.0, y=Ly / Lx, z=0.35)

# ---- 1) superficie 3D interattiva, footprint = rettangolo ------------------
fig = go.Figure(go.Surface(x=Xg, y=Yg, z=W, colorscale="Turbo",
                           colorbar=dict(title="W [mm]")))
fig.update_layout(
    title=f"Forma d'onda nel rettangolo {Lx:.0f}×{Ly:.0f} m — 6 armoniche sovrapposte",
    scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
               zaxis_title="W(x,y) [mm]", aspectmode="manual", aspectratio=asp,
               zaxis=dict(range=[-zlim, zlim])))
fig.write_html(os.path.join(OUT, "forma_donda_3d.html"), include_plotlyjs="cdn")

# statico png
figm = plt.figure(figsize=(10, 7)); ax = figm.add_subplot(111, projection="3d")
ax.plot_surface(Xg, Yg, W, cmap="turbo", linewidth=0, antialiased=True)
ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("W [mm]")
ax.set_box_aspect((Lx, Ly, 0.4 * max(Lx, Ly)))
ax.set_title(f"Forma d'onda nel rettangolo {Lx:.0f}×{Ly:.0f} m (6 armoniche)")
figm.savefig(os.path.join(OUT, "forma_donda_3d.png"), dpi=130, bbox_inches="tight")
plt.close(figm)

# ---- 2) animazione: sovrapposizione progressiva delle 6 armoniche ----------
frames = [go.Frame(data=[go.Surface(x=Xg, y=Yg, z=cum[k], colorscale="Turbo",
                                    cmin=-zlim, cmax=zlim)],
                   name=str(k + 1),
                   layout=go.Layout(title=f"Sovrapposizione armoniche 1..{k+1} di 6"))
          for k in range(nS)]
figd = go.Figure(data=[go.Surface(x=Xg, y=Yg, z=cum[0], colorscale="Turbo",
                                  cmin=-zlim, cmax=zlim, colorbar=dict(title="W [mm]"))],
                 frames=frames)
figd.update_layout(
    title="Costruzione della forma d'onda: armonica 1 di 6",
    scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
               zaxis_title="W [mm]", aspectmode="manual", aspectratio=asp,
               zaxis=dict(range=[-zlim, zlim])),
    updatemenus=[dict(type="buttons", x=0.05, y=0.05, showactive=False,
                      buttons=[dict(label="▶ Play", method="animate",
                                    args=[None, dict(frame=dict(duration=600, redraw=True),
                                                     fromcurrent=True)]),
                               dict(label="⏸", method="animate",
                                    args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                       mode="immediate")])])],
    sliders=[dict(currentvalue=dict(prefix="armoniche 1.."),
                  steps=[dict(method="animate", label=str(k + 1),
                              args=[[str(k + 1)], dict(mode="immediate",
                                                       frame=dict(duration=0, redraw=True))])
                         for k in range(nS)])])
figd.write_html(os.path.join(OUT, "forma_donda_dinamico.html"), include_plotlyjs="cdn")

# ---- 3) densita': punti di variazione di frequenza lungo x -----------------
ds = xf[1] - xf[0]
an = hilbert(W, axis=1)
inst_f = np.gradient(np.unwrap(np.angle(an), axis=1), ds, axis=1) / (2 * np.pi)
dfreq = np.abs(np.gradient(inst_f, ds, axis=1))
thr = np.percentile(dfreq, 95)
yy, xx = np.where(dfreq >= thr)
fig3 = go.Figure(go.Scatter3d(
    x=xf[xx], y=y[yy], z=W[yy, xx], mode="markers",
    marker=dict(size=3, color=dfreq[yy, xx], colorscale="Inferno",
                colorbar=dict(title="|Δfreq| densità"))))
fig3.update_layout(
    title="Variazione di frequenza (=> densità) sulla forma d'onda del rettangolo",
    scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]", zaxis_title="W [mm]",
               aspectmode="manual", aspectratio=asp))
fig3.write_html(os.path.join(OUT, "forma_donda_densita.html"), include_plotlyjs="cdn")

fig4, ax = plt.subplots(figsize=(9, 5))
im = ax.imshow(dfreq, cmap="inferno", origin="lower", aspect="auto",
               extent=[0, Lx, 0, Ly])
ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]")
ax.set_title("Variazione di frequenza |Δfreq| (proxy densità) nel rettangolo")
fig4.colorbar(im, ax=ax, label="|Δfreq| [cicli/m²]")
fig4.savefig(os.path.join(OUT, "forma_donda_densita.png"), dpi=130, bbox_inches="tight")
plt.close(fig4)

print("salvati in", OUT + "/:",
      "forma_donda_3d.html/.png, forma_donda_dinamico.html, "
      "forma_donda_densita.html/.png, forma_donda.npz")

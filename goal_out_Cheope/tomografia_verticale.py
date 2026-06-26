#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tomografia_verticale.py — stile Figura 5 (pag.9) di 2208.00811v1.
Base = rettangolo geografico 268x276 m (il "green plane", superficie/medium boundary).
Su una griglia di pixel della base, sviluppa VERSO IL BASSO (profondita') la forma
d'onda data dall'ANTI-TRASFORMATA DI FOURIER delle 6 armoniche: i 6 strati .tiff
sono i coefficienti spettrali k=1..6, irfft -> onda stazionaria lungo la profondita'
(la "molla" oscillante ancorata a ogni pixel).
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
ZTOP, ZBOT = 210.0, -200.0   # finestra quota reale s.l.m. [m]
ZDEPTH = ZTOP - ZBOT         # estensione verticale = 400 m (sviluppo onde)
NZ = 256                # campioni in profondita'
CONF = dict(scrollZoom=True, displayModeBar=True)   # zoom abilitato
STEP_X = 6              # downsampling griglia molle (range)
STEP_Y = 2              # downsampling griglia molle (azimuth)
WIGGLE = 14.0          # ampiezza laterale max delle onde [m] (solo resa grafica)

d = np.load(os.path.join(OUT, "forma_donda.npz"), allow_pickle=True)
D = d["D"]                         # (6, ny, nx) strati di spostamento [mm]
x_px = d["x_px"]; y = d["y"]; Lx = float(d["Lx"]); Ly = float(d["Ly"])
nS, ny, nx = D.shape
# quota reale del terreno (DEM altopiano) alle coordinate reali del box
z0 = np.load(os.path.join(OUT, "dem.npz"))["z0"].astype(float)   # (ny, nx) [m s.l.m.]
# + geometria REALE della piramide di Chefren (quote da ricerca web):
#   base 215.5 m, altezza attuale 136.4 m, base quadrata allineata ai cardinali,
#   apice al centro del box (il box e' centrato su Chefren). Il DEM 90 m non la
#   risolve, quindi la sovrappongo esplicitamente all'altopiano reale.
H_PYR, BASE = 136.4, 215.5
half = BASE / 2.0
xc, yc = Lx / 2.0, Ly / 2.0
XX, YY = np.meshgrid(x_px, y)
cheb = np.maximum(np.abs(XX - xc), np.abs(YY - yc))      # base quadrata (Chebyshev)
pyr = np.clip(H_PYR * (1.0 - cheb / half), 0.0, None)    # profilo piramide [m]
z0 = z0 + pyr
print(f"Quota reale: altopiano {z0.min():.0f} m + piramide Chefren -> apice {z0.max():.0f} m s.l.m.")

# fase interferometrica reale per scena (per dare una fase distinta a ogni armonica
# -> evita l'antisimmetria a meta' profondita' che raddoppiava i valori)
box = np.load("box.npz", allow_pickle=True)
import re as _re
vv = box["vv"]; nomi = [str(x) for x in box["vv_scene"]]
_date = [_re.search(r"-vv-(\d{8})t", n).group(1) for n in nomi]
vv = vv[np.argsort(_date)]
phi = np.angle(vv * np.conj(vv.mean(0))[None])     # (nS, ny, nx) fase per scena
# coefficiente armonico complesso: ampiezza = |spostamento| [mm], fase = fase interfer.
coeff = (np.abs(D) * np.exp(1j * phi)).astype(complex)   # (nS, ny, nx)

# asse profondita' verso il basso (z negativo)
z = np.linspace(0, ZDEPTH, NZ)
nf = NZ // 2 + 1

# ---- ANTI-TRASFORMATA: i 6 strati = coefficienti armonici complessi k=1..6 -
# per ogni pixel costruisci lo spettro (solo bin 1..6) e fai irfft -> w(z)
Cpix = coeff.reshape(nS, ny * nx)                   # (6, npix) complesso
spec = np.zeros((ny * nx, nf), complex)
spec[:, 1:1 + nS] = Cpix.T                          # ampiezza+fase reali per armonica
w = np.fft.irfft(spec, n=NZ, axis=1) * NZ          # (npix, NZ) forma d'onda in profondita'
w = w.reshape(ny, nx, NZ).astype(np.float32)
wmax = float(np.abs(w).max())
print(f"base {Lx:.0f}x{Ly:.0f} m, profondita' 0..{ZDEPTH:.0f} m; "
      f"anti-FT 6 armoniche -> w {w.shape}, |w|max={wmax:.1f}")

# ---- griglia di "molle" sviluppate verso il basso --------------------------
ii = np.arange(0, ny, STEP_Y); jj = np.arange(0, nx, STEP_X)
sc = WIGGLE / wmax
X, Y, Z, Cc = [], [], [], []
for i in ii:
    for j in jj:
        wp = w[i, j]; gap = np.array([np.nan])
        X.append(x_px[j] + sc * wp); X.append(gap)
        Y.append(np.full(NZ, y[i])); Y.append(gap)
        Z.append(z0[i, j] - z);     Z.append(gap)     # parte dalla quota reale
        Cc.append(np.abs(wp));      Cc.append(gap)
X = np.concatenate(X); Y = np.concatenate(Y); Z = np.concatenate(Z); Cc = np.concatenate(Cc)

# piano di superficie alla QUOTA REALE (DEM), colorato per quota [m s.l.m.]
Xs, Ys = np.meshgrid(x_px, y)
plane = go.Surface(x=Xs, y=Ys, z=z0, surfacecolor=z0,
                   colorscale="Earth", opacity=0.7,
                   colorbar=dict(title="quota s.l.m. [m]", x=-0.08),
                   name="superficie reale (DEM)")
springs = go.Scatter3d(x=X, y=Y, z=Z, mode="lines",
                       line=dict(width=3, color=Cc, colorscale="Turbo",
                                 colorbar=dict(title="|w| profondità")),
                       connectgaps=False, name="onde (anti-FT 6 armoniche)")

# punti dell'altezza del suolo (DEM) in rosso
ground = go.Scatter3d(x=Xs.ravel(), y=Ys.ravel(), z=z0.ravel(), mode="markers",
                      marker=dict(size=2.5, color="red"), name="altezza suolo (DEM)")
fig = go.Figure([plane, springs, ground])
fig.update_layout(
    title="Tomografia verticale stile Fig.5 — onde dalla quota reale "
          f"(altopiano + piramide Chefren, {z0.min():.0f}–{z0.max():.0f} m s.l.m.)",
    scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
               zaxis_title="quota reale s.l.m. [m]", zaxis=dict(range=[ZBOT, ZTOP]),
               aspectmode="manual",
               aspectratio=dict(x=1.0, y=Ly / Lx, z=ZDEPTH / Lx)))
fig.write_html(os.path.join(OUT, "tomografia_verticale.html"), include_plotlyjs="cdn",
               config=CONF)

# B-scan: sezione verticale (range x profondita') su una linea azimuth centrale
i0 = ny // 2
B = w[i0].T                      # (NZ, nx)
figb, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(B, cmap="seismic", aspect="auto", origin="upper",
               extent=[0, Lx, ZDEPTH, 0], vmin=-wmax, vmax=wmax)
ax.set_xlabel("x E-W [m]"); ax.set_ylabel("profondità [m]")
ax.set_title(f"B-scan: anti-FT delle 6 armoniche sviluppata in profondità "
             f"(linea y={y[i0]:.0f} m)")
figb.colorbar(im, ax=ax, label="ampiezza onda w")
figb.savefig(os.path.join(OUT, "tomografia_bscan.png"), dpi=130, bbox_inches="tight")
plt.close(figb)

# anteprima statica 3D delle molle
figm = plt.figure(figsize=(9, 8)); ax = figm.add_subplot(111, projection="3d")
ax.plot_surface(Xs, Ys, z0, cmap="terrain", alpha=0.4, linewidth=0)
ax.scatter(Xs.ravel(), Ys.ravel(), z0.ravel(), c="red", s=4, depthshade=False)
for i in ii:
    for j in jj:
        ax.plot(x_px[j] + sc * w[i, j], np.full(NZ, y[i]), z0[i, j] - z, lw=0.6)
ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZDEPTH))
ax.set_title("Onde anti-FT dalla quota reale (finestra +200/−200 m s.l.m.)")
figm.savefig(os.path.join(OUT, "tomografia_verticale.png"), dpi=130, bbox_inches="tight")
plt.close(figm)

np.savez_compressed(os.path.join(OUT, "tomografia_verticale.npz"),
                    w=w, z=z, x=x_px, y=y, Lx=Lx, Ly=Ly, zdepth=ZDEPTH, z0=z0)
print("salvati in", OUT + "/: tomografia_verticale.html/.png, tomografia_bscan.png, .npz")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
goal_pipeline.py — esegue gli step 3..6 del goal.txt a partire da box.npz
(il chip-stack complesso ritagliato dallo step 2, skill sentinel1-slc-reader).

 3) array 3D degli spostamenti (DInSAR LOS): (n_strati, ny, nx)  [mm]
 4) array 3D di funzioni di Fourier: per ogni punto (x,y) i valori dei
    n_strati sono usati come ampiezze delle armoniche di una funzione di
    Fourier "stirata" su un range di 1000 m  ->  (ny, nx, Ns)
 5) grafico 3D dinamico (plotly) delle funzioni di Fourier
 6) grafico 3D dei punti dove varia la frequenza (=> densita') delle funzioni,
    via frequenza istantanea (trasformata di Hilbert).
"""
import sys, re, os, argparse
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
NPZ = "box.npz"
RANGE_M = 1000.0     # "range di 1000 metri" dello step 4
NS = 400             # campioni della funzione di Fourier sul range di 1000 m
OUTDIR = "."         # impostato in main()
SAT = ""             # filtro sensore (es. "s1d"); vuoto = tutti


def P(name):         # path di output nella cartella scelta
    return os.path.join(OUTDIR, name)


# ------------------------------------------------------------------ step 3
def step3_spostamenti():
    d = np.load(NPZ, allow_pickle=True)
    vv = d["vv"]                                  # (nS, ny, nx) complex
    nomi = [str(x) for x in d["vv_scene"]]
    inc = float(d["incid_mid"])
    f_em = float(d["f_em"]) if "f_em" in d.files else 5.405e9
    dR = float(d["dR"]); dA = float(d["dA"])
    lam = C / f_em
    K = lam / (4 * np.pi) * 1000.0                # mm per radiante
    gr = dR / np.sin(np.radians(inc))             # ground-range pixel [m]

    if SAT:                                       # tieni un solo sensore
        keep = [i for i, n in enumerate(nomi) if SAT.lower() in n.lower()]
        vv = vv[keep]; nomi = [nomi[i] for i in keep]
        print(f"    filtro sensore '{SAT}': {len(keep)} scene")

    # ordina per data, riferimento = scena piu' antica
    date = [re.search(r"-vv-(\d{8})t", n).group(1) for n in nomi]
    d0 = np.datetime64(f"{date[0][:4]}-{date[0][4:6]}-{date[0][6:8]}")
    giorni = np.array([(np.datetime64(f"{x[:4]}-{x[4:6]}-{x[6:8]}") - d0) /
                       np.timedelta64(1, "D") for x in date], float)
    order = np.argsort(giorni)
    vv = vv[order]; date = [date[i] for i in order]; giorni = giorni[order]

    ref = vv[0]
    phi = np.angle(vv * np.conj(ref)[None])       # (nS, ny, nx) wrapped
    disp = (-K * phi).astype(np.float32)          # mm  -> strati di spostamento
    coh = np.abs(np.exp(1j * phi).mean(0))        # coerenza temporale per pixel

    nS, ny, nx = disp.shape
    x = np.arange(nx) * gr                         # asse x [m] (range a terra)
    y = np.arange(ny) * dA                         # asse y [m] (azimuth)

    np.savez_compressed(P("array3d_spostamenti.npz"),
                        disp=disp, x=x, y=y, dates=np.array(date),
                        giorni=giorni, coh=coh.astype(np.float32),
                        gr=gr, dA=dA, lam=lam)
    print(f"[3] array spostamenti {disp.shape} (strati, y, x)  "
          f"disp [{disp.min():.1f},{disp.max():.1f}] mm  coh medio {coh.mean():.2f}")
    print(f"    estensione box: x 0..{x[-1]:.0f} m, y 0..{y[-1]:.0f} m")

    # figura riassuntiva degli strati
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for k, ax in enumerate(axes.ravel()):
        if k < nS:
            im = ax.imshow(disp[k], cmap="RdBu_r", vmin=-15, vmax=15,
                           aspect="auto", origin="lower")
            ax.set_title(f"strato {k}  ({date[k]})", fontsize=9)
            ax.set_xlabel("x (range) px"); ax.set_ylabel("y (az) px")
        else:
            ax.axis("off")
    fig.colorbar(im, ax=axes, shrink=0.7, label="spostamento LOS [mm]")
    fig.suptitle("Step 3 — strati 3D degli spostamenti (DInSAR)")
    fig.savefig(P("step3_strati_spostamenti.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    return disp, x, y, coh


# ------------------------------------------------------------------ step 4
def step4_fourier(disp, x, y):
    nS, ny, nx = disp.shape
    s = np.linspace(0.0, RANGE_M, NS)             # dominio "stirato" 0..1000 m
    # base armonica: k-esimo strato -> armonica k con periodo fondamentale 1000 m
    k = np.arange(1, nS + 1)
    base = np.sin(2 * np.pi * np.outer(k, s) / RANGE_M)   # (nS, NS)

    # f(x,y; s) = somma_k  disp_k(x,y) * sin(2 pi k s / 1000)
    fourier = np.tensordot(disp.transpose(1, 2, 0), base, axes=([2], [0]))
    fourier = fourier.astype(np.float32)          # (ny, nx, NS)

    np.savez_compressed(P("array3d_fourier.npz"),
                        fourier=fourier, s=s, x=x, y=y)
    print(f"[4] array Fourier {fourier.shape} (y, x, s)  "
          f"s 0..{RANGE_M:.0f} m in {NS} campioni, {nS} armoniche")
    return fourier, s


# ------------------------------------------------------------------ step 5
def step5_grafico_dinamico(fourier, s, x, y):
    ny, nx, ns = fourier.shape
    smin, smax = float(fourier.min()), float(fourier.max())
    # una superficie per ogni colonna range j: Z su (azimuth y) x (s),
    # animata mentre j scorre su tutto il box  -> grafico 3D dinamico.
    Sg, Yg = np.meshgrid(s, y)                     # (ny, ns)

    def surf(j):
        return go.Surface(x=Sg, y=Yg, z=fourier[:, j, :], cmin=smin, cmax=smax,
                          colorscale="Viridis", colorbar=dict(title="f(s) [mm]"))

    j0 = nx // 2
    frames = [go.Frame(data=[surf(j)], name=str(j),
                       layout=go.Layout(title=f"Funzioni di Fourier — colonna range j={j} "
                                              f"(x={x[j]:.0f} m)"))
              for j in range(nx)]
    fig = go.Figure(data=[surf(j0)], frames=frames)
    fig.update_layout(
        title=f"Step 5 — funzioni di Fourier 3D (dinamico) — colonna range j={j0}",
        scene=dict(xaxis_title="s — range stirato [m] (0..1000)",
                   yaxis_title="y azimuth [m]",
                   zaxis_title="f(s) [mm]",
                   zaxis=dict(range=[smin, smax])),
        updatemenus=[dict(type="buttons", showactive=False, x=0.05, y=0.05,
                          buttons=[dict(label="▶ Play", method="animate",
                                        args=[None, dict(frame=dict(duration=120, redraw=True),
                                                         fromcurrent=True)]),
                                   dict(label="⏸ Pause", method="animate",
                                        args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                           mode="immediate")])])],
        sliders=[dict(active=j0, currentvalue=dict(prefix="colonna range j = "),
                      steps=[dict(method="animate", label=str(j),
                                  args=[[str(j)], dict(mode="immediate",
                                                       frame=dict(duration=0, redraw=True))])
                             for j in range(nx)])])
    fig.write_html(P("step5_fourier_3d_dinamico.html"), include_plotlyjs="cdn")
    print("[5] salvato step5_fourier_3d_dinamico.html (interattivo + animato)")


# ------------------------------------------------------------------ step 6
def step6_variazioni_frequenza(fourier, s, x, y, coh):
    ny, nx, ns = fourier.shape
    ds = s[1] - s[0]
    # frequenza istantanea per ogni funzione di Fourier (Hilbert)
    an = hilbert(fourier, axis=2)
    inst_phase = np.unwrap(np.angle(an), axis=2)
    inst_freq = np.gradient(inst_phase, ds, axis=2) / (2 * np.pi)   # cicli/m
    # variazione di frequenza = |d(freq)/ds|  -> proxy di variazione di densita'
    dfreq = np.abs(np.gradient(inst_freq, ds, axis=2)).astype(np.float32)

    # soglia: punti con variazione marcata (percentile alto)
    thr = np.percentile(dfreq, 96.0)
    yy, xx, ss = np.where(dfreq >= thr)
    val = dfreq[yy, xx, ss]
    # mappa su coordinate fisiche: x range [m], y azimuth [m], s [m]
    px = x[xx]; py = y[yy]; ps = s[ss]
    print(f"[6] soglia variazione freq = {thr:.3e} cicli/m^2 ; "
          f"{len(val)} punti marcati su {dfreq.size}")

    np.savez_compressed(P("variazioni_frequenza.npz"),
                        inst_freq=inst_freq.astype(np.float32), dfreq=dfreq,
                        thr=thr, px=px, py=py, ps=ps, val=val)

    fig = go.Figure(data=[go.Scatter3d(
        x=ps, y=py, z=px, mode="markers",
        marker=dict(size=3, color=val, colorscale="Inferno", opacity=0.8,
                    colorbar=dict(title="|Δfreq| (densità)")))])
    fig.update_layout(
        title="Step 6 — punti di variazione di frequenza (=> densità) delle funzioni di Fourier",
        scene=dict(xaxis_title="s — range stirato [m]",
                   yaxis_title="y azimuth [m]",
                   zaxis_title="x range [m]"))
    fig.write_html(P("step6_variazioni_frequenza.html"), include_plotlyjs="cdn")

    # PNG statico riassuntivo: mappa della variazione media di freq per pixel
    var_map = dfreq.mean(axis=2)
    fig2, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(var_map, cmap="inferno", aspect="auto", origin="lower",
                   extent=[x[0], x[-1], y[0], y[-1]])
    ax.set_xlabel("x range [m]"); ax.set_ylabel("y azimuth [m]")
    ax.set_title("Step 6 — variazione media di frequenza per punto (proxy densità)")
    fig2.colorbar(im, ax=ax, label="media |Δfreq| [cicli/m²]")
    fig2.savefig(P("step6_mappa_densita.png"), dpi=130, bbox_inches="tight")
    plt.close(fig2)
    print("[6] salvato step6_variazioni_frequenza.html + step6_mappa_densita.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sat", default="", help="tieni un solo sensore, es. s1d")
    ap.add_argument("--outdir", default=".", help="cartella di output")
    a = ap.parse_args()
    SAT = a.sat
    OUTDIR = a.outdir
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"== pipeline (sat='{SAT or 'tutti'}', outdir='{OUTDIR}') ==")
    disp, x, y, coh = step3_spostamenti()
    fourier, s = step4_fourier(disp, x, y)
    step5_grafico_dinamico(fourier, s, x, y)
    step6_variazioni_frequenza(fourier, s, x, y, coh)
    print("\nPipeline goal.txt completata.")

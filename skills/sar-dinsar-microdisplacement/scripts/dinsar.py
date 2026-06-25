#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dinsar.py — per-pixel LOS micro-displacement (mm) from a Sentinel-1 SLC chip-stack .npz
via differential interferometry (DInSAR). Reference = earliest scene.

  phi = angle(z_s · conj(z_ref));  d = -(lambda/4pi)·phi  [mm];  fringe 2pi = lambda/2.

Outputs: per-pixel/date CSV, per-pixel velocity+coherence CSV, waveform PNG, 3-D cloud PNG.
The .npz is the one written by the sentinel1-slc-reader skill (keys: vv, dR, dA,
incid_mid, vv_scene). Use --sat to keep a single satellite (cleaner phase).

  python dinsar.py --npz box.npz --sat s1d --out disp
"""
import sys, re, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default="box_chefren_slc.npz")
    ap.add_argument("--sat", default="", help="keep only this sensor, e.g. s1d (default: all)")
    ap.add_argument("--out", default="microdisplacement")
    a = ap.parse_args()

    d = np.load(a.npz, allow_pickle=True)
    vv = d["vv"]                                   # (n_scene, az, rng) complex
    nomi = [str(x) for x in d["vv_scene"]]
    inc = float(d["incid_mid"])
    f_em = float(d["f_em"]) if "f_em" in d.files else 5.405e9
    lam = 299_792_458.0 / f_em
    K = lam / (4 * np.pi) * 1000.0                 # mm per radian

    if a.sat:
        keep = [i for i, n in enumerate(nomi) if a.sat.lower() in n.lower()]
        vv = vv[keep]; nomi = [nomi[i] for i in keep]
    nS, H, W = vv.shape

    def data_di(nm):
        return re.search(r"-vv-(\d{8})t", nm).group(1)
    date = [data_di(n) for n in nomi]
    d0 = np.datetime64(f"{date[0][:4]}-{date[0][4:6]}-{date[0][6:8]}")
    giorni = np.array([(np.datetime64(f"{x[:4]}-{x[4:6]}-{x[6:8]}") - d0) /
                       np.timedelta64(1, "D") for x in date], float)
    order = np.argsort(giorni)
    vv, nomi, date, giorni = vv[order], [nomi[i] for i in order], \
        [date[i] for i in order], giorni[order]
    print(f"{nS} scenes (ref={date[0]}):", date)

    ref = vv[0]
    ifg = vv * np.conj(ref)[None]                  # (nS, H, W)
    phi = np.angle(ifg)
    disp = -K * phi                                # mm, wrapped
    # per-pixel temporal phase coherence: |mean_s exp(j phi_s)| (1=stable, 0=noise)
    coh = np.abs(np.exp(1j * phi).mean(0))

    # linear velocity per pixel
    G = np.vstack([giorni, np.ones_like(giorni)]).T
    flat = disp.reshape(nS, -1)
    sol, *_ = np.linalg.lstsq(G, flat, rcond=None)
    vel = sol[0].reshape(H, W) * 365.0             # mm/yr
    resid = flat - G @ sol
    rmse = np.sqrt((resid ** 2).mean(0)).reshape(H, W)

    np.savetxt(f"{a.out}_pixel.csv",
               np.column_stack([np.repeat(np.arange(H), W), np.tile(np.arange(W), H),
                                vel.ravel(), coh.ravel(), rmse.ravel()]),
               delimiter=",", header="az,rng,vel_mm_yr,coherence,rmse_mm", comments="")
    print(f"disp/date: min={disp.min():.2f} max={disp.max():.2f} std={disp.std():.2f} mm | "
          f"mean coherence={coh.mean():.2f}")

    # waveforms (coherent pixels highlighted)
    fig, ax = plt.subplots(figsize=(9, 5))
    mask = coh.ravel() > 0.6
    for p in range(H * W):
        ax.plot(giorni, flat[:, p], color=("crimson" if mask[p] else "0.8"),
                lw=(1.0 if mask[p] else 0.3), alpha=(0.8 if mask[p] else 0.3))
    ax.set_xlabel("days from reference"); ax.set_ylabel("LOS displacement [mm]")
    ax.set_title(f"Micro-displacement waveforms — {nS} scenes "
                 f"({mask.sum()} coherent>0.6, fringe={lam/2*1000:.1f} mm)")
    fig.tight_layout(); fig.savefig(f"{a.out}_waveforms.png", dpi=140); plt.close(fig)

    # 3-D cloud: cumulative displacement of the last scene
    from mpl_toolkits.mplot3d import Axes3D  # noqa
    fig = plt.figure(figsize=(9, 7)); ax = fig.add_subplot(111, projection="3d")
    yy, xx = np.mgrid[0:H, 0:W]
    cum = disp[-1]
    sc = ax.scatter(xx.ravel(), yy.ravel(), cum.ravel(), c=coh.ravel(),
                    cmap="viridis", s=10)
    ax.set_xlabel("range px"); ax.set_ylabel("azimuth px"); ax.set_zlabel("cum. disp [mm]")
    ax.set_title(f"Cumulative LOS displacement ({date[0]}→{date[-1]})")
    fig.colorbar(sc, ax=ax, shrink=0.6, label="coherence")
    fig.tight_layout(); fig.savefig(f"{a.out}_3d.png", dpi=140); plt.close(fig)
    print(f"Saved: {a.out}_pixel.csv, {a.out}_waveforms.png, {a.out}_3d.png")


if __name__ == "__main__":
    main()

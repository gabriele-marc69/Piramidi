#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tomographic_images.py — build the tomographic reflectivity volume V[az,rng,z] from a
Sentinel-1 SLC stack over a ground box and render the publication-style figures
(Biondi & Malanga, RS 2022, 14, 5231):

  sezione_range_quota.png    vertical B-scan reflectivity vs (ground range, depth)
  sezione_azimuth_quota.png  vertical B-scan reflectivity vs (azimuth, depth)
  slice_orizzontali.png      horizontal slices of V at increasing depth
  mappa_quota.png            dominant-scatterer height (argmax|h|) + peak intensity
  volume_tomografico.npy     full V for 3-D scatterer clouds

Usage:
  python tomographic_images.py --stack stack_slc --outdir tomo_img \
      --nw 29 58 38.0 N 31 7 45.4 E --se 29 58 29.0 N 31 7 55.4 E \
      --zmax 8.5 --klook 12 --ovs 4
Defaults target the Khafre (Chefren) pyramid box. Stay within z_amb (printed) or depths alias.
"""
import os, sys, glob, re, argparse
import numpy as np
import rasterio
from rasterio.windows import Window
from rasterio.transform import GCPTransformer
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sar_tomo as T
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def dms2dec(d, m, s, h):
    v = float(d) + float(m) / 60 + float(s) / 3600
    return -v if str(h).upper() in ("S", "W") else v


def annotation(tif, stack):
    m = re.search(r"(s1[abcd])-iw(\d)-slc-vv-(\d{8})", os.path.basename(tif).lower())
    s, sw, dt = m.groups()
    for x in glob.glob(os.path.join(stack, "*.annotation.xml")):
        xb = os.path.basename(x).lower()
        if s in xb and f"iw{sw}" in xb and dt in xb:
            return x


def leggi_complesso(ds, win):
    a = ds.read(1, window=win)
    if np.iscomplexobj(a):
        return a.astype(np.complex64)
    if a.dtype.names and len(a.dtype.names) == 2:
        n0, n1 = a.dtype.names
        return a[n0].astype(np.float32) + 1j * a[n1].astype(np.float32)
    return a.astype(np.complex64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stack", default="stack_slc")
    ap.add_argument("--outdir", default="tomo_img")
    ap.add_argument("--nw", nargs=6, default=["29", "58", "38.0", "N", "31", "7"],
                    help="DMS NW: deg min sec hemiLat deg min ... (see --se)")
    ap.add_argument("--nw2", nargs=2, default=["45.4", "E"])
    ap.add_argument("--se", nargs=6, default=["29", "58", "29.0", "N", "31", "7"])
    ap.add_argument("--se2", nargs=2, default=["55.4", "E"])
    ap.add_argument("--zmax", type=float, default=8.5)
    ap.add_argument("--klook", type=int, default=12)
    ap.add_argument("--ovs", type=float, default=4.0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    lat_nw = dms2dec(a.nw[0], a.nw[1], a.nw[2], a.nw[3])
    lon_nw = dms2dec(a.nw[4], a.nw[5], a.nw2[0], a.nw2[1])
    lat_se = dms2dec(a.se[0], a.se[1], a.se[2], a.se[3])
    lon_se = dms2dec(a.se[4], a.se[5], a.se2[0], a.se2[1])
    latN, latS = max(lat_nw, lat_se), min(lat_nw, lat_se)
    lonW, lonE = min(lon_nw, lon_se), max(lon_nw, lon_se)
    corners = [(lonW, latN), (lonE, latN), (lonE, latS), (lonW, latS)]

    tiffs = sorted(glob.glob(os.path.join(a.stack, "*-iw[1-3]-slc-vv-*.tiff"))) or \
            sorted(glob.glob(os.path.join(a.stack, "*slc-vv-*.tiff")))
    if not tiffs:
        sys.exit(f"No VV SLC tiff in {a.stack}")

    chips, ref = [], None
    for t in tiffs:
        with rasterio.open(t) as ds:
            if not ds.gcps[0]:
                continue
            with GCPTransformer(ds.gcps[0]) as tf:
                rc = [tf.rowcol(lon, lat) for lon, lat in corners]
            rows = [r for r, c in rc]; cols = [c for r, c in rc]
            if not all(0 <= r < ds.height and 0 <= c < ds.width for r, c in rc):
                continue
            r0 = max(0, int(np.floor(min(rows)))); c0 = max(0, int(np.floor(min(cols))))
            r1 = min(ds.height, int(np.ceil(max(rows))) + 1)
            c1 = min(ds.width, int(np.ceil(max(cols))) + 1)
            chips.append(leggi_complesso(ds, Window(c0, r0, c1 - c0, r1 - r0)))
        if ref is None:
            ref = {"r0": r0, "c0": c0, "tif": t}
    if not chips:
        sys.exit("No scene covers the box.")
    H = min(c.shape[0] for c in chips); W = min(c.shape[1] for c in chips)
    print(f"{len(chips)} scene | box {W}(range) x {H}(azimuth) px")

    g = T.leggi_geometria(annotation(ref["tif"], a.stack))
    geo = T.GeometriaSAR()
    geo.f_em = g["f_em"]; geo.V = g["V"]; geo.theta_deg = g["incid_mid"]
    geo.R0 = g["R_near"] + (ref["c0"] + W / 2.0) * g["dR"]; geo.B_doppler = 1.0 / g["dt_az"]
    dz = geo.risoluzione_tomografica() / a.ovs
    z = np.arange(0.0, a.zmax + dz, dz)
    A = T.steering_matrix(geo, z, a.klook)
    z_amb = T.altezza_ambiguita(geo, a.klook)
    print(f"z 0..{a.zmax} m, {len(z)} samples (dz={dz:.2f} m) | z_amb={z_amb:.1f} m")
    if a.zmax > z_amb * 1.01:
        print(f"  WARNING: zmax > z_amb: depths beyond {z_amb:.1f} m alias!")

    V = np.zeros((H, W, len(z))); n = np.zeros((H, W))
    for s, chip in enumerate(chips):
        chip = chip[:H, :W]
        subs = [None if np.abs(chip[:, c]).max() == 0
                else T.sotto_aperture_doppler(chip[:, c], n_sub=a.klook + 1) for c in range(W)]
        for c in range(W):
            if subs[c] is None:
                continue
            for i in range(H):
                Y = np.array([T.pixel_tracking(subs[c][m], subs[c][m + 1], i)
                              for m in range(a.klook)], complex)
                if np.abs(Y).max() == 0:
                    continue
                V[i, c] += np.abs(A.conj().T @ Y); n[i, c] += 1
        print(f"  scene {s + 1}/{len(chips)}")
    V /= np.maximum(n[..., None], 1)
    Vn = V / (V.max() + 1e-12)
    np.save(os.path.join(a.outdir, "volume_tomografico.npy"), V)

    dA = g["dA"]; gr = g["dR"] / np.sin(np.radians(g["incid_mid"]))
    ew = np.arange(W) * gr; ns = np.arange(H) * dA
    db = lambda x: 10 * np.log10(np.maximum(x, 1e-6))

    ir = H // 2
    fig, ax = plt.subplots(figsize=(9, 4))
    im = ax.imshow(db(Vn[ir].T), origin="lower", aspect="auto", cmap="inferno",
                   extent=[0, ew[-1], 0, a.zmax])
    ax.set_xlabel("ground range [m]"); ax.set_ylabel("depth z [m]")
    ax.set_title(f"Vertical tomographic section (range-depth) — azimuth row {ir}")
    fig.colorbar(im, ax=ax, label="reflectivity [dB]"); fig.tight_layout()
    fig.savefig(f"{a.outdir}/sezione_range_quota.png", dpi=150); plt.close(fig)

    ic = W // 2
    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(db(Vn[:, ic].T), origin="lower", aspect="auto", cmap="inferno",
                   extent=[0, ns[-1], 0, a.zmax])
    ax.set_xlabel("azimuth [m]"); ax.set_ylabel("depth z [m]")
    ax.set_title(f"Vertical tomographic section (azimuth-depth) — range col {ic}")
    fig.colorbar(im, ax=ax, label="reflectivity [dB]"); fig.tight_layout()
    fig.savefig(f"{a.outdir}/sezione_azimuth_quota.png", dpi=150); plt.close(fig)

    quote = [0.0, a.zmax * 0.25, a.zmax * 0.5, a.zmax * 0.75]
    fig, axs = plt.subplots(1, len(quote), figsize=(15, 3.6))
    for k, zz in enumerate(quote):
        iz = int(np.argmin(np.abs(z - zz)))
        im = axs[k].imshow(db(Vn[:, :, iz]), origin="upper", aspect="auto", cmap="inferno",
                           extent=[0, ew[-1], ns[-1], 0])
        axs[k].set_title(f"slice z={z[iz]:.1f} m"); axs[k].set_xlabel("range [m]")
        if k == 0:
            axs[k].set_ylabel("azimuth [m]")
    fig.colorbar(im, ax=axs, label="reflectivity [dB]", fraction=0.02)
    fig.suptitle("Horizontal slices of the tomographic volume")
    fig.savefig(f"{a.outdir}/slice_orizzontali.png", dpi=150); plt.close(fig)

    quota = z[np.argmax(V, axis=2)]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    im0 = axs[0].imshow(quota, origin="upper", aspect="auto", cmap="turbo",
                        extent=[0, ew[-1], ns[-1], 0])
    axs[0].set_title("Tomographic height (max |h(z)|)")
    axs[0].set_xlabel("range [m]"); axs[0].set_ylabel("azimuth [m]")
    fig.colorbar(im0, ax=axs[0], label="height [m]")
    im1 = axs[1].imshow(db(Vn.max(2)), origin="upper", aspect="auto", cmap="inferno",
                        extent=[0, ew[-1], ns[-1], 0])
    axs[1].set_title("Peak intensity"); axs[1].set_xlabel("range [m]")
    fig.colorbar(im1, ax=axs[1], label="[dB]")
    fig.tight_layout(); fig.savefig(f"{a.outdir}/mappa_quota.png", dpi=150); plt.close(fig)

    print(f"\nSaved in {a.outdir}/: sezione_range_quota.png, sezione_azimuth_quota.png, "
          "slice_orizzontali.png, mappa_quota.png, volume_tomografico.npy")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_box.py — geolocate a lat/lon box on every Sentinel-1 SLC TIFF of a stack via the
product GCPs and extract a co-registered complex chip stack (azimuth x range).

Pairs each *-iw?-slc-<pol>-YYYYMMDD-*.tiff with its *.annotation.xml, reads geometry,
maps the box corners to pixels with rasterio GCPTransformer, windowed-reads the complex
chip (never the full scene), crops all scenes to the common minimum size, stacks them.

Saves an .npz (complex `vv` array + geometry scalars incl. f_em) and a per-pixel CSV.

  python extract_box.py --stack stack_slc --pol vv --out box.npz --nw 29 58 38.0 N 31 7 45.4 E --se 29 58 29.0 N 31 7 55.4 E
Defaults: Khafre (Chefren) pyramid box, Giza.
"""
import os, sys, glob, re, argparse
import xml.etree.ElementTree as ET
import numpy as np
import rasterio
from rasterio.windows import Window
from rasterio.transform import GCPTransformer
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
C = 299_792_458.0


def dms2dec(d, m, s, h):
    v = float(d) + float(m) / 60 + float(s) / 3600
    return -v if str(h).upper() in ("S", "W") else v


def trova_annotation(tif, stack, pol):
    m = re.search(rf"(s1[abcd])-iw(\d)-slc-{pol}-(\d{{8}})", os.path.basename(tif).lower())
    if not m:
        return None
    s, sw, dt = m.groups()
    for x in glob.glob(os.path.join(stack, "*.annotation.xml")):
        xb = os.path.basename(x).lower()
        if s in xb and f"iw{sw}" in xb and dt in xb:
            return x


def leggi_geometria(ann):
    r = ET.parse(ann).getroot()
    f = lambda p: float(r.find(p).text)
    g = {"f_em": f(".//radarFrequency"), "dt_az": f(".//azimuthTimeInterval"),
         "slantRangeTime": f(".//imageInformation/slantRangeTime"),
         "incid_mid": f(".//incidenceAngleMidSwath"),
         "dR": f(".//rangePixelSpacing"), "dA": f(".//azimuthPixelSpacing")}
    g["R_near"] = C * g["slantRangeTime"] / 2.0
    return g


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
    ap.add_argument("--pol", default="vv")
    ap.add_argument("--out", default="box.npz")
    ap.add_argument("--nw", nargs=8, default=["29", "58", "38.0", "N", "31", "7", "45.4", "E"])
    ap.add_argument("--se", nargs=8, default=["29", "58", "29.0", "N", "31", "7", "55.4", "E"])
    ap.add_argument("--csv", action="store_true", help="also write per-pixel CSV")
    a = ap.parse_args()

    lat_nw = dms2dec(*a.nw[0:3], a.nw[3]); lon_nw = dms2dec(*a.nw[4:7], a.nw[7])
    lat_se = dms2dec(*a.se[0:3], a.se[3]); lon_se = dms2dec(*a.se[4:7], a.se[7])
    latN, latS = max(lat_nw, lat_se), min(lat_nw, lat_se)
    lonW, lonE = min(lon_nw, lon_se), max(lon_nw, lon_se)
    corners = [(lonW, latN), (lonE, latN), (lonE, latS), (lonW, latS)]

    tiffs = sorted(glob.glob(os.path.join(a.stack, f"*-iw[1-3]-slc-{a.pol}-*.tiff"))) or \
            sorted(glob.glob(os.path.join(a.stack, f"*slc-{a.pol}-*.tiff")))
    if not tiffs:
        sys.exit(f"No {a.pol.upper()} SLC tiff in {a.stack}")

    chips, nomi, geom = [], [], None
    for tif in tiffs:
        ann = trova_annotation(tif, a.stack, a.pol)
        if ann is None:
            print(f"  ({os.path.basename(tif)}: no annotation, skip)"); continue
        g = leggi_geometria(ann)
        with rasterio.open(tif) as ds:
            if not ds.gcps[0]:
                print(f"  ({os.path.basename(tif)}: no GCP, skip)"); continue
            with GCPTransformer(ds.gcps[0]) as tf:
                rc = [tf.rowcol(lon, lat) for lon, lat in corners]
            if not all(0 <= r < ds.height and 0 <= c < ds.width for r, c in rc):
                print(f"  ({os.path.basename(tif)}: box outside raster, skip)"); continue
            rows = [r for r, c in rc]; cols = [c for r, c in rc]
            r0 = max(0, int(np.floor(min(rows)))); c0 = max(0, int(np.floor(min(cols))))
            r1 = min(ds.height, int(np.ceil(max(rows))) + 1)
            c1 = min(ds.width, int(np.ceil(max(cols))) + 1)
            chip = leggi_complesso(ds, Window(c0, r0, c1 - c0, r1 - r0))
        validi = float(np.mean(np.abs(chip) > 0) * 100)
        print(f"  >> {os.path.basename(tif)} {chip.shape[1]}x{chip.shape[0]}px "
              f"(range x az) valid {validi:.0f}%")
        chips.append(chip); nomi.append(os.path.basename(tif)); geom = g
    if not chips:
        sys.exit("No scene covers the box.")

    H = min(c.shape[0] for c in chips); W = min(c.shape[1] for c in chips)
    arr = np.stack([c[:H, :W] for c in chips], 0)
    print(f"stacked: {arr.shape} (n_scene, azimuth, range)")
    np.savez_compressed(a.out, vv=arr, vv_scene=np.array(nomi),
                        corners_lonlat=np.array(corners),
                        dR=geom["dR"], dA=geom["dA"], incid_mid=geom["incid_mid"],
                        f_em=geom["f_em"], R_near=geom["R_near"], dt_az=geom["dt_az"])
    print(f"saved {a.out}")

    if a.csv:
        gr = geom["dR"] / np.sin(np.radians(geom["incid_mid"]))
        out = a.out.rsplit(".", 1)[0] + "_pixels.csv"
        with open(out, "w", encoding="utf-8") as f:
            f.write("scene_idx,scene,az_idx,rng_idx,azimuth_m,range_m,real,imag,amp,phase_rad\n")
            for s in range(arr.shape[0]):
                for i in range(H):
                    for j in range(W):
                        z = arr[s, i, j]
                        f.write(f"{s},{nomi[s]},{i},{j},{i*geom['dA']:.2f},{j*gr:.2f},"
                                f"{z.real:.6g},{z.imag:.6g},{abs(z):.6g},{np.angle(z):.4f}\n")
        print(f"saved {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_dem.py — scarica la quota reale del terreno (DEM, Copernicus GLO-90 via
Open-Meteo) alle coordinate reali del box e la salva sulla griglia dei pixel
(ny x nx) -> rettangolo/dem.npz  (z0 [m s.l.m.], lon, lat).
"""
import sys, os, json, time, urllib.request, urllib.parse
import numpy as np
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = "rettangolo"
box = np.load("box.npz", allow_pickle=True)
ny, nx = box["vv"].shape[1:]
c = box["corners_lonlat"]; lonW, latN = c[0]; lonE, latS = c[2]

# griglia fine = pixel reali; griglia grossolana per il download (DEM ~90 m)
j = np.arange(nx); i = np.arange(ny)
lon = lonW + (j / (nx - 1)) * (lonE - lonW)        # E-W
lat = latS + (i / (ny - 1)) * (latN - latS)        # S->N

GC = 10                                            # griglia grossolana GCxGC <=100
clon = np.linspace(lonW, lonE, GC)
clat = np.linspace(latS, latN, GC)
CLON, CLAT = np.meshgrid(clon, clat)
qs = ",".join(f"{v:.6f}" for v in CLAT.ravel())
qo = ",".join(f"{v:.6f}" for v in CLON.ravel())
url = "https://api.open-meteo.com/v1/elevation?latitude=%s&longitude=%s" % (qs, qo)
print(f"download DEM griglia {GC}x{GC} ({GC*GC} punti) ...")
ce = None
for attempt in range(6):
    try:
        with urllib.request.urlopen(url, timeout=40) as r:
            ce = np.array(json.load(r)["elevation"]).reshape(GC, GC)
        break
    except Exception as ex:
        print(f"  retry {attempt} ({ex})"); time.sleep(5 + 5 * attempt)
if ce is None:
    sys.exit("download DEM fallito")

# interpolazione bilineare della griglia grossolana sui pixel reali
from scipy.interpolate import RegularGridInterpolator
itp = RegularGridInterpolator((clat, clon), ce, bounds_error=False, fill_value=None)
LAT, LON = np.meshgrid(lat, lon, indexing="ij")
z0 = itp(np.stack([LAT.ravel(), LON.ravel()], -1)).reshape(ny, nx)
os.makedirs(OUT, exist_ok=True)
np.savez_compressed(os.path.join(OUT, "dem.npz"), z0=z0, lon=lon, lat=lat,
                    LON=LON, LAT=LAT)
print(f"DEM: quota {z0.min():.0f}..{z0.max():.0f} m s.l.m. (media {z0.mean():.0f}); "
      f"salvato {OUT}/dem.npz")

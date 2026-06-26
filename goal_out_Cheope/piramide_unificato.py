#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramide_unificato.py
=====================
UN SOLO programma che esegue l'intera catena del goal, dai 6 step separati
(extract_box + goal_pipeline + forma_donda + tomografia_verticale +
variazioni_onde + fetch_dem) uniti in un'unica pipeline guidata dalle coordinate.

Step (allineati al goal):
  1) cerca i file .tiff (prodotti Sentinel-1 SLC) che CONTENGONO le coordinate
     passate al programma -> almeno 12 (catalogo Copernicus Data Space, OData).
  2) scarica i .tiff dei prodotti trovati (estrae dal pacchetto SAFE il/i
     measurement TIFF della polarizzazione scelta + annotation.xml).
  3) estrae l'AREA delle coordinate da ogni .tiff via i GCP del prodotto ->
     stack complesso co-registrato (n_scene, azimuth, range) -> box.npz.
  4) estrae la COMPONENTE VERTICALE dello spostamento dei singoli pixel
     (DInSAR LOS -> verticale) per le aree estratte -> array a 12 strati
     (n_strati, ny, nx).
  5) genera il grafico 3D dell'array: i 12 numeri di ogni pixel sono i
     coefficienti armonici di una forma d'onda sviluppata in profondita'
     (anti-trasformata di Fourier) -> superficie reale + "molle" 3D.
  6) estrae i punti dove le forme d'onda variano in frequenza (=> densita',
     Biondi-Malanga) e li disegna in 3D tenendo conto delle ALTEZZE REALI dei
     punti rispetto al livello 0 (DEM Copernicus + geometria del bersaglio).

Esempio:
  python piramide_unificato.py --nw 29 58 38.0 N 31 7 45.4 E \  
                               --se 29 58 29.0 N 31 7 55.4 E \
                               --download            # richiede credenziali CDSE
  python piramide_unificato.py --steps 3-6           # solo elaborazione (dati locali)

Credenziali download (step 2), via variabili d'ambiente:
  CDSE_USER / CDSE_PASS   (account gratuito su dataspace.copernicus.eu)
"""
import os, sys, re, glob, json, time, zipfile, argparse
import urllib.request, urllib.parse
import numpy as np

# matplotlib/scipy/plotly importati lazy negli step che li usano, cosi' gli step
# 1-2 (ricerca/download) girano anche su una macchina senza quei pacchetti.

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

C = 299_792_458.0
DBG = False


def log(*a):
    print(*a, flush=True)


def dbg(*a):
    if DBG:
        print("   [debug]", *a, flush=True)


def dms2dec(d, m, s, h):
    v = float(d) + float(m) / 60 + float(s) / 3600
    return -v if str(h).upper() in ("S", "W") else v


# ===================================================================== STEP 1
def step1_cerca(aoi, start, end, max_records, pol, outdir):
    """Cerca i prodotti Sentinel-1 SLC che contengono l'AOI (OData CDSE).
    Ritorna la lista [{Name, Id, ContentLength, Start}] ordinata per data."""
    lonW, lonE, latS, latN = aoi
    poly = (f"POLYGON(({lonW} {latS},{lonE} {latS},{lonE} {latN},"
            f"{lonW} {latN},{lonW} {latS}))")
    flt = ("Collection/Name eq 'SENTINEL-1' "
           f"and OData.CSC.Intersects(area=geography'SRID=4326;{poly}') "
           "and contains(Name,'SLC') "
           f"and ContentDate/Start gt {start}T00:00:00.000Z "
           f"and ContentDate/Start lt {end}T23:59:59.000Z")
    q = {"$filter": flt, "$orderby": "ContentDate/Start asc",
         "$top": str(max_records)}
    url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?" + \
          urllib.parse.urlencode(q)
    dbg("query:", url)
    log(f"[1] ricerca SLC che contengono l'AOI "
        f"(lon {lonW:.5f}..{lonE:.5f}, lat {latS:.5f}..{latN:.5f}), {start}..{end}")
    try:
        with urllib.request.urlopen(urllib.request.Request(url), timeout=90) as r:
            data = json.load(r)
    except Exception as ex:
        log(f"    ! ricerca catalogo fallita ({ex}); proseguo con i dati locali")
        return []
    prods = [{"Name": f["Name"], "Id": f["Id"],
              "ContentLength": f.get("ContentLength", 0),
              "Start": f.get("ContentDate", {}).get("Start", "")}
             for f in data.get("value", [])]
    log(f"[1] trovati {len(prods)} prodotti SLC che contengono le coordinate "
        f"{'(>=12 OK)' if len(prods) >= 12 else '(ATTENZIONE: < 12)'}")
    for p in prods[:20]:
        log(f"      - {p['Name']}  ({p['ContentLength']/1e9:.1f} GB)")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "prodotti_trovati.json"), "w", encoding="utf-8") as f:
        json.dump(prods, f, indent=2)
    log(f"    elenco salvato in {os.path.join(outdir, 'prodotti_trovati.json')}")
    return prods


# ===================================================================== STEP 2
def _cdse_token(user, pwd):
    url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    body = urllib.parse.urlencode({"client_id": "cdse-public", "grant_type": "password",
                                   "username": user, "password": pwd}).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=body), timeout=60) as r:
        return json.load(r)["access_token"]


def _scarica_ripristinabile(url, dest, get_token, expected=0, n_try=6):
    """Scarica `url` in `dest` riprendendo da dove si era interrotto.

    Il file viene scritto in `dest + '.part'` in append: a ogni tentativo si
    legge la dimensione gia' scaricata e si chiede al server solo il resto con
    l'header Range. Tollera errori transitori (SSL/timeout/reset) ritentando
    con backoff e ricavando ogni volta un token fresco da `get_token()`
    (i token CDSE scadono e i download da ~8 GB durano a lungo).
    A download completato `.part` viene rinominato in `dest`.
    Ritorna True se il file e' stato scaricato per intero."""
    part = dest + ".part"
    for att in range(1, n_try + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if expected and have >= expected:        # gia' completo da un giro precedente
            break
        headers = {"Authorization": f"Bearer {get_token()}"}
        if have:
            headers["Range"] = f"bytes={have}-"
            log(f"    riprendo da {have/1e9:.2f} GB"
                f"{f' / {expected/1e9:.1f} GB' if expected else ''} (tentativo {att}/{n_try})")
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=1800) as r:
                # se il server ignora il Range (200 invece di 206) si riparte da capo
                mode = "ab" if (have and r.status == 206) else "wb"
                if mode == "wb":
                    have = 0
                with open(part, mode) as fo:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        fo.write(chunk)
                        have += len(chunk)
            if not expected or have >= expected:
                break
            log(f"    ! ricevuti {have/1e9:.2f}/{expected/1e9:.1f} GB; ritento")
        except Exception as ex:
            wait = min(60, 5 * att)
            log(f"    ! interrotto ({ex}); ritento tra {wait}s "
                f"({att}/{n_try}, {os.path.getsize(part)/1e9 if os.path.exists(part) else 0:.2f} GB salvati)")
            time.sleep(wait)
    else:
        log(f"    ! download non completato dopo {n_try} tentativi"); return False
    if expected and os.path.getsize(part) < expected:
        log("    ! download incompleto"); return False
    os.replace(part, dest)
    return True


def step2_scarica(prods, n_want, pol, stackdir, user=None, pwd=None):
    """Scarica i prodotti e ne estrae i measurement TIFF (+ annotation) della
    polarizzazione richiesta nello stack locale. Credenziali CDSE: dai parametri
    --cdse-user/--cdse-pass oppure dalle variabili d'ambiente CDSE_USER/CDSE_PASS.
    Se mancano le credenziali (o la lista e' vuota), salta e usa quel che c'e'."""
    os.makedirs(stackdir, exist_ok=True)
    have = sorted(glob.glob(os.path.join(stackdir, f"*slc-{pol}-*.tiff")))
    if have:
        log(f"[2] gia' presenti {len(have)} TIFF '{pol}' in {stackdir} (riuso)")
    user = user or os.environ.get("CDSE_USER")
    pwd = pwd or os.environ.get("CDSE_PASS")
    if not prods:
        log("[2] nessun prodotto dalla ricerca: salto il download")
        return have
    if not (user and pwd):
        log("[2] credenziali CDSE assenti: salto il download.")
        log("    -> i SAFE pesano ~8 GB l'uno. Per scaricarli davvero, passa")
        log("       --cdse-user <user> --cdse-pass <pass>  (o le env CDSE_USER/CDSE_PASS)")
        return have
    try:
        _cdse_token(user, pwd)            # verifica subito le credenziali
    except Exception as ex:
        log(f"[2] login CDSE fallito ({ex}); salto il download"); return have
    get_token = lambda: _cdse_token(user, pwd)   # token fresco a ogni (ri)tentativo

    got = list(have)
    for p in prods:
        if len([g for g in got]) >= n_want:
            break
        safe_zip = os.path.join(stackdir, p["Name"].replace(".SAFE", "") + ".zip")
        url = (f"https://catalogue.dataspace.copernicus.eu/odata/v1/"
               f"Products({p['Id']})/$value")
        log(f"[2] scarico {p['Name']} ({p['ContentLength']/1e9:.1f} GB) ...")
        if not _scarica_ripristinabile(url, safe_zip, get_token,
                                       expected=int(p.get("ContentLength", 0) or 0)):
            log("    ! download fallito; continuo"); continue
        # estrai dal SAFE i measurement tiff + annotation della polarizzazione
        try:
            with zipfile.ZipFile(safe_zip) as z:
                for nm in z.namelist():
                    low = nm.lower()
                    if f"-{pol}-" in low and low.endswith(".tiff") and "/measurement/" in low:
                        dst = os.path.join(stackdir, os.path.basename(low))
                        with z.open(nm) as zi, open(dst, "wb") as fo:
                            fo.write(zi.read())
                        got.append(dst)
                    if f"-{pol}-" in low and low.endswith(".xml") and "/annotation/" in low \
                            and "/calibration/" not in low:
                        dst = os.path.join(stackdir, os.path.basename(low).replace(
                            ".xml", ".annotation.xml"))
                        with z.open(nm) as zi, open(dst, "wb") as fo:
                            fo.write(zi.read())
            os.remove(safe_zip)
        except Exception as ex:
            log(f"    ! estrazione SAFE fallita ({ex})")
    log(f"[2] TIFF '{pol}' disponibili nello stack: {len(got)}")
    return got


# ===================================================================== STEP 3
def _trova_annotation(tif, stack, pol):
    m = re.search(rf"(s1[abcd])-iw(\d)-slc-{pol}-(\d{{8}})",
                  os.path.basename(tif).lower())
    if not m:
        return None
    s, sw, dt = m.groups()
    for x in glob.glob(os.path.join(stack, "*.annotation.xml")):
        xb = os.path.basename(x).lower()
        if s in xb and f"iw{sw}" in xb and dt in xb:
            return x


def _leggi_geometria(ann):
    import xml.etree.ElementTree as ET
    r = ET.parse(ann).getroot()
    f = lambda p: float(r.find(p).text)
    g = {"f_em": f(".//radarFrequency"), "dt_az": f(".//azimuthTimeInterval"),
         "slantRangeTime": f(".//imageInformation/slantRangeTime"),
         "incid_mid": f(".//incidenceAngleMidSwath"),
         "dR": f(".//rangePixelSpacing"), "dA": f(".//azimuthPixelSpacing")}
    g["R_near"] = C * g["slantRangeTime"] / 2.0
    return g


def _leggi_complesso(ds, win):
    a = ds.read(1, window=win)
    if np.iscomplexobj(a):
        return a.astype(np.complex64)
    if a.dtype.names and len(a.dtype.names) == 2:
        n0, n1 = a.dtype.names
        return a[n0].astype(np.float32) + 1j * a[n1].astype(np.float32)
    return a.astype(np.complex64)


def step3_estrai_box(stackdir, pol, aoi, boxpath, fallback_box):
    """Geolocalizza l'AOI su ogni TIFF via i GCP e ritaglia lo stack complesso."""
    lonW, lonE, latS, latN = aoi
    corners = [(lonW, latN), (lonE, latN), (lonE, latS), (lonW, latS)]
    tiffs = sorted(glob.glob(os.path.join(stackdir, f"*-iw[1-3]-slc-{pol}-*.tiff"))) or \
            sorted(glob.glob(os.path.join(stackdir, f"*slc-{pol}-*.tiff")))
    if not tiffs:
        if os.path.exists(fallback_box):
            log(f"[3] nessun TIFF in {stackdir}: uso il box gia' ritagliato "
                f"{fallback_box} (fallback)")
            d = np.load(fallback_box, allow_pickle=True)
            np.savez_compressed(boxpath, **{k: d[k] for k in d.files})
            return boxpath
        sys.exit(f"[3] nessun TIFF in {stackdir} e nessun fallback {fallback_box}")

    import rasterio
    from rasterio.windows import Window
    from rasterio.transform import GCPTransformer
    chips, nomi, geom = [], [], None
    for tif in tiffs:
        ann = _trova_annotation(tif, stackdir, pol)
        if ann is None:
            log(f"    ({os.path.basename(tif)}: annotation mancante, skip)"); continue
        g = _leggi_geometria(ann)
        with rasterio.open(tif) as ds:
            if not ds.gcps[0]:
                log(f"    ({os.path.basename(tif)}: GCP assenti, skip)"); continue
            with GCPTransformer(ds.gcps[0]) as tf:
                rc = [tf.rowcol(lon, lat) for lon, lat in corners]
            if not all(0 <= r < ds.height and 0 <= c < ds.width for r, c in rc):
                log(f"    ({os.path.basename(tif)}: AOI fuori dal raster, skip)"); continue
            rows = [r for r, c in rc]; cols = [c for r, c in rc]
            r0 = max(0, int(np.floor(min(rows)))); c0 = max(0, int(np.floor(min(cols))))
            r1 = min(ds.height, int(np.ceil(max(rows))) + 1)
            c1 = min(ds.width, int(np.ceil(max(cols))) + 1)
            chip = _leggi_complesso(ds, Window(c0, r0, c1 - c0, r1 - r0))
        valid = float(np.mean(np.abs(chip) > 0) * 100)
        log(f"    >> {os.path.basename(tif)} {chip.shape[1]}x{chip.shape[0]}px "
            f"(rng x az) valid {valid:.0f}%")
        chips.append(chip); nomi.append(os.path.basename(tif)); geom = g
    if not chips:
        sys.exit("[3] nessuna scena copre l'AOI")
    H = min(c.shape[0] for c in chips); W = min(c.shape[1] for c in chips)
    arr = np.stack([c[:H, :W] for c in chips], 0)
    np.savez_compressed(boxpath, vv=arr, vv_scene=np.array(nomi),
                        corners_lonlat=np.array(corners),
                        dR=geom["dR"], dA=geom["dA"], incid_mid=geom["incid_mid"],
                        f_em=geom["f_em"], R_near=geom["R_near"], dt_az=geom["dt_az"])
    log(f"[3] stack ritagliato {arr.shape} (n_scene, az, rng) -> {boxpath}")
    return boxpath


# ===================================================================== STEP 4
def _data_scena(nome):
    m = re.search(r"-(?:vv|vh|hh|hv)-(\d{8})t", nome.lower()) or \
        re.search(r"(\d{8})t\d{6}", nome.lower())
    return m.group(1) if m else "00000000"


def step4_array_12(boxpath, n_layers, outdir):
    """Costruisce l'array a n_layers strati della COMPONENTE VERTICALE dello
    spostamento. Ogni strato = un interferogramma (coppia di scene) a baseline
    temporale crescente: phi = arg(s_j * conj(s_i)); LOS = -lambda/(4pi)*phi;
    verticale = LOS / cos(theta_inc). Con 12+ scene bastano le coppie sequenziali;
    con poche scene si usano le coppie a baseline minima fino ad arrivare a 12."""
    d = np.load(boxpath, allow_pickle=True)
    vv = d["vv"]
    nomi = [str(x) for x in d["vv_scene"]]
    inc = float(d["incid_mid"]); f_em = float(d["f_em"])
    dR = float(d["dR"]); dA = float(d["dA"])
    lam = C / f_em
    K = lam / (4 * np.pi) * 1000.0                 # mm per radiante
    cth = np.cos(np.radians(inc))                  # proiezione LOS -> verticale
    gr = dR / np.sin(np.radians(inc))              # passo ground-range [m]

    date = [_data_scena(n) for n in nomi]
    order = np.argsort(date)
    vv = vv[order]; nomi = [nomi[i] for i in order]; date = [date[i] for i in order]
    nS, ny, nx = vv.shape
    dnum = np.array([np.datetime64(f"{x[:4]}-{x[4:6]}-{x[6:8]}") for x in date])

    # coppie (i<j) ordinate per baseline temporale, poi per data
    pairs = [(i, j, int((dnum[j] - dnum[i]) / np.timedelta64(1, "D")))
             for i in range(nS) for j in range(i + 1, nS)]
    pairs.sort(key=lambda t: (t[2], t[0]))
    if nS - 1 >= n_layers:                          # abbastanza scene: sequenziali
        pairs = [(i, i + 1, int((dnum[i+1]-dnum[i])/np.timedelta64(1, "D")))
                 for i in range(nS - 1)][:n_layers]
    else:
        pairs = pairs[:n_layers]
    if len(pairs) < n_layers:
        log(f"[4] ATTENZIONE: solo {len(pairs)} coppie disponibili (< {n_layers})")

    disp = np.empty((len(pairs), ny, nx), np.float32)   # componente verticale [mm]
    phi_st = np.empty((len(pairs), ny, nx), np.float32) # fase interferometrica [rad]
    etich = []
    for s, (i, j, bt) in enumerate(pairs):
        phi = np.angle(vv[j] * np.conj(vv[i]))
        phi_st[s] = phi
        disp[s] = (-K * phi / cth).astype(np.float32)
        etich.append(f"{date[i]}~{date[j]} ({bt}d)")
    coh = np.abs(np.exp(1j * phi_st).mean(0))

    x = np.arange(nx) * gr                          # asse range a terra [m]
    y = np.arange(ny) * dA                          # asse azimuth [m]
    os.makedirs(outdir, exist_ok=True)
    np.savez_compressed(os.path.join(outdir, "array3d_12strati.npz"),
                        disp=disp, phi=phi_st, x=x, y=y, coh=coh.astype(np.float32),
                        etich=np.array(etich), inc=inc, lam=lam, gr=gr, dA=dA)
    log(f"[4] array componente verticale {disp.shape} (strati, y, x); "
        f"vert [{disp.min():.1f},{disp.max():.1f}] mm; coh media {coh.mean():.2f}")
    dbg("strati:", etich)

    # figura riassuntiva degli strati
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nrow = int(np.ceil(len(pairs) / 4))
    fig, axes = plt.subplots(nrow, 4, figsize=(15, 3.2 * nrow))
    for k, ax in enumerate(np.atleast_1d(axes).ravel()):
        if k < len(pairs):
            im = ax.imshow(disp[k], cmap="RdBu_r", vmin=-20, vmax=20,
                           aspect="auto", origin="lower")
            ax.set_title(etich[k], fontsize=8)
        else:
            ax.axis("off")
    fig.colorbar(im, ax=axes, shrink=0.7, label="spostamento verticale [mm]")
    fig.suptitle(f"Step 4 — array a {len(pairs)} strati (componente verticale)")
    fig.savefig(os.path.join(outdir, "step4_strati.png"), dpi=120, bbox_inches="tight")
    plt.close(fig)
    return disp, phi_st, x, y, coh


# ----------------------------------------------------------- DEM (altezze reali)
def _fetch_dem(aoi, ny, nx, outdir, pyramid):
    """Quota reale del terreno (DEM Copernicus via Open-Meteo) sulla griglia
    pixel (ny,nx). Opzionale overlay della geometria della piramide di Chefren
    se l'AOI e' centrata sul plateau di Giza."""
    cache = os.path.join(outdir, "dem.npz")
    lonW, lonE, latS, latN = aoi
    j = np.arange(nx); i = np.arange(ny)
    lon = lonW + (j / max(nx - 1, 1)) * (lonE - lonW)
    lat = latS + (i / max(ny - 1, 1)) * (latN - latS)
    if os.path.exists(cache):
        # la cache contiene il DEM NUDO (senza overlay piramide)
        z0 = np.load(cache)["z0"].astype(float)
        log(f"[DEM] riuso {cache}: quota terreno {z0.min():.0f}..{z0.max():.0f} m s.l.m.")
    else:
        GC = 10
        CLAT, CLON = np.meshgrid(np.linspace(latS, latN, GC),
                                 np.linspace(lonW, lonE, GC), indexing="ij")
        qs = ",".join(f"{v:.6f}" for v in CLAT.ravel())
        qo = ",".join(f"{v:.6f}" for v in CLON.ravel())
        url = f"https://api.open-meteo.com/v1/elevation?latitude={qs}&longitude={qo}"
        ce = None
        for att in range(5):
            try:
                with urllib.request.urlopen(url, timeout=40) as r:
                    ce = np.array(json.load(r)["elevation"], float).reshape(GC, GC)
                break
            except Exception as ex:
                log(f"[DEM] retry {att} ({ex})"); time.sleep(3 + 3 * att)
        if ce is None:
            log("[DEM] download fallito: uso quota 0 (livello mare)")
            z0 = np.zeros((ny, nx))
        else:
            from scipy.interpolate import RegularGridInterpolator
            itp = RegularGridInterpolator(
                (np.linspace(latS, latN, GC), np.linspace(lonW, lonE, GC)), ce,
                bounds_error=False, fill_value=None)
            LAT, LON = np.meshgrid(lat, lon, indexing="ij")
            z0 = itp(np.stack([LAT.ravel(), LON.ravel()], -1)).reshape(ny, nx)
        log(f"[DEM] quota {z0.min():.0f}..{z0.max():.0f} m s.l.m. (media {z0.mean():.0f})")
    # salva SEMPRE il DEM nudo in cache (l'overlay piramide e' ricalcolato ogni volta)
    np.savez_compressed(cache, z0=z0, lon=lon, lat=lat)

    if pyramid:
        # geometria reale piramide di Chefren (base 215.5 m, altezza attuale 136.4 m)
        latm = np.radians((latN + latS) / 2)
        Lx = (lonE - lonW) * 111320 * np.cos(latm)
        Ly = (latN - latS) * 110540
        x_px = np.linspace(0, Lx, nx); y_px = np.linspace(0, Ly, ny)
        H_PYR, BASE = 136.4, 215.5; half = BASE / 2.0
        XX, YY = np.meshgrid(x_px, y_px)
        cheb = np.maximum(np.abs(XX - Lx / 2), np.abs(YY - Ly / 2))
        z0 = z0 + np.clip(H_PYR * (1.0 - cheb / half), 0.0, None)
        log(f"[DEM] overlay piramide Chefren -> apice {z0.max():.0f} m s.l.m.")
    return z0


# ===================================================================== STEP 5
def step5_grafico_onde(disp, phi, x, y, aoi, outdir, pyramid,
                       NZ=256, ZTOP=210.0, ZBOT=-200.0, STEP_X=4, STEP_Y=2,
                       WIGGLE=14.0):
    """Grafico 3D dell'array: i numeri dei n_layers strati di ogni pixel sono i
    coefficienti armonici complessi (ampiezza=|vert|, fase=fase interferometrica)
    di una forma d'onda sviluppata in profondita' via anti-trasformata di Fourier.
    Le onde partono dalla QUOTA REALE del suolo (DEM)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go

    nS, ny, nx = disp.shape
    latm = np.radians((aoi[3] + aoi[2]) / 2)
    Lx = (aoi[1] - aoi[0]) * 111320 * np.cos(latm)
    Ly = (aoi[3] - aoi[2]) * 110540
    x_m = np.linspace(0, Lx, nx); y_m = np.linspace(0, Ly, ny)
    z0 = _fetch_dem(aoi, ny, nx, outdir, pyramid)

    ZD = ZTOP - ZBOT
    z = np.linspace(0, ZD, NZ); nf = NZ // 2 + 1
    coeff = (np.abs(disp) * np.exp(1j * phi)).astype(complex)   # (nS, ny, nx)
    spec = np.zeros((ny * nx, nf), complex)
    spec[:, 1:1 + nS] = coeff.reshape(nS, ny * nx).T
    w = (np.fft.irfft(spec, n=NZ, axis=1) * NZ).reshape(ny, nx, NZ).astype(np.float32)
    wmax = float(np.abs(w).max()) or 1.0
    log(f"[5] forme d'onda: anti-FT di {nS} armoniche -> w {w.shape}, |w|max={wmax:.1f}; "
        f"base {Lx:.0f}x{Ly:.0f} m, profondita' {ZD:.0f} m")

    np.savez_compressed(os.path.join(outdir, "forme_donda.npz"),
                        w=w, z=z, x=x_m, y=y_m, z0=z0, Lx=Lx, Ly=Ly, zdepth=ZD)

    # "molle" sviluppate verso il basso dalla quota reale
    ii = np.arange(0, ny, STEP_Y); jj = np.arange(0, nx, STEP_X)
    sc = WIGGLE / wmax
    X, Y, Z, Cc = [], [], [], []
    gap = np.array([np.nan])
    for i in ii:
        for j in jj:
            wp = w[i, j]
            X += [x_m[j] + sc * wp, gap]; Y += [np.full(NZ, y_m[i]), gap]
            Z += [z0[i, j] - z, gap]; Cc += [np.abs(wp), gap]
    X = np.concatenate(X); Y = np.concatenate(Y); Z = np.concatenate(Z); Cc = np.concatenate(Cc)

    Xs, Ys = np.meshgrid(x_m, y_m)
    # altezza/profondita' del suolo rispetto al LIVELLO 0 (s.l.m.)
    h_suolo = np.where(z0 >= 0, z0, 0.0)        # altezza dal livello 0 [m]
    p_suolo = np.where(z0 < 0, -z0, 0.0)        # profondita' dal livello 0 [m]
    txt = np.array([[f"x={x_m[j]:.0f} m, y={y_m[i]:.0f} m<br>"
                     f"quota {z0[i,j]:+.0f} m s.l.m.<br>"
                     f"altezza dal livello 0: {h_suolo[i,j]:.0f} m<br>"
                     f"profondita' dal livello 0: {p_suolo[i,j]:.0f} m"
                     for j in range(nx)] for i in range(ny)]).ravel()
    plane = go.Surface(x=Xs, y=Ys, z=z0, surfacecolor=z0, colorscale="Earth",
                       opacity=0.7, colorbar=dict(title="quota s.l.m. [m]", x=-0.08),
                       name="superficie reale (DEM)")
    # piano di riferimento a quota 0 (livello del mare)
    zero_plane = go.Surface(x=Xs, y=Ys, z=np.zeros_like(z0),
                            showscale=False, opacity=0.25,
                            colorscale=[[0, "gray"], [1, "gray"]],
                            name="livello 0 (s.l.m.)")
    springs = go.Scatter3d(x=X, y=Y, z=Z, mode="lines", connectgaps=False,
                           line=dict(width=3, color=Cc, colorscale="Turbo",
                                     colorbar=dict(title="|w|")),
                           name=f"forme d'onda ({nS} armoniche)")
    ground = go.Scatter3d(x=Xs.ravel(), y=Ys.ravel(), z=z0.ravel(), mode="markers",
                          marker=dict(size=1.2, color="red"), name="suolo (DEM)",
                          text=txt, hoverinfo="text")
    fig = go.Figure([plane, zero_plane, springs, ground])
    log(f"[5] altezza max del suolo dal livello 0: {h_suolo.max():.0f} m; "
        f"onde sviluppate fino a quota {ZBOT:.0f} m s.l.m. "
        f"(= {abs(ZBOT):.0f} m di profondita' sotto il livello 0)")
    fig.update_layout(
        title=f"Step 5 — grafico 3D dell'array a {nS} strati: forme d'onda dalla "
              f"quota reale ({z0.min():.0f}-{z0.max():.0f} m s.l.m.), livello 0 = mare",
        scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
                   zaxis_title="quota reale s.l.m. [m]", zaxis=dict(range=[ZBOT, ZTOP]),
                   aspectmode="manual",
                   aspectratio=dict(x=1.0, y=Ly / Lx, z=ZD / Lx)))
    fig.write_html(os.path.join(outdir, "step5_forme_donda_3d.html"),
                   include_plotlyjs="cdn",
                   config=dict(scrollZoom=True, displayModeBar=True))

    # anteprima statica + B-scan
    figm = plt.figure(figsize=(9, 8)); ax = figm.add_subplot(111, projection="3d")
    ax.plot_surface(Xs, Ys, z0, cmap="terrain", alpha=0.4, linewidth=0)
    for i in ii:
        for j in jj:
            ax.plot(x_m[j] + sc * w[i, j], np.full(NZ, y_m[i]), z0[i, j] - z, lw=0.6)
    ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
    ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZD))
    ax.set_title(f"Step 5 — forme d'onda dalla quota reale ({nS} armoniche)")
    figm.savefig(os.path.join(outdir, "step5_forme_donda_3d.png"), dpi=120,
                 bbox_inches="tight")
    plt.close(figm)

    i0 = ny // 2
    figb, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(w[i0].T, cmap="seismic", aspect="auto", origin="upper",
                   extent=[0, Lx, ZD, 0], vmin=-wmax, vmax=wmax)
    ax.set_xlabel("x E-W [m]"); ax.set_ylabel("profondita' [m]")
    ax.set_title(f"Step 5 — B-scan forme d'onda (linea y={y_m[i0]:.0f} m)")
    figb.colorbar(im, ax=ax, label="ampiezza w")
    figb.savefig(os.path.join(outdir, "step5_bscan.png"), dpi=120, bbox_inches="tight")
    plt.close(figb)
    log(f"[5] salvati step5_forme_donda_3d.html/.png, step5_bscan.png")
    return w, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT


# ===================================================================== STEP 6
def step6_variazioni(w, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT, outdir,
                     DMIN=0.035, DMAX=0.06):
    """Punti dove le forme d'onda variano in FREQUENZA (=> densita'): frequenza
    istantanea lungo la profondita' (Hilbert) -> |d(freq)/dz|; i punti nella
    banda densita' sono disegnati in 3D alle ALTEZZE REALI (quota DEM - profondita')."""
    from scipy.signal import hilbert
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go

    ny, nx, nz = w.shape
    dz = z[1] - z[0]
    an = hilbert(w, axis=2)
    inst_f = np.gradient(np.unwrap(np.angle(an), axis=2), dz, axis=2) / (2 * np.pi)
    dfreq = np.abs(np.gradient(inst_f, dz, axis=2))
    ii, jj, kk = np.where((dfreq > DMIN) & (dfreq <= DMAX))
    px = x_m[jj]; py = y_m[ii]; pz = z0[ii, jj] - z[kk]; pv = dfreq[ii, jj, kk]
    m = (pz >= ZBOT) & (pz <= ZTOP)
    px, py, pz, pv = px[m], py[m], pz[m], pv[m]
    # altezza dal livello 0 (quota>0) e profondita' dal livello 0 (quota<0)
    altezza = np.where(pz >= 0, pz, 0.0)
    profondita = np.where(pz < 0, -pz, 0.0)
    log(f"[6] banda densita' {DMIN}<|Δfreq|<={DMAX}: {len(px)} punti su {dfreq.size} "
        f"(max {dfreq.max():.4f})")
    n_sopra = int((pz >= 0).sum()); n_sotto = int((pz < 0).sum())
    log(f"[6] rispetto al livello 0: {n_sopra} punti sopra (altezza max "
        f"{altezza.max():.0f} m), {n_sotto} punti sotto (profondita' max "
        f"{profondita.max():.0f} m)")
    np.savez_compressed(os.path.join(outdir, "variazioni_frequenza.npz"),
                        px=px, py=py, pz=pz, pv=pv,
                        altezza_dal_livello0=altezza, profondita_dal_livello0=profondita,
                        dmin=DMIN, dmax=DMAX)

    # hover: quota s.l.m. + altezza/profondita' dal livello 0 per ogni punto
    txt = []
    for a, b, c, d in zip(px, py, pz, pv):
        rel = (f"altezza dal livello 0: {c:.0f} m" if c >= 0
               else f"profondita' dal livello 0: {-c:.0f} m")
        txt.append(f"x={a:.0f} m, y={b:.0f} m<br>quota {c:+.0f} m s.l.m.<br>"
                   f"{rel}<br>|Δfreq|={d:.4f}")

    Xg, Yg = np.meshgrid(x_m, y_m)
    ground = go.Scatter3d(x=Xg.ravel(), y=Yg.ravel(), z=z0.ravel(), mode="markers",
                          marker=dict(size=1.2, color="red"), name="suolo (DEM)")
    # piano di riferimento a quota 0 (livello del mare)
    zero_plane = go.Surface(x=Xg, y=Yg, z=np.zeros_like(Xg),
                            showscale=False, opacity=0.25,
                            colorscale=[[0, "gray"], [1, "gray"]],
                            name="livello 0 (s.l.m.)")
    fig = go.Figure([go.Scatter3d(
        x=px, y=py, z=pz, mode="markers", name="variazioni freq.",
        text=txt, hoverinfo="text",
        marker=dict(size=1.2, color=pv, colorscale="Inferno", opacity=0.85,
                    colorbar=dict(title="|Δfreq| (densita')"))), zero_plane, ground])
    fig.update_layout(
        title=f"Step 6 — punti di variazione di frequenza (=> densita') alle altezze "
              f"reali — base {Lx:.0f}×{Ly:.0f} m (livello 0 = mare)",
        scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
                   zaxis_title="quota s.l.m. [m] (+ = altezza, − = profondita' dal livello 0)",
                   zaxis=dict(range=[ZBOT, ZTOP]),
                   aspectmode="manual",
                   aspectratio=dict(x=1.0, y=Ly / Lx, z=(ZTOP - ZBOT) / Lx)))
    fig.write_html(os.path.join(outdir, "step6_variazioni_3d.html"),
                   include_plotlyjs="cdn",
                   config=dict(scrollZoom=True, displayModeBar=True))

    figm = plt.figure(figsize=(8, 9)); ax = figm.add_subplot(111, projection="3d")
    ax.scatter(Xg.ravel(), Yg.ravel(), z0.ravel(), c="red", s=2, depthshade=False)
    sc = ax.scatter(px, py, pz, c=pv, cmap="inferno", s=1.5)
    ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
    ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZTOP - ZBOT))
    ax.set_title("Step 6 — variazioni di frequenza alle altezze reali")
    figm.colorbar(sc, ax=ax, shrink=0.5, label="|Δfreq|")
    figm.savefig(os.path.join(outdir, "step6_variazioni_3d.png"), dpi=120,
                 bbox_inches="tight")
    plt.close(figm)
    log(f"[6] salvati step6_variazioni_3d.html/.png")


# ===================================================================== main
def parse_steps(s):
    if "-" in s:
        a, b = s.split("-"); return set(range(int(a), int(b) + 1))
    return {int(x) for x in s.split(",")}


def main():
    global DBG
    ap = argparse.ArgumentParser(description="Pipeline unica piramide (step 1-6)")
    ap.add_argument("--nw", nargs=4, metavar=("D", "M", "S", "H"),
                    default=["29", "58", "48.0", "N"], help="angolo NW lat (DMS)") 
    ap.add_argument("--nw-lon", nargs=4, metavar=("D", "M", "S", "H"),
                    default=["31", "7", "58.4", "E"], help="angolo NW lon (DMS)")  
    ap.add_argument("--se", nargs=4, metavar=("D", "M", "S", "H"),
                    default=["29", "58", "41.0", "N"], help="angolo SE lat (DMS)")
    ap.add_argument("--se-lon", nargs=4, metavar=("D", "M", "S", "H"),
                    default=["31", "8", "07.7", "E"], help="angolo SE lon (DMS)")
    ap.add_argument("--pol", default="vv")
    ap.add_argument("--layers", type=int, default=12, help="numero di strati (>=12)")
    ap.add_argument("--start", default="2025-01-01", help="data inizio ricerca YYYY-MM-DD")
    ap.add_argument("--end", default="2026-06-25", help="data fine ricerca YYYY-MM-DD")
    ap.add_argument("--max-search", type=int, default=40)
    ap.add_argument("--stack", default=None, help="cartella TIFF scaricati")
    ap.add_argument("--outdir", default=None, help="cartella output")
    ap.add_argument("--download", action="store_true", help="scarica davvero (servono credenziali CDSE)")
    ap.add_argument("--cdse-user", default=None,
                    help="username/email CDSE (altrimenti env CDSE_USER)")
    ap.add_argument("--cdse-pass", default=None,
                    help="password CDSE (altrimenti env CDSE_PASS)")
    ap.add_argument("--steps", default="1-6", help="es. 1-6 oppure 3,4,5")
    ap.add_argument("--dmin", type=float, default=0.035,
                    help="variazione di frequenza MINIMA della banda densita' (step 6)")
    ap.add_argument("--dmax", type=float, default=0.06,
                    help="variazione di frequenza MASSIMA della banda densita' (step 6)")
    ap.add_argument("--no-pyramid", action="store_true", help="non sovrapporre la piramide al DEM")
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args()
    DBG = a.debug

    here = os.path.dirname(os.path.abspath(__file__))
    outdir = a.outdir or os.path.join(here, "unificato")
    stackdir = a.stack or os.path.join(here, "stack_slc")
    fallback_box = os.path.join(here, "box.npz")
    boxpath = os.path.join(outdir, "box.npz")
    os.makedirs(outdir, exist_ok=True)
    pyramid = not a.no_pyramid
    steps = parse_steps(a.steps)

    lat_nw = dms2dec(*a.nw[:3], a.nw[3]); lon_nw = dms2dec(*a.nw_lon[:3], a.nw_lon[3])
    lat_se = dms2dec(*a.se[:3], a.se[3]); lon_se = dms2dec(*a.se_lon[:3], a.se_lon[3])
    latN, latS = max(lat_nw, lat_se), min(lat_nw, lat_se)
    lonW, lonE = min(lon_nw, lon_se), max(lon_nw, lon_se)
    aoi = (lonW, lonE, latS, latN)

    log("=" * 70)
    log(f"PIPELINE UNIFICATA — AOI lon {lonW:.5f}..{lonE:.5f}, lat {latS:.5f}..{latN:.5f}")
    log(f"step {sorted(steps)} | pol {a.pol} | strati {a.layers} | outdir {outdir}")
    log("=" * 70)

    prods = []
    if 1 in steps:
        prods = step1_cerca(aoi, a.start, a.end, a.max_search, a.pol, outdir)
    if 2 in steps:
        # senza --download non si scaricano gli ~8 GB/prodotto: si riusa lo stack
        step2_scarica(prods if a.download else [], a.layers, a.pol, stackdir,
                      user=a.cdse_user, pwd=a.cdse_pass)
    if 3 in steps:
        step3_estrai_box(stackdir, a.pol, aoi, boxpath, fallback_box)
    elif not os.path.exists(boxpath) and os.path.exists(fallback_box):
        boxpath = fallback_box

    if steps & {4, 5, 6}:
        disp, phi, x, y, coh = step4_array_12(boxpath, a.layers, outdir)
        if steps & {5, 6}:
            (w, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT) = \
                step5_grafico_onde(disp, phi, x, y, aoi, outdir, pyramid)
            if 6 in steps:
                step6_variazioni(w, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT, outdir,
                                 DMIN=a.dmin, DMAX=a.dmax)

    log("=" * 70)
    log(f"FATTO. Output in: {outdir}")
    log("=" * 70)


if __name__ == "__main__":
    main()

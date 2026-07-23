#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramide_acustica_vh.py
========================
Variante "acustica", SOLO canale VH, della pipeline di piramide_unificato.py.

*** COSA FA ***
Usa ESCLUSIVAMENTE la polarizzazione incrociata VH degli stessi TIFF Sentinel-1
gia' scaricati da piramide_unificato.py (stack_slc/), li tratta come un segnale
di VIBRAZIONE riflesso dalla struttura interna e ne ricava:
  1) un array a N strati della componente verticale dello spostamento VH
     (DInSAR classica: fase interferometrica tra coppie di acquisizioni ->
     spostamento verticale [mm]); N = TUTTE le scene VH disponibili nello
     stack meno una (--layers None, default: usa tutti i file VH presenti,
     ordinati per data, un layer per coppia di date consecutive -> nessuna
     scena esclusa), oppure un numero fisso con --layers;
  2) quei N strati reinterpretati come CAMPIONI di un'onda acustica composita
     (ampiezza = |spostamento|, fase = fase interferometrica) -> anti-trasformata
     di Fourier -> forma d'onda sviluppata in profondita' + inviluppo/frequenza
     istantanea di Hilbert (stessa matematica dello step 5 di
     piramide_unificato.py, qui ristretta a VH e presentata col lessico
     "vibrazione acustica");
  3) la tomografia Doppler VERA del paper Biondi & Malanga (steering matrix,
     h(z) = A^H*Y), anch'essa ristretta al canale VH;
  4) le variazioni PIENO<->VUOTO (tetto/fondo) con una selezione COERENTE
     per pixel (variazioni_coerenti_vh): a differenza dello step 6 standard
     di piramide_unificato.py, che seleziona indipendentemente le N
     variazioni piu' forti per TIPO (tetto/fondo) su tutto il volume e puo'
     lasciare un tetto (PIENO->VUOTO) senza il suo fondo (VUOTO->PIENO) o
     viceversa (verificato sui dati VH reali: 77% dei pixel selezionati
     risultava orfano), qui si selezionano SOLO coppie complete
     tetto+fondo dello STESSO pixel: ogni "vuoto" disegnato ha sempre sia
     il tetto sia il fondo -> output in step6_variazioni_3d_coerenti.html /
     variazioni_frequenza_coerenti.npz;
  5) i relativi grafici 3D (Plotly) e un ricontrollo automatico finale, che
     verifica ANCHE punto-per-punto la coerenza tetto/fondo (non solo che
     i conteggi globali coincidano).

*** PERCHE' "ACUSTICA" NON E' UNA LICENZA POETICA A CASO ***
Il metodo Biondi-Malanga (skill biondi-malanga-sar-tomography /
sar-doppler-tomography) NON pretende che l'onda radar (elettromagnetica)
attraversi la piramide: usa il micro-Doppler indotto dalla MICRO-VIBRAZIONE
della struttura (vento, traffico, microsismi) come se il SAR fosse un
RICEVITORE di onde acustiche/sismiche, non un illuminatore penetrante. Nel
codice della skill (sar_tomo.py: GeometriaSAR.lambda_sonic, v_sonic, f_sonic)
questa e' gia' l'interpretazione "sonica" ufficiale del paper. Questo script
resta fedele a quella stessa interpretazione, solo:
  - ristretto al canale VH (polarizzazione incrociata: piu' sensibile alla
    dispersione di volume/struttura che alla riflessione speculare pura di VV);
  - con vocabolario esplicitamente "acustico" nei log/output.

*** COSA NON E' (disclaimer onesto, eredita quelli di piramide_unificato.py) ***
- Il dato di partenza (ampiezza/fase VH) resta backscatter radar, NON una
  registrazione acustica reale: la reinterpretazione "onda/vibrazione" e'
  un'ANALOGIA MATEMATICA (stessi coefficienti di Fourier, stessa inversione
  h(z)=A^H*Y), non una misura fisica di suono.
- Con lo stack Sentinel-1 IW aperto i "look" tomografici sono pochi ->
  altezza di ambiguita' z_amb piccola (skill sar-doppler-tomography):
  profondita' oltre z_amb vanno in ALIAS, la profondita' assoluta a scala
  piramide NON e' recuperabile da questo stack.
- Usare i risultati per ipotesi esplorative, MAI come prova di strutture
  interne non documentate.

*** RIUSO DI CODICE (nessuna duplicazione della pipeline) ***
Importa direttamente da piramide_unificato.py (stesso repo) gli step
3 (estrazione box) / 4 (array a N strati) / 5 (grafico onde 3D) / 6
(variazioni PIENO/VUOTO 3D), forzando la polarizzazione a 'vh'. Per la
tomografia Doppler vera scrive una lettura-stack dedicata al VH (la skill
scripts/tomographic_images.py e' hardcoded su '-vv-': qui NON viene
modificata, si duplica solo il ritaglio del box + il loop di inversione,
usando le funzioni di libreria GENERICHE di
skills/sar-doppler-tomography/scripts/sar_tomo.py, che non sanno nulla di
polarizzazione).

Esempi:
  python piramide_acustica_vh.py --pyramid kefren
  python piramide_acustica_vh.py --pyramid cheope --steps 3-6 --skip-vera-tomografia
"""
import os, sys, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import piramide_unificato as pu  # noqa: E402  (import dopo sys.path.insert, voluto)

_SAR_TOMO_DIR = os.path.join(HERE, "skills", "sar-doppler-tomography", "scripts")
sys.path.insert(0, _SAR_TOMO_DIR)
import sar_tomo as T  # noqa: E402


def log(*a):
    print("[acustica-vh]", *a)


# Passo di profondita' della forma d'onda anti-FT (step 5-6): l'onda generata
# dalla trasformata di Fourier viene campionata ogni DZ_ONDA_DEFAULT metri
# (default 0.33 m). NZ (numero di campioni) viene derivato da questo passo e
# dalla finestra z_amb del dato sorgente (NZ = z_amb/dz + 1), non piu' un
# numero di campioni fisso: sovrascrivibile con --dz-onda o direttamente con
# --nz (che ha precedenza).
DZ_ONDA_DEFAULT = 0.33

# Soglia PIENO/VUOTO di default per lo step 6: MEDIA della frequenza
# istantanea su TUTTI i campioni della forma d'onda (ogni dz, non solo il
# picco per pixel): campione con frequenza SOPRA la media = VUOTO, SOTTO la
# media = PIENO. Sovrascrivibile con --soglia esplicita.


# ===================================================================== self-test
def self_test_tomografia():
    """Valida la catena di inversione h(z)=A^H*Y su riflettori SINTETICI a
    profondita' note (12/41/78 m) PRIMA di applicarla al canale VH reale:
    se l'inversione non recupera i riflettori iniettati, non ha senso
    interpretarne l'output come "acustica" della piramide (vedi CLAUDE.md:
    verificare le modifiche alla catena di tomografia col forward model
    sintetico)."""
    from scipy.signal import find_peaks
    geo = T.GeometriaSAR()
    k = 48
    Y = T.genera_Y_sintetico(geo, k, snr_db=20.0)
    z = np.linspace(0, 100, 400)
    h = np.abs(T.tomografia_beamforming(T.steering_matrix(geo, z, k), Y))
    pk, _ = find_peaks(h / h.max(), height=0.3, distance=8)
    trovati = np.round(z[pk], 1).tolist()
    attesi = [12.0, 41.0, 78.0]
    ok = len(pk) >= 3 and all(any(abs(t - a) <= 2.0 for t in trovati) for a in attesi)
    log(f"self-test inversione tomografica: attesi {attesi} m, trovati {trovati} m "
        f"-> {'OK' if ok else 'FALLITO'}")
    if not ok:
        sys.exit("[acustica-vh] self-test della tomografia FALLITO: la catena "
                  "h(z)=A^H*Y non recupera i riflettori sintetici iniettati. "
                  "Interrompo (--skip-self-test per bypassare, sconsigliato).")
    return True


# ===================================================================== tomografia VH
def tomografia_vh(stackdir, aoi_dms, outdir, zmax=8.5, klook=12, ovs=4.0):
    """Tomografia Doppler VERA (Biondi & Malanga) SOLO sul canale VH: stessa
    catena di skills/sar-doppler-tomography/scripts/tomographic_images.py
    (che legge '-vv-'), riscritta qui per leggere '-vh-' SENZA toccare lo
    script della skill. Ritorna un dict coi risultati oppure False/None se
    non ci sono TIFF/annotation VH utilizzabili."""
    import glob, re
    import rasterio
    from rasterio.windows import Window
    from rasterio.transform import GCPTransformer
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nw, nw_lon, se, se_lon = aoi_dms
    lat_nw = pu.dms2dec(*nw[:3], nw[3]); lon_nw = pu.dms2dec(*nw_lon[:3], nw_lon[3])
    lat_se = pu.dms2dec(*se[:3], se[3]); lon_se = pu.dms2dec(*se_lon[:3], se_lon[3])
    latN, latS = max(lat_nw, lat_se), min(lat_nw, lat_se)
    lonW, lonE = min(lon_nw, lon_se), max(lon_nw, lon_se)
    corners = [(lonW, latN), (lonE, latN), (lonE, latS), (lonW, latS)]

    tiffs = sorted(glob.glob(os.path.join(stackdir, "*-iw[1-3]-slc-vh-*.tiff"))) or \
            sorted(glob.glob(os.path.join(stackdir, "*slc-vh-*.tiff")))
    if not tiffs:
        log(f"nessun TIFF VH in {stackdir}: salto la tomografia Doppler vera")
        return None

    def leggi_complesso(ds, win):
        a = ds.read(1, window=win)
        if np.iscomplexobj(a):
            return a.astype(np.complex64)
        if a.dtype.names and len(a.dtype.names) == 2:
            n0, n1 = a.dtype.names
            return a[n0].astype(np.float32) + 1j * a[n1].astype(np.float32)
        return a.astype(np.complex64)

    def annotation_vh(tif):
        m = re.search(r"(s1[abcd])-iw(\d)-slc-vh-(\d{8})", os.path.basename(tif).lower())
        if not m:
            return None
        s, sw, dt = m.groups()
        for x in glob.glob(os.path.join(stackdir, "*.annotation.xml")):
            xb = os.path.basename(x).lower()
            if s in xb and f"iw{sw}" in xb and dt in xb:
                return x
        return None

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
            chip = leggi_complesso(ds, Window(c0, r0, c1 - c0, r1 - r0))
        # un sub-swath (es. iw1/iw3) puo' passare il controllo sui bounds del
        # raster pur avendo l'AOI fuori dal proprio swath illuminato: il chip
        # e' quasi tutto zeri (stesso controllo di _estrai_chips_pol nella
        # pipeline principale). Scartarlo qui e' NECESSARIO non solo per
        # pulizia del log ma per correttezza: se il PRIMO chip valido preso
        # come riferimento geometrico ('ref') fosse uno di questi, l'intera
        # tomografia userebbe l'incidenza/range del sub-swath SBAGLIATO.
        valid = float(np.mean(np.abs(chip) > 0) * 100)
        if valid < pu.MIN_VALID_PCT:
            log(f"  {os.path.basename(t)}: copertura {valid:.0f}% < "
                f"{pu.MIN_VALID_PCT:.0f}% (AOI fuori da questo sub-swath), skip")
            continue
        chips.append(chip)
        if ref is None:
            ref = {"r0": r0, "c0": c0, "tif": t}
    if not chips:
        log("nessuna scena VH copre il box: salto la tomografia Doppler vera")
        return None
    H = min(c.shape[0] for c in chips); W = min(c.shape[1] for c in chips)
    log(f"tomografia VH: {len(chips)} scene | box {W}(range) x {H}(azimuth) px")

    ann = annotation_vh(ref["tif"])
    if ann is None:
        log("annotation.xml VH non trovato: salto la tomografia Doppler vera")
        return None
    g = T.leggi_geometria(ann)
    geo = T.GeometriaSAR()
    geo.f_em = g["f_em"]; geo.V = g["V"]; geo.theta_deg = g["incid_mid"]
    geo.R0 = g["R_near"] + (ref["c0"] + W / 2.0) * g["dR"]; geo.B_doppler = 1.0 / g["dt_az"]
    # apertura sintetica REALE del dato sorgente (non il default 84 km del paper,
    # che e' il COSMO-SkyMed Spotlight della loro Table 1): stessa identica
    # formula/branching di piramide_unificato._profondita_massima_fourier, cosi'
    # la tomografia Doppler vera usa la STESSA apertura dell'analisi esplorativa
    # di Fourier invece di ereditare in silenzio quella, molto piu' fine, del
    # paper -- vedi sar_tomo.apertura_sintetica_em
    if geo.f_em > 8e9:                              # X-band: COSMO-SkyMed
        pass                                         # apertura di default (84 km) gia' corretta
    else:                                            # C-band: Sentinel-1
        iw_src = "iw" in os.path.basename(ref["tif"]).lower()
        res_az = 22.0 if iw_src else 5.0
        geo.apertura_orbitale = T.apertura_sintetica_em(geo.f_em, geo.R0, res_az)
        log(f"apertura sintetica reale (Sentinel-1 {'IW' if iw_src else 'SM'}, "
            f"C-band): A={geo.apertura_orbitale/1000:.2f} km (non il default 84 km "
            f"COSMO-SkyMed Spotlight del paper)")
    dz = geo.risoluzione_tomografica() / ovs
    z = np.arange(0.0, zmax + dz, dz)
    A = T.steering_matrix(geo, z, klook)
    z_amb = T.altezza_ambiguita(geo, klook)
    log(f"tomografia VH: z 0..{zmax} m, {len(z)} campioni (dz={dz:.2f} m) | "
        f"z_amb={z_amb:.1f} m (oltre = alias, coerente con la skill "
        f"sar-doppler-tomography)")
    if zmax > z_amb * 1.01:
        log(f"ATTENZIONE: zmax > z_amb: profondita' oltre {z_amb:.1f} m in ALIAS")

    V = np.zeros((H, W, len(z))); V2 = np.zeros((H, W, len(z))); n = np.zeros((H, W))
    for s, chip in enumerate(chips):
        chip = chip[:H, :W]
        subs = [None if np.abs(chip[:, c]).max() == 0
                else T.sotto_aperture_doppler(chip[:, c], n_sub=klook + 1) for c in range(W)]
        for c in range(W):
            if subs[c] is None:
                continue
            for i in range(H):
                Y = np.array([T.pixel_tracking(subs[c][m], subs[c][m + 1], i)
                              for m in range(klook)], complex)
                if np.abs(Y).max() == 0:
                    continue
                h = np.abs(A.conj().T @ Y)
                V[i, c] += h; V2[i, c] += h ** 2; n[i, c] += 1
        log(f"  scena VH {s + 1}/{len(chips)} elaborata")
    n3 = np.maximum(n[..., None], 1)
    V /= n3
    Vstd = np.sqrt(np.maximum(V2 / n3 - V ** 2, 0.0))
    Vn = V / (V.max() + 1e-12)

    tomo_out = os.path.join(outdir, "tomografia_vera_vh")
    os.makedirs(tomo_out, exist_ok=True)
    np.save(os.path.join(tomo_out, "volume_acustico_vh.npy"), V)
    np.save(os.path.join(tomo_out, "volume_acustico_vh_std.npy"), Vstd)

    dA = g["dA"]; gr = g["dR"] / np.sin(np.radians(g["incid_mid"]))
    ew = np.arange(W) * gr; ns = np.arange(H) * dA
    db = lambda x: 10 * np.log10(np.maximum(x, 1e-6))

    ir = H // 2
    fig, ax = plt.subplots(figsize=(9, 4))
    im = ax.imshow(db(Vn[ir].T), origin="lower", aspect="auto", cmap="inferno",
                   extent=[0, ew[-1], 0, zmax])
    ax.set_xlabel("ground range [m]"); ax.set_ylabel("profondita' z [m]")
    ax.set_title(f"Tomografia acustica VH (range-profondita') — riga azimuth {ir}")
    fig.colorbar(im, ax=ax, label="riflettivita' [dB]"); fig.tight_layout()
    fig.savefig(os.path.join(tomo_out, "sezione_range_quota_vh.png"), dpi=150); plt.close(fig)

    ipk = np.argmax(V, axis=2)
    quota = z[ipk]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    im0 = axs[0].imshow(quota, origin="upper", aspect="auto", cmap="turbo",
                        extent=[0, ew[-1], ns[-1], 0])
    axs[0].set_title("Quota del riflettore acustico dominante (max |h(z)|)")
    axs[0].set_xlabel("range [m]"); axs[0].set_ylabel("azimuth [m]")
    fig.colorbar(im0, ax=axs[0], label="profondita' [m]")
    im1 = axs[1].imshow(db(Vn.max(2)), origin="upper", aspect="auto", cmap="inferno",
                        extent=[0, ew[-1], ns[-1], 0])
    axs[1].set_title("Intensita' di picco"); axs[1].set_xlabel("range [m]")
    fig.colorbar(im1, ax=axs[1], label="[dB]")
    fig.tight_layout(); fig.savefig(os.path.join(tomo_out, "mappa_quota_vh.png"), dpi=150)
    plt.close(fig)

    log(f"tomografia acustica VH salvata in {tomo_out}/ "
        f"(sezione_range_quota_vh.png, mappa_quota_vh.png, volume_acustico_vh.npy)")
    return dict(V=V, Vstd=Vstd, z=z, z_amb=z_amb, n_scene=len(chips))


# ===================================================================== array 3D grezzo
def grafico_array_3d_vh(disp, x, y, etich, outdir):
    """Grafico 3D Plotly dell'array RAW a N strati dello step 4 (TUTTE le
    coppie interferometriche VH disponibili, non solo quelle sequenziali:
    vedi n_layers in main()): un PIANO ORIZZONTALE per ogni strato, colorato
    dallo spostamento verticale DInSAR [mm] (stessa scala RdBu_r di
    step4_strati.png). A differenza dello step 5 (che reinterpreta gli N
    strati come coefficienti di Fourier di una forma d'onda sviluppata in
    profondita' via anti-trasformata), qui i piani restano il DATO GREZZO
    cosi' come esce dalla DInSAR classica (fase interferometrica -> fase/
    cos(incidenza) -> spostamento verticale): nessuna trasformata, nessuna
    reinterpretazione.
    L'asse verticale e' la BASELINE TEMPORALE REALE [giorni] tra le due
    acquisizioni di ogni coppia (letta dall'etichetta "AAAAMMGG~AAAAMMGG
    (Nd)" prodotta da _interferogrammi_pol), NON un indice arbitrario di
    posizione nella lista: coppie con baseline diverse (specie ora che si
    usano TUTTE le combinazioni N*(N-1)/2 e non solo le N-1 sequenziali)
    hanno un significato fisico diverso (piu' baseline = piu' tempo per la
    decorrelazione/il rumore, meno baseline = misura piu' pulita ma piu'
    ravvicinata) -- l'asse deve rifletterlo, non un ordine di lista."""
    import re
    import plotly.graph_objects as go
    n_layers, ny, nx = disp.shape
    Xg, Yg = np.meshgrid(x, y)
    vmax = float(np.percentile(np.abs(disp), 98)) or 1.0
    baseline_gg = []
    for e in etich:
        m = re.search(r"\((\d+)d\)", str(e))
        baseline_gg.append(int(m.group(1)) if m else None)
    if any(b is None for b in baseline_gg):
        log("ATTENZIONE: baseline temporale non estraibile da alcune etichette "
            "(formato inatteso), uso l'indice di strato come fallback per quei piani")
        baseline_gg = [b if b is not None else k for k, b in enumerate(baseline_gg)]
    traces = []
    for k in range(n_layers):
        bt = baseline_gg[k]
        traces.append(go.Surface(
            x=Xg, y=Yg, z=np.full_like(Xg, float(bt)),
            surfacecolor=disp[k], cmin=-vmax, cmax=vmax, colorscale="RdBu_r",
            showscale=(k == 0),
            colorbar=dict(title="spost.<br>verticale<br>[mm]") if k == 0 else None,
            name=str(etich[k]), hovertext=f"{etich[k]} (VH)", hoverinfo="text",
            opacity=0.97))
    fig = go.Figure(traces)
    pu._tema_plotly(fig)
    ratio_y = (float(y.max()) / float(x.max())) if x.size and x.max() else 1.0
    bt_span = max(max(baseline_gg) - min(baseline_gg), 1)
    fig.update_layout(scene=dict(
        xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
        zaxis_title="baseline temporale della coppia [giorni]",
        aspectmode="manual",
        aspectratio=dict(x=1.0, y=ratio_y, z=max(0.4, min(2.0, bt_span / 60.0)))))
    path = os.path.join(outdir, "step4_array_3d_vh.html")
    pu._scrivi_pagina_html(
        fig, path,
        titolo="Step 4 — array 3D grezzo a N strati (VH, tutte le coppie)",
        sottotitolo="Ogni piano e' un interferogramma DInSAR tra due acquisizioni "
                    "VH (spostamento verticale, DATO GREZZO: nessuna trasformata "
                    "di Fourier, a differenza dello step 5). L'array usa TUTTE le "
                    "coppie interferometriche possibili tra le scene VH disponibili "
                    "nello stack (non solo quelle sequenziali), per il massimo "
                    "numero di strati; l'asse verticale e' la baseline temporale "
                    "REALE [giorni] di ciascuna coppia, non un indice di lista. "
                    "Trascina per ruotare, rotella per lo zoom.",
        badges=[("strati (VH, tutte le coppie possibili)", str(n_layers)),
                ("pixel per strato (y-az x x-rng)", f"{ny}x{nx}"),
                ("baseline min-max", f"{min(baseline_gg)}-{max(baseline_gg)} gg"),
                ("spost. verticale max (98° perc.)", f"±{vmax:.1f} mm")],
        nota="Diverso dallo step 5 (forme d'onda 3D): qui NON c'e' alcuna "
             "reinterpretazione Fourier/onda -- sono gli interferogrammi DInSAR "
             "cosi' come escono dalla pipeline (fase -> spostamento verticale), "
             "un piano per coppia. L'asse verticale e' la baseline temporale "
             "reale (giorni tra le due date), NON una quota o profondita' fisica: "
             "coppie a baseline lunga sono piu' soggette a decorrelazione/rumore "
             "di quelle a baseline corta, e la posizione sull'asse lo rende "
             "visibile invece di appiattire tutto su un indice arbitrario.")
    log(f"array 3D grezzo (VH, {n_layers} strati/coppie, baseline "
        f"{min(baseline_gg)}-{max(baseline_gg)} gg) salvato in {path}")
    return path


# ===================================================================== coerenza PIENO/VUOTO
def variazioni_coerenti_vh(freq_pol, z, x_m, y_m, z0, outdir, soglia, isteresi,
                           max_variazioni, riferimenti, ZBOT, ZTOP, Lx, Ly):
    """Ricostruisce le variazioni PIENO<->VUOTO con lo STESSO trigger di
    Schmitt di step6_variazioni (piramide_unificato.py) ma con una selezione
    finale che garantisce la COERENZA per pixel: un "tetto" (PIENO->VUOTO)
    disegnato ha SEMPRE il proprio "fondo" (VUOTO->PIENO) nello stesso
    pixel, e viceversa. step6_variazioni seleziona indipendentemente le N
    variazioni piu' forti per TIPO su tutto il volume: verificato sui dati
    VH reali, questo lascia il 77% delle colonne con un punto orfano (tetto
    senza fondo o fondo senza tetto) -> incoerente. Qui si selezionano SOLO
    coppie complete (stesso pixel), ordinate per la forza combinata della
    coppia."""
    import plotly.graph_objects as go

    pols = list(freq_pol)
    nz = z.size
    z_mid = (z[:-1] + z[1:]) / 2.0 if len(z) > 1 else z
    tutti = np.concatenate([freq_pol[p].ravel() for p in pols])
    q0 = float((tutti <= soglia).mean())
    s_basso = float(np.quantile(tutti, max(0.0, q0 - isteresi / 2)))
    s_alto = float(np.quantile(tutti, min(1.0, q0 + isteresi / 2)))

    coppie = []   # (salto_coppia, pol, i, j, k_tetto, k_fondo)
    aperti = 0    # eventi SENZA partner nella finestra di profondita' (bordo)
    for p in pols:
        fr = freq_pol[p]
        stato_pieno = fr[:, :, 0] < soglia
        pendente = {}   # (i,j) -> (k_tetto, salto_tetto) in attesa del fondo
        for k in range(1, nz):
            diventa_vuoto = stato_pieno & (fr[:, :, k] >= s_alto)
            diventa_pieno = (~stato_pieno) & (fr[:, :, k] <= s_basso)
            salto_k = np.abs(fr[:, :, k] - fr[:, :, k - 1])
            for i, j in zip(*np.nonzero(diventa_vuoto)):
                pendente[(i, j)] = (k, float(salto_k[i, j]))
            for i, j in zip(*np.nonzero(diventa_pieno)):
                if (i, j) in pendente:
                    k_t, salto_t = pendente.pop((i, j))
                    coppie.append((max(salto_t, float(salto_k[i, j])), p, i, j, k_t, k))
                else:
                    aperti += 1   # fondo senza tetto precedente (colonna gia' VUOTO al bordo z=0)
            stato_pieno = np.where(diventa_vuoto, False,
                                   np.where(diventa_pieno, True, stato_pieno))
        aperti += len(pendente)   # tetti il cui vuoto arriva al fondo della finestra, senza fondo

    log(f"coerenza PIENO/VUOTO: ricostruite {len(coppie)} coppie COMPLETE "
        f"tetto+fondo (stesso pixel) + {aperti} eventi singoli SCARTATI "
        f"(vuoto che tocca il bordo superiore/inferiore della finestra di "
        f"profondita', nessun partner disponibile)")

    max_coppie = max(max_variazioni // 2, 1)
    coppie.sort(key=lambda c: -c[0])
    coppie = coppie[:max_coppie]

    tx_, ty_, tz_, tv_, tpol_, ttipo_, salto_ = [], [], [], [], [], [], []
    for salto, p, i, j, k_t, k_f in coppie:
        fr = freq_pol[p]
        for k, tipo in ((k_t, "tetto"), (k_f, "fondo")):
            kv = k if tipo == "tetto" else k - 1
            tx_.append(x_m[j]); ty_.append(y_m[i]); tz_.append(z0[i, j] - z_mid[k - 1])
            tv_.append(float(fr[i, j, kv])); tpol_.append(p); ttipo_.append(tipo)
            salto_.append(salto)
    tx = np.array(tx_); ty = np.array(ty_); tz = np.array(tz_); tv = np.array(tv_)
    tpol = np.array(tpol_, dtype=object); ttipo = np.array(ttipo_, dtype=object)
    salto = np.array(salto_)

    n_tetti = int((ttipo == "tetto").sum()); n_fondi = int((ttipo == "fondo").sum())
    if len(coppie):
        assert n_tetti == n_fondi == len(coppie), \
            "coerenza PIENO/VUOTO rotta: numero tetti != numero fondi"
    log(f"coerenza PIENO/VUOTO: selezione COERENTE = {len(coppie)} coppie -> "
        f"{tx.size} punti ({n_tetti} tetti + {n_fondi} fondi, per costruzione "
        f"SEMPRE in numero uguale e appaiati sullo stesso pixel)")

    np.savez_compressed(os.path.join(outdir, "variazioni_frequenza_coerenti.npz"),
                        trans_x=tx, trans_y=ty, trans_z=tz, trans_freq_vuoto=tv,
                        trans_salto=salto, trans_pol=tpol.astype(str),
                        trans_tipo=ttipo.astype(str), soglia=soglia, isteresi=isteresi,
                        n_coppie=len(coppie), n_scartati=aperti)

    COL_TETTO, COL_FONDO = "#0d0d12", "#b45309"
    Xg, Yg = np.meshgrid(x_m, y_m)
    ground = go.Scatter3d(x=Xg.ravel(), y=Yg.ravel(), z=z0.ravel(), mode="markers",
                          marker=dict(size=1.2, color="red"), name="suolo (DEM)")
    zero_plane = go.Surface(x=Xg, y=Yg, z=np.zeros_like(Xg), showscale=False,
                            opacity=0.25, colorscale=[[0, "gray"], [1, "gray"]],
                            name="livello 0 (s.l.m.)")
    extra_traces = pu._tracce_riferimenti(riferimenti) if riferimenti else []
    # NIENTE linee verticali tetto->fondo: appesantivano il grafico e a densita'
    # alta formavano una "tenda" grigia che nascondeva i punti. L'appaiamento
    # resta garantito dalla selezione (solo coppie complete per pixel) ed e'
    # verificabile in variazioni_frequenza_coerenti.npz; nel grafico i due punti
    # della coppia condividono la stessa verticale (stesso x,y).
    tracce = []
    for p in pols:
        for tipo, colore, etich in (("tetto", COL_TETTO, "tetto (PIENO->VUOTO)"),
                                    ("fondo", COL_FONDO, "fondo (VUOTO->PIENO)")):
            m = (tpol == p) & (ttipo == tipo)
            tracce.append(go.Scatter3d(
                x=tx[m], y=ty[m], z=tz[m], mode="markers",
                name=f"{etich} {p.upper()} ({int(m.sum())} punti, coppie coerenti)",
                marker=dict(size=2, color=colore, opacity=0.95)))
    fig = go.Figure([*tracce, zero_plane, ground, *extra_traces])
    pu._tema_plotly(fig)
    fig.update_layout(
        scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
                   zaxis_title="quota s.l.m. [m]", zaxis=dict(range=[ZBOT, ZTOP]),
                   aspectmode="manual", aspectratio=dict(x=1.0, y=Ly / Lx, z=(ZTOP - ZBOT) / Lx)))
    pu._scrivi_pagina_html(
        fig, os.path.join(outdir, "step6_variazioni_3d_coerenti.html"),
        titolo="Step 6 (coerente) — variazioni PIENO/VUOTO accoppiate VH",
        sottotitolo="Selezione COERENTE per pixel: ogni tetto disegnato ha il proprio "
                    "fondo sulla stessa verticale (stesso pixel; nessuna linea di "
                    "collegamento, l'appaiamento e' nel .npz). Trascina per ruotare, "
                    "rotella per lo zoom.",
        badges=[("coppie tetto+fondo", str(len(coppie))),
                ("punti totali", str(tx.size)),
                ("eventi scartati (bordo finestra, senza partner)", str(aperti))],
        nota="A differenza dello step 6 standard (che seleziona indipendentemente "
             "le N variazioni piu' forti per TIPO su tutto il volume, e puo' "
             "lasciare un tetto senza il suo fondo o viceversa: verificato sui "
             "dati VH reali, il 77% delle colonne selezionate risultava orfano), "
             "qui si selezionano SOLO coppie complete tetto+fondo dello STESSO "
             "pixel: ogni vuoto disegnato ha sempre sia il tetto che il fondo.")
    return dict(n_coppie=len(coppie), n_scartati=aperti, n_punti=int(tx.size))


# ===================================================================== PIENO/VUOTO alternativo per-scena
def pieno_vuoto_comune_vh(boxpath, x_m, y_m, z0, Lx, Ly, ZTOP, ZBOT, outdir,
                          riferimenti, freq_pol=None, z=None, soglia_freq=None,
                          soglia_ampiezza=None):
    """Metodo ALTERNATIVO di classificazione PIENO/VUOTO, indipendente dalla
    forma d'onda di Fourier degli step 5-6: per OGNI scena VH grezza (le
    n_scene acquisizioni ritagliate in box_vh.npz, non le coppie
    interferometriche) classifica ogni pixel PIENO (ampiezza di backscatter
    >= soglia, riflettore forte) o VUOTO (ampiezza < soglia, riflettore
    debole), con UNA SOGLIA GLOBALE (mediana dell'ampiezza raggruppata su
    TUTTE le scene insieme, cosi' il confronto tra scene e' fisicamente
    comparabile: una soglia per-scena forzerebbe artificialmente ~lo stesso
    split ogni volta, mascherando le differenze reali). Poi confronta le
    n_scene mappe booleane PIXEL PER PIXEL: tiene SOLO i pixel dove TUTTE le
    scene sono concordi (sempre PIENO o sempre VUOTO) -> punti "in comune";
    scarta i pixel che cambiano classificazione da una scena all'altra
    (instabilita' temporale, rumore): questo da' un insieme di colonne (x,y)
    "affidabili" -- ma da SOLO l'ampiezza (una scena = un'immagine 2D) non
    esiste alcuna profondita' propria, quindi di per se' resterebbe un punto
    per colonna, sempre in superficie.
    Se vengono passati freq_pol/z/soglia_freq (la stessa forma d'onda
    Hilbert gia' calcolata e validata dagli step 5-6, sviluppata in
    profondita' da anti-FT delle coppie interferometriche), le colonne
    affidabili vengono ESTESE in profondita': per ogni colonna concorde
    PIENO (ampiezza) si tengono i singoli campioni z della sua forma d'onda
    ANCH'ESSI classificati PIENO dalla frequenza istantanea (idem VUOTO),
    cioe' i punti dove le DUE classificazioni INDIPENDENTI (ampiezza
    per-scena grezza e frequenza Hilbert della forma d'onda combinata)
    CONCORDANO, non solo in superficie ma a ogni profondita' sviluppata.
    Senza questi argomenti la funzione ricade sul solo punto di superficie
    per colonna (comportamento precedente).
    NOTA METODOLOGICA: l'ampiezza di backscatter radar grezzo (singola
    acquisizione) e la frequenza istantanea di Hilbert (forma d'onda
    sviluppata in profondita') sono due segnali DIVERSI derivati dagli stessi
    dati: qui si usa il PRIMO come filtro di affidabilita' sulla colonna e il
    SECONDO per la profondita' dei punti. Analisi ESPLORATIVA alternativa,
    non sostituisce e non e' piu' affidabile della tomografia Doppler vera
    del paper Biondi-Malanga."""
    import plotly.graph_objects as go
    d = np.load(boxpath, allow_pickle=True)
    if "vh" not in d.files:
        log("pieno/vuoto per-scena: box senza 'vh', salto")
        return None
    vh = d["vh"]  # (n_scene, ny, nx) complesso
    nomi = [str(n) for n in d["vh_scene"]]
    amp = np.abs(vh)
    n_scene, ny, nx = amp.shape
    if soglia_ampiezza is None:
        soglia_ampiezza = float(np.median(amp))
    pieno_s = amp >= soglia_ampiezza          # (n_scene, ny, nx) booleano
    vuoto_s = ~pieno_s
    pieno_comune = pieno_s.all(axis=0)        # concorde PIENO in TUTTE le scene
    vuoto_comune = vuoto_s.all(axis=0)        # concorde VUOTO in TUTTE le scene
    n_p = int(pieno_comune.sum()); n_v = int(vuoto_comune.sum())
    n_tot = ny * nx
    log(f"pieno/vuoto per-scena VH: {n_scene} scene ({', '.join(nomi[:3])}"
        f"{'...' if n_scene > 3 else ''}), soglia ampiezza globale (mediana "
        f"pooled su tutte le scene) = {soglia_ampiezza:.4g}")
    log(f"confronto pixel per pixel tra le {n_scene} mappe booleane -> "
        f"{n_p} pixel SEMPRE PIENO + {n_v} pixel SEMPRE VUOTO = "
        f"{n_p + n_v}/{n_tot} punti in comune ({100 * (n_p + n_v) / n_tot:.0f}%), "
        f"{n_tot - n_p - n_v} pixel scartati (classificazione non concorde "
        f"tra le scene)")

    # estensione in profondita': se disponibile la forma d'onda Hilbert
    # (freq_pol/z/soglia_freq, la stessa di step 5-6), per ogni colonna
    # concorde in ampiezza si tengono i campioni z ANCHE concordi in
    # frequenza istantanea -> punti non solo in superficie
    con_profondita = (freq_pol is not None and z is not None
                      and soglia_freq is not None and "vh" in freq_pol)
    if con_profondita:
        fr = freq_pol["vh"]                 # (ny, nx, nz)
        pieno_voxel = fr < soglia_freq       # stessa convenzione di step6/coerenza
        vuoto_voxel = ~pieno_voxel
        mask_p3 = pieno_voxel & pieno_comune[:, :, None]
        mask_v3 = vuoto_voxel & vuoto_comune[:, :, None]
        ip, jp, kp = np.nonzero(mask_p3)
        iv, jv, kv = np.nonzero(mask_v3)
        px = x_m[jp]; py = y_m[ip]; pz = z0[ip, jp] - z[kp]
        vx = x_m[jv]; vy = y_m[iv]; vz = z0[iv, jv] - z[kv]
        log(f"estensione in profondita' (frequenza Hilbert, soglia {soglia_freq:.4g} "
            f"1/m, stessa forma d'onda di step 5-6): sulle {n_p} colonne PIENO "
            f"in comune -> {px.size} punti anche PIENO in frequenza a varie "
            f"profondita'; sulle {n_v} colonne VUOTO in comune -> {vx.size} "
            f"punti anche VUOTO in frequenza")
    else:
        Xg0, Yg0 = np.meshgrid(x_m, y_m)
        px = Xg0[pieno_comune]; py = Yg0[pieno_comune]; pz = z0[pieno_comune]
        vx = Xg0[vuoto_comune]; vy = Yg0[vuoto_comune]; vz = z0[vuoto_comune]

    np.savez_compressed(os.path.join(outdir, "pieno_vuoto_comune_vh.npz"),
                        pieno_comune=pieno_comune, vuoto_comune=vuoto_comune,
                        pieno_per_scena=pieno_s, soglia_ampiezza=soglia_ampiezza,
                        vh_scene=np.array(nomi), x=x_m, y=y_m, z0=z0,
                        pieno_x=px, pieno_y=py, pieno_z=pz,
                        vuoto_x=vx, vuoto_y=vy, vuoto_z=vz,
                        con_profondita=con_profondita,
                        soglia_freq=(soglia_freq if soglia_freq is not None else np.nan))

    COL_PIENO, COL_VUOTO = "#0d0d12", "#b45309"
    Xg, Yg = np.meshgrid(x_m, y_m)
    tracce = [
        go.Scatter3d(x=px, y=py, z=pz, mode="markers",
                     name=f"PIENO (comune a tutte le scene): {px.size} punti",
                     marker=dict(size=2.2, color=COL_PIENO, opacity=0.9)),
        go.Scatter3d(x=vx, y=vy, z=vz, mode="markers",
                     name=f"VUOTO (comune a tutte le scene): {vx.size} punti",
                     marker=dict(size=2.2, color=COL_VUOTO, opacity=0.9)),
    ]
    zero_plane = go.Surface(x=Xg, y=Yg, z=np.zeros_like(Xg), showscale=False,
                            opacity=0.2, colorscale=[[0, "gray"], [1, "gray"]],
                            name="livello 0 (s.l.m.)")
    ground = go.Scatter3d(x=Xg.ravel(), y=Yg.ravel(), z=z0.ravel(), mode="markers",
                          marker=dict(size=1.2, color="red"), name="suolo (DEM)")
    extra_traces = pu._tracce_riferimenti(riferimenti) if riferimenti else []
    fig = go.Figure([*tracce, zero_plane, ground, *extra_traces])
    pu._tema_plotly(fig)
    zmin_plot = float(min(z0.min(), pz.min() if pz.size else z0.min(),
                          vz.min() if vz.size else z0.min())) - 10
    fig.update_layout(scene=dict(
        xaxis_title="x E-W [m]", yaxis_title="y N-S [m]", zaxis_title="quota s.l.m. [m]",
        zaxis=dict(range=[zmin_plot, ZTOP]),
        aspectmode="manual", aspectratio=dict(x=1.0, y=Ly / Lx, z=0.6)))
    path = os.path.join(outdir, "step_pieno_vuoto_comune_vh.html")
    sottotitolo = (
        f"Classificazione INDIPENDENTE per ognuna delle {n_scene} scene VH grezze "
        "(ampiezza di backscatter, soglia globale = mediana pooled su tutte le "
        "scene), confrontate pixel per pixel: le colonne dove tutte le scene "
        "concordano (sempre PIENO o sempre VUOTO) sono poi sviluppate in "
        "profondita' con la forma d'onda Hilbert di step 5-6, tenendo SOLO i "
        "campioni dove anche la frequenza istantanea concorda -- non solo in "
        "superficie. Drappeggio DEM incluso. Trascina per ruotare, rotella per lo zoom."
        if con_profondita else
        f"Classificazione INDIPENDENTE per ognuna delle {n_scene} scene VH grezze "
        "(ampiezza di backscatter, soglia globale = mediana pooled su tutte le "
        "scene), confrontate pixel per pixel: disegnati SOLO i pixel dove tutte "
        "le scene concordano (sempre PIENO o sempre VUOTO). Nessuna profondita' "
        "disponibile (freq_pol/z non passati): punti in superficie (DEM). "
        "Trascina per ruotare, rotella per lo zoom.")
    pu._scrivi_pagina_html(
        fig, path,
        titolo="Metodo alternativo — PIENO/VUOTO per-scena, punti in comune",
        sottotitolo=sottotitolo,
        badges=[("scene VH usate", str(n_scene)),
                ("soglia ampiezza (mediana pooled)", f"{soglia_ampiezza:.3g}"),
                ("colonne PIENO in comune (ampiezza)", str(n_p)),
                ("colonne VUOTO in comune (ampiezza)", str(n_v)),
                ("colonne scartate (non concordi)", str(n_tot - n_p - n_v))]
               + ([("punti PIENO (anche in profondita')", str(px.size)),
                   ("punti VUOTO (anche in profondita')", str(vx.size)),
                   ("soglia frequenza (step 5-6)", f"{soglia_freq:.4g} 1/m")]
                  if con_profondita else []),
        nota="Metodo ESPLORATIVO alternativo a step 5-6: la SELEZIONE delle colonne "
             "(quali x,y fidarsi) viene dalla sola ampiezza di backscatter radar per "
             "singola acquisizione (riflettore forte/debole, NON fase/frequenza); "
             "la PROFONDITA' dei punti disegnati, quando disponibile, viene invece "
             "dalla stessa frequenza istantanea di Hilbert usata da step 5-6 -- sono "
             "due segnali indipendenti derivati dagli stessi dati, qui incrociati: "
             "un punto e' disegnato solo se ENTRAMBI i metodi concordano sul tipo "
             "(pieno/vuoto), a quella colonna e a quella profondita'. Non e' "
             "tomografia Doppler del paper, ne' evidenza di strutture interne.")
    log(f"grafico 3D 'punti in comune'{' (con profondita)' if con_profondita else ''} "
        f"salvato in {path}")
    return dict(n_scene=n_scene, n_pieno=n_p, n_vuoto=n_v, soglia=soglia_ampiezza)


# ===================================================================== ricontrollo finale
def verifica_finale(outdir, tomo_info, coerenza_info=None, pv_comune_info=None):
    """Ricontrollo automatico di fine pipeline: nessuna eccezione non e'
    sufficiente, si verifica che gli array numerici salvati non contengano
    NaN/Inf e che i file attesi esistano davvero, prima di dichiarare
    l'esecuzione affidabile."""
    problemi = []

    array_path = os.path.join(outdir, "array3d_12strati.npz")
    if os.path.exists(array_path):
        d = np.load(array_path, allow_pickle=True)
        if "disp" not in d.files:
            problemi.append(f"{array_path}: manca la chiave 'disp'")
        else:
            disp = d["disp"]
            finito = bool(np.isfinite(disp).all())
            if not finito:
                problemi.append(f"{array_path}: 'disp' contiene NaN/Inf")
            log(f"verifica: array a {disp.shape[0]} strati VH, "
                f"spostamento verticale in [{disp.min():.1f}, {disp.max():.1f}] mm, "
                f"finito={finito}")
    else:
        problemi.append(f"{array_path}: mancante")

    onde_path = os.path.join(outdir, "forme_donda.npz")
    if os.path.exists(onde_path):
        d = np.load(onde_path, allow_pickle=True)
        if "w" not in d.files:
            problemi.append(f"{onde_path}: manca la chiave 'w'")
        else:
            w = d["w"]
            finito = bool(np.isfinite(w).all())
            if not finito:
                problemi.append(f"{onde_path}: la forma d'onda 'w' contiene NaN/Inf")
            log(f"verifica: forma d'onda acustica VH {w.shape}, "
                f"|w|max={np.abs(w).max():.2f}, finito={finito}")
    else:
        problemi.append(f"{onde_path}: mancante (step 5 non eseguito?)")

    html_path = os.path.join(outdir, "step5_forme_donda_3d.html")
    if not os.path.exists(html_path):
        problemi.append(f"{html_path}: mancante (grafico 3D onde non generato)")

    array3d_html = os.path.join(outdir, "step4_array_3d_vh.html")
    if not os.path.exists(array3d_html):
        problemi.append(f"{array3d_html}: mancante (grafico 3D array grezzo non generato)")

    if tomo_info:
        V = tomo_info["V"]
        finito = bool(np.isfinite(V).all())
        if not finito:
            problemi.append("volume tomografico VH: contiene NaN/Inf")
        log(f"verifica: volume tomografico VH {V.shape}, finito={finito}, "
            f"z_amb={tomo_info['z_amb']:.1f} m, scene usate={tomo_info['n_scene']}")

    if coerenza_info is not None:
        coer_path = os.path.join(outdir, "variazioni_frequenza_coerenti.npz")
        if not os.path.exists(coer_path):
            problemi.append(f"{coer_path}: mancante")
        else:
            d = np.load(coer_path, allow_pickle=True)
            tipo = d["trans_tipo"]
            n_t = int((tipo == "tetto").sum()); n_f = int((tipo == "fondo").sum())
            if n_t != n_f:
                problemi.append(f"{coer_path}: {n_t} tetti != {n_f} fondi (dovrebbero "
                                f"essere sempre uguali: ogni coppia e' un tetto+un fondo)")
            # verifica REALE punto per punto: per ogni pixel selezionato, ogni
            # tetto deve avere esattamente un fondo allo STESSO pixel/pol e la
            # quota del fondo dev'essere sotto (piu' profonda del) quella del
            # tetto -- non ci si fida del solo conteggio globale
            px = d["trans_x"]; py = d["trans_y"]; pz = d["trans_z"]; ppol = d["trans_pol"]
            gruppi = {}
            for xx, yy, zz, tt, pp in zip(px, py, pz, tipo, ppol):
                gruppi.setdefault((round(float(xx), 3), round(float(yy), 3), str(pp)), []).append(
                    (float(zz), str(tt)))
            orfani = 0
            for key, ev in gruppi.items():
                ev.sort(key=lambda e: -e[0])   # dalla superficie verso il basso
                tipi = [t for _, t in ev]
                if tipi.count("tetto") != tipi.count("fondo") or \
                        any(tipi[i] == tipi[i + 1] for i in range(len(tipi) - 1)):
                    orfani += 1
            if orfani:
                problemi.append(f"{coer_path}: {orfani} pixel con punti NON coerenti "
                                f"(tetto senza fondo o viceversa, o due tetti/fondi "
                                f"consecutivi sullo stesso pixel)")
            log(f"verifica: coerenza PIENO/VUOTO VH -> {n_t} tetti + {n_f} fondi, "
                f"{len(gruppi)} pixel coinvolti, {orfani} incoerenti (atteso 0), "
                f"{coerenza_info['n_scartati']} eventi di bordo scartati (attesi, "
                f"non un errore: vuoto che tocca il limite della finestra di profondita')")

    if pv_comune_info is not None:
        pv_path = os.path.join(outdir, "pieno_vuoto_comune_vh.npz")
        if not os.path.exists(pv_path):
            problemi.append(f"{pv_path}: mancante")
        else:
            d = np.load(pv_path, allow_pickle=True)
            pc = d["pieno_comune"]; vc = d["vuoto_comune"]
            sovrapposti = int((pc & vc).sum())
            if sovrapposti:
                problemi.append(f"{pv_path}: {sovrapposti} pixel marcati SIA "
                                f"pieno_comune SIA vuoto_comune (dovrebbero "
                                f"essere insiemi disgiunti per costruzione)")
            attesi_p, attesi_v = pv_comune_info["n_pieno"], pv_comune_info["n_vuoto"]
            if int(pc.sum()) != attesi_p or int(vc.sum()) != attesi_v:
                problemi.append(f"{pv_path}: conteggi salvati ({int(pc.sum())} "
                                f"pieno, {int(vc.sum())} vuoto) diversi da quelli "
                                f"riportati a video ({attesi_p} pieno, {attesi_v} vuoto)")
            log(f"verifica: pieno/vuoto per-scena VH -> {int(pc.sum())} pixel "
                f"PIENO comune + {int(vc.sum())} pixel VUOTO comune, "
                f"{sovrapposti} sovrapposti (atteso 0)")
        pv_html = os.path.join(outdir, "step_pieno_vuoto_comune_vh.html")
        if not os.path.exists(pv_html):
            problemi.append(f"{pv_html}: mancante (grafico 3D punti in comune non generato)")

    if problemi:
        log("VERIFICA FINALE: PROBLEMI TROVATI:")
        for p in problemi:
            log("  -", p)
        sys.exit("[acustica-vh] verifica finale FALLITA: vedi problemi sopra.")
    log("VERIFICA FINALE: OK, nessun errore/NaN/Inf trovato negli output prodotti.")


# ===================================================================== main
def main():
    ap = argparse.ArgumentParser(
        description="Tomografia 'acustica' SOLO-VH della piramide: reinterpreta "
                    "il canale VH come vibrazione riflessa dalla struttura interna.")
    ap.add_argument("--pyramid", choices=sorted(pu.PRESETS), default=None,
                    help="box di default per piramide (cheope|kefren); "
                         "sovrascrivibile con --nw/--nw-lon/--se/--se-lon")
    ap.add_argument("--nw", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--nw-lon", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--se", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--se-lon", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--stack", default=None,
                    help="cartella TIFF (default: stack_slc condiviso del preset)")
    ap.add_argument("--outdir", default=None,
                    help="cartella output (default: <preset>/acustica_vh)")
    ap.add_argument("--layers", type=int, default=None,
                    help="numero di strati VH; default None = usa TUTTE le scene VH "
                         "disponibili nello stack, ordinate per data (N scene -> N-1 "
                         "interferogrammi sequenziali)")
    ap.add_argument("--steps", default="3-6",
                    help="3=estrazione box VH, 4=array a N strati, 5=grafico onde 3D, "
                         "6=variazioni PIENO/VUOTO 3D")
    ap.add_argument("--no-track-filter", action="store_true")
    ap.add_argument("--no-pyramid", action="store_true", help="non sovrapporre la piramide al DEM")
    ap.add_argument("--nz", type=int, default=None,
                    help="numero di campioni della forma d'onda anti-FT; default None = "
                         "derivato da --dz-onda e dalla finestra z_amb del dato "
                         "(NZ = z_amb/dz_onda + 1). Se passato esplicitamente ha "
                         "precedenza su --dz-onda.")
    ap.add_argument("--dz-onda", type=float, default=DZ_ONDA_DEFAULT,
                    help=f"passo di profondita' [m] a cui campionare l'onda generata "
                         f"dalla Fourier (default {DZ_ONDA_DEFAULT} m); usato per "
                         "derivare --nz quando --nz non e' esplicito")
    ap.add_argument("--step-x", type=int, default=None)
    ap.add_argument("--step-y", type=int, default=1)
    ap.add_argument("--soglia", type=float, default=None,
                    help="soglia PIENO/VUOTO esplicita [1/m]; default None = MEDIA della "
                         "frequenza istantanea su tutti i campioni della forma d'onda "
                         "(ogni --dz-onda metri): frequenza sopra la media -> VUOTO, "
                         "sotto la media -> PIENO")
    ap.add_argument("--isteresi", type=float, default=0.10)
    ap.add_argument("--max-variazioni", type=int, default=1000)
    ap.add_argument("--zmax", type=float, default=200.0,
                    help="profondita' massima [m] per la tomografia vera. NOTA: con lo "
                         "stack Sentinel-1 IW aperto e l'apertura sintetica REALE del dato "
                         "(~1.1 km/burst, non gli 84 km COSMO-SkyMed Spotlight del paper) "
                         "z_amb e' tipicamente ~660 m con klook=12 ma la risoluzione per "
                         "cella e' grossolana (~24 m con ovs=4): oltre z_amb i valori vanno "
                         "in ALIAS (si ripiegano periodicamente), lo script lo segnala a "
                         "video ma calcola comunque fino a --zmax su richiesta esplicita")
    ap.add_argument("--klook", type=int, default=12,
                    help="numero di sub-aperture Doppler (look tomografici)")
    ap.add_argument("--ovs", type=float, default=4.0, help="oversampling in profondita'")
    ap.add_argument("--skip-vera-tomografia", action="store_true",
                    help="salta la tomografia Doppler vera (restano gli step 3-6)")
    ap.add_argument("--skip-self-test", action="store_true",
                    help="salta la validazione sintetica dell'inversione (sconsigliato)")
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args()
    pu.DBG = a.debug

    if a.nw is None or a.se is None or a.nw_lon is None or a.se_lon is None:
        if a.pyramid is None:
            sys.exit("Serve --pyramid {cheope,kefren} oppure tutti e 4 gli angoli "
                      "--nw/--nw-lon/--se/--se-lon espliciti.")
        preset = pu.PRESETS[a.pyramid]
        a.nw = a.nw or list(preset["nw"]); a.nw_lon = a.nw_lon or list(preset["nw_lon"])
        a.se = a.se or list(preset["se"]); a.se_lon = a.se_lon or list(preset["se_lon"])

    preset_dir = pu.PRESETS[a.pyramid]["outdir_name"] if a.pyramid else "goal_out"
    base_dir = os.path.join(HERE, preset_dir)
    outdir = a.outdir or os.path.join(base_dir, "acustica_vh")
    stack_base = (pu.PRESETS[a.pyramid].get("stack_outdir_name", preset_dir)
                  if a.pyramid else preset_dir)
    stackdir = a.stack or os.path.join(HERE, stack_base, "stack_slc")
    os.makedirs(outdir, exist_ok=True)
    boxpath = os.path.join(outdir, "box_vh.npz")
    # nessun fallback: questo script deve lavorare SOLO su dati VH reali
    # ritagliati ora, mai riusare box.npz misti VV+VH di altre run
    fallback_box = os.path.join(outdir, "_nessun_fallback_vh.npz")

    pyr_geom = None
    if not a.no_pyramid and a.pyramid:
        _p = pu.PRESETS[a.pyramid]
        pyr_geom = dict(base_m=_p["base_m"], h_m=_p["h_m"], label=_p["label"])

    lat_nw = pu.dms2dec(*a.nw[:3], a.nw[3]); lon_nw = pu.dms2dec(*a.nw_lon[:3], a.nw_lon[3])
    lat_se = pu.dms2dec(*a.se[:3], a.se[3]); lon_se = pu.dms2dec(*a.se_lon[:3], a.se_lon[3])
    latN, latS = max(lat_nw, lat_se), min(lat_nw, lat_se)
    lonW, lonE = min(lon_nw, lon_se), max(lon_nw, lon_se)
    aoi = (lonW, lonE, latS, latN)

    steps = pu.parse_steps(a.steps)

    log("=" * 70)
    log("PIRAMIDE ACUSTICA VH — canale VH reinterpretato come vibrazione")
    log("acustica riflessa dalla struttura interna (SOLO polarizzazione VH)")
    log(f"piramide: {a.pyramid or 'custom'} | strati: {a.layers} | outdir: {outdir}")
    log("NOTA METODOLOGICA: il canale VH resta backscatter radar (EM, "
        "polarizzazione incrociata); trattarlo come 'onda acustica' e' "
        "un'analogia matematica (stessi coefficienti di Fourier / stessa "
        "inversione h(z)=A^H*Y del metodo Biondi-Malanga, che gia' usa "
        "lambda_sonic per il micro-Doppler), NON una registrazione di suono "
        "reale. Vedi docstring in cima a questo file e CLAUDE.md.")
    log("=" * 70)

    if not a.skip_self_test:
        self_test_tomografia()

    if 3 in steps:
        pu.step3_estrai_box(stackdir, "vh", aoi, boxpath, fallback_box,
                            track_filter=not a.no_track_filter)
    elif not os.path.exists(boxpath):
        sys.exit(f"[acustica-vh] box VH non trovato ({boxpath}) e lo step 3 "
                  f"non e' incluso in --steps: aggiungi 3.")

    tomo_info = None
    coerenza_info = None
    pv_comune_info = None
    if steps & {4, 5, 6}:
        n_layers = a.layers
        if n_layers is None:
            d_box = np.load(boxpath, allow_pickle=True)
            if "vh_scene" not in d_box.files:
                sys.exit(f"[acustica-vh] box {boxpath} non contiene 'vh_scene'")
            n_scenes_vh = len(d_box["vh_scene"])
            # tutte le COPPIE possibili (combinazioni N*(N-1)/2), non solo le
            # N-1 sequenziali: piu' scene VH nello stack -> piu' strati nel
            # grafico, migliore precisione statistica (piu' campioni per la
            # sintesi Fourier dello step 5, gli errori casuali di coppie
            # diverse tendono a mediarsi invece di pesare su pochi strati)
            n_layers = max(n_scenes_vh * (n_scenes_vh - 1) // 2, 1)
            log(f"--layers non specificato: uso TUTTE le {n_scenes_vh} scene VH "
                f"disponibili in {stackdir} e TUTTE le {n_layers} coppie "
                f"interferometriche possibili tra loro (combinazioni N*(N-1)/2, "
                f"non solo le {max(n_scenes_vh - 1, 1)} sequenziali) per il "
                f"numero massimo di strati e la miglior precisione")
        per_pol, x, y, coh = pu.step4_array_12(boxpath, n_layers, outdir)
        if "vh" not in per_pol:
            sys.exit("[acustica-vh] il box non contiene la polarizzazione VH")
        d_arr = np.load(os.path.join(outdir, "array3d_12strati.npz"), allow_pickle=True)
        grafico_array_3d_vh(per_pol["vh"][0], x, y, d_arr["etich"], outdir)
        if steps & {5, 6}:
            k_strati = per_pol["vh"][0].shape[0]
            prof_info = pu._profondita_massima_fourier(boxpath, k_strati)
            nz_arg = a.nz
            if nz_arg is None:
                z_amb_pre = prof_info[1]
                nz_arg = max(int(round(z_amb_pre / a.dz_onda)) + 1, 2)
                log(f"--nz non specificato: onda campionata ogni --dz-onda="
                    f"{a.dz_onda:.2f} m sulla finestra z_amb={z_amb_pre:.0f} m "
                    f"-> NZ={nz_arg} campioni (dz effettivo "
                    f"{z_amb_pre / (nz_arg - 1):.3f} m)")
            (w_pol, freq_pol, env_pol, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT, riferimenti) = \
                pu.step5_grafico_onde(per_pol, x, y, aoi, outdir, pyr_geom,
                                      NZ=nz_arg, STEP_X=a.step_x, STEP_Y=a.step_y,
                                      pyramid=a.pyramid, stackdir=stackdir,
                                      prof_info=prof_info)
            if 6 in steps:
                # soglia PIENO/VUOTO RIDEFINITA: MEDIA della frequenza
                # istantanea su TUTTI i campioni della forma d'onda (ogni dz,
                # non solo il picco per pixel come nella versione precedente):
                # campione con frequenza sopra la media -> VUOTO, sotto la
                # media -> PIENO
                soglia_step6 = a.soglia
                if soglia_step6 is None:
                    dz_eff = float(z[1] - z[0]) if len(z) > 1 else float("nan")
                    tutti_freq = np.concatenate([freq_pol[p].ravel() for p in freq_pol])
                    soglia_step6 = float(np.mean(tutti_freq))
                    log(f"soglia PIENO/VUOTO ridefinita: MEDIA della frequenza "
                        f"istantanea su tutti i {tutti_freq.size} campioni della "
                        f"forma d'onda (ogni dz={dz_eff:.3f} m) = "
                        f"{soglia_step6:.4f} 1/m; frequenza > media -> VUOTO, "
                        f"< media -> PIENO; override con --soglia")
                pu.step6_variazioni(freq_pol, env_pol, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT,
                                    outdir, soglia=soglia_step6, riferimenti=riferimenti,
                                    prof_info=prof_info, isteresi=a.isteresi,
                                    max_variazioni=a.max_variazioni)
                # la selezione dello step 6 standard puo' lasciare punti orfani
                # (tetto senza fondo o viceversa): qui si ricostruisce una
                # selezione COERENTE, solo coppie complete per pixel
                d_var = np.load(os.path.join(outdir, "variazioni_frequenza.npz"),
                                allow_pickle=True)
                soglia_usata = float(d_var["soglia"])
                coerenza_info = variazioni_coerenti_vh(
                    freq_pol, z, x_m, y_m, z0, outdir, soglia_usata, a.isteresi,
                    a.max_variazioni, riferimenti, ZBOT, ZTOP, Lx, Ly)
                # metodo alternativo: PIENO/VUOTO per singola scena grezza
                # (ampiezza, non frequenza di Hilbert), punti in comune tra
                # le n_scene classificazioni indipendenti
                pv_comune_info = pieno_vuoto_comune_vh(
                    boxpath, x_m, y_m, z0, Lx, Ly, ZTOP, ZBOT, outdir, riferimenti,
                    freq_pol=freq_pol, z=z, soglia_freq=soglia_usata)

    if not a.skip_vera_tomografia:
        tomo_info = tomografia_vh(stackdir, (a.nw, a.nw_lon, a.se, a.se_lon), outdir,
                                  zmax=a.zmax, klook=a.klook, ovs=a.ovs)

    verifica_finale(outdir, tomo_info, coerenza_info, pv_comune_info)

    log("=" * 70)
    log(f"FATTO. Output in: {outdir}")
    log("=" * 70)


if __name__ == "__main__":
    main()

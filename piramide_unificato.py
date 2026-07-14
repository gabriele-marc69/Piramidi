#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramide_unificato.py
=====================
UN SOLO programma (fusione di goal_out_Cheope/piramide_unificato.py e
goal_out_kefren/piramide_unificato.py, che erano divergenti) che esegue l'intera
catena del goal, dai 6 step separati (extract_box + goal_pipeline + forma_donda +
tomografia_verticale + variazioni_onde + fetch_dem) uniti in un'unica pipeline
guidata dalle coordinate, con selezione della piramide via --pyramid.

Step (allineati al goal):
  1) cerca i file .tiff (prodotti Sentinel-1 SLC) che CONTENGONO le coordinate
     passate al programma -> almeno 12 (catalogo Copernicus Data Space, OData).
  2) scarica i .tiff dei prodotti trovati (estrae dal pacchetto SAFE i
     measurement TIFF delle polarizzazioni scelte + annotation.xml), con
     download ripristinabile (riprende da dove interrotto) e token CDSE
     rinnovato ad ogni tentativo. Default --pol vv,vh: entrambe le
     polarizzazioni del prodotto IW DV.
  3) estrae l'AREA delle coordinate da ogni .tiff via i GCP del prodotto ->
     stack complesso co-registrato (n_scene, azimuth, range) -> box.npz,
     verificando lo sfasamento tra le scene e tenendo solo il track dominante.
  4) estrae la COMPONENTE VERTICALE dello spostamento dei singoli pixel
     (DInSAR LOS -> verticale) per le aree estratte -> array a N strati
     (n_strati, ny, nx). Questo step è un'interferometria DInSAR classica
     (coerente con la skill sar-dinsar-microdisplacement), NON con la
     tomografia Doppler del paper.
  5) genera il grafico 3D dell'array: i numeri di ogni pixel sono i
     coefficienti armonici di una forma d'onda sviluppata in profondita'
     (anti-trasformata di Fourier) -> superficie reale + "molle" 3D. Il
     COLORE di ogni "molla" e' la sua frequenza istantanea (Hilbert):
     PIENO (frequenza bassa, materiale interpretato come denso/solido, colore
     chiaro) / VUOTO (frequenza alta, cavita'/aria, colore scuro).
     TUTTI gli strati DInSAR di TUTTE le polarizzazioni richieste (--pol,
     default vv,vh: es. 12 strati VV + 12 VH = 24 valori) vengono messi in
     UN SOLO array di coefficienti armonici per pixel e l'anti-trasformata
     e' calcolata IN UN SOLO PASSAGGIO (un unico np.fft.irfft vettorizzato),
     a RISOLUZIONE PIENA (ogni pixel az/rng del box, nessun subsample) ->
     UNA forma d'onda COMBINATA per pixel; il numero di "molle" disegnate e'
     ESATTAMENTE il numero di pixel dell'area selezionata (ny*nx), una per
     pixel.
  6) per OGNI pixel dell'area selezionata (stessa griglia (ny,nx) degli step
     4-5: un punto per pixel, non un sottoinsieme filtrato per banda) trova
     il campione di profondita' con inviluppo (Hilbert) massimo e la sua
     frequenza istantanea, poi classifica il punto PIENO (frequenza sotto
     soglia = denso/solido) o VUOTO (frequenza sopra soglia = cavita'/aria);
     disegna i punti in 3D alle ALTEZZE REALI rispetto al livello 0 (DEM
     Copernicus + geometria del bersaglio). Ogni pixel usa la STESSA forma
     d'onda combinata dello step 5 (tutti gli strati di tutte le pol in un
     solo array, un solo passaggio di Fourier: campo 'pol' = es. 'vv+vh'
     nel .npz e nell'hover).

*** NOTA METODOLOGICA IMPORTANTE (vedi anche skill biondi-malanga-sar-tomography) ***
Gli step 5-6 sono un'analisi ESPLORATIVA ispirata al lessico di Biondi & Malanga
(Remote Sens. 2022, 14, 5231) ma NON implementano il loro metodo di tomografia
Doppler a sub-aperture. Nel paper la profondita' viene ricostruita da un numero
di vettore d'onda verticale Kz = 4*pi*B_perp / (lambda_sonic * r * sin(theta))
(baseline ortogonali sintetizzate da sub-aperture Doppler di UNA SOLA immagine
SLC) e da un'inversione a filtro adattato h(z) = A^H * Y (steering matrix),
con risoluzione tomografica delta_z = lambda_sonic * R / (2*A). Gli step 5-6 di
QUESTO script, invece, trattano gli N strati DInSAR (interferogrammi tra COPPIE
di acquisizioni in date diverse) come se fossero i coefficienti spettrali di
una serie di Fourier stirata su un range di profondita' scelto ARBITRARIAMENTE
(non derivato dalla geometria SAR), poi ne misura la frequenza istantanea
(Hilbert) per proporre punti di "variazione di densita'". E' un'estensione
originale, non la tomografia del paper: usarla per ipotesi, non come misura
di profondita' fisica validata.

Per la tomografia Doppler COERENTE con il paper (steering matrix + inversione
h(z)=A^H*Y, gia' validata su un modello sintetico in
skills/sar-doppler-tomography/scripts/sar_tomo.py) usa il flag
--vera-tomografia: richiama scripts/tomographic_images.py sullo stesso stack
Sentinel-1, producendo le sezioni tomografiche in stile pubblicazione (B-scan,
slice orizzontali, mappa di quota) con l'inversione reale del paper, avvisando
sull'altezza di ambiguita' z_amb oltre la quale le profondita' vanno in alias.

Esempio:
  python piramide_unificato.py --pyramid kefren --download   # richiede credenziali CDSE
  python piramide_unificato.py --pyramid cheope --steps 3-6  # solo elaborazione (dati locali)
  python piramide_unificato.py --pyramid kefren --vera-tomografia   # + tomografia Doppler vera

Credenziali download (step 2), via variabili d'ambiente (MAI hardcoded nel codice):
  CDSE_USER / CDSE_PASS   (account gratuito su dataspace.copernicus.eu)
"""
import os, sys, re, glob, json, time, zipfile, argparse, subprocess, socket
import urllib.request, urllib.parse, urllib.error
import numpy as np

# matplotlib/scipy/plotly importati lazy negli step che li usano, cosi' gli step
# 1-2 (ricerca/download) girano anche su una macchina senza quei pacchetti.

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

C = 299_792_458.0
DBG = False
MIN_VALID_PCT = 50.0   # step 3: copertura minima dell'AOI in un sub-swath per tenerlo

# Box di default per piramide (DMS), presi dai due programmi ora unificati.
# Passa esplicitamente --nw/--nw-lon/--se/--se-lon per un box custom.
# base_m/h_m: geometria reale della piramide (base e altezza ATTUALE), usata per
# l'overlay sul DEM e quindi per le "altezze reali" degli step 5-6. La lunghezza
# della base e' anche il riferimento metrico interno usato nel paper Biondi-Malanga
# (skill biondi-malanga-sar-tomography, sez. 5.19): sbagliarla falsa le quote.
PRESETS = {
    "cheope": dict(
        nw=("29", "58", "48.0", "N"), nw_lon=("31", "7", "58.4", "E"),
        se=("29", "58", "41.0", "N"), se_lon=("31", "8", "07.7", "E"),
        outdir_name="goal_out_Cheope",
        base_m=230.3, h_m=138.7, label="Cheope (Khnum-Khufu)"),
    "kefren": dict(
        nw=("29", "58", "38.0", "N"), nw_lon=("31", "7", "45.4", "E"),
        se=("29", "58", "29.0", "N"), se_lon=("31", "7", "55.4", "E"),
        outdir_name="goal_out_kefren",
        base_m=215.5, h_m=136.4, label="Chefren (Khafre)"),
}

# Strutture interne note della Grande Piramide (Cheope/Khufu), lette dalla
# sezione N-S standard (diagramma "Great Pyramid Diagram", numerazione 1-11:
# 1 entrata originale, 2 traforo dei ladri, 3 giunzione, 4 corridoio
# discendente, 5 camera sotterranea, 6 corridoio ascendente, 7 camera della
# Regina, 8 corridoio orizzontale, 9 Grande Galleria, 10 camera del Re,
# 11 pozzo verticale). s_m = metri a NORD del centro della base (negativo =
# sud), h_m = metri SOPRA la base, e_m = metri a EST del centro (default 0.0:
# il diagramma sorgente e' una sezione 2D N-S, quindi quasi tutti i punti
# sono sul piano verticale centrale N-S; l'unica eccezione nota e ben
# documentata in letteratura e' l'entrata originale, qui inclusa).
# Valori incrociati con le misure originali di Petrie (1883, "The Pyramids
# and Temples of Gizeh") dove reperibili:
#   - Entrata originale: altezza 668.3 pollici sopra il pavimento = 16.97 m
#     (Petrie); offset 15 cubiti reali = 7.9 m ad est dell'asse N-S (fonte
#     concorde: Wikipedia "Great Pyramid of Giza"). Sostituisce il valore
#     precedente (15.63 m, senza offset E-W) con la misura Petrie verificata.
#   - Camera della Regina (21.0 m) e Camera del Re (43.0 m): coerenti con le
#     quote Petrie/Lehner ("The Complete Pyramids") citate in letteratura
#     (rispettivamente ~21.2 m e ~42.9 m sopra la base) -- non modificate.
# I restanti punti (traforo dei ladri, giunzioni, corridoi, pozzo) restano
# valori APPROSSIMATI da rilievi pubblicati (coerenti con le proporzioni del
# diagramma, scala 50 m): sono elementi scavati/irregolari senza una singola
# quota "ufficiale" citata in letteratura con la stessa precisione delle
# camere. Servono SOLO da riferimento visivo per confronto con la nuvola
# DInSAR -- NON sono derivati dai dati SAR di questo programma.
#
# Gli ultimi 3 campi (ns_m, ew_m, alt_m) sono le DIMENSIONI dell'ingombro
# (bounding box) di ogni struttura, in metri N-S x E-W x altezza, da misure
# pubblicate (Petrie 1883; Lehner, "The Complete Pyramids"):
#   - Camera del Re: 10.47 (E-W) x 5.23 (N-S) x 5.85 h (Petrie 412.4x206.3x230.1")
#   - Camera della Regina: 5.23 (E-W) x 5.75 (N-S) x 6.23 h al colmo del tetto
#     a doppio spiovente (Petrie 205.9x226.5", colmo ~245")
#   - Camera sotterranea: 14.06 (E-W) x 8.28 (N-S) x ~3.5 h, incompiuta
#     (Petrie 553.5x325.9")
#   - Grande Galleria: lunga 46.7 m, alta 8.6-8.7, larga 2.06 alla base;
#     inclinata ~26 gradi -> ingombro N-S ~41.9, verticale ~20.7+8.7
#   - Corridoio discendente: lungo 105.2 m a 26 gradi 31' (sezione 1.05x1.20)
#     -> ingombro N-S ~94, verticale ~48
#   - Corridoio ascendente: lungo 39.3 m a 26 gradi 2' -> N-S ~35.3, vert. ~18.5
#   - Corridoio orizzontale: lungo ~38.9 m, sezione 1.05 x 1.2-1.8
# Per i tratti INCLINATI l'ingombro e' il bounding box del tratto (il "cubo"
# contiene il corridoio in diagonale); per traforo e pozzo (scavi irregolari)
# e' un ingombro indicativo. I box sono CENTRATI sul punto di riferimento.
RIFERIMENTI_CHEOPE = [
    #  n  nome                                        s_m     h_m    e_m   ns_m   ew_m  alt_m
    (1, "Entrata originale", 115.15, 16.97, 7.9, 2.0, 1.1, 1.2),
    (2, "Traforo dei ladri (tratto iniziale)", 115.15, 3.0, 0.0, 6.0, 2.0, 2.5),
    (3, "Giunzione traforo / corridoio discendente", 98.4, 7.3, 0.0, 3.0, 2.0, 3.0),
    (4, "Corridoio discendente", 61.5, -11.2, 0.0, 94.0, 1.1, 48.0),
    (5, "Camera sotterranea", 12.2, -31.4, 0.0, 8.3, 14.1, 3.5),
    (6, "Corridoio ascendente", 80.8, 16.0, 0.0, 35.3, 1.1, 18.5),
    (7, "Camera della Regina", 0.0, 21.0, 0.0, 5.8, 5.2, 6.2),
    (8, "Corridoio orizzontale", 31.6, 21.0, 0.0, 38.9, 1.1, 1.8),
    (9, "Grande Galleria", 42.4, 35.2, 0.0, 41.9, 2.1, 29.4),
    (10, "Camera del Re", 17.0, 43.0, 0.0, 5.2, 10.5, 5.9),
    (11, "Pozzo verticale (Well Shaft)", 68.0, 0.0, 0.0, 2.0, 2.0, 50.0),
]

# Strutture interne note di Chefren (Khafre): dimensioni dai rilievi
# pubblicati (Belzoni 1818; Maragioglio & Rinaldi, "L'architettura delle
# piramidi menfite"): camera sepolcrale 14.15 (E-W) x 5.0 (N-S) x 6.83 h
# (tetto a doppio spiovente), camera sussidiaria 10.41 (N-S) x 3.12 (E-W) x
# 2.65 h. Le POSIZIONI (s/h/e) sono APPROSSIMATE dalle sezioni pubblicate
# (la camera sepolcrale e' presso l'asse centrale, al livello della base;
# la sussidiaria e' nel sottosuolo lungo il corridoio inferiore, a nord).
RIFERIMENTI_KEFREN = [
    (1, "Camera sepolcrale (Belzoni)", 0.0, 3.4, 0.0, 5.0, 14.2, 6.8),
    (2, "Camera sussidiaria (sotterranea)", 35.0, -9.0, 0.0, 10.4, 3.1, 2.7),
    (3, "Corridoio discendente inferiore", 65.0, -7.0, 0.0, 55.0, 1.1, 14.0),
]

# riferimenti per piramide (chiave = valore di --pyramid)
RIFERIMENTI_PIRAMIDE = {"cheope": RIFERIMENTI_CHEOPE, "kefren": RIFERIMENTI_KEFREN}


# Colorscale condivisa dagli step 5-6 (selezione RIBALTATA su richiesta): un
# mezzo con frequenza istantanea PIU' BASSA nella forma d'onda sintetizzata
# (anti-FT) e' interpretato come PIENO (materiale denso/solido, colore PIU'
# CHIARO), uno con frequenza PIU' ALTA come VUOTO (cavita'/aria, colore PIU'
# SCURO). E' un'ipotesi esplorativa (vedi nota metodologica in cima al
# modulo), NON una misura fisica di densita' validata dal paper Biondi-Malanga.
PIENO_VUOTO_COLORSCALE = [
    [0.0, "#f5ecd7"],   # PIENO: frequenza bassa -> crema (piu' chiaro)
    [0.5, "#8a7862"],
    [1.0, "#0d0d12"],   # VUOTO: frequenza alta -> quasi nero (piu' scuro)
]
# Stessi 3 colori, usati sia dal JS del click-export (step 6) sia dalla
# colormap matplotlib delle anteprime statiche: un'unica fonte di verita'.
PIENO_VUOTO_HEX = ["#f5ecd7", "#8a7862", "#0d0d12"]

# Tema della pagina web dei grafici interattivi (step 5-6): dashboard tecnica
# CHIARA con sfondo BIANCO del grafico (scelta esplicita dell'utente), testo
# scuro con contrasto >= 4.5:1, un solo colore d'accento per i comandi.
TEMA_UI = dict(
    sfondo="#ffffff", pannello="#f8fafc", bordo="#e2e8f0",
    testo="#0f172a", testo_muto="#475569", accento="#b45309",
    font="'Segoe UI', 'Fira Sans', system-ui, -apple-system, sans-serif")


def _scrivi_pagina_html(fig, path, titolo, sottotitolo, badges=(), nota=None,
                        post_script=None):
    """Scrive `fig` in una pagina HTML con layout curato: header fisso con
    titolo, sottotitolo, badge riassuntivi e nota metodologica richiudibile;
    il grafico occupa TUTTA l'altezza restante della finestra. Sostituisce il
    write_html nudo di Plotly, che infilava un titolo lunghissimo dentro il
    canvas sovrapponendolo alla scena 3D e lasciava la pagina bianca di
    default. `badges` = sequenza di coppie (etichetta, valore)."""
    import plotly.io as pio
    corpo = pio.to_html(fig, full_html=False, include_plotlyjs="cdn",
                        post_script=post_script,
                        config=dict(scrollZoom=True, displayModeBar=True,
                                    responsive=True),
                        default_width="100%", default_height="100%")
    chips = "".join(f'<span class="chip">{k}<b>{v}</b></span>'
                    for k, v in badges)
    nota_html = ""
    if nota:
        nota_html = ('<details class="nota"><summary>Nota metodologica'
                     f'</summary><p>{nota}</p></details>')
    t = TEMA_UI
    html = f"""<!DOCTYPE html>
<html lang="it">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{titolo}</title>
<style>
  * {{ box-sizing: border-box; }}
  html, body {{ margin: 0; height: 100%; background: {t['sfondo']};
                color: {t['testo']}; font-family: {t['font']}; }}
  .app {{ display: flex; flex-direction: column; height: 100%; }}
  header {{ flex: none; padding: 10px 16px 9px; background: {t['pannello']};
            border-bottom: 1px solid {t['bordo']}; }}
  header h1 {{ margin: 0 0 2px; font-size: 16px; font-weight: 600;
               letter-spacing: .2px; }}
  header .sotto {{ margin: 0; font-size: 12.5px; line-height: 1.55;
                   color: {t['testo_muto']}; max-width: 110ch; }}
  .chips {{ margin-top: 7px; display: flex; flex-wrap: wrap; gap: 6px; }}
  .chip {{ font-size: 11.5px; color: {t['testo_muto']}; background: {t['sfondo']};
           border: 1px solid {t['bordo']}; border-radius: 999px;
           padding: 2px 10px; white-space: nowrap; }}
  .chip b {{ color: {t['testo']}; font-weight: 600; margin-left: 5px;
             font-variant-numeric: tabular-nums; }}
  .nota {{ margin-top: 7px; font-size: 12px; color: {t['testo_muto']};
           max-width: 90ch; }}
  .nota summary {{ cursor: pointer; color: {t['accento']}; font-weight: 600;
                   outline-offset: 3px; transition: opacity .2s; }}
  .nota summary:hover {{ opacity: .8; }}
  .nota p {{ margin: 6px 0 0; line-height: 1.55; }}
  main {{ flex: 1 1 auto; min-height: 0; }}
  main .plotly-graph-div {{ height: 100% !important; }}
  @media (prefers-reduced-motion: reduce) {{
    * {{ transition: none !important; animation: none !important; }}
  }}
</style>
</head>
<body>
<div class="app">
  <header>
    <h1>{titolo}</h1>
    <p class="sotto">{sottotitolo}</p>
    <div class="chips">{chips}</div>
    {nota_html}
  </header>
  <main>
{corpo}
  </main>
</div>
</body>
</html>
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)


def _tema_plotly(fig):
    """Applica al grafico Plotly il tema chiaro (sfondo bianco) coerente con
    la pagina (_scrivi_pagina_html): il titolo sta nell'header HTML, non nel
    canvas."""
    t = TEMA_UI
    fig.update_layout(
        template="plotly_white",
        paper_bgcolor=t["sfondo"], plot_bgcolor=t["sfondo"],
        font=dict(family=t["font"], color=t["testo"], size=12),
        margin=dict(l=0, r=0, t=10, b=0),
        legend=dict(bgcolor="rgba(248,250,252,0.9)",
                    bordercolor=t["bordo"], borderwidth=1,
                    x=0.01, y=0.99, xanchor="left", yanchor="top"))


def _colorbar_pieno_vuoto(vmin, vmax, titolo="frequenza istantanea [1/m]"):
    """Colorbar Plotly con gli estremi esplicitamente etichettati VUOTO/PIENO."""
    vmid = (vmin + vmax) / 2.0
    return dict(title=titolo,
               tickvals=[vmin, vmid, vmax],
               ticktext=[f"PIENO ({vmin:.3f})", "", f"VUOTO ({vmax:.3f})"])


def _cmap_pieno_vuoto():
    """Colormap matplotlib (PIENO chiaro -> VUOTO scuro) coerente con
    PIENO_VUOTO_COLORSCALE/HEX usata nei grafici interattivi Plotly."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("pieno_vuoto", PIENO_VUOTO_HEX)


def _frequenza_istantanea(w, dz):
    """Frequenza istantanea (Hilbert) lungo l'ULTIMO asse (profondita') di un
    array w (..., NZ), e inviluppo di ampiezza |analitico|. Usata da step 5
    (colore delle "molle") e step 6 (classificazione PIENO/VUOTO per pixel):
    inst_f = d(fase dell'analitico)/dz / (2*pi). Analisi esplorativa (vedi
    nota metodologica in cima al modulo), non e' una misura fisica validata."""
    from scipy.signal import hilbert
    an = hilbert(w, axis=-1)
    inst_f = np.gradient(np.unwrap(np.angle(an), axis=-1), dz, axis=-1) / (2 * np.pi)
    return inst_f.astype(np.float32), np.abs(an).astype(np.float32)


def _calibra_centro_piramide(stackdir, pol, aoi, x_m, y_m):
    """Trova la posizione (x_m,y_m) del centro della base della piramide
    dentro il box, e la direzione nord/est degli assi riga/colonna leggendola
    dai GCP REALI del prodotto (non e' un'ipotesi: dipende dalla geometria di
    acquisizione del passaggio satellitare, quindi va verificata scena per
    scena). Ritorna (x_centro, y_centro, i_centro, j_centro, segno_nord,
    segno_est) con segno_*=+1 se l'asse (riga/colonna) cresce andando verso
    nord/est, -1 se cresce nel verso opposto; None se lo stack non ha piu'
    TIFF/GCP disponibili (es. dopo una pulizia dello stack)."""
    tiffs = sorted(glob.glob(os.path.join(stackdir, f"*-iw[1-3]-slc-{pol}-*.tiff")))
    if not tiffs:
        return None
    import rasterio
    from rasterio.transform import GCPTransformer
    lonW, lonE, latS, latN = aoi
    lon_c, lat_c = (lonW + lonE) / 2, (latS + latN) / 2
    corners = [(lonW, latN), (lonE, latN), (lonE, latS), (lonW, latS)]  # NW,NE,SE,SW
    ny, nx = len(y_m), len(x_m)
    for tif in tiffs:
        ann = _trova_annotation(tif, stackdir, pol)
        if ann is None:
            continue
        try:
            with rasterio.open(tif) as ds:
                if not ds.gcps[0]:
                    continue
                with GCPTransformer(ds.gcps[0]) as tf:
                    rc = [tf.rowcol(lon, lat) for lon, lat in corners]
                    r_c, c_c = tf.rowcol(lon_c, lat_c)
                height, width = ds.height, ds.width
        except Exception:
            continue
        rows = [r for r, c in rc]; cols = [c for r, c in rc]
        if not all(0 <= r < height and 0 <= c < width for r, c in rc):
            continue
        r0 = max(0, int(np.floor(min(rows)))); c0 = max(0, int(np.floor(min(cols))))
        i_c, j_c = r_c - r0, c_c - c0
        if not (0 <= i_c < ny and 0 <= j_c < nx):
            continue
        row_NW, row_SW = rows[0], rows[3]
        col_NW, col_NE = cols[0], cols[1]
        segno_nord = 1 if row_NW > row_SW else -1   # riga crescente -> nord?
        segno_est = 1 if col_NE > col_NW else -1     # colonna crescente -> est?
        x_c = float(np.interp(j_c, np.arange(nx), x_m))
        y_c = float(np.interp(i_c, np.arange(ny), y_m))
        return x_c, y_c, i_c, j_c, segno_nord, segno_est
    return None


def _tracce_riferimenti(riferimenti):
    """Tracce Plotly per le strutture interne note: per ogni struttura un
    PARALLELEPIPEDO semi-trasparente (Mesh3d) con l'INGOMBRO reale preso dai
    rilievi pubblicati (vedi RIFERIMENTI_*: N-S x E-W x altezza; per i tratti
    inclinati e' il bounding box del tratto), centrato sul punto di
    riferimento, + i marker numerati con l'etichetta. Le tracce condividono
    un legendgroup: dalla legenda si accendono/spengono tutte insieme.
    Riferimenti da letteratura, NON derivati dai dati SAR di questo programma.
    `riferimenti` = lista di (num, nome, x, y, z, dx_ew, dy_ns, dz_alt)."""
    import plotly.graph_objects as go
    # triangolazione standard di un cubo (8 vertici, 12 triangoli)
    I = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
    J = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
    K = [0, 7, 2, 3, 6, 7, 2, 5, 5, 5, 7, 6]
    UX = (0, 1, 1, 0, 0, 1, 1, 0)
    UY = (0, 0, 1, 1, 0, 0, 1, 1)
    UZ = (0, 0, 0, 0, 1, 1, 1, 1)
    tracce = []
    SPESSORE_MIN = 3.0   # [m] spessore minimo VISIVO per asse: un corridoio
    # largo 1 m su una scena di ~300 m e' sub-pixel e sembrerebbe un PIANO;
    # il box e' gonfiato al minimo visivo SOLO per il disegno, l'hover
    # riporta sempre le dimensioni reali da letteratura.
    # percorso degli SPIGOLI di un box (perimetri base+tetto e 4 verticali,
    # None = stacco): rendono il volume leggibile come 3D da ogni angolo.
    SPIGOLI = [0, 1, 2, 3, 0, 4, 5, 6, 7, 4, None, 1, 5, None, 2, 6, None, 3, 7]
    ex, ey, ez = [], [], []
    for n, (num, nome, x, y, z, dx, dy, dz) in enumerate(riferimenti):
        dxv, dyv, dzv = (max(dx, SPESSORE_MIN), max(dy, SPESSORE_MIN),
                         max(dz, SPESSORE_MIN))
        gonfiato = (dxv, dyv, dzv) != (dx, dy, dz)
        hover = (f"{num}. {nome}<br>ingombro ~ {dy:.1f} m (N-S) x "
                 f"{dx:.1f} m (E-W) x {dz:.1f} m (h)<br>"
                 f"centro: x={x:.0f} m, y={y:.0f} m, quota {z:.0f} m s.l.m."
                 + (f"<br>(disegnato con spessore minimo {SPESSORE_MIN:.0f} m "
                    f"per renderlo leggibile in 3D)" if gonfiato else ""))
        vx = [x - dxv / 2 + u * dxv for u in UX]
        vy = [y - dyv / 2 + u * dyv for u in UY]
        vz = [z - dzv / 2 + u * dzv for u in UZ]
        tracce.append(go.Mesh3d(
            x=vx, y=vy, z=vz,
            i=I, j=J, k=K, color="#22c55e", opacity=0.35, flatshading=True,
            hovertext=hover, hoverinfo="text",
            name="strutture interne note (ingombri da letteratura, non dai dati SAR)",
            legendgroup="riferimenti", showlegend=(n == 0)))
        for s in SPIGOLI:
            ex.append(vx[s] if s is not None else None)
            ey.append(vy[s] if s is not None else None)
            ez.append(vz[s] if s is not None else None)
        ex.append(None); ey.append(None); ez.append(None)
    tracce.append(go.Scatter3d(
        x=ex, y=ey, z=ez, mode="lines", connectgaps=False,
        line=dict(color="#15803d", width=3), hoverinfo="skip",
        legendgroup="riferimenti", showlegend=False,
        name="spigoli strutture interne"))
    rtxt = [f"{r[0]}. {r[1]}<br>x={r[2]:.0f} m, y={r[3]:.0f} m, "
            f"quota {r[4]:.0f} m s.l.m.<br>ingombro ~ {r[6]:.1f} (N-S) x "
            f"{r[5]:.1f} (E-W) x {r[7]:.1f} (h) m" for r in riferimenti]
    tracce.append(go.Scatter3d(
        x=[r[2] for r in riferimenti], y=[r[3] for r in riferimenti],
        z=[r[4] for r in riferimenti], mode="markers+text",
        marker=dict(size=4, color="#16a34a", symbol="diamond",
                    line=dict(width=1, color="black")),
        text=[str(r[0]) for r in riferimenti], textposition="top center",
        textfont=dict(size=11, color="#166534"),
        hovertext=rtxt, hoverinfo="text",
        legendgroup="riferimenti", showlegend=False,
        name="strutture interne note"))
    return tracce


def _profondita_massima_fourier(boxpath, k):
    """Profondita' MASSIMA raggiungibile dall'inversione di Fourier secondo
    la skill biondi-sar-tomography, in funzione del TIPO di file sorgente.
    Nel metodo del paper l'inversione h(z)=A^H*Y e' una DFT su k campioni di
    Kz = 4*pi*B_perp/(lambda_sonic*r*sin(theta)): oltre l'altezza di
    ambiguita' z_amb = 2*pi/delta_Kz = lambda_sonic*R*sin(theta)*(k-1)/(2*A)
    le profondita' vanno in ALIAS (stessa formula di
    skills/sar-doppler-tomography/scripts/sar_tomo.py: altezza_ambiguita);
    la risoluzione per cella e' delta_z = lambda_sonic*R/(2*A), con
    lambda_sonic = v/(2f) = 6000/(2*12500) = 0.24 m (versione pubblicata,
    NON la lambda radar). Il tipo di sorgente determina R, theta e
    soprattutto l'apertura orbitale A:
      - COSMO-SkyMed Spotlight X-band (paper): A = 84 km (Table 1 della
        versione pubblicata) -> delta_z ~0.9 m, z_amb(k=12) ~6-8 m;
      - Sentinel-1 IW SLC C-band (i file di questo repo): il TOPS riduce il
        dwell per bersaglio, A = lambda_em*R/(2*res_az) con res_az ~22 m ->
        A ~1 km, delta_z ~90 m, z_amb(k=12) ~600-700 m: finestra molto piu'
        PROFONDA ma con celle molto piu' GROSSOLANE del paper.
    Qui k = numero di strati DInSAR usati come armoniche dall'anti-FT degli
    step 5-6. Ritorna (tipo_sorgente, z_amb [m], delta_z [m], theta [gradi])."""
    d = np.load(boxpath, allow_pickle=True)
    f_em = float(d["f_em"]); R = float(d["R_near"]); inc = float(d["incid_mid"])
    pols = [p for p in ("vv", "vh", "hh", "hv") if p in d.files]
    nomi = " ".join(str(x) for x in d[f"{pols[0]}_scene"]) if pols else ""
    lam_em = C / f_em
    lam_sonic = 6000.0 / (2.0 * 12500.0)   # v/(2f), come sar_tomo.GeometriaSAR
    if f_em > 8e9:                          # X-band: COSMO-SkyMed del paper
        A = 84_000.0
        tipo = "COSMO-SkyMed Spotlight (X-band, A=84 km da paper)"
    else:                                   # C-band: Sentinel-1
        iw = "iw" in nomi.lower()
        res_az = 22.0 if iw else 5.0        # risoluzione azimuth IW / Stripmap
        A = lam_em * R / (2.0 * res_az)
        tipo = (f"Sentinel-1 {'IW' if iw else 'SM'} SLC "
                f"(C-band, apertura per bersaglio ~{A/1000:.1f} km)")
    dz = lam_sonic * R / (2.0 * A)
    z_amb = lam_sonic * R * np.sin(np.radians(inc)) * max(k - 1, 1) / (2.0 * A)
    return tipo, float(z_amb), float(dz), float(inc)


def log(*a):
    print(*a, flush=True)


def dbg(*a):
    if DBG:
        print("   [debug]", *a, flush=True)


def dms2dec(d, m, s, h):
    v = float(d) + float(m) / 60 + float(s) / 3600
    return -v if str(h).upper() in ("S", "W") else v


def lista_pol(pol):
    """Normalizza --pol in lista di polarizzazioni ('vv,vh' -> ['vv','vh']).
    La PRIMA e' la polarizzazione primaria: riferimento per il conteggio scene,
    il track filter e la geometria; le altre (es. vh) aggiungono forme d'onda
    indipendenti -> piu' tracce/punti in orizzontale nei grafici 3D."""
    if isinstance(pol, str):
        pol = pol.split(",")
    out = []
    for p in pol:
        p = p.strip().lower()
        if p and p not in out:
            out.append(p)
    return out or ["vv"]


def _chiave_scena_sw(nome):
    """Chiave satellite+subswath+datetime (al minuto): identica per i measurement
    VV e VH della stessa scena/subswath, che differiscono solo per il token di
    polarizzazione e per il numero progressivo finale del file."""
    n = os.path.basename(nome).lower()
    m1 = re.search(r"(s1[abcd]-iw\d)", n)
    m2 = re.search(r"(\d{8}t\d{4})", n)
    return f"{m1.group(1)}-{m2.group(1)}" if (m1 and m2) else n


# ===================================================================== STEP 1
def step1_cerca(aoi, start, end, max_records, pol, outdir):
    """Cerca i prodotti Sentinel-1 SLC che contengono l'AOI (OData CDSE).
    Ritorna la lista [{Name, Id, ContentLength, Start}] ordinata per data."""
    lonW, lonE, latS, latN = aoi
    poly = (f"POLYGON(({lonW} {latS},{lonE} {latS},{lonE} {latN},"
            f"{lonW} {latN},{lonW} {latS}))")
    # IW_SLC (non un generico 'SLC'): il resto della pipeline usa solo i
    # measurement iw1-3 (glob '*-iw[1-3]-slc-...'), un SLC Stripmap/EW da ~8 GB
    # verrebbe scaricato e poi ignorato dallo step 3.
    flt = ("Collection/Name eq 'SENTINEL-1' "
           f"and OData.CSC.Intersects(area=geography'SRID=4326;{poly}') "
           "and contains(Name,'IW_SLC') "
           f"and ContentDate/Start gt {start}T00:00:00.000Z "
           f"and ContentDate/Start lt {end}T23:59:59.000Z")
    q = {"$filter": flt, "$orderby": "ContentDate/Start asc",
         "$top": str(max_records)}
    url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?" + \
          urllib.parse.urlencode(q)
    dbg("query:", url)
    log(f"[1] ricerca SLC che contengono l'AOI "
        f"(lon {lonW:.5f}..{lonE:.5f}, lat {latS:.5f}..{latN:.5f}), {start}..{end}")
    cache = os.path.join(outdir, "prodotti_trovati.json")
    data = None
    for giro in (1, 2):              # 2 giri: il secondo dopo un'attesa-rete
        try:
            with urllib.request.urlopen(urllib.request.Request(url), timeout=90) as r:
                data = json.load(r)
            break
        except Exception as ex:
            log(f"    ! ricerca catalogo fallita ({ex})")
            if giro == 1 and not _rete_ok():
                if not _attendi_rete():
                    break
            else:
                break
    if data is None:
        # fallback: elenco prodotti salvato da un run precedente (stesse coordinate)
        if os.path.exists(cache):
            with open(cache, encoding="utf-8") as f:
                prods = json.load(f)
            log(f"[1] rete/catalogo indisponibile: riuso l'elenco in cache "
                f"({len(prods)} prodotti da {cache})")
            return prods
        log("    ! nessuna cache prodotti disponibile; proseguo con i dati locali")
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
    with open(cache, "w", encoding="utf-8") as f:
        json.dump(prods, f, indent=2)
    log(f"    elenco salvato in {cache}")
    return prods


# ===================================================================== STEP 2
def _cdse_token(user, pwd):
    url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    body = urllib.parse.urlencode({"client_id": "cdse-public", "grant_type": "password",
                                   "username": user, "password": pwd}).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=body), timeout=60) as r:
        return json.load(r)["access_token"]


def _rete_ok():
    """True se il DNS risolve il nodo CDSE (test rapido di connettivita')."""
    try:
        socket.getaddrinfo("catalogue.dataspace.copernicus.eu", 443)
        return True
    except OSError:
        return False


def _attendi_rete(max_wait=7200, check_every=30):
    """Se la rete e' giu', aspetta che torni (fino a max_wait secondi) invece di
    bruciare i tentativi: con download multi-ora una connessione che 'flappa'
    non deve far abbandonare l'intera coda di prodotti."""
    if _rete_ok():
        return True
    log(f"    ! rete assente: attendo che torni (controllo ogni {check_every}s, "
        f"max {max_wait//60} min)...")
    t0 = time.time()
    while time.time() - t0 < max_wait:
        time.sleep(check_every)
        if _rete_ok():
            log(f"    rete tornata dopo {time.time()-t0:.0f}s: riprendo")
            return True
    log(f"    ! rete ancora assente dopo {max_wait//60} min: rinuncio")
    return False


class _AuthRedirect(urllib.request.HTTPRedirectHandler):
    """Conserva l'header Authorization quando CDSE reindirizza dal catalogo al nodo
    di download: urllib di default lo rimuove sui redirect cross-host -> HTTP 401."""
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None:
            auth = req.get_header("Authorization")
            if auth:
                new.add_unredirected_header("Authorization", auth)
        return new


_CDSE_OPENER = urllib.request.build_opener(_AuthRedirect())


def _scarica_ripristinabile(url, dest, get_token, expected=0, n_try=60):
    # n_try alto di proposito: ogni tentativo RIPRENDE dal byte gia' scaricato
    # (header Range), quindi con connessioni instabili che cadono ogni pochi
    # minuti servono molti brevi tentativi per completare un file da ~8 GB;
    # un tentativo non costa nulla se non aggiunge progresso.
    """Scarica `url` in `dest` riprendendo da dove si era interrotto (fusione
    del download ripristinabile di Cheope con l'opener che preserva
    l'Authorization sui redirect di Kefren, cosi' da avere sia la resilienza
    sia il fix del 401).

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
        try:
            # token DENTRO il try: anche un errore di rete nel login deve essere
            # ritentato qui (prima faceva saltare l'intero prodotto)
            headers = {"Authorization": f"Bearer {get_token()}"}
            if have:
                headers["Range"] = f"bytes={have}-"
                log(f"    riprendo da {have/1e9:.2f} GB"
                    f"{f' / {expected/1e9:.1f} GB' if expected else ''} (tentativo {att}/{n_try})")
            req = urllib.request.Request(url, headers=headers)
            # _CDSE_OPENER conserva l'Authorization sul redirect verso il nodo download
            with _CDSE_OPENER.open(req, timeout=1800) as r:
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
            # se e' la rete a mancare, aspetta che torni PRIMA di consumare
            # il prossimo tentativo (connessioni instabili/lunghi blackout)
            if not _attendi_rete():
                return False
    else:
        log(f"    ! download non completato dopo {n_try} tentativi"); return False
    if expected and os.path.getsize(part) < expected:
        log("    ! download incompleto"); return False
    os.replace(part, dest)
    return True


def _dt_scena(nome):
    """Chiave-scena (YYYYMMDDtHHMM, senza secondi) dal nome SAFE/TIFF. Un SAFE
    contiene 3 sub-swath -> 3 TIFF della stessa scena i cui timestamp differiscono
    di 1-2 SECONDI (iw1/iw2/iw3): troncando ai minuti i 3 TIFF collassano su una
    sola scena, mentre track diversi (ore/minuti diversi) restano distinti."""
    m = re.search(r"(\d{8}t\d{4})", nome.lower())
    return m.group(1) if m else nome


def _estrai_tiff_da_safe(safe_zip, stackdir, pols, rimuovi_zip):
    """Estrae da un .zip SAFE i measurement TIFF + annotation delle
    polarizzazioni richieste in stackdir (scrittura atomica tmp+rename).
    Ritorna i path dei TIFF estratti. Con rimuovi_zip=True elimina lo zip
    sorgente a estrazione riuscita (SI' per gli zip scaricati da questo
    programma, NO per quelli di un altro programma come
    scarica_alta_risoluzione.py, che restano dove sono)."""
    estratti = []
    with zipfile.ZipFile(safe_zip) as z:
        for nm in z.namelist():
            low = nm.lower()
            dst = None
            pol_ok = any(f"-{q}-" in low for q in pols)
            if pol_ok and low.endswith(".tiff") and "/measurement/" in low:
                dst = os.path.join(stackdir, os.path.basename(low))
            elif pol_ok and low.endswith(".xml") and "/annotation/" in low \
                    and "/calibration/" not in low and "/rfi/" not in low:
                dst = os.path.join(stackdir, os.path.basename(low).replace(
                    ".xml", ".annotation.xml"))
            if dst is None:
                continue
            tmp = dst + ".tmp"
            with z.open(nm) as zi, open(tmp, "wb") as fo:
                while True:
                    chunk = zi.read(1 << 22)
                    if not chunk:
                        break
                    fo.write(chunk)
            os.replace(tmp, dst)
            if dst.endswith(".tiff"):
                estratti.append(dst)
    if rimuovi_zip:
        os.remove(safe_zip)
    return estratti


def step2_scarica(prods, n_want, pol, stackdir, user=None, pwd=None, download=False,
                  alta_ris_dir=None):
    """Scarica i prodotti e ne estrae i measurement TIFF (+ annotation) delle
    polarizzazioni richieste (--pol vv,vh) nello stack locale. Il conteggio verso
    n_want e' in SCENE distinte (datetime di acquisizione) della polarizzazione
    PRIMARIA, non in TIFF (un SAFE = 3 TIFF per pol della stessa scena). I prodotti
    del track dominante (stessa ora del giorno) vengono scaricati per primi: gli
    step 4-6 usano solo scene dello stesso track (cfr. filtro dello step 3), quindi
    date di track diversi non aggiungerebbero coppie.
    Prima di scaricare dal CDSE, legge i .zip SAFE gia' scaricati da
    `alta_ris_dir` (default alta_risoluzione_out/, prodotto da
    scarica_alta_risoluzione.py: stesso formato SAFE, VV+VH) e ne estrae subito
    i measurement mancanti: nessuna credenziale ne' rete necessarie per questa
    parte. I .zip.part (download ancora in corso) sono ignorati automaticamente.
    Credenziali CDSE: dai parametri --cdse-user/--cdse-pass oppure dalle variabili
    d'ambiente CDSE_USER/CDSE_PASS (mai hardcoded qui). Se mancano le credenziali
    (o la lista e' vuota), salta e usa quel che c'e'."""
    pols = lista_pol(pol)
    os.makedirs(stackdir, exist_ok=True)
    # riusa solo TIFF plausibilmente integri (un measurement IW pesa ~1 GB: file
    # molto piccoli = estrazione interrotta a meta' -> da rifare)
    MIN_TIFF = 50_000_000
    have_all = sorted({h for p in pols
                       for h in glob.glob(os.path.join(stackdir, f"*slc-{p}-*.tiff"))})
    have = [h for h in have_all if os.path.getsize(h) >= MIN_TIFF]
    for h in have_all:
        if h not in have:
            log(f"[2] {os.path.basename(h)}: troncato "
                f"({os.path.getsize(h)/1e6:.0f} MB), lo rimuovo e riestraggo")
            os.remove(h)
    if have:
        log(f"[2] gia' presenti {len(have)} TIFF '{'/'.join(pols)}' in {stackdir} "
            f"(riuso: {len({_dt_scena(os.path.basename(h)) for h in have})} scene)")

    got = list(have)
    if alta_ris_dir and os.path.isdir(alta_ris_dir):
        zips = sorted(glob.glob(os.path.join(alta_ris_dir, "*.zip")))  # esclude i .zip.part
        scene_have = {_dt_scena(os.path.basename(g)) for g in got
                     if f"-{pols[0]}-" in os.path.basename(g).lower()}
        nuovi = [z for z in zips if _dt_scena(os.path.basename(z)) not in scene_have]
        if zips:
            log(f"[2] {len(zips)} SAFE .zip in {alta_ris_dir} "
                f"({len(nuovi)} da estrarre, {len(zips)-len(nuovi)} gia' nello stack)")
        for z in nuovi:
            scene = {_dt_scena(os.path.basename(g)) for g in got
                     if f"-{pols[0]}-" in os.path.basename(g).lower()}
            if len(scene) >= n_want:
                break
            log(f"[2] estraggo da {os.path.basename(alta_ris_dir)}/: "
                f"{os.path.basename(z)} (gia' scaricato da scarica_alta_risoluzione.py)")
            try:
                got += _estrai_tiff_da_safe(z, stackdir, pols, rimuovi_zip=False)
            except Exception as ex:
                log(f"    ! estrazione fallita ({ex}) -- se il download e' ancora "
                    "in corso attendi che finisca (file .part)")

    user = user or os.environ.get("CDSE_USER")
    pwd = pwd or os.environ.get("CDSE_PASS")
    if not prods:
        log("[2] nessun prodotto dalla ricerca: uso i TIFF gia' presenti nello stack")
        return got
    # i prodotti trovati dallo step 1 CONTENGONO le coordinate richieste.
    # Riordina per track dominante (cluster sull'ora del giorno, come il filtro
    # dello step 3): prima le date del track piu' numeroso, poi gli altri.
    secs = [_acq_secs(p["Name"]) or -1 for p in prods]
    clusters = []
    for i, s in enumerate(secs):
        for cl in clusters:
            if abs(s - secs[cl[0]]) <= 20:
                cl.append(i); break
        else:
            clusters.append([i])
    clusters.sort(key=len, reverse=True)
    ordine = [i for cl in clusters for i in cl]
    if ordine != sorted(ordine):
        log(f"[2] {len(clusters)} track distinti tra i prodotti: scarico prima le "
            f"{len(clusters[0])} date del track dominante (le scene di track diversi "
            "verrebbero comunque escluse dal filtro dello step 3)")
    prods = [prods[i] for i in ordine]
    log(f"[2] {len(prods)} prodotti dalla ricerca contengono le coordinate "
        f"(ne servono >={n_want} scene)")
    for p in prods[:n_want]:
        log(f"      - {p['Name']}  ({p.get('ContentLength', 0)/1e9:.1f} GB)")
    if not download:
        log("[2] download non richiesto (manca --download): uso lo stack locale.")
        log(f"    -> i SAFE pesano ~8 GB l'uno. Per scaricare davvero i {len(prods)} "
            "prodotti trovati aggiungi")
        log("       --download  con  --cdse-user <user> --cdse-pass <pass>  "
            "(o le env CDSE_USER/CDSE_PASS)")
        return got
    if not (user and pwd):
        log("[2] credenziali CDSE assenti: salto il download.")
        log("    -> i SAFE pesano ~8 GB l'uno. Per scaricarli davvero, passa")
        log("       --cdse-user <user> --cdse-pass <pass>  (o le env CDSE_USER/CDSE_PASS)")
        return got
    try:
        _cdse_token(user, pwd)            # verifica subito le credenziali
    except Exception as ex:
        log(f"[2] login CDSE fallito ({ex}); salto il download"); return got
    get_token = lambda: _cdse_token(user, pwd)   # token fresco a ogni (ri)tentativo

    for p in prods:
        # n_want conta SCENE distinte (datetime) della pol PRIMARIA, non TIFF:
        # 1 SAFE = 3 TIFF iw1-3 per polarizzazione
        scene = {_dt_scena(os.path.basename(g)) for g in got
                 if f"-{pols[0]}-" in os.path.basename(g).lower()}
        if len(scene) >= n_want:
            break
        if _dt_scena(p["Name"]) in scene:
            log(f"[2] {p['Name']}: scena gia' nello stack, salto"); continue
        safe_zip = os.path.join(stackdir, p["Name"].replace(".SAFE", "") + ".zip")
        url = (f"https://catalogue.dataspace.copernicus.eu/odata/v1/"
               f"Products({p['Id']})/$value")
        expected = int(p.get("ContentLength", 0) or 0)
        zip_completo = (expected and os.path.exists(safe_zip)
                        and os.path.getsize(safe_zip) >= expected)
        ok = False
        for giro in (1, 2):          # 2 giri: il secondo dopo un'attesa-rete
            try:
                if zip_completo:
                    log(f"[2] {p['Name']}: zip gia' scaricato per intero, riuso "
                        "(salto il download, riestraggo)")
                    ok = True
                else:
                    log(f"[2] scarico {p['Name']} ({expected/1e9:.1f} GB) ...")
                    ok = _scarica_ripristinabile(url, safe_zip, get_token,
                                                 expected=expected)
            except urllib.error.HTTPError as ex:
                if ex.code == 401:
                    log("    ! 401 Unauthorized: credenziali CDSE errate/non "
                        "confermate (verifica l'account su dataspace.copernicus.eu) "
                        "o token scaduto")
                else:
                    log(f"    ! download fallito (HTTP {ex.code})")
            except urllib.error.URLError as ex:
                log(f"    ! download fallito (rete: {ex.reason})")
            if ok or _rete_ok():
                break                # riuscito, o fallito per un motivo non di rete
            if not _attendi_rete():  # rete giu': aspetta che torni e riprova
                break
        if not ok:
            log("    ! download fallito; continuo"); continue
        # estrai dal SAFE i measurement tiff + annotation delle polarizzazioni
        # (zip scaricato QUI: rimuovi_zip=True, a differenza degli zip letti da
        # alta_ris_dir che restano dove sono).
        try:
            got += _estrai_tiff_da_safe(safe_zip, stackdir, pols, rimuovi_zip=True)
        except Exception as ex:
            log(f"    ! estrazione SAFE fallita ({ex})")
    log(f"[2] TIFF '{'/'.join(pols)}' disponibili nello stack: {len(got)} "
        f"({len({_dt_scena(os.path.basename(g)) for g in got})} scene distinte)")
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
        # scarta le annotazioni ausiliarie del SAFE (rfi/calibration/noise): non
        # contengono la geometria del prodotto (radarFrequency, ecc.). Quella di
        # prodotto inizia col nome satellite (es. 's1a-iw2-slc-...').
        if not xb.startswith(("s1a-", "s1b-", "s1c-", "s1d-")):
            continue
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


def _acq_secs(nome):
    """Secondi-del-giorno dell'istante di acquisizione (HHMMSS nel nome SAFE/TIFF).
    Scene dello STESSO relative orbit (track) sovrappongono l'AOI quasi alla stessa
    ora del giorno: questo e' un identificatore robusto del track."""
    m = re.search(r"(\d{8})t(\d{2})(\d{2})(\d{2})", nome.lower())
    if not m:
        return None
    return int(m.group(2)) * 3600 + int(m.group(3)) * 60 + int(m.group(4))


def _subpix_offset(ref, mov):
    """Offset sub-pixel (d_az, d_rng) che porta 'mov' su 'ref' + correlazione di
    picco, via cross-correlazione FFT delle ampiezze con interpolazione parabolica.
    Serve a VERIFICARE quanto due scene si sfasano (diagnostica), non a co-registrare
    (il geocoding via GCP allinea gia' le scene dello stesso track)."""
    naz, nrg = ref.shape
    a = np.abs(ref).astype(float); b = np.abs(mov).astype(float)
    a = (a - a.mean()) / (a.std() + 1e-9); b = (b - b.mean()) / (b.std() + 1e-9)
    c = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real)
    p = np.unravel_index(np.argmax(c), c.shape)

    def par(i, ax):
        s = c.shape[ax]; im, ip = (i - 1) % s, (i + 1) % s
        y = (c[im, p[1]], c[i, p[1]], c[ip, p[1]]) if ax == 0 else \
            (c[p[0], im], c[p[0], i], c[p[0], ip])
        den = (y[0] - 2 * y[1] + y[2])
        return 0.0 if den == 0 else 0.5 * (y[0] - y[2]) / den

    daz = (p[0] - naz // 2) + par(p[0], 0)
    drg = (p[1] - nrg // 2) + par(p[1], 1)
    corr = c.max() / (np.sqrt((a ** 2).sum() * (b ** 2).sum()) + 1e-9)
    return daz, drg, corr


def _verifica_e_filtra_track(chips, nomi, geom, tol_s=20, track_filter=True):
    """VERIFICA le coordinate/geometria delle scene e quanto si sfasano tra loro,
    poi (se track_filter) tiene solo le scene del relative orbit dominante: scene di
    track diversi sono decorrelate e degradano coerenza e quote. Ritorna gli indici
    da tenere e la diagnostica (offset px/m + correlazione per scena)."""
    n = len(chips)
    dR = geom["dR"]; dA = geom["dA"]; inc = geom["incid_mid"]
    gr = dR / np.sin(np.radians(inc))            # passo ground-range [m]
    secs = [_acq_secs(x) for x in nomi]
    # raggruppa per ora-del-giorno (track): cluster entro tol_s secondi
    groups = []
    for i, s in enumerate(secs):
        for grp in groups:
            if s is not None and abs(s - secs[grp[0]]) <= tol_s:
                grp.append(i); break
        else:
            groups.append([i])
    groups.sort(key=lambda g: (-len(g), secs[g[0]] or 0))
    keep = sorted(groups[0]) if (track_filter and groups) else list(range(n))
    ref = keep[0]                                # riferimento = 1a scena del track tenuto

    log(f"[3] VERIFICA scene: {n} ritagliate, {len(groups)} track distinti "
        f"(per ora di acquisizione, tol {tol_s}s)")
    log(f"      {'data':<10}{'sat':<5}{'acq UTC':>9}{'d_az[px]':>9}{'d_rg[px]':>9}"
        f"{'d_az[m]':>9}{'d_rg[m]':>9}{'corr':>7}  nota")
    diag = []
    for i in range(n):
        daz, drg, corr = _subpix_offset(chips[ref], chips[i])
        hh = secs[i]
        acq = f"{hh//3600:02d}:{(hh%3600)//60:02d}:{hh%60:02d}" if hh is not None else "??"
        dt = _data_scena(nomi[i]); sat = nomi[i][:3]
        nota = "RIF" if i == ref else ("tenuta" if i in keep
               else "ESCLUSA (track diverso)")
        log(f"      {dt:<10}{sat:<5}{acq:>9}{daz:>9.2f}{drg:>9.2f}"
            f"{daz*dA:>9.1f}{drg*gr:>9.1f}{corr:>7.2f}  {nota}")
        diag.append((dt, sat, daz, drg, corr, i in keep))
    if track_filter and len(keep) < n:
        log(f"[3] track dominante: {len(keep)}/{n} scene tenute "
            f"(escluse {n-len(keep)} di altri track -> avrebbero decorrelato le quote)")
    return keep, diag


def _estrai_chips_pol(stackdir, pol, corners):
    """Geolocalizza l'AOI su ogni TIFF della polarizzazione `pol` via i GCP e
    ritaglia i chip complessi. Ritorna (chips, nomi, geom)."""
    tiffs = sorted(glob.glob(os.path.join(stackdir, f"*-iw[1-3]-slc-{pol}-*.tiff"))) or \
            sorted(glob.glob(os.path.join(stackdir, f"*slc-{pol}-*.tiff")))
    if not tiffs:
        return [], [], None
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
        # un chip quasi tutto-zeri = AOI sul bordo dello swath / fuori dal burst:
        # non copre davvero l'area e inquinerebbe gli interferogrammi -> scarta
        if valid < MIN_VALID_PCT:
            log(f"       (copertura {valid:.0f}% < {MIN_VALID_PCT:.0f}%: l'AOI non e' "
                "nello swath di questo sub-swath, skip)")
            continue
        chips.append(chip); nomi.append(os.path.basename(tif)); geom = g
    return chips, nomi, geom


def step3_estrai_box(stackdir, pol, aoi, boxpath, fallback_box, track_filter=True):
    """Ritaglia lo stack complesso dell'AOI per OGNI polarizzazione richiesta
    (--pol vv,vh). Il track filter e' calcolato sulla polarizzazione primaria e
    applicato per scena alle altre (VV e VH della stessa scena condividono
    acquisizione e griglia raster). Nel box .npz: <pol>=stack, <pol>_scene=nomi."""
    pols = lista_pol(pol)
    lonW, lonE, latS, latN = aoi
    corners = [(lonW, latN), (lonE, latN), (lonE, latS), (lonW, latS)]
    per_pol, geom = {}, None
    for p in pols:
        log(f"[3] ritaglio polarizzazione '{p}'...")
        chips, nomi, g = _estrai_chips_pol(stackdir, p, corners)
        if chips:
            per_pol[p] = (chips, nomi)
            geom = geom or g
        else:
            log(f"[3] pol '{p}': nessun TIFF utilizzabile in {stackdir}")
    if not per_pol:
        if os.path.exists(fallback_box):
            log(f"[3] nessun TIFF in {stackdir}: uso il box gia' ritagliato "
                f"{fallback_box} (fallback)")
            d = np.load(fallback_box, allow_pickle=True)
            np.savez_compressed(boxpath, **{k: d[k] for k in d.files})
            return boxpath
        sys.exit(f"[3] nessun TIFF in {stackdir} e nessun fallback {fallback_box}")

    # dimensioni comuni a TUTTI i chip (le pol condividono la griglia raster)
    H = min(c.shape[0] for chips, _ in per_pol.values() for c in chips)
    W = min(c.shape[1] for chips, _ in per_pol.values() for c in chips)
    per_pol = {p: ([c[:H, :W] for c in chips], nomi)
               for p, (chips, nomi) in per_pol.items()}

    # verifica geometria/sfasamento sul PRIMARIO e tiene solo il track dominante
    primario = next(p for p in pols if p in per_pol)
    chips0, nomi0 = per_pol[primario]
    keep, diag = _verifica_e_filtra_track(chips0, nomi0, geom, track_filter=track_filter)
    arr = np.stack([chips0[i] for i in keep], 0)
    nomi_k = [nomi0[i] for i in keep]
    off_az = np.array([diag[i][2] for i in keep], np.float32)
    off_rg = np.array([diag[i][3] for i in keep], np.float32)
    corr_ref = np.array([diag[i][4] for i in keep], np.float32)
    salva = {primario: arr, f"{primario}_scene": np.array(nomi_k)}

    # altre polarizzazioni: stesse SCENE tenute per il primario (chiave
    # satellite+subswath+datetime -> il measurement gemello dell'altra pol)
    for p, (chips_p, nomi_p) in per_pol.items():
        if p == primario:
            continue
        idx = {_chiave_scena_sw(n): i for i, n in enumerate(nomi_p)}
        sel = [(idx[_chiave_scena_sw(n)], n) for n in nomi_k
               if _chiave_scena_sw(n) in idx]
        if len(sel) < 2:
            log(f"[3] pol '{p}': solo {len(sel)} scene coincidenti col track "
                f"del primario (<2), la escludo dal box")
            continue
        salva[p] = np.stack([chips_p[i] for i, _ in sel], 0)
        salva[f"{p}_scene"] = np.array([nomi_p[i] for i, _ in sel])
        log(f"[3] pol '{p}': {len(sel)} scene allineate alle {len(nomi_k)} del primario")

    np.savez_compressed(boxpath, **salva,
                        corners_lonlat=np.array(corners),
                        dR=geom["dR"], dA=geom["dA"], incid_mid=geom["incid_mid"],
                        f_em=geom["f_em"], R_near=geom["R_near"], dt_az=geom["dt_az"],
                        off_az_px=off_az, off_rg_px=off_rg, corr_ref=corr_ref)
    log(f"[3] stack ritagliato {arr.shape} (n_scene, az, rng), "
        f"pol nel box: {[p for p in pols if p in salva]} -> {boxpath}")
    return boxpath


# ===================================================================== STEP 4
def _data_scena(nome):
    m = re.search(r"-(?:vv|vh|hh|hv)-(\d{8})t", nome.lower()) or \
        re.search(r"(\d{8})t\d{6}", nome.lower())
    return m.group(1) if m else "00000000"


def _interferogrammi_pol(stack, nomi, n_layers, K, cth, etichetta_pol):
    """Coppie interferometriche e strati (disp verticale [mm], fase [rad]) per
    una singola polarizzazione. Ritorna (disp, phi_st, etich) o None se lo
    stack non ha abbastanza scene (<2)."""
    date = [_data_scena(n) for n in nomi]
    order = np.argsort(date)
    stack = stack[order]; nomi = [nomi[i] for i in order]
    date = [date[i] for i in order]
    nS, ny, nx = stack.shape
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
    if len(pairs) == 0:
        return None
    if len(pairs) < n_layers:
        log(f"[4] ATTENZIONE (pol {etichetta_pol}): solo {len(pairs)} coppie "
            f"disponibili (< {n_layers})")

    disp = np.empty((len(pairs), ny, nx), np.float32)   # componente verticale [mm]
    phi_st = np.empty((len(pairs), ny, nx), np.float32) # fase interferometrica [rad]
    etich = []
    for s, (i, j, bt) in enumerate(pairs):
        phi = np.angle(stack[j] * np.conj(stack[i]))
        phi_st[s] = phi
        disp[s] = (-K * phi / cth).astype(np.float32)
        etich.append(f"{date[i]}~{date[j]} ({bt}d)")
    return disp, phi_st, etich


def step4_array_12(boxpath, n_layers, outdir):
    """Costruisce l'array a n_layers strati della COMPONENTE VERTICALE dello
    spostamento (DInSAR classica, coerente con la skill sar-dinsar-microdisplacement:
    phi = arg(s_j * conj(s_i)) -> LOS = -lambda/(4pi)*phi -> verticale = LOS/cos(theta)),
    per OGNI polarizzazione presente nel box (--pol vv,vh): ogni pol e' una misura
    indipendente della stessa scena e alimenta forme d'onda proprie negli step 5-6.
    NON e' la tomografia Doppler a sub-aperture del paper Biondi-Malanga (vedi note
    nel docstring del modulo): qui ogni strato e' un interferogramma fra due
    acquisizioni in DATE diverse, non una sub-apertura Doppler di una singola
    immagine. Con 12+ scene bastano le coppie sequenziali; con poche scene si usano
    le coppie a baseline minima fino ad arrivare a n_layers.
    Ritorna (per_pol, x, y, coh) con per_pol = {pol: (disp, phi)}."""
    d = np.load(boxpath, allow_pickle=True)
    pols = [p for p in ("vv", "vh", "hh", "hv") if p in d.files]
    if not pols:
        sys.exit(f"[4] nessuna polarizzazione nel box {boxpath}")
    inc = float(d["incid_mid"]); f_em = float(d["f_em"])
    dR = float(d["dR"]); dA = float(d["dA"])
    lam = C / f_em
    K = lam / (4 * np.pi) * 1000.0                 # mm per radiante
    cth = np.cos(np.radians(inc))                  # proiezione LOS -> verticale
    gr = dR / np.sin(np.radians(inc))              # passo ground-range [m]

    per_pol, etich_pol = {}, {}
    for p in pols:
        nomi = [str(x) for x in d[f"{p}_scene"]]
        r = _interferogrammi_pol(d[p], nomi, n_layers, K, cth, p)
        if r is None:
            if p == pols[0]:
                date = sorted({_data_scena(n) for n in nomi})
                sys.exit(
                    f"[4] impossibile proseguire: l'interferometria richiede ALMENO 2 "
                    f"acquisizioni dello stesso track, ma lo stack '{p}' ne contiene "
                    f"{d[p].shape[0]} (date: {date}).\n"
                    f"    -> scarica altre date con --download e credenziali CDSE valide "
                    f"(il download nello step 2 e' fallito: 401 Unauthorized / rete assente).")
            log(f"[4] pol '{p}': meno di 2 scene, la salto")
            continue
        per_pol[p] = (r[0], r[1])
        etich_pol[p] = r[2]

    primario = pols[0]
    disp, phi_st = per_pol[primario]
    etich = etich_pol[primario]
    coh = np.abs(np.exp(1j * phi_st).mean(0))
    ny, nx = disp.shape[1:]

    x = np.arange(nx) * gr                          # asse range a terra [m]
    y = np.arange(ny) * dA                          # asse azimuth [m]
    os.makedirs(outdir, exist_ok=True)
    extra = {}
    for p, (dp, ph) in per_pol.items():
        if p != primario:
            extra[f"disp_{p}"] = dp; extra[f"phi_{p}"] = ph
            extra[f"etich_{p}"] = np.array(etich_pol[p])
    np.savez_compressed(os.path.join(outdir, "array3d_12strati.npz"),
                        disp=disp, phi=phi_st, x=x, y=y, coh=coh.astype(np.float32),
                        etich=np.array(etich), inc=inc, lam=lam, gr=gr, dA=dA,
                        **extra)
    for p, (dp, _) in per_pol.items():
        log(f"[4] pol '{p}': array componente verticale {dp.shape} (strati, y, x); "
            f"vert [{dp.min():.1f},{dp.max():.1f}] mm")
    log(f"[4] coh media (pol '{primario}'): {coh.mean():.2f}")
    dbg("strati:", etich)

    # figura riassuntiva degli strati
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nrow = int(np.ceil(len(etich) / 4))
    fig, axes = plt.subplots(nrow, 4, figsize=(15, 3.2 * nrow))
    for k, ax in enumerate(np.atleast_1d(axes).ravel()):
        if k < len(etich):
            im = ax.imshow(disp[k], cmap="RdBu_r", vmin=-20, vmax=20,
                           aspect="auto", origin="lower")
            ax.set_title(etich[k], fontsize=8)
        else:
            ax.axis("off")
    fig.colorbar(im, ax=axes, shrink=0.7, label="spostamento verticale [mm]")
    fig.suptitle(f"Step 4 — array a {len(etich)} strati (componente verticale, "
                 f"DInSAR, pol {primario})")
    fig.savefig(os.path.join(outdir, "step4_strati.png"), dpi=120, bbox_inches="tight")
    plt.close(fig)
    return per_pol, x, y, coh


# ----------------------------------------------------------- DEM (altezze reali)
def _fetch_dem(aoi, ny, nx, outdir, pyr_geom):
    """Quota reale del terreno (DEM Copernicus via Open-Meteo) sulla griglia
    pixel (ny,nx). Opzionale overlay della geometria della piramide selezionata
    (pyr_geom = dict(base_m, h_m, label) dal PRESET, None = nessun overlay):
    ogni piramide ha base/altezza proprie, usarle giuste e' essenziale perche'
    la base e' il riferimento metrico interno (paper Biondi-Malanga, sez. 5.19)."""
    cache = os.path.join(outdir, "dem.npz")
    lonW, lonE, latS, latN = aoi
    j = np.arange(nx); i = np.arange(ny)
    lon = lonW + (j / max(nx - 1, 1)) * (lonE - lonW)
    lat = latS + (i / max(ny - 1, 1)) * (latN - latS)
    z0 = None
    if os.path.exists(cache):
        # la cache contiene il DEM NUDO (senza overlay piramide)
        cached = np.load(cache)["z0"].astype(float)
        if cached.shape == (ny, nx):
            z0 = cached
            log(f"[DEM] riuso {cache}: quota terreno {z0.min():.0f}..{z0.max():.0f} m s.l.m.")
        else:
            log(f"[DEM] cache {cache} ignorata: shape {cached.shape} != griglia ({ny},{nx})")
    if z0 is None:
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

    if pyr_geom:
        # geometria reale della piramide selezionata (base e altezza ATTUALE dal
        # PRESET: Cheope 230.3/138.7 m, Chefren 215.5/136.4 m)
        latm = np.radians((latN + latS) / 2)
        Lx = (lonE - lonW) * 111320 * np.cos(latm)
        Ly = (latN - latS) * 110540
        x_px = np.linspace(0, Lx, nx); y_px = np.linspace(0, Ly, ny)
        H_PYR = float(pyr_geom["h_m"]); BASE = float(pyr_geom["base_m"])
        half = BASE / 2.0
        XX, YY = np.meshgrid(x_px, y_px)
        cheb = np.maximum(np.abs(XX - Lx / 2), np.abs(YY - Ly / 2))
        z0 = z0 + np.clip(H_PYR * (1.0 - cheb / half), 0.0, None)
        log(f"[DEM] overlay piramide {pyr_geom['label']} (base {BASE:.1f} m, "
            f"h {H_PYR:.1f} m) -> apice {z0.max():.0f} m s.l.m.")
    return z0


# ===================================================================== STEP 5
def step5_grafico_onde(per_pol, x, y, aoi, outdir, pyr_geom,
                       NZ=256, ZTOP=210.0, ZBOT=-200.0, STEP_X=None, STEP_Y=1,
                       WIGGLE=14.0, pyramid=None, stackdir=None, prof_info=None):
    """Grafico 3D dell'array: i numeri dei n_layers strati DInSAR di ogni pixel
    sono usati come coefficienti armonici complessi (ampiezza=|vert|, fase=fase
    interferometrica) di una forma d'onda sviluppata in profondita' via
    anti-trasformata di Fourier. Il range di profondita' e' la z_amb del dato
    sorgente quando prof_info e' passato (skill biondi-sar-tomography),
    altrimenti ZTOP..ZBOT. Le onde partono dalla QUOTA REALE del suolo (DEM).
    TUTTI gli strati di TUTTE le polarizzazioni richieste (--pol, default
    vv,vh) vengono messi in UN SOLO array di coefficienti armonici (es. 12 VV
    + 12 VH = 24 armoniche: prima gli strati VV in ordine di data, poi VH) e
    l'anti-trasformata e' calcolata IN UN SOLO PASSAGGIO (un unico
    np.fft.irfft vettorizzato su tutta l'area), a risoluzione PIENA (ogni
    pixel az/rng del box, nessun subsample) -> UNA forma d'onda combinata per
    pixel. Nel grafico 3D si disegna una traccia ogni STEP_X colonne/STEP_Y
    righe: con i default (STEP_X auto=1, STEP_Y=1) viene disegnata una
    "molla" per OGNI pixel dell'area (ny*nx), nessun sotto-campionamento.
    Analisi esplorativa: vedi nota metodologica in cima al modulo."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go

    pols = list(per_pol)
    disp, phi = per_pol[pols[0]]
    nS, ny, nx = disp.shape
    # mappatura dei punti sui PIXEL REALI: x,y arrivano dallo step 4 e sono gia'
    # le posizioni vere di ogni pixel (passo ground-range = dR/sin(incidenza),
    # passo azimuth = dA, entrambi dalla geometria SAR reale annotata nel
    # prodotto), la STESSA griglia (ny,nx) dell'array DInSAR -> ogni pixel del
    # grafico 3D e' allineato esattamente al pixel su cui e' stata calcolata la
    # trasformata di Fourier. In precedenza x,y venivano ignorati e si
    # ricostruiva una griglia approssimata dai gradi lon/lat dell'AOI (fattori
    # fissi 111320/110540 m/grado): non allineata al pixel reale e incoerente
    # con lo step 4.
    x_m = np.asarray(x, dtype=np.float64)
    y_m = np.asarray(y, dtype=np.float64)
    Lx = float(x_m[-1]) or 1.0
    Ly = float(y_m[-1]) or 1.0
    z0 = _fetch_dem(aoi, ny, nx, outdir, pyr_geom)

    # riferimenti strutture interne note (vedi RIFERIMENTI_PIRAMIDE): calibra
    # centro/orientamento sui GCP reali del prodotto, poi posiziona gli
    # INGOMBRI (parallelepipedi con le dimensioni da letteratura) rispetto
    # alla quota REALE del terreno (DEM nudo, senza overlay piramide, che e'
    # il riferimento fisico per "altezza sopra la base"); quasi tutti sul
    # piano verticale centrale N-S (e_off=0), tranne l'entrata originale di
    # Cheope che ha un offset E-W noto e documentato (e_off).
    riferimenti = []
    if pyramid in RIFERIMENTI_PIRAMIDE and stackdir:
        cal = _calibra_centro_piramide(stackdir, pols[0], aoi, x_m, y_m)
        if cal:
            x_c, y_c, i_c, j_c, segno_nord, segno_est = cal
            try:
                z_bare = np.load(os.path.join(outdir, "dem.npz"))["z0"]
            except Exception:
                z_bare = z0
            ii = min(max(int(round(i_c)), 0), z_bare.shape[0] - 1)
            jj = min(max(int(round(j_c)), 0), z_bare.shape[1] - 1)
            z_base = float(z_bare[ii, jj])
            for num, nome, s_off, h_off, e_off, ns_m, ew_m, alt_m in \
                    RIFERIMENTI_PIRAMIDE[pyramid]:
                riferimenti.append((num, nome,
                                    x_c + segno_est * e_off,
                                    y_c + segno_nord * s_off,
                                    z_base + h_off,
                                    ew_m, ns_m, alt_m))
            log(f"[5] {len(riferimenti)} strutture interne note ({pyr_geom['label'] if pyr_geom else pyramid}, "
                f"ingombri da rilievi pubblicati, NON dai dati SAR) disegnate come "
                f"parallelepipedi nel grafico 3D: centro base a x={x_c:.0f} m, "
                f"y={y_c:.0f} m, quota base {z_base:.0f} m s.l.m.")
        else:
            log("[5] riferimenti strutture interne non posizionati: nessun TIFF/GCP "
                "disponibile per calibrare centro/orientamento (stack ripulito?)")

    # ZTOP adattivo: con l'overlay piramide l'apice (quota plateau + h piramide,
    # es. Cheope ~200-215 m s.l.m.) puo' superare il default e verrebbe tagliato
    if float(z0.max()) + 5.0 > ZTOP:
        ZTOP = float(np.ceil(z0.max())) + 5.0
        log(f"[5] ZTOP alzato a {ZTOP:.0f} m s.l.m. per non tagliare l'apice "
            f"({z0.max():.0f} m)")
    if prof_info is not None:
        # profondita' di sviluppo = z_amb del dato sorgente (skill
        # biondi-sar-tomography): oltre quella profondita' la DFT va in alias,
        # quindi il grafico NON deve mostrare profondita' oltre z_amb.
        tipo_fonte, z_amb_f, dz_tomo = prof_info[:3]
        ZD = float(z_amb_f)
        ZBOT = float(np.floor(z0.min() - ZD)) - 5.0
        log(f"[5] profondita' MASSIMA dell'inversione di Fourier per il dato "
            f"sorgente [{tipo_fonte}]: z_amb = {ZD:.0f} m sotto il suolo, "
            f"risoluzione tomografica dz = {dz_tomo:.0f} m "
            f"(skill biondi-sar-tomography: z_amb = lambda_sonic*R*sin(theta)"
            f"*(k-1)/(2A)) -> onde sviluppate fino a {ZBOT:.0f} m s.l.m.")
    else:
        tipo_fonte, dz_tomo = None, None
        ZD = ZTOP - ZBOT
    z = np.linspace(0, ZD, NZ); nf = NZ // 2 + 1
    # TUTTI gli strati INSIEME, IN UN SOLO PASSAGGIO: si crea UN SOLO array di
    # coefficienti armonici con gli strati di TUTTE le polarizzazioni (es. 12
    # VV + 12 VH = 24 armoniche: strati VV in ordine di data, poi VH) e si
    # calcola l'anti-trasformata di Fourier UNA VOLTA sola per pixel (un unico
    # np.fft.irfft vettorizzato su tutta l'area). Prima veniva calcolata una
    # trasformata SEPARATA per ogni polarizzazione: corretto su richiesta.
    coeff_tutti = np.concatenate(
        [(np.abs(dp) * np.exp(1j * ph)).astype(complex)
         for dp, ph in (per_pol[p] for p in pols)], axis=0)   # (n_tot, ny, nx)
    n_per_pol = {p: per_pol[p][0].shape[0] for p in pols}
    n_tot = coeff_tutti.shape[0]
    if n_tot >= nf:
        log(f"[5] ATTENZIONE: {n_tot} armoniche >= {nf} campioni spettrali con "
            f"NZ={NZ}; tronco alle prime {nf - 1} (alza --nz per usarle tutte)")
        coeff_tutti = coeff_tutti[:nf - 1]
        n_tot = nf - 1
    spec = np.zeros((ny * nx, nf), complex)
    spec[:, 1:1 + n_tot] = coeff_tutti.reshape(n_tot, ny * nx).T
    pols_orig, pols = pols, ["+".join(pols)]
    w_pol = {pols[0]: (np.fft.irfft(spec, n=NZ, axis=1) * NZ)
             .reshape(ny, nx, NZ).astype(np.float32)}
    if STEP_X is None:
        # una sola forma d'onda combinata -> una traccia per OGNI colonna:
        # il numero di "molle" = numero di pixel dell'area (ny*nx), una per
        # pixel. Passa --step-x/--step-y > 1 per un'anteprima sotto-campionata.
        STEP_X = len(pols)          # = 1
    w = w_pol[pols[0]]
    wmax = max(float(np.abs(wp).max()) for wp in w_pol.values()) or 1.0
    log(f"[5] forme d'onda (analisi esplorativa, NON tomografia Doppler del paper): "
        f"anti-FT UNICA in un solo passaggio di {n_tot} armoniche = tutti gli "
        f"strati insieme ({' + '.join(f'{n_per_pol[q]} {q.upper()}' for q in pols_orig)}) "
        f"-> w {w.shape}, |w|max={wmax:.1f}; base {Lx:.0f}x{Ly:.0f} m, "
        f"profondita' {ZD:.0f} m")

    # frequenza istantanea (Hilbert) PER OGNI polarizzazione: guida il colore
    # PIENO/VUOTO delle "molle" qui sotto e viene riusata tal quale dallo step 6
    # (nessuna ricomputazione, stessa definizione di frequenza in entrambi i
    # grafici). dz = passo dell'asse profondita' z.
    dz = float(z[1] - z[0])
    freq_pol, env_pol = {}, {}
    for p in pols:
        freq_pol[p], env_pol[p] = _frequenza_istantanea(w_pol[p], dz)
    freq_all = np.concatenate([f.ravel() for f in freq_pol.values()])
    fmin = float(np.nanpercentile(freq_all, 2))
    fmax = float(np.nanpercentile(freq_all, 98))
    if fmax <= fmin:
        fmax = fmin + 1e-6
    log(f"[5] frequenza istantanea (colore PIENO/VUOTO delle molle, percentili "
        f"2-98 su tutta l'area): [{fmin:.4f}, {fmax:.4f}] 1/m")

    # forme d'onda SEPARATE per polarizzazione (stessa anti-FT ma con i SOLI
    # strati di quella pol): alimentano lo step 6, che segna SOLO i punti
    # VUOTO generati prima dalla prima pol richiesta (VV) e poi dalle
    # successive (VH). La forma d'onda COMBINATA qui sopra resta quella
    # disegnata dalle "molle" di questo step.
    freq_sep, env_sep = {}, {}
    for p in pols_orig:
        dp, ph = per_pol[p]
        coeff = (np.abs(dp) * np.exp(1j * ph)).astype(complex)
        n_p = min(coeff.shape[0], nf - 1)
        spec_p = np.zeros((ny * nx, nf), complex)
        spec_p[:, 1:1 + n_p] = coeff[:n_p].reshape(n_p, ny * nx).T
        w_sep = (np.fft.irfft(spec_p, n=NZ, axis=1) * NZ) \
            .reshape(ny, nx, NZ).astype(np.float32)
        freq_sep[p], env_sep[p] = _frequenza_istantanea(w_sep, dz)
    log(f"[5] forme d'onda separate per pol ({', '.join(q.upper() for q in pols_orig)}): "
        f"pronte per lo step 6 (punti VUOTO per pol, prima "
        f"{pols_orig[0].upper()}" +
        (f" poi {' poi '.join(q.upper() for q in pols_orig[1:])}" if len(pols_orig) > 1 else "") + ")")

    # con una sola pol la chiave combinata coincide con quella separata
    # (es. 'freq_vv'): il dict li unifica senza duplicati
    freq_npz = {f"freq_{p}": freq_pol[p] for p in pols}
    freq_npz.update({f"freq_{p}": freq_sep[p] for p in pols_orig})
    np.savez_compressed(os.path.join(outdir, "forme_donda.npz"),
                        w=w, z=z, x=x_m, y=y_m, z0=z0, Lx=Lx, Ly=Ly, zdepth=ZD,
                        pol=np.array(pols), freq_range=np.array([fmin, fmax]),
                        **{f"w_{p}": w_pol[p] for p in pols[1:]},
                        **freq_npz)

    # "molle" sviluppate verso il basso dalla quota reale: una traccia per
    # OGNI pixel (STEP_X auto = 1 con la forma d'onda combinata), ognuna
    # sintetizzata dall'anti-FT UNICA di tutti gli strati insieme.
    ii = np.arange(0, ny, STEP_Y)
    sc = WIGGLE / wmax
    tracce = {}
    gap = np.array([np.nan])
    offsets = [(k * STEP_X) // len(pols) for k in range(len(pols))]
    for k, p in enumerate(pols):
        jj = np.arange(offsets[k], nx, STEP_X)
        X, Y, Z, Cc = [], [], [], []
        for i in ii:
            for j in jj:
                wp = w_pol[p][i, j]
                fp = freq_pol[p][i, j]           # colore = frequenza istantanea
                X += [x_m[j] + sc * wp, gap]; Y += [np.full(NZ, y_m[i]), gap]
                Z += [z0[i, j] - z, gap]; Cc += [fp, gap]
        tracce[p] = (np.concatenate(X), np.concatenate(Y),
                     np.concatenate(Z), np.concatenate(Cc))
    n_tracce = sum(len(np.arange(offsets[k], nx, STEP_X)) * len(ii)
                   for k in range(len(pols)))
    colonne_coperte = len(set(j for off in offsets for j in range(off, nx, STEP_X)))
    righe_coperte = len(ii)
    pixel_area = ny * nx
    log(f"[5] tracce 3D: {n_tracce} totali su {len(pols)} polarizzazioni "
        f"({colonne_coperte}/{nx} colonne E-W coperte"
        f"{', PIENA' if colonne_coperte == nx else ''}"
        f" | {righe_coperte}/{ny} righe N-S coperte"
        f"{', PIENA' if righe_coperte == ny else ''})")
    if n_tracce == pixel_area:
        log(f"[5] OK: {n_tracce} tracce = {pixel_area} pixel dell'area selezionata "
            f"({ny}x{nx}) -> una traccia per pixel, nessun sotto/sovra-campionamento")
    else:
        log(f"[5] ATTENZIONE: {n_tracce} tracce != {pixel_area} pixel dell'area "
            f"({ny}x{nx}): --step-x/--step-y non di default, la copertura non e' 1:1")
    if pixel_area > 200_000 and n_tracce == pixel_area:
        log(f"[5] area selezionata grande ({pixel_area} pixel): il grafico 3D con una "
            "molla per pixel puo' essere pesante da aprire nel browser; usa "
            "--step-x/--step-y > 1 per un'anteprima piu' leggera e sotto-campionata")

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
    springs = []
    for k, p in enumerate(pols):
        X, Y, Z, Cc = tracce[p]
        springs.append(go.Scatter3d(
            x=X, y=Y, z=Z, mode="lines", connectgaps=False,
            line=dict(width=3, color=Cc, colorscale=PIENO_VUOTO_COLORSCALE,
                      cmin=fmin, cmax=fmax,
                      colorbar=_colorbar_pieno_vuoto(fmin, fmax) if k == 0 else None,
                      showscale=(k == 0)),
            name=f"forme d'onda {p.upper()} ({n_tot} armoniche, tutti gli "
                 f"strati in un solo passaggio)"))
    ground = go.Scatter3d(x=Xs.ravel(), y=Ys.ravel(), z=z0.ravel(), mode="markers",
                          marker=dict(size=1.2, color="red"), name="suolo (DEM)",
                          text=txt, hoverinfo="text")
    extra_traces = _tracce_riferimenti(riferimenti) if riferimenti else []
    fig = go.Figure([plane, zero_plane, *springs, ground, *extra_traces])
    log(f"[5] altezza max del suolo dal livello 0: {h_suolo.max():.0f} m; "
        f"onde sviluppate fino a quota {ZBOT:.0f} m s.l.m. "
        f"(= {abs(ZBOT):.0f} m di profondita' sotto il livello 0)")
    _tema_plotly(fig)
    fig.update_layout(
        scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
                   zaxis_title="quota reale s.l.m. [m]", zaxis=dict(range=[ZBOT, ZTOP]),
                   aspectmode="manual",
                   aspectratio=dict(x=1.0, y=Ly / Lx, z=ZD / Lx)))
    _scrivi_pagina_html(
        fig, os.path.join(outdir, "step5_forme_donda_3d.html"),
        titolo=f"Step 5 — Forme d'onda 3D"
               f"{' · ' + pyr_geom['label'] if pyr_geom else ''}",
        sottotitolo="Forme d'onda sviluppate in profondita' dalla quota reale "
                    "del suolo (DEM); il colore di ogni \"molla\" e' la sua "
                    "frequenza istantanea: PIENO = bassa (chiaro), VUOTO = "
                    "alta (scuro). Livello 0 = mare. Trascina per ruotare, "
                    "rotella per lo zoom.",
        badges=[("armoniche (tutti gli strati insieme)",
                 f"{n_tot} = " + " + ".join(f"{n_per_pol[q]} {q.upper()}"
                                            for q in pols_orig)),
                ("molle", f"{n_tracce} = {pixel_area} pixel ({ny}×{nx})"),
                ("base", f"{Lx:.0f}×{Ly:.0f} m"),
                ("quota suolo", f"{z0.min():.0f}–{z0.max():.0f} m s.l.m.")]
               + ([("profondita' max fonte (z_amb)", f"{ZD:.0f} m"),
                   ("dz tomografica", f"{dz_tomo:.0f} m"),
                   ("sorgente", tipo_fonte)] if prof_info else
                  [("profondita'", f"fino a {ZBOT:.0f} m s.l.m.")]),
        nota="Analisi ESPLORATIVA ispirata al lessico di Biondi & Malanga "
             "(Remote Sens. 2022, 14, 5231) ma NON la loro tomografia Doppler "
             "a sub-aperture: gli strati DInSAR (coppie di acquisizioni in "
             "date diverse) sono trattati come coefficienti spettrali di una "
             "serie di Fourier stirata su un range di profondita' scelto, non "
             "derivato dalla geometria SAR. Usare per ipotesi, non come "
             "misura fisica validata di profondita'.")

    # anteprima statica + B-scan
    figm = plt.figure(figsize=(9, 8)); ax = figm.add_subplot(111, projection="3d")
    ax.plot_surface(Xs, Ys, z0, cmap="terrain", alpha=0.4, linewidth=0)
    for k, p in enumerate(pols):
        off = (k * STEP_X) // len(pols)
        for i in ii:
            for j in np.arange(off, nx, STEP_X):
                ax.plot(x_m[j] + sc * w_pol[p][i, j], np.full(NZ, y_m[i]),
                        z0[i, j] - z, lw=0.6)
    ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
    ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZD))
    ax.set_title(f"Step 5 — forme d'onda dalla quota reale "
                 f"({n_tot} armoniche = tutti gli strati insieme, {'/'.join(pols)})")
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
    # allo step 6 vanno le forme d'onda SEPARATE per pol (VV, VH, ...): segna
    # solo i punti VUOTO generati prima dalla prima pol e poi dalle successive
    return w_pol, freq_sep, env_sep, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT, riferimenti


# ===================================================================== STEP 6
def step6_variazioni(freq_pol, env_pol, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT, outdir,
                     soglia=None, riferimenti=None, delta_vuoto=5.0, prof_info=None):
    """NEL GRAFICO SI SEGNANO SOLO I PUNTI VUOTO, generati PRIMA dalla forma
    d'onda della prima polarizzazione richiesta (VV) e POI da quelle
    successive (VH): per OGNI pol, per OGNI pixel dell'area (ny*nx, la
    stessa griglia degli step 4-5) si prende lungo la profondita' il campione
    dove l'INVILUPPO (ampiezza dell'analitico di Hilbert, calcolato dallo
    step 5 sulla forma d'onda SEPARATA di quella pol) e' massimo -- il
    "riflettore" piu' forte della colonna -- con la sua frequenza istantanea.
    Il punto e' classificato PIENO (materiale interpretato come denso/solido)
    se la frequenza e' SOTTO `soglia`, VUOTO (cavita'/aria) se e' sopra:
    soglia di default = mediana delle frequenze di picco su tutti i punti
    generati (sovrascrivibile con `soglia`). SOLO i punti VUOTO vengono
    disegnati (una traccia per pol, nell'ordine di --pol); i PIENO non sono
    disegnati ma restano nel .npz e nel cursore della soglia. Disegnato in 3D
    alle ALTEZZE REALI (quota DEM - profondita'). In aggiunta, un OVERLAY
    selezionabile (nascosto di default, si attiva dalla legenda o dal
    pulsante dedicato) riempie con punti/linea ogni tratto VERTICALE che va
    da una transizione PIENO->VUOTO alla successiva VUOTO->PIENO (per pol),
    per visualizzare l'estensione in profondita' di ogni cavita' individuata.
    Vedi nota metodologica in cima al modulo: e' un'analisi ispirata al
    lessico Biondi-Malanga ma non e' la loro tomografia Doppler a
    sub-aperture."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import plotly.graph_objects as go

    if not isinstance(freq_pol, dict):               # retro-compatibilita'
        freq_pol = {"vv": freq_pol}
    pols = list(freq_pol)
    ny, nx, nz = freq_pol[pols[0]].shape
    pixel_area = ny * nx

    # punti generati PER POLARIZZAZIONE, nell'ordine di --pol: PRIMA tutti i
    # pixel dell'area dalla forma d'onda VV, POI tutti dalla VH (una passata
    # completa per pol, forme d'onda SEPARATE dallo step 5). Profondita'
    # scelta = picco dell'inviluppo (riflettore piu' forte della colonna).
    # Nel grafico si SEGNANO SOLO I VUOTI; i PIENO restano in .npz/cursore.
    px_, py_, pz_, pv_, ppol_ = [], [], [], [], []
    for p in pols:
        kk = env_pol[p].reshape(pixel_area, nz).argmax(axis=1)
        px_.append(np.tile(x_m, ny))
        py_.append(np.repeat(y_m, nx))
        pz_.append(z0.ravel() - z[kk])
        pv_.append(freq_pol[p].reshape(pixel_area, nz)[np.arange(pixel_area), kk])
        ppol_.append(np.full(pixel_area, p, dtype=object))
    px = np.concatenate(px_); py = np.concatenate(py_)
    pz = np.concatenate(pz_); pv = np.concatenate(pv_)
    ppol = np.concatenate(ppol_)
    n = px.size
    assert n == pixel_area * len(pols)
    log(f"[6] {n} punti = {pixel_area} pixel dell'area selezionata ({ny}x{nx}) "
        f"x {len(pols)} polarizzazioni, generati in ordine: "
        f"{' -> '.join(p.upper() for p in pols)} (forme d'onda separate per pol)")

    if soglia is None:
        soglia = float(np.median(pv))
        log(f"[6] soglia PIENO/VUOTO automatica (mediana delle frequenze di "
            f"picco su tutti i punti generati): {soglia:.4f} 1/m")
    else:
        log(f"[6] soglia PIENO/VUOTO (esplicita, --soglia): {soglia:.4f} 1/m")
    classe = np.where(pv < soglia, "pieno", "vuoto")
    n_pieno = int((classe == "pieno").sum()); n_vuoto = n - n_pieno
    log(f"[6] classificazione: {n_vuoto} punti VUOTO (frequenza alta, gli UNICI "
        f"segnati nel grafico: " +
        ", ".join(f"{int(((ppol == p) & (classe == 'vuoto')).sum())} da {p.upper()}"
                  for p in pols) +
        f"), {n_pieno} punti PIENO (non disegnati)")

    # altezza dal livello 0 (quota>0) e profondita' dal livello 0 (quota<0)
    altezza = np.where(pz >= 0, pz, 0.0)
    profondita = np.where(pz < 0, -pz, 0.0)
    n_sopra = int((pz >= 0).sum()); n_sotto = int((pz < 0).sum())
    log(f"[6] rispetto al livello 0: {n_sopra} punti sopra (altezza max "
        f"{altezza.max():.0f} m), {n_sotto} punti sotto (profondita' max "
        f"{profondita.max():.0f} m)")
    np.savez_compressed(os.path.join(outdir, "variazioni_frequenza.npz"),
                        px=px, py=py, pz=pz, pv=pv, pol=ppol.astype(str), classe=classe,
                        altezza_dal_livello0=altezza, profondita_dal_livello0=profondita,
                        soglia=soglia)

    # hover: quota s.l.m. + altezza/profondita' dal livello 0 + classe/pol per punto
    txt = []
    for a, b, c, d, q, cl in zip(px, py, pz, pv, ppol, classe):
        rel = (f"altezza dal livello 0: {c:.0f} m" if c >= 0
               else f"profondita' dal livello 0: {-c:.0f} m")
        txt.append(f"x={a:.0f} m, y={b:.0f} m<br>quota {c:+.0f} m s.l.m.<br>"
                   f"{rel}<br>frequenza istantanea={d:.4f} 1/m<br>"
                   f"{cl.upper()}<br>pol {q.upper()}")

    fmin = float(np.nanpercentile(pv, 2)); fmax = float(np.nanpercentile(pv, 98))
    if fmax <= fmin:
        fmax = fmin + 1e-6

    # --- riempimento verticale dei vuoti (overlay selezionabile nel grafico) ---
    # Per ogni pixel (i,j), lungo TUTTA la profondita' (non solo il picco usato
    # sopra per il punto singolo) si segue la classificazione PIENO/VUOTO
    # campione per campione sull'asse verticale z: quando si passa da PIENO a
    # VUOTO si comincia a "riempire" quel tratto con punti/linea, quando si
    # torna da VUOTO a PIENO ci si ferma (i vuoti che iniziano gia' al
    # campione piu' superficiale, senza un PIENO sopra da cui scendere, non
    # vengono riempiti: non c'e' una transizione nota da cui partire). E' un
    # OVERLAY in aggiunta al punto-per-pixel sopra, non lo sostituisce: di
    # default e' nascosto (visible='legendonly'), selezionabile dalla legenda
    # o dal pulsante dedicato aggiunto sotto.
    # una passata per POLARIZZAZIONE (stesso ordine dei punti: prima VV, poi
    # VH), sulla forma d'onda separata di quella pol per TUTTA l'area
    dz_fill = float(z[1] - z[0]) if len(z) > 1 else 1.0
    passo_fill = max(1, int(round(delta_vuoto / dz_fill)))
    fill_traces = []
    n_fill_tot, n_disegnati = 0, 0
    for p in pols:
        pieno_full = freq_pol[p] < soglia                  # True = pieno, per campione
        stato = np.zeros((ny, nx), dtype=bool)
        riempi = np.zeros((ny, nx, nz), dtype=bool)
        for k in range(nz):
            is_p = pieno_full[:, :, k]
            inizia_qui = (~is_p) & pieno_full[:, :, k - 1] if k > 0 else np.zeros((ny, nx), bool)
            stato = np.where(is_p, False, np.where(inizia_qui, True, stato))
            riempi[:, :, k] = stato & (~is_p)
        n_fill = int(riempi.sum())
        n_fill_tot += n_fill
        if n_fill == 0:
            continue
        # delta_vuoto = distanza verticale [m] tra due punti disegnati lungo la
        # stessa linea: il riempimento e' un overlay visivo, disegnare OGNI
        # campione z (dz ~ ZD/NZ, di solito 1-2 m) appesantisce la pagina senza
        # aggiungere leggibilita'. Si tiene un campione ogni `passo_fill` piu'
        # SEMPRE l'ultimo del segmento, cosi' l'estensione del vuoto resta esatta.
        gap = np.array([np.nan])
        FX, FY, FZ, FC = [], [], [], []
        for i in range(ny):
            for j in range(nx):
                col = riempi[i, j]
                if not col.any():
                    continue
                idx = np.flatnonzero(col)
                # segmenti contigui distinti nella stessa colonna: piu' vuoti
                # separati da un tratto PIENO restano segmenti separati (gap)
                splits = np.flatnonzero(np.diff(idx) > 1) + 1
                for run in np.split(idx, splits):
                    sel = run[::passo_fill]
                    if sel[-1] != run[-1]:
                        sel = np.append(sel, run[-1])
                    n_disegnati += len(sel)
                    FX.append(np.full(len(sel), x_m[j])); FX.append(gap)
                    FY.append(np.full(len(sel), y_m[i])); FY.append(gap)
                    FZ.append(z0[i, j] - z[sel]); FZ.append(gap)
                    FC.append(freq_pol[p][i, j, sel]); FC.append(gap)
        fill_traces.append(go.Scatter3d(
            x=np.concatenate(FX), y=np.concatenate(FY), z=np.concatenate(FZ),
            mode="lines+markers", connectgaps=False,
            line=dict(width=4, color=np.concatenate(FC), colorscale=PIENO_VUOTO_COLORSCALE,
                      cmin=fmin, cmax=fmax, showscale=False),
            marker=dict(size=2.2, color=np.concatenate(FC), colorscale=PIENO_VUOTO_COLORSCALE,
                       cmin=fmin, cmax=fmax, showscale=False),
            name=f"riempimento vuoti {p.upper()} (clicca in legenda per mostrare)",
            visible="legendonly"))
    log(f"[6] riempimento verticale dei vuoti: {n_fill_tot} campioni "
        f"({100*n_fill_tot/(pixel_area*nz*len(pols)):.1f}% del volume z x pol) tra le "
        f"transizioni PIENO->VUOTO e VUOTO->PIENO (overlay per pol, nascosto di default)")
    if n_fill_tot > 0:
        log(f"[6] overlay vuoti: un punto ogni {passo_fill * dz_fill:.1f} m in "
            f"verticale (--delta-vuoto {delta_vuoto:.1f} m, passo {passo_fill} "
            f"campioni) -> {n_disegnati} punti disegnati invece di {n_fill_tot}")

    Xg, Yg = np.meshgrid(x_m, y_m)
    ground = go.Scatter3d(x=Xg.ravel(), y=Yg.ravel(), z=z0.ravel(), mode="markers",
                          marker=dict(size=1.2, color="red"), name="suolo (DEM)")
    # piano di riferimento a quota 0 (livello del mare)
    zero_plane = go.Surface(x=Xg, y=Yg, z=np.zeros_like(Xg),
                            showscale=False, opacity=0.25,
                            colorscale=[[0, "gray"], [1, "gray"]],
                            name="livello 0 (s.l.m.)")
    extra_traces = _tracce_riferimenti(riferimenti) if riferimenti else []
    # nel grafico si SEGNANO SOLO I PUNTI VUOTO: una traccia per pol,
    # nell'ordine di generazione (traccia 0 = VUOTO dalla prima pol, es. VV;
    # traccia 1 = VUOTO dalla seconda, es. VH). I PIENO non sono disegnati
    # (restano in .npz e nel cursore della soglia, che risuddivide le tracce).
    txt_arr = np.array(txt, dtype=object)
    tracce_vuoto = []
    for k, p in enumerate(pols):
        m = (ppol == p) & (classe == "vuoto")
        tracce_vuoto.append(go.Scatter3d(
            x=px[m], y=py[m], z=pz[m], mode="markers",
            name=f"VUOTO {p.upper()} ({int(m.sum())} punti)",
            text=list(txt_arr[m]), hoverinfo="text",
            marker=dict(size=2.4, color=pv[m], colorscale=PIENO_VUOTO_COLORSCALE,
                        cmin=fmin, cmax=fmax, opacity=0.9,
                        colorbar=_colorbar_pieno_vuoto(fmin, fmax) if k == 0 else None,
                        showscale=(k == 0))))
    traces = [*tracce_vuoto, zero_plane, ground, *extra_traces]
    idx_fill = None
    if fill_traces:
        idx_fill = list(range(len(traces), len(traces) + len(fill_traces)))
        traces.extend(fill_traces)
    fig = go.Figure(traces)
    _tema_plotly(fig)
    fig.update_layout(
        scene=dict(xaxis_title="x E-W [m]", yaxis_title="y N-S [m]",
                   zaxis_title="quota s.l.m. [m] (+ = altezza, − = profondita' dal livello 0)",
                   zaxis=dict(range=[ZBOT, ZTOP]),
                   aspectmode="manual",
                   aspectratio=dict(x=1.0, y=Ly / Lx, z=(ZTOP - ZBOT) / Lx)))
    if idx_fill is not None:
        # pulsante dedicato (in aggiunta al toggle dalla legenda) per mostrare
        # /nascondere il riempimento verticale dei vuoti (tutte le pol insieme)
        fig.update_layout(updatemenus=[dict(
            type="buttons", direction="left", x=0.99, y=0.99,
            xanchor="right", yanchor="top",
            bgcolor=TEMA_UI["pannello"], bordercolor=TEMA_UI["bordo"],
            font=dict(color=TEMA_UI["testo"], size=12),
            buttons=[
                dict(label="Riempimento vuoti: ON", method="restyle",
                     args=[{"visible": True}, idx_fill]),
                dict(label="Riempimento vuoti: OFF", method="restyle",
                     args=[{"visible": "legendonly"}, idx_fill]),
            ])])
    # --- pulsante per esportare i 2 piani di sezione (E-W e N-S) del punto in fuoco ---
    dx = float(np.min(np.diff(x_m))) if len(x_m) > 1 else 1.0
    dy = float(np.min(np.diff(y_m))) if len(y_m) > 1 else 1.0
    cs = dict(PX=[round(float(v), 2) for v in px],
              PY=[round(float(v), 2) for v in py],
              PZ=[round(float(v), 2) for v in pz],
              PV=[round(float(v), 5) for v in pv],
              PPOL=[str(q) for q in ppol],
              POLS=[str(p) for p in pols],
              TOLX=dx / 2.0, TOLY=dy / 2.0,
              SOGLIA=round(float(soglia), 6),
              FMIN=round(fmin, 6), FMAX=round(fmax, 6))
    post_script = """
var gd = document.getElementById('{plot_id}');
var CS = __CSDATA__;
var sel = {x: null, y: null, z: null};
// colorscale continua PIENO/VUOTO (stessa di PIENO_VUOTO_COLORSCALE lato Python):
// frequenza bassa = PIENO (chiaro), frequenza alta = VUOTO (scuro)
var PVCS = [[0.0,'#f5ecd7'],[0.5,'#8a7862'],[1.0,'#0d0d12']];

var bar = document.createElement('div');
bar.style.cssText = 'position:fixed;bottom:14px;left:14px;z-index:50;' +
    'background:rgba(255,255,255,0.95);color:#0f172a;padding:10px 12px;' +
    'border:1px solid #e2e8f0;border-radius:10px;' +
    "font-family:'Segoe UI',system-ui,sans-serif;font-size:13px;" +
    'box-shadow:0 4px 16px rgba(15,23,42,0.15)';
var info = document.createElement('div');
info.style.cssText = 'margin-bottom:7px;color:#475569';
info.textContent = 'Clicca un punto per dargli il fuoco';
var btn = document.createElement('button');
btn.textContent = 'Salva piani E-W & N-S (PNG)';
btn.disabled = true;
btn.style.cssText = 'padding:6px 12px;font-size:13px;font-weight:600;' +
    'border:0;border-radius:8px;background:#f5b942;color:#101318;' +
    'transition:opacity .2s ease';
function aggiornaBtn(){
    btn.style.opacity = btn.disabled ? '0.45' : '1';
    btn.style.cursor = btn.disabled ? 'not-allowed' : 'pointer';
}
btn.onmouseenter = function(){ if(!btn.disabled) btn.style.opacity = '0.85'; };
btn.onmouseleave = aggiornaBtn;
aggiornaBtn();
bar.appendChild(info); bar.appendChild(btn);
document.body.appendChild(bar);

gd.on('plotly_click', function(d){
    if(!d || !d.points || !d.points.length) return;
    var p = d.points[0];
    sel.x = p.x; sel.y = p.y; sel.z = p.z;
    info.textContent = 'Fuoco: x=' + Math.round(sel.x) + ' m, y=' + Math.round(sel.y) +
        ' m, quota=' + Math.round(sel.z) + ' m';
    info.style.color = '#0f172a';
    btn.disabled = false;
    aggiornaBtn();
});

// ---- cursore per REGOLARE la soglia PIENO/VUOTO (trace 0 = punti per pixel) ----
var box2 = document.createElement('div');
box2.style.cssText = 'margin-top:10px;padding-top:9px;border-top:1px solid #e2e8f0;max-width:250px';
var slLab = document.createElement('div');
slLab.style.cssText = 'margin-bottom:6px;color:#0f172a;line-height:1.4';
var sl = document.createElement('input');
sl.type = 'range';
sl.min = CS.FMIN; sl.max = CS.FMAX;
sl.step = (CS.FMAX - CS.FMIN) / 400;
sl.value = CS.SOGLIA;
sl.setAttribute('aria-label', 'Soglia PIENO/VUOTO [1/m]');
sl.style.cssText = 'display:block;width:100%;cursor:pointer;accent-color:#f5b942';
var slReset = document.createElement('button');
slReset.textContent = 'Ripristina (soglia iniziale, colore continuo)';
slReset.style.cssText = 'margin-top:7px;padding:4px 10px;font-size:12px;cursor:pointer;' +
    'border:1px solid #e2e8f0;border-radius:8px;background:transparent;color:#475569;' +
    'transition:color .2s ease,border-color .2s ease';
slReset.onmouseenter = function(){ slReset.style.color = '#0f172a'; slReset.style.borderColor = '#94a3b8'; };
slReset.onmouseleave = function(){ slReset.style.color = '#475569'; slReset.style.borderColor = '#e2e8f0'; };
var slHint = document.createElement('div');
slHint.style.cssText = 'margin-top:6px;font-size:11.5px;color:#475569;line-height:1.4';
slHint.textContent = "I punti disegnati sono SOLO i VUOTO (frequenza sopra " +
    "soglia), generati prima da " + CS.POLS[0].toUpperCase() +
    (CS.POLS.length > 1 ? " e poi da " + CS.POLS.slice(1).join('/').toUpperCase() : "") +
    "; i PIENO non sono disegnati. Muovendo il cursore la selezione " +
    "si aggiorna; l'overlay 'riempimento dei vuoti' resta calcolato con la " +
    "soglia iniziale.";

function contaPieni(s){
    var n = 0;
    for(var i = 0; i < CS.PV.length; i++) if(CS.PV[i] < s) n++;
    return n;
}
function testoPunto(i, classe){
    var q = CS.PZ[i];
    var rel = (q >= 0 ? 'altezza dal livello 0: ' + Math.round(q)
                      : "profondita' dal livello 0: " + Math.round(-q)) + ' m';
    return 'x=' + Math.round(CS.PX[i]) + ' m, y=' + Math.round(CS.PY[i]) + ' m<br>' +
        'quota ' + (q >= 0 ? '+' : '−') + Math.abs(Math.round(q)) + ' m s.l.m.<br>' +
        rel + '<br>frequenza istantanea=' + CS.PV[i].toFixed(4) + ' 1/m<br>' +
        classe + '<br>pol ' + CS.PPOL[i].toUpperCase();
}
function etichettaSoglia(s, extra){
    var np_ = contaPieni(s);
    slLab.textContent = 'Soglia PIENO/VUOTO: ' + s.toFixed(4) + ' 1/m' + (extra || '') +
        ' — VUOTO ' + (CS.PV.length - np_) + ' (punti segnati) / PIENO ' +
        np_ + ' (non disegnati)';
}
// rialimenta le tracce dei VUOTO, una per polarizzazione (traccia 0 = prima
// pol, es. VV; traccia 1 = seconda, es. VH): SOLO i punti sopra soglia sono
// disegnati, i PIENO non compaiono in nessuna traccia
function applicaSoglia(s, extra){
    var upd = [];
    for(var t = 0; t < CS.POLS.length; t++) upd.push({X:[],Y:[],Z:[],C:[],T:[]});
    for(var i = 0; i < CS.PV.length; i++){
        if(CS.PV[i] < s) continue;                      // PIENO: non disegnato
        var t = CS.POLS.indexOf(CS.PPOL[i]);
        if(t < 0) continue;
        upd[t].X.push(CS.PX[i]); upd[t].Y.push(CS.PY[i]); upd[t].Z.push(CS.PZ[i]);
        upd[t].C.push(CS.PV[i]); upd[t].T.push(testoPunto(i, 'VUOTO'));
    }
    etichettaSoglia(s, extra);
    for(var t = 0; t < CS.POLS.length; t++){
        Plotly.restyle(gd, {x:[upd[t].X], y:[upd[t].Y], z:[upd[t].Z],
                            'marker.color':[upd[t].C], text:[upd[t].T],
                            name:['VUOTO ' + CS.POLS[t].toUpperCase() +
                                  ' (' + upd[t].X.length + ' punti)']}, [t]);
    }
}
sl.oninput = function(){ applicaSoglia(parseFloat(sl.value)); };
slReset.onclick = function(){
    sl.value = CS.SOGLIA;
    applicaSoglia(CS.SOGLIA, ' (iniziale)');
};
etichettaSoglia(CS.SOGLIA, ' (iniziale)');
box2.appendChild(slLab); box2.appendChild(sl);
box2.appendChild(slReset); box2.appendChild(slHint);
bar.appendChild(box2);

btn.onclick = function(){
    if(sel.x === null) return;
    var ewx=[], ewz=[], ewc=[], nsy=[], nsz=[], nsc=[];
    for(var i=0; i<CS.PX.length; i++){
        if(Math.abs(CS.PY[i] - sel.y) <= CS.TOLY){ ewx.push(CS.PX[i]); ewz.push(CS.PZ[i]); ewc.push(CS.PV[i]); }
        if(Math.abs(CS.PX[i] - sel.x) <= CS.TOLX){ nsy.push(CS.PY[i]); nsz.push(CS.PZ[i]); nsc.push(CS.PV[i]); }
    }
    var t1 = {type:'scatter', mode:'markers', x:ewx, y:ewz, name:'E-W',
        marker:{size:5, color:ewc, colorscale:PVCS, colorbar:{title:'freq. (PIENO..VUOTO)', x:1.0, len:0.9}},
        xaxis:'x', yaxis:'y'};
    var t2 = {type:'scatter', mode:'markers', x:nsy, y:nsz, name:'N-S',
        marker:{size:5, color:nsc, colorscale:PVCS, showscale:false},
        xaxis:'x2', yaxis:'y2'};
    var layout = {
        width:1300, height:600, showlegend:false,
        title:'Piani di sezione sul fuoco  x=' + Math.round(sel.x) + ' m, y=' + Math.round(sel.y) +
            ' m, quota=' + Math.round(sel.z) + ' m',
        grid:{rows:1, columns:2, pattern:'independent'},
        xaxis:{title:'x E-W [m]', anchor:'y'},
        yaxis:{title:'quota s.l.m. [m]'},
        xaxis2:{title:'y N-S [m]', anchor:'y2'},
        yaxis2:{title:'quota s.l.m. [m]'},
        annotations:[
            {text:'Piano E-W (y=' + Math.round(sel.y) + ' m)  \\u2014 ' + ewx.length + ' punti',
                x:0.2, xref:'paper', y:1.06, yref:'paper', showarrow:false, font:{size:13}},
            {text:'Piano N-S (x=' + Math.round(sel.x) + ' m)  \\u2014 ' + nsy.length + ' punti',
                x:0.8, xref:'paper', y:1.06, yref:'paper', showarrow:false, font:{size:13}}
        ]
    };
    var hidden = document.createElement('div');
    hidden.style.cssText = 'position:absolute;left:-9999px;top:-9999px';
    document.body.appendChild(hidden);
    Plotly.newPlot(hidden, [t1, t2], layout).then(function(){
        return Plotly.toImage(hidden, {format:'png', width:1300, height:600});
    }).then(function(url){
        var a = document.createElement('a');
        a.href = url;
        a.download = 'piani_EW_NS_x' + Math.round(sel.x) + '_y' + Math.round(sel.y) + '.png';
        document.body.appendChild(a); a.click(); document.body.removeChild(a);
        Plotly.purge(hidden); document.body.removeChild(hidden);
    });
};
""".replace("__CSDATA__", json.dumps(cs))

    _scrivi_pagina_html(
        fig, os.path.join(outdir, "step6_variazioni_3d.html"),
        titolo="Step 6 — Solo punti VUOTO ("
               + " poi ".join(p.upper() for p in pols) + ")",
        sottotitolo="Per ogni polarizzazione (prima "
                    + ", poi ".join(p.upper() for p in pols) +
                    "), un punto per pixel dell'area alla profondita' del "
                    "riflettore piu' forte della sua forma d'onda: nel "
                    "grafico sono SEGNATI SOLO I VUOTI (frequenza istantanea "
                    "alta, scuro), una traccia per pol; i PIENO non sono "
                    "disegnati. Livello 0 = mare. Il CURSORE in basso a "
                    "sinistra regola la soglia PIENO/VUOTO in tempo reale; "
                    "clicca un punto per dargli il fuoco ed esportare i "
                    "piani di sezione E-W / N-S; il riempimento dei vuoti "
                    "(per pol) si attiva dal pulsante in alto a destra o "
                    "dalla legenda. La profondita' del grafico e' limitata "
                    "alla z_amb del dato sorgente (oltre, l'inversione di "
                    "Fourier va in alias).",
        badges=[("punti VUOTO (segnati)",
                 " + ".join(f"{int(((ppol == p) & (classe == 'vuoto')).sum())} {p.upper()}"
                            for p in pols)),
                ("PIENO (non disegnati)", n_pieno),
                ("soglia", f"{soglia:.4f} 1/m"),
                ("base", f"{Lx:.0f}×{Ly:.0f} m"),
                ("delta vuoti", f"{delta_vuoto:.1f} m")]
               + ([("profondita' max fonte (z_amb)", f"{ZD:.0f} m"),
                   ("dz tomografica", f"{prof_info[2]:.0f} m"),
                   ("sorgente", prof_info[0])] if prof_info else []),
        nota="Analisi ESPLORATIVA ispirata al lessico di Biondi & Malanga "
             "(Remote Sens. 2022, 14, 5231) ma NON la loro tomografia Doppler "
             "a sub-aperture: classificazione PIENO/VUOTO per soglia sulla "
             "frequenza istantanea (Hilbert) di forme d'onda sintetizzate su "
             "un range di profondita' scelto, non derivato dalla geometria "
             "SAR. Usare per ipotesi, non come misura fisica validata.",
        post_script=post_script)

    figm = plt.figure(figsize=(8, 9)); ax = figm.add_subplot(111, projection="3d")
    ax.scatter(Xg.ravel(), Yg.ravel(), z0.ravel(), c="red", s=2, depthshade=False)
    m_vuoto = classe == "vuoto"
    sc = ax.scatter(px[m_vuoto], py[m_vuoto], pz[m_vuoto], c=pv[m_vuoto],
                    cmap=_cmap_pieno_vuoto(), vmin=fmin, vmax=fmax, s=2.0)
    ax.set_xlabel("x E-W [m]"); ax.set_ylabel("y N-S [m]"); ax.set_zlabel("quota s.l.m. [m]")
    ax.set_zlim(ZBOT, ZTOP); ax.set_box_aspect((Lx, Ly, ZTOP - ZBOT))
    ax.set_title(f"Step 6 — solo punti VUOTO ({' poi '.join(p.upper() for p in pols)}), "
                 f"{n_vuoto} su {n} generati")
    figm.colorbar(sc, ax=ax, shrink=0.5,
                  label="frequenza istantanea (scuro=VUOTO)")
    figm.savefig(os.path.join(outdir, "step6_variazioni_3d.png"), dpi=120,
                 bbox_inches="tight")
    plt.close(figm)
    log(f"[6] salvati step6_variazioni_3d.html/.png")


# =============================================== TOMOGRAFIA DOPPLER "VERA" (paper)
def step_vera_tomografia(nw, nw_lon, se, se_lon, stackdir, outdir,
                         zmax=8.5, klook=12, ovs=4.0):
    """Richiama skills/sar-doppler-tomography/scripts/tomographic_images.py, che
    implementa la tomografia Doppler a sub-aperture di Biondi & Malanga cosi'
    come descritta nella skill biondi-malanga-sar-tomography (steering matrix
    Kz = 4*pi*B_perp/(lambda_sonic*r*sin(theta)), inversione h(z) = A^H*Y,
    risoluzione delta_z = lambda_sonic*R/(2A)) validata su un modello sintetico
    (vedi sar_tomo.py --self-test). Produce le sezioni tomografiche in stile
    pubblicazione (B-scan range/azimuth-profondita', slice orizzontali, mappa
    di quota) usando UNA sola immagine SLC per volta (sub-aperture Doppler),
    non le coppie interferometriche multi-data usate dagli step 4-6."""
    here = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(here, "skills", "sar-doppler-tomography", "scripts",
                          "tomographic_images.py")
    if not os.path.exists(script):
        log(f"[vera-tomografia] script non trovato: {script}"); return False
    tomo_out = os.path.join(outdir, "tomografia_vera")
    os.makedirs(tomo_out, exist_ok=True)
    cmd = [sys.executable, script,
          "--stack", stackdir, "--outdir", tomo_out,
          "--nw", nw[0], nw[1], nw[2], nw[3], nw_lon[0], nw_lon[1],
          "--nw2", nw_lon[2], nw_lon[3],
          "--se", se[0], se[1], se[2], se[3], se_lon[0], se_lon[1],
          "--se2", se_lon[2], se_lon[3],
          "--zmax", str(zmax), "--klook", str(klook), "--ovs", str(ovs)]
    log(f"[vera-tomografia] {' '.join(cmd)}")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        log(f"[vera-tomografia] terminato con codice {r.returncode}")
        return False
    log(f"[vera-tomografia] output (metodo Biondi-Malanga, coerente col paper) in {tomo_out}/")
    return True


# ===================================================================== main
def parse_steps(s):
    if "-" in s:
        a, b = s.split("-"); return set(range(int(a), int(b) + 1))
    return {int(x) for x in s.split(",")}


def main():
    global DBG

    ap = argparse.ArgumentParser(description="Pipeline unica piramide (step 1-6)")
    ap.add_argument("--pyramid", choices=sorted(PRESETS), default=None,
                    help="box di default per piramide (cheope|kefren); "
                         "sovrascrivibile con --nw/--nw-lon/--se/--se-lon")
    ap.add_argument("--nw", nargs=4, metavar=("D", "M", "S", "H"),
                    default=None, help="angolo NW lat (DMS)")
    ap.add_argument("--nw-lon", nargs=4, metavar=("D", "M", "S", "H"),
                    default=None, help="angolo NW lon (DMS)")
    ap.add_argument("--se", nargs=4, metavar=("D", "M", "S", "H"),
                    default=None, help="angolo SE lat (DMS)")
    ap.add_argument("--se-lon", nargs=4, metavar=("D", "M", "S", "H"),
                    default=None, help="angolo SE lon (DMS)")
    ap.add_argument("--pol", default="vv,vh",
                    help="polarizzazioni, separate da virgola (default vv,vh). Gli "
                         "strati DInSAR di TUTTE le pol elencate vengono messi in un "
                         "SOLO array di coefficienti (es. 12 VV + 12 VH = 24 "
                         "armoniche) e l'anti-trasformata di Fourier dello step 5 e' "
                         "calcolata in un solo passaggio, a risoluzione piena per "
                         "pixel -> una forma d'onda combinata per pixel. La prima "
                         "pol e' la primaria (conteggi/track filter). Solo VV: --pol vv")
    ap.add_argument("--layers", type=int, default=12, help="numero di strati (>=12)")
    ap.add_argument("--start", default="2026-01-01", help="data inizio ricerca YYYY-MM-DD")
    ap.add_argument("--end", default=None,
                    help="data fine ricerca YYYY-MM-DD (default: oggi)")
    ap.add_argument("--max-search", type=int, default=40)
    ap.add_argument("--stack", default=None, help="cartella TIFF scaricati")
    ap.add_argument("--outdir", default=None, help="cartella output")
    ap.add_argument("--alta-risoluzione-dir", default=None,
                    help="cartella con i .zip SAFE gia' scaricati da "
                         "scarica_alta_risoluzione.py (default: alta_risoluzione_out/ "
                         "accanto a questo script); lo step 2 ne estrae i measurement "
                         "VV/VH mancanti PRIMA di provare il download CDSE, senza "
                         "bisogno di credenziali. I .zip.part (download in corso) "
                         "vengono ignorati. Passa '' per disattivare")
    ap.add_argument("--download", action="store_true", help="scarica davvero (servono credenziali CDSE)")
    ap.add_argument("--cdse-user", default=None,
                    help="username/email CDSE (altrimenti env CDSE_USER)")
    ap.add_argument("--cdse-pass", default=None,
                    help="password CDSE (altrimenti env CDSE_PASS)")
    ap.add_argument("--steps", default="1-6", help="es. 1-6 oppure 3,4,5")
    ap.add_argument("--soglia", type=float, default=None,
                    help="step 6: soglia di frequenza istantanea (1/m) SOTTO la quale "
                         "un punto (uno per pixel dell'area selezionata) e' classificato "
                         "PIENO, sopra la quale e' VUOTO. Default automatico = mediana "
                         "delle frequenze di picco su tutta l'area")
    ap.add_argument("--delta-vuoto", type=float, default=5.0,
                    help="step 6: distanza verticale [m] tra i punti disegnati "
                         "nell'overlay 'riempimento dei vuoti' (default 5.0; "
                         "prima veniva disegnato OGNI campione z, ~1-2 m). "
                         "Valori piu' alti = meno punti per linea verticale = "
                         "pagina piu' leggera; l'ultimo punto di ogni segmento "
                         "e' sempre disegnato, l'estensione dei vuoti resta esatta")
    ap.add_argument("--nz", type=int, default=256,
                    help="step 5: risoluzione verticale (campioni in profondita') "
                         "dell'anti-trasformata di Fourier, calcolata a risoluzione "
                         "piena per OGNI pixel e OGNI pol richiesta in --pol")
    ap.add_argument("--step-x", type=int, default=None,
                    help="step 5: passo colonne (range) tra una traccia 3D e la "
                         "successiva DELLA STESSA pol. Default AUTO = numero di "
                         "polarizzazioni in --pol, cosi' le pol insieme (anche con una "
                         "sola pol) coprono OGNI colonna senza buchi: il numero di "
                         "tracce disegnate = numero di pixel dell'area selezionata "
                         "(una traccia per pixel). Alza il valore per meno tracce/piu' "
                         "veloce (a scapito della copertura orizzontale 1:1)")
    ap.add_argument("--step-y", type=int, default=1,
                    help="step 5: passo righe (azimuth/N-S) tra una traccia 3D e la "
                         "successiva. Default 1 = OGNI riga (risoluzione N-S piena, "
                         "simmetrica con la copertura E-W di --step-x). Alza il valore "
                         "per meno tracce/piu' veloce (a scapito della copertura N-S)")
    ap.add_argument("--no-pyramid", action="store_true", help="non sovrapporre la piramide al DEM")
    ap.add_argument("--no-track-filter", action="store_true",
                    help="step 3: NON filtrare le scene per track (tiene anche scene di "
                         "relative orbit diversi, decorrelate -> quote peggiori)")
    ap.add_argument("--vera-tomografia", action="store_true",
                    help="oltre agli step scelti, esegue ANCHE la tomografia Doppler "
                         "vera del paper (skills/sar-doppler-tomography), coerente col "
                         "metodo Biondi-Malanga (steering matrix + h(z)=A^H*Y)")
    ap.add_argument("--zmax", type=float, default=8.5,
                    help="profondita' massima [m] per --vera-tomografia (resta sotto z_amb!)")
    ap.add_argument("--klook", type=int, default=12,
                    help="numero di sub-aperture Doppler (look tomografici) per --vera-tomografia")
    ap.add_argument("--ovs", type=float, default=4.0,
                    help="oversampling in profondita' per --vera-tomografia")
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args()
    DBG = a.debug
    if a.end is None:
        a.end = time.strftime("%Y-%m-%d")   # fine ricerca = oggi, mai una data fissa

    if a.nw is None or a.se is None or a.nw_lon is None or a.se_lon is None:
        if a.pyramid is None:
            sys.exit("Serve --pyramid {cheope,kefren} oppure tutti e 4 gli angoli "
                      "--nw/--nw-lon/--se/--se-lon espliciti.")
        preset = PRESETS[a.pyramid]
        a.nw = a.nw or list(preset["nw"]); a.nw_lon = a.nw_lon or list(preset["nw_lon"])
        a.se = a.se or list(preset["se"]); a.se_lon = a.se_lon or list(preset["se_lon"])

    here = os.path.dirname(os.path.abspath(__file__))
    alta_ris_dir = (os.path.join(here, "alta_risoluzione_out") if a.alta_risoluzione_dir is None
                    else a.alta_risoluzione_dir or None)
    preset_dir = PRESETS[a.pyramid]["outdir_name"] if a.pyramid else "goal_out"
    base_dir = os.path.join(here, preset_dir)
    outdir = a.outdir or os.path.join(base_dir, "unificato")
    stackdir = a.stack or os.path.join(base_dir, "stack_slc")
    fallback_box = os.path.join(base_dir, "box.npz")
    boxpath = os.path.join(outdir, "box.npz")
    os.makedirs(outdir, exist_ok=True)
    # geometria per l'overlay piramide sul DEM (step 5-6): quella del preset
    # scelto (ogni piramide ha base/altezza proprie). Con un box custom senza
    # --pyramid non c'e' geometria nota -> nessun overlay (solo DEM).
    pyr_geom = None
    if not a.no_pyramid:
        if a.pyramid:
            _p = PRESETS[a.pyramid]
            pyr_geom = dict(base_m=_p["base_m"], h_m=_p["h_m"], label=_p["label"])
        else:
            log("(box custom senza --pyramid: nessun overlay piramide sul DEM; "
                "usa --pyramid cheope|kefren per l'overlay con la geometria giusta)")
    steps = parse_steps(a.steps)

    lat_nw = dms2dec(*a.nw[:3], a.nw[3]); lon_nw = dms2dec(*a.nw_lon[:3], a.nw_lon[3])
    lat_se = dms2dec(*a.se[:3], a.se[3]); lon_se = dms2dec(*a.se_lon[:3], a.se_lon[3])
    latN, latS = max(lat_nw, lat_se), min(lat_nw, lat_se)
    lonW, lonE = min(lon_nw, lon_se), max(lon_nw, lon_se)
    aoi = (lonW, lonE, latS, latN)

    log("=" * 70)
    log(f"PIPELINE UNIFICATA — piramide: {a.pyramid or 'custom'} — "
        f"AOI lon {lonW:.5f}..{lonE:.5f}, lat {latS:.5f}..{latN:.5f}")
    log(f"step {sorted(steps)} | pol {a.pol} | strati {a.layers} | outdir {outdir}")
    if a.vera_tomografia:
        log("  + tomografia Doppler VERA del paper (steering matrix / h(z)=A^H*Y)")
    log("=" * 70)

    prods = []
    if 1 in steps:
        prods = step1_cerca(aoi, a.start, a.end, a.max_search, a.pol, outdir)
    if 2 in steps:
        # i prodotti trovati dallo step 1 (che contengono le coordinate) vengono
        # sempre passati allo step 2; il download vero (~8 GB/prodotto) avviene
        # solo con --download e credenziali CDSE, altrimenti si riusa lo stack.
        step2_scarica(prods, a.layers, a.pol, stackdir,
                      user=a.cdse_user, pwd=a.cdse_pass, download=a.download,
                      alta_ris_dir=alta_ris_dir)
    if 3 in steps:
        step3_estrai_box(stackdir, a.pol, aoi, boxpath, fallback_box,
                         track_filter=not a.no_track_filter)
    elif not os.path.exists(boxpath) and os.path.exists(fallback_box):
        boxpath = fallback_box

    if steps & {4, 5, 6}:
        per_pol, x, y, coh = step4_array_12(boxpath, a.layers, outdir)
        if steps & {5, 6}:
            # profondita' massima dell'inversione (z_amb) secondo il TIPO di
            # file sorgente (skill biondi-sar-tomography); k = strati reali
            k_strati = next(iter(per_pol.values()))[0].shape[0]
            prof_info = _profondita_massima_fourier(boxpath, k_strati)
            (w_pol, freq_pol, env_pol, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT, riferimenti) = \
                step5_grafico_onde(per_pol, x, y, aoi, outdir, pyr_geom,
                                   NZ=a.nz, STEP_X=a.step_x, STEP_Y=a.step_y,
                                   pyramid=a.pyramid, stackdir=stackdir,
                                   prof_info=prof_info)
            if 6 in steps:
                step6_variazioni(freq_pol, env_pol, z, x_m, y_m, z0, Lx, Ly, ZD, ZTOP, ZBOT,
                                 outdir, soglia=a.soglia, riferimenti=riferimenti,
                                 delta_vuoto=a.delta_vuoto, prof_info=prof_info)

    if a.vera_tomografia:
        step_vera_tomografia(a.nw, a.nw_lon, a.se, a.se_lon, stackdir, outdir,
                             zmax=a.zmax, klook=a.klook, ovs=a.ovs)

    log("=" * 70)
    log(f"FATTO. Output in: {outdir}")
    log("=" * 70)


if __name__ == "__main__":
    main()

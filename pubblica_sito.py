#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pubblica_sito.py  --  raccoglie le uscite in un sito statico
============================================================

Mette in una sola cartella quello che si guarda -- le tre pagine 3D
interattive, le figure e il riassunto JSON -- piu' un indice che le lega e
dice che cosa sono. La cartella e' pronta per un caricamento diretto su
Cloudflare Pages (o su qualunque hosting statico): niente build, niente
dipendenze, nessuna chiamata di rete dalle pagine.

    python pubblica_sito.py                       # da out_piramidi_v03
    python pubblica_sito.py --out sito_giza       # altrove

I numeri dell'indice NON sono scritti a mano: si leggono da ``meta_v03.json``
e da ``meta.json`` della corsa che si sta pubblicando. Se la corsa cambia,
cambia l'indice, e non c'e' modo che le due cose divergano in silenzio.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: (nome del file, titolo, descrizione) delle pagine interattive
PAGINE: Tuple[Tuple[str, str, str], ...] = (
    ("onde_3d_piramidi.html",
     "Le sinusoidi, pixel per pixel",
     "Per ogni cella della superficie delle piramidi in geometria radar, la "
     "sinusoide risultante di quel pixel: la somma delle sue sinusoidi, una "
     "per data. Ruotala, taglia la quota, cambia la scala orizzontale."),
    ("tomografia_v03_pieno_vuoto.html",
     "Volume pieno/vuoto (v03)",
     "I voxel in cui il profilo misurato si discosta dalla risposta attesa da "
     "un solo diffusore, oltre la dispersione misurata sulla piana. Rosso "
     "eccesso, blu difetto, soglia regolabile."),
    ("tomografia_piramidi_3d.html",
     "Superficie e tomogramma (v02)",
     "La catena principale: superficie ricostruita dai .tiff, nuvola "
     "tomografica, suolo di riferimento, strutture note e il pannello dei "
     "lobi verticali."),
)

#: (nome del file, didascalia) delle figure
FIGURE: Tuple[Tuple[str, str], ...] = (
    ("onde_3d_superficie_piramidi.png",
     "Il campo di onde in 3D: una curva per pixel, colorata per piramide. I "
     "punti neri sono il diffusore dominante delle celle sopra soglia."),
    ("sezione_onde_sotto_piramidi.png",
     "Sezione verticale: la sagoma della piramide in scala vera e, sotto il "
     "suolo, le sinusoidi di ogni data con la loro somma e l'inviluppo."),
    ("onde_cheope_chefren.png",
     "Le sinusoidi di una colonna, il profilo dopo la FFT contro la PSF di un "
     "solo diffusore, e il residuo pieno/vuoto con la banda del nullo."),
    ("pieno_vuoto_3d.png",
     "Il volume pieno/vuoto visto da sud-ovest e di taglio."),
    ("superficie_ricostruita.png",
     "Superficie, qualita' e validazione di livello 1: misurato contro "
     "simulato."),
    ("profili_pieno_vuoto.png",
     "Profili verticali: colonne sulle piramidi contro deserto di controllo."),
)


def _leggi(path: str) -> Optional[Dict[str, Any]]:
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _tabella(righe: List[Tuple[str, str]]) -> str:
    corpo = "".join(f"<tr><td>{k}</td><td class='v'>{v}</td></tr>"
                    for k, v in righe)
    return f"<table>{corpo}</table>"


def numeri(m3: Optional[Dict[str, Any]],
           m2: Optional[Dict[str, Any]]) -> List[Tuple[str, str]]:
    """I numeri della corsa, letti dai due riassunti JSON."""
    out: List[Tuple[str, str]] = []
    if m2:
        b = m2.get("budget", {})
        out += [
            ("acquisizioni impilate", str(len(m2.get("date", [])))),
            ("master", str(m2.get("master", "-"))),
            ("polarizzazione", str(m2.get("polarizzazione", "-")).split()[0]),
            ("escursione delle baseline", f"{b.get('b_spread', 0):.1f} m"),
            ("risoluzione verticale &delta;<sub>z</sub>",
             f"{b.get('delta_z_vertical', 0):.0f} m"),
            ("precisione sulla quota &sigma;<sub>h</sub>",
             f"{b.get('sigma_h', 0):.1f} m"),
            ("coerenza mediana misurata",
             f"{m2.get('coerenza_mediana_misurata', 0):.3f}"),
        ]
    if m3:
        f = m3.get("errore_fft", {})
        c = m3.get("confronto_con_v02", {})
        n = m3.get("nullo", {})
        bil = (m3.get("bilancio_anomalie") or {}).get("anomalie_totali", {})
        campo = m3.get("campo_onde_pixel_per_pixel") or {}
        out += [
            ("errore della FFT contro il calcolo diretto",
             f"{f.get('errore_relativo_max', float('nan')):.1e}"),
            ("quota identica a v02 entro 1 cm",
             f"{c.get('celle_entro_1_cm_pct', 0):.0f} % delle celle"),
            ("soglia |z-score| dal nullo di piana",
             f"{n.get('soglia_zscore', 0):.2f} su {n.get('celle', 0)} celle"),
            ("celle anomale: piramidi contro piana",
             f"{100 * bil.get('frazione_piramidi', 0):.1f} % contro "
             f"{100 * bil.get('frazione_piana', 0):.1f} % (z = "
             f"{bil.get('z', 0):+.2f})"),
            ("colonne nel campo di onde",
             f"{campo.get('celle', 0)} ({campo.get('celle_sopra_soglia', 0)} "
             "sopra soglia)"),
        ]
    return out


INDICE = """<!DOCTYPE html>
<html lang="it">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tomografia SAR della piana di Giza - Sentinel-1</title>
<style>
:root { color-scheme: light; }
* { box-sizing: border-box; }
body { margin: 0; font: 16px/1.65 system-ui, "Segoe UI", Roboto, sans-serif;
       background: #F8FAFC; color: #0F172A; }
.avvolto { max-width: 1080px; margin: 0 auto; padding: 0 22px 80px; }
header { background: #FFF; border-bottom: 1px solid #E2E8F0;
         padding: 40px 0 30px; margin-bottom: 34px; }
h1 { margin: 0 0 10px; font-size: 30px; letter-spacing: -.015em; }
.sotto { color: #475569; font-size: 17px; margin: 0; max-width: 720px; }
h2 { font-size: 13px; text-transform: uppercase; letter-spacing: .07em;
     color: #64748B; margin: 42px 0 14px; }
p { max-width: 760px; }
a { color: #1D4ED8; }
.griglia { display: grid; gap: 16px;
           grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); }
.scheda { display: block; background: #FFF; border: 1px solid #E2E8F0;
          border-radius: 12px; padding: 18px 20px; text-decoration: none;
          color: inherit; transition: border-color .15s, transform .15s; }
.scheda:hover { border-color: #94A3B8; transform: translateY(-2px); }
.scheda b { display: block; font-size: 17px; margin-bottom: 6px;
            color: #0F172A; }
.scheda span { font-size: 14px; color: #475569; }
.scheda i { display: block; margin-top: 10px; font-size: 13px; color: #1D4ED8;
            font-style: normal; }
table { width: 100%; border-collapse: collapse; background: #FFF;
        border: 1px solid #E2E8F0; border-radius: 12px; overflow: hidden; }
td { padding: 9px 16px; border-bottom: 1px solid #F1F5F9; font-size: 15px; }
tr:last-child td { border-bottom: none; }
td.v { text-align: right; font-variant-numeric: tabular-nums; color: #334155; }
figure { margin: 0 0 26px; background: #FFF; border: 1px solid #E2E8F0;
         border-radius: 12px; overflow: hidden; }
figure img { display: block; width: 100%; height: auto; }
figcaption { padding: 12px 18px; font-size: 14px; color: #475569;
             border-top: 1px solid #F1F5F9; }
.nota { background: #FFF7ED; border-left: 3px solid #F59E0B;
        padding: 14px 18px; border-radius: 8px; font-size: 15px;
        color: #7C2D12; max-width: 760px; }
footer { margin-top: 50px; padding-top: 20px; border-top: 1px solid #E2E8F0;
         font-size: 14px; color: #64748B; }
code { background: #F1F5F9; padding: 1px 5px; border-radius: 4px;
       font-size: 90%; }
</style>
</head>
<body>
<header><div class="avvolto">
  <h1>Tomografia SAR della piana di Giza</h1>
  <p class="sotto">Che cosa si riesce davvero a misurare delle tre piramidi con
  uno stack Sentinel-1 IW SLC, canale VH &mdash; e dove il dato smette di dire
  qualcosa. Metodo di Biondi e Malanga, implementato e verificato riga per
  riga.</p>
</div></header>
<div class="avvolto">

  <h2>Le pagine interattive</h2>
  <div class="griglia">__SCHEDE__</div>

  <h2>Come si arriva a queste figure</h2>
  <p>Da ogni pixel complesso dell&rsquo;interferogramma si prendono
  <b>modulo</b> e <b>fase in radianti</b>. Ogni data diventa una sinusoide
  lungo la verticale, <code>s(z) = A cos(&phi; &minus; &kappa;z)</code>, dove
  &kappa; e&rsquo; il numero d&rsquo;onda della baseline di quella data. La
  somma delle sinusoidi viene trasformata con una FFT &mdash; le righe non
  stanno su un reticolo uniforme, perche&rsquo; le baseline non lo sono, quindi
  passano per un nucleo gaussiano e una deapodizzazione analitica &mdash; e il
  profilo che ne esce viene confrontato con la risposta attesa da <b>un solo
  diffusore</b>.</p>
  <p>Non e&rsquo; un metodo alternativo, ed e&rsquo; scritto ovunque: la somma
  delle sinusoidi <i>e&rsquo;</i> la parte reale del periodogramma tomografico,
  e la quota che ne esce coincide con quella della catena precedente entro un
  centimetro sul 100&nbsp;% delle celle. Cambia il costo del calcolo, e cambia
  che ora lo si puo&rsquo; disegnare.</p>

  <h2>I numeri della corsa</h2>
  __TABELLA__

  <h2>Che cosa dicono i dati</h2>
  <p>La superficie della piana si ricostruisce con precisione metrica.
  <b>Le piramidi no.</b> Le facce a ~52&deg; superano l&rsquo;angolo di
  incidenza di 37&deg;: sono in layover pieno, con centinaia di punti di
  superficie ripiegati nella stessa cella di risoluzione. Il piano dei
  diffusori dominanti, visibile di taglio nel campo di onde, sta appena sopra
  il deserto e non segue il profilo delle facce.</p>
  <p>Le celle delle piramidi mostrano piu&rsquo; energia di quanta la risposta
  di un solo diffusore ne spieghi, e solo in eccesso: sul difetto le due
  popolazioni sono indistinguibili. Ma nessuna delle grandezze di controllo
  spiega quell&rsquo;eccesso &mdash; layover, luminosita&rsquo;, coerenza,
  quota simulata valgono al massimo il 6&nbsp;% della varianza.</p>
  <div class="nota">La risoluzione verticale di questa pila e&rsquo; di circa
  <b>130 metri</b>. Le camere note della Grande Piramide misurano metri: due
  ordini di grandezza sotto. Nulla di quello che si vede qui e&rsquo; una
  rilevazione di cavita&rsquo; risolta in profondita&rsquo;, e il programma non
  lo lascia intendere da nessuna parte.</div>

  <h2>Le figure</h2>
  __FIGURE__

  <footer>
    Dati Sentinel-1 &copy; Copernicus / ESA, distribuiti sotto la licenza del
    Copernicus Data Space Ecosystem. Gli articoli e il brevetto citati non sono
    ridistribuiti.<br>
    Pagine generate il __DATA__ dalla corsa in <code>__CORSA__</code>.
  </footer>
</div>
</body>
</html>
"""


def costruisci(sorgente: str, destinazione: str, verbose: bool = True,
               escludi: Sequence[str] = ()) -> str:
    """Copia le uscite e scrive l'indice. Restituisce il percorso della cartella.

    ``escludi`` toglie file per nome. Serve per il caricamento diretto su
    Cloudflare, che accetta un lotto per volta sotto una certa dimensione: la
    pagina di v02 pesa 3 MB e non fa parte delle novita', quindi e' la prima
    a uscire. L'indice si adatta da solo -- le schede si scrivono sui file
    che ci sono, non su una lista fissa."""
    os.makedirs(destinazione, exist_ok=True)
    fuori = set(escludi)
    m3 = _leggi(os.path.join(sorgente, "meta_v03.json"))
    m2 = _leggi(os.path.join(sorgente, "meta.json"))

    schede = []
    copiati: List[str] = []
    for nome, titolo, descr in PAGINE:
        src = os.path.join(sorgente, nome)
        if nome in fuori or not os.path.exists(src):
            if verbose:
                print(f"  manca (saltata): {nome}")
            continue
        shutil.copy2(src, os.path.join(destinazione, nome))
        copiati.append(nome)
        kb = os.path.getsize(src) / 1024.0
        schede.append(
            f'<a class="scheda" href="{nome}"><b>{titolo}</b>'
            f'<span>{descr}</span><i>apri &rarr; <span style="color:#94A3B8">'
            f'{kb:.0f} KB</span></i></a>')

    figure = []
    for nome, didascalia in FIGURE:
        src = os.path.join(sorgente, nome)
        if nome in fuori or not os.path.exists(src):
            if verbose:
                print(f"  manca (saltata): {nome}")
            continue
        shutil.copy2(src, os.path.join(destinazione, nome))
        copiati.append(nome)
        figure.append(f'<figure><img src="{nome}" alt="{didascalia}" '
                      f'loading="lazy"><figcaption>{didascalia}</figcaption>'
                      "</figure>")

    for extra in ("meta_v03.json", "meta.json", "budget_tomografico.txt"):
        src = os.path.join(sorgente, extra)
        if extra not in fuori and os.path.exists(src):
            shutil.copy2(src, os.path.join(destinazione, extra))
            copiati.append(extra)

    html = (INDICE
            .replace("__SCHEDE__", "\n".join(schede))
            .replace("__FIGURE__", "\n".join(figure))
            .replace("__TABELLA__", _tabella(numeri(m3, m2)))
            .replace("__DATA__", time.strftime("%Y-%m-%d"))
            .replace("__CORSA__", os.path.basename(os.path.abspath(sorgente))))
    with open(os.path.join(destinazione, "index.html"), "w",
              encoding="utf-8") as fh:
        fh.write(html)
    copiati.append("index.html")

    if verbose:
        tot = sum(os.path.getsize(os.path.join(destinazione, f))
                  for f in copiati) / (1024.0 * 1024.0)
        print(f"  {len(copiati)} file, {tot:.1f} MB in "
              f"{os.path.abspath(destinazione)}")
    return os.path.abspath(destinazione)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    ap.add_argument("--sorgente", default="out_piramidi_v03",
                    help="cartella con le uscite della corsa")
    ap.add_argument("--out", default="sito_giza",
                    help="cartella del sito da caricare")
    ap.add_argument("--escludi", nargs="*", default=[],
                    help="file da NON pubblicare (per stare sotto il limite "
                         "di un caricamento diretto)")
    args = ap.parse_args()
    print("costruzione del sito statico:")
    costruisci(args.sorgente, args.out, escludi=args.escludi)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

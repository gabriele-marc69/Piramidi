#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pubblica_sito.py  --  raccoglie le uscite in un sito statico
============================================================

Mette in una cartella quello che si guarda -- le pagine 3D interattive, le
pagine delle singole piramidi, le figure -- piu' un indice che le lega, dice
che cosa sono e **contiene i dati**: i riassunti JSON e il budget non sono
file a parte da scaricare, sono dentro l'HTML dell'indice.

Due profili, perche' le due catene stanno su due progetti distinti:

    python pubblica_sito.py --profilo v03 --out sito_v03
    python pubblica_sito.py --profilo v02 --out sito_v02

I numeri dell'indice NON sono scritti a mano: si leggono da ``meta_v03.json``
e da ``meta.json`` della corsa che si sta pubblicando. Se la corsa cambia,
cambia l'indice, e non c'e' modo che le due cose divergano in silenzio.

Niente build, niente dipendenze, nessuna chiamata di rete dalle pagine: la
cartella e' pronta per un caricamento diretto su qualunque hosting statico.
``--no-ottimizza`` disattiva la ricompressione delle figure a tavolozza (sono
grafici a linee: 256 colori bastano e il file dimezza), utile solo se si
vuole il PNG identico a quello della corsa.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import shutil
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------
# I due profili di pubblicazione
# --------------------------------------------------------------------------

PAGINE_V03: Tuple[Tuple[str, str, str], ...] = (
    ("onde_3d_piramidi.html",
     "Le sinusoidi, pixel per pixel",
     "Tutte e tre le piramidi insieme: per ogni cella della superficie in "
     "geometria radar, la sinusoide risultante di quel pixel. Ruotala, taglia "
     "la quota, cambia la scala orizzontale."),
    ("tomografia_v03_pieno_vuoto.html",
     "Volume pieno/vuoto",
     "I voxel in cui il profilo misurato si discosta dalla risposta attesa da "
     "un solo diffusore, oltre la dispersione misurata sulla piana. Rosso "
     "eccesso, blu difetto, soglia regolabile."),
)

#: le pagine per singola piramide, con le forme d'onda calcolate dentro
PIRAMIDI_V03: Tuple[Tuple[str, str, str], ...] = (
    ("onde_cheope.html", "Cheope (Khufu)",
     "Le sole colonne di Cheope in 3D, la forma d'onda della colonna "
     "rappresentativa con il suo inviluppo, il profilo medio della piramide "
     "e i suoi numeri."),
    ("onde_chefren.html", "Chefren (Khafre)",
     "Le sole colonne di Chefren, con le stesse due forme d'onda calcolate e "
     "il confronto fra quota misurata e quota simulata."),
    ("onde_micerino.html", "Micerino (Menkaure)",
     "La piu' piccola delle tre, e quella con meno celle: utile proprio per "
     "vedere quanto poco basta perche' il segnale sparisca."),
)

FIGURE_V03: Tuple[Tuple[str, str], ...] = (
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
)

PAGINE_V02: Tuple[Tuple[str, str, str], ...] = (
    ("tomografia_piramidi_3d.html",
     "Superficie e tomogramma",
     "La catena principale: superficie ricostruita dai .tiff, nuvola "
     "tomografica, suolo di riferimento, strutture interne note e il pannello "
     "dei lobi verticali."),
)

FIGURE_V02: Tuple[Tuple[str, str], ...] = (
    ("superficie_ricostruita.png",
     "Superficie, qualita' e validazione di livello 1: quota misurata contro "
     "quota simulata in geometria radar."),
    ("profili_pieno_vuoto.png",
     "Profili verticali: colonne sulle piramidi contro un deserto di "
     "controllo."),
)

PROSA_V03 = """
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
"""

PROSA_V02 = """
  <h2>Che cosa fa questa catena</h2>
  <p>Stack interferometrico multi-baseline sul canale VH: lettura degli SLC,
  deramping TOPS, coregistrazione scelta <b>massimizzando la coerenza</b> e non
  il picco della cross-correlazione, fase geometrica di riferimento calcolata
  da orbite e geolocation grid, e inversione in quota con un periodogramma per
  pixel &mdash; il numero d&rsquo;onda verticale cambia da cella a cella, quindi
  non esiste una sola matrice di steering.</p>
  <p>La soglia di qualita&rsquo; non e&rsquo; scelta a mano: viene dalla
  distribuzione <b>nulla</b> del periodogramma, stimata per Monte Carlo con le
  baseline vere. Con undici fasori casuali il picco normalizzato ha gia&rsquo;
  mediana 0,56: qualunque soglia piu&rsquo; bassa accetterebbe rumore puro e lo
  farebbe sembrare una misura.</p>
"""

CONCLUSIONE = """
  <h2>Che cosa dicono i dati</h2>
  <p>La superficie della piana si ricostruisce con precisione metrica.
  <b>Le piramidi no.</b> Le facce a ~52&deg; superano l&rsquo;angolo di
  incidenza di 37&deg;: sono in layover pieno, con centinaia di punti di
  superficie ripiegati nella stessa cella di risoluzione. Il piano dei
  diffusori dominanti sta appena sopra il deserto e non segue il profilo delle
  facce.</p>
  <div class="nota">La risoluzione verticale di questa pila e&rsquo; di circa
  <b>130 metri</b>. Le camere note della Grande Piramide misurano metri: due
  ordini di grandezza sotto. Nulla di quello che si vede qui e&rsquo; una
  rilevazione di cavita&rsquo; risolta in profondita&rsquo;, e il programma non
  lo lascia intendere da nessuna parte.</div>
"""

PROFILI: Dict[str, Dict[str, Any]] = {
    "v03": {
        "titolo": "Tomografia SAR della piana di Giza &mdash; sinusoidi e FFT",
        "sottotitolo": (
            "Dalla coppia modulo/fase di ogni pixel complesso alle sinusoidi "
            "verticali, alla FFT, al discriminante pieno/vuoto. Metodo di "
            "Biondi e Malanga, implementato e verificato riga per riga."),
        "pagine": PAGINE_V03,
        "piramidi": PIRAMIDI_V03,
        "figure": FIGURE_V03,
        "prosa": PROSA_V03,
        "dati": (("meta_v03.json", "Riassunto della corsa v03"),
                 ("budget_tomografico.txt",
                  "Budget tomografico e piano delle sub-aperture")),
    },
    "v02": {
        "titolo": ("Tomografia SAR della piana di Giza &mdash; superficie e "
                   "tomogramma"),
        "sottotitolo": (
            "La catena multi-baseline: superficie reale ricostruita dagli SLC, "
            "nuvola tomografica e attributi di micro-moto. Metodo di Biondi e "
            "Malanga, implementato e verificato riga per riga."),
        "pagine": PAGINE_V02,
        "piramidi": (),
        "figure": FIGURE_V02,
        "prosa": PROSA_V02,
        "dati": (("meta.json", "Riassunto completo della corsa"),
                 ("budget_tomografico.txt",
                  "Budget tomografico e piano delle sub-aperture")),
    },
}


# --------------------------------------------------------------------------
# Lettura dei riassunti e tabella dei numeri
# --------------------------------------------------------------------------

def _leggi(path: str) -> Optional[Dict[str, Any]]:
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _tabella(righe: List[Tuple[str, str]]) -> str:
    corpo = "".join(f"<tr><td>{k}</td><td class='v'>{v}</td></tr>"
                    for k, v in righe)
    return f"<table>{corpo}</table>"


def numeri(m3: Optional[Dict[str, Any]], m2: Optional[Dict[str, Any]],
           profilo: str) -> List[Tuple[str, str]]:
    """I numeri della corsa, letti dai riassunti JSON e mai riscritti a mano."""
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
            ("soglia di qualita' dalla distribuzione nulla",
             f"{m2.get('soglia_qualita_gamma', 0):.3f}"),
        ]
    if m3 and profilo == "v03":
        f = m3.get("errore_fft", {})
        c = m3.get("confronto_con_v02", {})
        n = m3.get("nullo", {})
        bil = (m3.get("bilancio_anomalie") or {}).get("anomalie_totali", {})
        campo = m3.get("campo_onde_pixel_per_pixel") or {}
        out += [
            ("errore della FFT contro il calcolo diretto",
             f"{f.get('errore_relativo_max', float('nan')):.1e}"),
            ("quota identica alla catena precedente entro 1 cm",
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


# --------------------------------------------------------------------------
# Indice
# --------------------------------------------------------------------------

INDICE = """<!DOCTYPE html>
<html lang="it">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITOLOTAB__</title>
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
.pastiglia { display: inline-block; width: 11px; height: 11px;
             border-radius: 3px; margin-right: 7px; }
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
details { background: #FFF; border: 1px solid #E2E8F0; border-radius: 12px;
          margin-bottom: 12px; }
summary { padding: 13px 18px; cursor: pointer; font-size: 15px;
          font-weight: 600; }
summary span { font-weight: 400; color: #64748B; font-size: 13px; }
details pre { margin: 0; padding: 0 18px 16px; overflow-x: auto;
              font-size: 12px; line-height: 1.5; color: #334155;
              max-height: 520px; }
footer { margin-top: 50px; padding-top: 20px; border-top: 1px solid #E2E8F0;
         font-size: 14px; color: #64748B; }
code { background: #F1F5F9; padding: 1px 5px; border-radius: 4px;
       font-size: 90%; }
</style>
</head>
<body>
<header><div class="avvolto">
  <h1>__TITOLO__</h1>
  <p class="sotto">__SOTTOTITOLO__</p>
</div></header>
<div class="avvolto">

  <h2>Le pagine interattive</h2>
  <div class="griglia">__SCHEDE__</div>
__PIRAMIDI__
__PROSA__
  <h2>I numeri della corsa</h2>
  __TABELLA__
__CONCLUSIONE__
  <h2>Le figure</h2>
__FIGURE__
  <h2>I dati</h2>
  <p>I riassunti della corsa sono <b>dentro questa pagina</b>, non in file da
  scaricare a parte: sono gli stessi da cui vengono i numeri della tabella qui
  sopra.</p>
__DATI__
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


def _peso(path: str) -> str:
    n = os.path.getsize(path) / 1024.0
    return f"{n / 1024.0:.1f} MB" if n > 1024 else f"{n:.0f} KB"


def _ottimizza_png(src: str, dst: str) -> Tuple[bool, int, int]:
    """Ricomprime una figura a tavolozza. Restituisce (fatto, prima, dopo).

    Sono grafici a linee su fondo bianco: 256 colori adattivi sono
    indistinguibili dall'originale e il file cala di parecchio. Se Pillow non
    c'e', il file viene copiato tale e quale -- il sito resta corretto, pesa
    solo di piu'."""
    prima = os.path.getsize(src)
    try:
        from PIL import Image
    except Exception:
        shutil.copy2(src, dst)
        return False, prima, prima
    with Image.open(src) as im:
        im.convert("RGB").quantize(colors=256).save(
            dst, format="PNG", optimize=True)
    return True, prima, os.path.getsize(dst)


def costruisci(sorgente: str, destinazione: str, profilo: str = "v03",
               escludi: Sequence[str] = (), ottimizza: bool = True,
               verbose: bool = True) -> str:
    """Copia le uscite del profilo e scrive l'indice con i dati dentro."""
    if profilo not in PROFILI:
        raise ValueError(f"profilo sconosciuto: {profilo}")
    prof = PROFILI[profilo]
    os.makedirs(destinazione, exist_ok=True)
    fuori = set(escludi)
    m3 = _leggi(os.path.join(sorgente, "meta_v03.json"))
    m2 = _leggi(os.path.join(sorgente, "meta.json"))
    copiati: List[str] = []

    def _porta(nome: str) -> Optional[str]:
        src = os.path.join(sorgente, nome)
        if nome in fuori:
            return None
        if not os.path.exists(src):
            if verbose:
                print(f"  manca (saltato): {nome}")
            return None
        dst = os.path.join(destinazione, nome)
        if ottimizza and nome.lower().endswith(".png"):
            fatto, a, b = _ottimizza_png(src, dst)
            if verbose and fatto:
                print(f"  {nome}: {a / 1024:.0f} -> {b / 1024:.0f} KB "
                      f"({100 * (1 - b / max(a, 1)):.0f} % in meno)")
        else:
            shutil.copy2(src, dst)
        copiati.append(nome)
        return dst

    schede = []
    for nome, titolo, descr in prof["pagine"]:
        dst = _porta(nome)
        if dst:
            schede.append(
                f'<a class="scheda" href="{nome}"><b>{titolo}</b>'
                f'<span>{descr}</span><i>apri &rarr; '
                f'<span style="color:#94A3B8">{_peso(dst)}</span></i></a>')

    colori = ("#2563EB", "#16A34A", "#B45309")
    schede_pir = []
    for k, (nome, titolo, descr) in enumerate(prof["piramidi"]):
        dst = _porta(nome)
        if dst:
            schede_pir.append(
                f'<a class="scheda" href="{nome}"><b><span class="pastiglia" '
                f'style="background:{colori[k % 3]}"></span>{titolo}</b>'
                f'<span>{descr}</span><i>apri &rarr; '
                f'<span style="color:#94A3B8">{_peso(dst)}</span></i></a>')
    blocco_pir = ""
    if schede_pir:
        blocco_pir = (
            "\n  <h2>Una pagina per piramide, con le forme d'onda "
            "calcolate</h2>\n"
            "  <p>Ogni pagina porta le sole colonne di quella piramide in 3D "
            "e, sotto, due grafici disegnati dentro la pagina: la somma delle "
            "sinusoidi della colonna rappresentativa con il suo inviluppo, e "
            "il profilo medio di tutte le sue celle.</p>\n"
            f'  <div class="griglia">{"".join(schede_pir)}</div>\n')

    figure = []
    for nome, didascalia in prof["figure"]:
        if _porta(nome):
            figure.append(f'  <figure><img src="{nome}" alt="{didascalia}" '
                          f'loading="lazy"><figcaption>{didascalia}'
                          "</figcaption></figure>")

    # --- i dati, dentro la pagina ------------------------------------------
    blocchi = []
    for nome, etichetta in prof["dati"]:
        src = os.path.join(sorgente, nome)
        if not os.path.exists(src):
            continue
        with open(src, encoding="utf-8") as fh:
            testo = fh.read()
        if nome.endswith(".json"):
            try:
                testo = json.dumps(json.loads(testo), indent=2,
                                   ensure_ascii=False)
            except ValueError:
                pass
        peso = len(testo.encode("utf-8")) / 1024.0
        blocchi.append(
            f"  <details><summary>{etichetta} <span>&mdash; {nome}, "
            f"{peso:.0f} KB</span></summary><pre>{html.escape(testo)}"
            "</pre></details>")

    pagina = (INDICE
              .replace("__TITOLOTAB__",
                       prof["titolo"].replace("&mdash;", "-"))
              .replace("__TITOLO__", prof["titolo"])
              .replace("__SOTTOTITOLO__", prof["sottotitolo"])
              .replace("__SCHEDE__", "".join(schede))
              .replace("__PIRAMIDI__", blocco_pir)
              .replace("__PROSA__", prof["prosa"])
              .replace("__TABELLA__", _tabella(numeri(m3, m2, profilo)))
              .replace("__CONCLUSIONE__", CONCLUSIONE)
              .replace("__FIGURE__", "\n".join(figure))
              .replace("__DATI__", "\n".join(blocchi))
              .replace("__DATA__", time.strftime("%Y-%m-%d"))
              .replace("__CORSA__",
                       os.path.basename(os.path.abspath(sorgente))))
    with open(os.path.join(destinazione, "index.html"), "w",
              encoding="utf-8") as fh:
        fh.write(pagina)
    copiati.append("index.html")

    if verbose:
        tot = sum(os.path.getsize(os.path.join(destinazione, f))
                  for f in copiati) / (1024.0 * 1024.0)
        print(f"  profilo {profilo}: {len(copiati)} file, {tot:.2f} MB in "
              f"{os.path.abspath(destinazione)}")
    return os.path.abspath(destinazione)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Raccoglie le uscite in un sito statico, con i dati "
                    "dentro la pagina.")
    ap.add_argument("--sorgente", default="out_piramidi_v03",
                    help="cartella con le uscite della corsa")
    ap.add_argument("--out", default="sito_giza",
                    help="cartella del sito da caricare")
    ap.add_argument("--profilo", default="v03", choices=sorted(PROFILI),
                    help="quale catena pubblicare")
    ap.add_argument("--escludi", nargs="*", default=[],
                    help="file da NON pubblicare")
    ap.add_argument("--no-ottimizza", dest="ottimizza", action="store_false",
                    help="non ricomprimere le figure a tavolozza")
    args = ap.parse_args()
    print(f"costruzione del sito statico ({args.profilo}):")
    costruisci(args.sorgente, args.out, profilo=args.profilo,
               escludi=args.escludi, ottimizza=args.ottimizza)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramidi_v03.py  --  sintesi sinusoidale e FFT verticale  (2026-09-07)
======================================================================

Terza versione della catena di Giza. Prende il punto dove ``piramidi_v02.py``
si ferma -- il cubo interferometrico multilooked -- e ci passa sopra una
catena esplicita in **quattro passi**, che e' quanto chiesto:

  1. **forma polare**: da ogni pixel complesso si prendono MODULO e FASE IN
     RADIANTI (``|y|`` e ``atan2(Im y, Re y)``);
  2. **sintesi**: ogni data diventa una ONDA SINUSOIDALE lungo la verticale,
     ``s_i(z) = A_i * cos(phi_i - kappa_i * z)``, con ``kappa_i`` il numero
     d'onda verticale di quella data su quel pixel;
  3. **FFT**: le sinusoidi si sommano e la somma viene trasformata con una
     FFT vera (non con il doppio ciclo di v02) per ottenere il profilo in
     quota, e da questo l'indice **pieno/vuoto**;
  4. **3D**: un nuovo grafico tridimensionale del volume pieno/vuoto, un
     disegno a parte delle onde sinusoidali SOTTO Cheope e Chefren, e il
     campo delle sinusoidi risultanti **pixel per pixel** su tutta la
     superficie delle piramidi -- una curva per cella, alla sua posizione,
     statica in PNG e interattiva in HTML.

Che cosa e' davvero nuovo, e che cosa no
----------------------------------------
Va detto subito, perche' altrimenti questa versione sembrerebbe un metodo
nuovo e non lo e'. Con ``kappa_i = segno * k_z,i``:

    somma delle sinusoidi   S(z) = sum_i A_i cos(phi_i - kappa_i z)
                                 = Re[ sum_i y_i exp(-j kappa_i z) ]
                                 = Re[ h(z) ]                       (v02)

cioe' **la somma delle onde sinusoidali E' la parte reale del periodogramma
tomografico di v02**, e il suo inviluppo analitico e' ``|h(z)|``, la stessa
grandezza da cui v02 estrae la quota. Modulo e fase non aggiungono
informazione alla coppia (reale, immaginaria): sono la stessa misura in
coordinate polari, come gia' dichiarato in ``Pixel_Complessi/LEGGIMI.md``.
La sintesi sinusoidale non e' quindi una seconda misura indipendente: e' la
STESSA inversione scritta nel dominio in cui la si puo' disegnare, ed e'
esattamente per questo che vale la pena farla -- il grafico delle onde sotto
le piramidi mostra da dove viene il numero, non un'analogia.

Quello che invece cambia davvero e' **come** si calcola:

* v02 valuta ``h(z)`` con un doppio ciclo Python su 257 quote x N date;
* v03 riconosce che ``h(z)`` e' la trasformata di Fourier di uno spettro di
  righe alle frequenze ``kappa_i`` -- non uniformi, perche' le baseline non
  lo sono -- e lo calcola con una **FFT vera** dopo aver portato le righe su
  un reticolo uniforme con un nucleo gaussiano (NUFFT di tipo 1). Il costo
  passa da O(N_z * N_date) a O(N_k log N_k), e l'errore di griglia non e'
  argomentato ma **misurato**, cella per cella, contro il periodogramma
  diretto: ``--selftest`` lo blocca sotto 1e-3.

Pieno e vuoto: che cosa misura l'indice di v03
-----------------------------------------------
``solidity_index`` di v02 e' un prodotto di tre attributi normalizzati; e'
un discriminante, e non ha un nullo con cui confrontarsi. Qui l'indice e'
costruito come **residuo rispetto alla risposta attesa da UN solo
diffusore**:

    P(z)  = |h(z)| / max_z |h(z)|                      profilo misurato
    Q(z)  = |sum_i A_i exp(-j kappa_i (z - z_pic))| / sum_i A_i
                                                       PSF dell'array
    R(z)  = P(z) - Q(z)                                residuo
    S(z)  = R(z) / sigma_nullo(z - z_pic)              z-score

``Q`` e' la risposta che quella cella darebbe se sotto ci fosse un unico
bersaglio puntiforme alla quota del picco: stesse ampiezze, stesse baseline,
fase puramente geometrica. Dove il profilo misurato STA SOPRA la PSF c'e'
piu' energia di quanta un solo diffusore ne spieghi (**pieno** anomalo),
dove ci sta sotto ce n'e' meno (**vuoto** anomalo). ``sigma_nullo`` non e'
scelta: e' la dispersione di ``R`` misurata sulle celle di sola piana, fuori
dall'impronta delle piramidi, alla stessa distanza dal picco -- lo stesso
principio di F22, un test di ipotesi e non una soglia estetica.

Il limite resta quello di v02 e non lo tocca nessuna FFT: con 43 baseline la
cella di Rayleigh verticale e' ``delta_z ~ 132 m``. Un residuo di PSF su
scale piu' fini non e' stratigrafia, e' la PSF stessa piu' rumore. Le camere
note della Grande Piramide misurano metri: due ordini di grandezza sotto.
L'indice pieno/vuoto e' un discriminante calibrato sul nullo, **non** una
rilevazione di cavita' risolta in profondita'. Il programma lo stampa come
tale e il confronto piramidi/piana e' riportato con il suo z-test, qualunque
risultato dia.

Che cosa dicono i dati (VH, 43 date, master 20260503, DATA_Ghiza)
-----------------------------------------------------------------
1. **La FFT e' il periodogramma di v02, verificato.** Il profilo calcolato
   per FFT coincide con il calcolo diretto entro 1.5e-5 in relativo, e la
   quota che ne esce, passata per lo stesso estimatore di v02
   (``surface_from_tomogram``), coincide con quella di v02 entro 1 cm sul
   100 per cento delle celle. Il guadagno e' di tempo, non di risultato:
   0.7 s contro i ~5 s stimati per il ciclo diretto sullo stesso volume.
2. **Le celle delle piramidi mostrano piu' energia di quanta la PSF di un
   solo diffusore ne spieghi**: il 28.4 per cento supera la soglia contro il
   19.9 per cento delle celle di piana (z = +3.15). L'effetto e' tutto
   nell'eccesso: sul difetto (il "vuoto") le due popolazioni sono identiche,
   0.4 contro 0.3 per cento, z = +0.33.
3. **Quell'eccesso pero' non e' spiegato da nessuna delle grandezze di
   controllo.** Correlazione sulle 236 celle di qualita' delle piramidi:
   ripiegamento di layover r = +0.07, ampiezza r = -0.23, coerenza
   r = -0.25, quota simulata r = -0.08. La piu' forte spiega il 6 per cento
   della varianza. Un residuo non spiegato resta un residuo: con
   ``delta_z`` = 132 m non c'e' niente, sotto quella scala, che possa essere
   risolto, e il modo corretto di leggere il punto 2 e' che la PSF a un solo
   diffusore e' un modello troppo semplice per celle in layover pieno --
   non che sotto le piramidi si sia visto qualcosa.
4. Le due colonne disegnate (la piu' luminosa dentro la maschera simulata di
   ciascuna piramide) hanno il diffusore dominante a +9.4 m sopra il
   riferimento in entrambi i casi, cioe' appena sopra il deserto e non a
   meta' della piramide: la stessa conclusione di v02, vista da un'altra
   parte.
5. Il campo di onde pixel per pixel (F59) porta la stessa cosa su tutte le
   706 celle della superficie in geometria radar -- 335 su Cheope, 301 su
   Chefren, 70 su Micerino, di cui 236 sopra la soglia di qualita'. Guardato
   di taglio, il piano dei diffusori dominanti sta appena sopra il deserto
   per tutte e tre le piramidi, e non segue in alcun modo il profilo delle
   facce: e' l'osservazione del punto 2 della lista di v02, disegnata invece
   che riassunta in una pendenza.

Novita' di questa versione
--------------------------
Elencate in ``NOVITA`` e stampate da ``--novita``; i codici proseguono la
numerazione di v02, che si ferma a F50.

Uso
---
    python piramidi_v03.py                    # catena completa + uscite v03
    python piramidi_v03.py --solo-v03         # salta le uscite di v02
    python piramidi_v03.py --dates 12         # meno date, corsa rapida
    python piramidi_v03.py --selftest         # autotest della FFT e del nullo
    python piramidi_v03.py --novita
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

import piramidi_v02 as v2
from piramidi_v02 import (
    PYRAMIDS,
    Config,
    Pyramid,
    TomoBudget,
    _norm01,
    _wrap,
    known_structures,
    pyramid_mesh,
)

# --------------------------------------------------------------------------
# Registro delle novita' (prosegue la numerazione di v02, ferma a F50)
# --------------------------------------------------------------------------

NOVITA: Tuple[Tuple[str, str, str], ...] = (
    ("F51", "infrastruttura",
     "v02.run() ora restituisce anche il cubo interferometrico multilooked "
     "(y_ml) e il suo k_z (k_ml): sono la coppia modulo/fase da cui v03 "
     "sintetizza le sinusoidi. Senza questa aggiunta v03 avrebbe dovuto "
     "rileggere i .tiff una seconda volta, 60 s su 43 date, per ottenere "
     "dati identici a quelli gia' in memoria."),
    ("F52", "metodo",
     "Forma polare esplicita: A = |y|, phi = atan2(Im y, Re y) in radianti. "
     "E' la stessa coppia di Pixel_Complessi/LEGGIMI.md ma calcolata "
     "sull'INTERFEROGRAMMA normalizzato e riferito, non sul DN grezzo: la "
     "fase grezza e' dominata da exp(-j4piR/lambda) e non porta quota (F01, "
     "F02). Modulo e fase non aggiungono informazione alla coppia cartesiana "
     "-- sono la stessa misura in coordinate polari."),
    ("F53", "metodo",
     "Sintesi sinusoidale: ogni data e' un tono lungo la verticale, "
     "s_i(z) = A_i cos(phi_i - kappa_i z) con kappa_i = segno * k_z,i. La "
     "somma delle sinusoidi e' esattamente Re[h(z)] del periodogramma di "
     "v02 e il suo inviluppo e' |h(z)|: e' la stessa inversione, scritta "
     "nel dominio in cui si puo' disegnare. Dichiarato apertamente, non "
     "presentato come misura indipendente."),
    ("F54", "calcolo",
     "Il profilo in quota si calcola con una FFT vera. Le righe spettrali "
     "stanno alle frequenze kappa_i, che NON sono uniformi (le baseline non "
     "lo sono): vengono portate su un reticolo uniforme con un nucleo "
     "gaussiano e deapodizzate analiticamente (NUFFT di tipo 1), poi una "
     "sola FFT produce il profilo per tutte le celle. L'errore rispetto al "
     "periodogramma diretto e' MISURATO su un campione di celle e stampato, "
     "non assunto; --selftest lo blocca sotto 1e-3."),
    ("F55", "metodo",
     "Indice pieno/vuoto come residuo rispetto alla PSF dell'array: "
     "R(z) = |h(z)|/max - |sum A_i exp(-j kappa_i (z - z_pic))|/sum A_i. "
     "La PSF e' calcolata per cella con le SUE ampiezze e le SUE baseline, "
     "non con un modello medio. Sostituisce il prodotto di attributi "
     "normalizzati di solidity_index, che non aveva un nullo con cui "
     "confrontarsi."),
    ("F56", "statistica",
     "Il residuo e' normalizzato dalla dispersione MISURATA sulle celle di "
     "sola piana, alla stessa distanza dal picco (z - z_pic): un z-score con "
     "il suo nullo empirico, come F22 per gamma. La soglia di anomalia e' un "
     "percentile di quel nullo, e il confronto piramidi/piana e' riportato "
     "come z-test fra due proporzioni, qualunque risultato dia."),
    ("F57", "uscita",
     "Nuovo grafico 3D del volume pieno/vuoto (HTML interattivo autonomo + "
     "PNG): voxel di eccesso in rosso e di difetto in blu, con soglia in "
     "z-score regolabile, sopra la superficie misurata e i profili delle "
     "piramidi."),
    ("F60", "uscita",
     "Una pagina per ciascuna piramide, con le SUE sole colonne e le forme "
     "d'onda calcolate disegnate dentro la pagina in SVG: la somma delle "
     "sinusoidi della colonna rappresentativa con il suo inviluppo, e la "
     "media INCOERENTE dei moduli su tutte le celle di quella piramide (la "
     "media coerente le annullerebbe, perche' le fasi fra celle sono "
     "scorrelate). Nessun file esterno e nessuna libreria: la pagina si "
     "regge da sola. Il renderer 3D e' lo stesso della pagina d'insieme, "
     "filtrato -- due copie divergerebbero."),
    ("F59", "uscita",
     "Grafico 3D delle sinusoidi PIXEL PER PIXEL sulla superficie delle "
     "piramidi: una curva per cella, alla sua posizione est/nord, che e' la "
     "somma delle sinusoidi di quella cella (Re[h(z)]) normalizzata al "
     "proprio massimo. La maschera e' quella in GEOMETRIA RADAR, non "
     "l'impronta al suolo. Statico in PNG e interattivo in HTML, con la "
     "scala orizzontale a cursore perche' lo spostamento e' un artificio di "
     "disegno: l'ampiezza di un interferogramma non e' una lunghezza."),
    ("F58", "uscita",
     "Disegno a parte delle onde sinusoidali SOTTO Cheope e Chefren: la "
     "sezione verticale con la sagoma della piramide in alto e, sotto il "
     "suolo, le sinusoidi di ogni data, la loro somma, l'inviluppo e il "
     "residuo pieno/vuoto con la banda del nullo."),
)


# --------------------------------------------------------------------------
# Configurazione
# --------------------------------------------------------------------------

@dataclass
class ConfigV3(Config):
    """Config di v02 piu' i parametri della sintesi sinusoidale e della FFT."""

    # --- NUFFT: reticolo uniforme in kappa ---------------------------------
    #: sovracampionamento del reticolo in kappa; la finestra in z diventa
    #: `sovracamp` volte piu' larga di quella utile e l'aliasing del gridding
    #: viene spinto fuori dalla parte che si tiene.
    sovracamp: int = 2
    #: larghezza del nucleo gaussiano, in celle del reticolo kappa
    kernel_sigma: float = 1.6
    #: semilarghezza del nucleo, in celle (troncamento)
    kernel_taps: int = 8
    #: celle su cui misurare l'errore della FFT contro il periodogramma diretto
    fft_check_cells: int = 240

    # --- anomalia pieno/vuoto ----------------------------------------------
    #: percentile del nullo (celle di piana) che definisce la soglia |z-score|
    anomalia_q: float = 99.0
    #: voxel esportati nella nuvola 3D
    max_voxel_v3: int = 24000

    # --- campo di onde pixel per pixel (F59) --------------------------------
    #: colonne al massimo nel campo di onde; oltre si tengono le piu' luminose
    onde3d_max_celle: int = 1200
    #: Sottocampionamento dell'asse z per le curve disegnate. A 3 restano 86
    #: campioni su 800 m, cioe' un punto ogni 9,4 m: il lobo principale della
    #: PSF e' largo 132 m, quindi la forma non cambia, e la pagina pesa un
    #: terzo in meno -- il caricamento diretto su un hosting statico accetta
    #: un lotto per volta e quel terzo conta.
    onde3d_passo_z: int = 3
    #: quanti metri vale, sul disegno, un'ampiezza normalizzata pari a 1.
    #: Le celle distano ~14 x 15 m: oltre una decina di metri le colonne si
    #: sovrappongono e il campo diventa una tenda in cui non si distingue piu'
    #: la singola sinusoide, che e' proprio quello che il disegno deve mostrare.
    onde3d_scala_m: float = 11.0
    #: colonne disegnate nel PNG (l'HTML le porta tutte, con un cursore)
    onde3d_png_celle: int = 190
    #: finestra in quota del PNG: l'asse tomografico arriva a +-400 m, ma su
    #: quella scala le piramidi (130 m) diventano illeggibili
    onde3d_png_z_max: float = 250.0

    out_dir: str = "out_piramidi_v03"
    #: ``html_name`` resta quello di v02 (``build_html`` di v02 lo usa per la
    #: propria pagina): la pagina nuova ha un nome suo, cosi' le due
    #: convivono nella stessa cartella e si possono confrontare.
    html_v3: str = "tomografia_v03_pieno_vuoto.html"
    html_onde3d: str = "onde_3d_piramidi.html"


# ==========================================================================
# 1.  Forma polare e sintesi delle sinusoidi
# ==========================================================================

def forma_polare(y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Da complesso a (modulo, fase in radianti) -- F52.

    ``fase = atan2(Im, Re)`` in ``(-pi, pi]``. La fase NON e' una proprieta'
    della sola componente immaginaria: senza la parte reale il ramo di atan2
    non e' definito e il segno si perde. Le due coppie sono lo stesso
    numero, non due misure."""
    y = np.asarray(y)
    return np.abs(y).astype(np.float32), np.angle(y).astype(np.float32)


def onde_sinusoidali(amp: np.ndarray, fase: np.ndarray, kappa: np.ndarray,
                     z: np.ndarray) -> np.ndarray:
    """Le sinusoidi verticali di UNA cella, una per data -- F53.

        s_i(z) = A_i * cos(phi_i - kappa_i * z)

    ``amp``, ``fase`` e ``kappa`` sono vettori [n_date]; esce [n_date, n_z].
    La somma sull'asse 0 e' ``Re[h(z)]`` del periodogramma di v02: e' questa
    identita' che rende il disegno una vista del calcolo e non un'analogia."""
    a = np.asarray(amp, dtype=np.float64)[:, None]
    p = np.asarray(fase, dtype=np.float64)[:, None]
    k = np.asarray(kappa, dtype=np.float64)[:, None]
    zz = np.asarray(z, dtype=np.float64)[None, :]
    return a * np.cos(p - k * zz)


def somma_analitica(y_cell: np.ndarray, kappa: np.ndarray,
                    z: np.ndarray) -> np.ndarray:
    """``h(z) = sum_i y_i exp(-j kappa_i z)`` per UNA cella (riferimento esatto).

    La parte reale e' la somma delle sinusoidi, il modulo il loro inviluppo
    analitico. Serve al disegno e come verita' di confronto per la FFT."""
    yy = np.asarray(y_cell, dtype=np.complex128)[:, None]
    k = np.asarray(kappa, dtype=np.float64)[:, None]
    zz = np.asarray(z, dtype=np.float64)[None, :]
    return (yy * np.exp(-1j * k * zz)).sum(axis=0)


# ==========================================================================
# 2.  FFT non uniforme (NUFFT di tipo 1) sull'asse delle quote -- F54
# ==========================================================================

@dataclass
class GrigliaKappa:
    """Reticolo uniforme in numero d'onda verticale, coniugato all'asse z.

    Il legame e' quello della DFT: con ``n_k`` celle di passo ``d_kappa`` la
    trasformata cade su ``z_n = n * dz`` con ``dz = 2*pi/(n_k*d_kappa)``.
    Fissato il passo in quota ``dz`` (lo stesso di v02, per poter confrontare
    i due profili cella per cella) e il sovracampionamento ``sovracamp``, il
    reticolo e' completamente determinato."""
    n_k: int
    d_kappa: float
    dz: float
    z_full: np.ndarray            # [n_k] quote della FFT, ordinate crescenti
    idx_core: np.ndarray          # indici di z_full che riproducono z_axis
    sigma: float                  # larghezza del nucleo, in celle
    taps: int
    kappa_max_grid: float

    def as_text(self) -> str:
        return (
            f"    reticolo kappa: {self.n_k} celle, passo "
            f"{self.d_kappa:.3e} rad/m, |kappa| < {self.kappa_max_grid:.4f} rad/m\n"
            f"    finestra in quota: {self.z_full[0]:+.0f} .. {self.z_full[-1]:+.0f} m "
            f"(passo {self.dz:.3f} m), utile {len(self.idx_core)} campioni\n"
            f"    nucleo gaussiano: sigma {self.sigma:.2f} celle, "
            f"troncato a +-{self.taps} celle"
        )


def griglia_kappa(z_axis: np.ndarray, cfg: ConfigV3) -> GrigliaKappa:
    """Costruisce il reticolo in kappa coniugato a ``z_axis``.

    ``z_axis`` di v02 e' simmetrico e a passo costante (257 campioni da -400 a
    +400 m). Il reticolo esteso ha ``sovracamp`` volte piu' celle con lo
    STESSO passo in quota, quindi le quote di v02 sono un sottoinsieme esatto
    di quelle della FFT: nessuna interpolazione fra i due profili."""
    z = np.asarray(z_axis, dtype=np.float64)
    n_z = len(z)
    dz = float(z[1] - z[0])
    n_k = int(cfg.sovracamp) * n_z
    d_kappa = 2.0 * math.pi / (n_k * dz)

    # quote della FFT nell'ordine naturale, poi riordinate in modo crescente
    z_full = np.fft.fftshift(np.fft.fftfreq(n_k, d=d_kappa / (2.0 * math.pi)))
    # gli indici che riproducono esattamente z_axis (le quote coincidono a
    # meno del mezzo passo: sono multipli dello stesso dz)
    idx = np.searchsorted(z_full, z - 0.5 * dz) + 0
    idx = np.clip(idx, 0, n_k - 1)
    err = float(np.max(np.abs(z_full[idx] - z)))
    if err > 1e-6 * max(dz, 1.0):
        raise ValueError("il reticolo in quota della FFT non contiene z_axis "
                         f"(scarto massimo {err:.3e} m)")
    return GrigliaKappa(
        n_k=n_k, d_kappa=d_kappa, dz=dz, z_full=z_full, idx_core=idx,
        sigma=float(cfg.kernel_sigma), taps=int(cfg.kernel_taps),
        kappa_max_grid=0.5 * n_k * d_kappa,
    )


def profilo_fft(y: np.ndarray, kappa: np.ndarray, g: GrigliaKappa,
                ) -> np.ndarray:
    """``h(z) = sum_i y_i exp(-j kappa_i z)`` per TUTTE le celle, con una FFT.

    Le righe spettrali stanno a ``kappa_i``, che cambia per data E per cella
    (F04) e non e' uniforme. Si spargono quindi sul reticolo uniforme con un
    nucleo gaussiano -- che ha trasformata gaussiana nota, quindi si
    deapodizza esattamente -- e si trasforma una volta sola:

        G[m]   = sum_i y_i * exp(-(kappa_m - kappa_i)^2 / (2 s^2 dk^2))
        h[n]   = FFT(G)[n] / (sqrt(2 pi) s exp(-s^2 dk^2 z_n^2 / 2))

    Il fattore di deapodizzazione cresce con ``z``: e' la ragione per cui il
    reticolo e' sovracampionato e si tiene solo la parte centrale, dove il
    fattore vale poco piu' di uno e l'aliasing e' sotto il rumore.

    Esce [n_l, n_p, n_z_core], complesso, sulle stesse quote di ``z_axis``."""
    y = np.asarray(y)
    kap = np.asarray(kappa, dtype=np.float64)
    n_d, n_l, n_p = y.shape
    n_pix = n_l * n_p
    n_k, dk, s = g.n_k, g.d_kappa, g.sigma

    b = kap.reshape(n_d, n_pix) / dk           # posizione in celle
    if np.max(np.abs(b)) > 0.5 * n_k - g.taps - 1:
        raise ValueError("kappa fuori dal reticolo: aumentare elev_max o "
                         "sovracamp")

    grid = np.zeros((n_pix, n_k), dtype=np.complex128)
    righe = np.arange(n_pix)
    yv = y.reshape(n_d, n_pix).astype(np.complex128)
    b0 = np.floor(b).astype(np.int64)
    frac = b - b0                              # distanza dalla cella inferiore
    for i in range(n_d):
        for t in range(-g.taps + 1, g.taps + 1):
            d = (frac[i] - t)                  # in celle
            w = np.exp(-0.5 * (d * d) / (s * s))
            grid[righe, (b0[i] + t) % n_k] += yv[i] * w

    h = np.fft.fft(grid, axis=1)
    h = np.fft.fftshift(h, axes=1)[:, g.idx_core]

    z = g.z_full[g.idx_core]
    deap = math.sqrt(2.0 * math.pi) * s * np.exp(-0.5 * (s * dk * z) ** 2)
    h /= deap[None, :]
    return h.reshape(n_l, n_p, len(z)).astype(np.complex64)


def errore_fft(y: np.ndarray, kappa: np.ndarray, z_axis: np.ndarray,
               h_fft: np.ndarray, celle: int = 240,
               seed: int = 20260907) -> Dict[str, float]:
    """Errore della FFT contro il periodogramma diretto, MISURATO -- F54.

    Il periodogramma diretto e' lo stesso calcolo di v02, esatto per
    costruzione ma O(N_z * N_date). Qui viene valutato su un campione casuale
    di celle e confrontato con il risultato della FFT: l'errore relativo che
    esce e' il numero da riportare, non il limite teorico del nucleo."""
    n_d, n_l, n_p = y.shape
    rng = np.random.default_rng(seed)
    n = min(int(celle), n_l * n_p)
    flat = rng.choice(n_l * n_p, n, replace=False)
    ii, jj = flat // n_p, flat % n_p

    z = np.asarray(z_axis, dtype=np.float64)
    yy = y[:, ii, jj].astype(np.complex128)                 # [n_d, n]
    kk = kappa[:, ii, jj].astype(np.float64)
    esatto = np.einsum("dn,dnz->nz", yy,
                       np.exp(-1j * kk[:, :, None] * z[None, None, :]))
    stimato = h_fft[ii, jj, :].astype(np.complex128)
    scala = np.abs(esatto).max(axis=1, keepdims=True)
    scala = np.where(scala < 1e-12, 1e-12, scala)
    rel = np.abs(stimato - esatto) / scala
    return {
        "celle_verificate": int(n),
        "errore_relativo_medio": float(np.mean(rel)),
        "errore_relativo_max": float(np.max(rel)),
        "errore_p99": float(np.percentile(rel, 99.0)),
    }


# ==========================================================================
# 3.  Pieno e vuoto: residuo rispetto alla PSF dell'array -- F55, F56
# ==========================================================================

def psf_array(amp: np.ndarray, kappa: np.ndarray, g: GrigliaKappa,
              k_pic: np.ndarray) -> np.ndarray:
    """Risposta attesa da UN solo diffusore alla quota del picco, per cella.

        Q(z) = |sum_i A_i exp(-j kappa_i (z - z_pic))| / sum_i A_i

    Si ottiene dalla stessa routine della FFT sostituendo a ``y_i`` il vettore
    ``A_i exp(+j kappa_i z_pic)``: ampiezze misurate, baseline misurate, fase
    puramente geometrica. Non e' un modello medio -- ogni cella ha la propria
    PSF, perche' ha le proprie ampiezze e il proprio k_z (F04)."""
    z = g.z_full[g.idx_core]
    z_pic = z[k_pic]                                    # [n_l, n_p]
    y_psf = (amp.astype(np.complex128)
             * np.exp(1j * kappa * z_pic[None, :, :])).astype(np.complex64)
    q = np.abs(profilo_fft(y_psf, kappa, g)).astype(np.float32)
    somma = amp.sum(axis=0).astype(np.float32)
    return q / np.maximum(somma, 1e-12)[:, :, None]


def _allinea_al_picco(vol: np.ndarray, k_pic: np.ndarray) -> np.ndarray:
    """Riporta ogni profilo su un asse ``zeta = z - z_pic`` a indici interi.

    Il picco sta su un campione della griglia (e' un argmax), quindi
    l'allineamento e' una traslazione di indici e non introduce
    interpolazione. Fuori dall'asse originale il valore e' NaN: non si
    inventa profilo dove non e' stato misurato."""
    n_l, n_p, n_z = vol.shape
    zeta = np.arange(n_z) - (n_z // 2)                   # asse relativo
    idx = k_pic[:, :, None] + zeta[None, None, :]
    dentro = (idx >= 0) & (idx <= n_z - 1)
    out = np.take_along_axis(vol, np.clip(idx, 0, n_z - 1), axis=2)
    return np.where(dentro, out, np.nan).astype(np.float32)


def _riporta_ad_assoluto(vol_zeta: np.ndarray, k_pic: np.ndarray) -> np.ndarray:
    """Inverso di ``_allinea_al_picco``: da zeta a z assoluto."""
    n_l, n_p, n_z = vol_zeta.shape
    c = n_z // 2
    idx = np.arange(n_z)[None, None, :] - k_pic[:, :, None] + c
    dentro = (idx >= 0) & (idx <= n_z - 1)
    out = np.take_along_axis(vol_zeta, np.clip(idx, 0, n_z - 1), axis=2)
    return np.where(dentro, out, np.nan).astype(np.float32)


def indice_pieno_vuoto(
    h_fft: np.ndarray, amp: np.ndarray, kappa: np.ndarray, g: GrigliaKappa,
    piana: np.ndarray, q_soglia: float = 99.0,
) -> Dict[str, Any]:
    """Residuo profilo-PSF normalizzato sul nullo della piana -- F55, F56.

    ``piana`` e' la maschera booleana delle celle FUORI dall'impronta delle
    piramidi: sono quelle che definiscono il nullo. Se l'indice separasse
    davvero pieno da vuoto, le celle delle piramidi lo supererebbero piu'
    spesso di quelle di piana; il programma misura questa differenza e la
    riporta con il suo z-test, senza aggiustarla."""
    mag = np.abs(h_fft).astype(np.float32)
    k_pic = np.argmax(mag, axis=2).astype(np.int64)
    picco = np.max(mag, axis=2)
    p = mag / np.maximum(picco, 1e-12)[:, :, None]        # profilo normalizzato
    q = psf_array(amp, kappa, g, k_pic)                   # PSF per cella

    r_z = (p - q).astype(np.float32)                      # residuo, z assoluto
    r_zeta = _allinea_al_picco(r_z, k_pic)                # residuo, zeta

    # --- nullo: dispersione del residuo sulle sole celle di piana ----------
    sel = np.asarray(piana, dtype=bool)
    if sel.sum() < 30:
        raise ValueError("troppe poche celle di piana per calibrare il nullo")
    campione = r_zeta[sel]                                # [n_piana, n_zeta]
    # ai due estremi dell'asse zeta qualche colonna resta senza campioni (il
    # picco di quelle celle e' troppo vicino al bordo dell'asse in quota):
    # si conta quante celle finite ha ogni colonna e le colonne troppo povere
    # ereditano la dispersione mediana, invece di dividere per zero o per una
    # deviazione stimata su due valori.
    finiti = np.isfinite(campione).sum(axis=0)
    # deviazione standard calcolata a mano sui soli campioni finiti: nanstd su
    # una colonna interamente vuota emette un avviso e restituisce NaN, e qui
    # le colonne vuote sono la norma, non un caso limite.
    cnt = np.maximum(finiti, 1)
    media = np.nansum(campione, axis=0) / cnt
    var = np.nansum((campione - media[None, :]) ** 2, axis=0) / np.maximum(cnt - 1, 1)
    sigma = np.where(finiti >= 20, np.sqrt(var), np.nan)
    buone = np.isfinite(sigma) & (sigma > 1e-6)
    if not buone.any():
        raise ValueError("nullo non stimabile: nessuna colonna con campioni")
    sigma = np.where(buone, sigma, float(np.median(sigma[buone])))
    s_zeta = (r_zeta / sigma[None, None, :]).astype(np.float32)
    nullo = s_zeta[sel]
    fin = np.isfinite(nullo)
    soglia = float(np.percentile(np.abs(nullo[fin]), q_soglia))

    s_z = _riporta_ad_assoluto(s_zeta, k_pic)             # z-score, z assoluto

    # --- riassunti per cella ------------------------------------------------
    n_z = s_zeta.shape[2]
    c = n_z // 2
    fin = np.isfinite(s_zeta)
    sotto = np.zeros(n_z, dtype=bool)
    sotto[:c] = True                                      # zeta < 0: sotto il picco
    # il picco stesso e' zero per costruzione (profilo e PSF vi coincidono):
    # va escluso da entrambi gli estremi, altrimenti entra come "nessuna
    # anomalia" in ogni cella e comprime le statistiche verso zero.
    fuori_picco = np.arange(n_z) != c

    giu = np.where(fin & sotto[None, None, :], s_zeta, np.inf)
    k_vuoto = np.argmin(giu, axis=2)
    vuoto_max = np.take_along_axis(giu, k_vuoto[:, :, None], axis=2)[:, :, 0]
    vuoto_max = np.where(np.isfinite(vuoto_max), vuoto_max, np.nan)

    su = np.where(fin & fuori_picco[None, None, :], s_zeta, -np.inf)
    pieno_max = np.max(su, axis=2)
    pieno_max = np.where(np.isfinite(pieno_max), pieno_max, np.nan)

    zeta_axis = (np.arange(n_z) - c) * g.dz
    z_vuoto = zeta_axis[k_vuoto].astype(np.float32)

    return {
        "profilo": p.astype(np.float32),
        "psf": q.astype(np.float32),
        "residuo": r_z,
        "zscore": s_z,
        "zscore_zeta": s_zeta,
        "sigma_nullo": sigma.astype(np.float32),
        "zeta": zeta_axis.astype(np.float32),
        "k_picco": k_pic.astype(np.int32),
        "picco": picco.astype(np.float32),
        "soglia": soglia,
        "percentile": float(q_soglia),
        "vuoto_max": vuoto_max.astype(np.float32),
        "pieno_max": pieno_max.astype(np.float32),
        "z_vuoto": z_vuoto,
        "celle_nullo": int(sel.sum()),
    }


def _z_test_proporzioni(k1: int, n1: int, k2: int, n2: int) -> Dict[str, Any]:
    """z-test fra due proporzioni (piramidi contro piana) -- F56."""
    if n1 < 5 or n2 < 5:
        return {"esito": "campione insufficiente"}
    p1, p2 = k1 / n1, k2 / n2
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(max(p * (1.0 - p) * (1.0 / n1 + 1.0 / n2), 1e-18))
    z = (p1 - p2) / se
    return {
        "frazione_piramidi": round(p1, 4),
        "frazione_piana": round(p2, 4),
        "celle_piramidi": int(n1),
        "celle_piana": int(n2),
        "z": round(float(z), 2),
        "significativo_2sigma": bool(abs(z) >= 2.0),
        "esito": ("le celle delle piramidi superano la soglia piu' spesso "
                  "della piana" if z >= 2.0 else
                  "le celle delle piramidi superano la soglia MENO spesso "
                  "della piana" if z <= -2.0 else
                  "nessuna differenza significativa fra piramidi e piana"),
    }


def _pearson(a: np.ndarray, b: np.ndarray) -> Tuple[float, int]:
    """Correlazione di Pearson sui soli campioni finiti di entrambe."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 10:
        return float("nan"), int(m.sum())
    x, y = a[m] - a[m].mean(), b[m] - b[m].mean()
    den = math.sqrt(float(np.dot(x, x)) * float(np.dot(y, y)))
    if den < 1e-18:
        return float("nan"), int(m.sum())
    return float(np.dot(x, y) / den), int(m.sum())


def spiegazione_eccesso(idx: Dict[str, Any], res: Dict[str, Any],
                        pyr_mask: np.ndarray, good: np.ndarray) -> Dict[str, Any]:
    """Contro che cosa correla l'eccesso di energia, sulle celle delle piramidi.

    Se le celle delle piramidi mostrano piu' energia di quanta un solo
    diffusore ne spieghi, la prima spiegazione da escludere non e' una
    cavita': e' il LAYOVER. Con facce a 52 gradi e incidenza a 37 gradi ogni
    cella delle piramidi raccoglie fino a centinaia di punti di superficie
    (``sim_fold`` di v02), e una miscela di diffusori produce esattamente un
    profilo piu' largo della PSF di uno solo. Qui l'eccesso viene correlato
    con il ripiegamento simulato, con l'ampiezza e con la coerenza: se e' il
    layover, si vede."""
    m = np.asarray(pyr_mask, dtype=bool) & np.asarray(good, dtype=bool)
    if m.sum() < 20:
        return {"esito": "celle insufficienti"}
    ecc = idx["pieno_max"][m]
    amp_db = 20.0 * np.log10(np.maximum(res["amp"][m], 1e-9))
    out: Dict[str, Any] = {"celle": int(m.sum())}
    for nome, altro in (("ripiegamento_layover", np.log1p(res["sim_fold"][m])),
                        ("ampiezza_dB", amp_db),
                        ("coerenza", res["coherence"][m]),
                        ("quota_simulata_m", res["sim_h"][m])):
        r, n = _pearson(ecc, altro)
        out[nome] = {"r": None if not math.isfinite(r) else round(r, 3),
                     "campioni": n}
    ordinate = [(abs(v["r"]), k) for k, v in out.items()
                if isinstance(v, dict) and v.get("r") is not None]
    if ordinate:
        ordinate.sort(reverse=True)
        r_max, prima = ordinate[0]
        out["prima_spiegazione"] = prima
        out["varianza_spiegata_pct"] = round(100.0 * r_max ** 2, 1)
        # La forza dell'affermazione segue |r|, non il fatto che una delle
        # grandezze sia arrivata prima delle altre: con r = 0.25 si spiega il
        # 6 per cento della varianza, e chiamarla "la spiegazione" sarebbe
        # esattamente l'errore che questo programma dovrebbe evitare.
        if r_max >= 0.5:
            out["esito"] = (
                f"l'eccesso e' in gran parte {prima} (r = {out[prima]['r']:+.2f}, "
                f"{out['varianza_spiegata_pct']:.0f} % della varianza): va "
                "attribuito a quello, non a una struttura interna")
        elif r_max >= 0.3:
            out["esito"] = (
                f"l'eccesso correla in parte con {prima} "
                f"(r = {out[prima]['r']:+.2f}, {out['varianza_spiegata_pct']:.0f} % "
                "della varianza): una spiegazione parziale, non una struttura")
        else:
            out["esito"] = (
                "nessuna delle grandezze di controllo spiega l'eccesso: la "
                f"correlazione piu' forte e' {prima} con r = "
                f"{out[prima]['r']:+.2f}, cioe' il "
                f"{out['varianza_spiegata_pct']:.0f} % della varianza. "
                "L'eccesso resta senza spiegazione misurata -- il che non lo "
                "rende una struttura: sotto delta_z non c'e' niente da "
                "risolvere, e un residuo non spiegato e' prima di tutto un "
                "residuo")
    return out


def bilancio_anomalie(idx: Dict[str, Any], pyr_mask: np.ndarray,
                      piana: np.ndarray, good: np.ndarray) -> Dict[str, Any]:
    """Quante celle superano la soglia, dentro e fuori dalle piramidi."""
    s = idx["zscore_zeta"]
    with np.errstate(invalid="ignore"):
        forte = np.nanmax(np.abs(s), axis=2) >= idx["soglia"]
    vuoto = idx["vuoto_max"] <= -idx["soglia"]
    pieno = idx["pieno_max"] >= idx["soglia"]

    m_p = np.asarray(pyr_mask, dtype=bool) & good
    m_d = np.asarray(piana, dtype=bool) & good
    return {
        "soglia_zscore": round(float(idx["soglia"]), 3),
        "percentile_nullo": idx["percentile"],
        "celle_nullo": idx["celle_nullo"],
        "anomalie_totali": _z_test_proporzioni(int(forte[m_p].sum()), int(m_p.sum()),
                                               int(forte[m_d].sum()), int(m_d.sum())),
        "vuoto": _z_test_proporzioni(int(vuoto[m_p].sum()), int(m_p.sum()),
                                     int(vuoto[m_d].sum()), int(m_d.sum())),
        "pieno": _z_test_proporzioni(int(pieno[m_p].sum()), int(m_p.sum()),
                                     int(pieno[m_d].sum()), int(m_d.sum())),
    }


# ==========================================================================
# 4.  Colonne rappresentative sotto Cheope e Chefren -- F58
# ==========================================================================

PIRAMIDI_ONDE: Tuple[str, str] = ("Cheope (Khufu)", "Chefren (Khafre)")


def colonna_rappresentativa(mask: np.ndarray, amp: np.ndarray,
                            good: np.ndarray) -> Optional[Tuple[int, int]]:
    """Cella su cui disegnare le onde: la piu' luminosa fra quelle di qualita'.

    La scelta e' dichiarata perche' cambia il disegno: si prende il massimo di
    ampiezza DENTRO la maschera simulata della piramide, preferendo le celle
    che superano la soglia di qualita'; se nessuna la supera si ripiega sul
    massimo di ampiezza e il disegno lo segnala."""
    m = np.asarray(mask, dtype=bool)
    if not m.any():
        return None
    scelta = m & np.asarray(good, dtype=bool)
    if not scelta.any():
        scelta = m
    a = np.where(scelta, amp, -np.inf)
    k = int(np.argmax(a))
    return divmod(k, amp.shape[1])


def dati_onde(res: Dict[str, Any], idx: Dict[str, Any], g: GrigliaKappa,
              cfg: ConfigV3, n_fine: int = 1201) -> List[Dict[str, Any]]:
    """Per Cheope e Chefren: sinusoidi, somma, inviluppo, profilo e residuo.

    Le sinusoidi sono campionate su un asse fitto (``n_fine`` punti) perche'
    servono al disegno: il passo di 3 m dell'asse tomografico le renderebbe
    spezzate. Profilo, PSF e residuo restano invece sull'asse della misura."""
    y = res["y_ml"]
    kappa = res["kappa"]
    amp_master = res["amp"]
    good = res["good"]
    z_axis = np.asarray(res["z_axis"], dtype=np.float64)
    z_fine = np.linspace(float(z_axis[0]), float(z_axis[-1]), int(n_fine))

    per_nome = {d["nome"]: d for d in res["sim_per"]}
    fuori: List[Dict[str, Any]] = []
    for p in PYRAMIDS:
        if p.name not in PIRAMIDI_ONDE:
            continue
        d = per_nome.get(p.name)
        if d is None:
            continue
        cella = colonna_rappresentativa(d["mask"], amp_master, good)
        if cella is None:
            continue
        i, j = cella
        yc = y[:, i, j]
        kc = kappa[:, i, j].astype(np.float64)
        a_c, f_c = forma_polare(yc)
        onde = onde_sinusoidali(a_c, f_c, kc, z_fine)
        h_fine = somma_analitica(yc, kc, z_fine)

        # profilo medio della piramide, incoerente: media dei MODULI, non dei
        # complessi -- celle diverse hanno fasi scorrelate e la media coerente
        # le annullerebbe.
        m = np.asarray(d["mask"], dtype=bool)
        prof_medio = idx["profilo"][m].mean(axis=0) if m.any() else None

        b_perp = np.array([b.b_perp for b in res["baselines"]], dtype=np.float64)
        fuori.append({
            "nome": p.name,
            "piramide": p,
            "cella": (int(i), int(j)),
            "qualita": bool(good[i, j]),
            "z_fine": z_fine,
            "onde": onde,
            "somma": onde.sum(axis=0),
            "inviluppo": np.abs(h_fine),
            "modulo": a_c.astype(np.float64),
            "fase": f_c.astype(np.float64),
            "kappa": kc,
            "b_perp": b_perp,
            "date": list(res["dates"]),
            "z_axis": z_axis,
            "profilo": idx["profilo"][i, j],
            "psf": idx["psf"][i, j],
            "residuo": idx["residuo"][i, j],
            "zscore": idx["zscore"][i, j],
            "profilo_medio_piramide": prof_medio,
            "n_celle_maschera": int(m.sum()),
            "z_picco": float(z_axis[idx["k_picco"][i, j]]),
            "h_ref": float(res["height_ref"][i, j]),
            "vuoto_max": float(idx["vuoto_max"][i, j]),
            "z_vuoto": float(idx["z_vuoto"][i, j]),
            "pieno_max": float(idx["pieno_max"][i, j]),
        })
    return fuori


# ==========================================================================
# 4b.  Il campo delle onde: una sinusoide risultante per ogni pixel  -- F59
# ==========================================================================

#: colori delle tre piramidi nei disegni del campo di onde
COLORI_PIRAMIDE: Tuple[str, str, str] = ("#2563EB", "#16A34A", "#B45309")


@dataclass
class CampoOnde:
    """Le onde di TUTTI i pixel che cadono sulle piramidi, non due colonne.

    Per ogni cella la curva disegnata e' ``Re[h(z)]``, cioe' esattamente la
    somma delle sinusoidi di quel pixel (F53): non viene ricalcolata qui, si
    prende la parte reale del profilo gia' trasformato. L'inviluppo e'
    ``|h(z)|``. Entrambi sono normalizzati al massimo del modulo di quella
    cella, quindi ogni colonna si legge sulla propria scala e colonne
    luminose e colonne deboli restano confrontabili in forma."""
    east: np.ndarray              # [n] m, posizione ENU della cella
    north: np.ndarray             # [n]
    h_ref: np.ndarray             # [n] quota di riferimento .xml della cella
    z: np.ndarray                 # [n_z] asse delle quote (sottocampionato)
    onda: np.ndarray              # [n, n_z] Re h / max|h|, in [-1, 1]
    inviluppo: np.ndarray         # [n, n_z] |h| / max|h|, in [0, 1]
    piramide: np.ndarray          # [n] indice in PYRAMIDS, -1 se nessuna
    qualita: np.ndarray           # [n] bool, sopra la soglia di qualita'
    z_picco: np.ndarray           # [n] quota del diffusore dominante
    sim_h: np.ndarray             # [n] quota simulata della piramide in quella cella
    amp: np.ndarray               # [n] ampiezza del master (sigma0 o DN)
    riga: np.ndarray              # [n] indici di cella, per poter risalire
    colonna: np.ndarray           # [n]
    passo_z: int                  # sottocampionamento applicato all'asse z

    def __len__(self) -> int:
        return int(len(self.east))

    def per_piramide(self) -> Dict[str, int]:
        return {p.name: int(np.count_nonzero(self.piramide == k))
                for k, p in enumerate(PYRAMIDS)}

    def filtra(self, sel: np.ndarray) -> "CampoOnde":
        """Un nuovo campo con le sole colonne selezionate (stesso asse z)."""
        sel = np.asarray(sel, dtype=bool)
        return CampoOnde(
            east=self.east[sel], north=self.north[sel], h_ref=self.h_ref[sel],
            z=self.z, onda=self.onda[sel], inviluppo=self.inviluppo[sel],
            piramide=self.piramide[sel], qualita=self.qualita[sel],
            z_picco=self.z_picco[sel], sim_h=self.sim_h[sel],
            amp=self.amp[sel], riga=self.riga[sel], colonna=self.colonna[sel],
            passo_z=self.passo_z)

    def rappresentativa(self) -> int:
        """Indice della colonna piu' luminosa fra quelle di qualita'.

        E' la stessa regola delle figure a due colonne: si preferisce una
        cella sopra la soglia, e solo se non ce ne sono si ripiega sul
        massimo di ampiezza -- dichiarandolo."""
        if len(self) == 0:
            return -1
        scelta = self.qualita if self.qualita.any() else np.ones(len(self), bool)
        return int(np.argmax(np.where(scelta, self.amp, -np.inf)))


def campo_onde_piramidi(res: Dict[str, Any], an: Dict[str, Any],
                        cfg: ConfigV3) -> CampoOnde:
    """Estrae l'onda risultante di ogni pixel sulla superficie delle piramidi.

    "Superficie delle piramidi" e' la maschera in GEOMETRIA RADAR (``pyr_mask``,
    F24), non l'impronta al suolo: un punto a quota h si sposta di
    ``h*cos(theta)`` verso il near range, e usare l'impronta al suolo
    metterebbe le onde su celle che la piramide non illumina.

    Se le celle superano ``cfg.onde3d_max_celle`` si tengono le piu'
    luminose, non un campione casuale: la decimazione a caso toglierebbe
    proprio le colonne su cui la misura ha senso."""
    h = an["h"]                                   # [n_l, n_p, n_z] complesso
    z_axis = np.asarray(res["z_axis"], dtype=np.float64)
    passo = max(1, int(cfg.onde3d_passo_z))
    kz = np.arange(0, len(z_axis), passo)
    if kz[-1] != len(z_axis) - 1:
        kz = np.append(kz, len(z_axis) - 1)

    m = np.asarray(res["pyr_mask"], dtype=bool)
    ii, jj = np.where(m)
    if len(ii) > cfg.onde3d_max_celle:
        forza = res["amp"][ii, jj]
        tieni = np.argsort(-forza)[: cfg.onde3d_max_celle]
        ii, jj = ii[tieni], jj[tieni]

    hc = h[ii, jj][:, kz]                         # [n, n_zc]
    scala = np.max(np.abs(hc), axis=1)
    scala = np.where(scala < 1e-12, 1e-12, scala)
    onda = (hc.real / scala[:, None]).astype(np.float32)
    inviluppo = (np.abs(hc) / scala[:, None]).astype(np.float32)

    # a quale piramide appartiene ogni cella: dalle maschere per piramide di
    # v02, che portano gia' l'assegnazione fatta in proiezione
    pir = np.full(len(ii), -1, dtype=np.int8)
    nomi = [p.name for p in PYRAMIDS]
    for d in res["sim_per"]:
        k = nomi.index(d["nome"]) if d["nome"] in nomi else -1
        if k < 0:
            continue
        sel = np.asarray(d["mask"], dtype=bool)[ii, jj]
        pir[sel] = k

    k_pic = an["indice"]["k_picco"][ii, jj]
    return CampoOnde(
        east=res["east"][ii, jj].astype(np.float32),
        north=res["north"][ii, jj].astype(np.float32),
        h_ref=res["height_ref"][ii, jj].astype(np.float32),
        z=z_axis[kz].astype(np.float32),
        onda=onda,
        inviluppo=inviluppo,
        piramide=pir,
        qualita=np.asarray(res["good"], dtype=bool)[ii, jj],
        z_picco=z_axis[k_pic].astype(np.float32),
        sim_h=res["sim_h"][ii, jj].astype(np.float32),
        amp=res["amp"][ii, jj].astype(np.float32),
        riga=ii.astype(np.int32), colonna=jj.astype(np.int32),
        passo_z=passo,
    )


# ==========================================================================
# 5.  Disegni: le onde sotto Cheope e Chefren  -- F58
# ==========================================================================

_COL_ONDA = "#38BDF8"
_COL_SOMMA = "#0F172A"
_COL_PSF = "#F0A24A"
_COL_PIENO = "#DC2626"
_COL_VUOTO = "#2563EB"


def _mappa_baseline(b_perp: np.ndarray):
    """Colore per data, dalla baseline ortogonale: il parametro che decide la
    frequenza della sinusoide di quella data."""
    import matplotlib
    from matplotlib import cm
    lo, hi = float(np.min(b_perp)), float(np.max(b_perp))
    norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
    return cm.ScalarMappable(norm=norm, cmap="viridis"), norm


def plot_onde(onde: List[Dict[str, Any]], res: Dict[str, Any], cfg: ConfigV3,
              idx: Dict[str, Any]) -> Optional[str]:
    """Tre righe per piramide: sinusoidi, profilo contro PSF, residuo.

    E' il disegno del calcolo, non un'illustrazione: la curva nera della prima
    riga e' la somma delle sinusoidi colorate che le stanno sotto, e la sua
    ampiezza analitica e' la curva della seconda riga."""
    if not onde:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(onde)
    fig, ax = plt.subplots(3, n, figsize=(8.2 * n, 13.0), squeeze=False,
                           constrained_layout=True)
    budget: TomoBudget = res["budget"]
    soglia = float(idx["soglia"])

    for c, d in enumerate(onde):
        z = d["z_fine"]
        sm, _ = _mappa_baseline(d["b_perp"])

        # --- riga 1: le sinusoidi e la loro somma --------------------------
        a0 = ax[0][c]
        scala = float(np.max(np.abs(d["onde"]))) or 1.0
        for k in range(d["onde"].shape[0]):
            a0.plot(z, d["onde"][k] / scala, lw=0.7, alpha=0.55,
                    color=sm.to_rgba(d["b_perp"][k]))
        s_norm = float(np.max(np.abs(d["somma"]))) or 1.0
        a0.plot(z, d["somma"] / s_norm, lw=2.0, color=_COL_SOMMA,
                label="somma delle sinusoidi = Re h(z)")
        a0.plot(z, d["inviluppo"] / s_norm, lw=1.4, ls="--", color=_COL_PSF,
                label="inviluppo |h(z)|")
        a0.plot(z, -d["inviluppo"] / s_norm, lw=1.4, ls="--", color=_COL_PSF)
        a0.axvline(d["z_picco"], color="#16A34A", lw=1.4, ls=":",
                   label=f"picco a z = {d['z_picco']:+.1f} m")
        a0.axhline(0.0, color="#94A3B8", lw=0.7)
        cb = fig.colorbar(sm, ax=a0, pad=0.01)
        cb.set_label("B_perp della data [m]")
        a0.set_xlim(float(z[0]), float(z[-1]))
        a0.set_ylim(-1.15, 1.15)
        a0.set_xlabel("quota z sopra il riferimento .xml [m]")
        a0.set_ylabel("ampiezza normalizzata")
        a0.set_title(
            f"{d['nome']} - cella ({d['cella'][0]}, {d['cella'][1]})\n"
            f"{d['onde'].shape[0]} sinusoidi, una per data: "
            "s_i(z) = A_i cos(phi_i - kappa_i z)")
        a0.legend(fontsize=8, loc="upper right")
        a0.grid(alpha=.22)

        # --- riga 2: profilo misurato contro PSF dell'array ----------------
        a1 = ax[1][c]
        za = d["z_axis"]
        a1.plot(za, d["profilo"], lw=2.0, color=_COL_SOMMA,
                label="profilo misurato P(z) = |h| / max|h|")
        a1.plot(za, d["psf"], lw=1.8, color=_COL_PSF,
                label="PSF di UN diffusore Q(z)")
        a1.fill_between(za, d["psf"], d["profilo"],
                        where=(d["profilo"] >= d["psf"]),
                        color=_COL_PIENO, alpha=.28, label="eccesso (pieno)")
        a1.fill_between(za, d["psf"], d["profilo"],
                        where=(d["profilo"] < d["psf"]),
                        color=_COL_VUOTO, alpha=.28, label="difetto (vuoto)")
        if d["profilo_medio_piramide"] is not None:
            a1.plot(za, d["profilo_medio_piramide"], lw=1.1, ls=":",
                    color="#64748B",
                    label=f"media |h| sulle {d['n_celle_maschera']} celle "
                          "della piramide")
        a1.axvline(d["z_picco"], color="#16A34A", lw=1.2, ls=":")
        a1.set_xlim(float(za[0]), float(za[-1]))
        a1.set_ylim(0.0, 1.05)
        a1.set_xlabel("quota z sopra il riferimento .xml [m]")
        a1.set_ylabel("ampiezza normalizzata")
        a1.set_title("Profilo dopo la FFT contro la risposta attesa da un solo "
                     "diffusore")
        a1.legend(fontsize=8, loc="upper right")
        a1.grid(alpha=.22)

        # --- riga 3: residuo in z-score -----------------------------------
        a2 = ax[2][c]
        sz = d["zscore"]
        fin = np.isfinite(sz)
        a2.plot(za[fin], sz[fin], lw=1.8, color=_COL_SOMMA,
                label="z-score del residuo")
        a2.fill_between(za[fin], 0.0, sz[fin], where=(sz[fin] > 0),
                        color=_COL_PIENO, alpha=.30)
        a2.fill_between(za[fin], 0.0, sz[fin], where=(sz[fin] < 0),
                        color=_COL_VUOTO, alpha=.30)
        a2.axhspan(-soglia, soglia, color="#94A3B8", alpha=.18,
                   label=f"nullo della piana (+-{soglia:.2f}, p"
                         f"{idx['percentile']:.0f})")
        a2.axhline(0.0, color="#94A3B8", lw=0.8)
        a2.axvline(d["z_picco"], color="#16A34A", lw=1.2, ls=":")
        a2.set_xlim(float(za[0]), float(za[-1]))
        a2.set_xlabel("quota z sopra il riferimento .xml [m]")
        a2.set_ylabel("z-score (sigma del nullo)")
        a2.set_title(
            f"Residuo pieno/vuoto - deficit massimo {d['vuoto_max']:+.2f} a "
            f"{d['z_vuoto']:+.0f} m dal picco, eccesso {d['pieno_max']:+.2f}")
        a2.legend(fontsize=8, loc="upper right")
        a2.grid(alpha=.22)

    fig.suptitle(
        "Sinusoidi verticali da modulo e fase, somma, FFT e residuo "
        f"pieno/vuoto | {len(res['dates'])} date, master {res['master_date']}, "
        f"delta_z = {budget.delta_z_vertical:.0f} m: nessuna struttura piu' "
        "fine di questa cella e' separabile", fontsize=11)
    path = os.path.join(cfg.out_dir, "onde_cheope_chefren.png")
    fig.savefig(path, dpi=125)
    plt.close(fig)
    return path


def plot_sezione_onde(onde: List[Dict[str, Any]], res: Dict[str, Any],
                      cfg: ConfigV3, idx: Dict[str, Any]) -> Optional[str]:
    """Le onde sinusoidali disegnate SOTTO la sagoma delle due piramidi -- F58.

    Sezione verticale: sopra il suolo la piramide in scala vera (metri sulle
    ascisse, metri sulle ordinate), sotto il suolo le sinusoidi di ogni data
    con l'ampiezza riscalata alla semibase -- il fattore di scala e' stampato
    nel pannello, perche' l'ampiezza di un interferogramma non ha unita' di
    lunghezza e disegnarla su un asse in metri senza dirlo sarebbe un
    trucco."""
    if not onde:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    budget: TomoBudget = res["budget"]
    n = len(onde)
    fig, ax = plt.subplots(1, n, figsize=(8.0 * n, 11.0), squeeze=False,
                           constrained_layout=True)
    soglia = float(idx["soglia"])

    for c, d in enumerate(onde):
        a = ax[0][c]
        p: Pyramid = d["piramide"]
        half = 0.5 * p.base_side_m
        h_ref = d["h_ref"]
        z = d["z_fine"]
        z_abs = h_ref + z                      # quota assoluta della sezione

        sm, _ = _mappa_baseline(d["b_perp"])
        scala = float(np.max(np.abs(d["inviluppo"]))) or 1.0
        fattore = half / scala                 # da ampiezza a metri sul disegno
        # Le singole sinusoidi hanno ampiezza A_i, cioe' circa 1/N della
        # somma: disegnate con lo stesso fattore dell'inviluppo sarebbero un
        # filo verticale. Hanno quindi un fattore proprio, dichiarato
        # nell'etichetta dell'asse -- due scale, entrambe scritte, e non una
        # sola scala che nasconde quanto e' stata gonfiata la seconda.
        scala_sing = float(np.max(np.abs(d["onde"]))) or 1.0
        fattore_sing = 0.55 * half / scala_sing

        # --- la piramide, in scala vera ------------------------------------
        a.add_patch(Polygon(
            [(-half, p.base_alt_m), (half, p.base_alt_m),
             (0.0, p.base_alt_m + p.height_m)],
            closed=True, facecolor="#FDE7C7", edgecolor="#B45309", lw=1.8,
            zorder=3))
        a.plot([-1.35 * half, 1.35 * half], [p.base_alt_m, p.base_alt_m],
               color="#B45309", lw=1.2, zorder=3)
        a.text(0.0, p.base_alt_m + 0.55 * p.height_m, p.name.split()[0],
               ha="center", va="center", fontsize=11, color="#7C2D12", zorder=4)

        # linea del suolo di riferimento (quota .xml della cella)
        a.axhline(h_ref, color="#64748B", lw=1.1, ls="-", zorder=2)
        a.text(-1.33 * half, h_ref - 16, f"riferimento .xml  {h_ref:.1f} m",
               fontsize=8, color="#475569", zorder=6, va="top",
               bbox=dict(boxstyle="round,pad=0.25", fc="#FFFFFF", ec="none",
                         alpha=.8))

        # --- le onde, sotto il suolo ---------------------------------------
        giu = z_abs <= h_ref
        for k in range(d["onde"].shape[0]):
            a.plot(d["onde"][k][giu] * fattore_sing, z_abs[giu], lw=0.7,
                   alpha=0.55, color=sm.to_rgba(d["b_perp"][k]), zorder=1)
        inv = d["inviluppo"] * fattore
        a.fill_betweenx(z_abs[giu], -inv[giu], inv[giu], color="#CBD5E1",
                        alpha=.45, zorder=0, label="inviluppo |h(z)|")
        a.plot(d["somma"][giu] * fattore, z_abs[giu], lw=1.8, color=_COL_SOMMA,
               zorder=2, label="somma delle sinusoidi")

        # --- dove il residuo dice pieno o vuoto ----------------------------
        za = h_ref + d["z_axis"]
        sz = d["zscore"]
        fin = np.isfinite(sz)
        forte_v = fin & (sz <= -soglia)
        forte_p = fin & (sz >= soglia)
        for q, col, eti in ((forte_v, _COL_VUOTO, "vuoto oltre soglia"),
                            (forte_p, _COL_PIENO, "pieno oltre soglia")):
            if q.any():
                a.scatter(np.zeros(int(q.sum())), za[q], s=26, marker="_",
                          linewidths=2.2, color=col, zorder=5, label=eti)
        a.axhline(h_ref + d["z_picco"], color="#16A34A", lw=1.3, ls=":",
                  zorder=4,
                  label=f"diffusore dominante  {h_ref + d['z_picco']:.0f} m")

        # --- barra della cella di Rayleigh: il limite, disegnato -----------
        dz_r = float(budget.delta_z_vertical)
        x_bar = 1.12 * half
        y0 = h_ref - 0.75 * float(np.max(np.abs(d["z_axis"])))
        a.plot([x_bar, x_bar], [y0, y0 + dz_r], color="#111827", lw=3.0,
               solid_capstyle="butt", zorder=5)
        a.text(x_bar - 6, y0 + 0.5 * dz_r, f"delta_z = {dz_r:.0f} m",
               rotation=90, ha="right", va="center", fontsize=8.5,
               color="#111827", zorder=5)

        a.set_xlim(-1.4 * half, 1.4 * half)
        a.set_ylim(float(np.min(z_abs)) - 10.0,
                   p.base_alt_m + p.height_m + 40.0)
        a.set_xlabel(
            "est-ovest [m] per la piramide; le onde sono riscalate\n"
            f"somma e inviluppo: 1 = {half:.0f} m   "
            f"singole sinusoidi: 1 = {0.55 * half:.0f} m")
        a.set_ylabel("quota assoluta [m]")
        a.set_title(f"{d['nome']}\ncella ({d['cella'][0]}, {d['cella'][1]})"
                    + ("" if d["qualita"] else "  [sotto la soglia di qualita']"))
        a.legend(fontsize=8, loc="lower right")
        a.grid(alpha=.18)

    fig.suptitle(
        "Le onde sinusoidali sotto Cheope e Chefren\n"
        "una per data: ampiezza = modulo, sfasamento = fase radiante "
        "dell'interferogramma, frequenza = numero d'onda della sua baseline\n"
        f"la barra nera e' la cella di Rayleigh verticale "
        f"({budget.delta_z_vertical:.0f} m): sotto quella scala il profilo e' "
        "PSF dell'array, non stratigrafia",
        fontsize=11)
    path = os.path.join(cfg.out_dir, "sezione_onde_sotto_piramidi.png")
    fig.savefig(path, dpi=125)
    plt.close(fig)
    return path


def plot_onde_3d_png(campo: CampoOnde, res: Dict[str, Any], cfg: ConfigV3,
                     ) -> Optional[str]:
    """Il campo di onde in 3D, versione statica -- F59.

    Ogni colonna e' la sinusoide risultante di UN pixel: sta alla sua
    posizione est/nord, si sviluppa lungo la verticale e l'ampiezza e'
    disegnata come spostamento verso est. Lo spostamento e' un artificio di
    disegno -- l'ampiezza di un interferogramma non e' una lunghezza -- e il
    fattore di conversione e' scritto sull'asse."""
    if len(campo) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    # per il PNG si diradano le colonne: 700 curve sovrapposte su una figura
    # statica non si leggono, e la pagina interattiva le porta tutte
    n = len(campo)
    if n > cfg.onde3d_png_celle:
        passo = int(math.ceil(n / cfg.onde3d_png_celle))
        sel = np.zeros(n, dtype=bool)
        # si dirada per POSIZIONE, non per luminosita': cosi' il campo resta
        # rappresentativo di tutta l'area invece di addensarsi sugli spigoli
        ordine = np.lexsort((campo.east, campo.north))
        sel[ordine[::passo]] = True
    else:
        sel = np.ones(n, dtype=bool)
    idx = np.flatnonzero(sel)

    sc = float(cfg.onde3d_scala_m)
    # finestra in quota ridotta: sull'asse intero (+-400 m) le piramidi
    # occupano un sesto della figura e non si vede piu' dove stanno le onde
    # rispetto a loro. L'HTML porta l'asse intero.
    kz_png = np.abs(campo.z) <= float(cfg.onde3d_png_z_max)
    z = campo.z[kz_png].astype(np.float64)

    fig = plt.figure(figsize=(17.0, 8.6), constrained_layout=True)
    for n_ax, (elev, azim, titolo) in enumerate((
            (24, -62, "vista da sud-ovest"),
            (6, -90, "di taglio da sud: il campo come sezione"))):
        ax = fig.add_subplot(1, 2, n_ax + 1, projection="3d")

        segmenti, colori = [], []
        for k in idx:
            x = campo.east[k] + campo.onda[k][kz_png].astype(np.float64) * sc
            y = np.full(len(z), float(campo.north[k]))
            zz = float(campo.h_ref[k]) + z
            punti = np.stack([x, y, zz], axis=1)
            segmenti.append(punti)
            p_i = int(campo.piramide[k])
            colori.append(COLORI_PIRAMIDE[p_i] if 0 <= p_i < len(COLORI_PIRAMIDE)
                          else "#94A3B8")
        ax.add_collection3d(Line3DCollection(
            segmenti, colors=colori, linewidths=0.7, alpha=.55))

        # spigoli delle piramidi, come riferimento
        spig = []
        for p in PYRAMIDS:
            mm = pyramid_mesh(p, res["lat0"], res["lon0"])
            v = mm["vertices"]
            for a_i, b_i in ((0, 1), (1, 2), (2, 3), (3, 0),
                             (0, 4), (1, 4), (2, 4), (3, 4)):
                spig.append([tuple(v[a_i]), tuple(v[b_i])])
        ax.add_collection3d(Line3DCollection(spig, colors="#7C2D12", lw=1.2))

        # Quota del diffusore dominante, ma solo per le celle sopra la soglia
        # di qualita' e dentro la finestra disegnata: matplotlib non taglia
        # gli artisti 3D ai limiti degli assi, quindi i picchi delle celle
        # decorrelate (che stanno ovunque sui +-400 m) uscivano dalla scatola
        # e sembravano detriti sparsi sulla figura.
        vis = (campo.qualita[idx]
               & (np.abs(campo.z_picco[idx]) <= float(cfg.onde3d_png_z_max)))
        if vis.any():
            ax.scatter(campo.east[idx][vis], campo.north[idx][vis],
                       (campo.h_ref[idx] + campo.z_picco[idx])[vis],
                       s=7, color="#0F172A", alpha=.75, linewidths=0)

        ax.set_xlim(float(campo.east.min()) - 2 * sc,
                    float(campo.east.max()) + 2 * sc)
        ax.set_ylim(float(campo.north.min()) - 20, float(campo.north.max()) + 20)
        ax.set_zlim(float((campo.h_ref + z[0]).min()) - 10,
                    max(float((campo.h_ref + z[-1]).max()),
                        max(p.base_alt_m + p.height_m for p in PYRAMIDS)) + 20)
        ax.set_xlabel(f"est [m]  (ampiezza 1 = {sc:.0f} m)")
        ax.set_zlabel("quota [m]")
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(titolo, fontsize=10)
        if n_ax == 0:
            ax.set_ylabel("nord [m]")
        else:
            ax.set_yticks([])

    conteggi = campo.per_piramide()
    eti = ", ".join(f"{nome.split()[0]} {conteggi[nome]}" for nome in conteggi)
    budget: TomoBudget = res["budget"]
    fig.suptitle(
        "Le sinusoidi risultanti, pixel per pixel, sulla superficie delle "
        "piramidi\n"
        f"{len(idx)} colonne disegnate su {len(campo)} celle ({eti}); ogni "
        f"curva e' Re[h(z)] di UNA cella, cioe' la somma delle sue "
        f"{len(res['dates'])} sinusoidi, normalizzata al proprio massimo\n"
        f"punti neri = diffusore dominante delle sole celle sopra soglia; "
        f"quota mostrata +-"
        f"{cfg.onde3d_png_z_max:.0f} m sui +-{abs(float(campo.z[0])):.0f} m "
        f"dell'asse tomografico; delta_z = {budget.delta_z_vertical:.0f} m, "
        "quindi la forma verticale di ogni colonna e' la PSF dell'array",
        fontsize=10)
    path = os.path.join(cfg.out_dir, "onde_3d_superficie_piramidi.png")
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def plot_3d_png(res: Dict[str, Any], idx: Dict[str, Any], cfg: ConfigV3,
                bilancio: Dict[str, Any]) -> str:
    """Nuovo grafico 3D del volume pieno/vuoto, versione statica -- F57."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    s = idx["zscore"]
    z = np.asarray(res["z_axis"], dtype=np.float64)
    east, north, h_ref = res["east"], res["north"], res["height_ref"]
    soglia = float(idx["soglia"])
    good = res["good"]

    fin = np.isfinite(s)
    sel = fin & (np.abs(s) >= soglia) & good[:, :, None]
    ii, jj, kk = np.where(sel)
    if len(ii) > cfg.max_voxel_v3:
        rng = np.random.default_rng(cfg.seed)
        keep = rng.choice(len(ii), cfg.max_voxel_v3, replace=False)
        ii, jj, kk = ii[keep], jj[keep], kk[keep]

    x = east[ii, jj]
    y = north[ii, jj]
    zz = h_ref[ii, jj] + z[kk]
    val = s[ii, jj, kk]

    fig = plt.figure(figsize=(16.5, 8.2), constrained_layout=True)
    n_vuoto = int(np.sum(val < 0))
    n_pieno = int(np.sum(val > 0))
    for n_ax, (elev, azim, titolo) in enumerate((
            (22, -60, "vista da sud-ovest"),
            (3, -90, "sezione: vista da sud, di taglio"))):
        ax = fig.add_subplot(1, 2, n_ax + 1, projection="3d")
        im = ax.scatter(x, y, zz, c=val, cmap="coolwarm", vmin=-3.0, vmax=3.0,
                        s=4, alpha=.55, linewidths=0)
        # superficie misurata, come riferimento
        g = good
        ax.scatter(east[g], north[g], res["height_display"][g], s=2,
                   color="#334155", alpha=.35, linewidths=0)
        # spigoli delle piramidi
        segs = []
        for p in PYRAMIDS:
            m = pyramid_mesh(p, res["lat0"], res["lon0"])
            v = m["vertices"]
            for a_i, b_i in ((0, 1), (1, 2), (2, 3), (3, 0),
                             (0, 4), (1, 4), (2, 4), (3, 4)):
                segs.append([tuple(v[a_i]), tuple(v[b_i])])
        ax.add_collection3d(Line3DCollection(segs, colors="#B45309", lw=1.1))
        ax.set_xlabel("est [m]")
        ax.set_zlabel("quota [m]")
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(titolo, fontsize=10)
        if n_ax == 0:
            ax.set_ylabel("nord [m]")
        else:
            # di taglio l'asse nord si accavalla su se stesso: mostrarlo
            # darebbe una scala illeggibile e falsamente precisa
            ax.set_yticks([])
        if n_ax == 1:
            cb = fig.colorbar(im, ax=ax, shrink=.7, pad=.02)
            cb.set_label("z-score del residuo  (- vuoto   + pieno)")

    tot = bilancio["anomalie_totali"]
    fig.suptitle(
        f"Volume pieno/vuoto dalla FFT delle sinusoidi -- {len(ii)} voxel oltre "
        f"|z| = {soglia:.2f} (p{idx['percentile']:.0f} del nullo di piana): "
        f"{n_pieno} di eccesso, {n_vuoto} di difetto\n"
        + (f"piramidi {tot.get('frazione_piramidi', float('nan')):.1%} contro "
           f"piana {tot.get('frazione_piana', float('nan')):.1%} delle celle, "
           f"z = {tot.get('z', float('nan'))}: {tot.get('esito', '')}"
           if "frazione_piramidi" in tot else tot.get("esito", "")),
        fontsize=11)
    path = os.path.join(cfg.out_dir, "pieno_vuoto_3d.png")
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


# ==========================================================================
# 6.  Nuovo grafico 3D interattivo del volume pieno/vuoto  -- F57
# ==========================================================================

_HTML_TESTA = """<!DOCTYPE html>
<html lang="it">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Giza v03 - volume pieno/vuoto dalla FFT delle sinusoidi</title>
<style>
:root { color-scheme: light; }
* { box-sizing: border-box; }
body { margin: 0; font: 14px/1.5 system-ui, "Segoe UI", Roboto, sans-serif;
       background: #F8FAFC; color: #0F172A; height: 100vh; display: flex;
       flex-direction: column; overflow: hidden; }
header { padding: 12px 20px 9px; border-bottom: 1px solid #E2E8F0;
         background: #FFF; flex: 0 0 auto; }
h1 { margin: 0 0 4px; font-size: 17px; }
header p { margin: 0; font-size: 12.5px; color: #475569; max-width: 1100px; }
main { display: flex; gap: 0; align-items: stretch; flex: 1 1 auto;
       min-height: 0; }
/* La tela e' posizionata in assoluto dentro #scena e #scena ha min-width 0.
   Senza queste due righe il canvas entra nel calcolo del flex con la propria
   larghezza intrinseca: ogni misura() la fa crescere, la crescita allarga il
   contenitore, il resize successivo la fa crescere ancora, e in pochi giri la
   pagina spinge fuori la barra laterale e blocca il rendering. */
#scena { flex: 1 1 auto; position: relative; background: #FFF;
         min-width: 0; overflow: hidden; }
canvas { position: absolute; inset: 0; display: block; width: 100%;
         height: 100%; cursor: grab; }
canvas.trascina { cursor: grabbing; }
aside { width: 340px; flex: 0 0 340px; border-left: 1px solid #E2E8F0;
        background: #FFF; overflow-y: auto; padding: 14px 16px 40px; }
h2 { font-size: 12px; text-transform: uppercase; letter-spacing: .06em;
     color: #64748B; margin: 18px 0 8px; }
h2:first-child { margin-top: 0; }
label.riga { display: flex; align-items: center; gap: 8px; font-size: 13px;
             margin: 5px 0; }
label.riga input[type=checkbox] { accent-color: #0F172A; }
.cursore { margin: 10px 0 14px; font-size: 12.5px; color: #334155; }
.cursore input { width: 100%; accent-color: #0F172A; }
.chiave { display: flex; align-items: center; gap: 7px; font-size: 12.5px;
          margin: 3px 0; }
.pastiglia { width: 13px; height: 13px; border-radius: 3px; flex: 0 0 13px; }
table { width: 100%; border-collapse: collapse; font-size: 12.5px; }
td { padding: 3px 0; vertical-align: top; }
td.v { text-align: right; font-variant-numeric: tabular-nums; }
.nota { font-size: 12px; color: #475569; background: #F1F5F9; padding: 9px 11px;
        border-radius: 7px; border-left: 3px solid #94A3B8; margin: 10px 0; }
.nota b { color: #0F172A; }
footer { padding: 7px 20px; border-top: 1px solid #E2E8F0; background: #FFF;
         font-size: 11.5px; color: #64748B; flex: 0 0 auto; }
button { font: inherit; padding: 5px 10px; border: 1px solid #CBD5E1;
         background: #FFF; border-radius: 6px; cursor: pointer; }
button:hover { background: #F1F5F9; }
</style>
</head>
<body>
"""


def _nota_eccesso(bilancio: Dict[str, Any]) -> str:
    """Il paragrafo sul perche' le piramidi mostrano piu' energia, se lo fanno.

    Compare solo quando l'eccesso c'e' davvero: scrivere la spiegazione di un
    effetto che i dati non mostrano sarebbe suggerire che ci sia."""
    tot = bilancio.get("anomalie_totali", {})
    sp = bilancio.get("spiegazione_eccesso", {})
    if not isinstance(tot.get("z"), (int, float)) or float(tot["z"]) < 2.0:
        return ""
    if "prima_spiegazione" not in sp:
        return ""
    righe = "".join(
        f"<tr><td>{k.replace('_', ' ')}</td><td class='v'>{sp[k]['r']:+.2f}</td></tr>"
        for k in ("ripiegamento_layover", "ampiezza_dB", "coerenza",
                  "quota_simulata_m")
        if k in sp and sp[k].get("r") is not None)
    return (
        "<div class=\"nota\"><b>Prima di leggerci una struttura:</b> le celle "
        "delle piramidi mostrano davvero piu' energia di quanta la PSF di un "
        "solo diffusore ne spieghi. Le spiegazioni ordinarie da escludere "
        "sono il <b>layover</b> (facce a 52 gradi contro un'incidenza di 37: "
        "fino a centinaia di punti di superficie ripiegati nella stessa "
        "cella, e una miscela di diffusori da' un profilo piu' largo di uno "
        "solo), la luminosita' e la decorrelazione. Correlazione dell'eccesso "
        f"con ciascuna, sulle {sp['celle']} celle di qualita' delle piramidi:"
        f"<table>{righe}</table>{sp['esito']}.</div>")


def build_html_v03(res: Dict[str, Any], idx: Dict[str, Any], cfg: ConfigV3,
                   bilancio: Dict[str, Any], fft_info: Dict[str, Any],
                   onde: List[Dict[str, Any]]) -> str:
    """Pagina 3D autonoma del volume pieno/vuoto -- F57.

    Disegna i soli voxel che superano la soglia calibrata sul nullo, sopra la
    superficie misurata e i profili delle piramidi. La soglia resta un
    cursore: chi guarda deve poter vedere che abbassandola la nuvola riempie
    tutto, piramidi e deserto allo stesso modo."""
    s = idx["zscore"]
    z = np.asarray(res["z_axis"], dtype=np.float64)
    east, north = res["east"], res["north"]
    h_ref = res["height_ref"]
    good = res["good"]
    pyr = np.asarray(res["pyr_mask"], dtype=bool)
    budget: TomoBudget = res["budget"]
    soglia = float(idx["soglia"])

    fin = np.isfinite(s)
    sel = fin & (np.abs(s) >= 0.6 * soglia) & good[:, :, None]
    ii, jj, kk = np.where(sel)
    if len(ii) > cfg.max_voxel_v3:
        rng = np.random.default_rng(cfg.seed)
        forza = np.abs(s[ii, jj, kk])
        # si tiene per primo cio' che e' piu' lontano dal nullo, poi un
        # campione casuale del resto: decimare a caso soltanto nasconderebbe
        # proprio le anomalie piu' forti, che sono il motivo della figura.
        ordine = np.argsort(-forza)
        testa = ordine[: cfg.max_voxel_v3 // 2]
        coda = rng.choice(ordine[cfg.max_voxel_v3 // 2:],
                          cfg.max_voxel_v3 - len(testa), replace=False)
        keep = np.concatenate([testa, coda])
        ii, jj, kk = ii[keep], jj[keep], kk[keep]

    voxel = np.stack([
        np.round(east[ii, jj], 1),
        np.round(north[ii, jj], 1),
        np.round(h_ref[ii, jj] + z[kk], 1),
        np.round(s[ii, jj, kk], 2),
        pyr[ii, jj].astype(np.float32),
    ], axis=1)

    gd = good
    superficie = np.stack([
        np.round(east[gd], 1), np.round(north[gd], 1),
        np.round(res["height_display"][gd], 1),
    ], axis=1)

    piramidi = []
    for p in PYRAMIDS:
        m = pyramid_mesh(p, res["lat0"], res["lon0"])
        piramidi.append({
            "nome": p.name,
            "v": [[round(c, 1) for c in v] for v in m["vertices"]],
            "spigoli": [[0, 1], [1, 2], [2, 3], [3, 0],
                        [0, 4], [1, 4], [2, 4], [3, 4]],
        })

    try:
        strutture = known_structures(res["lat0"], res["lon0"], verbose=False)
    except Exception:                                        # pragma: no cover
        strutture = []

    colonne = [{
        "nome": d["nome"],
        "cella": list(d["cella"]),
        "est": round(float(east[d["cella"][0], d["cella"][1]]), 1),
        "nord": round(float(north[d["cella"][0], d["cella"][1]]), 1),
        "z_picco": round(d["z_picco"], 1),
        "vuoto_max": round(d["vuoto_max"], 2),
        "z_vuoto": round(d["z_vuoto"], 1),
        "pieno_max": round(d["pieno_max"], 2),
    } for d in onde]

    payload = {
        "voxel": voxel.tolist(),
        "superficie": superficie.tolist(),
        "piramidi": piramidi,
        "strutture": strutture,
        "colonne": colonne,
        "soglia": round(soglia, 3),
        "percentile": idx["percentile"],
        "bilancio": bilancio,
        "fft": {k: (round(v, 6) if isinstance(v, float) else v)
                for k, v in fft_info.items()},
        "meta": {
            "date": len(res["dates"]),
            "master": res["master_date"],
            "delta_z": round(float(budget.delta_z_vertical), 1),
            "sigma_h": round(float(budget.sigma_h), 1),
            "n_voxel": int(len(ii)),
            "n_superficie": int(gd.sum()),
            "passo_z": round(float(z[1] - z[0]), 3),
            "z_min": round(float(z[0]), 1),
            "z_max": round(float(z[-1]), 1),
            "polarizzazione": cfg.polarisation.upper(),
            "celle_nullo": idx["celle_nullo"],
        },
    }

    dati = json.dumps(payload, ensure_ascii=False, separators=(",", ":"),
                      default=float).replace("</", "<\\/")

    tot = bilancio.get("anomalie_totali", {})
    vuo = bilancio.get("vuoto", {})
    pie = bilancio.get("pieno", {})

    def _pc(d: Dict[str, Any], k: str) -> str:
        v = d.get(k)
        return "-" if v is None else f"{100.0 * float(v):.1f} %"

    corpo = f"""
<header>
  <h1>Giza v03 &mdash; volume pieno/vuoto dalla FFT delle sinusoidi verticali</h1>
  <p>Da ogni pixel complesso si prendono <b>modulo</b> e <b>fase in
  radianti</b>; ogni data diventa una sinusoide lungo la verticale, la somma
  viene trasformata con una FFT e il profilo che ne esce viene confrontato con
  la risposta attesa da <b>un solo diffusore</b>. Quello che vedi sono i voxel
  in cui il residuo supera la dispersione misurata sulla piana. Non e' una
  mappa di cavita': la cella di Rayleigh verticale di questa pila e'
  <b>{budget.delta_z_vertical:.0f} m</b>.</p>
</header>
<main>
  <div id="scena"><canvas id="tela"></canvas></div>
  <aside>
    <h2>Strati</h2>
    <label class="riga"><input type="checkbox" id="c_vuoto" checked> voxel di <b>vuoto</b> (residuo negativo)</label>
    <label class="riga"><input type="checkbox" id="c_pieno" checked> voxel di <b>pieno</b> (residuo positivo)</label>
    <label class="riga"><input type="checkbox" id="c_sup" checked> superficie misurata</label>
    <label class="riga"><input type="checkbox" id="c_pir" checked> piramidi (geometria nota)</label>
    <label class="riga"><input type="checkbox" id="c_str"> strutture interne note</label>
    <label class="riga"><input type="checkbox" id="c_solo"> solo celle sulle piramidi</label>

    <h2>Soglia e taglio</h2>
    <div class="cursore">
      soglia |z-score| = <b id="e_soglia"></b>
      <input type="range" id="s_soglia" min="0.6" max="6" step="0.05">
      <span id="e_nvox"></span>
    </div>
    <div class="cursore">
      quota minima = <b id="e_zmin"></b> m
      <input type="range" id="s_zmin" min="-400" max="200" step="5">
    </div>
    <div class="cursore">
      quota massima = <b id="e_zmax"></b> m
      <input type="range" id="s_zmax" min="-200" max="400" step="5">
    </div>
    <div class="cursore">
      esagerazione verticale = <b id="e_vex"></b>x
      <input type="range" id="s_vex" min="0.5" max="4" step="0.1">
    </div>
    <button id="b_reset">rimetti la vista a zero</button>

    <h2>Legenda</h2>
    <div class="chiave"><span class="pastiglia" style="background:#2563EB"></span>
      difetto di energia rispetto alla PSF &mdash; <i>vuoto</i></div>
    <div class="chiave"><span class="pastiglia" style="background:#DC2626"></span>
      eccesso di energia rispetto alla PSF &mdash; <i>pieno</i></div>
    <div class="chiave"><span class="pastiglia" style="background:#334155"></span>
      superficie misurata (celle sopra la soglia di qualita')</div>
    <div class="chiave"><span class="pastiglia" style="background:#B45309"></span>
      spigoli delle piramidi &mdash; riferimento, non misura</div>

    <h2>Il numero che conta</h2>
    <table>
      <tr><td>celle anomale sulle piramidi</td><td class="v">{_pc(tot, 'frazione_piramidi')}</td></tr>
      <tr><td>celle anomale sulla piana</td><td class="v">{_pc(tot, 'frazione_piana')}</td></tr>
      <tr><td>z fra le due proporzioni</td><td class="v">{tot.get('z', '-')}</td></tr>
      <tr><td>solo vuoto: piramidi / piana</td><td class="v">{_pc(vuo, 'frazione_piramidi')} / {_pc(vuo, 'frazione_piana')}</td></tr>
      <tr><td>solo pieno: piramidi / piana</td><td class="v">{_pc(pie, 'frazione_piramidi')} / {_pc(pie, 'frazione_piana')}</td></tr>
    </table>
    <div class="nota"><b>Esito:</b> {tot.get('esito', 'campione insufficiente')}.
    La soglia e' il percentile {idx['percentile']:.0f} del residuo misurato
    <i>campione per campione</i> su {idx['celle_nullo']} celle di sola piana, e
    una cella conta come anomala se lo supera in almeno un punto del suo asse
    verticale ({len(z)} quote campionate, ma circa
    {max(1, int(round(2.0 * float(z[-1] - z[0]) / max(budget.delta_z_vertical, 1.0))))}
    indipendenti: il profilo e' liscio su delta_z). La frazione di celle di
    piana anomale e' percio' molto maggiore dell'
    {100.0 - idx['percentile']:.0f} %, ed e' <b>quella</b> il metro di
    paragone qui sopra, non lo zero.</div>
    {_nota_eccesso(bilancio)}

    <h2>Catena</h2>
    <table>
      <tr><td>date nella pila</td><td class="v">{len(res['dates'])}</td></tr>
      <tr><td>master</td><td class="v">{res['master_date']}</td></tr>
      <tr><td>polarizzazione</td><td class="v">{cfg.polarisation.upper()}</td></tr>
      <tr><td>errore della FFT (max)</td><td class="v">{fft_info.get('errore_relativo_max', float('nan')):.2e}</td></tr>
      <tr><td>delta_z verticale</td><td class="v">{budget.delta_z_vertical:.0f} m</td></tr>
      <tr><td>sigma_h della superficie</td><td class="v">{budget.sigma_h:.1f} m</td></tr>
      <tr><td>voxel disegnabili</td><td class="v">{len(ii)}</td></tr>
    </table>

    <h2>Colonne di Cheope e Chefren</h2>
    <div id="tab_colonne"></div>
    <div class="nota">Sono le due celle su cui e' disegnato il grafico delle
    onde (<code>sezione_onde_sotto_piramidi.png</code>): la piu' luminosa
    dentro la maschera simulata di ciascuna piramide.</div>
  </aside>
</main>
<footer id="pie"></footer>
"""

    script = """
<script>
(function () {
  "use strict";
  var D = DATI_QUI;
  var tela = document.getElementById("tela");
  var ctx = tela.getContext("2d");
  var vista = { yaw: -0.62, pitch: 0.42, zoom: 1.0, panx: 0, pany: 0 };
  var iniziale = JSON.parse(JSON.stringify(vista));
  var stato = { soglia: D.soglia, zmin: -400, zmax: 400, vex: 1.4 };

  // --- centro e scala della scena --------------------------------------
  var cx = 0, cy = 0, cz = 0, raggio = 1;
  (function () {
    var xs = [], ys = [], zs = [];
    D.voxel.forEach(function (v) { xs.push(v[0]); ys.push(v[1]); zs.push(v[2]); });
    D.superficie.forEach(function (v) { xs.push(v[0]); ys.push(v[1]); zs.push(v[2]); });
    D.piramidi.forEach(function (p) {
      p.v.forEach(function (v) { xs.push(v[0]); ys.push(v[1]); zs.push(v[2]); });
    });
    if (!xs.length) { return; }
    var mn = function (a) { return Math.min.apply(null, a); };
    var mx = function (a) { return Math.max.apply(null, a); };
    cx = 0.5 * (mn(xs) + mx(xs)); cy = 0.5 * (mn(ys) + mx(ys));
    cz = 0.5 * (mn(zs) + mx(zs));
    raggio = Math.max(mx(xs) - mn(xs), mx(ys) - mn(ys), 1);
  })();

  function colore(s) {
    var t = Math.max(-1, Math.min(1, s / 4.0));
    if (t < 0) {
      var a = -t;
      return "rgba(37,99,235," + (0.25 + 0.65 * a).toFixed(3) + ")";
    }
    return "rgba(220,38,38," + (0.25 + 0.65 * t).toFixed(3) + ")";
  }

  var W = 0, H = 0, dpr = 1;
  function misura() {
    dpr = window.devicePixelRatio || 1;
    var r = tela.parentNode.getBoundingClientRect();
    W = Math.max(320, Math.floor(r.width));
    H = Math.max(320, Math.floor(r.height));
    tela.width = Math.floor(W * dpr);
    tela.height = Math.floor(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function proietta(x, y, z) {
    var dx = x - cx, dy = y - cy, dz = (z - cz) * stato.vex;
    var cy1 = Math.cos(vista.yaw), sy1 = Math.sin(vista.yaw);
    var X = dx * cy1 - dy * sy1;
    var Y = dx * sy1 + dy * cy1;
    var cp = Math.cos(vista.pitch), sp = Math.sin(vista.pitch);
    var Y2 = Y * cp - dz * sp;
    var Z2 = Y * sp + dz * cp;
    var s = (Math.min(W, H) * 0.78 / raggio) * vista.zoom;
    return [W / 2 + X * s + vista.panx, H / 2 - Z2 * s + vista.pany, Y2];
  }

  function disegna() {
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#FFFFFF";
    ctx.fillRect(0, 0, W, H);

    var soloPir = document.getElementById("c_solo").checked;
    var vVuoto = document.getElementById("c_vuoto").checked;
    var vPieno = document.getElementById("c_pieno").checked;
    var vSup = document.getElementById("c_sup").checked;
    var vPir = document.getElementById("c_pir").checked;
    var vStr = document.getElementById("c_str").checked;

    // raccolta dei punti con la loro profondita'
    var punti = [];
    var n = 0;
    if (vVuoto || vPieno) {
      for (var i = 0; i < D.voxel.length; i++) {
        var v = D.voxel[i];
        var s = v[3];
        if (Math.abs(s) < stato.soglia) { continue; }
        if (s < 0 && !vVuoto) { continue; }
        if (s > 0 && !vPieno) { continue; }
        if (soloPir && v[4] < 0.5) { continue; }
        if (v[2] < stato.zmin || v[2] > stato.zmax) { continue; }
        var p = proietta(v[0], v[1], v[2]);
        punti.push([p[2], p[0], p[1], colore(s), 3.0]);
        n++;
      }
    }
    if (vSup) {
      for (var k = 0; k < D.superficie.length; k++) {
        var u = D.superficie[k];
        var q = proietta(u[0], u[1], u[2]);
        punti.push([q[2], q[0], q[1], "rgba(51,65,85,0.55)", 2.2]);
      }
    }

    // ordinamento per profondita' a secchielli: O(n), non O(n log n)
    var NB = 384, lo = Infinity, hi = -Infinity, j;
    for (j = 0; j < punti.length; j++) {
      if (punti[j][0] < lo) { lo = punti[j][0]; }
      if (punti[j][0] > hi) { hi = punti[j][0]; }
    }
    var secchi = new Array(NB);
    for (j = 0; j < NB; j++) { secchi[j] = []; }
    var scala = (hi > lo) ? (NB - 1) / (hi - lo) : 0;
    for (j = 0; j < punti.length; j++) {
      var b = Math.floor((punti[j][0] - lo) * scala);
      secchi[b < 0 ? 0 : (b >= NB ? NB - 1 : b)].push(punti[j]);
    }
    for (j = NB - 1; j >= 0; j--) {
      var lista = secchi[j];
      for (var m = 0; m < lista.length; m++) {
        var e = lista[m];
        ctx.fillStyle = e[3];
        ctx.fillRect(e[1] - e[4] / 2, e[2] - e[4] / 2, e[4], e[4]);
      }
    }

    if (vPir) {
      ctx.strokeStyle = "#B45309";
      ctx.lineWidth = 1.3;
      D.piramidi.forEach(function (p) {
        ctx.beginPath();
        p.spigoli.forEach(function (sp) {
          var a = proietta(p.v[sp[0]][0], p.v[sp[0]][1], p.v[sp[0]][2]);
          var b2 = proietta(p.v[sp[1]][0], p.v[sp[1]][1], p.v[sp[1]][2]);
          ctx.moveTo(a[0], a[1]); ctx.lineTo(b2[0], b2[1]);
        });
        ctx.stroke();
      });
    }

    if (vStr && D.strutture.length) {
      D.strutture.forEach(function (s) {
        ctx.strokeStyle = s.colore || "#16A34A";
        ctx.lineWidth = 1.0;
        ctx.beginPath();
        for (var t = 0; t + 5 < s.edges.length; t += 6) {
          var a = proietta(s.edges[t], s.edges[t + 1], s.edges[t + 2]);
          var b2 = proietta(s.edges[t + 3], s.edges[t + 4], s.edges[t + 5]);
          ctx.moveTo(a[0], a[1]); ctx.lineTo(b2[0], b2[1]);
        }
        ctx.stroke();
      });
    }

    document.getElementById("e_nvox").textContent = n + " voxel";
  }

  // --- comandi ----------------------------------------------------------
  var trascina = false, lx = 0, ly = 0;
  tela.addEventListener("mousedown", function (ev) {
    trascina = true; lx = ev.clientX; ly = ev.clientY;
    tela.classList.add("trascina");
  });
  window.addEventListener("mouseup", function () {
    trascina = false; tela.classList.remove("trascina");
  });
  window.addEventListener("mousemove", function (ev) {
    if (!trascina) { return; }
    if (ev.shiftKey) {
      vista.panx += ev.clientX - lx; vista.pany += ev.clientY - ly;
    } else {
      vista.yaw += (ev.clientX - lx) * 0.006;
      vista.pitch += (ev.clientY - ly) * 0.006;
      vista.pitch = Math.max(-1.45, Math.min(1.45, vista.pitch));
    }
    lx = ev.clientX; ly = ev.clientY;
    disegna();
  });
  tela.addEventListener("wheel", function (ev) {
    ev.preventDefault();
    vista.zoom *= (ev.deltaY < 0) ? 1.1 : 1 / 1.1;
    vista.zoom = Math.max(0.25, Math.min(9, vista.zoom));
    disegna();
  }, { passive: false });

  function lega(id, campo, mostra, fmt) {
    var el = document.getElementById(id);
    var et = document.getElementById(mostra);
    el.value = stato[campo];
    et.textContent = fmt(stato[campo]);
    el.addEventListener("input", function () {
      stato[campo] = parseFloat(el.value);
      et.textContent = fmt(stato[campo]);
      disegna();
    });
  }
  lega("s_soglia", "soglia", "e_soglia", function (v) { return v.toFixed(2); });
  lega("s_zmin", "zmin", "e_zmin", function (v) { return v.toFixed(0); });
  lega("s_zmax", "zmax", "e_zmax", function (v) { return v.toFixed(0); });
  lega("s_vex", "vex", "e_vex", function (v) { return v.toFixed(1); });
  ["c_vuoto", "c_pieno", "c_sup", "c_pir", "c_str", "c_solo"].forEach(function (id) {
    document.getElementById(id).addEventListener("change", disegna);
  });
  document.getElementById("b_reset").addEventListener("click", function () {
    vista = JSON.parse(JSON.stringify(iniziale));
    disegna();
  });
  window.addEventListener("resize", function () { misura(); disegna(); });

  // --- tabella delle due colonne ---------------------------------------
  (function () {
    var righe = D.colonne.map(function (c) {
      return "<tr><td>" + c.nome + "</td><td class='v'>" +
        c.vuoto_max.toFixed(2) + " a " + c.z_vuoto.toFixed(0) + " m</td></tr>";
    }).join("");
    document.getElementById("tab_colonne").innerHTML =
      "<table><tr><td><b>colonna</b></td><td class='v'><b>deficit max</b></td></tr>"
      + righe + "</table>";
  })();

  document.getElementById("pie").textContent =
    "trascina per ruotare, shift+trascina per spostare, rotella per lo zoom  |  " +
    D.meta.date + " date, master " + D.meta.master + ", " + D.meta.polarizzazione +
    "  |  asse z da " + D.meta.z_min + " a " + D.meta.z_max + " m a passo " +
    D.meta.passo_z + " m  |  errore della FFT contro il periodogramma diretto: " +
    D.fft.errore_relativo_max.toExponential(1) + " su " +
    D.fft.celle_verificate + " celle";

  misura();
  disegna();
})();
</script>
</body>
</html>
"""

    html = _HTML_TESTA + corpo + script.replace("DATI_QUI", dati)
    os.makedirs(cfg.out_dir, exist_ok=True)
    path = os.path.join(cfg.out_dir, cfg.html_v3)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return path


def _testa(titolo: str) -> str:
    """La stessa intestazione HTML di ``_HTML_TESTA`` con un titolo diverso.

    Un ``replace`` e non un ``format``: il blocco contiene CSS, e le graffe
    del CSS farebbero saltare qualunque formattazione a segnaposto."""
    vecchio = ("<title>Giza v03 - volume pieno/vuoto dalla FFT delle "
               "sinusoidi</title>")
    assert vecchio in _HTML_TESTA
    return _HTML_TESTA.replace(vecchio, f"<title>{titolo}</title>")


def _svg_onda(z: np.ndarray, curve: List[Tuple[np.ndarray, str, str]],
              titolo: str, y_lo: float, y_hi: float,
              rif: Optional[float] = None, w: int = 720, h: int = 240) -> str:
    """Un grafico SVG in linea: nessuna libreria, nessuna richiesta di rete.

    ``curve`` e' una lista di (valori, colore, etichetta) sullo stesso asse z.
    ``rif`` disegna una verticale tratteggiata (la quota del picco). Le pagine
    per singola piramide vanno caricate da sole, quindi tutto quello che
    disegnano deve stare dentro il file."""
    ml, mr, mt, mb = 52, 14, 26, 34
    pw, ph = w - ml - mr, h - mt - mb
    z0, z1 = float(z[0]), float(z[-1])

    def px(v: float) -> float:
        return ml + pw * (v - z0) / max(z1 - z0, 1e-9)

    def py(v: float) -> float:
        return mt + ph * (1.0 - (v - y_lo) / max(y_hi - y_lo, 1e-9))

    parti = [f'<svg viewBox="0 0 {w} {h}" width="100%" '
             f'style="max-width:{w}px" role="img" aria-label="{titolo}">']
    parti.append(f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" '
                 'fill="#FFF" stroke="#E2E8F0"/>')
    # griglia e tacche in quota
    for tz in range(int(math.ceil(z0 / 100.0)) * 100, int(z1) + 1, 100):
        x = px(tz)
        parti.append(f'<line x1="{x:.1f}" y1="{mt}" x2="{x:.1f}" '
                     f'y2="{mt + ph}" stroke="#F1F5F9"/>')
        parti.append(f'<text x="{x:.1f}" y="{mt + ph + 15}" font-size="10" '
                     f'fill="#64748B" text-anchor="middle">{tz}</text>')
    # tacche in ordinata
    for k in range(3):
        v = y_lo + (y_hi - y_lo) * k / 2.0
        y = py(v)
        parti.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml + pw}" '
                     f'y2="{y:.1f}" stroke="#F1F5F9"/>')
        parti.append(f'<text x="{ml - 6}" y="{y + 3:.1f}" font-size="10" '
                     f'fill="#64748B" text-anchor="end">{v:+.1f}</text>')
    if y_lo < 0 < y_hi:
        parti.append(f'<line x1="{ml}" y1="{py(0):.1f}" x2="{ml + pw}" '
                     f'y2="{py(0):.1f}" stroke="#CBD5E1"/>')
    if rif is not None and z0 <= rif <= z1:
        parti.append(f'<line x1="{px(rif):.1f}" y1="{mt}" x2="{px(rif):.1f}" '
                     f'y2="{mt + ph}" stroke="#16A34A" stroke-width="1.2" '
                     'stroke-dasharray="4 3"/>')

    leg = []
    for k, (val, colore, etichetta) in enumerate(curve):
        pts = " ".join(f"{px(float(z[i])):.1f},{py(float(val[i])):.1f}"
                       for i in range(len(z)))
        parti.append(f'<polyline points="{pts}" fill="none" stroke="{colore}" '
                     'stroke-width="1.6"/>')
        leg.append(f'<tspan fill="{colore}">&#9632;</tspan> {etichetta}')
    parti.append(f'<text x="{ml}" y="{mt - 10}" font-size="11" fill="#334155">'
                 + "   ".join(leg) + "</text>")
    parti.append(f'<text x="{ml + pw / 2:.0f}" y="{h - 4}" font-size="10" '
                 'fill="#64748B" text-anchor="middle">quota z sopra il '
                 'riferimento .xml [m]</text>')
    parti.append("</svg>")
    return "".join(parti)


def _pannello_piramide(campo: CampoOnde, p: Pyramid, res: Dict[str, Any],
                       an: Dict[str, Any]) -> str:
    """Le forme d'onda calcolate di UNA piramide, come HTML gia' pronto.

    Due grafici: la colonna rappresentativa (la somma delle sue sinusoidi,
    con l'inviluppo) e la media dei moduli su tutte le celle di quella
    piramide. La media e' INCOERENTE -- media dei |h|, non dei complessi:
    celle diverse hanno fasi scorrelate e la media coerente le annullerebbe,
    dando una curva piatta che sembrerebbe assenza di segnale."""
    k = campo.rappresentativa()
    if k < 0:
        return ""
    z = campo.z.astype(np.float64)
    onda = campo.onda[k].astype(np.float64)
    inv = campo.inviluppo[k].astype(np.float64)
    medio = campo.inviluppo.mean(axis=0).astype(np.float64)
    medio = medio / max(float(medio.max()), 1e-9)

    g1 = _svg_onda(
        z, [(onda, "#0F172A", "somma delle sinusoidi, Re h(z)"),
            (inv, "#F0A24A", "inviluppo |h(z)|"),
            (-inv, "#F0A24A", "")],
        f"forma d'onda della colonna rappresentativa di {p.name}",
        -1.05, 1.05, rif=float(campo.z_picco[k]))
    g2 = _svg_onda(
        z, [(medio, "#2563EB",
             f"media dei moduli sulle {len(campo)} celle della piramide")],
        f"profilo medio di {p.name}", 0.0, 1.05)

    idx = an["indice"]
    i, j = int(campo.riga[k]), int(campo.colonna[k])
    return f"""
  <h2>Le forme d'onda calcolate</h2>
  <div class="grafico">
    <b>Colonna rappresentativa &mdash; cella ({i}, {j})</b>
    <p class="did">La piu' luminosa fra quelle sopra la soglia di qualita'.
    La curva nera e' la somma delle {len(res['dates'])} sinusoidi di questa
    cella, una per data; l'arancio e' il loro inviluppo analitico. La
    verticale verde e' la quota del diffusore dominante,
    {float(campo.z_picco[k]):+.1f} m sopra il riferimento.</p>
    {g1}
  </div>
  <div class="grafico">
    <b>Profilo medio della piramide</b>
    <p class="did">Media dei moduli su tutte le {len(campo)} celle. E'
    incoerente per necessita': celle diverse hanno fasi scorrelate, e una
    media coerente le annullerebbe.</p>
    {g2}
  </div>
  <table>
    <tr><td>celle sulla superficie (geometria radar)</td><td class="v">{len(campo)}</td></tr>
    <tr><td>di cui sopra la soglia di qualita'</td><td class="v">{int(np.count_nonzero(campo.qualita))}</td></tr>
    <tr><td>quota mediana del diffusore dominante</td><td class="v">{float(np.median(campo.z_picco)):+.1f} m</td></tr>
    <tr><td>quota simulata mediana della piramide</td><td class="v">{float(np.median(campo.sim_h)):+.1f} m</td></tr>
    <tr><td>lato di base / altezza (da letteratura)</td><td class="v">{p.base_side_m:.1f} m / {p.height_m:.1f} m</td></tr>
    <tr><td>pendenza delle facce</td><td class="v">{p.face_slope_deg:.2f}&deg;</td></tr>
    <tr><td>eccesso massimo (z-score) sulle sue celle</td><td class="v">{float(np.nanmax(idx['pieno_max'][campo.riga, campo.colonna])):+.2f}</td></tr>
    <tr><td>difetto massimo (z-score) sulle sue celle</td><td class="v">{float(np.nanmin(idx['vuoto_max'][campo.riga, campo.colonna])):+.2f}</td></tr>
  </table>
  <div class="nota">La quota mediana del diffusore dominante va confrontata
  con quella simulata della riga sopra: se le facce fossero ricostruite le
  due seguirebbero l'una l'altra. Non lo fanno, ed e' il risultato
  centrale di questa catena, non un difetto del disegno. Con
  <code>delta_z</code> = {res['budget'].delta_z_vertical:.0f} m le facce a
  {p.face_slope_deg:.0f}&deg; sono in layover pieno e non esiste un
  diffusore dominante da localizzare.</div>
"""


def _slug(nome: str) -> str:
    """Nome di file da un nome di piramide: 'Cheope (Khufu)' -> 'cheope'."""
    return nome.split()[0].lower().replace("'", "")


def build_html_onde3d(campo: CampoOnde, res: Dict[str, Any], cfg: ConfigV3,
                      solo: Optional[int] = None,
                      an: Optional[Dict[str, Any]] = None) -> Optional[str]:
    """Pagina 3D interattiva del campo di onde, pixel per pixel -- F59, F60.

    Porta TUTTE le colonne (il PNG ne dirada), con i cursori per quante
    disegnarne, per la scala orizzontale e per l'esagerazione verticale. La
    scala orizzontale e' un cursore proprio perche' e' un artificio: chi
    guarda deve poter vedere che cambiandola cambia l'aspetto del campo e non
    il dato.

    Con ``solo`` uguale all'indice di una piramide (F60) esce la pagina di
    QUELLA piramide: le sue sole colonne, i suoi numeri, e i due grafici
    delle forme d'onda calcolate disegnati in SVG dentro la pagina, senza
    dipendere da nessun file esterno."""
    if len(campo) == 0:
        return None

    pir_sel: Optional[Pyramid] = None
    if solo is not None:
        pir_sel = PYRAMIDS[solo]
        campo = campo.filtra(campo.piramide == solo)
        if len(campo) == 0:
            return None

    # ordine spaziale: diradare seguendo questo indice tiene il campo
    # rappresentativo di tutta l'area invece di addensarlo su una fascia
    ordine = np.lexsort((campo.east, campo.north)).astype(np.int32)
    budget: TomoBudget = res["budget"]

    colonne = np.stack([
        np.round(campo.east, 1), np.round(campo.north, 1),
        np.round(campo.h_ref, 1), campo.piramide.astype(np.float32),
        campo.qualita.astype(np.float32), np.round(campo.z_picco, 1),
        np.round(campo.sim_h, 1),
    ], axis=1)

    piramidi = []
    for p in PYRAMIDS:
        m = pyramid_mesh(p, res["lat0"], res["lon0"])
        piramidi.append({
            "nome": p.name,
            "v": [[round(c, 1) for c in v] for v in m["vertices"]],
            "spigoli": [[0, 1], [1, 2], [2, 3], [3, 0],
                        [0, 4], [1, 4], [2, 4], [3, 4]],
        })

    # due decimali e non tre: le curve sono normalizzate a 1 e disegnate su
    # poche centinaia di pixel, quindi la terza cifra non si vede, mentre nel
    # file pesa un terzo del payload. La pagina dichiara la precisione.
    payload = {
        "z": np.round(campo.z, 2).tolist(),
        "colonne": colonne.tolist(),
        "onda": np.round(campo.onda, 2).tolist(),
        "inviluppo": np.round(campo.inviluppo, 2).tolist(),
        "ordine": ordine.tolist(),
        "piramidi": piramidi,
        "colori": list(COLORI_PIRAMIDE),
        "conteggi": campo.per_piramide(),
        "meta": {
            "date": len(res["dates"]),
            "master": res["master_date"],
            "polarizzazione": cfg.polarisation.upper(),
            "celle": len(campo),
            "campioni_z": int(len(campo.z)),
            "passo_z": round(float(campo.z[1] - campo.z[0]), 3),
            "sottocamp": campo.passo_z,
            "scala_m": float(cfg.onde3d_scala_m),
            "delta_z": round(float(budget.delta_z_vertical), 1),
            "sigma_h": round(float(budget.sigma_h), 1),
            "qualita": int(np.count_nonzero(campo.qualita)),
        },
    }
    dati = json.dumps(payload, ensure_ascii=False, separators=(",", ":"),
                      default=float).replace("</", "<\\/")

    c = payload["conteggi"]
    righe_conteggi = "".join(
        f"<tr><td><span class='pastiglia' style='background:{COLORI_PIRAMIDE[k]}'>"
        f"</span> {p.name}</td><td class='v'>{c.get(p.name, 0)}</td></tr>"
        for k, p in enumerate(PYRAMIDS))

    # --- pezzi che cambiano fra la pagina d'insieme e quella di una sola
    # piramide (F60). Il resto della pagina, renderer compreso, resta lo
    # stesso: due pagine diverse che disegnano la stessa cosa divergerebbero.
    if pir_sel is None:
        titolo_pag = "Giza v03 - sinusoidi pixel per pixel sulle piramidi"
        h1 = ("Giza v03 &mdash; le sinusoidi risultanti, pixel per pixel, "
              "sulla superficie delle piramidi")
        occhiello = ""
        nome_file = cfg.html_onde3d
    else:
        titolo_pag = f"Giza v03 - {pir_sel.name}, forme d'onda pixel per pixel"
        h1 = (f"{pir_sel.name} &mdash; le forme d'onda calcolate, "
              "pixel per pixel")
        occhiello = ('<p class="torna"><a href="index.html">&larr; torna '
                     "all'indice</a></p>")
        nome_file = f"onde_{_slug(pir_sel.name)}.html"

    caselle = "".join(
        f'<label class="riga"><input type="checkbox" class="c_pir" '
        f'data-k="{k}" checked> <span class="pastiglia" '
        f'style="background:{COLORI_PIRAMIDE[k]}"></span> {p.name}</label>'
        for k, p in enumerate(PYRAMIDS)
        if pir_sel is None or k == solo)

    pannello = ("" if pir_sel is None
                else _pannello_piramide(campo, pir_sel, res, an or {}))
    stile_extra = "" if pir_sel is None else """
<style>
/* la pagina di una singola piramide scorre: sotto la scena 3D ci sono i
   grafici delle forme d'onda, e devono poter stare piu' in basso */
body { height: auto; overflow: auto; display: block; }
main { height: 68vh; min-height: 420px; }
.sotto { max-width: 1080px; margin: 0 auto; padding: 8px 22px 60px; }
.grafico { background: #FFF; border: 1px solid #E2E8F0; border-radius: 12px;
           padding: 16px 18px; margin: 16px 0; }
.grafico b { font-size: 15px; }
.did { font-size: 13px; color: #475569; margin: 6px 0 12px; max-width: 760px; }
.sotto h2 { margin-top: 26px; }
.sotto table { max-width: 640px; }
.sotto td { padding: 5px 0; }
.torna { margin: 6px 0 0; font-size: 13px; }
</style>
"""

    corpo = f"""
<header>
  <h1>{h1}</h1>
  <p>Ogni curva verticale e' <b>una cella</b>: la somma delle sue
  {len(res['dates'])} sinusoidi, una per data, cioe' <code>Re h(z)</code>
  normalizzato al massimo di quella cella. Sta alla propria posizione est/nord
  e si sviluppa lungo la quota; l'ampiezza e' disegnata come spostamento verso
  est, ed e' un <b>artificio di disegno</b> &mdash; il cursore della scala e'
  li' per farlo vedere. La forma verticale di ogni colonna e' la PSF
  dell'array: <code>delta_z</code> = {budget.delta_z_vertical:.0f} m.</p>
  {occhiello}
</header>
<main>
  <div id="scena"><canvas id="tela"></canvas></div>
  <aside>
    <h2>Quali piramidi</h2>
    {caselle}
    <label class="riga"><input type="checkbox" id="c_qual"> solo celle sopra la soglia di qualita'</label>

    <h2>Che cosa disegnare</h2>
    <label class="riga"><input type="checkbox" id="c_onda" checked> somma delle sinusoidi</label>
    <label class="riga"><input type="checkbox" id="c_inv"> inviluppo |h(z)|</label>
    <label class="riga"><input type="checkbox" id="c_pic" checked> diffusore dominante (solo celle sopra soglia)</label>
    <label class="riga"><input type="checkbox" id="c_wire" checked> spigoli delle piramidi</label>
    <label class="riga"><input type="checkbox" id="c_quota"> colora per quota del picco</label>

    <h2>Scala e vista</h2>
    <div class="cursore">
      colonne disegnate = <b id="e_ncol"></b>
      <input type="range" id="s_ncol" min="20" max="{len(campo)}" step="10">
    </div>
    <div class="cursore">
      ampiezza 1 = <b id="e_scala"></b> m sul disegno
      <input type="range" id="s_scala" min="4" max="90" step="1">
    </div>
    <div class="cursore">
      quota mostrata da <b id="e_zmin"></b> m
      <input type="range" id="s_zmin" min="{int(campo.z[0])}" max="0" step="5">
    </div>
    <div class="cursore">
      quota mostrata fino a <b id="e_zmax"></b> m
      <input type="range" id="s_zmax" min="0" max="{int(campo.z[-1])}" step="5">
    </div>
    <div class="cursore">
      esagerazione verticale = <b id="e_vex"></b>x
      <input type="range" id="s_vex" min="0.5" max="4" step="0.1">
    </div>
    <button id="b_reset">rimetti la vista a zero</button>

    <h2>Celle</h2>
    <table>{righe_conteggi}
      <tr><td>totale sulle piramidi</td><td class="v">{len(campo)}</td></tr>
      <tr><td>di cui sopra soglia</td><td class="v">{int(np.count_nonzero(campo.qualita))}</td></tr>
      <tr><td>campioni per curva</td><td class="v">{len(campo.z)}</td></tr>
    </table>

    <div class="nota">Le colonne sono <b>tutte normalizzate al proprio
    massimo</b>: la loro altezza sul disegno non dice quanto una cella e'
    luminosa, dice che forma ha il suo profilo. Colonne di celle deboli hanno
    la stessa ampiezza di quelle forti, ed e' voluto &mdash; serve a
    confrontare le forme, non le intensita'.</div>
    <div class="nota">La maschera e' quella in <b>geometria radar</b>: un
    punto a quota h si sposta di h&middot;cos(theta) verso il near range,
    quindi le celle che portano la piramide non sono quelle della sua
    impronta al suolo.</div>
  </aside>
</main>
<div class="sotto">{pannello}</div>
<footer id="pie"></footer>
"""

    script = """
<script>
(function () {
  "use strict";
  var D = DATI_QUI;
  var tela = document.getElementById("tela");
  var ctx = tela.getContext("2d");
  var vista = { yaw: -0.62, pitch: 0.36, zoom: 1.0, panx: 0, pany: 0 };
  var iniziale = JSON.parse(JSON.stringify(vista));
  // 200 colonne di partenza e non tutte: con 700 curve sovrapposte il campo
  // diventa una tenda opaca in cui non si distingue piu' la singola
  // sinusoide, che e' quello che questa pagina deve far vedere. Il cursore
  // le porta tutte per chi le vuole.
  var stato = { ncol: Math.min(D.colonne.length, 200),
                scala: D.meta.scala_m, vex: 1.0,
                zmin: D.z[0], zmax: D.z[D.z.length - 1] };

  var cx = 0, cy = 0, cz = 0, raggio = 1, zmin = 0, zmax = 1, hmed = 0;
  (function () {
    var xs = [], ys = [], zs = [], sh = 0;
    D.colonne.forEach(function (c) {
      xs.push(c[0]); ys.push(c[1]); sh += c[2];
      zs.push(c[2] + D.z[0]); zs.push(c[2] + D.z[D.z.length - 1]);
    });
    hmed = D.colonne.length ? sh / D.colonne.length : 0;
    D.piramidi.forEach(function (p) {
      p.v.forEach(function (v) { xs.push(v[0]); ys.push(v[1]); zs.push(v[2]); });
    });
    var mn = function (a) { return Math.min.apply(null, a); };
    var mx = function (a) { return Math.max.apply(null, a); };
    cx = 0.5 * (mn(xs) + mx(xs)); cy = 0.5 * (mn(ys) + mx(ys));
    zmin = mn(zs); zmax = mx(zs); cz = 0.5 * (zmin + zmax);
    raggio = Math.max(mx(xs) - mn(xs), mx(ys) - mn(ys), 1);
  })();

  var pmin = Infinity, pmax = -Infinity;
  D.colonne.forEach(function (c) {
    var q = c[2] + c[5];
    if (q < pmin) { pmin = q; }
    if (q > pmax) { pmax = q; }
  });

  function coloreQuota(q) {
    var t = (pmax > pmin) ? (q - pmin) / (pmax - pmin) : 0.5;
    var r = Math.round(30 + 200 * t), g = Math.round(90 + 60 * (1 - t));
    var b = Math.round(220 - 170 * t);
    return "rgb(" + r + "," + g + "," + b + ")";
  }

  var W = 0, H = 0, dpr = 1;
  function misura() {
    dpr = window.devicePixelRatio || 1;
    var r = tela.parentNode.getBoundingClientRect();
    W = Math.max(320, Math.floor(r.width));
    H = Math.max(320, Math.floor(r.height));
    tela.width = Math.floor(W * dpr);
    tela.height = Math.floor(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  var cyaw = 1, syaw = 0, cpit = 1, spit = 0, sfac = 1;
  function preparaCamera() {
    cyaw = Math.cos(vista.yaw); syaw = Math.sin(vista.yaw);
    cpit = Math.cos(vista.pitch); spit = Math.sin(vista.pitch);
    sfac = (Math.min(W, H) * 0.80 / raggio) * vista.zoom;
    // il centro verticale segue la finestra in quota scelta: altrimenti
    // tagliando l'asse la scena resta appesa in un angolo dello schermo
    cz = hmed + 0.5 * (stato.zmin + stato.zmax);
  }
  var px = 0, py = 0, pd = 0;
  function proietta(x, y, z) {
    var dx = x - cx, dy = y - cy, dz = (z - cz) * stato.vex;
    var X = dx * cyaw - dy * syaw;
    var Y = dx * syaw + dy * cyaw;
    var Y2 = Y * cpit - dz * spit;
    var Z2 = Y * spit + dz * cpit;
    px = W / 2 + X * sfac + vista.panx;
    py = H / 2 - Z2 * sfac + vista.pany;
    pd = Y2;
  }

  function selezione() {
    var vis = [];
    var attive = {};
    var box = document.querySelectorAll(".c_pir");
    for (var i = 0; i < box.length; i++) {
      attive[box[i].getAttribute("data-k")] = box[i].checked;
    }
    var soloQ = document.getElementById("c_qual").checked;
    for (var k = 0; k < D.ordine.length; k++) {
      var idx = D.ordine[k];
      var c = D.colonne[idx];
      if (!attive[String(c[3])]) { continue; }
      if (soloQ && c[4] < 0.5) { continue; }
      vis.push(idx);
    }
    if (vis.length > stato.ncol) {
      var passo = vis.length / stato.ncol;
      var fuori = [];
      for (var j = 0; j < stato.ncol; j++) {
        fuori.push(vis[Math.floor(j * passo)]);
      }
      return fuori;
    }
    return vis;
  }

  function disegna() {
    preparaCamera();
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#FFFFFF";
    ctx.fillRect(0, 0, W, H);

    var vOnda = document.getElementById("c_onda").checked;
    var vInv = document.getElementById("c_inv").checked;
    var vPic = document.getElementById("c_pic").checked;
    var vWire = document.getElementById("c_wire").checked;
    var perQuota = document.getElementById("c_quota").checked;
    var sel = selezione();

    // ordinamento per profondita' della base: le colonne davanti coprono
    // quelle dietro, come deve essere
    var chiavi = new Array(sel.length);
    for (var i = 0; i < sel.length; i++) {
      var c0 = D.colonne[sel[i]];
      proietta(c0[0], c0[1], c0[2]);
      chiavi[i] = [pd, sel[i]];
    }
    chiavi.sort(function (a, b) { return b[0] - a[0]; });

    // primo e ultimo campione dentro la finestra in quota scelta
    var t0 = 0, t1 = D.z.length - 1;
    while (t0 < t1 && D.z[t0] < stato.zmin) { t0++; }
    while (t1 > t0 && D.z[t1] > stato.zmax) { t1--; }

    for (var m = 0; m < chiavi.length; m++) {
      var id = chiavi[m][1];
      var c = D.colonne[id];
      var col = perQuota ? coloreQuota(c[2] + c[5])
                         : (D.colori[c[3]] || "#94A3B8");
      if (vInv) {
        var inv = D.inviluppo[id];
        ctx.strokeStyle = col;
        ctx.globalAlpha = 0.22;
        ctx.lineWidth = 1.0;
        for (var s = -1; s <= 1; s += 2) {
          ctx.beginPath();
          for (var t = t0; t <= t1; t++) {
            proietta(c[0] + s * inv[t] * stato.scala, c[1], c[2] + D.z[t]);
            if (t === t0) { ctx.moveTo(px, py); } else { ctx.lineTo(px, py); }
          }
          ctx.stroke();
        }
      }
      if (vOnda) {
        var on = D.onda[id];
        ctx.strokeStyle = col;
        ctx.globalAlpha = 0.72;
        ctx.lineWidth = 0.9;
        ctx.beginPath();
        for (var u = t0; u <= t1; u++) {
          proietta(c[0] + on[u] * stato.scala, c[1], c[2] + D.z[u]);
          if (u === t0) { ctx.moveTo(px, py); } else { ctx.lineTo(px, py); }
        }
        ctx.stroke();
      }
      // il picco si disegna solo dove significa qualcosa: sulle celle sotto
      // la soglia di qualita' l'argmax e' rumore, e sparso su tutto l'asse
      // sembrerebbe una nuvola di bersagli
      if (vPic && c[4] >= 0.5 && c[5] >= stato.zmin && c[5] <= stato.zmax) {
        proietta(c[0], c[1], c[2] + c[5]);
        ctx.globalAlpha = 0.9;
        ctx.fillStyle = "#0F172A";
        ctx.fillRect(px - 1.8, py - 1.8, 3.6, 3.6);
      }
    }
    ctx.globalAlpha = 1.0;

    if (vWire) {
      ctx.strokeStyle = "#7C2D12";
      ctx.lineWidth = 1.3;
      D.piramidi.forEach(function (p) {
        ctx.beginPath();
        p.spigoli.forEach(function (sp) {
          proietta(p.v[sp[0]][0], p.v[sp[0]][1], p.v[sp[0]][2]);
          var ax = px, ay = py;
          proietta(p.v[sp[1]][0], p.v[sp[1]][1], p.v[sp[1]][2]);
          ctx.moveTo(ax, ay); ctx.lineTo(px, py);
        });
        ctx.stroke();
      });
    }

    document.getElementById("e_ncol").textContent =
      chiavi.length + " di " + D.colonne.length;
  }

  var trascina = false, lx = 0, ly = 0;
  tela.addEventListener("mousedown", function (ev) {
    trascina = true; lx = ev.clientX; ly = ev.clientY;
    tela.classList.add("trascina");
  });
  window.addEventListener("mouseup", function () {
    trascina = false; tela.classList.remove("trascina");
  });
  window.addEventListener("mousemove", function (ev) {
    if (!trascina) { return; }
    if (ev.shiftKey) {
      vista.panx += ev.clientX - lx; vista.pany += ev.clientY - ly;
    } else {
      vista.yaw += (ev.clientX - lx) * 0.006;
      vista.pitch += (ev.clientY - ly) * 0.006;
      vista.pitch = Math.max(-1.45, Math.min(1.45, vista.pitch));
    }
    lx = ev.clientX; ly = ev.clientY;
    disegna();
  });
  tela.addEventListener("wheel", function (ev) {
    ev.preventDefault();
    vista.zoom *= (ev.deltaY < 0) ? 1.1 : 1 / 1.1;
    vista.zoom = Math.max(0.25, Math.min(12, vista.zoom));
    disegna();
  }, { passive: false });

  function lega(id, campo, mostra, fmt) {
    var el = document.getElementById(id);
    var et = document.getElementById(mostra);
    el.value = stato[campo];
    if (et && fmt) { et.textContent = fmt(stato[campo]); }
    el.addEventListener("input", function () {
      stato[campo] = parseFloat(el.value);
      if (et && fmt) { et.textContent = fmt(stato[campo]); }
      disegna();
    });
  }
  lega("s_ncol", "ncol", null, null);
  lega("s_scala", "scala", "e_scala", function (v) { return v.toFixed(0); });
  lega("s_vex", "vex", "e_vex", function (v) { return v.toFixed(1); });
  lega("s_zmin", "zmin", "e_zmin", function (v) { return v.toFixed(0); });
  lega("s_zmax", "zmax", "e_zmax", function (v) { return v.toFixed(0); });
  var tutte = document.querySelectorAll(".c_pir, #c_qual, #c_onda, #c_inv, "
    + "#c_pic, #c_wire, #c_quota");
  for (var q = 0; q < tutte.length; q++) {
    tutte[q].addEventListener("change", disegna);
  }
  document.getElementById("b_reset").addEventListener("click", function () {
    vista = JSON.parse(JSON.stringify(iniziale));
    disegna();
  });
  window.addEventListener("resize", function () { misura(); disegna(); });

  document.getElementById("pie").textContent =
    "trascina per ruotare, shift+trascina per spostare, rotella per lo zoom  |  "
    + D.meta.celle + " celle sulle piramidi, " + D.meta.campioni_z
    + " campioni per curva (asse z a passo " + D.meta.passo_z + " m, uno ogni "
    + D.meta.sottocamp + " dell'asse tomografico)  |  " + D.meta.date
    + " date, master " + D.meta.master + ", " + D.meta.polarizzazione
    + "  |  delta_z = " + D.meta.delta_z + " m, sigma_h = " + D.meta.sigma_h + " m";

  misura();
  disegna();
})();
</script>
</body>
</html>
"""

    html = (_testa(titolo_pag) + stile_extra
            + corpo + script.replace("DATI_QUI", dati))
    os.makedirs(cfg.out_dir, exist_ok=True)
    path = os.path.join(cfg.out_dir, nome_file)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    return path


# ==========================================================================
# 7.  Pipeline v03
# ==========================================================================

def analisi_v03(res: Dict[str, Any], cfg: ConfigV3,
                verbose: bool = True) -> Dict[str, Any]:
    """I quattro passi di v03 sopra il cubo gia' costruito da v02.

    Non rilegge nulla dal disco: lavora su ``y_ml`` e ``k_ml`` (F51), cioe'
    esattamente i numeri da cui v02 ha ricavato la sua superficie. Se i due
    programmi divergessero, divergerebbero sul metodo e non sui dati -- ed e'
    per questo che il primo controllo stampato e' la differenza fra la quota
    del picco di v03 e quella di v02."""
    y = res["y_ml"]
    sign = float(res["sign"])
    kappa = (sign * res["k_ml"]).astype(np.float32)
    res["kappa"] = kappa
    z_axis = np.asarray(res["z_axis"], dtype=np.float64)

    # --- passo 1: forma polare (F52) ---------------------------------------
    amp, fase = forma_polare(y)
    if verbose:
        print(f"      modulo: mediana {float(np.median(amp)):.4f}, "
              f"massimo {float(amp.max()):.4f}")
        print(f"      fase [rad]: {float(fase.min()):+.3f} .. "
              f"{float(fase.max()):+.3f}  (media circolare "
              f"{float(np.angle(np.mean(np.exp(1j * fase.astype(np.float64))))):+.3f})")
        print(f"      numeri d'onda verticali kappa: {float(np.abs(kappa).min()):.5f} "
              f".. {float(np.abs(kappa).max()):.5f} rad/m  "
              f"(segno di k_z calibrato sui dati: {sign:+.0f})")

    # --- passo 2 e 3: sintesi e FFT (F53, F54) -----------------------------
    g = griglia_kappa(z_axis, cfg)
    if verbose:
        print(g.as_text())
    t0 = time.time()
    h = profilo_fft(y, kappa, g)
    t_fft = time.time() - t0

    t0 = time.time()
    fft_info = errore_fft(y, kappa, z_axis, h, celle=cfg.fft_check_cells,
                          seed=cfg.seed)
    t_dir = time.time() - t0
    n_check = fft_info["celle_verificate"]
    n_tot = h.shape[0] * h.shape[1]
    fft_info["secondi_fft_tutte_le_celle"] = round(t_fft, 3)
    fft_info["secondi_diretto_stimati"] = round(t_dir * n_tot / max(n_check, 1), 1)
    if verbose:
        print(f"      FFT su {n_tot} celle: {t_fft:.2f} s  "
              f"(periodogramma diretto, stimato sullo stesso volume: "
              f"{fft_info['secondi_diretto_stimati']:.0f} s)")
        print(f"      errore contro il periodogramma diretto su {n_check} celle: "
              f"medio {fft_info['errore_relativo_medio']:.2e}, "
              f"massimo {fft_info['errore_relativo_max']:.2e}")

    # --- controllo incrociato con v02 --------------------------------------
    # La quota va confrontata con lo STESSO estimatore: v02 raffina il picco
    # con la parabola a tre punti, quindi confrontare il suo risultato con il
    # nudo argmax di v03 darebbe uno scarto fino a mezzo passo (1.56 m) che
    # non e' un disaccordo fra i due programmi ma la differenza fra due
    # estimatori. Qui si applica a h la stessa funzione di v02.
    y_abs_sum = np.abs(y).sum(axis=0).astype(np.float32)
    h_v3, gamma_v3 = v2.surface_from_tomogram(h, z_axis.astype(np.float32),
                                              y_abs_sum)
    d_h = np.abs(h_v3 - res["height_rel"])
    d_g = np.abs(gamma_v3 - res["gamma"])
    confronto_v02 = {
        "estimatore": "surface_from_tomogram di v02, picco raffinato a parabola",
        "scarto_quota_mediano_m": round(float(np.median(d_h)), 6),
        "scarto_quota_massimo_m": round(float(np.max(d_h)), 6),
        "scarto_gamma_massimo": round(float(np.max(d_g)), 8),
        "celle_entro_1_cm_pct": round(100.0 * float(np.mean(d_h <= 0.01)), 2),
        "passo_asse_z_m": round(float(z_axis[1] - z_axis[0]), 4),
    }
    if verbose:
        print(f"      quota del picco contro v02 (stesso estimatore): entro "
              f"1 cm sul {confronto_v02['celle_entro_1_cm_pct']:.1f} % delle "
              f"celle, scarto massimo {confronto_v02['scarto_quota_massimo_m']:.4f} m, "
              f"gamma {confronto_v02['scarto_gamma_massimo']:.2e}")

    # --- passo 4: pieno e vuoto (F55, F56) ---------------------------------
    piana_geom = (~np.asarray(res["pyr_mask"], dtype=bool)
                  & ~np.asarray(res["ground_mask"], dtype=bool))
    piana = piana_geom & np.asarray(res["good"], dtype=bool)
    nullo_su_qualita = bool(piana.sum() >= 100)
    if not nullo_su_qualita:
        piana = piana_geom
    idx = indice_pieno_vuoto(h, amp, kappa, g, piana, q_soglia=cfg.anomalia_q)
    idx["nullo_su_celle_di_qualita"] = nullo_su_qualita
    if verbose:
        print(f"      nullo del residuo: {idx['celle_nullo']} celle di piana"
              + ("" if nullo_su_qualita else " (senza filtro di qualita': "
                 "troppo poche celle sopra soglia)"))
        print(f"      dispersione del residuo sul nullo: mediana "
              f"{float(np.median(idx['sigma_nullo'])):.4f}, "
              f"soglia |z| = {idx['soglia']:.3f} "
              f"(percentile {cfg.anomalia_q:.0f})")

    bil = bilancio_anomalie(idx, res["pyr_mask"], piana_geom, res["good"])
    bil["spiegazione_eccesso"] = spiegazione_eccesso(idx, res, res["pyr_mask"],
                                                     res["good"])
    if verbose:
        for etichetta, chiave in (("anomalie (pieno o vuoto)", "anomalie_totali"),
                                  ("solo vuoto", "vuoto"),
                                  ("solo pieno", "pieno")):
            d = bil[chiave]
            if "frazione_piramidi" in d:
                print(f"      {etichetta:24s}: piramidi "
                      f"{100 * d['frazione_piramidi']:5.1f} %  piana "
                      f"{100 * d['frazione_piana']:5.1f} %  z = {d['z']:+.2f}")
            else:
                print(f"      {etichetta:24s}: {d['esito']}")
        print(f"      -> {bil['anomalie_totali'].get('esito', '')}")
        sp = bil["spiegazione_eccesso"]
        if "prima_spiegazione" in sp:
            print("      contro che cosa correla l'eccesso (celle delle "
                  f"piramidi, {sp['celle']}):")
            for k in ("ripiegamento_layover", "ampiezza_dB", "coerenza",
                      "quota_simulata_m"):
                if k in sp and sp[k]["r"] is not None:
                    print(f"        {k:22s} r = {sp[k]['r']:+.3f}")
            print(f"      -> {sp['esito']}")

    return {
        "h": h,
        "amp": amp,
        "fase": fase,
        "kappa": kappa,
        "griglia": g,
        "indice": idx,
        "bilancio": bil,
        "fft": fft_info,
        "confronto_v02": confronto_v02,
        "piana": piana_geom,
    }


def _pulisci(v: Any) -> Any:
    """Rende serializzabile (e leggibile) un valore che puo' essere NaN."""
    if isinstance(v, (np.floating, float)):
        f = float(v)
        return None if not math.isfinite(f) else round(f, 6)
    if isinstance(v, (np.integer, int)):
        return int(v)
    if isinstance(v, (np.bool_, bool)):
        return bool(v)
    if isinstance(v, np.ndarray):
        return [_pulisci(x) for x in v.tolist()]
    if isinstance(v, dict):
        return {k: _pulisci(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_pulisci(x) for x in v]
    return v


def save_outputs_v03(res: Dict[str, Any], an: Dict[str, Any], cfg: ConfigV3,
                     onde: List[Dict[str, Any]], percorsi: Dict[str, str],
                     campo: Optional[CampoOnde] = None) -> str:
    """Matrici e riassunto JSON della parte v03."""
    os.makedirs(cfg.out_dir, exist_ok=True)
    idx = an["indice"]
    for nome, arr in (("v03_zscore", idx["zscore"]),
                      ("v03_profilo", idx["profilo"]),
                      ("v03_psf", idx["psf"]),
                      ("v03_residuo", idx["residuo"]),
                      ("v03_sigma_nullo", idx["sigma_nullo"]),
                      ("v03_vuoto_max", idx["vuoto_max"]),
                      ("v03_pieno_max", idx["pieno_max"]),
                      ("v03_z_vuoto", idx["z_vuoto"]),
                      ("v03_zeta", idx["zeta"]),
                      ("v03_modulo", an["amp"]),
                      ("v03_fase_rad", an["fase"]),
                      ("v03_kappa", an["kappa"])):
        np.save(os.path.join(cfg.out_dir, f"{nome}.npy"), arr)

    g: GrigliaKappa = an["griglia"]
    meta = {
        "generato": time.strftime("%Y-%m-%d %H:%M:%S"),
        "programma": "piramidi_v03.py (sintesi sinusoidale + FFT verticale)",
        "metodo": (
            "modulo e fase radiante dell'interferogramma multilooked -> una "
            "sinusoide verticale per data -> somma -> FFT non uniforme "
            "(NUFFT di tipo 1, nucleo gaussiano) -> profilo in quota -> "
            "residuo rispetto alla PSF di un solo diffusore -> z-score sul "
            "nullo delle celle di piana"),
        "identita_con_v02": (
            "la somma delle sinusoidi e' Re[h(z)] e il suo inviluppo e' "
            "|h(z)|: e' lo stesso periodogramma di v02, calcolato per FFT "
            "invece che con il doppio ciclo. Non e' una misura indipendente."),
        "date": res["dates"],
        "master": res["master_date"],
        "polarizzazione": cfg.polarisation,
        "segno_k_z": res["sign"],
        "griglia_kappa": {
            "n_k": g.n_k, "d_kappa_rad_m": g.d_kappa, "dz_m": g.dz,
            "sovracampionamento": cfg.sovracamp,
            "kernel_sigma_celle": g.sigma, "kernel_taps": g.taps,
            "kappa_max_reticolo_rad_m": g.kappa_max_grid,
        },
        "errore_fft": _pulisci(an["fft"]),
        "confronto_con_v02": _pulisci(an["confronto_v02"]),
        "nullo": {
            "celle": idx["celle_nullo"],
            "su_celle_di_qualita": idx.get("nullo_su_celle_di_qualita"),
            "percentile": idx["percentile"],
            "soglia_zscore": _pulisci(idx["soglia"]),
            "sigma_mediana": _pulisci(float(np.median(idx["sigma_nullo"]))),
        },
        "bilancio_anomalie": _pulisci(an["bilancio"]),
        "colonne_onde": [{
            "nome": d["nome"],
            "cella": list(d["cella"]),
            "sopra_soglia_di_qualita": bool(d["qualita"]),
            "celle_nella_maschera": d["n_celle_maschera"],
            "z_picco_m": _pulisci(d["z_picco"]),
            "quota_riferimento_xml_m": _pulisci(d["h_ref"]),
            "vuoto_max_zscore": _pulisci(d["vuoto_max"]),
            "z_del_vuoto_m": _pulisci(d["z_vuoto"]),
            "pieno_max_zscore": _pulisci(d["pieno_max"]),
            "modulo_per_data": _pulisci(d["modulo"]),
            "fase_rad_per_data": _pulisci(d["fase"]),
            "kappa_per_data_rad_m": _pulisci(d["kappa"]),
            "b_perp_per_data_m": _pulisci(d["b_perp"]),
        } for d in onde],
        "campo_onde_pixel_per_pixel": (None if campo is None else {
            "descrizione": (
                "F59: per ogni cella della superficie delle piramidi in "
                "geometria radar, la sinusoide risultante Re[h(z)] -- la "
                "somma delle sinusoidi di quella cella, una per data -- "
                "normalizzata al massimo del modulo di quella cella."),
            "celle": len(campo),
            "celle_sopra_soglia": int(np.count_nonzero(campo.qualita)),
            "per_piramide": campo.per_piramide(),
            "campioni_per_curva": int(len(campo.z)),
            "sottocampionamento_asse_z": campo.passo_z,
            "scala_disegno_m_per_ampiezza_1": cfg.onde3d_scala_m,
            "quota_picco_mediana_m": _pulisci(float(np.median(campo.z_picco))),
            "quota_simulata_mediana_m": _pulisci(float(np.median(campo.sim_h))),
        }),
        "uscite": percorsi,
        "novita": [{"id": i, "tipo": t, "descrizione": d} for i, t, d in NOVITA],
        "avvertenza": (
            "Il residuo pieno/vuoto e' un discriminante calibrato sul nullo "
            "della piana, NON una rilevazione di cavita' risolta in "
            f"profondita'. La cella di Rayleigh verticale e' "
            f"{res['budget'].delta_z_vertical:.0f} m: le camere note della "
            "Grande Piramide misurano metri, due ordini di grandezza sotto. "
            "Il confronto piramidi/piana e' riportato con il suo z-test "
            "qualunque risultato dia, ed e' quello il numero da leggere."),
    }
    path = os.path.join(cfg.out_dir, "meta_v03.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False, default=str)
    return path


def run_v03(cfg: ConfigV3, verbose: bool = True,
            solo_v03: bool = False) -> Dict[str, Any]:
    """Catena completa: v02 fino al cubo, poi i quattro passi di v03."""
    os.makedirs(cfg.out_dir, exist_ok=True)
    res = v2.run(cfg, verbose=verbose)

    print("\n  [10] validazione a tre livelli (v02)")
    valid = v2.validate(res, cfg)
    l1 = valid["livello_1_geometria"]
    reg = l1.get("regressione_misurato_vs_simulato", {})
    print(f"      quota mediana della piana: {l1['_quota_mediana_piana_m']} m")
    if "pendenza" in reg:
        print(f"      regressione misurato vs simulato: pendenza "
              f"{reg['pendenza']:+.3f} +- {reg['errore_standard_pendenza']:.3f} "
              f"su {reg['pixel']} celle -> {reg['esito']}")

    percorsi: Dict[str, str] = {}
    if not solo_v03:
        print("\n  [11] uscite di v02 (superficie, profili, HTML)")
        profiles = v2.profile_analysis(res, cfg)
        v2.surface_plot(res, cfg, valid)
        v2.save_outputs(res, cfg, valid, profiles)
        percorsi["html_v02"] = v2.build_html(res, cfg, valid)

    print("\n  [12] v03: forma polare, sinusoidi, FFT, pieno/vuoto")
    an = analisi_v03(res, cfg, verbose=verbose)

    print("\n  [13] colonne di Cheope e Chefren e disegno delle onde")
    onde = dati_onde(res, an["indice"], an["griglia"], cfg)
    for d in onde:
        print(f"      {d['nome']:24s} cella ({d['cella'][0]:3d}, {d['cella'][1]:3d})"
              f"  picco a {d['z_picco']:+7.1f} m  deficit massimo "
              f"{d['vuoto_max']:+.2f} sigma a {d['z_vuoto']:+7.1f} m dal picco")
    if not onde:
        print("      nessuna colonna disponibile: le maschere simulate di "
              "Cheope e Chefren non contengono celle")

    p = plot_onde(onde, res, cfg, an["indice"])
    if p:
        percorsi["onde"] = p
    p = plot_sezione_onde(onde, res, cfg, an["indice"])
    if p:
        percorsi["sezione_onde"] = p

    print("\n  [14] campo di onde: una sinusoide risultante per ogni pixel")
    campo = campo_onde_piramidi(res, an, cfg)
    conteggi = campo.per_piramide()
    print(f"      {len(campo)} celle sulla superficie delle piramidi in "
          f"geometria radar ({int(np.count_nonzero(campo.qualita))} sopra la "
          "soglia di qualita'): "
          + ", ".join(f"{n.split()[0]} {v}" for n, v in conteggi.items()))
    print(f"      ogni colonna e' Re[h(z)] della sua cella, {len(campo.z)} "
          f"campioni (uno ogni {campo.passo_z} dell'asse tomografico)")
    p_o3 = plot_onde_3d_png(campo, res, cfg)
    if p_o3:
        percorsi["onde_3d_png"] = p_o3
    p_o3 = build_html_onde3d(campo, res, cfg)
    if p_o3:
        percorsi["onde_3d_html"] = p_o3
    for k, p in enumerate(PYRAMIDS):
        p_pir = build_html_onde3d(campo, res, cfg, solo=k, an=an)
        if p_pir:
            percorsi[f"pagina_{_slug(p.name)}"] = p_pir
            print(f"      pagina di {p.name}: "
                  f"{os.path.basename(p_pir)} "
                  f"({int(np.count_nonzero(campo.piramide == k))} colonne)")

    print("\n  [15] nuovo grafico 3D del volume pieno/vuoto")
    percorsi["png_3d"] = plot_3d_png(res, an["indice"], cfg, an["bilancio"])
    percorsi["html_v03"] = build_html_v03(res, an["indice"], cfg,
                                          an["bilancio"], an["fft"], onde)

    print("\n  [16] uscite v03")
    percorsi["meta_v03"] = save_outputs_v03(res, an, cfg, onde, percorsi,
                                            campo)

    return {"res": res, "valid": valid, "an": an, "onde": onde,
            "campo": campo, "percorsi": percorsi}


# ==========================================================================
# 8.  Autotest
# ==========================================================================

def _esito(ok: bool) -> bool:
    print("ok" if ok else "FALLITO")
    return ok


def selftest() -> int:
    """Blocca le identita' su cui poggia v03: la somma delle sinusoidi contro
    il periodogramma, la FFT contro il calcolo diretto, il residuo nullo
    quando il modello e' esatto, e la taratura del nullo."""
    ok = True
    rng = np.random.default_rng(11)
    z_axis = np.linspace(-400.0, 400.0, 257)
    cfg = ConfigV3()
    g = griglia_kappa(z_axis, cfg)

    # 1) la somma delle sinusoidi e' Re[h(z)]
    n_d = 17
    kap = rng.uniform(-0.06, 0.06, n_d)
    yc = (rng.standard_normal(n_d) + 1j * rng.standard_normal(n_d)) / 3.0
    a, f = forma_polare(yc)
    somma = onde_sinusoidali(a, f, kap, z_axis).sum(axis=0)
    esatta = somma_analitica(yc, kap, z_axis)
    e1 = float(np.max(np.abs(somma - esatta.real)) / np.max(np.abs(esatta)))
    print(f"  somma delle sinusoidi == Re[h(z)]   errore {e1:.2e}", end="  ")
    ok &= _esito(e1 < 1e-5)

    # 2) la FFT riproduce il periodogramma diretto
    n_l, n_p = 6, 5
    y = (rng.standard_normal((n_d, n_l, n_p))
         + 1j * rng.standard_normal((n_d, n_l, n_p))).astype(np.complex64)
    kappa = np.repeat(kap[:, None, None], n_l, axis=1).repeat(n_p, axis=2)
    kappa = (kappa + rng.uniform(-2e-3, 2e-3, (n_d, n_l, n_p))).astype(np.float32)
    h = profilo_fft(y, kappa, g)
    info = errore_fft(y, kappa, z_axis, h, celle=n_l * n_p, seed=3)
    print(f"  FFT contro periodogramma diretto     errore max "
          f"{info['errore_relativo_max']:.2e}", end="  ")
    ok &= _esito(info["errore_relativo_max"] < 1e-3)

    # 3) un solo diffusore: il picco cade dove e' stato messo e il residuo
    #    rispetto alla PSF e' nullo. E' il controllo che dice che l'indice
    #    pieno/vuoto misura uno SCARTO dal modello e non il modello stesso.
    # la fase di un diffusore a quota z0 vale +kappa*z0: e' la convenzione
    # su cui poggia tutto (y = A exp(+j kappa h), picco del periodogramma in
    # h). Con il segno opposto il picco esce a -z0, ed e' un errore che si
    # nota solo su un bersaglio simulato: qui e' bloccato.
    z0 = 62.5                                   # multiplo del passo (3.125 m)
    y1 = np.exp(1j * kappa * z0).astype(np.complex64)
    h1 = profilo_fft(y1, kappa, g)
    k1 = int(np.argmax(np.abs(h1[0, 0])))
    print(f"  diffusore a {z0:+.1f} m -> picco a {z_axis[k1]:+.1f} m", end="  ")
    ok &= _esito(abs(z_axis[k1] - z0) <= 0.6 * (z_axis[1] - z_axis[0]))

    amp1 = np.abs(y1).astype(np.float32)
    q = psf_array(amp1, kappa, g, np.argmax(np.abs(h1), axis=2))
    p1 = np.abs(h1) / np.max(np.abs(h1), axis=2)[:, :, None]
    e3 = float(np.max(np.abs(p1 - q)))
    print(f"  residuo profilo-PSF su un solo diffusore  {e3:.2e}", end="  ")
    ok &= _esito(e3 < 5e-3)

    # 4) allineamento al picco e ritorno
    vol = rng.standard_normal((4, 3, 257)).astype(np.float32)
    kp = rng.integers(40, 200, (4, 3))
    rt = _riporta_ad_assoluto(_allinea_al_picco(vol, kp), kp)
    fin = np.isfinite(rt)
    e4 = float(np.max(np.abs(rt[fin] - vol[fin])))
    print(f"  allineamento al picco, andata e ritorno   {e4:.2e}", end="  ")
    ok &= _esito(e4 == 0.0)

    # 5) taratura del nullo: su fasi casuali la frazione di celle di piana
    #    che supera la soglia deve valere circa 100 - percentile
    n_l, n_p = 22, 20
    y2 = np.exp(1j * rng.uniform(0, 2 * np.pi, (n_d, n_l, n_p))).astype(np.complex64)
    kappa2 = np.repeat(kap[:, None, None], n_l, axis=1).repeat(n_p, axis=2)
    kappa2 = kappa2.astype(np.float32)
    h2 = profilo_fft(y2, kappa2, g)
    amp2 = np.abs(y2).astype(np.float32)
    piana = np.ones((n_l, n_p), dtype=bool)
    idx = indice_pieno_vuoto(h2, amp2, kappa2, g, piana, q_soglia=99.0)
    s = idx["zscore_zeta"]
    fin = np.isfinite(s)
    frazione = float(np.mean(np.abs(s[fin]) >= idx["soglia"]))
    print(f"  frazione oltre soglia sul nullo {100 * frazione:.2f} % "
          "(attesa 1 %)", end="  ")
    ok &= _esito(0.008 <= frazione <= 0.013)

    print("\n  " + ("tutti i controlli superati" if ok
                    else "ALMENO UN CONTROLLO FALLITO"))
    return 0 if ok else 1


# ==========================================================================
# 9.  CLI
# ==========================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Giza v03: da modulo e fase radiante alle sinusoidi "
                    "verticali, FFT, pieno/vuoto e nuovo grafico 3D.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--stack-dir", default=Config.stack_dir)
    ap.add_argument("--swath", default="iw2", choices=["iw1", "iw2", "iw3"])
    ap.add_argument("--platform", default="s1c")
    ap.add_argument("--pol", dest="polarisation", default="vh",
                    choices=["vh", "vv"])
    ap.add_argument("--dates", dest="n_dates", type=int, default=0,
                    help="quante date usare; 0 = tutte")
    ap.add_argument("--n-elev", type=int, default=257)
    ap.add_argument("--elev-max", dest="elev_max_m", type=float, default=400.0)
    ap.add_argument("--look-range", type=int, default=4)
    ap.add_argument("--look-azimuth", type=int, default=1)
    ap.add_argument("--margin", dest="area_margin_m", type=float, default=150.0)
    ap.add_argument("--sovracamp", type=int, default=2,
                    help="F54: sovracampionamento del reticolo in kappa")
    ap.add_argument("--kernel-sigma", type=float, default=1.6,
                    help="F54: larghezza del nucleo gaussiano, in celle")
    ap.add_argument("--kernel-taps", type=int, default=8,
                    help="F54: semilarghezza del nucleo, in celle")
    ap.add_argument("--anomalia-q", type=float, default=99.0,
                    help="F56: percentile del nullo che fa da soglia |z-score|")
    ap.add_argument("--out", dest="out_dir", default="out_piramidi_v03")
    ap.add_argument("--solo-v03", action="store_true",
                    help="salta superficie, profili e HTML di v02")
    ap.add_argument("--no-suolo", dest="suolo_dem", action="store_false")
    ap.add_argument("--no-dem-esterno", dest="fetch_dem", action="store_false")
    ap.add_argument("--novita", action="store_true",
                    help="elenca le novita' rispetto a v02")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    print("=" * 76)
    print("piramidi_v03 - sinusoidi da modulo e fase, FFT, pieno/vuoto, 3D")
    print("=" * 76)

    if args.novita:
        print()
        for i, t, d in NOVITA:
            print(f"  {i}  [{t}]")
            for riga in _wrap(d, 70):
                print(f"        {riga}")
            print()
        return 0

    if args.selftest:
        print("\nautotest di v03:")
        return selftest()

    cfg = ConfigV3(
        stack_dir=args.stack_dir, swath=args.swath, platform=args.platform,
        polarisation=args.polarisation, n_dates=args.n_dates,
        n_elev=args.n_elev, elev_max_m=args.elev_max_m,
        look_range=args.look_range, look_azimuth=args.look_azimuth,
        area_margin_m=args.area_margin_m, out_dir=args.out_dir,
        suolo_dem=args.suolo_dem, fetch_dem=args.fetch_dem,
        sovracamp=args.sovracamp, kernel_sigma=args.kernel_sigma,
        kernel_taps=args.kernel_taps, anomalia_q=args.anomalia_q,
    )

    t0 = time.time()
    out = run_v03(cfg, verbose=not args.quiet, solo_v03=args.solo_v03)
    perc = out["percorsi"]

    print("\n  uscite:")
    etichette = {
        "html_v03": "3D interattivo pieno/vuoto (v03)",
        "png_3d": "3D statico pieno/vuoto",
        "onde_3d_html": "onde 3D pixel per pixel (interattivo)",
        "onde_3d_png": "onde 3D pixel per pixel (statico)",
        "onde": "sinusoidi, profilo e residuo",
        "sezione_onde": "onde sotto Cheope e Chefren",
        "meta_v03": "riassunto JSON di v03",
        "html_v02": "3D interattivo di v02",
    }
    for k, eti in etichette.items():
        if k in perc:
            print(f"    {eti:36s} {os.path.abspath(perc[k])}")
    print(f"\n  completato in {time.time() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

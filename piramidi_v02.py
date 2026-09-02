#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramidi_v02.py  --  revisione critica del 2026-08-28
=====================================================

Tomografia multi-baseline della piana di Giza da stack Sentinel-1 IW, canale
**VH**, con ricostruzione della **superficie reale** e uscita 3D interattiva.

Dati usati
----------
* geometria, orbite, geocodifica e quota di riferimento  -> file ``*.annotation.xml``
* misura (ampiezza e fase, pixel per pixel)              -> file ``*.tiff``

Che cosa e' cambiato
--------------------
Le correzioni sono elencate in ``FIXES`` e stampate da ``--fixes``. Le
sette piu' pesanti:

* **F42** -- ogni data veniva LETTA alla stessa finestra assoluta del master e
  riallineata dopo con ``apply_shift()``, che trasla in modo CIRCOLARE. Ma
  l'inquadramento del prodotto cambia da una missione all'altra: sulla stessa
  traccia relativa 58, il bersaglio di Giza cade 4510 linee piu' indietro nei
  prodotti S1A che nei S1C. Con un ritaglio di 89 linee, quello "slave" era
  deserto a 63 km di distanza, ripiegato su se stesso e dirampato con i
  parametri del burst sbagliato: 16 date su 28 erano rumore, e la pila
  multi-missione -- proprio quella che porta l'escursione di baseline da 163 a
  273 m -- non era utilizzabile. Ora l'offset intero entra nella finestra di
  lettura e nell'indice di burst; alla rampa di fase resta il residuo
  sub-pixel.
* **F01/F02** -- il cubo conteneva SLC grezze e la fase di terra piatta veniva
  tolta con una regressione su ``np.unwrap`` applicato a uno scatter mascherato
  e appiattito. Su dati decorrelati quell'unwrap e' un random walk: iniettava
  uno schermo di fase casuale, e l'inversione tomografica non poteva che
  restituire rumore. Ora il cubo e' interferometrico e la fase di riferimento
  e' calcolata analiticamente da orbite e geolocation grid.
* **F15** -- risoluzione e precisione erano confuse. ``delta_z`` (243 m) e' la
  separazione di Rayleigh fra DUE diffusori nella stessa cella; ``sigma_h``
  (15 m) e' l'errore sulla quota di UNO. E' la seconda che rende disegnabile
  una superficie, e la versione precedente non la calcolava.
* **F22** -- la soglia di qualita' era fissata a mano a 0.35, sotto la mediana
  della distribuzione NULLA del periodogramma (0.56 con 11 baseline). Accettava
  il 100 per cento dei pixel, rumore compreso.
* **F28/F29** -- si processava un rettangolo di 1.8 x 2.2 km su tutta la piana
  (50 000 celle) e la vista 3D ridisegnava 8 000 poligoni ordinati a ogni
  fotogramma, con la stringa di colore ricostruita ogni volta; con
  prefers-reduced-motion il ciclo di disegno non partiva nemmeno e nessun
  comando aveva effetto. Ora si elabora la sola area delle piramidi in
  geometria radar (1.6 x 1.3 km, 8 400 celle, 14 s invece di 54) e il disegno
  e' a richiesta con colori precalcolati e ordinamento O(n): 1.4 ms per
  fotogramma.
* **F30** -- la superficie misurata era disegnata come maglia di poligoni: ogni
  faccia univa quattro nodi e riempiva il vuoto fra di essi, cosi' un campo in
  cui solo il 23 per cento delle celle supera la soglia usciva sullo schermo
  come una superficie continua. Ora e' una nuvola di punti, un punto per nodo
  misurato, e i nodi sotto soglia restano punti piccoli e spenti. Le uniche
  maglie rimaste sono i riferimenti dichiarati (piramidi ideali, superficie
  simulata, superficie .xml e strutture note), che non sono misure.
* **F31** -- il canale di coerenza del micro-moto era |media(dev)|/media(|dev|)
  su un vettore appena demediato: zero per costruzione, su ogni pixel, sempre.
  Ora e' la CONCENTRAZIONE spettrale della traccia (quota di energia nella riga
  dominante, 1 per un tono puro e 0.091 per rumore bianco con N_D = 12), e con
  essa viene riportata la frequenza meccanica di quella riga.

Che cosa dicono i dati (VH, 11 date, master 20260310)
------------------------------------------------------
1. La catena funziona: il 23 per cento delle celle supera la soglia calibrata
   sul nullo (contro l'1 per cento atteso per caso), e la piana esce piatta a
   ~63.5 m contro i 63.1-64.2 m della quota di riferimento dello .xml. Quella
   quota di riferimento e' la bilineare LOCALE (F36) fra i quattro nodi VERI
   della geolocation grid dello .xml che racchiudono il ritaglio: la grid ha
   solo 231 nodi su tutta la scena (spaziatura ~1500 x 1300 pixel) contro i
   ~90 x 370 pixel del ritaglio, quindi dentro il ritaglio la bilineare esce
   quasi piatta per limite di risoluzione del riferimento, non perche' il
   suolo vero (o il calcolo) lo sia -- sull'intera scena quella stessa grid
   porta rilievo fino a centinaia di metri (F35). Una spline bicubica globale
   sugli stessi 231 nodi (usata prima di F36) dava un valore diverso di 6 m
   proprio al centro di questo ritaglio, per curvatura importata da nodi a
   decine di km di distanza: la bilineare locale usa solo i quattro nodi
   reali della cella che contiene il ritaglio.
2. **Le piramidi non vengono ricostruite.** La regressione fra quota misurata e
   quota simulata in geometria radar da' pendenza Theil-Sen -0.009 con IC95
   [-0.058, +0.039] contro l'1.0 atteso, su 168 celle. Le mediane per fascia
   restano entro pochi metri da zero per quote simulate da 5 a 110 m: nessuna
   dipendenza dalla quota simulata, non un profilo.
3. La causa e' geometrica, non di processing. Le facce a ~52 gradi superano
   l'incidenza di 37 gradi: sono in layover pieno, con fino a 705 punti di
   superficie ripiegati nella stessa cella. Con ``delta_z`` = 243 m non esiste
   un diffusore dominante da localizzare, e il periodogramma restituisce il
   centro di fase della miscela, che sta appena sopra il deserto.
4. Controprova cambiando canale (VV): stessa conclusione, pendenza +0.038.
5. **Il canale VH e' la scelta peggiore per questi bersagli** (F26): il picco su
   Cheope vale +10.5 dB in VH contro +22.2 dB in VV. Il ritorno delle facce a
   gradoni e' co-polarizzato. VH resta il default perche' e' quanto richiesto,
   ma il costo e' misurato e dichiarato, non giustificato a parole.

Coerenza con le fonti (Biondi & Malanga 2022, WO 2024/008365, arXiv 2206.09200)
------------------------------------------------------------------------------
Quello che segue le fonti alla lettera:

* banda di guardia ``B_DL = B_cD/2`` sempre sottratta, master focalizzato su
  ``B_sub = B_cD - B_DL`` (ch03, ch11)  ->  ``Config.guard_fraction = 0.5``;
* ``B_shift`` come selettore della frequenza meccanica e ``N_D`` come frequenza
  di campionamento della vibrazione, non come parametro di calcolo (ch11)  ->
  ``--b-shift`` e ``--nd``, con il vincolo ``B_shift < B_DL/3`` verificato e
  stampato in ``MMPlan.as_text``;
* schema a 11 blocchi (ch13): FFT2 diretta calcolata UNA volta fuori dal ciclo,
  dentro il ciclo solo i due filtri, le due IFFT2 e il tracker  ->
  ``micro_motion_energy``;
* oscillatore linearizzato a 2 gradi di liberta' ``r(t) = (a cos, b sin)
  e^{-lt/2}`` (ch12): {a, b} sono esattamente gli shift che il coregistratore
  misura, quindi l'energia e' la deviazione standard del vettore COMPLESSO
  demediato, non di ``|shift|`` (F09);
* ``delta_z = lambda*R/(2*A)`` con lambda ACUSTICA, mai quella radar (ch04,
  ch12, ch14)  ->  ``MMPlan.delta_z_acoustic``, con v = 6000 m/s come nel
  paper di Giza;
* ``K_z = 4*pi*B_perp/(lambda*r_i*sin(theta))`` per pixel e per data (F04);
* protocollo di validazione a tre livelli (ch15) e distinzione fra risultato
  misurato e interpretazione (ch08-ch09).

Le divergenze, tutte deliberate e dichiarate:

1. **La profondita' qui viene dalle baseline orbitali, non dal micro-moto.**
   Le fonti lavorano su spotlight COSMO-SkyMed con ``B_cD`` ~ 22 kHz e indagano
   a 12.5 kHz; il TOPS di Sentinel-1 da' ``B_cD`` = 313 Hz, 143 ms di
   illuminazione e una apertura sintetica di 959 m, quindi ``delta_z``
   acustica = 68 km (contro i 36 m del Vesuvio e gli 0.92 m di Giza). Il banco
   di sub-aperture resta, ma produce un ATTRIBUTO di superficie: i blocchi
   9-11 delle fonti (focalizzazione FFT2 -> mappa tomografica -> filtro e
   geocodifica) non sono applicati alla traccia di vibrazione, e il programma
   non lo nasconde.
2. **Passo di marcia.** Le fonti danno ``passo = (B_cD - B_DL)/N_D``; con quel
   passo il master arriva esattamente al bordo dello spettro e lo slave, tenuto
   a ``+B_shift``, lo supera. Qui il passo e' ``(B_cD - B_sub - B_shift)/(N_D-1)``,
   che e' la stessa marcia con lo slave dentro banda (F06, verificato in
   ``--selftest``).
3. **Livelli 2 e 3 della validazione sono interni.** Le fonti confrontano la
   misura grezza con un sismografo co-locato e la struttura con una modalita'
   diversa; qui non esistono ne' l'uno ne' l'altra, e i due livelli sono
   controlli di autoconsistenza (dispersione sulla piana contro ``sigma_h``
   teorica; separabilita' piramide/piana e bilancio di layover). Vanno pesati
   come tali. La controprova disponibile e' il cambio di parametro di ch14,
   qui fatto sulla polarizzazione (VH contro VV).

Avvertenza epistemica (ch08-ch09 delle fonti): la superficie ricostruita e la
separabilita' delle piramidi dal deserto sono misure, con il loro errore.
L'indice di solidita' e' un discriminante multi-attributo dichiarato come tale,
NON una rilevazione di cavita' risolta in profondita'. Le camere note misurano
metri: due ordini di grandezza sotto la risoluzione di questi dati. E il
brevetto WO 2024/008365 e' una domanda pubblicata con rapporto di ricerca di
categoria X su tutte e 10 le rivendicazioni: va citato come divulgazione di un
metodo, mai come brevetto concesso (ch10).
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import time
import urllib.request
from dataclasses import dataclass, asdict, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import xml.etree.ElementTree as ET

from piramidi_v01 import (
    C_LIGHT,
    DEFAULT_STACK,
    Chip,
    Config as ConfigV1,
    S1Annotation,
    StackEntry,
    _azimuth_band_mask,
    _batch_subpixel_shift,
    _block2_fft2,
    _block34_bandpass,
    _block56_ifft2,
    _txt,
    parse_annotation,
    tops_deramp,
)

# --------------------------------------------------------------------------
# Registro delle correzioni
# --------------------------------------------------------------------------

FIXES: Tuple[Tuple[str, str, str], ...] = (
    ("F01", "critico",
     "Il cubo tomografico ora contiene l'interferogramma y_i = s_i * conj(s_m) "
     "normalizzato, non le SLC grezze. Con le SLC grezze la fase di ogni data e' "
     "dominata dal termine di propagazione exp(-j4piR_i/lambda), che varia di "
     "migliaia di cicli fra date e sommerge completamente il termine k_z*z su cui "
     "si regge l'inversione."),
    ("F02", "critico",
     "La fase di terra piatta e di topografia di riferimento e' rimossa "
     "analiticamente da orbita + geolocation grid (delta_R geometrico), non con "
     "una regressione su np.unwrap di uno scatter mascherato e appiattito, che su "
     "dati decorrelati e' un random walk e iniettava uno schermo di fase casuale."),
    ("F03", "alto",
     "Interpolazione orbitale di Lagrange all'ordine 8 al posto di np.interp "
     "lineare. I vettori di stato distano 10 s: l'interpolazione lineare sbaglia "
     "di metri, e l'escursione di baseline utile e' di appena ~163 m."),
    ("F04", "alto",
     "k_z calcolato con lo slant range e l'angolo di incidenza PROPRI di ciascuna "
     "acquisizione e di ciascun pixel (K_z = 4*pi*B_perp/(lambda*r_i*sin(theta))), "
     "non con un R medio scalare."),
    ("F05", "alto",
     "Il segno di k_z non e' piu' assunto: viene calibrato sui dati verificando "
     "che le piramidi risultino SOPRA la piana e non sotto, e il margine della "
     "decisione viene riportato (validazione di livello 1, ch15)."),
    ("F06", "alto",
     "Le sub-bande Doppler master/slave restano dentro lo spettro. Prima l'ultima "
     "sub-apertura slave sconfinava a 0.583*B_cD contro un limite di 0.5*B_cD, "
     "quindi le ultime marce del banco leggevano spettro vuoto."),
    ("F07", "alto",
     "Il micro-moto si calcola su un chip esteso in azimuth (1024 linee dentro il "
     "burst) invece che sulle ~185 linee del chip tomografico. Con 185 linee la "
     "griglia Doppler ha 7.8 Hz per bin su 313 Hz totali: le 12 sub-aperture "
     "distavano meno di due bin l'una dall'altra e la decomposizione era priva di "
     "significato."),
    ("F08", "medio",
     "Finestratura di Hann sulle sub-bande al posto del rect, che generava lobi "
     "laterali sinc in azimuth e sporcava il pixel tracking."),
    ("F09", "alto",
     "L'energia di micro-moto e' la deviazione standard del vettore di shift "
     "COMPLESSO demediato, non std(|shift|). Con std(|shift|) una vibrazione a "
     "modulo costante -- il caso dell'oscillatore a 2 gradi di liberta' delle "
     "fonti, eq. (20) -- restituiva energia nulla."),
    ("F10", "medio",
     "La griglia del micro-moto e' riportata al reticolo pieno con "
     "map_coordinates sulle coordinate reali dei nodi, non con zoom: zoom "
     "introduceva un offset di mezza finestra (~8 px) e poteva restituire una "
     "forma diversa da quella del cubo, rompendo il broadcast."),
    ("F11", "medio",
     "I profili verticali mediano in POTENZA e convertono in dB alla fine. "
     "Mediare direttamente i dB introduce un bias sistematico."),
    ("F12", "medio",
     "La colonna di deserto di controllo e' scelta dentro la copertura effettiva "
     "e verificata non vuota. Prima poteva cadere fuori dal chip, restituire zeri "
     "e produrre una separabilita' falsa."),
    ("F13", "medio",
     "Il master e' il 'supermaster' che minimizza il costo di decorrelazione "
     "combinato (baseline + temporale), non semplicemente la prima data."),
    ("F14", "alto",
     "Geocodifica per interpolazione bilineare LOCALE (F36) sulla geolocation "
     "grid, che e' un reticolo regolare 11x21. Prima si usava griddata + "
     "nan_to_num(media): i punti fuori dall'inviluppo venivano riempiti con la "
     "media, deformando la geometria."),
    ("F15", "concettuale",
     "Separazione fra RISOLUZIONE e PRECISIONE. delta_z ~ 227 m e' la separazione "
     "di Rayleigh fra due diffusori nella stessa cella; la precisione sulla quota "
     "di un diffusore dominante e' sigma_h = lambda*R*sin(theta)*sigma_phi / "
     "(4*pi*sigma_Bperp) e con 11 baseline vale alcune decine di metri. La "
     "versione precedente concludeva 'non risolvibile' confondendo le due cose, e "
     "per questo non disegnava alcuna superficie."),
    ("F16", "concettuale",
     "Layover e foreshortening sono calcolati e riportati. Con incidenza 39.4 "
     "gradi (angolo di vista ~34.8) e facce a 52 gradi la faccia rivolta al "
     "sensore e' in layover pieno: e' una proprieta' della geometria, non un "
     "artefatto, e va dichiarata prima di interpretare (ch14)."),
    ("F17", "nuovo",
     "Superficie reale ricostruita pixel per pixel: quota di riferimento dalla "
     "geolocation grid dello .xml piu' quota residua stimata dal periodogramma "
     "multi-baseline sui .tiff VH, con maschera di qualita' sulla coerenza di "
     "fit. Prima il 3D mostrava solo piramidi geometriche ideali disegnate da "
     "costanti scritte a mano."),
    ("F18", "nuovo",
     "Vista 3D con molti gradi di liberta': imbardata, beccheggio, rollio, "
     "traslazione, zoom, prospettiva/ortografica, esagerazione verticale, doppia "
     "soglia in quota, soglia di intensita', piano di sezione, attributo di "
     "colore selezionabile, illuminazione orientabile, livelli commutabili, "
     "controlli da tastiera e rispetto di prefers-reduced-motion."),
    ("F19", "basso",
     "meta.json viene scritto una volta sola invece che due, e l'analisi dei "
     "profili non e' piu' annidata dentro la funzione di salvataggio."),
    ("F20", "basso",
     "Ordinamento painter esplicito dal piu' lontano al piu' vicino, sia per i "
     "voxel sia per le facce della superficie. Prima l'ordine era ascendente su "
     "una profondita' di segno non verificato."),
    ("F21", "basso",
     "Autotest interno (--selftest) che blocca le convenzioni di segno di "
     "apply_shift, del tracker sub-pixel e dell'interpolatore orbitale."),
    ("F22", "critico",
     "La soglia di qualita' sulla superficie e' calibrata sulla DISTRIBUZIONE "
     "NULLA del periodogramma, non fissata a mano. Con 11 baseline e fasi "
     "puramente casuali il picco normalizzato ha gia' mediana 0.56 e 99mo "
     "percentile 0.80: una soglia scritta a mano a 0.35 accettava il 100 per "
     "cento dei pixel, rumore compreso, e faceva sembrare misurato cio' che non "
     "lo era."),
    ("F23", "alto",
     "Stimatore di coerenza corretto. Poiche' y_i = s_i*conj(s_m)/|s_m|, la "
     "media di y_i NON e' la media di s_i*conj(s_m): la coerenza va ricostruita "
     "ripesando per l'ampiezza del master, altrimenti e' uno stimatore diverso e "
     "distorto."),
    ("F24", "critico",
     "L'impronta delle piramidi e' calcolata nella GEOMETRIA RADAR, proiettando "
     "in avanti la superficie 3D attraverso range-Doppler, non come impronta al "
     "suolo. Un punto a 138 m di quota si sposta di h*cos(theta) in slant range, "
     "cioe' ~48 pixel: la maschera al suolo NON copriva i pixel dove l'apice "
     "cade davvero, e quindi calibrazione del segno e validazione guardavano i "
     "pixel sbagliati."),
    ("F28", "nuovo",
     "L'area processata e' stretta sulle sole piramidi in GEOMETRIA RADAR piu' "
     "un margine di 150 m, invece di un rettangolo di 1.8 x 2.2 km su tutta la "
     "piana. Il ritaglio comprende lo spostamento in slant dovuto alla quota, "
     "altrimenti gli apici resterebbero fuori. Le celle da elaborare scendono "
     "da ~50 000 a ~8 000: il calcolo e' piu' rapido, il rendering 3D e' "
     "fluido, e la statistica non e' piu' diluita dal deserto circostante. "
     "--full-scene ripristina il comportamento precedente."),
    ("F29", "critico",
     "La vista 3D non rispondeva ai comandi. Due cause: con "
     "prefers-reduced-motion attivo il ciclo di disegno non partiva affatto e "
     "la pagina restava all'unico fotogramma iniziale, quindi nessun cursore "
     "aveva effetto; e senza quell'impostazione il disegno rifaceva a ogni "
     "fotogramma 8 000 poligoni ordinati, con la stringa di colore ricostruita "
     "ogni volta, a pochi fotogrammi al secondo. Ora il disegno e' a richiesta "
     "su flag di aggiornamento, i vertici sono proiettati una volta sola per "
     "fotogramma invece di quattro volte per poligono, i colori vengono da una "
     "tabella precalcolata, e durante il trascinamento si disegna a "
     "risoluzione ridotta."),
    ("F27", "medio",
     "Despicatura dichiarata della superficie per il rendering: i nodi che si "
     "discostano dalla mediana locale piu' di 3 sigma_h sono errori grossolani "
     "di aggancio su un lobo laterale, valgono centinaia di metri e da soli "
     "dominavano la scala delle quote, rendendo illeggibile una superficie che "
     "per il 98 per cento sta in poche decine di metri. La quota grezza resta "
     "nei .npy; il conteggio delle sostituzioni e' riportato."),
    ("F26", "concettuale",
     "Il canale VH e' documentato come una SCELTA, con la sua misura di costo. "
     "La versione precedente lo giustificava come 'cross-pol, sensibile allo "
     "scattering di volume e alla vibrazione del suolo'. La misura dice il "
     "contrario: nella zona di layover le piramidi emergono di +1.7 dB in VH e "
     "di +8.4 dB in VV, e il picco su Cheope vale +10.5 dB in VH contro +22.2 dB "
     "in VV. Il ritorno delle facce a gradoni e' co-polarizzato; il cross-pol "
     "butta via oltre 10 dB proprio sui bersagli da misurare. Il programma "
     "resta su VH come richiesto, ma misura e riporta il costo, e --pol vv "
     "permette il confronto diretto."),
    ("F25", "alto",
     "La validazione di livello 1 e' una REGRESSIONE fra quota misurata e quota "
     "simulata in geometria radar su tutti i pixel di qualita', non il "
     "percentile 90 della quota dentro l'impronta al suolo. Con le facce in "
     "layover il percentile 90 misurava la coda del rumore e restituiva rilievi "
     "di 300 m dove ce ne sono 138."),
    ("F30", "alto",
     "La superficie misurata non e' piu' disegnata come maglia di poligoni ne' "
     "come campo continuo di celle: e' una NUVOLA DI PUNTI, un punto per nodo "
     "misurato, sia nella vista 3D sia nel pannello matplotlib. Il poligono "
     "univa quattro nodi qualsiasi e riempiva il vuoto fra di essi, quindi "
     "mostrava come superficie continua un campo in cui solo il 23 per cento "
     "delle celle supera la soglia calibrata sul nullo: il 77 per cento del "
     "colore sullo schermo era interpolazione, non misura. I nodi sotto soglia "
     "restano visibili come punti piccoli e spenti, e non passano piu' per "
     "misura solo perche' un vicino affidabile li tirava dentro un poligono; "
     "il filtro sulle celle di qualita' li toglie del tutto. Le sole maglie "
     "rimaste sono i riferimenti dichiarati - piramidi ideali e superficie "
     "simulata - che non sono misure. Il disegno resta sul budget di F29: "
     "i dischetti vengono da un atlante precalcolato, uno per colore della "
     "tabella, e in movimento la nuvola si decima di un fattore due per riga "
     "e colonna."),
    ("F31", "critico",
     "Il canale di coerenza del micro-moto era vuoto. Era calcolato come "
     "|media(dev)| / media(|dev|) con dev = shifts - media(shifts): il "
     "numeratore e' la media di un vettore appena demediato, cioe' zero per "
     "costruzione, su ogni pixel e per qualsiasi dato. Il canale esisteva, "
     "veniva salvato in mm_coh.npy e non conteneva nulla. Ora la domanda e' "
     "quella del modello a 2 gradi di liberta' (ch12): la traccia e' un tono "
     "singolo o rumore bianco? Si risponde con lo spettro della traccia - che "
     "e' il blocco 9 delle fonti applicato alla singola marcia - e la misura e' "
     "la CONCENTRAZIONE spettrale, cioe' la quota di energia nella riga "
     "dominante esclusa la continua: 1 per un tono puro, 1/(N_D-1) = 0.091 per "
     "rumore bianco con N_D = 12. Insieme alla concentrazione viene ora "
     "riportata la frequenza meccanica della riga dominante, coerente con N_D "
     "come frequenza di campionamento della vibrazione (ch11), e la frazione "
     "di celle la cui riga cade dentro la banda osservabile dichiarata dal "
     "banco."),
    ("F32", "nuovo",
     "Tre livelli nuovi nella scena 3D. (a) SUPERFICIE DI RIFERIMENTO: le quote "
     "dei pixel lette dalla geolocation grid degli .annotation.xml e "
     "interpolate con la bilineare locale di F14/F36, disegnate come fondale a "
     "filo di ferro. Non e' un DEM del plateau: sono 231 nodi su tutto lo swath, e "
     "sul ritaglio delle piramidi l'escursione totale e' di poco piu' di due "
     "metri; il valore e' scritto nella barra laterale perche' chi guarda non "
     "la scambi per topografia. (b) PUNTI DI MICRO-MOTO: un punto per pixel la "
     "cui traccia di vibrazione supera la soglia di concentrazione spettrale di "
     "F31, colorato sull'energia. Poggiano sulla superficie .xml perche' il "
     "micro-moto e' un attributo di superficie: il banco di sub-aperture non "
     "gli assegna una profondita', e disegnarlo a quote diverse suggerirebbe "
     "una tomografia che questi dati non fanno. (c) STRUTTURE INTERNE NOTE di "
     "Cheope e Chefren, prese da piramide_cheope_3d.py e piramide_kefren_3d.py, "
     "che restano la sorgente unica di quelle quote. Sono un riferimento "
     "archeologico e sono disegnate come tale, a filo di ferro e mai come "
     "punti: con delta_z verticale di 242 m e camere che misurano metri, "
     "nessuna di esse e' rilevabile da questi dati, e la sovrapposizione serve "
     "a dare la scala di quel divario, non a suggerire una rilevazione."),
    ("F33", "nuovo",
     "Quattro correzioni alla scena 3D (goal utente 2026-08-31). (a) La "
     "superficie SIMULATA (il layover geometrico di simulate_pyramids_radar) "
     "non e' piu' disegnata: restava l'ultima delle maglie di riferimento "
     "che poteva essere confusa con una ricostruzione. La simulazione resta "
     "internamente per la validazione di livello 1 (regressione "
     "misurato-vs-simulato, F25) e per pyr_mask: non e' un risultato, e' un "
     "termine di paragone, e ora si vede solo li'. (b) La superficie di "
     "riferimento .xml e i punti di micro-moto che vi poggiano hanno un "
     "proprio reticolo (xml_ref), l'intera griglia multilooked n_l x n_p, "
     "mai decimato da surface_max_nodes: prima condividevano il reticolo "
     "rado della nuvola misurata e coincidevano con l'intera griglia solo "
     "per la dimensione di questo dataset, non per costruzione. (c) Le "
     "piramidi ideali e le strutture interne note condividono esattamente "
     "centro, angolo di rotazione (azimuth_deg, oggi 0 per tutte e tre) e "
     "quota di base, e sono proiettate dalla stessa project() rigida della "
     "vista: restano solidali sotto qualunque imbardata, beccheggio o "
     "rollio, verificato a vista ruotando la scena. (d) I voxel sono ora "
     "colorati per indice di PIENO/VUOTO (solidity_index, F31: energia del "
     "periodogramma x coerenza x (1 - micro-moto), quindi il micro-moto "
     "entra nel colore), non piu' per la sola ampiezza; il valore esportato "
     "e' ristirato per percentili solo ai fini del colore, il dato "
     "scientifico in solidity.npy resta la formula dichiarata. Resta un "
     "discriminante multi-attributo, non una rilevazione di cavita' risolta "
     "in profondita' (avvertenza invariata, ch08-ch09)."),
    ("F34", "alto",
     "La superficie .xml (F33) non si vedeva come pixel. Era disegnata come "
     "maglia di quadrilateri semitrasparenti: su un terreno quasi piatto "
     "(poco piu' di 2 m di escursione sul ritaglio) il riempimento si "
     "fondeva otticamente con la griglia di sfondo a 100 m, e i singoli "
     "pixel del reticolo xml_ref non erano piu' riconoscibili come tali, "
     "anche se il dato sotto era gia' corretto e completo (F33). Ora e' un "
     "punto per pixel, come la nuvola misurata (F30), colorato sulla "
     "propria quota con una rampa grigio-acciaio dedicata (mai la rampa "
     "calda della misura, cosi' il riferimento resta riconoscibile), "
     "disegnato con lo stesso atlante precalcolato del resto della scena "
     "(F29) per restare fluido sugli ~8 000 pixel del ritaglio anche in "
     "movimento."),
    ("F35", "chiarimento",
     "F34 documentava il quasi-piano di xml_ref (~2 m di escursione) come un "
     "fatto del ritaglio, senza spiegare la causa: i punti della geolocation "
     "grid dello .xml sono posizioni fisiche di pixel con una quota propria, "
     "NON un piano -- sull'intera scena il rilievo vero arriva a centinaia di "
     "metri (raw_gcp_nodes). Ma quella grid ha solo 11 x 21 = 231 nodi su "
     "tutta la scena, spaziati circa 1500 linee x 1300 pixel, mentre il "
     "ritaglio delle piramidi ne occupa una piccola frazione di UNO: nessun "
     "nodo reale cade abbastanza vicino da portare rilievo dentro la spline "
     "locale (height_ref/xml_ref), che quindi esce correttamente quasi piatta "
     "-- limite di risoluzione del riferimento .xml, non errore di calcolo. "
     "Ora i nodi grezzi (non interpolati) sono esposti nel payload e in un "
     "pannello dedicato della UI, con la distanza del nodo piu' vicino e il "
     "rilievo vero sulla scena, cosi' il limite e' ispezionabile invece che "
     "implicito."),
    ("F36", "critico",
     "F35 spiegava perche' xml_ref esce quasi piatta, ma Geocoder usava ancora "
     "una spline bicubica GLOBALE (kx=ky=3) su tutti e 231 i nodi della scena: "
     "verificato sui quattro nodi VERI che racchiudono il ritaglio delle "
     "piramidi (linee 7540/9048, pixel 1304/2608, quote 80.0/51.0/64.0/64.0 m "
     "-- 29 m di escursione fra loro), la spline dava 70.0 m al centro del "
     "ritaglio mentre la bilineare fra quei quattro nodi da' 63.7 m: 6 m di "
     "differenza per curvatura importata da nodi a decine di km, non dal "
     "terreno locale. Il ritaglio sta interamente dentro UNA cella del "
     "reticolo (F35): la geocodifica ora usa RegularGridInterpolator "
     "bilineare, che per un punto in una cella usa esclusivamente i quattro "
     "nodi reali di quella cella -- i valori letti dallo .xml corrispondente, "
     "non un adattamento globale. L'estrapolazione lineare oltre il bordo del "
     "reticolo (fill_value=None) sostituisce la vecchia estrapolazione della "
     "spline, con lo stesso comportamento definito ovunque."),
    ("F37", "richiesta",
     "Rimosso lo strato 3D 'Superficie .xml' (drawRef/buildRefColors/XCIDX): "
     "su questo ritaglio la quota di riferimento varia di poco piu' di 1 m "
     "(F35/F36), quindi disegnata come nuvola propria non aggiungeva nulla "
     "che il pannello 'Nodi grezzi della geolocation grid' non mostrasse gia' "
     "in forma numerica ispezionabile. NOTA sui file .xml: non esiste un tag "
     "<position> con la quota dei singoli pixel -- i 17 nodi <position> "
     "dell'annotation.xml sono i VETTORI DI STATO ORBITALI del satellite "
     "(coordinate ECEF x/y/z), non posizioni al suolo. L'unica quota per "
     "pixel nello .xml e' <geolocationGridPoint><height>, 231 nodi per "
     "scena: e' quella che height_ref/xml_ref usano (F36), presa dal file "
     ".annotation.xml del MASTER -- l'unico geometricamente corretto da "
     "usare per l'intero stack, perche' tutte le date sono coregistrate sul "
     "suo reticolo di pixel."),
    ("F38", "richiesta",
     "Aggiunto il layer 'suolo (DEM)', portato da _fetch_dem() di "
     "piramide_acustica_vh.py / piramide_unificato.py (Piramid_V3): un DEM "
     "ESTERNO vero (Copernicus via l'API di elevazione Open-Meteo), campionato "
     "su un reticolo fitto SOLO nell'area del ritaglio e interpolato "
     "(RegularGridInterpolator bilineare) sulla stessa griglia (east, north) "
     "di height_ref. Risolve il rilievo che la geolocation grid dello .xml "
     "non puo' dare a questa scala (F35/F36): verificato, 18..94 m di "
     "escursione su un campione dello stesso ritaglio, contro gli ~1 m della "
     "bilineare .xml. E' un RIFERIMENTO esterno (fetch_dem_suolo(), cache in "
     "dem_suolo.npz, disattivabile con cfg.fetch_dem=False), non entra in "
     "nessun calcolo tomografico: la geometria .xml resta l'unica usata per "
     "fase, k_z e geocodifica. CORRETTO da F39: vedi sotto."),
    ("F39", "richiesta",
     "Il 'suolo (DEM)' di F38 non seguiva il profilo delle piramidi e non era "
     "sul datum dei dati. Tre difetti, tutti verificati sui dati. (a) Il DEM "
     "esterno le piramidi non le contiene: interrogando l'API di elevazione "
     "esattamente sui tre apici noti restituisce 66 / 76 / 78 m, cioe' il "
     "plateau nudo -- l'apice di Cheope, che sta a ~198 m, usciva 130 m piu' "
     "in basso e il 'suolo' passava DENTRO le piramidi invece di seguirne le "
     "facce. Non era il passo del reticolo: anche interrogando il singolo "
     "punto dell'apice il valore resta 66 m. (b) Il livello veniva da un "
     "dataset con un datum verticale suo: 68.6 m di mediana sul ritaglio e "
     "75.5 m sotto Cheope, contro i 63.1..64.2 m delle quote .xml, cioe' fino "
     "a ~12 m di scarto rispetto alla superficie misurata, che poggia sul "
     "datum .xml. (c) Se lo scaricamento falliva il layer spariva del tutto. "
     "Ora ground_dem_suolo() lo costruisce in tre pezzi dichiarati: DATUM = "
     "le quote lette dai file di stack_slc (<geolocationGridPoint><height> "
     "degli annotation.xml, bilineare locale del master, F36), lette da TUTTE "
     "le date come controllo di consistenza (plateau_heights_from_stack(): 11 "
     "date VH IW2, tutte 64.0 m sul nodo grezzo piu' vicino, che cade a "
     "2.2-3.2 km); RILIEVO = il DEM esterno riportato su quel datum togliendo "
     "il suo scarto mediano sull'impronta delle piramidi, e senza rete resta "
     "la bilineare .xml, cosi' il layer c'e' comunque; PIRAMIDI = "
     "pyramid_profile_enu(), z = base + altezza * (1 - max(|dx|,|dy|)/semilato) "
     "sul reticolo ENU, con la base presa dal terreno sotto ciascuna impronta "
     "(non dal base_alt_m di letteratura) perche' il layer sia continuo con il "
     "terreno intorno, composto con np.maximum. Corretto anche il reticolo di "
     "query: l'API accetta 100 coordinate per richiesta e prima se ne mandava "
     "una sola, quindi il reticolo era bloccato a 10x10 (~170 m di passo); ora "
     "la richiesta e' spezzata in blocchi e cfg.dem_grid vale 24 (~70 m). "
     "Interruttori: cfg.suolo_dem (il layer), cfg.fetch_dem (il rilievo "
     "esterno), cfg.dem_pyramids (il profilo delle piramidi)."),
    ("F41", "richiesta",
     "La nuvola di voxel mette il 46.4% dei punti SOPRA la superficie e il "
     "23.4% sopra l'apice di Cheope, e la pagina non diceva perche'. Non e' "
     "un difetto del disegno, sono due fatti sovrapposti. (1) L'asse z e' un "
     "intervallo di RICERCA simmetrico attorno al datum (+-elev_max = 400 m): "
     "meta' sta in aria per costruzione e nulla nell'inversione obbliga la "
     "soluzione a stare sottoterra. (2) Dove finisce l'energia lo decide la "
     "PSF verticale della pila, e con 8 baseline quella PSF non discrimina. "
     "Misurato sul ritaglio delle piramidi da vertical_lobe_profile(): picco "
     "a -12.5 m, cioe' sulla superficie, corretto; lobo principale a -3 dB "
     "largo 150 m (da -87.5 a +62.5 m); ma il LOBO LATERALE PEGGIORE sta a "
     "-1.16 dB, poco piu' di un dB sotto il picco, e ai bordi dell'asse la "
     "curva vale ancora -1.9 / -2.7 dB, per 5.07 dB di contrasto su tutti gli "
     "800 m dell'asse. Un picco che i suoi lobi quasi eguagliano non "
     "localizza nulla in profondita'. Aggravante di resa: l'asse e' "
     "campionato a 3.125 m, 91 volte piu' fitto della delta_z di Rayleigh "
     "(283.5 m), quindi la nuvola sembra dettagliata per sovracampionamento. "
     "build_html disegna ora quel profilo nel pannello 'Lobi verticali' "
     "accanto ai comandi dei voxel, con la riga dei -3 dB, la barra di "
     "delta_z in scala e le due frazioni di voxel sopra il suolo. La nuvola "
     "NON viene ritagliata a z<=0: nasconderne meta' darebbe un'immagine "
     "piu' credibile e piu' falsa. La correzione vera e' piu' baseline, non "
     "un filtro sul disegno."),
    ("F42", "critico",
     "Ogni data veniva LETTA alla stessa finestra assoluta del master e "
     "riallineata dopo con apply_shift(), che e' una traslazione CIRCOLARE via "
     "rampa di fase. Vale per frazioni di pixel, non per l'inquadramento del "
     "prodotto: sulla stessa traccia relativa 58 il bersaglio di Giza cade a "
     "-4510 linee nei prodotti S1A e a +30 nei S1D rispetto al master S1C, e "
     "fino a -108 pixel in range. Con 89 linee di ritaglio, per 16 date su 28 "
     "lo 'slave' era un pezzo di deserto a 63 km di distanza, ripiegato su se "
     "stesso dalla traslazione circolare e per giunta dirampato con i "
     "parametri del burst sbagliato: interferogrammi di puro rumore, e la "
     "pila multi-missione (che e' esattamente cio' che allarga l'escursione "
     "di baseline da 163 a 273 m) non era utilizzabile. Anche restando su una "
     "sola missione l'offset in range di alcune decine di pixel faceva "
     "avvolgere una colonna su otto del ritaglio. Ora l'offset INTERO entra "
     "nella finestra di lettura e nell'indice di burst passato al deramping "
     "TOPS; alla rampa di fase resta il solo residuo sub-pixel. Se il "
     "ritaglio spostato scavalcherebbe il bordo del burst o del prodotto la "
     "data viene scartata e detta a voce alta, perche' il deramping TOPS non "
     "e' definito a cavallo di due burst."),
    ("F43", "alto",
     "_omogenea() separava le acquisizioni per il solo verso di passaggio. "
     "Due tracce ascendenti DIVERSE sono incompatibili quanto un'ascendente e "
     "una discendente: guardano lo stesso punto da posizioni orbitali "
     "distanti centinaia di chilometri, quindi baseline fuori scala e "
     "coerenza nulla. Il criterio ora e' verso + orbita relativa, letta dal "
     "manifest.safe accanto all'annotation (con ricaduta sul solo verso nel "
     "layout piatto, dove il manifest non c'e')."),
    ("F44", "medio",
     "Il default di --gamma-min era 0.35, cioe' proprio la costante che F22 "
     "aveva tolto, e main() lo passava sempre alla Config (che ha default "
     "0.0). Restava innocuo solo perche' la soglia calibrata sul nullo e' piu' "
     "alta di 0.35 con le pile finora usate; con una pila abbastanza numerosa "
     "la soglia nulla scende e il pavimento tornerebbe a decidere lui. Il "
     "default e' ora 0.0, coerente con Config."),
    ("F45", "medio",
     "ground_dem_suolo() usava `m_lat` nella riga di resoconto del reticolo "
     "DEM, ma in quello scope il nome non esiste: NameError. Il ramo scatta "
     "solo al PRIMO scaricamento riuscito (nessuna cache .npz) e con verbose, "
     "quindi bastava avere gia' il file per non vederlo mai. Ora usa la "
     "costante M_PER_GRADO_LAT. Trovato con mypy, non a occhio."),
)


# --------------------------------------------------------------------------
# Geometria nota delle piramidi -- SOLO come riferimento di validazione
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Pyramid:
    """Geometria nota da letteratura archeologica.

    Serve ESCLUSIVAMENTE come verita' di riferimento per la validazione di
    livello 1 (geometria) e come sovrapposizione nel rendering. Non entra in
    nessun calcolo che produce la superficie: quella viene dai .tiff.
    """
    name: str
    lat: float
    lon: float
    base_side_m: float
    height_m: float
    base_alt_m: float
    face_slope_deg: float
    azimuth_deg: float = 0.0


PYRAMIDS: Tuple[Pyramid, ...] = (
    Pyramid("Cheope (Khufu)",       29.979235, 31.134202, 230.33, 138.5, 60.0, 51.84),
    Pyramid("Chefren (Khafre)",     29.976111, 31.130833, 215.25, 136.4, 70.0, 53.13),
    Pyramid("Micerino (Menkaure)",  29.972500, 31.128056, 102.20,  61.0, 66.0, 51.34),
)


@dataclass
class Config:
    stack_dir: str = DEFAULT_STACK
    swath: str = "iw2"
    polarisation: str = "vh"
    platform: str = "s1c"

    n_dates: int = 11

    # --- finestra di analisi ------------------------------------------------
    # F28: per default si processa SOLO l'area delle piramidi in geometria
    # radar, piu' un anello di margine per la calibrazione e il controllo.
    area_margin_m: float = 150.0
    full_scene: bool = False          # True = vecchio rettangolo su tutta la piana
    chip_pad_range_m: float = 900.0   # usati solo con full_scene
    chip_pad_azim_m: float = 1100.0
    mm_lines: int = 1024              # F07: linee azimutali per il micro-moto

    # --- asse di elevazione -------------------------------------------------
    n_elev: int = 257
    elev_max_m: float = 400.0

    # --- micro-moto ---------------------------------------------------------
    n_d: int = 12
    guard_fraction: float = 0.5       # B_DL / B_cD, imposto dalle fonti
    b_shift_fraction: float = 0.09    # B_shift / B_cD; deve stare sotto guard/3
    corr_window: int = 32

    # --- multilooking -------------------------------------------------------
    look_range: int = 4               # 4 * 2.33 m / sin(theta) ~ 14.7 m a terra
    look_azimuth: int = 1             # 1 * 13.955 m -> celle quasi quadrate
    coh_window: int = 7

    # --- qualita' -----------------------------------------------------------
    gamma_null_q: float = 99.0        # F22: percentile della distribuzione nulla
    gamma_null_trials: int = 20000
    gamma_min: float = 0.0            # eventuale pavimento aggiuntivo
    surface_max_nodes: int = 12000    # nodi esportati nella nuvola 3D (F30)

    out_dir: str = "out_piramidi_v02"
    html_name: str = "tomografia_piramidi_3d.html"
    seed: int = 20260828

    # --- F38/F39: suolo (DEM) ------------------------------------------------
    # Layer di riferimento visivo ("suolo (DEM)"): la geometria .xml resta
    # l'unica usata nel calcolo tomografico (fase, k_z, geocodifica).
    # F39: il layer e' terreno + PROFILO DELLE PIRAMIDI, con il datum preso
    # dalle quote lette negli annotation.xml dello stack (ground_dem_suolo()).
    # --- F40: radiometria dai .xml di calibrazione e rumore ----------------
    # I prodotti .SAFE portano calibration-*.xml (sigmaNought) e noise-*.xml
    # (NESZ). Applicarli converte i conteggi DN in sigma0 e toglie il rumore
    # termico dall'intensita'. La fase non cambia -- si divide per un reale
    # positivo -- quindi baseline, k_z e quote restano identici: cambia solo
    # la radiometria, che in VH cross-pol e' proprio dove sta il problema.
    calibrazione: bool = True

    suolo_dem: bool = True            # produce il layer (anche senza rete)
    fetch_dem: bool = True            # aggiunge il rilievo del DEM esterno
    dem_pyramids: bool = True         # F39: il suolo segue le piramidi
    dem_grid: int = 24                # nodi per lato del reticolo di query API
    dem_cache_name: str = "dem_suolo.npz"


# ==========================================================================
# 1.  Orbita: interpolazione di Lagrange e stato a Doppler zero
# ==========================================================================

def ecef_from_llh(lat, lon, h=0.0) -> np.ndarray:
    """Da lat/lon/quota WGS84 a coordinate cartesiane geocentriche.

    Vettorizzata: accetta scalari o array e restituisce (..., 3)."""
    a, f = 6378137.0, 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    la = np.radians(np.asarray(lat, dtype=np.float64))
    lo = np.radians(np.asarray(lon, dtype=np.float64))
    hh = np.asarray(h, dtype=np.float64)
    n = a / np.sqrt(1.0 - e2 * np.sin(la) ** 2)
    return np.stack([
        (n + hh) * np.cos(la) * np.cos(lo),
        (n + hh) * np.cos(la) * np.sin(lo),
        (n * (1.0 - e2) + hh) * np.sin(la),
    ], axis=-1)


@dataclass
class Orbit:
    """Vettori di stato del sensore, con interpolazione di Lagrange (F03)."""
    t: np.ndarray                 # [n] secondi dal primo vettore di stato
    pos: np.ndarray               # [n, 3] m, ECEF
    vel: np.ndarray               # [n, 3] m/s, ECEF
    t0: datetime

    order: int = 8

    def _lagrange(self, tq: np.ndarray, vals: np.ndarray) -> np.ndarray:
        """Interpolazione di Lagrange di ordine `order` sui campioni piu' vicini.

        I vettori di stato Sentinel-1 distano 10 s: su quell'intervallo l'orbita
        non e' affatto lineare e np.interp introduce errori metrici, dello stesso
        ordine di grandezza delle baseline che stiamo misurando."""
        tq = np.atleast_1d(np.asarray(tq, dtype=np.float64))
        n = len(self.t)
        k = min(self.order, n)
        # indice del campione immediatamente precedente
        idx = np.clip(np.searchsorted(self.t, tq) - k // 2, 0, n - k)
        sel = idx[:, None] + np.arange(k)[None, :]          # [nq, k]
        ts = self.t[sel]                                    # [nq, k]
        # pesi di Lagrange
        diff = tq[:, None] - ts                             # [nq, k]
        w = np.ones_like(diff)
        for j in range(k):
            num = np.delete(diff, j, axis=1)
            den = ts[:, j:j + 1] - np.delete(ts, j, axis=1)
            w[:, j] = np.prod(num / den, axis=1)
        return np.einsum("qk,qkc->qc", w, vals[sel])

    def state(self, tq) -> Tuple[np.ndarray, np.ndarray]:
        return self._lagrange(tq, self.pos), self._lagrange(tq, self.vel)

    def zero_doppler_time(self, targets: np.ndarray, n_coarse: int = 6000) -> np.ndarray:
        """Istante di Doppler zero (P(t) - T) . V(t) = 0 per ogni target.

        Ricerca grossolana su griglia + interpolazione lineare del passaggio per
        zero: f(t) e' regolare e quasi lineare sull'intervallo di 4 ms della
        griglia, quindi la radice interpolata e' esatta ben oltre il necessario."""
        tg = np.linspace(self.t[0], self.t[-1], n_coarse)
        pg, vg = self.state(tg)                                    # [n_coarse, 3]
        tgt = np.atleast_2d(targets)                               # [m, 3]
        # f[m, n_coarse]
        f = np.einsum("nc,nc->n", pg, vg)[None, :] \
            - tgt @ vg.T
        sign = np.signbit(f)
        cross = np.diff(sign.astype(np.int8), axis=1) != 0
        out = np.empty(len(tgt), dtype=np.float64)
        for i in range(len(tgt)):
            j = np.flatnonzero(cross[i])
            if len(j) == 0:
                out[i] = tg[np.argmin(np.abs(f[i]))]
                continue
            j = j[np.argmin(np.abs(f[i, j]))]
            f0, f1 = f[i, j], f[i, j + 1]
            out[i] = tg[j] + (tg[j + 1] - tg[j]) * f0 / (f0 - f1)
        return out

    def slant_range(self, targets: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Slant range a Doppler zero, con posizione e velocita' del sensore."""
        tz = self.zero_doppler_time(targets)
        p, v = self.state(tz)
        r = np.linalg.norm(p - np.atleast_2d(targets), axis=1)
        return r, p, v


def read_orbit(path: str) -> Orbit:
    root = ET.parse(path).getroot()
    svs = root.findall(".//orbitList/orbit") or root.findall(".//orbit")
    if not svs:
        raise ValueError(f"nessun vettore di stato in {path}")
    times = [datetime.strptime(_txt(s, "time")[:26], "%Y-%m-%dT%H:%M:%S.%f")
             for s in svs]
    pos = np.array([[float(_txt(s, f"position/{c}")) for c in "xyz"] for s in svs])
    vel = np.array([[float(_txt(s, f"velocity/{c}")) for c in "xyz"] for s in svs])
    t0 = times[0]
    ts = np.array([(t - t0).total_seconds() for t in times])
    order = int(np.clip(len(ts), 2, 8))
    return Orbit(t=ts, pos=pos, vel=vel, t0=t0, order=order)


# ==========================================================================
# 2.  Geocodifica: interpolazione BILINEARE LOCALE sulla geolocation grid
#     regolare  (F14, corretto F36)
# ==========================================================================

class Geocoder:
    """La geolocation grid di un annotation IW e' un reticolo REGOLARE
    (11 linee x 21 pixel = 231 nodi) che copre l'INTERA scena (~250 x 250 km),
    con nodi reali spaziati circa 1500 linee x 1300 pixel.

    F36: prima veniva interpolata con una spline bicubica GLOBALE (kx=ky=3) su
    tutti e 231 i nodi. Una spline globale e' influenzata anche dai nodi
    lontani: verificato sui dati, agli angoli VERI che racchiudono il ritaglio
    delle piramidi (7540/9048 linee x 1304/2608 pixel, quota 80.0 / 51.0 /
    64.0 / 64.0 m -- un'escursione di 29 m) la spline dava 70.0 m al centro
    del ritaglio, mentre la bilineare LOCALE fra quei quattro nodi veri da'
    63.7 m: 6 m di differenza dovuti a curvatura importata da nodi a decine di
    km di distanza, non al terreno locale. Per un punto che sta dentro UNA
    sola cella del reticolo (il caso di questo intero programma: F35) la
    bilineare sui quattro nodi reali che la racchiudono e' la sola
    interpolazione che usa esclusivamente i valori di coordinate letti dallo
    .xml corrispondenti a quell'area, come richiesto. Resta definita ovunque
    (fill_value=None estrapola linearmente oltre il bordo del reticolo, cosi'
    come la spline non restituiva NaN fuori dall'inviluppo convesso)."""

    def __init__(self, ann: S1Annotation):
        from scipy.interpolate import RegularGridInterpolator

        root = ET.parse(ann.path).getroot()
        gp = root.findall(".//geolocationGrid//geolocationGridPoint")
        height = np.array([float(_txt(p, "height")) for p in gp])

        lines = np.unique(ann.geo_line)
        pixels = np.unique(ann.geo_pixel)
        n_l, n_p = len(lines), len(pixels)
        if n_l * n_p != len(ann.geo_line):
            raise ValueError("geolocation grid non regolare: "
                             f"{n_l}x{n_p} != {len(ann.geo_line)}")

        li = np.searchsorted(lines, ann.geo_line)
        pi = np.searchsorted(pixels, ann.geo_pixel)

        def grid(v: np.ndarray) -> np.ndarray:
            g = np.zeros((n_l, n_p))
            g[li, pi] = v
            return g

        def interp(g: np.ndarray) -> "RegularGridInterpolator":
            return RegularGridInterpolator(
                (lines, pixels), g, method="linear",
                bounds_error=False, fill_value=None)

        self.lines, self.pixels = lines, pixels
        self._lat = interp(grid(ann.geo_lat))
        self._lon = interp(grid(ann.geo_lon))
        self._hgt = interp(grid(height))
        self._inc = interp(grid(ann.geo_incidence))
        self.height_nodes = grid(height)
        # F35: i nodi GREZZI del reticolo (gli stessi usati sopra per la
        # bilineare locale), per poterli esporre anche come punti fisici veri
        # -- vedi raw_gcp_nodes() piu' sotto.
        self.lat_nodes = grid(ann.geo_lat)
        self.lon_nodes = grid(ann.geo_lon)
        self.ann = ann

    def _ev(self, interp, line: np.ndarray, pixel: np.ndarray) -> np.ndarray:
        line = np.asarray(line, dtype=np.float64)
        pixel = np.asarray(pixel, dtype=np.float64)
        shape = line.shape
        pts = np.column_stack([line.ravel(), pixel.ravel()])
        return interp(pts).reshape(shape)

    def llh(self, line, pixel) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (self._ev(self._lat, line, pixel),
                self._ev(self._lon, line, pixel),
                self._ev(self._hgt, line, pixel))

    def incidence(self, line, pixel) -> np.ndarray:
        return self._ev(self._inc, line, pixel)

    def latlon_to_line_pixel(self, lat: float, lon: float,
                             iters: int = 40) -> Tuple[float, float]:
        """Inversione della geocodifica per Newton 2D sulla spline.

        Molto piu' accurata dell'interpolazione lineare su triangoli usata prima,
        e sempre definita anche vicino al bordo del reticolo."""
        l = float(np.mean(self.lines))
        p = float(np.mean(self.pixels))
        dl = max(1.0, (self.lines[-1] - self.lines[0]) * 1e-4)
        dp = max(1.0, (self.pixels[-1] - self.pixels[0]) * 1e-4)
        for _ in range(iters):
            la, lo, _ = self.llh(l, p)
            f = np.array([float(la) - lat, float(lo) - lon])
            if np.hypot(*f) < 1e-11:
                break
            la_l, lo_l, _ = self.llh(l + dl, p)
            la_p, lo_p, _ = self.llh(l, p + dp)
            j = np.array([[(float(la_l) - float(la)) / dl, (float(la_p) - float(la)) / dp],
                          [(float(lo_l) - float(lo)) / dl, (float(lo_p) - float(lo)) / dp]])
            try:
                step = np.linalg.solve(j, f)
            except np.linalg.LinAlgError:
                break
            l -= float(step[0])
            p -= float(step[1])
            l = float(np.clip(l, self.lines[0], self.lines[-1]))
            p = float(np.clip(p, self.pixels[0], self.pixels[-1]))
        return l, p


# ==========================================================================
# 3.  Baseline e budget
# ==========================================================================

@dataclass
class Baseline:
    date: str
    b_perp: float
    b_par: float
    b_temp: float                 # giorni rispetto al master
    slant_range: float


def compute_baselines(
    orbits: Sequence[Tuple[str, Orbit]], target: np.ndarray, master_idx: int,
) -> List[Baseline]:
    """Baseline ortogonale, parallela e temporale rispetto al master.

    B_perp e' la componente della baseline ortogonale alla linea di vista e
    all'asse di volo: e' l'apertura sintetica della tomografia multi-baseline.
    B_par non porta informazione di elevazione."""
    tgt = np.atleast_2d(target)
    states = []
    for date, orb in orbits:
        rng, pos_s, vel_s = orb.slant_range(tgt)
        states.append((date, pos_s[0], vel_s[0], float(rng[0]), orb.t0))

    _, p_m, v_m, _, t0_m = states[master_idx]
    los = (target - p_m)
    los = los / np.linalg.norm(los)
    normal = np.cross(los, v_m / np.linalg.norm(v_m))
    normal /= np.linalg.norm(normal)

    out: List[Baseline] = []
    for date, p, _v, r, t0 in states:
        b = p - p_m
        b_par = float(np.dot(b, los))
        b_perp = float(np.dot(b - b_par * los, normal))
        out.append(Baseline(date, b_perp, b_par,
                            (t0 - t0_m).total_seconds() / 86400.0, r))
    return out


def pick_supermaster(baselines: Sequence[Baseline]) -> int:
    """F13: master che minimizza il costo di decorrelazione combinato.

    Il costo e' la somma su tutte le altre date di (B_perp/B_crit)^2 +
    (dt/tau)^2. La scelta non cambia l'escursione totale di baseline -- quella e'
    invariante -- ma riduce la decorrelazione media, quindi migliora la coerenza
    con cui si stima la quota."""
    bp = np.array([b.b_perp for b in baselines])
    dt = np.array([b.b_temp for b in baselines])
    b_crit, tau = 5000.0, 60.0
    cost = [np.sum(((bp - bp[i]) / b_crit) ** 2 + ((dt - dt[i]) / tau) ** 2)
            for i in range(len(baselines))]
    return int(np.argmin(cost))


@dataclass
class TomoBudget:
    wavelength: float
    slant_range: float
    incidence_deg: float
    b_spread: float
    b_std: float
    n_baselines: int
    delta_z_slant: float
    delta_z_vertical: float
    ambiguity_height: float
    target_height: float
    gamma_typ: float
    sigma_h: float                     # F15: precisione, non risoluzione

    @property
    def cells_over_target(self) -> float:
        return self.target_height / self.delta_z_vertical

    @property
    def resolves_interior(self) -> bool:
        """Servono ~3 celle sull'altezza per separare strutture interne."""
        return self.cells_over_target >= 3.0

    @property
    def surface_measurable(self) -> bool:
        """La superficie e' misurabile se la PRECISIONE sulla quota di un
        diffusore dominante e' molto piu' fine dell'altezza del target."""
        return self.sigma_h <= self.target_height / 4.0

    def as_text(self) -> str:
        return "\n".join([
            "-" * 76,
            "BUDGET TOMOGRAFICO MULTI-BASELINE",
            "-" * 76,
            f"  lunghezza d'onda radar                = {self.wavelength * 100:10.2f} cm",
            f"  slant range                      R    = {self.slant_range / 1000:10.2f} km",
            f"  angolo di incidenza              th   = {self.incidence_deg:10.2f} gradi",
            f"  baseline impilate                     = {self.n_baselines:10d}",
            f"  escursione baseline ortogonali        = {self.b_spread:10.1f} m",
            f"  deviazione standard delle baseline    = {self.b_std:10.1f} m",
            "",
            "  RISOLUZIONE  (separazione di Rayleigh fra due diffusori)",
            f"    delta_z = lambda*R/(2*B)  in slant  = {self.delta_z_slant:10.1f} m",
            f"    proiettata in verticale             = {self.delta_z_vertical:10.1f} m",
            f"    altezza di ambiguita'               = {self.ambiguity_height:10.1f} m",
            f"    celle sull'altezza del target       = {self.cells_over_target:10.2f}",
            f"    >> struttura interna separabile?      "
            f"{'SI' if self.resolves_interior else 'NO':>10s}",
            "",
            "  PRECISIONE  (quota di UN diffusore dominante)  -- F15",
            f"    coerenza tipica assunta       gamma = {self.gamma_typ:10.2f}",
            f"    sigma_h = lambda*R*sin(th)*sigma_phi/(4*pi*sigma_B)",
            f"                                        = {self.sigma_h:10.1f} m",
            f"    altezza del target                  = {self.target_height:10.1f} m",
            f"    >> superficie reale misurabile?       "
            f"{'SI' if self.surface_measurable else 'NO':>10s}",
            "",
            "  Le due righe non si contraddicono: separare due diffusori dentro la",
            "  stessa cella e stimare la quota di quello dominante sono problemi",
            "  diversi. Confonderli e' l'errore concettuale corretto in F15.",
            "-" * 76,
        ])


def compute_tomo_budget(
    baselines: Sequence[Baseline], ann: S1Annotation, incidence_deg: float,
    target_height: float, gamma_typ: float = 0.6,
) -> TomoBudget:
    bp = np.array([b.b_perp for b in baselines])
    spread = float(bp.max() - bp.min())
    b_std = float(bp.std(ddof=1)) if len(bp) > 1 else 0.0
    r = float(np.mean([b.slant_range for b in baselines]))
    lam = ann.wavelength
    sin_t = math.sin(math.radians(incidence_deg))

    dz_slant = lam * r / (2.0 * spread) if spread > 0 else float("inf")
    dz_vert = dz_slant / sin_t

    gaps = np.diff(np.sort(bp))
    gaps = gaps[gaps > 1e-6]
    med_gap = float(np.median(gaps)) if len(gaps) else spread
    amb = lam * r / (2.0 * med_gap) / sin_t if med_gap > 0 else float("inf")

    # F15 -- limite di Cramer-Rao sulla quota di un diffusore dominante.
    # sigma_phi per interferogramma multilook, poi propagata dalla pendenza
    # dphi/dh = 4*pi*B_perp/(lambda*R*sin(theta)).
    n = max(len(bp), 2)
    g2 = float(np.clip(gamma_typ, 1e-3, 0.999)) ** 2
    sigma_phi = math.sqrt((1.0 - g2) / (2.0 * n * g2))
    sigma_h = (lam * r * sin_t * sigma_phi / (4.0 * math.pi * max(b_std, 1e-6))) \
        if b_std > 0 else float("inf")

    return TomoBudget(
        wavelength=lam, slant_range=r, incidence_deg=incidence_deg,
        b_spread=spread, b_std=b_std, n_baselines=len(baselines),
        delta_z_slant=dz_slant, delta_z_vertical=dz_vert,
        ambiguity_height=amb, target_height=target_height,
        gamma_typ=gamma_typ, sigma_h=sigma_h,
    )


def look_angle_deg(incidence_deg: float, sat_height_m: float = 693_000.0,
                   earth_radius_m: float = 6_371_000.0) -> float:
    """Angolo di vista (off-nadir) dall'angolo di incidenza al suolo.

    NON e' 90 - incidenza: quella e' l'elevazione locale. Con Terra sferica vale
    sin(look) = R_e * sin(incidenza) / (R_e + h_sat)."""
    s = earth_radius_m * math.sin(math.radians(incidence_deg)) / \
        (earth_radius_m + sat_height_m)
    return math.degrees(math.asin(min(max(s, -1.0), 1.0)))


def layover_report(incidence_deg: float,
                   sat_height_m: float = 693_000.0) -> Dict[str, Any]:
    """F16 -- diagnosi geometrica prima di qualunque interpretazione.

    Criteri corretti, con theta = angolo di INCIDENZA al suolo e alpha =
    pendenza della faccia:
        alpha < theta          -> foreshortening (faccia rivolta al sensore)
        alpha > theta          -> layover
        alpha > 90 - theta     -> ombra (faccia opposta al sensore)

    Il confronto va fatto con l'incidenza, non con l'angolo di vista ne' con
    90 - incidenza: 90 - theta e' la soglia dell'OMBRA, non del layover."""
    th = incidence_deg
    look = look_angle_deg(incidence_deg, sat_height_m)
    rows = []
    for p in PYRAMIDS:
        a = p.face_slope_deg
        near = "layover" if a > th else ("foreshortening" if a > 0 else "piano")
        far = "ombra" if a > (90.0 - th) else "illuminata"
        # la faccia vicina si comprime in slant range in questo intervallo
        w_ground = p.base_side_m / 2.0
        slant_extent = abs(w_ground * math.sin(math.radians(th))
                           - p.height_m * math.cos(math.radians(th)))
        rows.append({
            "piramide": p.name,
            "pendenza_faccia_deg": a,
            "incidenza_locale_faccia_vicina_deg": round(th - a, 2),
            "incidenza_locale_faccia_lontana_deg": round(th + a, 2),
            "faccia_vicina": near,
            "faccia_lontana": far,
            "estensione_slant_faccia_vicina_m": round(slant_extent, 1),
            "compressione": round(w_ground / max(slant_extent, 1e-3), 1),
        })
    return {
        "angolo_di_incidenza_deg": round(th, 2),
        "angolo_di_vista_off_nadir_deg": round(look, 2),
        "soglia_layover_deg": round(th, 2),
        "soglia_ombra_deg": round(90.0 - th, 2),
        "piramidi": rows,
        "nota": ("Tutte e tre le facce rivolte al sensore superano l'angolo di "
                 "incidenza: sono in layover pieno, e le facce opposte sono in "
                 "ombra. E' geometria, non un artefatto del processing. Le fonti "
                 "raccomandano di preferire il lato in layover: comprime piu' "
                 "diffusori nella stessa cella di range ma li separa lungo l'asse "
                 "di elevazione, che e' esattamente l'asse ricostruito qui."),
    }


# ==========================================================================
# 4.  Lettura dei chip e coregistrazione
# ==========================================================================

@dataclass
class ChipWindow:
    l0: int
    l1: int
    p0: int
    p1: int

    @property
    def n_l(self) -> int:
        return self.l1 - self.l0 + 1

    @property
    def n_p(self) -> int:
        return self.p1 - self.p0 + 1


#: metri per grado di latitudine, e per grado di longitudine all'equatore.
#: Erano due letterali ripetuti in quattordici punti del file; ora stanno qui e
#: le due funzioni sotto sono gli unici modi previsti per usarli.
M_PER_GRADO_LAT = 111_132.0
M_PER_GRADO_LON_EQ = 111_320.0


def _enu_offset(lat, lon, lat0: float, lon0: float):
    """Offset locale (est, nord) in metri di (lat, lon) rispetto a (lat0, lon0).

    Piano tangente con i due fattori metrici standard. Queste due costanti
    comparivano scritte a mano in QUATTORDICI punti del file, ogni volta con un
    contorno diverso: bastava sbagliarne una copia perche' due geometrie si
    spostassero l'una rispetto all'altra senza che nulla desse errore, ed e'
    esattamente il guasto che l'autotest F39 esiste per intercettare. Funziona
    identica su scalari e su array perche' non converte nulla: e' aritmetica."""
    m_lat = M_PER_GRADO_LAT
    m_lon = M_PER_GRADO_LON_EQ * math.cos(math.radians(lat0))
    return (lon - lon0) * m_lon, (lat - lat0) * m_lat


def _latlon_from_enu(east, north, lat0: float, lon0: float):
    """L'inversa esatta di _enu_offset(): (est, nord) -> (lat, lon).

    Stessi due fattori, quindi le due funzioni restano coerenti per
    costruzione; separarle era il modo migliore per farle divergere."""
    m_lat = M_PER_GRADO_LAT
    m_lon = M_PER_GRADO_LON_EQ * math.cos(math.radians(lat0))
    return lat0 + north / m_lat, lon0 + east / m_lon


def pyramid_extent_radar(ann: S1Annotation, geo: Geocoder
                         ) -> Tuple[float, float, float, float]:
    """Estensione in (linea, pixel) occupata dalle piramidi in geometria radar.

    F28: comprende lo spostamento in slant range dovuto alla quota
    (h*cos(theta)/dr, ~47 pixel per l'apice di Cheope). Delimitare l'area sui
    soli centri al suolo lascerebbe fuori proprio le celle in cui cadono gli
    apici."""
    l_min = p_min = float("inf")
    l_max = p_max = float("-inf")
    for p in PYRAMIDS:
        half = p.base_side_m / 2.0
        inv = _local_affine_inverse(geo, p.lat, p.lon, 1.5 * half)
        # i quattro spigoli di base piu' l'apice
        pts = [(-half, -half, p.base_alt_m), (half, -half, p.base_alt_m),
               (half, half, p.base_alt_m), (-half, half, p.base_alt_m),
               (0.0, 0.0, p.base_alt_m + p.height_m)]
        a = math.radians(p.azimuth_deg)
        ca, sa = math.cos(a), math.sin(a)
        for dx, dy, h in pts:
            x = dx * ca - dy * sa
            y = dx * sa + dy * ca
            lat, lon = _latlon_from_enu(x, y, p.lat, p.lon)
            l, q = inv(np.array([lat]), np.array([lon]))
            l, q = float(l[0]), float(q[0])
            href = float(geo.llh(l, q)[2])
            th = math.radians(float(geo.incidence(l, q)))
            q -= (h - href) * math.cos(th) / ann.range_pixel_spacing
            l_min, l_max = min(l_min, l), max(l_max, l)
            p_min, p_max = min(p_min, q), max(p_max, q)
    return l_min, l_max, p_min, p_max


def target_window(ann: S1Annotation, geo: Geocoder, cfg: Config,
                  n_lines: Optional[int] = None) -> Tuple[ChipWindow, int]:
    """Finestra di analisi, confinata al burst.

    F28: per default la finestra e' stretta sulle PIRAMIDI in geometria radar
    piu' un margine di ``cfg.area_margin_m``, non un rettangolo di 1.8 x 2.2 km
    attorno alla piana. Il margine serve a due cose e a nient'altro: dare alla
    coregistrazione e alla calibrazione del piano di riferimento un anello di
    terreno stabile fuori dalle impronte, e permettere una colonna di controllo.
    Tutto il resto della scena non entra nel calcolo.

    Il confinamento al burst e' obbligatorio: il deramping TOPS e la
    decomposizione in sub-aperture sono definiti solo dentro un burst."""
    l_lo, l_hi, p_lo, p_hi = pyramid_extent_radar(ann, geo)

    if cfg.full_scene:                       # comportamento storico, opzionale
        lp = [geo.latlon_to_line_pixel(p.lat, p.lon) for p in PYRAMIDS]
        l_lo, l_hi = min(x[0] for x in lp), max(x[0] for x in lp)
        p_lo, p_hi = min(x[1] for x in lp), max(x[1] for x in lp)
        pad_l = cfg.chip_pad_azim_m / ann.azimuth_pixel_spacing
        pad_p = cfg.chip_pad_range_m / ann.range_pixel_spacing
    else:
        pad_l = cfg.area_margin_m / ann.azimuth_pixel_spacing
        pad_p = cfg.area_margin_m / ann.range_pixel_spacing

    burst = ann.burst_of_line(0.5 * (l_lo + l_hi))
    b_first = burst * ann.lines_per_burst
    b_last = b_first + ann.lines_per_burst - 1

    if n_lines is None:
        l0 = int(math.floor(l_lo - pad_l))
        l1 = int(math.ceil(l_hi + pad_l))
    else:                                   # F07: finestra azimutale estesa
        c = int(round(0.5 * (l_lo + l_hi)))
        l0 = c - n_lines // 2
        l1 = l0 + n_lines - 1

    l0 = max(l0, b_first)
    l1 = min(l1, b_last)
    p0 = max(int(math.floor(p_lo - pad_p)), 0)
    p1 = min(int(math.ceil(p_hi + pad_p)), ann.n_samples - 1)

    # il reticolo multilooked deve avere un numero intero di blocchi
    n_p = p1 - p0 + 1
    p1 -= n_p % max(cfg.look_range, 1)
    n_l = l1 - l0 + 1
    l1 -= n_l % max(cfg.look_azimuth, 1)
    return ChipWindow(l0, l1, p0, p1), burst


@dataclass
class LutRadiometrica:
    """Una tabella di calibrazione o di rumore, sul reticolo sparso dello .xml.

    Sentinel-1 distribuisce entrambe come vettori campionati: alcune decine di
    righe in azimuth, qualche centinaio di colonne in range. Fra i nodi si
    interpola bilinearmente sulle coordinate ASSOLUTE del prodotto, perche' e'
    a quelle che si riferiscono e non al ritaglio."""

    lines: np.ndarray                  # [n_l] righe dei nodi
    pixels: np.ndarray                 # [n_p] colonne dei nodi
    valori: np.ndarray                 # [n_l, n_p]

    def valuta(self, l0: int, p0: int, n_l: int, n_p: int) -> np.ndarray:
        """La tabella interpolata sul ritaglio [l0:l0+n_l, p0:p0+n_p]."""
        ll = np.arange(l0, l0 + n_l, dtype=np.float64)
        pp = np.arange(p0, p0 + n_p, dtype=np.float64)
        # interpolazione separabile: prima in range su ogni riga di nodi, poi
        # in azimuth fra le righe. Fuori dal reticolo si estende il bordo.
        per_riga = np.empty((len(self.lines), n_p), dtype=np.float64)
        for i in range(len(self.lines)):
            per_riga[i] = np.interp(pp, self.pixels, self.valori[i])
        if len(self.lines) == 1:
            return np.repeat(per_riga, n_l, axis=0).astype(np.float32)
        fuori = np.clip(np.searchsorted(self.lines, ll) - 1,
                        0, len(self.lines) - 2)
        l_lo = self.lines[fuori]
        l_hi = self.lines[fuori + 1]
        w = np.clip((ll - l_lo) / np.maximum(l_hi - l_lo, 1e-9), 0.0, 1.0)
        out = (per_riga[fuori] * (1.0 - w[:, None])
               + per_riga[fuori + 1] * w[:, None])
        return out.astype(np.float32)


def _vettori(root: ET.Element, tag_lista: str, tag_pixel: str,
             tag_valori: str) -> Optional[LutRadiometrica]:
    """Raccoglie i vettori omonimi di uno .xml in una LutRadiometrica."""
    nodi = root.findall(f".//{tag_lista}")
    if not nodi:
        return None
    lines: List[float] = []
    pixels: Optional[np.ndarray] = None
    righe: List[np.ndarray] = []
    for v in nodi:
        t_val = v.findtext(tag_valori)
        t_pix = v.findtext(tag_pixel)
        if not t_val or not t_pix:
            continue
        val = np.array([float(x) for x in t_val.split()], dtype=np.float64)
        pix = np.array([float(x) for x in t_pix.split()], dtype=np.float64)
        if pixels is None:
            pixels = pix
        elif len(pix) != len(pixels):
            val = np.interp(pixels, pix, val)
        lines.append(float(_txt(v, "line").split()[0]))
        righe.append(val)
    if pixels is None or not righe:
        return None
    ordine = np.argsort(lines)
    return LutRadiometrica(lines=np.array(lines)[ordine], pixels=pixels,
                           valori=np.array(righe)[ordine])


def leggi_calibrazione(path: str) -> Optional[LutRadiometrica]:
    """La tabella sigmaNought dal calibration-*.xml.

    Il DN complesso dello SLC si converte in ampiezza calibrata dividendo per
    questo fattore: sigma0 = |DN / A_sigma|^2. E' una divisione per un numero
    REALE e positivo, quindi la FASE non cambia: baseline, k_z e quote non ne
    risentono, ne risente solo la radiometria."""
    try:
        root = ET.parse(path).getroot()
    except Exception:                                   # pragma: no cover
        return None
    return _vettori(root, "calibrationVector", "pixel", "sigmaNought")


def leggi_rumore(path: str) -> Optional[LutRadiometrica]:
    """Il rumore termico dal noise-*.xml, in unita' di DN^2.

    Il prodotto porta un profilo in range (noiseRangeVector) e un profilo di
    modulazione in azimuth (noiseAzimuthVector); qui si usa il primo, che e'
    quello che descrive la rampa del NESZ attraverso lo swath, ed e' il termine
    dominante. Serve per sottrarre il pavimento di rumore dall'intensita': in
    VH cross-pol sul deserto il segnale ci sta sopra di pochissimo."""
    try:
        root = ET.parse(path).getroot()
    except Exception:                                   # pragma: no cover
        return None
    lut = _vettori(root, "noiseRangeVector", "pixel", "noiseRangeLut")
    if lut is None:                                     # prodotti IPF vecchi
        lut = _vettori(root, "noiseVector", "pixel", "noiseLut")
    return lut


def read_window(entry: StackEntry, win: ChipWindow, burst: int,
                calibra: bool = False) -> Chip:
    """Legge il ritaglio dello SLC, opzionalmente calibrato in ampiezza.

    Con `calibra` il DN viene diviso per sigmaNought, cosi' che |dato|^2 sia
    sigma0 e non un numero di conteggi arbitrario. La fase resta intatta."""
    import rasterio
    from rasterio.windows import Window

    with rasterio.open(entry.tiff) as ds:
        arr = ds.read(1, window=Window(win.p0, win.l0, win.n_p, win.n_l))
    dato = np.ascontiguousarray(arr.astype(np.complex64))

    if calibra and entry.calibration:
        cal = leggi_calibrazione(entry.calibration)
        if cal is not None:
            a = cal.valuta(win.l0, win.p0, win.n_l, win.n_p)
            dato = (dato / np.maximum(a, 1e-9)).astype(np.complex64)

    return Chip(data=dato, line0=win.l0, pixel0=win.p0,
                date=entry.date, burst=burst)


def mappa_rumore(entry: StackEntry, win: ChipWindow) -> Optional[np.ndarray]:
    """Il pavimento di rumore sul ritaglio, nelle stesse unita' di sigma0.

    Il noise-*.xml da' il rumore in DN^2; per confrontarlo con un'intensita'
    calibrata va diviso per sigmaNought^2, gli stessi fattori usati da
    read_window(calibra=True). Il risultato e' il NESZ pixel per pixel."""
    if not entry.noise or not entry.calibration:
        return None
    noi = leggi_rumore(entry.noise)
    cal = leggi_calibrazione(entry.calibration)
    if noi is None or cal is None:
        return None
    n = noi.valuta(win.l0, win.p0, win.n_l, win.n_p)
    a = cal.valuta(win.l0, win.p0, win.n_l, win.n_p)
    return (n / np.maximum(a, 1e-9) ** 2).astype(np.float32)


def apply_shift(img: np.ndarray, shift: complex) -> np.ndarray:
    """Trasla il contenuto di `img` di +shift (righe, colonne) via rampa di fase.

    Convenzione bloccata dall'autotest: apply_shift(f, d)(x) == f(x - d)."""
    n_l, n_p = img.shape
    fl = np.fft.fftfreq(n_l)[:, None]
    fp = np.fft.fftfreq(n_p)[None, :]
    ramp = np.exp(-2j * np.pi * (fl * shift.real + fp * shift.imag))
    return np.fft.ifft2(np.fft.fft2(img) * ramp).astype(np.complex64)


def estimate_shift(master: np.ndarray, slave: np.ndarray) -> Tuple[complex, float]:
    """Shift globale per cross-correlazione di fase, con peak-to-median ratio.

    Su VH cross-pol nel deserto la decorrelazione a 12 giorni e' forte e il picco
    puo' essere spurio: il PSR permette di riconoscerlo e ricadere sulla
    geolocalizzazione, che qui e' piu' affidabile perche' non dipende dalla
    coerenza radiometrica."""
    m = master - master.mean()
    s = slave - slave.mean()
    cross = np.fft.fft2(m) * np.conj(np.fft.fft2(s))
    mag = np.abs(cross)
    mag[mag == 0] = 1.0
    corr = np.abs(np.fft.ifft2(cross / mag))
    psr = float(corr.max() / max(np.median(corr), 1e-12))
    shift = complex(_batch_subpixel_shift(master[None], slave[None], refine="dft")[0])
    return shift, psr


# ==========================================================================
# 5.  Fase geometrica di riferimento (F02)
# ==========================================================================

def _node_grid(geo: Geocoder, win: ChipWindow, n_nodes: Tuple[int, int]):
    """Reticolo di nodi sul chip e relativi bersagli ECEF.

    Il preambolo era scritto due volte, identico, in reference_range_difference
    e in perp_baseline_field: stesso reticolo, stessa geolocation grid, stessa
    conversione in ECEF. Cambiava solo cosa se ne faceva dopo."""
    nl, np_ = n_nodes
    ln = np.linspace(win.l0, win.l1, nl)
    pn = np.linspace(win.p0, win.p1, np_)
    gl, gp = np.meshgrid(ln, pn, indexing="ij")
    lat, lon, h = geo.llh(gl.ravel(), gp.ravel())
    return ln, pn, ecef_from_llh(lat, lon, h)


def _spline_to_pixels(ln: np.ndarray, pn: np.ndarray, val: np.ndarray,
                      win: ChipWindow, dtype) -> np.ndarray:
    """Spline bicubica sui nodi, valutata su ogni pixel del chip.

    Anche questa coda era duplicata. Le grandezze interpolate (delta_R in metri
    e B_perp in metri) sono geometriche e lisce, quindi la spline non ha nulla
    a che vedere con l'avvolgimento della fase."""
    from scipy.interpolate import RectBivariateSpline

    nl, np_ = len(ln), len(pn)
    spl = RectBivariateSpline(ln, pn, val.reshape(nl, np_),
                              kx=min(3, nl - 1), ky=min(3, np_ - 1))
    ll = np.arange(win.l0, win.l1 + 1, dtype=np.float64)
    pp = np.arange(win.p0, win.p1 + 1, dtype=np.float64)
    return spl(ll, pp).astype(dtype)


def reference_range_difference(
    orb_s: Orbit, orb_m: Orbit, geo: Geocoder, win: ChipWindow,
    n_nodes: Tuple[int, int] = (16, 24),
) -> np.ndarray:
    """delta_R(l,p) = R_slave(l,p) - R_master(l,p) sulla superficie di riferimento.

    E' la fase di terra piatta PIU' la topografia di riferimento, calcolata
    esattamente dalla geometria (orbite + geolocation grid dello .xml) invece che
    stimata con una regressione sulla fase avvolta. Il calcolo esatto e' fatto su
    un reticolo di nodi e interpolato con spline bicubica: delta_R e' una
    quantita' geometrica liscia in metri, quindi l'interpolazione non ha nulla a
    che vedere con l'avvolgimento della fase."""
    ln, pn, tgt = _node_grid(geo, win, n_nodes)
    r_s, _, _ = orb_s.slant_range(tgt)
    r_m, _, _ = orb_m.slant_range(tgt)
    return _spline_to_pixels(ln, pn, r_s - r_m, win, np.float64)


def perp_baseline_field(
    orb_s: Orbit, orb_m: Orbit, geo: Geocoder, win: ChipWindow,
    n_nodes: Tuple[int, int] = (8, 12),
) -> np.ndarray:
    """B_perp(l,p): la baseline ortogonale varia attraverso il chip.

    Usarne un solo valore scalare (come faceva la versione precedente) e' una
    approssimazione accettabile su 3 km, ma calcolarlo per pixel costa poco e
    toglie un termine di errore sistematico dal k_z."""
    ln, pn, tgt = _node_grid(geo, win, n_nodes)

    _, p_m, v_m = orb_m.slant_range(tgt)
    _, p_s, _ = orb_s.slant_range(tgt)

    los = tgt - p_m
    los /= np.linalg.norm(los, axis=1, keepdims=True)
    vh = v_m / np.linalg.norm(v_m, axis=1, keepdims=True)
    nrm = np.cross(los, vh)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)

    b = p_s - p_m
    b_par = np.einsum("ij,ij->i", b, los)
    b_perp = np.einsum("ij,ij->i", b - b_par[:, None] * los, nrm)

    return _spline_to_pixels(ln, pn, b_perp, win, np.float32)


# ==========================================================================
# 6.  Cubo interferometrico
# ==========================================================================

@dataclass
class InterfCube:
    y: np.ndarray                  # complex64 [n_date, n_l, n_p]  (interferogrammi)
    amp_master: np.ndarray         # float32 [n_l, n_p]
    k_z: np.ndarray                # float32 [n_date, n_l, n_p]  rad/m
    dates: List[str]
    baselines: List[Baseline]
    win: ChipWindow
    incidence: np.ndarray          # float32 [n_l, n_p] gradi
    ann: S1Annotation
    geo: Geocoder
    master_date: str
    coreg: List[Dict[str, Any]] = field(default_factory=list)
    flat_mask: Optional[np.ndarray] = None
    #: NESZ del master sul ritaglio, in unita' di sigma0 (None se il prodotto
    #: non porta i .xml di calibrazione e rumore)
    nesz: Optional[np.ndarray] = None
    #: True se le ampiezze sono sigma0 calibrato e non conteggi DN
    calibrato: bool = False


def _pyramid_footprint_mask(geo: Geocoder, win: ChipWindow, scale: float = 1.0) -> np.ndarray:
    """Maschera booleana dei pixel che cadono sull'impronta delle piramidi."""
    ll = np.arange(win.l0, win.l1 + 1, dtype=np.float64)
    pp = np.arange(win.p0, win.p1 + 1, dtype=np.float64)
    gl, gp = np.meshgrid(ll, pp, indexing="ij")
    lat, lon, _ = geo.llh(gl.ravel(), gp.ravel())
    lat = lat.reshape(gl.shape)
    lon = lon.reshape(gl.shape)

    m = np.zeros(gl.shape, dtype=bool)
    for p in PYRAMIDS:
        dx, dy = _enu_offset(lat, lon, p.lat, p.lon)
        half = 0.5 * p.base_side_m * scale
        m |= (np.abs(dx) <= half) & (np.abs(dy) <= half)
    return m


def build_interf_cube(
    entries: Sequence[StackEntry], cfg: Config, verbose: bool = True
) -> InterfCube:
    """Legge, dirampa, coregistra e riferisce interferometricamente lo stack VH.

    Il cubo che esce contiene y_i = s_i * conj(s_m) / |s_m|, cioe' interferogrammi
    normalizzati (F01), con la fase geometrica di riferimento gia' rimossa (F02)
    e la costante e la rampa residue calibrate sulla piana (fuori dalle impronte
    delle piramidi). Solo a quel punto la fase di y_i vale k_z * h e l'inversione
    tomografica ha senso."""
    anns = [parse_annotation(e.annotation) for e in entries]
    orbits = [read_orbit(e.annotation) for e in entries]
    geos = [Geocoder(a) for a in anns]

    tgt0 = ecef_from_llh(PYRAMIDS[0].lat, PYRAMIDS[0].lon, PYRAMIDS[0].base_alt_m)
    bl_probe = compute_baselines(list(zip([e.date for e in entries], orbits)), tgt0, 0)
    mi = pick_supermaster(bl_probe)
    baselines = compute_baselines(list(zip([e.date for e in entries], orbits)), tgt0, mi)
    if verbose:
        print(f"  supermaster selezionato: {entries[mi].date} "
              f"(F13: minimo costo di decorrelazione combinato)")

    ann_m, orb_m, geo_m = anns[mi], orbits[mi], geos[mi]
    win, burst = target_window(ann_m, geo_m, cfg)
    n_l, n_p = win.n_l, win.n_p

    # --- F42: finestra dello slave presa dove il bersaglio sta DAVVERO -----
    # Ogni prodotto ha il proprio inquadramento lungo l'orbita. Sulla stessa
    # traccia relativa, il bersaglio cade a linee molto diverse a seconda della
    # missione: fra S1C e S1A su Giza la differenza misurata e' di ~4510 linee
    # (63 km di volo) e fino a un centinaio di pixel in range. La versione
    # precedente leggeva TUTTE le date alla stessa finestra assoluta del master
    # e correggeva dopo con apply_shift(), che e' una traslazione CIRCOLARE via
    # rampa di fase: valida per frazioni di pixel, non per migliaia di linee.
    # Con un offset del genere lo slave era un pezzo di deserto a 63 km di
    # distanza, ripiegato su se stesso, e per giunta dirampato con i parametri
    # del burst sbagliato. L'offset intero ora entra nella finestra di
    # lettura (e nell'indice di burst dello slave, che il deramping TOPS usa);
    # alla rampa di fase resta solo il residuo sub-pixel.
    l_m, p_m = geo_m.latlon_to_line_pixel(PYRAMIDS[0].lat, PYRAMIDS[0].lon)
    finestre: List[Tuple[ChipWindow, int, complex]] = []
    tenute: List[int] = []
    fuori: List[str] = []
    for i in range(len(entries)):
        if i == mi:
            finestre.append((win, burst, 0j))
            tenute.append(i)
            continue
        l_s, p_s = geos[i].latlon_to_line_pixel(PYRAMIDS[0].lat, PYRAMIDS[0].lon)
        d_l, d_p = l_s - l_m, p_s - p_m
        off_l, off_p = int(round(d_l)), int(round(d_p))
        w_s = ChipWindow(win.l0 + off_l, win.l1 + off_l,
                         win.p0 + off_p, win.p1 + off_p)
        b_s = anns[i].burst_of_line(0.5 * (w_s.l0 + w_s.l1))
        b_first = b_s * anns[i].lines_per_burst
        b_last = b_first + anns[i].lines_per_burst - 1
        if (w_s.l0 < b_first or w_s.l1 > b_last
                or w_s.p0 < 0 or w_s.p1 > anns[i].n_samples - 1):
            # il ritaglio scavalcherebbe il bordo del burst o del prodotto:
            # il deramping TOPS non e' definito a cavallo di due burst.
            fuori.append(f"{entries[i].date} (offset {off_l:+d} linee, "
                         f"{off_p:+d} pixel)")
            continue
        finestre.append((w_s, b_s, complex(d_l - off_l, d_p - off_p)))
        tenute.append(i)

    if len(tenute) < len(entries):
        if verbose:
            print(f"  F42: {len(entries) - len(tenute)} acquisizioni scartate, il "
                  f"ritaglio cade fuori dal burst: " + ", ".join(fuori))
        entries = [entries[i] for i in tenute]
        anns = [anns[i] for i in tenute]
        orbits = [orbits[i] for i in tenute]
        geos = [geos[i] for i in tenute]
        mi = tenute.index(mi)
        ann_m, orb_m, geo_m = anns[mi], orbits[mi], geos[mi]
        baselines = compute_baselines(
            list(zip([e.date for e in entries], orbits)), tgt0, mi)
    if verbose:
        dl_max = max(abs(w.l0 - win.l0) for w, _b, _s in finestre)
        dp_max = max(abs(w.p0 - win.p0) for w, _b, _s in finestre)
        print(f"  F42: ritaglio riposizionato per data (offset massimo "
              f"{dl_max} linee, {dp_max} pixel); alla rampa di fase resta "
              f"solo il residuo sub-pixel")

    # Rete di sicurezza sulla baseline critica: oltre di essa lo spostamento
    # spettrale in range supera l'intera banda del chirp e fra le due
    # acquisizioni non resta nessuna frequenza in comune. La coerenza e' zero
    # per costruzione, non per decorrelazione temporale.
    r0 = 0.5 * C_LIGHT * ann_m.slant_range_time
    th0 = math.radians(float(np.median(geo_m.incidence(
        np.array([0.5 * (win.l0 + win.l1)]), np.array([0.5 * (win.p0 + win.p1)])))))
    b_crit = (ann_m.range_bandwidth * r0 * math.tan(th0)
              * ann_m.wavelength / C_LIGHT)
    oltre = [b for b in baselines if abs(b.b_perp) > b_crit]
    if oltre and verbose:
        print(f"  ATTENZIONE: {len(oltre)} acquisizioni oltre la baseline "
              f"critica ({b_crit / 1000.0:.1f} km): "
              + ", ".join(f"{b.date} ({b.b_perp / 1000.0:+.0f} km)" for b in oltre)
              + " -- non sono interferometriche.")

    calibra = cfg.calibrazione and bool(entries[mi].calibration)
    chip_m = read_window(entries[mi], win, burst, calibra=calibra)
    img_m = tops_deramp(chip_m, ann_m)
    amp_m = np.abs(img_m).astype(np.float32)

    # --- calibrazione radiometrica e rumore termico (F40) -----------------
    # In VH cross-pol sul deserto il segnale sta pochi dB sopra il NESZ: una
    # parte dell'ampiezza della piana NON e' retrodiffusione ma rumore
    # termico. Sottrarlo in INTENSITA' (non in ampiezza) e' l'unico modo
    # corretto, perche' le potenze si sommano e le ampiezze no.
    # Il rumore NON va sottratto qui: amp_m normalizza gli interferogrammi
    # (yi /= amp_m) e azzerare le celle sotto il pavimento produrrebbe
    # divisioni per zero che distruggono la stima di coerenza -- misurato:
    # la coerenza mediana crollava da 0.281 a 0.083. La sottrazione si fa a
    # valle del multilooking, dove l'intensita' media e' stabile.
    nesz = mappa_rumore(entries[mi], win) if calibra else None
    if nesz is not None:
        if verbose:
            def _db(x: float) -> float:
                return 10.0 * math.log10(max(x, 1e-30))
            n_med = float(np.median(nesz))
            s_med = float(np.median(amp_m.astype(np.float64) ** 2))
            print(f"      radiometria (F40): sigma0 calibrato dal "
                  f"calibration-*.xml; NESZ mediano {_db(n_med):.1f} dB, "
                  f"segnale {_db(s_med):.1f} dB, SNR {_db(s_med / max(n_med, 1e-30)):.1f} dB")
            print(f"      celle entro 3 dB dal rumore: "
                  f"{100.0 * float(np.mean(amp_m.astype(np.float64) ** 2 < 2.0 * nesz)):.1f} "
                  f"% -- il pavimento viene sottratto dopo il multilooking")

    ll = np.arange(win.l0, win.l1 + 1, dtype=np.float64)
    pp = np.arange(win.p0, win.p1 + 1, dtype=np.float64)
    gl, gp = np.meshgrid(ll, pp, indexing="ij")
    inc = geo_m.incidence(gl, gp).astype(np.float32)
    sin_t = np.sin(np.radians(inc)).astype(np.float32)

    # slant range per pixel del master, dalla geometria del prodotto
    r_pix = (ann_m.slant_range_near
             + (pp - 0.0) * ann_m.range_pixel_spacing)[None, :].astype(np.float32)
    r_pix = np.broadcast_to(r_pix, (n_l, n_p)).astype(np.float32)

    # la piana di riferimento: tutto tranne le impronte delle piramidi allargate
    flat_mask = ~_pyramid_footprint_mask(geo_m, win, scale=1.6)

    lam = ann_m.wavelength
    y = np.zeros((len(entries), n_l, n_p), dtype=np.complex64)
    k_z = np.zeros((len(entries), n_l, n_p), dtype=np.float32)
    coreg: List[Dict[str, Any]] = []

    for i, entry in enumerate(entries):
        ann, orb, geo = anns[i], orbits[i], geos[i]
        win_i, burst_i, geo_shift = finestre[i]

        if i == mi:
            img = img_m
            shift, src, psr = 0j, "master", float("inf")
        else:
            # F42: la finestra e' gia' quella dello slave, quindi il chip
            # contiene lo stesso pezzo di terreno del master e il residuo da
            # correggere e' sub-pixel.
            chip = read_window(entry, win_i, burst_i, calibra=calibra)
            img = tops_deramp(chip, ann)

            corr_shift, psr = estimate_shift(img_m, img)
            if abs(corr_shift - geo_shift) <= 2.0 and psr >= 8.0:
                shift, src = corr_shift, "corr"
            else:
                shift, src = geo_shift, "geo"
            img = apply_shift(img, -shift)

        # --- F02: fase geometrica di riferimento -------------------------
        d_r = reference_range_difference(orb, orb_m, geo_m, win)
        phi_ref = (-4.0 * np.pi / lam) * d_r
        yi = (img * np.conj(img_m)).astype(np.complex64)
        yi *= np.exp(-1j * phi_ref).astype(np.complex64)

        # --- normalizzazione: |y| = |s_i|, fase = fase differenziale ------
        yi /= np.maximum(amp_m, 1e-6)

        # --- F04: k_z per pixel ------------------------------------------
        bpf = perp_baseline_field(orb, orb_m, geo_m, win)
        k_z[i] = (4.0 * np.pi * bpf / (lam * r_pix * sin_t)).astype(np.float32)

        y[i] = yi
        coreg.append({
            "data": entry.date,
            "b_perp_m": round(baselines[i].b_perp, 2),
            "b_temp_giorni": round(baselines[i].b_temp, 1),
            "shift_px": [round(shift.real, 3), round(shift.imag, 3)],
            "offset_finestra_px": [win_i.l0 - win.l0, win_i.p0 - win.p0],
            "burst": burst_i,
            "sorgente_shift": src,
            "psr": None if math.isinf(psr) else round(psr, 2),
        })
        if verbose:
            print(f"    [{entry.date}] B_perp={baselines[i].b_perp:+8.2f} m  "
                  f"dt={baselines[i].b_temp:+6.0f} d  "
                  f"finestra=({win_i.l0 - win.l0:+6d},{win_i.p0 - win.p0:+5d}) px "
                  f"burst {burst_i}  "
                  f"residuo=({shift.real:+6.3f},{shift.imag:+6.3f}) px [{src}]"
                  + ("" if math.isinf(psr) else f"  PSR={psr:5.1f}"))

    cube = InterfCube(
        y=y, amp_master=amp_m, k_z=k_z, nesz=nesz, calibrato=calibra,
        dates=[e.date for e in entries], baselines=baselines, win=win,
        incidence=inc, ann=ann_m, geo=geo_m, master_date=entries[mi].date,
        coreg=coreg, flat_mask=flat_mask,
    )
    calibrate_reference_plane(cube, verbose=verbose)
    return cube


def calibrate_reference_plane(cube: InterfCube, verbose: bool = True) -> None:
    """Azzera costante e rampa residue di ciascuna data SULLA PIANA.

    Restano sempre, dopo la rimozione geometrica, un offset di fase per data
    (atmosfera, errore d'orbita, fase del diffusore) e una rampa lentissima
    (errore di baseline). Si stimano sui soli pixel della piana -- esclusa
    l'impronta allargata delle piramidi -- cosi' che la stima non possa
    assorbire il rilievo che vogliamo misurare. La rampa e' cercata come picco
    spettrale del complesso, mai per unwrap della fase avvolta."""
    n_d, n_l, n_p = cube.y.shape
    mask = cube.flat_mask
    if mask is None or mask.sum() < 100:
        return

    yy, xx = np.mgrid[0:n_l, 0:n_p].astype(np.float32)
    for i in range(n_d):
        if cube.dates[i] == cube.master_date:
            continue
        yi = cube.y[i]

        # rampa: picco della FFT2 del complesso sui soli pixel di piana
        z = np.where(mask, yi, 0)
        sp = np.fft.fft2(z)
        k = int(np.argmax(np.abs(sp)))
        kl, kp = k // n_p, k % n_p
        fl = np.fft.fftfreq(n_l)[kl]
        fp = np.fft.fftfreq(n_p)[kp]
        yi = yi * np.exp(-2j * np.pi * (fl * yy + fp * xx)).astype(np.complex64)

        # costante: fase media complessa sulla piana
        c = np.mean(yi[mask])
        if abs(c) > 1e-12:
            yi = yi * np.exp(-1j * np.angle(c)).astype(np.complex64)
        cube.y[i] = yi.astype(np.complex64)

    if verbose:
        print("  piano di riferimento calibrato su "
              f"{int(mask.sum())} pixel di piana (F02)")


# ==========================================================================
# 7.  Multilooking, inversione, superficie
# ==========================================================================

def multilook(a: np.ndarray, sl: int, sp: int) -> np.ndarray:
    """Media a blocchi sugli ultimi due assi spaziali (assi -2, -1 di un 2D,
    oppure assi 1 e 2 di un cubo [n, l, p])."""
    if sl == 1 and sp == 1:
        return a
    if a.ndim == 2:
        n_l, n_p = a.shape
        n_l2, n_p2 = n_l // sl, n_p // sp
        return a[:n_l2 * sl, :n_p2 * sp].reshape(n_l2, sl, n_p2, sp).mean(axis=(1, 3))
    n_d, n_l, n_p = a.shape
    n_l2, n_p2 = n_l // sl, n_p // sp
    return a[:, :n_l2 * sl, :n_p2 * sp].reshape(n_d, n_l2, sl, n_p2, sp).mean(axis=(2, 4))


def tomographic_periodogram(
    y: np.ndarray, k_z: np.ndarray, z_axis: np.ndarray, sign: float = 1.0,
) -> np.ndarray:
    """h(z) = sum_i y_i * exp(-j*s*k_z_i*z)  --  beamforming a filtro adattato.

    E' l'inversione h = A^H Y delle fonti, con la differenza che k_z varia per
    pixel (F04): il prodotto va quindi fatto pixel per pixel invece che con una
    singola matrice di steering condivisa."""
    n_d, n_l, n_p = y.shape
    out = np.zeros((n_l, n_p, len(z_axis)), dtype=np.complex64)
    for zi, z in enumerate(z_axis):
        acc = np.zeros((n_l, n_p), dtype=np.complex64)
        for i in range(n_d):
            acc += y[i] * np.exp(-1j * sign * k_z[i] * z)
        out[:, :, zi] = acc
    return out


def surface_from_tomogram(
    tomo: np.ndarray, z_axis: np.ndarray, y_abs_sum: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Quota del diffusore dominante e coerenza di fit, pixel per pixel.

    quota  = argmax_z |h(z)|, raffinato con parabola a tre punti sul picco
    gamma  = |h(z_pic)| / sum_i |y_i|   in [0, 1]

    gamma e' la frazione di energia che si somma in fase alla quota stimata:
    e' la misura di quanto la stima e' credibile, e diventa la maschera di
    qualita' della superficie."""
    mag = np.abs(tomo)
    k = np.argmax(mag, axis=2)
    n_l, n_p, n_z = mag.shape
    ii, jj = np.meshgrid(np.arange(n_l), np.arange(n_p), indexing="ij")

    km1 = np.clip(k - 1, 0, n_z - 1)
    kp1 = np.clip(k + 1, 0, n_z - 1)
    c0 = mag[ii, jj, k]
    cm = mag[ii, jj, km1]
    cp = mag[ii, jj, kp1]
    den = cm - 2.0 * c0 + cp
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)
    frac = np.clip(0.5 * (cm - cp) / den, -1.0, 1.0)
    frac = np.where((k == 0) | (k == n_z - 1), 0.0, frac)

    dz = float(z_axis[1] - z_axis[0])
    height = (z_axis[k] + frac * dz).astype(np.float32)
    gamma = (c0 / np.maximum(y_abs_sum, 1e-9)).astype(np.float32)
    return height, np.clip(gamma, 0.0, 1.0)


def null_gamma_threshold(
    k_z_typ: np.ndarray, z_axis: np.ndarray, q: float = 99.0,
    trials: int = 20000, seed: int = 20260828,
) -> Tuple[float, Dict[str, float]]:
    """F22 -- soglia di qualita' calibrata sulla distribuzione NULLA.

    Con N baseline e fasi puramente casuali, il picco del periodogramma
    normalizzato non vale zero: e' il massimo di N fasori casuali su tutto
    l'asse z, e con N = 11 ha gia' mediana ~0.56. Una soglia scelta a mano piu'
    bassa di quel valore accetta rumore puro e fa sembrare misurato cio' che non
    lo e'. La soglia corretta e' un percentile alto di questa distribuzione:
    e' un test di ipotesi, non una preferenza estetica."""
    rng = np.random.default_rng(seed)
    n = len(k_z_typ)
    a = np.exp(-1j * np.outer(k_z_typ, z_axis))               # [N, n_z]
    ph = rng.uniform(0.0, 2.0 * np.pi, (trials, n))
    peaks = np.abs(np.exp(1j * ph) @ a).max(axis=1) / n
    thr = float(np.percentile(peaks, q))
    return thr, {
        "mediana_nulla": round(float(np.median(peaks)), 4),
        "p90_nulla": round(float(np.percentile(peaks, 90)), 4),
        "p99_nulla": round(float(np.percentile(peaks, 99)), 4),
        "p99_9_nulla": round(float(np.percentile(peaks, 99.9)), 4),
        "soglia_usata": round(thr, 4),
        "percentile_usato": q,
        "prove": trials,
    }


def _local_affine_inverse(geo: Geocoder, lat0: float, lon0: float, span_m: float):
    """Mappa affine locale (lat, lon) -> (linea, pixel), calibrata su 9 ancore.

    L'inversione esatta e' un Newton sulla spline: corretta ma troppo lenta per
    decine di migliaia di punti. Su poche centinaia di metri la mappa e' affine
    con errore molto inferiore al pixel, e le ancore la vincolano ai minimi
    quadrati."""
    d = np.linspace(-span_m, span_m, 3)
    a_lat, a_lon, a_l, a_p = [], [], [], []
    for dy in d:
        for dx in d:
            la, lo = _latlon_from_enu(dx, dy, lat0, lon0)
            l, q = geo.latlon_to_line_pixel(la, lo)
            a_lat.append(la)
            a_lon.append(lo)
            a_l.append(l)
            a_p.append(q)
    a = np.column_stack([np.array(a_lat) - lat0, np.array(a_lon) - lon0,
                         np.ones(len(a_lat))])
    cl, *_ = np.linalg.lstsq(a, np.array(a_l), rcond=None)
    cp, *_ = np.linalg.lstsq(a, np.array(a_p), rcond=None)

    def inv(lat, lon):
        lat = np.asarray(lat, dtype=np.float64)
        lon = np.asarray(lon, dtype=np.float64)
        dl = lat - lat0
        do = lon - lon0
        return (cl[0] * dl + cl[1] * do + cl[2],
                cp[0] * dl + cp[1] * do + cp[2])
    return inv


def simulate_pyramids_radar(
    geo: Geocoder, ann: S1Annotation, win: ChipWindow,
    incidence: np.ndarray, n_side: int = 220,
) -> Dict[str, Any]:
    """F24 -- proiezione in avanti delle piramidi nella GEOMETRIA RADAR.

    Un punto a quota h sopra la superficie di riferimento ha slant range piu'
    corto di h*cos(theta), cioe' si sposta verso il near range di
    h*cos(theta)/spaziatura_range pixel. Con h = 138 m, theta = 37 gradi e
    spaziatura 2.33 m sono ~47 pixel: usare l'impronta AL SUOLO come maschera
    delle piramidi guarda pixel che non contengono la piramide, e mette la
    calibrazione del segno e la validazione fuori bersaglio.

    Restituisce, sul reticolo del chip:
      sim_h      quota massima (sopra il riferimento) che cade in ciascuna cella
      sim_mask   celle raggiunte da almeno un punto di superficie
      sim_fold   quanti punti di superficie distinti vi si ripiegano (layover)
      sim_name   indice della piramide dominante (-1 se nessuna)"""
    n_l, n_p = win.n_l, win.n_p
    sim_h = np.full((n_l, n_p), -1e9, dtype=np.float32)
    sim_fold = np.zeros((n_l, n_p), dtype=np.int32)
    sim_name = np.full((n_l, n_p), -1, dtype=np.int8)

    cos_t = np.cos(np.radians(np.median(incidence)))
    dr = ann.range_pixel_spacing

    for pi, p in enumerate(PYRAMIDS):
        # reticolo sulla base, quota della superficie piramidale (piramide retta)
        u = np.linspace(-1.0, 1.0, n_side)
        uu, vv = np.meshgrid(u, u, indexing="ij")
        frac = 1.0 - np.maximum(np.abs(uu), np.abs(vv))       # 1 all'apice, 0 alla base
        h_srf = p.base_alt_m + p.height_m * frac
        half = p.base_side_m / 2.0
        a = math.radians(p.azimuth_deg)
        ca, sa = math.cos(a), math.sin(a)
        dx = uu * half * ca - vv * half * sa
        dy = uu * half * sa + vv * half * ca

        lat, lon = _latlon_from_enu(dx, dy, p.lat, p.lon)

        # (lat, lon) -> (linea, pixel) alla quota di RIFERIMENTO.
        # L'inversione di Newton costa troppo per decine di migliaia di punti:
        # su un'area di poche centinaia di metri la mappa e' affine con errore
        # ben sotto il pixel, quindi la si calibra su 9 ancore e si applica
        # vettorialmente.
        inv = _local_affine_inverse(geo, p.lat, p.lon, 1.5 * half)
        ll, pp = inv(lat.ravel(), lon.ravel())
        href = geo.llh(ll, pp)[2]

        dh = h_srf.ravel() - href
        px = pp - dh * cos_t / dr
        li = np.round(ll - win.l0).astype(int)
        pj = np.round(px - win.p0).astype(int)

        ok = (li >= 0) & (li < n_l) & (pj >= 0) & (pj < n_p)
        li, pj, dh = li[ok], pj[ok], dh[ok]
        np.add.at(sim_fold, (li, pj), 1)
        # quota massima per cella
        order = np.argsort(dh)
        sim_h[li[order], pj[order]] = np.maximum(
            sim_h[li[order], pj[order]], dh[order].astype(np.float32))
        sim_name[li, pj] = pi

    mask = sim_h > -1e8
    sim_h = np.where(mask, sim_h, 0.0).astype(np.float32)
    per = [{"nome": PYRAMIDS[k].name,
            "mask": (sim_name == k) & mask,
            "h": np.where((sim_name == k) & mask, sim_h, 0.0).astype(np.float32)}
           for k in range(len(PYRAMIDS))]
    return {"sim_h": sim_h, "sim_mask": mask, "sim_fold": sim_fold,
            "sim_name": sim_name, "per": per}


def despike_surface(
    h: np.ndarray, good: np.ndarray, sigma_h: float,
    window: int = 5, n_sigma: float = 3.0,
) -> Tuple[np.ndarray, int]:
    """F27 -- rimozione dichiarata degli errori grossolani di quota.

    Anche sopra la soglia di qualita' resta una coda di pixel il cui
    periodogramma ha agganciato un lobo laterale invece del lobo principale:
    sono errori grossolani, non rumore gaussiano, e valgono centinaia di metri.
    Se restano dentro dominano la scala del rendering e rendono illeggibile una
    superficie che per il 98 per cento sta in poche decine di metri.

    Il criterio e' esplicito: un nodo viene sostituito con la mediana locale se
    se ne discosta di piu' di n_sigma volte la precisione teorica sigma_h, o se
    non supera la soglia di qualita'. Il numero di sostituzioni viene riportato:
    la quota grezza resta nei .npy, quella ripulita va solo nel rendering."""
    from scipy.ndimage import median_filter

    med = median_filter(h.astype(np.float32), size=window, mode="nearest")
    bad = (~good) | (np.abs(h - med) > n_sigma * max(sigma_h, 1e-6))
    return np.where(bad, med, h).astype(np.float32), int((good & bad).sum())


def calibrate_kz_sign(
    y: np.ndarray, k_z: np.ndarray, z_axis: np.ndarray,
    pyr_mask: np.ndarray, y_abs_sum: np.ndarray, verbose: bool = True,
) -> Tuple[float, Dict[str, Any]]:
    """F05 -- il segno di k_z e' deciso dai dati, non assunto.

    La convenzione di segno dipende da come sono orientati LOS e normale, e un
    segno sbagliato ribalta la topografia. Il test e' la validazione di livello 1
    (ch15): le piramidi devono risultare SOPRA la piana. Si prova entrambi i
    segni e si riporta il margine, cosi' la decisione resta ispezionabile."""
    res: Dict[str, Any] = {}
    best, best_margin = 1.0, -np.inf
    for s in (1.0, -1.0):
        tomo = tomographic_periodogram(y, k_z, z_axis, sign=s)
        h, g = surface_from_tomogram(tomo, z_axis, y_abs_sum)
        w = g > np.percentile(g, 60)
        on = pyr_mask & w
        off = (~pyr_mask) & w
        if on.sum() < 20 or off.sum() < 20:
            margin = -np.inf
            h_on = h_off = float("nan")
        else:
            h_on = float(np.median(h[on]))
            h_off = float(np.median(h[off]))
            margin = h_on - h_off
        res[f"segno_{'+1' if s > 0 else '-1'}"] = {
            "quota_mediana_piramidi_m": round(h_on, 2) if np.isfinite(h_on) else None,
            "quota_mediana_piana_m": round(h_off, 2) if np.isfinite(h_off) else None,
            "margine_m": round(margin, 2) if np.isfinite(margin) else None,
            "pixel_su_piramide": int(on.sum()),
            "pixel_su_piana": int(off.sum()),
        }
        if margin > best_margin:
            best, best_margin = s, margin
    res["segno_scelto"] = int(best)
    res["criterio"] = ("le piramidi devono risultare sopra la piana; si sceglie il "
                       "segno che massimizza la mediana(quota sulle impronte) - "
                       "mediana(quota sulla piana)")
    res["margine_scelto_m"] = round(float(best_margin), 2) if np.isfinite(best_margin) else None
    if verbose:
        print(f"  segno di k_z calibrato: {int(best):+d}  "
              f"(margine piramidi-piana = {best_margin:+.1f} m)")
    return best, res


def stack_coherence(y: np.ndarray, amp_m: np.ndarray, window: int = 7,
                    master_date: str = "", dates: Optional[Sequence[str]] = None
                    ) -> np.ndarray:
    """F23 -- coerenza interferometrica di ciascuna data contro il master.

    Attenzione allo stimatore: il cubo contiene y_i = s_i*conj(s_m)/|s_m|, quindi
    <y_i> NON e' <s_i conj(s_m)> e usarlo direttamente da uno stimatore diverso
    e distorto. Ricostruendo il numeratore vero:

        <s_i conj(s_m)>  =  < y_i * |s_m| >
        <|s_i|^2>        =  < |y_i|^2 >
        <|s_m|^2>        =  < |s_m|^2 >

    e la coerenza e' |<s_i conj(s_m)>| / sqrt(<|s_i|^2> <|s_m|^2>)."""
    from scipy.ndimage import uniform_filter

    n_d = y.shape[0]
    a_m = amp_m.astype(np.float32)
    den_m = uniform_filter(a_m ** 2, window)
    acc = np.zeros(y.shape[1:], dtype=np.float32)
    cnt = 0
    for i in range(n_d):
        if dates is not None and master_date and dates[i] == master_date:
            continue
        yi = y[i] * a_m                                   # = s_i * conj(s_m)
        num = uniform_filter(yi.real.astype(np.float32), window) \
            + 1j * uniform_filter(yi.imag.astype(np.float32), window)
        den_i = uniform_filter((np.abs(y[i]) ** 2).astype(np.float32), window)
        acc += (np.abs(num) / np.maximum(np.sqrt(den_i * den_m), 1e-9)).astype(np.float32)
        cnt += 1
    return (acc / max(cnt, 1)).astype(np.float32)


# ==========================================================================
# 8.  Micro-moto Doppler  (F06 - F10)
# ==========================================================================

@dataclass
class MMPlan:
    b_cd: float
    b_dl: float
    b_sub: float
    b_shift: float
    step: float
    n_d: int
    df_bin: float
    bins_per_sub: float
    k_a: float
    k_t: float
    t_illum: float
    dt_step: float
    dt_shift: float
    t_window: float
    f_sample: float
    f_max_obs: float
    f_min_obs: float
    aperture_m: float
    lambda_sound_m: float
    delta_z_acoustic: float

    @property
    def usable_band(self) -> bool:
        """La marcia ha senso solo se la finestra osservata contiene almeno un
        ciclo della frequenza selezionata da B_shift. Vale per
        B_shift < B_DL / 3."""
        return self.f_max_obs >= self.f_min_obs

    def as_text(self) -> str:
        if self.f_max_obs > 2.0 * self.f_min_obs:
            verdict = "una banda di frequenze meccaniche"
        elif self.usable_band:
            verdict = "poco piu' di un bin di frequenza"
        else:
            verdict = ("NULLA: la finestra osservata non contiene un ciclo intero "
                       "della frequenza selezionata (serve B_shift < B_DL/3)")
        return "\n".join([
            "-" * 76,
            "BANCO DI SUB-APERTURE DOPPLER  (micro-moto)",
            "-" * 76,
            f"  banda Doppler totale        B_cD    = {self.b_cd:9.2f} Hz",
            f"  banda di guardia            B_DL    = {self.b_dl:9.2f} Hz",
            f"  larghezza sub-banda         B_sub   = {self.b_sub:9.2f} Hz",
            f"  separazione master/slave    B_shift = {self.b_shift:9.2f} Hz",
            f"  passo di marcia                     = {self.step:9.2f} Hz",
            f"  sub-aperture                N_D     = {self.n_d:9d}",
            f"  risoluzione griglia Doppler         = {self.df_bin:9.3f} Hz/bin",
            f"  bin per sub-banda                   = {self.bins_per_sub:9.1f}",
            "",
            "  MAPPATURA FREQUENZA -> TEMPO (dal Doppler rate dell'annotation)",
            f"  azimuth FM rate             k_a     = {self.k_a:9.1f} Hz/s",
            f"  rate TOPS combinata         k_t     = {self.k_t:9.1f} Hz/s",
            f"  tempo di illuminazione del target   = {self.t_illum * 1e3:9.1f} ms",
            f"  passo temporale fra sub-aperture    = {self.dt_step * 1e3:9.2f} ms",
            f"  ritardo master-slave                = {self.dt_shift * 1e3:9.2f} ms",
            f"  finestra osservata (marcia completa)= {self.t_window * 1e3:9.2f} ms",
            "",
            f"  campionamento meccanico     f_s     = {self.f_sample:9.2f} Hz",
            f"  frequenza osservabile max (Nyquist) = {self.f_max_obs:9.2f} Hz",
            f"  frequenza osservabile min (finestra)= {self.f_min_obs:9.2f} Hz",
            f"  >> il banco risolve {verdict}",
            "",
            "  PERCHE' LA TOMOGRAFIA ACUSTICA DELLE FONTI QUI NON E' APPLICABILE",
            f"  apertura sintetica          L_sa    = {self.aperture_m:9.0f} m",
            f"  lambda acustica a {self.f_max_obs:5.1f} Hz (v=6000 m/s) = {self.lambda_sound_m:9.0f} m",
            f"  delta_z acustica = lambda*R/(2*L_sa) = {self.delta_z_acoustic / 1000:8.1f} km",
            "",
            "  Limite dei dati: il TOPS di Sentinel-1 illumina un bersaglio per",
            f"  {self.t_illum * 1e3:.0f} ms soltanto, quindi la traccia di vibrazione dura "
            f"{self.t_window * 1e3:.0f} ms",
            "  e la frequenza meccanica accessibile e' di decine di Hz. Le fonti",
            "  lavorano su spotlight con B_cD ~ 22 kHz e indagano a 12.5 kHz: due",
            "  ordini di grandezza in piu'. E' questa la ragione per cui qui la",
            "  profondita' viene dalle baseline orbitali e non dal micro-moto, che",
            "  resta un attributo di superficie e nulla di piu'.",
            "-" * 76,
        ])


def azimuth_rates(ann: S1Annotation, burst_idx: int, pixel: float) -> Tuple[float, float]:
    """(k_a, k_t): azimuth FM rate e rate TOPS combinata al pixel dato.

    Sono gli stessi coefficienti che il deramping TOPS usa. k_a mappa la
    frequenza Doppler sul tempo di osservazione del bersaglio; k_t e' la rampa
    che il deramping rimuove."""
    b = ann.bursts[min(burst_idx, len(ann.bursts) - 1)]
    tau = ann.slant_range_time + pixel / ann.range_sampling_rate
    k_a = sum(c * (tau - b.fm_rate_t0) ** i for i, c in enumerate(b.fm_rate_poly))
    k_psi = math.radians(ann.azimuth_steering_rate_deg)
    k_s = 2.0 * ann.orbit_velocity * k_psi / ann.wavelength
    den = k_a - k_s
    k_t = k_a * k_s / (den if abs(den) > 1e-6 else 1e-6)
    return float(k_a), float(k_t)


def plan_subapertures(ann: S1Annotation, n_lines: int, cfg: Config,
                      burst_idx: int = 0, pixel: float = 0.0) -> MMPlan:
    """F06 -- banco di sub-aperture che resta DENTRO lo spettro.

    Regole dalle fonti: banda di guardia B_DL = B_cD/2 sempre sottratta; il
    master e' focalizzato su B_sub = B_cD - B_DL; master e slave sono tenuti a
    distanza rigida B_shift e marciano insieme. Perche' anche lo slave resti
    dentro lo spettro, la corsa disponibile e' B_cD - B_sub - B_shift, e il passo
    e' quella corsa divisa per N_D - 1. La versione precedente usava
    passo = B_sub/N_D indipendentemente da B_shift, e le ultime sub-aperture
    slave leggevano spettro vuoto.

    La conversione frequenza -> tempo usa il Doppler rate VERO letto
    dall'annotation, non una stima grossolana dalla durata del burst: sono due
    quantita' diverse di oltre un ordine di grandezza."""
    b_cd = ann.azimuth_bandwidth
    b_dl = b_cd * cfg.guard_fraction
    b_sub = b_cd - b_dl
    b_shift = b_cd * cfg.b_shift_fraction
    travel = b_cd - b_sub - b_shift
    if travel <= 0:
        b_shift = 0.5 * (b_cd - b_sub)
        travel = b_cd - b_sub - b_shift
    step = travel / max(cfg.n_d - 1, 1)

    df_bin = ann.prf / n_lines
    k_a, k_t = azimuth_rates(ann, burst_idx, pixel)
    ka = abs(k_a)
    t_illum = b_cd / max(ka, 1e-9)
    dt_step = step / max(ka, 1e-9)
    dt_shift = b_shift / max(ka, 1e-9)
    t_window = travel / max(ka, 1e-9)

    f_max_obs = 1.0 / (2.0 * max(dt_shift, 1e-12))
    # apertura sintetica effettivamente disponibile e risoluzione che ne
    # deriverebbe se si volesse fare tomografia ACUSTICA come nelle fonti
    v_ground = 6700.0
    v_sound = 6000.0
    aperture = v_ground * t_illum
    lam_sound = v_sound / max(f_max_obs, 1e-9)
    r_ref = 850_000.0
    dz_ac = lam_sound * r_ref / (2.0 * max(aperture, 1e-9))

    return MMPlan(
        b_cd=b_cd, b_dl=b_dl, b_sub=b_sub, b_shift=b_shift, step=step,
        n_d=cfg.n_d, df_bin=df_bin, bins_per_sub=b_sub / max(df_bin, 1e-9),
        k_a=k_a, k_t=k_t, t_illum=t_illum,
        dt_step=dt_step, dt_shift=dt_shift, t_window=t_window,
        f_sample=1.0 / max(dt_step, 1e-12),
        f_max_obs=f_max_obs,
        f_min_obs=1.0 / max(t_window, 1e-12),
        aperture_m=aperture, lambda_sound_m=lam_sound, delta_z_acoustic=dz_ac,
    )


def _tapered_band(n_lines: int, prf: float, f_lo: float, f_hi: float) -> np.ndarray:
    """F08 -- maschera di banda con finestra di Hann invece del rect.

    Il rect produce lobi laterali sinc in azimuth che sporcano il pixel tracking;
    la Hann li abbatte a costo di un modesto allargamento del lobo principale."""
    freqs = np.fft.fftfreq(n_lines, d=1.0 / prf)
    inside = (freqs >= f_lo) & (freqs < f_hi)
    m = np.zeros(n_lines, dtype=np.float32)
    if not inside.any():
        m[np.argmin(np.abs(freqs - 0.5 * (f_lo + f_hi)))] = 1.0
        return m
    u = (freqs[inside] - f_lo) / max(f_hi - f_lo, 1e-9)
    m[inside] = (0.5 - 0.5 * np.cos(2.0 * np.pi * u)).astype(np.float32)
    return m


def micro_motion_energy(
    img: np.ndarray, ann: S1Annotation, plan: MMPlan, cfg: Config,
    out_shape: Tuple[int, int], row_offset: int, verbose: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Energia di micro-moto per pixel dal banco di sub-aperture Doppler.

    Restituisce (energia, coerenza_della_traccia, frequenza_dominante)
    riportate sul reticolo del chip tomografico. `img` e' il chip esteso in
    azimuth (F07); `row_offset` e' la riga del chip esteso che corrisponde alla
    prima riga del chip tomografico.
    """
    from scipy.ndimage import map_coordinates

    spec = _block2_fft2(img)
    n_l, n_p = img.shape
    half = cfg.corr_window // 2
    stride = max(cfg.corr_window // 2, 1)
    ls = np.arange(half, n_l - half, stride)
    ps = np.arange(half, n_p - half, stride)
    if len(ls) == 0 or len(ps) == 0:
        z = np.zeros(out_shape, dtype=np.float32)
        return z, z.copy(), z.copy()
    gl, gp = np.meshgrid(ls, ps, indexing="ij")
    fl, fp = gl.ravel(), gp.ravel()

    shifts = np.zeros((len(fl), plan.n_d), dtype=np.complex64)
    f0 = -plan.b_cd / 2.0
    for k in range(plan.n_d):
        lo_m = f0 + k * plan.step
        mask_m = _tapered_band(n_l, ann.prf, lo_m, lo_m + plan.b_sub)
        mask_s = _tapered_band(n_l, ann.prf, lo_m + plan.b_shift,
                               lo_m + plan.b_sub + plan.b_shift)
        master = _block56_ifft2(_block34_bandpass(spec, mask_m))
        slave = _block56_ifft2(_block34_bandpass(spec, mask_s))

        mw = np.empty((len(fl), cfg.corr_window, cfg.corr_window), dtype=np.complex64)
        sw = np.empty_like(mw)
        for i, (l, p) in enumerate(zip(fl, fp)):
            mw[i] = master[l - half:l + half, p - half:p + half]
            sw[i] = slave[l - half:l + half, p - half:p + half]
        shifts[:, k] = _batch_subpixel_shift(mw, sw)
        if verbose:
            print(f"      sub-apertura {k + 1}/{plan.n_d}"
                  f"  [{lo_m:+7.1f} .. {lo_m + plan.b_sub:+7.1f}] Hz", flush=True)

    # F09: energia = deviazione standard del vettore COMPLESSO demediato.
    # std(|shift|) annulla una vibrazione a modulo costante, che e' proprio il
    # caso dell'oscillatore a 2 gradi di liberta' r(t) = (a cos, b sin) e^{-lt/2}.
    dev = shifts - shifts.mean(axis=1, keepdims=True)
    energy = np.sqrt((np.abs(dev) ** 2).mean(axis=1)).astype(np.float32)

    # F31: la "stabilita' della traccia" era |media(dev)| / media(|dev|), cioe'
    # il modulo della media di un vettore APPENA demediato: zero per
    # costruzione, su ogni pixel, sempre. Il canale esisteva ma non conteneva
    # nulla. La domanda giusta e' quella del modello a 2 gradi di liberta'
    # (ch12): la traccia e' un tono singolo o rumore bianco? Si risponde con lo
    # spettro della traccia, che e' anche il blocco 9 delle fonti applicato in
    # piccolo. Concentrazione = quota di energia nella riga dominante (esclusa
    # la continua, gia' tolta): 1 per un tono puro, ~1/(N_D-1) per rumore.
    F = np.fft.fft(dev, axis=1)
    P = (np.abs(F) ** 2)
    P[:, 0] = 0.0                       # la continua e' stata rimossa a monte
    tot = P.sum(axis=1)
    kpk = P.argmax(axis=1)
    conc = (P.max(axis=1) / np.maximum(tot, 1e-30)).astype(np.float32)
    conc[tot <= 0] = 0.0
    # frequenza meccanica della riga dominante, dal passo temporale della
    # marcia: N_D e' la frequenza di campionamento della vibrazione (ch11).
    f_axis = np.fft.fftfreq(plan.n_d, d=max(plan.dt_step, 1e-12))
    freq = np.abs(f_axis[kpk]).astype(np.float32)

    energy = energy.reshape(gl.shape)
    conc = conc.reshape(gl.shape)
    freq = freq.reshape(gl.shape)

    # F10: riporto sul reticolo del chip tomografico con map_coordinates sulle
    # coordinate REALI dei nodi, non con zoom (che spostava di mezza finestra).
    tl = np.arange(out_shape[0], dtype=np.float64) + row_offset
    tp = np.arange(out_shape[1], dtype=np.float64)
    ci = np.interp(tl, ls.astype(np.float64), np.arange(len(ls), dtype=np.float64))
    cj = np.interp(tp, ps.astype(np.float64), np.arange(len(ps), dtype=np.float64))
    ii, jj = np.meshgrid(ci, cj, indexing="ij")
    coords = np.stack([ii, jj])
    e = map_coordinates(energy, coords, order=1, mode="nearest").astype(np.float32)
    c = map_coordinates(conc, coords, order=1, mode="nearest").astype(np.float32)
    # la frequenza e' una etichetta discreta: interpolarla creerebbe righe che
    # il banco non puo' produrre, quindi si campiona al nodo piu' vicino
    f = map_coordinates(freq, coords, order=0, mode="nearest").astype(np.float32)
    return e, c, f


# ==========================================================================
# 9.  Attributi derivati e geometria locale
# ==========================================================================

def _norm01(a: np.ndarray, lo_p: float = 2.0, hi_p: float = 98.0) -> np.ndarray:
    a = np.nan_to_num(a.astype(np.float32))
    lo, hi = np.percentile(a, lo_p), np.percentile(a, hi_p)
    return np.clip((a - lo) / max(float(hi - lo), 1e-9), 0.0, 1.0)


def solidity_index(tomo_mag, coherence, mm_energy) -> np.ndarray:
    """Discriminante multi-attributo pieno/vuoto, dichiarato come euristica.

        indice = norm(energia) * norm(coerenza) * (1 - norm(micro_moto))

    Alto = ritorno forte, diffusore stabile, poca vibrazione anomala. NON e' una
    rilevazione di cavita' risolta in profondita': la risoluzione verticale non
    lo consente, e le fonti stesse chiedono di tenere separata la misura
    dall'interpretazione."""
    e = _norm01(tomo_mag)
    c = _norm01(coherence)[..., None]
    m = _norm01(mm_energy)[..., None]
    return (e * c * (1.0 - m)).astype(np.float32)


def vertical_lobe_profile(
    tomo_mag: np.ndarray, z_axis: np.ndarray, pyr_mask: np.ndarray,
    delta_z_vertical: float,
) -> Dict[str, Any]:
    """Profilo medio dell'energia lungo z: la PSF verticale della pila (F41).

    Risponde alla domanda "perche' i voxel stanno anche sopra il suolo". Con
    poche baseline la cella di Rayleigh verticale e' piu' larga dell'oggetto:
    la risposta di UN diffusore di superficie si spalma lungo tutto l'asse di
    ricerca e i lobi laterali finiscono in aria. Questo profilo e' quella
    risposta, MISURATA sui dati e non argomentata. Se il fascio non scende mai
    a -3 dB dentro l'asse non esiste un lobo principale separabile, e ogni
    struttura verticale nella nuvola e' PSF dell'array, non stratigrafia.
    """
    z = np.asarray(z_axis, dtype=np.float64)
    msk = np.asarray(pyr_mask, dtype=bool)

    def _curva(sel: np.ndarray) -> np.ndarray:
        if not sel.any():
            return np.full(len(z), -99.0)
        p = tomo_mag[sel].mean(axis=0).astype(np.float64)
        return 20.0 * np.log10(np.maximum(p, 1e-12) / max(float(p.max()), 1e-12))

    d_pyr = _curva(msk)
    d_all = _curva(np.ones(msk.shape, dtype=bool))

    # lobo principale: il tratto CONTIGUO sopra -3 dB che contiene il picco.
    k_pk = int(np.argmax(d_pyr))
    lo = hi = k_pk
    while lo > 0 and d_pyr[lo - 1] >= -3.0:
        lo -= 1
    while hi < len(z) - 1 and d_pyr[hi + 1] >= -3.0:
        hi += 1
    troncato = bool(lo == 0 or hi == len(z) - 1)

    # lobo laterale peggiore FUORI dal lobo principale; se il lobo principale
    # arriva ai bordi dell'asse non esiste un "fuori" e il campo resta nullo.
    fuori = np.ones(len(z), dtype=bool)
    fuori[lo:hi + 1] = False
    lobo_lat = float(d_pyr[fuori].max()) if fuori.any() else None

    passo = float(z[1] - z[0]) if len(z) > 1 else 0.0
    return {
        "z": np.round(z, 1).tolist(),
        "pyr": np.round(d_pyr, 2).tolist(),
        "tutto": np.round(d_all, 2).tolist(),
        "z_picco": round(float(z[k_pk]), 1),
        "lobo_lo": round(float(z[lo]), 1),
        "lobo_hi": round(float(z[hi]), 1),
        "lobo_larghezza": round(float(z[hi] - z[lo]), 1),
        "lobo_troncato": troncato,
        "lobo_laterale_db": None if lobo_lat is None else round(lobo_lat, 2),
        "contrasto_db": round(float(d_pyr.max() - d_pyr.min()), 2),
        "bordi_db": [round(float(d_pyr[0]), 2), round(float(d_pyr[-1]), 2)],
        "passo_z": round(passo, 3),
        "delta_z": round(float(delta_z_vertical), 1),
        "sovracamp": (round(float(delta_z_vertical / passo), 0)
                      if passo > 0 else None),
    }


def local_enu(geo: Geocoder, win: ChipWindow, sl: int, sp: int,
              shape: Tuple[int, int]
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Tuple[float, float]]:
    """Coordinate locali est/nord [m] del reticolo multilooked.

    F14: le coordinate sono quelle del CENTRO di ciascun blocco di multilooking,
    non del primo campione, e la geocodifica passa dalla bilineare locale
    sulla geolocation grid dello .xml (F36)."""
    n_l, n_p = shape
    ll = win.l0 + (np.arange(n_l) * sl + (sl - 1) / 2.0)
    pp = win.p0 + (np.arange(n_p) * sp + (sp - 1) / 2.0)
    gl, gp = np.meshgrid(ll, pp, indexing="ij")
    lat, lon, h_ref = geo.llh(gl.ravel(), gp.ravel())
    lat = lat.reshape(gl.shape)
    lon = lon.reshape(gl.shape)
    h_ref = h_ref.reshape(gl.shape).astype(np.float32)

    lat0, lon0 = PYRAMIDS[0].lat, PYRAMIDS[0].lon
    e_m, n_m = _enu_offset(lat, lon, lat0, lon0)
    east = e_m.astype(np.float32)
    north = n_m.astype(np.float32)
    return east, north, h_ref, (lat0, lon0)


def raw_gcp_nodes(geo: Geocoder, lat0: float, lon0: float) -> Dict[str, Any]:
    """I nodi VERI della geolocation grid dello .xml -- non la spline bicubica.

    F35: "height_ref"/"xml_ref" (local_enu, sopra) sono la spline valutata
    fittamente sul reticolo multilooked delle piramidi: SEMBRANO un campo
    continuo pixel-per-pixel perche' la valutazione e' fitta, ma l'informazione
    che portano non lo e'. La geolocation grid di un annotation IW ha 11 x 21 =
    231 nodi su TUTTA la scena (~250 x 250 km), quindi una spaziatura di
    centinaia di linee x centinaia di pixel -- vedi i valori reali sotto,
    stimati da geo.lines/geo.pixels. La finestra delle piramidi e' larga una
    piccola frazione di UNA sola cella di quel reticolo: dentro quella cella la
    spline non ha nessun dato aggiuntivo da interpolare, quindi il risultato è
    correttamente quasi piatto (F34, ~2 m di escursione). Non è un errore di
    calcolo: è il limite di risoluzione del riferimento .xml.

    Questa funzione restituisce gli stessi nodi SENZA passare dalla spline: la
    posizione fisica di ciascun pixel-nodo del reticolo, con la sua quota
    propria. Sull'intera scena il rilievo vero è tutt'altro che un piano (fino
    a diverse centinaia di metri di escursione, vedi il nodo piu' vicino
    stampato in console): la spline dentro il ritaglio non lo mostra solo
    perche' nessun nodo reale cade li' vicino."""
    lat = geo.lat_nodes.astype(np.float64)
    lon = geo.lon_nodes.astype(np.float64)
    h = geo.height_nodes.astype(np.float64)
    east, north = _enu_offset(lat, lon, lat0, lon0)
    ll, pp = np.meshgrid(geo.lines, geo.pixels, indexing="ij")
    return {
        "n_l": int(lat.shape[0]), "n_p": int(lat.shape[1]),
        "line": ll.ravel().astype(int).tolist(),
        "pixel": pp.ravel().astype(int).tolist(),
        "lat": np.round(lat, 6).ravel().tolist(),
        "lon": np.round(lon, 6).ravel().tolist(),
        "h": np.round(h, 1).ravel().tolist(),
        "east": np.round(east, 1).ravel().tolist(),
        "north": np.round(north, 1).ravel().tolist(),
    }


def plateau_heights_from_stack(cfg: Config, geo: Geocoder,
                               verbose: bool = True) -> Dict[str, Any]:
    """Quota del terreno alle tre piramidi LETTA DAI FILE DELLO STACK.

    F39. Le uniche quote presenti nei file di ``stack_slc`` sono i
    ``<geolocationGridPoint><height>`` degli ``*.annotation.xml`` (F37: i
    ``<position>`` sono vettori di stato ORBITALI, non punti al suolo).
    Questa funzione le legge da OGNI data dello stack, non solo dal master:

    * il valore usato e' la bilineare LOCALE del master (F36), l'unica
      geometricamente corretta per l'intero stack, perche' tutte le date sono
      coregistrate sul reticolo di pixel del master;
    * le altre date servono da controllo: se le quote lette cambiassero da una
      data all'altra il datum non sarebbe unico e il layer del suolo non
      sarebbe confrontabile con la superficie misurata.

    Verificato sui dati (11 date VH, IW2 di ``goal_out_kefren/stack_slc``):
    tutte e 11 danno lo stesso nodo grezzo a 64.0 m (il piu' vicino cade a
    2.2-3.2 km dal ritaglio) e la bilineare locale del master da' 63.1..64.2 m.
    Quello e' il livello del PLATEAU: la grid non ha nodi abbastanza vicini per
    contenere le piramidi (F35), che vanno quindi aggiunte a parte da
    ``pyramid_profile_enu()``."""
    per = []
    for p in PYRAMIDS:
        l, q = geo.latlon_to_line_pixel(p.lat, p.lon)
        per.append({"nome": p.name,
                    "h_xml_master": float(geo.llh(l, q)[2]),
                    "base_alt_letteratura": float(p.base_alt_m)})

    # controllo di consistenza su tutte le date dello stack
    per_data: List[Dict[str, Any]] = []
    # unico posto che NON usa _enu_offset, e di proposito: qui il fattore di
    # longitudine e' quello del centro della scena (PYRAMIDS[0]) mentre lo
    # scostamento e' misurato da CIASCUNA piramide, quindi i due riferimenti
    # non coincidono e la sostituzione cambierebbe il risultato.
    m_lon = M_PER_GRADO_LON_EQ * math.cos(math.radians(PYRAMIDS[0].lat))
    # verbose=False: l'avviso sulle tracce (F43) l'ha gia' dato run(); qui si
    # rileggono soltanto le quote dei nodi.
    for entry in discover_stack(cfg, verbose=False):
        try:
            root = ET.parse(entry.annotation).getroot()
            gp = root.findall(".//geolocationGrid//geolocationGridPoint")
            gl = np.array([float(_txt(g, "latitude")) for g in gp])
            go = np.array([float(_txt(g, "longitude")) for g in gp])
            gh = np.array([float(_txt(g, "height")) for g in gp])
        except Exception as ex:                                # pragma: no cover
            if verbose:
                print(f"      [suolo] {os.path.basename(entry.annotation)}: "
                      f"lettura fallita ({ex})")
            continue
        row: Dict[str, Any] = {"data": entry.date, "h": [], "d_km": []}
        for p in PYRAMIDS:
            dist = np.hypot((gl - p.lat) * M_PER_GRADO_LAT, (go - p.lon) * m_lon)
            k = int(np.argmin(dist))
            row["h"].append(float(gh[k]))
            row["d_km"].append(float(dist[k] / 1000.0))
        per_data.append(row)

    hh = np.array([r["h"] for r in per_data]) if per_data else None
    if verbose:
        print(f"      quote lette dai {len(per_data)} annotation.xml di "
              f"{cfg.stack_dir}")
        for i, d in enumerate(per):
            if hh is not None:
                nodo = (f"nodo grezzo {hh[:, i].min():.1f}..{hh[:, i].max():.1f} m "
                        f"a {np.mean([r['d_km'][i] for r in per_data]):.1f} km")
            else:
                nodo = "nodo grezzo non disponibile"
            print(f"        {d['nome']:<22} bilineare .xml {d['h_xml_master']:6.2f} m"
                  f"  ({nodo}; letteratura {d['base_alt_letteratura']:.0f} m)")
        if hh is not None and float(hh.max() - hh.min()) > 0.5:
            print("        ATTENZIONE: le date non concordano sulla quota di "
                  "riferimento: il layer del suolo non sta su un datum unico")
    return {"per_piramide": per, "per_data": per_data, "n_date": len(per_data),
            "concordi": bool(hh is not None and float(hh.max() - hh.min()) <= 0.5)}


def _pyramid_local_xy(east: np.ndarray, north: np.ndarray, p: Pyramid,
                      lat0: float, lon0: float) -> Tuple[np.ndarray, np.ndarray]:
    """(east, north) riportati nel sistema della base di ``p``.

    Rotazione INVERSA di ``azimuth_deg``, cioe' l'inversa esatta della
    convenzione con cui ``pyramid_mesh()`` costruisce i quattro spigoli."""
    cx, cy = _enu_offset(p.lat, p.lon, lat0, lon0)
    a = math.radians(p.azimuth_deg)
    ca, sa = math.cos(a), math.sin(a)
    de = np.asarray(east, dtype=np.float64) - cx
    dn = np.asarray(north, dtype=np.float64) - cy
    return de * ca + dn * sa, -de * sa + dn * ca


def pyramid_footprints_enu(east: np.ndarray, north: np.ndarray, lat0: float,
                           lon0: float, scale: float = 1.0) -> List[np.ndarray]:
    """Impronte al suolo delle tre piramidi sul reticolo ENU (una maschera
    per piramide)."""
    out = []
    for p in PYRAMIDS:
        dx, dy = _pyramid_local_xy(east, north, p, lat0, lon0)
        half = scale * p.base_side_m / 2.0
        out.append((np.abs(dx) <= half) & (np.abs(dy) <= half))
    return out


def pyramid_profile_enu(east: np.ndarray, north: np.ndarray, lat0: float,
                        lon0: float, base_alt: Sequence[float]) -> np.ndarray:
    """Quota della SUPERFICIE delle tre piramidi sul reticolo (east, north).

    F39. Piramide retta a base quadrata: dentro l'impronta, con (dx, dy)
    riportati nel sistema della base,

        z = base + altezza * (1 - max(|dx|, |dy|) / semilato)

    che vale ``base`` sul perimetro e ``base + altezza`` all'apice. E' la
    stessa formula gia' usata da ``simulate_pyramids_radar()``, li' per la
    proiezione in geometria radar, qui AL SUOLO, perche' questo e' un layer di
    terreno e non una simulazione di ritorno radar.

    ``base_alt[k]`` e' la quota del terreno sotto la piramide k, presa dal
    terreno stesso (ancorato alle quote lette dai file dello stack, vedi
    ``plateau_heights_from_stack()``), non dal valore di letteratura
    ``Pyramid.base_alt_m``: il layer deve essere continuo con il terreno che
    lo circonda, altrimenti le piramidi galleggiano o affondano.

    Fuori da ogni impronta restituisce ``-inf``, cosi' il chiamante compone
    terreno e piramidi con un ``np.maximum``. Il PERIMETRO appartiene alla
    piramide (``frac >= 0``, con la tolleranza che serve perche' i quattro
    spigoli calcolati da ``pyramid_mesh()`` cadano dentro nonostante
    l'arrotondamento): la stessa convenzione di ``pyramid_footprints_enu()``,
    altrimenti maschera e profilo non coinciderebbero sul bordo."""
    z = np.full(np.shape(east), -np.inf, dtype=np.float64)
    for k, p in enumerate(PYRAMIDS):
        dx, dy = _pyramid_local_xy(east, north, p, lat0, lon0)
        half = p.base_side_m / 2.0
        frac = 1.0 - np.maximum(np.abs(dx), np.abs(dy)) / half
        inside = frac >= -1e-9
        zk = float(base_alt[k]) + p.height_m * np.maximum(frac, 0.0)
        z = np.where(inside, np.maximum(z, zk), z)
    return z


def _open_meteo_elevation(gclat: np.ndarray, gclon: np.ndarray,
                          chunk: int = 100) -> Optional[np.ndarray]:
    """Elevazioni Copernicus dall'API Open-Meteo su un reticolo lat x lon.

    F39: l'API accetta al massimo 100 coordinate per richiesta. La versione
    precedente ne mandava una sola, quindi il reticolo non poteva superare i
    10x10 nodi: su un ritaglio di 1.6 x 1.3 km sono ~170 m di passo, piu'
    grossolani del terreno che si voleva mostrare. Qui la richiesta e'
    spezzata in blocchi, cosi' il passo lo decide ``cfg.dem_grid``."""
    CLAT, CLON = np.meshgrid(gclat, gclon, indexing="ij")
    fl, fo = CLAT.ravel(), CLON.ravel()
    out = np.empty(fl.size, dtype=np.float64)
    for s in range(0, fl.size, chunk):
        qs = ",".join(f"{v:.6f}" for v in fl[s:s + chunk])
        qo = ",".join(f"{v:.6f}" for v in fo[s:s + chunk])
        url = f"https://api.open-meteo.com/v1/elevation?latitude={qs}&longitude={qo}"
        ok = False
        for att in range(5):
            try:
                with urllib.request.urlopen(url, timeout=40) as r:
                    out[s:s + chunk] = np.array(json.load(r)["elevation"],
                                                dtype=np.float64)
                ok = True
                break
            except Exception as ex:
                print(f"      [DEM] blocco {s // chunk} tentativo {att} fallito ({ex})")
                time.sleep(2 + 2 * att)
        if not ok:
            return None
    return out.reshape(CLAT.shape)


def ground_dem_suolo(lat0: float, lon0: float, east: np.ndarray, north: np.ndarray,
                     h_xml: np.ndarray, out_dir: str, gc: int = 24,
                     cache_name: str = "dem_suolo.npz", use_external: bool = True,
                     add_pyramids: bool = True, verbose: bool = True
                     ) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Layer "suolo (DEM)": il terreno del ritaglio, PIRAMIDI COMPRESE.

    F39 -- che cosa era sbagliato. Il layer (F38) era il solo DEM esterno
    (Copernicus via l'API di elevazione Open-Meteo). Verificato interrogando
    l'API sui tre apici noti: restituisce 66 / 76 / 78 m, cioe' il PLATEAU
    NUDO. Le piramidi in quel dataset non ci sono: l'apice di Cheope, che sta
    a ~198 m, usciva 130 m piu' in basso, e il "suolo" invece di seguire il
    profilo delle piramidi passava dentro. Non era un problema di passo del
    reticolo: anche interrogando l'API esattamente sull'apice il valore resta
    66 m.

    Come e' costruito ora, in tre pezzi dichiarati:

    1. **datum**: le quote lette dai file di ``stack_slc``
       (``<geolocationGridPoint><height>`` degli ``*.annotation.xml``,
       bilineare locale del master, F36), 63.1..64.2 m sul ritaglio e concordi
       su tutte le date (``plateau_heights_from_stack()``).
    2. **rilievo locale**: il DEM esterno, se scaricabile, RIPORTATO su quel
       datum -- gli si toglie il suo scarto mediano dalle quote .xml
       sull'impronta delle piramidi. Resta il rilievo del plateau e della
       scarpata, ma il livello e' quello dei file dello stack e non quello di
       un dataset con un datum verticale diverso. Senza rete il terreno resta
       la bilineare .xml e il layer c'e' lo stesso (prima spariva).
    3. **profilo delle piramidi** (``pyramid_profile_enu()``), appoggiato sul
       terreno del punto 2: base = mediana del terreno sotto ciascuna
       impronta, apice = base + altezza nota. Composizione con ``np.maximum``,
       cosi' il suolo segue il terreno fuori dalle impronte e le facce delle
       piramidi dentro.

    Resta un RIFERIMENTO: non entra in nessun calcolo tomografico -- fase, k_z
    e geocodifica restano quelli della geometria .xml."""
    from scipy.interpolate import RegularGridInterpolator

    east = np.asarray(east, dtype=np.float64)
    north = np.asarray(north, dtype=np.float64)
    h_xml = np.asarray(h_xml, dtype=np.float64)
    info: Dict[str, Any] = {
        "datum_xml_m": [round(float(h_xml.min()), 2), round(float(h_xml.max()), 2)]}

    lat, lon = _latlon_from_enu(east, north, lat0, lon0)
    latS, latN = float(lat.min()), float(lat.max())
    lonW, lonE = float(lon.min()), float(lon.max())
    pad_lat = max(0.002, 0.15 * (latN - latS))
    pad_lon = max(0.002, 0.15 * (lonE - lonW))
    latS, latN = latS - pad_lat, latN + pad_lat
    lonW, lonE = lonW - pad_lon, lonE + pad_lon
    bbox = np.array([latS, latN, lonW, lonE])

    # --- 2. rilievo del terreno da DEM esterno (opzionale) -----------------
    terr = None
    if use_external:
        cache = os.path.join(out_dir, cache_name)
        if os.path.exists(cache):
            try:
                d = np.load(cache)
                if (tuple(int(v) for v in d["shape"]) == east.shape
                        and np.allclose(d["bbox"], bbox, atol=1e-4)
                        and int(d["gc"]) == int(gc)):
                    terr = d["z0"].astype(np.float64)
                    if verbose:
                        print(f"      [DEM] riuso {cache}: rilievo esterno "
                              f"{terr.min():.0f}..{terr.max():.0f} m")
            except Exception:
                terr = None
        if terr is None:
            gclat = np.linspace(latS, latN, gc)
            gclon = np.linspace(lonW, lonE, gc)
            ce = _open_meteo_elevation(gclat, gclon)
            if ce is None:
                if verbose:
                    print("      [DEM] scaricamento fallito: il terreno resta la "
                          "bilineare .xml, il layer del suolo c'e' lo stesso (F39)")
            else:
                itp = RegularGridInterpolator((gclat, gclon), ce,
                                              bounds_error=False, fill_value=None)
                terr = itp(np.stack([lat.ravel(), lon.ravel()],
                                    axis=-1)).reshape(east.shape)
                np.savez_compressed(cache, z0=terr, shape=np.array(east.shape),
                                    bbox=bbox, gc=np.array(int(gc)))
                if verbose:
                    # F45: era `m_lat`, che in questo scope non esiste. Il
                    # ramo scatta solo al primo scaricamento riuscito (senza
                    # cache) e con verbose: bastava avere gia' il .npz per non
                    # vederlo mai, ma li' e' un NameError, non un numero storto.
                    step_m = (latN - latS) * M_PER_GRADO_LAT / max(gc - 1, 1)
                    print(f"      [DEM] rilievo esterno (Copernicus via Open-Meteo): "
                          f"{terr.min():.0f}..{terr.max():.0f} m su reticolo "
                          f"{gc}x{gc}, passo ~{step_m:.0f} m")

    foot = pyramid_footprints_enu(east, north, lat0, lon0, 1.0)
    any_foot = np.zeros(east.shape, dtype=bool)
    for f in foot:
        any_foot |= f

    # --- 1. datum: le quote lette dai file dello stack ---------------------
    if terr is None:
        ground = h_xml.copy()
        info["rilievo_esterno"] = False
        info["bias_datum_m"] = 0.0
    else:
        ref = any_foot if any_foot.any() else np.ones(east.shape, dtype=bool)
        bias = float(np.median(terr[ref]) - np.median(h_xml[ref]))
        ground = terr - bias
        info["rilievo_esterno"] = True
        info["bias_datum_m"] = round(bias, 2)
        if verbose:
            print(f"      [suolo] rilievo esterno riportato sul datum .xml: "
                  f"{-bias:+.1f} m (mediana sull'impronta delle piramidi) -> "
                  f"terreno {ground.min():.0f}..{ground.max():.0f} m")
    info["terreno_m"] = [round(float(ground.min()), 1), round(float(ground.max()), 1)]

    # --- 3. profilo delle piramidi -----------------------------------------
    base_alt: List[float] = [
        float(np.median(ground[foot[k]])) if foot[k].any() else float(np.median(ground))
        for k in range(len(PYRAMIDS))]
    info["base_alt_m"] = [round(b, 2) for b in base_alt]
    info["apici_m"] = [round(b + p.height_m, 1) for b, p in zip(base_alt, PYRAMIDS)]
    info["piramidi"] = bool(add_pyramids)
    info["nodi_impronta"] = [int(f.sum()) for f in foot]
    # scarto della base ricavata dal terreno rispetto al valore di letteratura.
    # Non e' una correzione: e' il residuo dichiarato fra il datum dei file
    # dello stack (63.1..64.2 m sul ritaglio, quasi piatto per F35) e le quote
    # di base note delle singole piramidi, che quel datum non puo' risolvere.
    info["scarto_letteratura_m"] = [round(b - p.base_alt_m, 1)
                                    for b, p in zip(base_alt, PYRAMIDS)]

    if add_pyramids:
        zp = pyramid_profile_enu(east, north, lat0, lon0, base_alt)
        ground = np.maximum(ground, zp)
        if verbose:
            for k, p in enumerate(PYRAMIDS):
                print(f"      [suolo] {p.name:<22} base {base_alt[k]:6.1f} m   "
                      f"apice {base_alt[k] + p.height_m:6.1f} m   "
                      f"impronta {int(foot[k].sum())} nodi   "
                      f"(letteratura: base {p.base_alt_m:.0f} m, "
                      f"scarto {base_alt[k] - p.base_alt_m:+.1f} m)")
            print(f"      [suolo] con le piramidi: {ground.min():.0f}..{ground.max():.0f} m "
                  f"(terreno da solo {info['terreno_m'][0]:.0f}..{info['terreno_m'][1]:.0f} m)")

    info["suolo_m"] = [round(float(ground.min()), 1), round(float(ground.max()), 1)]
    return ground.astype(np.float64), info


def pyramid_mesh(p: Pyramid, lat0: float, lon0: float) -> Dict[str, Any]:
    """Vertici e facce della piramide ideale in ENU locale [m] -- riferimento."""
    cx, cy = _enu_offset(p.lat, p.lon, lat0, lon0)
    h = p.base_side_m / 2.0
    a = math.radians(p.azimuth_deg)
    ca, sa = math.cos(a), math.sin(a)
    corners = [[cx + dx * ca - dy * sa, cy + dx * sa + dy * ca, p.base_alt_m]
               for dx, dy in ((-h, -h), (h, -h), (h, h), (-h, h))]
    apex = [cx, cy, p.base_alt_m + p.height_m]
    return {
        "name": p.name, "vertices": corners + [apex], "apex": apex,
        "base_side": p.base_side_m, "height": p.height_m,
        "base_alt": p.base_alt_m, "slope": p.face_slope_deg,
    }


# --------------------------------------------------------------------------
# Strutture interne note  (F32)  --  RIFERIMENTO archeologico, non misura
# --------------------------------------------------------------------------
# La geometria viene da piramide_cheope_3d.py e piramide_kefren_3d.py, che
# sono la sorgente unica di quelle quote (rilievo Petrie 1883 per Cheope,
# Legon / Vyse & Perring / Lehner per Chefren). Qui non entra in nessun
# calcolo: sta nella scena 3D per far vedere DOVE stanno le camere note
# rispetto alla superficie ricostruita, ed e' disegnata come i riferimenti,
# a filo di ferro, mai come punti misurati. Con delta_z verticale di 242 m e
# camere che misurano metri, nessuna di queste strutture e' rilevabile da
# questi dati: il confronto serve a dare la scala dell'impossibilita', non a
# suggerire una rilevazione.

STRUCTURE_SOURCES: Tuple[Tuple[str, str], ...] = (
    ("piramide_cheope_3d", "Cheope (Khufu)"),
    ("piramide_kefren_3d", "Chefren (Khafre)"),
)


def _structure_edges(s: Dict[str, Any]) -> List[Tuple[Tuple[float, float, float],
                                                         Tuple[float, float, float]]]:
    """Spigoli di una struttura nel sistema locale della piramide.

    X = est, Y = nord, Z = quota sopra il piano di base. I tre tipi usati dai
    due file sorgente: ``box`` (parallelepipedo), ``cuspide`` (stesso ingombro
    ma con tetto a doppio spiovente, colmo lungo il lato maggiore) e
    ``corridoio`` (prisma a sezione rettangolare con asse lungo Y, inclinato
    fra i due estremi)."""
    def rect(pts):
        return [(pts[i], pts[(i + 1) % 4]) for i in range(4)]

    t = s["tipo"]
    if t in ("box", "cuspide"):
        x, y, z = s["x"], s["y"], s["z"]
        hx, hy = s["dx"] / 2.0, s["dy"] / 2.0
        z1 = z + s["dz"]
        low = [(x - hx, y - hy, z), (x + hx, y - hy, z),
               (x + hx, y + hy, z), (x - hx, y + hy, z)]
        top = [(a, b, z1) for a, b, _ in low]
        e = rect(low) + rect(top) + [(low[i], top[i]) for i in range(4)]
        if t == "cuspide":
            zc = z + s.get("dz_colmo", s["dz"])
            # il colmo corre lungo il lato maggiore
            if s["dx"] >= s["dy"]:
                ridge = [(x - hx, y, zc), (x + hx, y, zc)]
                rafters = [(ridge[0], top[0]), (ridge[0], top[3]),
                           (ridge[1], top[1]), (ridge[1], top[2])]
            else:
                ridge = [(x, y - hy, zc), (x, y + hy, zc)]
                rafters = [(ridge[0], top[0]), (ridge[0], top[1]),
                           (ridge[1], top[2]), (ridge[1], top[3])]
            e += [(ridge[0], ridge[1])] + rafters
        return e
    if t == "corridoio":
        x, w, h = s["x"], s["w"], s["h"]
        hx = w / 2.0
        ends = []
        for y, z in ((s["y0"], s["z0"]), (s["y1"], s["z1"])):
            ends.append([(x - hx, y, z), (x + hx, y, z),
                         (x + hx, y, z + h), (x - hx, y, z + h)])
        a, b = ends
        return rect(a) + rect(b) + [(a[i], b[i]) for i in range(4)]
    raise ValueError(f"tipo di struttura sconosciuto: {t}")


def known_structures(lat0: float, lon0: float,
                     verbose: bool = True) -> List[Dict[str, Any]]:
    """Strutture interne note delle piramidi, in ENU locale [m].

    Se uno dei due moduli manca il programma continua senza quel gruppo: la
    scena 3D perde un riferimento, non un risultato."""
    from matplotlib.colors import to_hex

    by_name = {p.name: p for p in PYRAMIDS}
    out: List[Dict[str, Any]] = []

    for mod_name, pyr_name in STRUCTURE_SOURCES:
        p = by_name.get(pyr_name)
        if p is None:
            continue
        try:
            mod = __import__(mod_name)
            strutture = list(mod.STRUTTURE)
        except Exception as exc:                       # pragma: no cover
            if verbose:
                print(f"      strutture: {mod_name} non disponibile ({exc})")
            continue

        cx, cy = _enu_offset(p.lat, p.lon, lat0, lon0)
        ang = math.radians(p.azimuth_deg)
        ca, sa = math.cos(ang), math.sin(ang)

        for s in strutture:
            try:
                edges = _structure_edges(s)
            except ValueError as exc:                  # pragma: no cover
                if verbose:
                    print(f"      strutture: {exc}")
                continue
            flat = []
            for a, b in edges:
                for vx, vy, vz in (a, b):
                    flat += [round(cx + vx * ca - vy * sa, 2),
                             round(cy + vx * sa + vy * ca, 2),
                             round(p.base_alt_m + vz, 2)]
            out.append({
                "pyr": pyr_name,
                "num": int(s["num"]),
                "nome": str(s["nome"]),
                "tipo": str(s["tipo"]),
                "colore": to_hex(s["colore"]),
                "certezza": str(s.get("certezza", "")),
                "edges": flat,
            })
        if verbose:
            print(f"      strutture note: {len(strutture)} da {mod_name} "
                  f"-> {pyr_name}")
    return out


# ==========================================================================
# 10.  Pipeline
# ==========================================================================

def polarisation_contrast(cfg: Config, verbose: bool = True) -> Dict[str, Any]:
    """F26 -- quanto emergono le piramidi in ciascuna polarizzazione disponibile.

    E' un controllo di coerenza fra metodi nel senso delle fonti: si misura, si
    riporta anche il disaccordo, e la scelta del canale smette di essere una
    giustificazione a parole e diventa un numero. Il contrasto e' calcolato
    nella zona di LAYOVER simulata, cioe' dove la piramide finisce davvero, non
    sull'impronta al suolo."""
    from scipy.ndimage import uniform_filter

    out: Dict[str, Any] = {}
    for pol in ("vh", "vv"):
        c = Config(**{**asdict(cfg), "polarisation": pol})
        try:
            # verbose=False: l'avviso sulle tracce l'ha gia' dato run()
            entries = discover_stack(c, verbose=False)
        except FileNotFoundError:
            continue
        e = entries[min(len(entries) // 2, len(entries) - 1)]
        ann = parse_annotation(e.annotation)
        geo = Geocoder(ann)
        win, _ = target_window(ann, geo, c)
        ll = np.arange(win.l0, win.l1 + 1.0)
        pp = np.arange(win.p0, win.p1 + 1.0)
        gl, gp = np.meshgrid(ll, pp, indexing="ij")
        sim = simulate_pyramids_radar(geo, ann, win, geo.incidence(gl, gp))

        # F40: qui si misurano dei dB, quindi conviene che siano dB di sigma0
        # e non di conteggi DN. Con i .xml del prodotto il contrasto diventa
        # una grandezza fisica confrontabile fra date, swath e polarizzazioni,
        # e togliendo il NESZ non si sta piu' misurando anche il rumore -- che
        # in VH sul deserto e' una frazione non trascurabile del fondo.
        chip = read_window(e, win, 0, calibra=c.calibrazione and bool(e.calibration))
        pw = uniform_filter(np.abs(chip.data) ** 2, (3, 9))
        if c.calibrazione:
            n_map = mappa_rumore(e, win)
            if n_map is not None:
                pw = np.maximum(pw - uniform_filter(n_map, (3, 9)), 0.0)

        m, sh = sim["sim_mask"], sim["sim_h"]
        bg = 10.0 * math.log10(float(pw[~m].mean()))
        bands = {}
        for lo, hi, lab in ((80.0, 1e9, "layover_alto_gt80m"),
                            (40.0, 80.0, "medio_40_80m"),
                            (0.0, 40.0, "basso_lt40m")):
            s = m & (sh >= lo) & (sh < hi)
            if s.sum() < 10:
                continue
            bands[lab] = round(10.0 * math.log10(float(pw[s].mean())) - bg, 2)
        peaks = {}
        for d in sim["per"]:
            s = d["mask"] & (sh >= 60.0)
            if s.sum() >= 5:
                peaks[d["nome"]] = round(
                    10.0 * math.log10(float(pw[s].max())) - bg, 2)
        out[pol] = {"data": e.date, "fondo_dB": round(bg, 2),
                    "contrasto_per_fascia_dB": bands, "picco_dB": peaks}
        if verbose:
            print(f"      {pol.upper()} ({e.date}): contrasto nella zona di layover "
                  + ", ".join(f"{k.split('_')[0]} {v:+.1f} dB"
                              for k, v in bands.items())
                  + (f"  |  picco Cheope {peaks.get(PYRAMIDS[0].name, float('nan')):+.1f} dB"
                     if PYRAMIDS[0].name in peaks else ""))
    if "vh" in out and "vv" in out:
        pk = PYRAMIDS[0].name
        a = out["vh"]["picco_dB"].get(pk)
        b = out["vv"]["picco_dB"].get(pk)
        if a is not None and b is not None:
            out["conclusione"] = (
                f"Il ritorno delle piramidi e' co-polarizzato: picco su Cheope "
                f"{b:+.1f} dB in VV contro {a:+.1f} dB in VH. Usare VH costa "
                f"{b - a:.1f} dB proprio sui bersagli da misurare. La scelta del "
                f"canale e' rispettata come richiesta, ma il costo e' misurato e "
                f"dichiarato, non giustificato a parole.")
    return out


def _data_da_nome(nome: str) -> str:
    m = re.search(r"-(\d{8})t", os.path.basename(nome))
    return m.group(1) if m else "?"


def _stack_piatto(cfg: Config) -> List[StackEntry]:
    """Layout storico: <platform>-<swath>-slc-<pol>-*.annotation.xml + .tiff."""
    pattern = os.path.join(
        cfg.stack_dir,
        f"{cfg.platform}-{cfg.swath}-slc-{cfg.polarisation}-*.annotation.xml",
    )
    out: List[StackEntry] = []
    for ann in sorted(glob.glob(pattern)):
        tif = ann.replace(".annotation.xml", ".tiff")
        if os.path.exists(tif):
            out.append(StackEntry(_data_da_nome(ann), ann, tif))
    return out


def _stack_safe(cfg: Config) -> List[StackEntry]:
    """Alberi .SAFE come li consegna CDSE, letti dove sono.

    Evita di dover costruire una cartella-ponte con copie o hard link: il
    prodotto Sentinel-1 ha gia' tutto al posto giusto, e da qui si raccolgono
    anche i due .xml di calibrazione e rumore che il layout piatto non porta.
    La piattaforma NON filtra: un S1A e un S1C sulla stessa orbita relativa
    appartengono alla stessa pila multi-baseline, ed e' proprio la diversita'
    di baseline che serve."""
    marca = f"-{cfg.swath}-slc-{cfg.polarisation}-".lower()
    out: List[StackEntry] = []
    for safe in sorted(glob.glob(os.path.join(cfg.stack_dir, "*.SAFE"))):
        ann = [a for a in glob.glob(os.path.join(safe, "annotation", "*.xml"))
               if marca in os.path.basename(a).lower()]
        if not ann:
            continue
        base = os.path.splitext(os.path.basename(ann[0]))[0]
        tif = os.path.join(safe, "measurement", base + ".tiff")
        if not os.path.exists(tif):
            continue                      # .tiff assente o ancora .part
        cal = os.path.join(safe, "annotation", "calibration",
                           "calibration-" + base + ".xml")
        noi = os.path.join(safe, "annotation", "calibration",
                           "noise-" + base + ".xml")
        out.append(StackEntry(
            _data_da_nome(base), ann[0], tif,
            calibration=cal if os.path.exists(cal) else None,
            noise=noi if os.path.exists(noi) else None))
    return out


def _passo_orbitale(annotation: str) -> str:
    """Ascending / Descending letto dall'annotation, senza parsarlo tutto."""
    try:
        root = ET.parse(annotation).getroot()
    except Exception:                                   # pragma: no cover
        return "?"
    t = root.findtext(".//generalAnnotation/productInformation/pass")
    return (t or "?").strip()


def _traccia(annotation: str) -> str:
    """F43 -- traccia di appartenenza: verso di passaggio PIU' orbita relativa.

    Il verso da solo non basta. Due tracce ascendenti diverse guardano lo
    stesso punto da posizioni orbitali distanti centinaia di chilometri: la
    baseline ortogonale e' fuori scala esattamente come fra ascendente e
    discendente, e la coerenza e' zero. Il numero di orbita relativa non sta
    nell'annotation ma nel manifest.safe che gli sta accanto nell'albero
    .SAFE; se non c'e' (layout piatto) si ricade sul solo verso, che e' quanto
    faceva la versione precedente."""
    passo = _passo_orbitale(annotation)
    manifest = os.path.join(os.path.dirname(os.path.dirname(annotation)),
                            "manifest.safe")
    if not os.path.exists(manifest):
        return passo
    try:
        with open(manifest, "r", encoding="utf-8", errors="replace") as fh:
            testo = fh.read()
    except OSError:                                     # pragma: no cover
        return passo
    m = re.search(r"relativeOrbitNumber[^>]*>\s*(\d+)\s*<", testo)
    return f"{passo}/traccia {m.group(1)}" if m else passo


def _omogenea(entries: List[StackEntry], verbose: bool = True) -> List[StackEntry]:
    """Tiene solo le acquisizioni geometricamente compatibili fra loro.

    Una pila interferometrica deve stare su UNA sola traccia: ascendente e
    discendente guardano lo stesso punto da lati opposti dell'orbita, e la
    loro baseline ortogonale vale centinaia di chilometri -- cinque ordini di
    grandezza oltre la baseline critica, quindi coerenza esattamente zero.
    F43: il criterio e' verso di passaggio PIU' orbita relativa, perche' due
    tracce dello stesso verso sono incompatibili quanto due versi diversi.
    Mescolarle non degrada il risultato, lo distrugge. Si tiene il gruppo piu'
    numeroso e si dice a voce alta cosa e' stato scartato."""
    if len(entries) < 2:
        return entries
    gruppi: Dict[str, List[StackEntry]] = {}
    for e in entries:
        gruppi.setdefault(_traccia(e.annotation), []).append(e)
    if len(gruppi) == 1:
        return entries
    tenuto = max(gruppi, key=lambda k: len(gruppi[k]))
    if verbose:
        scartati = {k: len(v) for k, v in gruppi.items() if k != tenuto}
        print(f"  ATTENZIONE: la cartella contiene passaggi di tracce diverse. "
              f"Tengo {len(gruppi[tenuto])} acquisizioni {tenuto}, scarto "
              + ", ".join(f"{n} {k}" for k, n in scartati.items())
              + " (baseline oltre la critica: coerenza nulla).")
    return gruppi[tenuto]


def discover_stack(cfg: Config, verbose: bool = True) -> List[StackEntry]:
    """Trova le acquisizioni, in uno qualsiasi dei due layout supportati."""
    entries = _stack_piatto(cfg) or _omogenea(_stack_safe(cfg), verbose)
    if not entries:
        raise FileNotFoundError(
            f"nessuna acquisizione {cfg.swath}/{cfg.polarisation} in "
            f"{cfg.stack_dir}\n"
            f"  layout piatto cercato: {cfg.platform}-{cfg.swath}-slc-"
            f"{cfg.polarisation}-*.annotation.xml\n"
            f"  alberi cercati:        *.SAFE/annotation/*-{cfg.swath}-slc-"
            f"{cfg.polarisation}-*.xml")
    entries.sort(key=lambda e: e.date)
    return entries


#: alias storico
discover_vh_stack = discover_stack


def run(cfg: Config, verbose: bool = True) -> Dict[str, Any]:
    t_start = time.time()
    entries = discover_stack(cfg)[: cfg.n_dates]

    print(f"\nstack VH: {len(entries)} acquisizioni "
          f"({entries[0].date} -> {entries[-1].date})")
    print("\n  [1] accumulo interferometrico (lettura .tiff, deramping TOPS, "
          "coregistrazione, fase geometrica di riferimento)")
    cube = build_interf_cube(entries, cfg, verbose=verbose)
    n_l, n_p = cube.y.shape[1:]
    inc_mid = float(np.median(cube.incidence))
    az_m = n_p_lines = None
    az_m = cube.win.n_l * cube.ann.azimuth_pixel_spacing
    gr_m = cube.win.n_p * cube.ann.range_pixel_spacing / math.sin(math.radians(inc_mid))
    print(f"      cubo: {cube.y.shape} [date, linee, pixel]  "
          f"incidenza mediana {inc_mid:.2f} gradi")
    print(f"      area processata (F28): {az_m:.0f} m azimuth x {gr_m:.0f} m ground "
          f"range, margine {cfg.area_margin_m:.0f} m attorno alle piramidi"
          + ("  [SCENA INTERA]" if cfg.full_scene else ""))

    budget = compute_tomo_budget(cube.baselines, cube.ann, inc_mid,
                                 PYRAMIDS[0].height_m)
    print(budget.as_text())

    lay = layover_report(inc_mid)
    print(f"  geometria (F16): incidenza {lay['angolo_di_incidenza_deg']} gradi, "
          f"vista off-nadir {lay['angolo_di_vista_off_nadir_deg']} gradi; "
          + "; ".join(f"{r['piramide'].split()[0]} {r['faccia_vicina']}"
                      for r in lay["piramidi"]))

    # --- multilooking degli interferogrammi --------------------------------
    sl, sp = cfg.look_azimuth, cfg.look_range
    y_ml = multilook(cube.y, sl, sp).astype(np.complex64)
    k_ml = multilook(cube.k_z, sl, sp).astype(np.float32)
    amp_ml = multilook(cube.amp_master.astype(np.float32), sl, sp)
    # F40: il pavimento di rumore si toglie QUI, in intensita' e sul reticolo
    # multilooked. Il multilooking media sl*sp celle: riduce la varianza del
    # rumore ma non la sua MEDIA, che resta il NESZ e va sottratta. Farlo su
    # una media di 4 celle invece che pixel per pixel evita gli zeri che
    # rovinerebbero i prodotti a valle.
    nesz_ml = None
    if cube.nesz is not None:
        nesz_ml = multilook(cube.nesz.astype(np.float32), sl, sp)
        i_pulita = np.maximum(amp_ml.astype(np.float64) ** 2 - nesz_ml, 0.0)
        amp_ml = np.sqrt(i_pulita).astype(np.float32)
    inc_ml = multilook(cube.incidence, sl, sp)
    shape: Tuple[int, int] = (y_ml.shape[1], y_ml.shape[2])
    print(f"\n  [2] multilooking {sl}x{sp} -> reticolo {shape} "
          f"(~{cfg.look_azimuth * cube.ann.azimuth_pixel_spacing:.1f} m azimuth x "
          f"{sp * cube.ann.range_pixel_spacing / math.sin(math.radians(inc_mid)):.1f} m "
          f"ground range)")

    y_abs_sum = np.abs(y_ml).sum(axis=0).astype(np.float32)

    # --- simulazione delle piramidi in geometria radar (F24) ---------------
    print("\n  [3] proiezione delle piramidi nella geometria radar")
    sim = simulate_pyramids_radar(cube.geo, cube.ann, cube.win, cube.incidence)
    pyr_mask = multilook(sim["sim_mask"].astype(np.float32), sl, sp) > 0.5
    sim_h = multilook(sim["sim_h"], sl, sp).astype(np.float32)
    sim_fold = multilook(sim["sim_fold"].astype(np.float32), sl, sp)
    sim_per = [{"nome": d["nome"],
                "mask": multilook(d["mask"].astype(np.float32), sl, sp) > 0.4,
                "h": multilook(d["h"], sl, sp).astype(np.float32)}
               for d in sim["per"]]
    ground_mask = multilook(
        _pyramid_footprint_mask(cube.geo, cube.win, 1.0).astype(np.float32),
        sl, sp) > 0.5
    dh_apex = float(sim["sim_h"].max())
    shift_px = dh_apex * math.cos(math.radians(inc_mid)) / cube.ann.range_pixel_spacing
    print(f"      apice piu' alto sopra il riferimento: {dh_apex:.1f} m -> "
          f"spostamento in slant {shift_px:.1f} pixel "
          f"({shift_px / sp:.1f} celle multilooked)")
    print(f"      celle raggiunte dalle piramidi: {int(pyr_mask.sum())} "
          f"(impronta al suolo: {int(ground_mask.sum())}, "
          f"sovrapposizione {100.0 * (pyr_mask & ground_mask).sum() / max(pyr_mask.sum(), 1):.0f} %)")
    sim_fold_max = float(sim["sim_fold"].max())
    print(f"      ripiegamento massimo per cella (layover): "
          f"{float(sim_fold_max):.0f} punti di superficie")

    print("\n      contrasto radiometrico delle piramidi per polarizzazione (F26):")
    pol_contrast = polarisation_contrast(cfg, verbose=True)
    if "conclusione" in pol_contrast:
        for line in _wrap(pol_contrast["conclusione"], 68):
            print(f"      {line}")

    # --- calibrazione del segno di k_z (F05) -------------------------------
    print("\n  [4] calibrazione del segno di k_z sui dati (validazione livello 1)")
    z_axis = np.linspace(-cfg.elev_max_m, cfg.elev_max_m, cfg.n_elev).astype(np.float32)
    high_mask = pyr_mask & (sim_h >= 0.5 * float(sim_h.max()))
    sign, sign_report = calibrate_kz_sign(y_ml, k_ml, z_axis, high_mask,
                                          y_abs_sum, verbose=verbose)

    # --- inversione tomografica --------------------------------------------
    print("\n  [5] inversione tomografica multi-baseline (periodogramma in quota)")
    tomo = tomographic_periodogram(y_ml, k_ml, z_axis, sign=sign)
    tomo_mag = np.abs(tomo).astype(np.float32)

    # --- superficie reale (F17) --------------------------------------------
    print("  [6] superficie reale: quota del diffusore dominante pixel per pixel")
    dh, gamma = surface_from_tomogram(tomo, z_axis, y_abs_sum)
    east, north, h_ref, (lat0, lon0) = local_enu(cube.geo, cube.win, sl, sp, shape)
    surf_abs = (h_ref + dh).astype(np.float32)

    # F35: i nodi grezzi della geolocation grid, per dimostrare che il quasi
    # piano di "height_ref" dentro il ritaglio non e' un errore di calcolo ma
    # il limite di risoluzione del riferimento .xml (vedi raw_gcp_nodes()).
    gcp_raw = raw_gcp_nodes(cube.geo, lat0, lon0)
    gcp_e = np.asarray(gcp_raw["east"])
    gcp_n = np.asarray(gcp_raw["north"])
    gcp_h = np.asarray(gcp_raw["h"])
    d_gcp = np.hypot(gcp_e - float(np.mean(east)), gcp_n - float(np.mean(north)))
    i_near = int(np.argmin(d_gcp))
    gcp_line_sp = int(np.diff(np.asarray(cube.geo.lines)).min())
    gcp_pix_sp = int(np.diff(np.asarray(cube.geo.pixels)).min())
    print(f"      reticolo GREZZO dello .xml (F35): {gcp_raw['n_l']} x "
          f"{gcp_raw['n_p']} nodi su tutta la scena, spaziatura "
          f"{gcp_line_sp} linee x {gcp_pix_sp} pixel -- il ritaglio delle "
          f"piramidi ({shape[0] * sl} x {shape[1] * sp} pixel) sta dentro UNA "
          f"frazione di una sola cella di quel reticolo.")
    print(f"      nodo .xml piu' vicino al ritaglio: a {d_gcp[i_near] / 1000:.1f} km, "
          f"quota {gcp_h[i_near]:.0f} m. Rilievo VERO sull'intera scena: "
          f"{gcp_h.min():.0f} .. {gcp_h.max():.0f} m (dev.std {gcp_h.std():.0f} m) "
          "-- i nodi grezzi non sono un piano, ma nessuno cade abbastanza "
          "vicino al ritaglio da farlo vedere nella spline locale.")

    # F38/F39: suolo (DEM). NON entra nel calcolo, solo un layer di riferimento
    # visivo. F39: il datum sono le quote lette dagli annotation.xml dello
    # stack, e il layer comprende il profilo delle piramidi.
    suolo_dem = None
    suolo_info: Dict[str, Any] = {}
    plateau = None
    if cfg.suolo_dem:
        print("  [6b] suolo (DEM): terreno + profilo delle piramidi (F39)")
        plateau = plateau_heights_from_stack(cfg, cube.geo, verbose=verbose)
        suolo_dem, suolo_info = ground_dem_suolo(
            lat0, lon0, east.astype(np.float64), north.astype(np.float64),
            h_ref.astype(np.float64), cfg.out_dir, gc=cfg.dem_grid,
            cache_name=cfg.dem_cache_name, use_external=cfg.fetch_dem,
            add_pyramids=cfg.dem_pyramids, verbose=verbose)
        print(f"      escursione sul ritaglio: {float(suolo_dem.max() - suolo_dem.min()):.0f} m "
              f"(contro {float(h_ref.max() - h_ref.min()):.1f} m della bilineare .xml, F36)")

    # F22: la soglia esce dalla distribuzione nulla, non da una costante
    k_typ = np.median(k_ml.reshape(k_ml.shape[0], -1), axis=1) * sign
    thr, null_stats = null_gamma_threshold(k_typ, z_axis, q=cfg.gamma_null_q,
                                           trials=cfg.gamma_null_trials,
                                           seed=cfg.seed)
    thr = max(thr, cfg.gamma_min)
    good = gamma >= thr
    print(f"      quota di riferimento dallo .xml: "
          f"{float(h_ref.min()):.1f} .. {float(h_ref.max()):.1f} m")
    print(f"      distribuzione NULLA del periodogramma (F22): mediana "
          f"{null_stats['mediana_nulla']}, p99 {null_stats['p99_nulla']}")
    print(f"      soglia di qualita' = {thr:.3f} (percentile {cfg.gamma_null_q} del nullo)")
    print(f"      pixel di qualita': {int(good.sum())} / {good.size} "
          f"({100.0 * good.mean():.1f} %)  -- gamma mediana misurata "
          f"{float(np.median(gamma)):.3f}")
    if float(np.median(gamma)) < null_stats["mediana_nulla"]:
        print("      NOTA: la gamma mediana e' sotto la mediana del nullo: la "
              "maggioranza dei pixel non porta informazione di quota.")

    # --- attributi ----------------------------------------------------------
    print("\n  [7] coerenza interferometrica")
    coh = stack_coherence(y_ml, amp_ml, window=cfg.coh_window,
                          master_date=cube.master_date, dates=cube.dates)
    gamma_meas = float(np.median(coh))
    budget_a_priori = budget
    budget = compute_tomo_budget(cube.baselines, cube.ann, inc_mid,
                                 PYRAMIDS[0].height_m, gamma_typ=gamma_meas)
    print(f"      coerenza mediana MISURATA = {gamma_meas:.3f} "
          f"(a priori si era assunto {budget_a_priori.gamma_typ:.2f})")
    print(f"      budget di precisione aggiornato: sigma_h = {budget.sigma_h:.1f} m "
          f"(era {budget_a_priori.sigma_h:.1f} m con la coerenza assunta)")
    if not budget.surface_measurable:
        print("      ATTENZIONE: con la coerenza misurata la superficie NON e' "
              "misurabile secondo il criterio sigma_h <= h/4. Il risultato va "
              "letto come limite superiore, non come misura.")

    surf_disp, n_spike = despike_surface(surf_abs, good, budget.sigma_h)
    print(f"      despicatura (F27): {n_spike} nodi di qualita' sostituiti con la "
          f"mediana locale ({100.0 * n_spike / max(int(good.sum()), 1):.1f} % dei "
          f"nodi sopra soglia)")

    print("  [8] micro-moto dalle sub-aperture Doppler (chip esteso, F07)")
    mm_win, mm_burst = target_window(cube.ann, cube.geo, cfg, n_lines=cfg.mm_lines)
    plan = plan_subapertures(cube.ann, mm_win.n_l, cfg, burst_idx=mm_burst,
                             pixel=0.5 * (mm_win.p0 + mm_win.p1))
    print(plan.as_text())
    mm_entry = next(e for e in entries if e.date == cube.master_date)
    mm_chip = read_window(mm_entry, mm_win, mm_burst)
    mm_img = tops_deramp(mm_chip, cube.ann)
    row_off = cube.win.l0 - mm_win.l0
    mm_full, mm_coh_full, mm_freq_full = micro_motion_energy(
        mm_img, cube.ann, plan, cfg, (n_l, n_p), row_off, verbose=verbose)
    mm = multilook(mm_full, sl, sp)
    mm_coh = multilook(mm_coh_full, sl, sp)
    # F31: la frequenza dominante e' una ETICHETTA su una griglia discreta di
    # righe (multipli di 1/t_window): mediarla su un blocco di multilooking
    # produrrebbe valori che il banco non puo' generare, quindi si prende il
    # campione centrale del blocco.
    mm_freq = mm_freq_full[(sl - 1) // 2::sl, (sp - 1) // 2::sp][:n_l // sl, :n_p // sp]
    mm_freq = np.ascontiguousarray(mm_freq)
    n_in_band = int(np.count_nonzero(
        (mm_freq >= plan.f_min_obs) & (mm_freq <= plan.f_max_obs)))
    print(f"      micro-moto: concentrazione spettrale mediana "
          f"{float(np.median(mm_coh)):.3f} (rumore bianco = {1.0 / (plan.n_d - 1):.3f}), "
          f"riga dominante entro [{plan.f_min_obs:.1f}, {plan.f_max_obs:.1f}] Hz "
          f"su {100.0 * n_in_band / max(mm_freq.size, 1):.0f}% delle celle")

    print("  [9] discriminante multi-attributo pieno/vuoto")
    solidity = solidity_index(tomo_mag, coh, mm)

    print(f"\n  volume finale: {tomo_mag.shape} [azimuth, range, elevazione]")

    return {
        "tomo_mag": tomo_mag,
        "solidity": solidity,
        "coherence": coh,
        "gamma": gamma,
        "height_rel": dh,
        "height_abs": surf_abs,
        "height_display": surf_disp,
        "n_spike": n_spike,
        "height_ref": h_ref,
        "gcp_raw": gcp_raw,
        "gcp_raw_near_km": round(float(d_gcp[i_near]) / 1000.0, 2),
        "gcp_raw_near_h": round(float(gcp_h[i_near]), 1),
        "gcp_raw_spacing_px": (gcp_line_sp, gcp_pix_sp),
        "suolo_dem": suolo_dem,
        "suolo_info": suolo_info,
        "suolo_plateau": plateau,
        "good": good,
        "amp": amp_ml,
        "mm": mm,
        "mm_coh": mm_coh,
        "mm_freq": mm_freq,
        "mm_band": (plan.f_min_obs, plan.f_max_obs),
        "mm_noise_floor": 1.0 / max(plan.n_d - 1, 1),
        "incidence": inc_ml,
        "z_axis": z_axis,
        "east": east,
        "north": north,
        "lat0": lat0,
        "lon0": lon0,
        "pyr_mask": pyr_mask,
        "area_m": (az_m, gr_m),
        "ground_mask": ground_mask,
        "sim_h": sim_h,
        "sim_fold": sim_fold,
        "sim_per": sim_per,
        "pol_contrast": pol_contrast,
        "gamma_thr": thr,
        "null_stats": null_stats,
        "budget": budget,
        "budget_a_priori": budget_a_priori,
        "gamma_misurata": gamma_meas,
        "layover": lay,
        "sign": sign,
        "sign_report": sign_report,
        "plan": plan,
        "baselines": cube.baselines,
        "coreg": cube.coreg,
        "dates": cube.dates,
        "master_date": cube.master_date,
        "cfg": cfg,
        "elapsed_s": time.time() - t_start,
    }


# ==========================================================================
# 11.  Validazione a tre livelli (ch15)
# ==========================================================================

def validate(res: Dict[str, Any], cfg: Config) -> Dict[str, Any]:
    """Protocollo a tre livelli: geometria, misura, struttura.

    Ogni livello puo' fallire senza gli altri; solo tutti e tre insieme
    vincolano la catena. I disaccordi vanno riportati quanto gli accordi."""
    h = res["height_abs"]
    g = res["gamma"]
    east, north = res["east"], res["north"]
    good = res["good"]
    lat0, lon0 = res["lat0"], res["lon0"]
    budget: TomoBudget = res["budget"]

    dh = res["height_rel"]
    sim_h = res["sim_h"]
    pyr = res["pyr_mask"]

    # --- livello 1: regressione misurato contro simulato in geometria radar --
    # F25: il confronto giusto non e' un percentile della quota dentro
    # l'impronta al suolo, ma la regressione fra la quota misurata e quella
    # simulata cella per cella nella geometria in cui il radar guarda davvero.
    lvl1: Dict[str, Any] = {}
    plane = good & (~pyr)
    h_plane = float(np.median(h[plane])) if plane.sum() > 20 else float("nan")
    dh_plane = float(np.median(dh[plane])) if plane.sum() > 20 else 0.0

    sel = good & pyr & (sim_h > 5.0)
    if sel.sum() >= 20:
        x = sim_h[sel].astype(np.float64)
        yv = (dh[sel] - dh_plane).astype(np.float64)
        slope, icept = np.polyfit(x, yv, 1)
        r = float(np.corrcoef(x, yv)[0, 1])
        resid = yv - (slope * x + icept)
        # significativita' della pendenza (errore standard OLS)
        se = float(np.std(resid, ddof=2) /
                   max(np.sqrt(np.sum((x - x.mean()) ** 2)), 1e-9))
        # stima robusta: pochi errori grossolani da centinaia di metri
        # trascinano una regressione ai minimi quadrati anche dopo la
        # despicatura, quindi la pendenza viene riportata anche alla Theil-Sen
        try:
            from scipy.stats import theilslopes
            ts = theilslopes(yv, x, 0.95)
            robusto = {"pendenza_theil_sen": round(float(ts[0]), 3),
                       "ic95": [round(float(ts[2]), 3), round(float(ts[3]), 3)]}
        except Exception:
            robusto = {}
        # mediane per fascia di quota simulata: dicono se il rilievo e' seguito
        # o se sopra la piramide c'e' solo un offset costante
        fasce = []
        for lo, hi in ((5, 30), (30, 60), (60, 90), (90, 1000)):
            b = (x >= lo) & (x < hi)
            if b.sum() >= 5:
                fasce.append({"simulato_m": [lo, hi if hi < 1000 else None],
                              "celle": int(b.sum()),
                              "misurato_mediano_m": round(float(np.median(yv[b])), 1)})
        lvl1["regressione_misurato_vs_simulato"] = {
            "pixel": int(sel.sum()),
            "pendenza": round(float(slope), 3),
            **robusto,
            "mediane_per_fascia": fasce,
            "pendenza_attesa": 1.0,
            "errore_standard_pendenza": round(se, 3),
            "t_pendenza": round(float(slope / max(se, 1e-9)), 1),
            "intercetta_m": round(float(icept), 1),
            "correlazione_r": round(r, 3),
            "rmse_m": round(float(np.std(resid)), 1),
            "esito": ("la quota misurata segue quella simulata"
                      if slope > 3 * se and r > 0.15 else
                      "nessuna dipendenza significativa dalla quota simulata"),
        }
    else:
        lvl1["regressione_misurato_vs_simulato"] = {
            "pixel": int(sel.sum()), "esito": "pixel di qualita' insufficienti"}

    for d in res["sim_per"]:
        m, sh = d["mask"], d["h"]
        if not m.any():
            lvl1[d["nome"]] = {"esito": "fuori copertura"}
            continue
        top = m & good & (sh >= 0.75 * float(sh[m].max()))
        base = m & good & (sh <= 0.25 * float(sh[m].max()))
        rec: Dict[str, Any] = {
            "celle_simulate": int(m.sum()),
            "celle_di_qualita": int((m & good).sum()),
            "celle_di_sommita": int(top.sum()),
        }
        if top.sum() >= 5:
            mis = float(np.median(dh[top]) - dh_plane)
            att = float(np.median(sh[top]))
            rec.update({
                "rilievo_misurato_m": round(mis, 1),
                "rilievo_simulato_m": round(att, 1),
                "errore_m": round(mis - att, 1),
                "entro_2_sigma_h": bool(abs(mis - att) <= 2.0 * budget.sigma_h),
            })
        if base.sum() >= 5:
            rec["rilievo_misurato_alla_base_m"] = round(
                float(np.median(dh[base]) - dh_plane), 1)
        lvl1[d["nome"]] = rec

    lvl1["_quota_mediana_piana_m"] = round(h_plane, 1) if np.isfinite(h_plane) else None
    lvl1["_sigma_h_attesa_m"] = round(budget.sigma_h, 1)
    lvl1["_soglia_qualita_gamma"] = round(float(res["gamma_thr"]), 3)
    lvl1["_distribuzione_nulla"] = res["null_stats"]
    lvl1["_riferimento"] = ("le quote note da letteratura archeologica entrano SOLO "
                            "nella simulazione di confronto, mai nella catena che "
                            "produce la superficie misurata")

    # --- livello 2: la misura contro se stessa ------------------------------
    #   dispersione della quota sulla piana, che dovrebbe essere piatta
    lvl2 = {
        "pixel_buoni_frazione": round(float(good.mean()), 4),
        "gamma_mediana": round(float(np.median(g)), 4),
        "gamma_p90": round(float(np.percentile(g, 90)), 4),
        "dispersione_quota_sulla_piana_m": (
            round(float(np.std(h[plane])), 1) if plane.sum() > 20 else None),
        "sigma_h_teorica_m": round(budget.sigma_h, 1),
        "nota": ("la dispersione misurata sulla piana e' il controllo empirico "
                 "della sigma_h teorica: se le due sono dello stesso ordine il "
                 "budget di precisione e' onesto"),
    }

    # --- livello 3: struttura contro geometria indipendente -----------------
    lvl3 = {
        "layover": res["layover"],
        "separabilita_piramide_vs_piana": None,
        "avvertenza_interna": (
            f"delta_z verticale {budget.delta_z_vertical:.0f} m contro "
            f"{budget.target_height:.0f} m di altezza: la piramide occupa "
            f"{budget.cells_over_target:.2f} celle di RISOLUZIONE. Nessuna "
            "struttura interna e' separabile con questi dati. Le camere note "
            "(Grande Galleria, camera del Re, camera della Regina) hanno "
            "dimensioni di alcuni metri, due ordini di grandezza sotto."),
    }
    if plane.sum() > 20:
        on = good & res["pyr_mask"]
        if on.sum() > 20:
            d = float(np.median(h[on]) - h_plane)
            s = float(np.sqrt(np.var(h[on]) / max(on.sum(), 1)
                              + np.var(h[plane]) / max(plane.sum(), 1)))
            lvl3["separabilita_piramide_vs_piana"] = {
                "delta_quota_mediana_m": round(d, 1),
                "errore_standard_m": round(s, 2),
                "z_score": round(d / max(s, 1e-6), 1),
                "pixel_su_piramide": int(on.sum()),
                "pixel_su_piana": int(plane.sum()),
            }

    return {
        "livello_1_geometria": lvl1,
        "livello_2_misura": lvl2,
        "livello_3_struttura": lvl3,
        "calibrazione_segno_k_z": res["sign_report"],
    }


def profile_analysis(res: Dict[str, Any], cfg: Config) -> Dict[str, Any]:
    """Profili verticali sotto ciascuna piramide contro una colonna di deserto.

    F11: la media e' fatta in POTENZA e convertita in dB alla fine.
    F12: la colonna di controllo e' verificata non vuota."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tomo, sol = res["tomo_mag"], res["solidity"]
    z, east, north = res["z_axis"], res["east"], res["north"]
    pw = np.maximum(tomo.astype(np.float64), 1e-12) ** 2
    lat0, lon0 = res["lat0"], res["lon0"]

    def column(e0, n0, half):
        m = (np.abs(east - e0) <= half) & (np.abs(north - n0) <= half)
        if m.sum() == 0:
            return None, None, 0
        p = pw[m].mean(axis=0)
        d = 10.0 * np.log10(np.maximum(p, 1e-12))
        return d, sol[m].mean(axis=0), int(m.sum())

    meshes = [pyramid_mesh(p, lat0, lon0) for p in PYRAMIDS]
    ref_db = None
    stats: Dict[str, Any] = {}

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), constrained_layout=True)
    curves = []
    for m, p in zip(meshes, PYRAMIDS):
        d, s, npx = column(m["apex"][0], m["apex"][1], p.base_side_m * 0.5)
        if d is None:
            stats[p.name] = {"pixel_nella_colonna": 0, "esito": "fuori copertura"}
            continue
        curves.append((p.name, d, s, npx))
        ref_db = d.max() if ref_db is None else max(ref_db, d.max())

    # F12: colonna di controllo dentro la copertura effettiva, verificata
    e_lo, e_hi = float(east.min()), float(east.max())
    n_lo, n_hi = float(north.min()), float(north.max())
    apex_e = [m["apex"][0] for m in meshes]
    apex_n = [m["apex"][1] for m in meshes]
    best = None
    for frac_e in np.linspace(0.08, 0.92, 15):
        for frac_n in np.linspace(0.08, 0.92, 15):
            e0 = e_lo + frac_e * (e_hi - e_lo)
            n0 = n_lo + frac_n * (n_hi - n_lo)
            dist = min(math.hypot(e0 - a, n0 - b) for a, b in zip(apex_e, apex_n))
            d, s, npx = column(e0, n0, 120.0)
            if npx >= 20 and (best is None or dist > best[0]):
                best = (dist, e0, n0, d, s, npx)
    if best is None:
        raise RuntimeError("nessuna colonna di deserto valida nel chip (F12)")
    _, bg_e, bg_n, d_bg, s_bg, npx_bg = best
    ref_db = d_bg.max() if ref_db is None else max(ref_db, d_bg.max())

    for name, d, s, npx in curves:
        axes[0].plot(d - ref_db, z, lw=2, label=f"{name} ({npx} px)")
        axes[1].plot(s, z, lw=2, label=name)
        stats[name] = {
            "pixel_nella_colonna": npx,
            "intensita_media_dB": round(float(10 * np.log10(
                np.maximum(np.mean(10 ** (d / 10)), 1e-12)) - ref_db), 2),
            "solidita_media": round(float(np.mean(s)), 4),
            "solidita_max": round(float(np.max(s)), 4),
            "quota_del_massimo_m": round(float(z[int(np.argmax(d))]), 1),
        }
    axes[0].plot(d_bg - ref_db, z, "--", color="#94A3B8", lw=1.8,
                 label=f"deserto di controllo ({npx_bg} px)")
    axes[1].plot(s_bg, z, "--", color="#94A3B8", lw=1.8, label="deserto di controllo")
    stats["deserto_controllo"] = {
        "pixel_nella_colonna": npx_bg,
        "posizione_enu_m": [round(bg_e, 1), round(bg_n, 1)],
        "intensita_media_dB": round(float(10 * np.log10(
            np.maximum(np.mean(10 ** (d_bg / 10)), 1e-12)) - ref_db), 2),
        "solidita_media": round(float(np.mean(s_bg)), 4),
    }

    # terzo pannello: istogramma delle quote misurate
    h = res["height_abs"][res["good"]]
    if h.size:
        axes[2].hist(h, bins=60, color="#16A34A", alpha=.75)
        axes[2].set_xlabel("quota assoluta misurata [m]")
        axes[2].set_ylabel("pixel")
        axes[2].set_title("Distribuzione delle quote (pixel coerenti)")
        for p in PYRAMIDS:
            axes[2].axvline(p.base_alt_m + p.height_m, color="#F0A24A", ls=":", lw=1.4)
        # F32: la banda delle quote di riferimento lette dagli .xml, per
        # vedere dove cade l'istogramma rispetto al terreno del prodotto
        hr = res["height_ref"]
        axes[2].axvspan(float(hr.min()), float(hr.max()), color="#94A3B8",
                        alpha=.35, zorder=0,
                        label=f"quote .xml ({hr.min():.1f}-{hr.max():.1f} m)")
        axes[2].legend(fontsize=8)
        axes[2].grid(alpha=.25)

    for ax, (xl, ti) in zip(axes[:2],
                            [("intensita' relativa [dB]", "Profilo verticale di intensita'"),
                             ("indice di solidita'", "Profilo verticale di solidita'")]):
        ax.axhline(0, color="#888", lw=.8, ls=":")
        ax.set_xlabel(xl)
        ax.set_ylabel("quota relativa [m]")
        ax.set_title(ti)
        ax.legend(fontsize=8)
        ax.grid(alpha=.25)
    fig.suptitle("Colonne verticali sotto le piramidi contro deserto di controllo  "
                 f"(risoluzione verticale {res['budget'].delta_z_vertical:.0f} m, "
                 f"precisione sulla quota {res['budget'].sigma_h:.0f} m)")
    fig.savefig(os.path.join(cfg.out_dir, "profili_pieno_vuoto.png"), dpi=140)
    plt.close(fig)

    sep = {}
    bgi = stats["deserto_controllo"]["intensita_media_dB"]
    bgs = stats["deserto_controllo"]["solidita_media"]
    for p in PYRAMIDS:
        st = stats.get(p.name, {})
        if "intensita_media_dB" not in st:
            continue
        sep[p.name] = {
            "delta_intensita_dB": round(st["intensita_media_dB"] - bgi, 2),
            "delta_solidita": round(st["solidita_media"] - bgs, 4),
        }
    stats["separabilita_vs_deserto"] = sep

    print("\n  --- colonne verticali: piramidi contro deserto ---")
    for p in PYRAMIDS:
        st = stats.get(p.name, {})
        if "intensita_media_dB" not in st:
            print(f"    {p.name:24s} (fuori copertura)")
            continue
        print(f"    {p.name:24s} int={st['intensita_media_dB']:7.2f} dB  "
              f"sol={st['solidita_media']:.4f}  picco a z={st['quota_del_massimo_m']:+7.1f} m")
    bg = stats["deserto_controllo"]
    print(f"    {'deserto di controllo':24s} int={bg['intensita_media_dB']:7.2f} dB  "
          f"sol={bg['solidita_media']:.4f}  ({bg['pixel_nella_colonna']} px)")
    return stats


def surface_plot(res: Dict[str, Any], cfg: Config,
                 valid: Dict[str, Any]) -> None:
    """Tre pannelli di ispezione: superficie, qualita', validazione di livello 1.

    I limiti di colore sono percentili robusti: con gli estremi la scala viene
    schiacciata dagli errori grossolani residui e il rilievo reale sparisce."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = np.where(res["good"], res["height_display"], np.nan)
    east, north = res["east"], res["north"]
    fin = np.isfinite(h)
    lo, hi = (np.percentile(h[fin], [2, 98]) if fin.any() else (0.0, 1.0))

    fig, ax = plt.subplots(1, 3, figsize=(19, 5.8), constrained_layout=True)

    # F30: un punto per cella misurata, non un campo continuo di quadrati.
    # pcolormesh disegna comunque una tessera piena per ogni cella e, con il
    # 23 per cento di celle sopra soglia, dava alla mappa l'aspetto di una
    # superficie campionata fitto. Il punto dice dove c'e' una misura.
    ax[0].set_facecolor("#F3F4F6")
    im = ax[0].scatter(east[fin], north[fin], c=h[fin], s=7, marker="o",
                       linewidths=0, cmap="terrain", vmin=lo, vmax=hi)
    fig.colorbar(im, ax=ax[0], label="quota assoluta [m]")
    ax[0].set_xlim(float(east.min()), float(east.max()))
    ax[0].set_ylim(float(north.min()), float(north.max()))
    ax[0].set_title(f"Superficie dai .tiff {cfg.polarisation.upper()} "
                    f"({int(fin.sum())} celle sopra soglia, scala p2-p98)")

    im2 = ax[1].pcolormesh(east, north, res["gamma"], shading="auto",
                           cmap="magma", vmin=0.4, vmax=1.0)
    fig.colorbar(im2, ax=ax[1], label="coerenza di fit")
    ax[1].set_title(f"Qualita' (soglia dal nullo = {res['gamma_thr']:.3f})")

    for a in ax[:2]:
        a.contour(east, north, res["pyr_mask"].astype(float), levels=[0.5],
                  colors="#38BDF8", linewidths=1.2)
        for p in PYRAMIDS:
            m = pyramid_mesh(p, res["lat0"], res["lon0"])
            xs = [v[0] for v in m["vertices"][:4]] + [m["vertices"][0][0]]
            ys = [v[1] for v in m["vertices"][:4]] + [m["vertices"][0][1]]
            a.plot(xs, ys, color="#F0A24A", lw=1.2)
        a.set_xlabel("est [m]")
        a.set_ylabel("nord [m]")
        a.set_aspect("equal")

    # pannello 3: la validazione di livello 1, misurato contro simulato
    sel = res["good"] & res["pyr_mask"] & (res["sim_h"] > 5.0)
    plane = res["good"] & (~res["pyr_mask"])
    dh0 = float(np.median(res["height_rel"][plane])) if plane.sum() > 20 else 0.0
    reg = valid["livello_1_geometria"]["regressione_misurato_vs_simulato"]
    if sel.sum() >= 20:
        x = res["sim_h"][sel]
        y = res["height_rel"][sel] - dh0
        ax[2].scatter(x, y, s=9, alpha=.45, color="#38BDF8", edgecolors="none",
                      label=f"celle di qualita' ({int(sel.sum())})")
        xr = np.linspace(0, float(x.max()), 10)
        ax[2].plot(xr, xr, "--", color="#F0A24A", lw=1.6, label="atteso (pendenza 1)")
        if "pendenza" in reg:
            ax[2].plot(xr, reg["pendenza"] * xr + reg["intercetta_m"], "-",
                       color="#16A34A", lw=2,
                       label=f"misurato (pendenza {reg['pendenza']:+.3f} "
                             f"+- {reg['errore_standard_pendenza']:.3f})")
        ax[2].axhline(0, color="#888", lw=.8, ls=":")
        ax[2].legend(fontsize=8, loc="upper left")
    ax[2].set_xlabel("quota simulata in geometria radar [m]")
    ax[2].set_ylabel("quota misurata sopra la piana [m]")
    ax[2].set_title("Validazione livello 1: misurato vs simulato")
    ax[2].grid(alpha=.25)

    fig.suptitle(
        f"{cfg.polarisation.upper()} - {len(res['dates'])} date, master "
        f"{res['master_date']} | sigma_h {res['budget'].sigma_h:.0f} m, "
        f"delta_z {res['budget'].delta_z_vertical:.0f} m, "
        f"{100 * float(res['good'].mean()):.0f}% celle sopra soglia")
    fig.savefig(os.path.join(cfg.out_dir, "superficie_ricostruita.png"), dpi=130)
    plt.close(fig)


# ==========================================================================
# 12.  Uscite
# ==========================================================================

def save_outputs(res: Dict[str, Any], cfg: Config,
                 valid: Dict[str, Any], profiles: Dict[str, Any]) -> None:
    os.makedirs(cfg.out_dir, exist_ok=True)
    for key in ("tomo_mag", "solidity", "coherence", "gamma", "z_axis",
                "height_abs", "height_display", "height_rel", "height_ref",
                "east", "north",
                "amp", "mm", "mm_coh", "mm_freq", "sim_h", "sim_fold"):
        np.save(os.path.join(cfg.out_dir, f"{key}.npy"), res[key])
    np.save(os.path.join(cfg.out_dir, "pyr_mask.npy"), res["pyr_mask"])
    np.save(os.path.join(cfg.out_dir, "good.npy"), res["good"])

    budget: TomoBudget = res["budget"]
    plan: MMPlan = res["plan"]
    meta = {
        "generato": time.strftime("%Y-%m-%d %H:%M:%S"),
        "programma": "piramidi_v02.py (revisione 2026-08-28)",
        "metodo": ("stack interferometrico multi-baseline VH + periodogramma in "
                   "quota + attributi di micro-moto Doppler"),
        "polarizzazione": "VH (cross-pol, scattering di volume)",
        "sorgente_geometria": "annotation.xml (orbite, geolocation grid, quota di riferimento)",
        "sorgente_misura": "tiff VH, pixel per pixel",
        "date": res["dates"],
        "master": res["master_date"],
        "coregistrazione": res["coreg"],
        "budget": asdict(budget),
        "budget_a_priori": asdict(res["budget_a_priori"]),
        "coerenza_mediana_misurata": round(float(res["gamma_misurata"]), 4),
        "soglia_qualita_gamma": round(float(res["gamma_thr"]), 4),
        "distribuzione_nulla_periodogramma": res["null_stats"],
        "nodi_despicati": int(res["n_spike"]),
        "contrasto_per_polarizzazione": res["pol_contrast"],
        "piano_sub_aperture": asdict(plan),
        "segno_k_z": res["sign"],
        "layover": res["layover"],
        "validazione": valid,
        "analisi_colonne": profiles,
        "geolocation_grid_xml": {
            "descrizione": (
                "F35/F36: i nodi VERI (non interpolati) della geolocation "
                "grid dello .xml. Non sono un piano: coprono tutta la scena "
                "con un rilievo reale fino a centinaia di metri. height_ref/"
                "xml_ref sono ora la bilineare LOCALE (F36) fra i quattro "
                "nodi reali della cella che contiene il ritaglio -- non una "
                "spline globale sui 231 nodi -- ma la spaziatura del "
                "reticolo e' enormemente piu' larga del ritaglio delle "
                "piramidi, quindi quella cella non ha rilievo locale da "
                "mostrare, e il risultato e' correttamente quasi piatto: "
                "e' un limite di risoluzione del riferimento, non un errore "
                "di calcolo."
            ),
            "spaziatura_reticolo_px": list(res["gcp_raw_spacing_px"]),
            "nodo_piu_vicino_al_ritaglio_km": res["gcp_raw_near_km"],
            "nodo_piu_vicino_al_ritaglio_quota_m": res["gcp_raw_near_h"],
            "rilievo_vero_su_tutta_la_scena_m": [
                round(float(min(res["gcp_raw"]["h"])), 1),
                round(float(max(res["gcp_raw"]["h"])), 1),
            ],
        },
        "correzioni_applicate": [
            {"id": i, "gravita": g, "descrizione": d} for i, g, d in FIXES
        ],
        "avvertenza": (
            "La superficie e' una misura: quota di riferimento dallo .xml piu' "
            "quota del diffusore dominante stimata dai .tiff. La sua precisione e' "
            f"~{budget.sigma_h:.0f} m. La RISOLUZIONE verticale e' "
            f"{budget.delta_z_vertical:.0f} m: nessuna struttura interna e' "
            "separabile. L'indice di solidita' e' un discriminante "
            "multi-attributo, non una rilevazione di cavita'."
        ),
    }
    with open(os.path.join(cfg.out_dir, "meta.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False, default=str)
    with open(os.path.join(cfg.out_dir, "budget_tomografico.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(budget.as_text() + "\n\n" + plan.as_text() + "\n")


def _decimate_to(n: int, target: int) -> np.ndarray:
    step = max(1, int(math.ceil(n / max(target, 1))))
    return np.arange(0, n, step)


def build_html(res: Dict[str, Any], cfg: Config,
               valid: Dict[str, Any]) -> str:
    """Pagina 3D autonoma: superficie reale + nuvola tomografica + riferimenti."""
    tomo = res["tomo_mag"]
    sol = res["solidity"]
    z = res["z_axis"]
    east, north = res["east"], res["north"]
    budget: TomoBudget = res["budget"]

    db = 20 * np.log10(np.maximum(tomo, 1e-12))
    db -= db.max()

    # --- superficie: reticolo decimato ------------------------------------
    n_l, n_p = res["height_abs"].shape
    target_side = int(math.sqrt(cfg.surface_max_nodes * n_l / max(n_p, 1)))
    si = _decimate_to(n_l, max(target_side, 8))
    sj = _decimate_to(n_p, max(int(cfg.surface_max_nodes / max(len(si), 1)), 8))
    gi, gj = np.meshgrid(si, sj, indexing="ij")

    h_s = res["height_display"][gi, gj]
    g_s = res["gamma"][gi, gj]
    good_s = res["good"][gi, gj]
    # normalizzata una volta sola: serve sia alla nuvola decimata sia al
    # reticolo xml_ref, e ricalcolarla due volte dava per forza lo stesso
    # risultato (i percentili sono quelli dell'intera matrice in entrambi i casi)
    mm_n = _norm01(res["mm"])
    # F30: i nodi non affidabili prendono la quota di riferimento SOLO per
    # avere un valore finito da proiettare e da usare nel calcolo delle
    # normali; nella nuvola escono come punti piccoli e spenti, e il pulsante
    # "Solo celle di qualita'" li toglie del tutto. Non riempiono piu' il
    # vuoto dentro un poligono facendosi passare per misura.
    h_s = np.where(good_s | np.isfinite(h_s), h_s, res["height_ref"][gi, gj])
    amp_db = 20 * np.log10(np.maximum(res["amp"][gi, gj], 1e-6))
    amp_n = (amp_db - np.percentile(amp_db, 2)) / \
        max(float(np.percentile(amp_db, 98) - np.percentile(amp_db, 2)), 1e-6)

    surface = {
        "n_l": int(len(si)), "n_p": int(len(sj)),
        "east": np.round(east[gi, gj], 1).ravel().tolist(),
        "north": np.round(north[gi, gj], 1).ravel().tolist(),
        "h": np.round(h_s, 2).ravel().tolist(),
        "gamma": np.round(g_s, 3).ravel().tolist(),
        "good": good_s.astype(np.int8).ravel().tolist(),
        "amp": np.round(np.clip(amp_n, 0, 1), 3).ravel().tolist(),
        "coh": np.round(res["coherence"][gi, gj], 3).ravel().tolist(),
        "mm": np.round(mm_n[gi, gj], 3).ravel().tolist(),
        # F31: concentrazione spettrale e riga dominante del micro-moto, sulla
        # stessa griglia (eventualmente decimata) della nuvola misurata: serve
        # solo a colorare quella nuvola con l'attributo "micro-moto", non e'
        # la rappresentazione della superficie .xml (quella e' in xml_ref).
        "mmc": np.round(res["mm_coh"][gi, gj], 3).ravel().tolist(),
        "mmf": np.round(res["mm_freq"][gi, gj], 1).ravel().tolist(),
        "simmask": res["pyr_mask"][gi, gj].astype(np.int8).ravel().tolist(),
    }

    # F33 (2026-08-31): "xml_ref" e' valutata PIXEL PER PIXEL, cioe' sull'intera
    # griglia multilooked (n_l x n_p, la stessa di height_ref e di ogni altro
    # attributo di questa pipeline), MAI decimata da cfg.surface_max_nodes.
    # Prima "href" veniva letta dalla stessa griglia rada della nuvola misurata
    # (F30): con questo dataset i due reticoli coincidevano per caso (90x93 <
    # 12000 nodi), ma il codice non lo garantiva. Ora la superficie .xml e i
    # punti di micro-moto che vi poggiano (F32) hanno il proprio reticolo,
    # indipendente e completo.
    #
    # F35: ATTENZIONE -- "pixel per pixel" qui vuol dire che la spline bicubica
    # e' CAMPIONATA una volta per ogni pixel del ritaglio, non che ogni pixel
    # porti un'informazione di quota indipendente. I nodi VERI della geolocation
    # grid dello .xml sono solo 231 su tutta la scena (raw_gcp_nodes, F35): il
    # ritaglio delle piramidi sta dentro una piccola frazione di UNA sola cella
    # di quel reticolo, quindi "xml_ref.h" e' localmente quasi un piano (F34,
    # ~2 m di escursione) per costruzione, non per errore. E' un'interpolazione
    # liscia di un riferimento rado, non una misura del profilo del suolo pixel
    # per pixel. Il payload porta anche "gcp_raw" (i nodi grezzi, non smussati,
    # con quota propria) proprio per rendere ispezionabile questo limite.
    xml_ref = {
        "n_l": int(n_l), "n_p": int(n_p),
        "east": np.round(east, 1).ravel().tolist(),
        "north": np.round(north, 1).ravel().tolist(),
        "h": np.round(res["height_ref"], 2).ravel().tolist(),
        "mask": res["pyr_mask"].astype(np.int8).ravel().tolist(),
        "mmc": np.round(res["mm_coh"], 3).ravel().tolist(),
        "mmf": np.round(res["mm_freq"], 1).ravel().tolist(),
        "mm": np.round(mm_n, 3).ravel().tolist(),
    }

    # --- nuvola tomografica ------------------------------------------------
    thr = float(np.percentile(db, 99.0))
    ii, jj, kk = np.where(db >= thr)
    if len(ii) > 20000:
        rng = np.random.default_rng(cfg.seed)
        keep = rng.choice(len(ii), 20000, replace=False)
        ii, jj, kk = ii[keep], jj[keep], kk[keep]

    pyr_mask = res["pyr_mask"]
    thr_default = float(np.median(db[ii, jj, kk])) if len(ii) else -8.0
    # F33: l'indice di solidita' e' il prodotto di tre termini normalizzati
    # 0..1 (F32), quindi si addensa vicino a zero; qui viene ristirato per
    # percentili SOLO per il colore dei voxel esportati (stesso trattamento
    # gia' usato per "amp" alla riga sopra), cosi' la rampa pieno/vuoto usa
    # davvero tutta la scala. Il valore scientifico resta res["solidity"],
    # non ristirato, usato dai profili verticali e dal JSON di validazione.
    sol_n = _norm01(sol[ii, jj, kk]) if len(ii) else np.zeros(0, dtype=np.float32)
    voxels = [[
        round(float(east[i, j]), 1), round(float(north[i, j]), 1),
        round(float(res["height_ref"][i, j] + z[k]), 1),
        round(float(db[i, j, k]), 2), round(float(sol_n[m]), 3),
        int(bool(pyr_mask[i, j])),
    ] for m, (i, j, k) in enumerate(zip(ii, jj, kk))]

    # F41: la PSF verticale della pila, cioe' la ragione per cui quella nuvola
    # ha punti anche sopra il suolo. Va calcolata QUI, dopo la selezione dei
    # voxel, perche' porta anche la frazione di punti DISEGNATI che sta sopra
    # il datum: e' il numero che il lettore vede davvero sullo schermo.
    lobi = vertical_lobe_profile(tomo, z, res["pyr_mask"],
                                 budget.delta_z_vertical)
    z_rel = z[kk] if len(kk) else np.zeros(0, dtype=np.float32)
    lobi["n_voxel"] = int(len(kk))
    lobi["frac_sopra"] = (round(float(np.mean(z_rel > 0.0)), 4)
                          if len(z_rel) else 0.0)
    apice = max(p.height_m for p in PYRAMIDS)
    lobi["quota_apice"] = round(float(apice), 1)
    lobi["frac_sopra_apice"] = (round(float(np.mean(z_rel > apice)), 4)
                                if len(z_rel) else 0.0)

    meshes = [pyramid_mesh(p, res["lat0"], res["lon0"]) for p in PYRAMIDS]

    hv = np.asarray(surface["h"], dtype=np.float64)
    apex = max(m["apex"][2] for m in meshes)
    base = min(m["base_alt"] for m in meshes)
    z_lo = float(min(np.percentile(hv, 0.5), base - 30.0))
    z_hi = float(max(np.percentile(hv, 99.5), apex + 30.0))
    # F39: il layer del suolo ora arriva fino agli apici (prima si fermava al
    # plateau), e il cursore delle quote taglia i punti fuori da [z_lo, z_hi]:
    # se non lo si tiene dentro l'intervallo, parte del suolo sparirebbe senza
    # che nulla lo segnali.
    _suolo = res.get("suolo_dem")
    if _suolo is not None:
        z_lo = min(z_lo, float(np.min(_suolo)) - 10.0)
        z_hi = max(z_hi, float(np.max(_suolo)) + 10.0)

    lvl1 = valid["livello_1_geometria"]
    val_rows = []
    for p in PYRAMIDS:
        d = lvl1.get(p.name, {})
        if "rilievo_misurato_m" in d:
            val_rows.append({
                "nome": p.name,
                "misurato": d.get("rilievo_misurato_m"),
                "atteso": d.get("rilievo_simulato_m"),
                "errore": d.get("errore_m"),
                "ok": bool(d.get("entro_2_sigma_h")),
            })
    reg = lvl1.get("regressione_misurato_vs_simulato", {})

    href_v = np.asarray(xml_ref["h"], dtype=np.float64)
    band = res["mm_band"]

    # F35: nodi grezzi della geolocation grid, con la distanza dal centro del
    # ritaglio gia' calcolata lato server -- vedi raw_gcp_nodes() e il
    # commento su xml_ref qui sopra.
    gcp_raw = dict(res["gcp_raw"])
    gcp_e = np.asarray(gcp_raw["east"])
    gcp_n = np.asarray(gcp_raw["north"])
    gcp_h_all = np.asarray(gcp_raw["h"], dtype=np.float64)
    cx_aoi = float(np.mean(east))
    cy_aoi = float(np.mean(north))
    gcp_raw["dist_km"] = np.round(
        np.hypot(gcp_e - cx_aoi, gcp_n - cy_aoi) / 1000.0, 1).tolist()

    # F38/F39: suolo (DEM) -- terreno sul datum delle quote lette dai file di
    # stack_slc, piu' il profilo delle piramidi. Vedi ground_dem_suolo().
    # None solo se cfg.suolo_dem=False.
    suolo_arr = res.get("suolo_dem")
    suolo_dem_payload = None
    if suolo_arr is not None:
        # est/nord non si ripetono: il suolo sta sullo STESSO reticolo di
        # xml_ref (l'intera griglia multilooked, mai decimata), quindi
        # serializzarne una seconda copia significava scrivere due volte le
        # stesse n_l*n_p coordinate nella pagina. La pagina le riusa da
        # xml_ref; n_l/n_p restano qui e il lettore JS verifica che coincidano.
        suolo_dem_payload = {
            "n_l": int(n_l), "n_p": int(n_p),
            "h": np.round(suolo_arr, 1).ravel().tolist(),
            "info": res.get("suolo_info") or {},
            "plateau": (res.get("suolo_plateau") or {}).get("per_piramide", []),
            "n_date_lette": (res.get("suolo_plateau") or {}).get("n_date", 0),
        }

    payload = {
        "surface": surface,
        "xml_ref": xml_ref,
        "suolo_dem": suolo_dem_payload,
        "gcp_raw": gcp_raw,
        "gcp_raw_near_km": res["gcp_raw_near_km"],
        "gcp_raw_near_h": res["gcp_raw_near_h"],
        "gcp_raw_spacing_px": list(res["gcp_raw_spacing_px"]),
        "gcp_raw_h_range": [round(float(gcp_h_all.min()), 1),
                            round(float(gcp_h_all.max()), 1)],
        "gcp_raw_h_std": round(float(gcp_h_all.std()), 1),
        "voxels": voxels,
        "lobi": lobi,
        "pyramids": meshes,
        "structures": known_structures(res["lat0"], res["lon0"], verbose=False),
        "mm": {
            "banda_hz": [round(float(band[0]), 1), round(float(band[1]), 1)],
            "rumore": round(float(res["mm_noise_floor"]), 3),
            "conc_mediana": round(float(np.median(res["mm_coh"])), 3),
            "righe": sorted({round(float(v), 1) for v in np.unique(res["mm_freq"])}),
        },
        "href_range": [round(float(href_v.min()), 1), round(float(href_v.max()), 1)],
        "thr_default": round(thr_default, 1),
        "pol": cfg.polarisation.upper(),
        "area_m": [round(float(res["area_m"][0]), 0), round(float(res["area_m"][1]), 0)],
        "n_cells": int(res["height_abs"].size),
        # scala z robusta: i percentili, non gli estremi, altrimenti pochi
        # outlier residui rendono inutilizzabili i cursori di quota
        "z_min": z_lo, "z_max": z_hi,
        "n_spike": int(res["n_spike"]),
        "dates": res["dates"],
        "master": res["master_date"],
        "baselines": [round(b.b_perp, 2) for b in res["baselines"]],
        "budget": {
            "delta_z_vertical": round(budget.delta_z_vertical, 1),
            "sigma_h": round(budget.sigma_h, 1),
            "b_spread": round(budget.b_spread, 1),
            "b_std": round(budget.b_std, 1),
            "n_baselines": budget.n_baselines,
            "target_height": budget.target_height,
            "cells": round(budget.cells_over_target, 2),
            "resolves": bool(budget.resolves_interior),
            "measurable": bool(budget.surface_measurable),
            "ambiguity": round(budget.ambiguity_height, 0),
        },
        "validation": val_rows,
        "regression": reg,
        "gamma_thr": round(float(res["gamma_thr"]), 3),
        "null_median": res["null_stats"]["mediana_nulla"],
        "good_frac": round(float(res["good"].mean()), 4),
        "plane_h": lvl1.get("_quota_mediana_piana_m"),
        "sign": int(res["sign"]),
        "layover": res["layover"],
        "n_fixes": len(FIXES),
    }

    html = _HTML_TEMPLATE.replace("__PAYLOAD__", json.dumps(
        payload, separators=(",", ":"), allow_nan=False, default=float))
    path = os.path.join(cfg.out_dir, cfg.html_name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"      superficie: {surface['n_l']}x{surface['n_p']} nodi; "
          f"voxel: {len(voxels)}")
    return path


_HTML_TEMPLATE = r"""<!doctype html>
<html lang="it">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tomografia 3D — Piana di Giza</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fira+Code:wght@400;500;600&family=Fira+Sans:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<style>
:root{
  --bg:#020617; --panel:#0E1223; --panel-2:#131a30; --line:#334155;
  --fg:#F8FAFC; --muted:#94A3B8; --accent:#16A34A; --accent-2:#38BDF8;
  --warn:#F59E0B; --bad:#DC2626; --ring:#F8FAFC;
  --sans:"Fira Sans",ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;
  --mono:"Fira Code",ui-monospace,SFMono-Regular,Consolas,monospace;
}
*{box-sizing:border-box}
html,body{height:100%}
body{margin:0;background:var(--bg);color:var(--fg);font-family:var(--sans);
     overflow:hidden;font-size:14px;line-height:1.5}
#wrap{display:flex;height:100vh}
#view{flex:1;position:relative;min-width:0}
canvas{display:block;width:100%;height:100%;cursor:grab;outline:none}
canvas:focus-visible{outline:2px solid var(--ring);outline-offset:-2px}
canvas.drag{cursor:grabbing}
/* F41: il grafico dei lobi e' un canvas di pannello, non la vista 3D:
   va sottratto alla regola width/height 100% qui sopra. */
#psf{width:100%;height:auto;cursor:default;border-radius:6px;
     background:#0A0F1F;border:1px solid var(--line);margin:2px 0 8px}
#side{width:340px;flex:0 0 340px;background:var(--panel);
      border-left:1px solid var(--line);padding:18px 18px 40px;
      overflow-y:auto;font-size:12.5px}
#side::-webkit-scrollbar{width:8px}
#side::-webkit-scrollbar-thumb{background:var(--line);border-radius:4px}
h1{font-size:16px;margin:0 0 2px;font-weight:600;letter-spacing:.01em}
.sub{color:var(--muted);font-size:11.5px;margin-bottom:14px}
.grp{margin-bottom:14px;padding-bottom:12px;border-bottom:1px solid var(--line)}
.grp:last-child{border-bottom:none}
.gt{color:var(--muted);font-size:10.5px;text-transform:uppercase;
    letter-spacing:.08em;font-weight:600;margin:0 0 8px}
label{display:block;color:var(--muted);font-size:11px;margin:9px 0 4px}
label:first-of-type{margin-top:0}
.val{float:right;color:var(--accent-2);font-family:var(--mono);font-size:11px}
input[type=range]{width:100%;accent-color:var(--accent);height:18px;cursor:pointer}
input[type=range]:focus-visible{outline:2px solid var(--ring);outline-offset:2px}
button{background:var(--panel-2);color:var(--fg);border:1px solid var(--line);
       border-radius:6px;padding:6px 10px;font-size:11.5px;font-family:var(--sans);
       cursor:pointer;margin:0 4px 4px 0;transition:background .18s,border-color .18s}
button:hover{border-color:var(--accent);background:#1a2338}
button:focus-visible{outline:2px solid var(--ring);outline-offset:2px}
button.on{background:var(--accent);color:#04140A;border-color:var(--accent);font-weight:600}
table{width:100%;border-collapse:collapse;font-size:11.5px}
td{padding:3px 0;color:var(--muted);vertical-align:top}
td:last-child{text-align:right;color:var(--fg);font-family:var(--mono);font-size:11px}
.note{background:#1a1206;border:1px solid #5a3f0f;border-radius:6px;padding:9px 10px;
      font-size:11px;line-height:1.55;color:#FDE9C8}
.ok{background:#04180d;border:1px solid #14532d;color:#BBF7D0}
.legend{display:flex;height:8px;border-radius:4px;overflow:hidden;margin:7px 0 3px}
.legend i{flex:1}
.lbl{display:flex;justify-content:space-between;color:var(--muted);font-size:10px;
     font-family:var(--mono)}
#hud{position:absolute;left:16px;top:14px;font-size:11px;color:var(--muted);
     pointer-events:none;line-height:1.65;font-family:var(--mono);
     text-shadow:0 0 10px rgba(0,0,0,.9)}
#hud b{color:var(--fg)}
#keys{position:absolute;left:16px;bottom:14px;font-size:10.5px;color:var(--muted);
      font-family:var(--mono);pointer-events:none;opacity:.8}
.pill{display:inline-block;padding:1px 6px;border-radius:999px;font-size:10px;
      font-family:var(--mono);border:1px solid var(--line)}
.pill.g{color:#BBF7D0;border-color:#14532d;background:#04180d}
.pill.r{color:#FECACA;border-color:#7f1d1d;background:#1b0606}
@media (prefers-reduced-motion: reduce){ *{transition:none !important} }
@media (max-width:900px){
  #wrap{flex-direction:column}
  #side{width:100%;flex:0 0 auto;max-height:46vh;border-left:none;
        border-top:1px solid var(--line)}
}
</style>
</head>
<body>
<div id="wrap">
  <div id="view">
    <canvas id="cv" tabindex="0" aria-label="Vista tridimensionale della superficie ricostruita"></canvas>
    <div id="hud"></div>
    <div id="keys">frecce ruota · shift+frecce trasla · +/− zoom · Q/E rollio · R reset · spazio rotazione</div>
  </div>
  <aside id="side">
    <h1>Tomografia 3D — Piana di Giza</h1>
    <div class="sub">Sentinel-1 · superficie ricostruita dai .tiff pixel per pixel,
      geometria dagli .xml. <b>Un punto per nodo misurato</b>, non una maglia:
      dove non c'è misura resta vuoto, e i nodi sotto soglia sono punti piccoli
      e spenti. La quota di ciascun punto e' <b>quota .xml (bilineare locale,
      F36) + quota relativa misurata dal periodogramma</b>: la geometria
      dello .xml resta il datum del calcolo tomografico (fase, k<sub>z</sub>,
      geocodifica), ma non e' piu' disegnata come strato a se' (rimosso, F37:
      sul ritaglio ha solo ~1 m di escursione, vedi pannello «Nodi grezzi»
      sotto per i valori .xml reali). Il <b>suolo (DEM)</b> in rosso e' il
      terreno, <b>piramidi comprese</b> (F39): il datum sono le quote lette
      dagli annotation.xml dello stack, il rilievo locale viene dal DEM
      esterno (Copernicus via Open-Meteo) riportato su quel datum, e sopra
      c'e' il profilo geometrico delle tre piramidi appoggiato al terreno —
      prima il DEM esterno da solo dava 66 m sull'apice di Cheope, cioe' il
      plateau nudo, e il «suolo» passava dentro le piramidi. Resta un
      RIFERIMENTO di contesto, mai usato nel
      calcolo. Sotto, a filo di ferro, gli altri <b>riferimenti</b>
      che misure non sono: le piramidi ideali e le camere note, solidali con
      la piramide sotto qualunque rotazione della vista. I <b>voxel</b> sono
      colorati per indice di pieno/vuoto, che incorpora il micro-moto.
      <b>Area ristretta alle sole piramidi</b> in geometria radar;
      <i>Solo piramidi</i> mostra anche il margine di calibrazione.</div>

    <div class="grp">
      <div class="gt">Vista</div>
      <button id="play" class="on">Rotazione</button>
      <button id="reset">Reset</button>
      <button id="proj">Prospettiva</button>
      <label>Imbardata <span class="val" id="yawV"></span></label>
      <input type="range" id="yaw" min="0" max="360" step="1" value="34">
      <label>Beccheggio <span class="val" id="pitV"></span></label>
      <input type="range" id="pit" min="-85" max="85" step="1" value="26">
      <label>Rollio <span class="val" id="rolV"></span></label>
      <input type="range" id="rol" min="-45" max="45" step="1" value="0">
      <label>Zoom <span class="val" id="zoomV"></span></label>
      <input type="range" id="zoom" min="0.3" max="6" step="0.05" value="1">
      <label>Esagerazione verticale <span class="val" id="exV"></span></label>
      <input type="range" id="ex" min="1" max="14" step="0.5" value="4">
    </div>

    <div class="grp">
      <div class="gt">Livelli</div>
      <button id="lSurf" class="on">Superficie misurata</button>
      <button id="lDem" class="on">Suolo (DEM)</button>
      <button id="lMM">Punti micro-moto</button>
      <button id="lStr" class="on">Strutture note</button>
      <button id="lVox">Voxel (pieno/vuoto)</button>
      <button id="lPyr" class="on">Piramidi ideali</button>
      <button id="lGrid" class="on">Griglia</button>
      <button id="lAxes" class="on">Assi</button>
      <button id="lOnlyPyr" class="on">Solo piramidi</button>
      <button id="lMask">Solo celle di qualita'</button>
    </div>

    <div class="grp">
      <div class="gt">Superficie</div>
      <label>Colore <span class="val" id="attrV">quota</span></label>
      <button class="attr on" data-a="h">Quota</button>
      <button class="attr" data-a="amp">Ampiezza</button>
      <button class="attr" data-a="gamma">Coerenza fit</button>
      <button class="attr" data-a="coh">Coerenza interf.</button>
      <button class="attr" data-a="mm">Micro-moto</button>
      <button class="attr" data-a="mmc">Conc. micro-moto</button>
      <div class="legend" id="lg"></div>
      <div class="lbl"><span id="lgA"></span><span id="lgB"></span></div>
      <label>Dimensione punto <span class="val" id="dotV"></span></label>
      <input type="range" id="dot" min="0.8" max="6" step="0.2" value="2.2"
             aria-label="Raggio in pixel dei punti della superficie">
      <label>Azimut del sole <span class="val" id="sazV"></span></label>
      <input type="range" id="saz" min="0" max="360" step="5" value="315">
      <label>Elevazione del sole <span class="val" id="selV"></span></label>
      <input type="range" id="sel" min="5" max="85" step="1" value="42">
      <label>Opacita' <span class="val" id="opV"></span></label>
      <input type="range" id="op" min="0.15" max="1" step="0.05" value="1">
    </div>

    <div class="grp">
      <div class="gt">Sezione e soglie</div>
      <label>Quota minima <span class="val" id="zloV"></span></label>
      <input type="range" id="zlo" min="0" max="100" step="1" value="0">
      <label>Quota massima <span class="val" id="zhiV"></span></label>
      <input type="range" id="zhi" min="0" max="100" step="1" value="100">
      <label>Piano di sezione <span class="val" id="secV"></span></label>
      <input type="range" id="sec" min="0" max="100" step="1" value="100">
      <button id="secAx">Sezione: Est</button>
      <label>Soglia voxel <span class="val" id="thrV"></span></label>
      <input type="range" id="thr" min="-40" max="0" step="0.5" value="-8"
             aria-label="Soglia di intensita dei voxel">
      <label>Dimensione voxel <span class="val" id="szV"></span></label>
      <input type="range" id="sz" min="1" max="6" step="0.5" value="2">
      <label>Colore voxel: pieno/vuoto (con micro-moto)</label>
      <div class="legend" id="lgVox"></div>
      <div class="lbl"><span>vuoto</span><span>pieno</span></div>
    </div>

    <div class="grp">
      <div class="gt">Lobi verticali — perché i voxel stanno anche in aria</div>
      <canvas id="psf" width="304" height="250"
              aria-label="Profilo verticale dell'energia tomografica in dB"></canvas>
      <div class="lbl" style="margin-bottom:6px">
        <span style="color:var(--accent-2)">&#9473; sole piramidi</span>
        <span style="color:#64748B">&#9473; tutto il ritaglio</span>
      </div>
      <table id="pt"></table>
      <div style="margin-top:6px;color:var(--muted);font-size:11px;line-height:1.5"
           id="psfnote"></div>
    </div>

    <div class="grp">
      <div class="gt">Micro-moto</div>
      <label>Soglia di concentrazione spettrale <span class="val" id="mmtV"></span></label>
      <input type="range" id="mmt" min="0" max="1" step="0.01" value="0.35"
             aria-label="Soglia sulla concentrazione spettrale della traccia">
      <table id="mt"></table>
    </div>

    <div class="grp">
      <div class="gt">Strutture interne note</div>
      <table id="st"></table>
    </div>

    <div class="grp">
      <div class="gt">Budget</div>
      <table id="bt"></table>
    </div>

    <div class="grp">
      <div class="gt">Validazione livello 1 — geometria</div>
      <table id="rt"></table>
      <table id="vt" style="margin-top:8px"></table>
    </div>

    <div class="grp">
      <div class="gt">Suolo (DEM): quote lette dallo stack — F39</div>
      <table id="suot"></table>
      <div style="margin-top:6px;color:var(--muted);font-size:11px;line-height:1.5">
        Il layer del suolo e' costruito in tre pezzi dichiarati: il DATUM sono
        le quote <code>&lt;geolocationGridPoint&gt;&lt;height&gt;</code> lette
        dagli <code>annotation.xml</code> dello stack (bilineare locale del
        master, F36, verificata su tutte le date); il RILIEVO viene dal DEM
        esterno riportato su quel datum; le PIRAMIDI sono il profilo
        geometrico appoggiato al terreno. Il DEM esterno da solo dava 66 m
        sull'apice di Cheope — il plateau nudo — e il suolo passava dentro le
        piramidi. Lo scarto fra la base ricavata dal terreno e quella di
        letteratura e' mostrato e non corretto: il datum dello .xml varia
        ~1 m sul ritaglio (F35) e non puo' risolvere le basi delle singole
        piramidi.
      </div>
    </div>

    <div class="grp">
      <div class="gt">Nodi grezzi della geolocation grid (.xml) — F35</div>
      <table id="gcpt"></table>
      <div style="margin-top:6px;color:var(--muted);font-size:11px;line-height:1.5">
        Il datum usato per ogni quota di questa pagina (quota .xml sopra,
        F36) e' la bilineare fra i quattro nodi VERI che racchiudono il
        ritaglio; questi sono TUTTI i nodi VERI della scena, senza
        interpolazione, con la loro quota propria. Sono lontani dal ritaglio
        perche' il reticolo dello .xml copre l'intera scena con appena 231
        nodi: la bilineare locale non ha rilievo da mostrare perche' nessun
        nodo reale cade li' vicino, non perche' il calcolo lo appiattisca.
      </div>
    </div>

    <div class="grp">
      <div id="note" class="note"></div>
    </div>
  </aside>
</div>
<script>
"use strict";
const D = __PAYLOAD__;
const S = D.surface;
const NL = S.n_l, NP = S.n_p, NN = NL * NP, NQ = (NL - 1) * (NP - 1);
const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
const idx = (i, j) => i * NP + j;
/* F33: la superficie di riferimento .xml (e i punti di micro-moto che vi
   poggiano) vive su un reticolo PROPRIO, l'intera griglia multilooked, mai
   decimato da cfg.surface_max_nodes come la nuvola misurata sopra. */
const XR = D.xml_ref;
const XNL = XR.n_l, XNP = XR.n_p, XNN = XNL * XNP;
const idxX = (i, j) => i * XNP + j;
/* F38/F39: "suolo (DEM)" -- il TERRENO, piramidi comprese. Datum = le quote
   lette dagli annotation.xml di stack_slc; rilievo = DEM esterno (Copernicus
   via Open-Meteo) riportato su quel datum; sopra, il profilo geometrico delle
   piramidi appoggiato al terreno (F39: il solo DEM esterno dava 66 m
   sull'apice di Cheope, cioe' il plateau nudo, e il suolo passava dentro le
   piramidi). Null solo se cfg.suolo_dem=False: in quel caso il layer e il suo
   pulsante restano semplicemente assenti. */
const DEM = D.suolo_dem;
const DNL = DEM ? DEM.n_l : 0, DNP = DEM ? DEM.n_p : 0, DNN = DNL * DNP;

/* =====================================================================
   Dati in array tipizzati: gli array JS generici costringono il motore a
   controllare il tipo a ogni accesso, e qui gli accessi sono centinaia di
   migliaia per fotogramma.
   ===================================================================== */
const EAST = Float32Array.from(S.east), NORTH = Float32Array.from(S.north);
const HGT  = Float32Array.from(S.h);
const GOOD = Uint8Array.from(S.good), SIMM = Uint8Array.from(S.simmask);
/* F32/F33: XHREF e' la quota dei pixel letta dalla geolocation grid degli
   .xml, pixel per pixel sul reticolo xml_ref; XMMC/XMMF/XMME sono la
   concentrazione spettrale, la riga dominante in Hz e l'energia del
   micro-moto sullo STESSO reticolo (F31). XHREF regge due livelli: la
   superficie di riferimento e la quota su cui poggiano i punti del
   micro-moto, che e' un fenomeno di superficie e non ha una profondita'
   propria. MMC/MMF/MME (sotto) restano sul reticolo della nuvola misurata,
   e servono solo a colorare quella nuvola con l'attributo "micro-moto". */
const XEAST = Float32Array.from(XR.east), XNORTH = Float32Array.from(XR.north);
const XHREF = Float32Array.from(XR.h),    XMASK  = Uint8Array.from(XR.mask);
const XMMC  = Float32Array.from(XR.mmc),  XMMF   = Float32Array.from(XR.mmf);
const XMME  = Float32Array.from(XR.mm);
const MMC  = Float32Array.from(S.mmc), MMF = Float32Array.from(S.mmf);
const MME  = Float32Array.from(S.mm);
const ATTR = {h:HGT, amp:Float32Array.from(S.amp), gamma:Float32Array.from(S.gamma),
              coh:Float32Array.from(S.coh), mm:MME, mmc:MMC};
const PX = new Float32Array(NN), PY = new Float32Array(NN), PZ = new Float32Array(NN);
const RX = new Float32Array(XNN), RY = new Float32Array(XNN), RZ = new Float32Array(XNN);
/* il suolo condivide il reticolo di xml_ref: le sue est/nord non
   viaggiano piu' nel payload, si riusano quelle. Se un giorno i due
   reticoli divergessero, DEM.east tornerebbe presente e vincerebbe. */
const DEAST = DEM ? (DEM.east ? Float32Array.from(DEM.east) : XEAST)
                  : new Float32Array(0);
const DNORTH = DEM ? (DEM.north ? Float32Array.from(DEM.north) : XNORTH)
                   : new Float32Array(0);
const DHGT  = DEM ? Float32Array.from(DEM.h) : new Float32Array(0);
const DX = new Float32Array(DNN), DY = new Float32Array(DNN), DZ = new Float32Array(DNN);

const cv = document.getElementById('cv'), ctx = cv.getContext('2d');
let W = 0, H = 0, dpr = Math.min(devicePixelRatio || 1, 2);
function resize(){
  const r = cv.parentElement.getBoundingClientRect();
  W = Math.max(1, r.width); H = Math.max(1, r.height);
  cv.width = W * dpr; cv.height = H * dpr; ctx.setTransform(dpr,0,0,dpr,0,0);
  spriteDirty = true;          // l'atlante e' alla risoluzione del dispositivo
  invalidate();
}
addEventListener('resize', resize);

/* ---------- stato: i gradi di liberta' della vista ---------- */
const V = {
  yaw:34*Math.PI/180, pitch:26*Math.PI/180, roll:0,
  zoom:1, panX:0, panY:0, zEx:4,
  persp:true, autoRot:!reduce,
  attr:'h', opacity:1, dot:2.2, sunAz:315, sunEl:42,
  zlo:0, zhi:1, sec:1, secAxis:'E',
  thr:(D.thr_default!=null?D.thr_default:-8), vsize:2,
  surf:true, vox:false, pyr:true, grid:true, axes:true,
  onlyGood:false, onlyPyr:true,
  /* F32: i tre livelli nuovi. mmThr e' una soglia sulla CONCENTRAZIONE
     spettrale, non sull'energia: filtra i pixel la cui traccia e' rumore
     bianco invece che una vibrazione a riga singola. */
  mmp:false, str:true, mmThr:0.35,
  dem:!!DEM   // F38: attivo solo se il DEM esterno e' stato scaricato
};

/* ---------- F29: disegno A RICHIESTA -----------------------------------
   Prima si ridisegnava a ogni fotogramma anche a scena ferma, e con
   prefers-reduced-motion il ciclo non partiva affatto: nessun cursore aveva
   effetto. Ora un flag dice quando la scena e' cambiata; a riposo il ciclo
   non fa nulla e il browser resta libero di rispondere. Il flag delle
   normali e' separato perche' l'illuminazione dipende solo da sole ed
   esagerazione verticale, non dal punto di vista.                          */
let dirty = true, shadeDirty = true, colorDirty = true, dragging = false;
function invalidate(){ dirty = true; }
function invalidateShade(){ shadeDirty = true; dirty = true; }
function invalidateColor(){ colorDirty = true; spriteDirty = true; dirty = true; }

/* ---------- centro e scala della scena ---------- */
let cxs=0, cys=0, zOff=0, span=1600;
let RANGE_E=[0,1], RANGE_N=[0,1];
{
  let se=0, sn=0, k=0, e0=1e18, e1=-1e18, n0=1e18, n1=-1e18;
  for(let i=0;i<NN;i++){
    if(EAST[i]<e0)e0=EAST[i]; if(EAST[i]>e1)e1=EAST[i];
    if(NORTH[i]<n0)n0=NORTH[i]; if(NORTH[i]>n1)n1=NORTH[i];
  }
  // il centro e' il baricentro delle sole celle delle piramidi: la vista
  // inquadra i monumenti, non il margine di calibrazione che sta attorno
  for(let i=0;i<NN;i++){ if(SIMM[i]){ se+=EAST[i]; sn+=NORTH[i]; k++; } }
  if(k<10){ se=0; sn=0; k=0; for(let i=0;i<NN;i++){se+=EAST[i];sn+=NORTH[i];k++;} }
  cxs=se/k; cys=sn/k;
  zOff = D.plane_h!=null ? D.plane_h : (D.z_min+D.z_max)/2;
  RANGE_E=[e0,e1]; RANGE_N=[n0,n1];
  span = Math.max(400, 1.18*Math.max(e1-e0, n1-n0));
}

/* ---------- proiezione ---------- */
let _cy=1,_sy=0,_cp=1,_sp=0,_cr=1,_sr=0,_scale=1,_pd=1;
function updateCam(){
  _cy=Math.cos(V.yaw); _sy=Math.sin(V.yaw);
  _cp=Math.cos(V.pitch); _sp=Math.sin(V.pitch);
  _cr=Math.cos(V.roll); _sr=Math.sin(V.roll);
  _scale=Math.min(W,H)/span*V.zoom;
  _pd=span*1.9;
}
/* i tre limiti di taglio dipendono solo dai cursori, quindi sono costanti
   dentro un fotogramma: erano ricalcolati in drawSurface, drawMM, drawVoxels,
   drawDem e hud, cioe' cinque volte le stesse due formule a ogni frame.
   Stessa scelta gia' fatta per la camera qui sopra. */
let _zlo=0,_zhi=0,_secLimit=0;
function updateCuts(){
  _zlo=D.z_min+(D.z_max-D.z_min)*V.zlo;
  _zhi=D.z_min+(D.z_max-D.z_min)*V.zhi;
  _secLimit = V.secAxis==='E'
      ? RANGE_E[0]+(RANGE_E[1]-RANGE_E[0])*V.sec
      : RANGE_N[0]+(RANGE_N[1]-RANGE_N[0])*V.sec;
}
function project(x,y,z){
  x-=cxs; y-=cys; z=(z-zOff)*V.zEx;
  const X=x*_cy-y*_sy, Y=x*_sy+y*_cy;
  const up=z*_cp-Y*_sp, depth=Y*_cp+z*_sp;
  const X2=X*_cr-up*_sr, U2=X*_sr+up*_cr;
  let s=_scale;
  if(V.persp) s*=_pd/Math.max(_pd+depth,_pd*0.12);
  return [W/2+X2*s+V.panX, H/2-U2*s+V.panY, depth];
}
/* proiezione di TUTTI i nodi in un colpo solo: con la maglia ogni poligono
   riproiettava i suoi quattro vertici, cioe' quattro volte lo stesso lavoro
   per ogni vertice interno; con la nuvola serve comunque una passata sola.
   Le cinque righe di algebra sono le stesse di project() ed e' una ripetizione
   VOLUTA, non una svista: chiamare project() qui dentro allocherebbe un array
   di tre elementi per ciascuno degli XNN/NN/DNN nodi a ogni fotogramma.
   F33: ex/nx/n parametrizzano il reticolo, cosi' la stessa funzione serve
   sia alla nuvola misurata (EAST/NORTH, NN) sia al reticolo xml_ref
   (XEAST/XNORTH, XNN), che ha una propria dimensione. */
function projectNodes(h, ex, nx, n, ox, oy, oz){
  for(let i=0;i<n;i++){
    const x=ex[i]-cxs, y=nx[i]-cys, z=(h[i]-zOff)*V.zEx;
    const X=x*_cy-y*_sy, Y=x*_sy+y*_cy;
    const up=z*_cp-Y*_sp, depth=Y*_cp+z*_sp;
    const X2=X*_cr-up*_sr, U2=X*_sr+up*_cr;
    let s=_scale;
    if(V.persp) s*=_pd/Math.max(_pd+depth,_pd*0.12);
    ox[i]=W/2+X2*s+V.panX; oy[i]=H/2-U2*s+V.panY; oz[i]=depth;
  }
}

/* ---------- rampe di colore e tabella precalcolata ---------- */
function mk(stops){ return t=>{
  t=t<0?0:(t>1?1:t);
  const f=t*(stops.length-1), i=Math.min(f|0,stops.length-2), g=f-i;
  const a=stops[i], b=stops[i+1];
  return [a[0]+(b[0]-a[0])*g|0, a[1]+(b[1]-a[1])*g|0, a[2]+(b[2]-a[2])*g|0];
};}
const RAMP = {
  h:    mk([[12,32,64],[24,88,110],[62,140,108],[168,178,96],[224,196,142],[248,246,242]]),
  amp:  mk([[8,6,28],[70,15,90],[158,38,100],[219,86,68],[248,158,54],[252,232,158]]),
  gamma:mk([[10,20,45],[22,101,52],[74,222,128],[220,252,231]]),
  coh:  mk([[10,20,45],[22,101,52],[74,222,128],[220,252,231]]),
  mm:   mk([[10,20,50],[30,120,180],[240,240,235],[220,120,50],[170,30,30]]),
  mmc:  mk([[14,16,34],[60,50,120],[150,70,150],[230,140,110],[252,232,158]]),
  /* F33: rampa pieno/vuoto per i voxel. Alto (pieno) = indice di solidita'
     alto = energia forte, coerente, poco micro-moto anomalo; basso (vuoto)
     = il contrario. E' un discriminante multi-attributo (F31), non una
     rilevazione di cavita' risolta in profondita'. */
  sol:  mk([[10,15,45],[35,70,130],[100,160,150],[215,190,110],[235,110,40]])
};
const ATTR_LABEL = {h:['quota bassa','quota alta'], amp:['debole','forte'],
  gamma:['stima incerta','stima solida'], coh:['decorrelato','coerente'],
  mm:['fermo','vibrante'], mmc:['rumore bianco','riga singola']};

const NC=40, NS=16;
const RANGE = {};
for(const k in ATTR){
  const a=Array.from(ATTR[k]).filter(Number.isFinite).sort((p,q)=>p-q);
  RANGE[k]= a.length ? [a[(a.length*0.02)|0], a[(a.length*0.98)|0]] : [0,1];
}
/* F29: la stringa "rgb(...)" costava una allocazione e un parsing CSS per
   ogni elemento di ogni fotogramma. Ora ne esistono NC*NS in tutto,
   costruite una volta, piu' una copia attenuata per i nodi sotto soglia. */
let COLTAB=null, CIDX=null;
function buildColors(){
  const ramp=RAMP[V.attr], arr=ATTR[V.attr];
  const lo=RANGE[V.attr][0], hi=RANGE[V.attr][1], inv=1/Math.max(hi-lo,1e-9);
  CIDX=new Uint8Array(NN);
  for(let i=0;i<NN;i++){
    let t=(arr[i]-lo)*inv; t=t<0?0:(t>1?1:t);
    CIDX[i]=Math.min(NC-1,(t*NC)|0);
  }
  COLTAB=[new Array(NC*NS), new Array(NC*NS)];
  for(let c=0;c<NC;c++){
    const col=ramp((c+0.5)/NC);
    for(let sIdx=0;sIdx<NS;sIdx++){
      const sh=0.34+0.66*((sIdx+0.5)/NS);
      COLTAB[0][c*NS+sIdx]='rgb('+(col[0]*sh|0)+','+(col[1]*sh|0)+','+(col[2]*sh|0)+')';
      COLTAB[1][c*NS+sIdx]='rgb('+(col[0]*sh*0.45|0)+','+(col[1]*sh*0.45|0)+','+(col[2]*sh*0.45|0)+')';
    }
  }
  colorDirty=false;
}

/* ---------- ombreggiatura: dipende solo da sole e esagerazione ---------- */
const SHADE=new Uint8Array(NQ);
function buildShade(){
  const az=V.sunAz*Math.PI/180, el=V.sunEl*Math.PI/180;
  const lx=Math.sin(az)*Math.cos(el), ly=Math.cos(az)*Math.cos(el), lz=Math.sin(el);
  for(let i=0;i<NL-1;i++){
    for(let j=0;j<NP-1;j++){
      const a=idx(i,j), b=idx(i,j+1), d=idx(i+1,j);
      const ux=EAST[b]-EAST[a], uy=NORTH[b]-NORTH[a], uz=(HGT[b]-HGT[a])*V.zEx;
      const vx=EAST[d]-EAST[a], vy=NORTH[d]-NORTH[a], vz=(HGT[d]-HGT[a])*V.zEx;
      let nx=uy*vz-uz*vy, ny=uz*vx-ux*vz, nz=ux*vy-uy*vx;
      const nn=Math.hypot(nx,ny,nz)||1; nx/=nn; ny/=nn; nz/=nn;
      if(nz<0){nx=-nx;ny=-ny;nz=-nz;}
      let lam=nx*lx+ny*ly+nz*lz; if(lam<0)lam=0;
      SHADE[i*(NP-1)+j]=Math.min(NS-1,(lam*NS)|0);
    }
  }
  shadeDirty=false;
}

/* ---------- atlante di dischetti (F30) ---------------------------------
   Con la maglia i poligoni erano poche migliaia e un fillStyle serviva per
   un'area grande; con la nuvola servono 8 000 dischetti per fotogramma, e
   costruire un arco e riempirlo uno per uno costava 96 ms. Qui il dischetto
   e' disegnato UNA volta per ciascuno dei NC*NS colori della tabella (per
   due volte: raggio pieno per i nodi sopra soglia, ridotto per quelli
   sotto), e il disegno diventa una drawImage per punto. L'atlante e' alla
   risoluzione del dispositivo, altrimenti su schermi HiDPI i punti
   risulterebbero sfocati.
   L'ultima riga dell'atlante (indice 2*NS) e' quella dei punti di
   micro-moto: colore dalla rampa mm, nessuna ombreggiatura (una vibrazione
   non ha una normale) e raggio maggiorato, cosi' il livello si distingue a
   colpo d'occhio dalla nuvola della superficie.                           */
let SPRITE=null, SP_CELL=0, SP_DEV=0, MM_ROW=0, spriteDirty=true;
function invalidateSprite(){ spriteDirty=true; DEM_SPRITE=null; dirty=true; }
function buildSprites(){
  if(colorDirty) buildColors();
  const r=V.dot, rMM=r*1.5;
  const cell=Math.ceil(2*Math.max(r,rMM)+2), dev=Math.ceil(cell*dpr);
  const cvs=document.createElement('canvas');
  cvs.width=NC*dev; cvs.height=(2*NS+1)*dev;
  const c=cvs.getContext('2d');
  const k=dev/cell;
  for(let g=0; g<2; g++){
    const rr=(g===0 ? r : Math.max(0.6, r*0.55))*k;
    for(let ci=0; ci<NC; ci++){
      for(let sIdx=0; sIdx<NS; sIdx++){
        c.fillStyle=COLTAB[g][ci*NS+sIdx];
        c.beginPath();
        c.arc(ci*dev+dev/2, (g*NS+sIdx)*dev+dev/2, rr, 0, 6.2832);
        c.fill();
      }
    }
  }
  MM_ROW=2*NS;
  for(let ci=0; ci<NC; ci++){
    const col=RAMP.mm((ci+0.5)/NC);
    c.fillStyle='rgb('+col[0]+','+col[1]+','+col[2]+')';
    c.beginPath();
    c.arc(ci*dev+dev/2, MM_ROW*dev+dev/2, rMM*k, 0, 6.2832);
    c.fill();
  }
  SPRITE=cvs; SP_CELL=cell; SP_DEV=dev; spriteDirty=false;
}

/* ---------- griglia e assi ---------- */
function drawGrid(){
  const step=100, n=Math.ceil(span/2/step);
  ctx.strokeStyle='rgba(148,163,184,.13)'; ctx.lineWidth=1;
  ctx.beginPath();
  for(let i=-n;i<=n;i++){
    let a=project(cxs+i*step,cys-n*step,zOff), b=project(cxs+i*step,cys+n*step,zOff);
    ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]);
    a=project(cxs-n*step,cys+i*step,zOff); b=project(cxs+n*step,cys+i*step,zOff);
    ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]);
  }
  ctx.stroke();
}
function drawAxes(){
  const L=200, o=project(cxs,cys,zOff);
  const ax=[[project(cxs+L,cys,zOff),'E','#38BDF8'],
            [project(cxs,cys+L,zOff),'N','#16A34A'],
            [project(cxs,cys,zOff+L/2),'Z','#F59E0B']];
  ctx.lineWidth=1.6; ctx.font='600 11px "Fira Code",monospace';
  for(const row of ax){
    const p=row[0], lab=row[1], col=row[2];
    ctx.strokeStyle=col; ctx.beginPath(); ctx.moveTo(o[0],o[1]); ctx.lineTo(p[0],p[1]); ctx.stroke();
    ctx.fillStyle=col; ctx.fillText(lab,p[0]+5,p[1]+4);
  }
}

/* ---------- superficie: nuvola di punti (F30) ----------------------------
   Un punto per NODO misurato, non una faccia per quadrupla di nodi. Il
   poligono riempiva lo spazio fra quattro nodi qualsiasi, e con il 23 per
   cento di celle sopra soglia quel riempimento era la maggior parte del
   colore sullo schermo: una superficie continua disegnata sopra un campo
   che continuo non e'. Il punto mostra dove c'e' una misura e lascia vuoto
   dove non ce n'e'. I nodi sotto soglia restano, ma piccoli e spenti.      */
function drawSurface(step){
  if(colorDirty) buildColors();
  if(shadeDirty) buildShade();
  const zlo=_zlo, zhi=_zhi, secLimit=_secLimit;
  /* F20/F29: per un campo di quote su reticolo regolare l'ordine painter si
     ottiene percorrendo righe e colonne dal lato lontano a quello vicino.
     E' esatto (la profondita' e' affine in (i,j)) e costa O(n) invece di un
     ordinamento su migliaia di elementi a ogni fotogramma. */
  const iRev = PZ[idx(NL-1,0)] > PZ[idx(0,0)];
  const jRev = PZ[idx(0,NP-1)] > PZ[idx(0,0)];
  if(spriteDirty) buildSprites();
  ctx.globalAlpha=V.opacity;
  const cell=SP_CELL, dev=SP_DEV, half=cell/2;
  let drawn=0;
  for(let ai=0; ai<NL; ai+=step){
    const i = iRev ? NL-1-ai : ai;
    const si = i<NL-1 ? i : Math.max(NL-2, 0);
    for(let aj=0; aj<NP; aj+=step){
      const j = jRev ? NP-1-aj : aj;
      const a=idx(i,j);
      if(V.onlyPyr && !SIMM[a]) continue;
      const ok=GOOD[a]!==0;
      if(V.onlyGood && !ok) continue;
      const h=HGT[a];
      if(!(h>=zlo && h<=zhi)) continue;
      if(V.secAxis==='E' ? EAST[a]>secLimit : NORTH[a]>secLimit) continue;
      const sj = j<NP-1 ? j : Math.max(NP-2, 0);
      const sh = SHADE[si*(NP-1)+sj];
      ctx.drawImage(SPRITE, CIDX[a]*dev, ((ok?0:NS)+sh)*dev, dev, dev,
                    PX[a]-half, PY[a]-half, cell, cell);
      drawn++;
    }
  }
  ctx.globalAlpha=1;
  return drawn;
}

/* F37: la superficie di riferimento .xml (F32/F33/F34/F35/F36) non e' piu'
   disegnata come strato a se' (drawRef, rimosso su richiesta esplicita: sul
   ritaglio delle piramidi coincide quasi esattamente col piano orizzontale a
   quota costante, quindi come layer visivo non aggiungeva informazione oltre
   al pannello «Nodi grezzi» sotto). I dati (XHREF/XEAST/XNORTH, dal .xml
   corrispondente al master, F36) restano: sono il datum su cui poggiano i
   punti di micro-moto qui sotto, e su cui e' costruita ogni quota di questa
   pagina (surface.h = height_ref + dh). */

/* ---------- punti di micro-moto (F32/F33) --------------------------------
   Un punto per pixel del reticolo xml_ref la cui traccia di vibrazione e'
   una riga singola e non rumore bianco: il filtro e' sulla CONCENTRAZIONE
   spettrale (F31), il colore sull'energia. I punti stanno sulla superficie
   di riferimento degli .xml perche' il micro-moto e' un attributo di
   SUPERFICIE - il banco di sub-aperture non gli assegna una profondita', e
   fingere il contrario sarebbe esattamente la confusione che le fonti
   chiedono di evitare.                                                    */
function drawMM(step){
  if(spriteDirty) buildSprites();
  const zlo=_zlo, zhi=_zhi, secLimit=_secLimit;
  const iRev = RZ[idxX(XNL-1,0)] > RZ[idxX(0,0)];
  const jRev = RZ[idxX(0,XNP-1)] > RZ[idxX(0,0)];
  const cell=SP_CELL, dev=SP_DEV, half=cell/2, row=MM_ROW*SP_DEV;
  let n=0;
  for(let ai=0; ai<XNL; ai+=step){
    const i = iRev ? XNL-1-ai : ai;
    for(let aj=0; aj<XNP; aj+=step){
      const j = jRev ? XNP-1-aj : aj;
      const a=idxX(i,j);
      if(XMMC[a] < V.mmThr) continue;
      if(V.onlyPyr && !XMASK[a]) continue;
      const h=XHREF[a];
      if(!(h>=zlo && h<=zhi)) continue;
      if(V.secAxis==='E' ? XEAST[a]>secLimit : XNORTH[a]>secLimit) continue;
      let t=XMME[a]; t=t<0?0:(t>1?1:t);
      ctx.drawImage(SPRITE, Math.min(NC-1,(t*NC)|0)*dev, row, dev, dev,
                    RX[a]-half, RY[a]-half, cell, cell);
      n++;
    }
  }
  return n;
}

/* ---------- strutture interne note (F32) --------------------------------
   Camere e corridoi di Cheope e Chefren da piramide_cheope_3d.py e
   piramide_kefren_3d.py. Sono un RIFERIMENTO archeologico: filo di ferro,
   mai punti, cosi' non possono essere scambiati per qualcosa che questi
   dati hanno rilevato. Ordine painter sulla profondita' del baricentro.  */
function drawStructures(){
  const ST=D.structures||[];
  if(!ST.length) return 0;
  const ord=[];
  for(let k=0;k<ST.length;k++){
    const e=ST[k].edges; let d=0;
    for(let q=0;q<e.length;q+=3) d+=project(e[q],e[q+1],e[q+2])[2];
    ord.push([k, d/(e.length/3)]);
  }
  ord.sort((a,b)=>b[1]-a[1]);
  ctx.lineWidth=1.1;
  for(const row of ord){
    const s=ST[row[0]], e=s.edges;
    ctx.strokeStyle=s.colore;
    ctx.beginPath();
    for(let q=0;q<e.length;q+=6){
      const p=project(e[q],e[q+1],e[q+2]), r=project(e[q+3],e[q+4],e[q+5]);
      ctx.moveTo(p[0],p[1]); ctx.lineTo(r[0],r[1]);
    }
    ctx.stroke();
  }
  return ST.length;
}

/* ---------- voxel: nuvola pieno/vuoto (F33) ------------------------------
   La soglia "thr" resta sull'AMPIEZZA (v[3], dB): e' il cancello di qualita'
   che decide se il voxel e' un ritorno tomografico abbastanza forte da
   essere disegnato. Il COLORE pero' e' l'indice di solidita' (v[4],
   solidity_index di F31/F32: energia x coerenza x (1 - micro-moto),
   ristirato per percentili solo per la resa - vedi build_html), cioe' il
   discriminante pieno/vuoto calcolato CON il micro-moto e disegnato lungo
   tutto l'asse verticale z della nuvola tomografica, come richiesto: alto =
   pieno (caldo), basso = vuoto (freddo). Resta un discriminante
   multi-attributo dichiarato come tale, non una rilevazione di cavita'
   risolta in profondita' (ch08-ch09 delle fonti). */
const NV=D.voxels.length;
const VXY=new Float32Array(NV*2), VZ=new Float32Array(NV), VSRC=new Int32Array(NV);
function drawVoxels(step){
  const zlo=_zlo, zhi=_zhi;
  let m=0;
  for(let k=0;k<NV;k+=step){
    const v=D.voxels[k];
    if(v[3]<V.thr||v[2]<zlo||v[2]>zhi) continue;
    if(V.onlyPyr && v[5]!==1) continue;
    const p=project(v[0],v[1],v[2]);
    VXY[2*m]=p[0]; VXY[2*m+1]=p[1]; VZ[m]=p[2]; VSRC[m]=k; m++;
  }
  const ord=new Array(m);
  for(let q=0;q<m;q++) ord[q]=q;
  ord.sort((x,y)=>VZ[y]-VZ[x]);
  let last=null;
  for(let q=0;q<m;q++){
    const w=ord[q], sol=D.voxels[VSRC[w]][4];
    const c=RAMP.sol(sol);
    const key='rgba('+c[0]+','+c[1]+','+c[2]+','+(0.22+0.65*sol).toFixed(2)+')';
    if(key!==last){ ctx.fillStyle=key; last=key; }
    ctx.beginPath(); ctx.arc(VXY[2*w],VXY[2*w+1],V.vsize,0,6.2832); ctx.fill();
  }
  return m;
}

/* ---------- piramidi ideali ---------- */
function drawPyramids(){
  ctx.font='500 11px "Fira Sans",sans-serif';
  for(const p of D.pyramids){
    const Vt=p.vertices.map(v=>project(v[0],v[1],v[2]));
    const faces=[[0,1,4],[1,2,4],[2,3,4],[3,0,4]];
    const ord=faces.map((f,i)=>[i,(Vt[f[0]][2]+Vt[f[1]][2]+Vt[f[2]][2])/3])
                   .sort((a,b)=>b[1]-a[1]);
    ctx.strokeStyle='rgba(245,158,11,.85)'; ctx.lineWidth=1.3;
    for(const row of ord){
      const f=faces[row[0]];
      ctx.beginPath(); ctx.moveTo(Vt[f[0]][0],Vt[f[0]][1]);
      ctx.lineTo(Vt[f[1]][0],Vt[f[1]][1]); ctx.lineTo(Vt[f[2]][0],Vt[f[2]][1]);
      ctx.closePath(); ctx.stroke();
    }
    const ap=project(p.apex[0],p.apex[1],p.apex[2]);
    ctx.fillStyle='rgba(248,250,252,.95)';
    ctx.fillText(p.name, ap[0]+8, ap[1]-5);
    ctx.fillStyle='rgba(148,163,184,.9)';
    ctx.font='400 10px "Fira Code",monospace';
    ctx.fillText(p.base_side.toFixed(0)+'m x '+p.height.toFixed(0)+'m', ap[0]+8, ap[1]+8);
    ctx.font='500 11px "Fira Sans",sans-serif';
  }
}

/* ---------- suolo (DEM, F38/F39) ----------------------------------------
   Come piramide_acustica_vh.py: nuvola di punti rossi, colore FISSO (non una
   rampa per quota, a differenza degli altri livelli) perche' e' esattamente
   lo stile del programma sorgente. F39: i punti seguono anche le facce delle
   piramidi, quindi la nuvola sale fino all'apice invece di fermarsi al
   plateau. E' un RIFERIMENTO, non una misura di questa pipeline: la geometria
   .xml resta l'unica usata nel calcolo. */
let DEM_SPRITE=null, DEM_CELL=0;
function buildDemSprite(){
  const r=Math.max(1.0, V.dot*0.65);
  const cell=Math.ceil(2*r+2), dev=Math.ceil(cell*dpr);
  const cvs=document.createElement('canvas');
  cvs.width=dev; cvs.height=dev;
  const c=cvs.getContext('2d');
  c.fillStyle='#dc2626';
  c.beginPath(); c.arc(dev/2,dev/2, r*(dev/cell), 0, 6.2832); c.fill();
  DEM_SPRITE=cvs; DEM_CELL=cell;
}
function drawDem(step){
  if(!DEM || !DNN) return 0;
  if(!DEM_SPRITE) buildDemSprite();
  const zlo=_zlo, zhi=_zhi, secLimit=_secLimit;
  const half=DEM_CELL/2;
  let n=0;
  for(let ai=0; ai<DNL; ai+=step){
    for(let aj=0; aj<DNP; aj+=step){
      const a=ai*DNP+aj;
      const h=DHGT[a];
      if(!(h>=zlo && h<=zhi)) continue;
      if(V.secAxis==='E' ? DEAST[a]>secLimit : DNORTH[a]>secLimit) continue;
      ctx.drawImage(DEM_SPRITE, DX[a]-half, DY[a]-half, DEM_CELL, DEM_CELL);
      n++;
    }
  }
  return n;
}

/* ---------- disegno completo ---------- */
let nDots=0, nVox=0, nMM=0, nDem=0, lastMs=0;
/* F30: quanti nodi la vista puo' arrivare a disegnare. La nuvola costa una
   drawImage per punto, quindi con tutta la piana in campo (8 000 punti) un
   fotogramma sta sui 40 ms: durante il trascinamento e la rotazione si
   decima di un fattore due per riga e colonna e si torna a pochi ms, mentre
   a scena ferma si disegna tutto. Il criterio guarda i nodi CANDIDATI, non
   quelli disegnati nel fotogramma precedente, altrimenti la densita'
   oscillerebbe fra un fotogramma e l'altro. */
let N_PYR=0; for(let i=0;i<NN;i++) if(SIMM[i]) N_PYR++;
function render(){
  const t0=performance.now();
  const cand = V.onlyPyr ? N_PYR : NN;
  const step = ((dragging||V.autoRot) && cand>4000) ? 2 : 1;
  updateCam();
  updateCuts();
  ctx.clearRect(0,0,W,H);
  ctx.fillStyle='#020617'; ctx.fillRect(0,0,W,H);
  if(V.grid) drawGrid();
  /* i nodi si proiettano una volta per campo di quote e per fotogramma:
     HGT per la superficie misurata (reticolo NN), XHREF per i punti di
     micro-moto che poggiano sul datum .xml (reticolo XNN, F33: proprio e
     completo, non quello della nuvola misurata; F37: non piu' disegnato
     come strato a se', vedi sopra). */
  if(V.mmp) projectNodes(XHREF,XEAST,XNORTH,XNN,RX,RY,RZ);
  if(V.surf) projectNodes(HGT,EAST,NORTH,NN,PX,PY,PZ);
  if(V.dem && DEM) projectNodes(DHGT,DEAST,DNORTH,DNN,DX,DY,DZ);
  nDem = V.dem ? drawDem(step) : 0; // il suolo va sotto a tutto il resto
  nMM = V.mmp ? drawMM(step) : 0;
  nDots = V.surf ? drawSurface(step) : 0;
  nVox = V.vox ? drawVoxels(step) : 0;
  if(V.pyr) drawPyramids();
  if(V.str) drawStructures();
  if(V.axes) drawAxes();
  lastMs = performance.now()-t0;
  hud();
  dirty=false;
}
function hud(){
  const zlo=_zlo, zhi=_zhi;
  const area = D.area_m ? (D.area_m[0].toFixed(0)+'x'+D.area_m[1].toFixed(0)+' m') : '-';
  document.getElementById('hud').innerHTML =
    'superficie <b>'+nDots+'</b> punti &nbsp; suolo DEM <b>'+nDem+'</b> &nbsp; micro-moto <b>'+nMM+'</b>'+
    ' &nbsp; voxel <b>'+nVox+'</b>'+
    ' &nbsp; <b>'+lastMs.toFixed(0)+' ms</b>/fotogramma<br>'+
    'quota <b>'+zlo.toFixed(0)+' .. '+zhi.toFixed(0)+' m</b> &nbsp; area <b>'+area+'</b><br>'+
    'date '+(D.pol||'VH')+' <b>'+D.dates.length+'</b> &nbsp; master <b>'+D.master+'</b><br>'+
    'sigma_h <b>'+D.budget.sigma_h+' m</b> &nbsp; delta_z <b>'+D.budget.delta_z_vertical+' m</b>';
}
function loop(){
  if(V.autoRot){ V.yaw+=0.0026; if(V.yaw>2*Math.PI)V.yaw-=2*Math.PI; syncYaw(); dirty=true; }
  if(dirty) render();
  requestAnimationFrame(loop);
}

/* ---------- interazione ---------- */
const $=id=>document.getElementById(id);
let drag=null, lx=0, ly=0;
cv.addEventListener('mousedown',e=>{
  drag = (e.shiftKey||e.button===1) ? 'pan' : 'rot';
  dragging=true; lx=e.clientX; ly=e.clientY; cv.classList.add('drag');
  V.autoRot=false; syncPlay(); cv.focus(); invalidate();
});
addEventListener('mouseup',()=>{
  if(drag){ drag=null; dragging=false; cv.classList.remove('drag'); invalidate(); }
});
addEventListener('mousemove',e=>{
  if(!drag) return;
  const dx=e.clientX-lx, dy=e.clientY-ly; lx=e.clientX; ly=e.clientY;
  if(drag==='pan'){ V.panX+=dx; V.panY+=dy; }
  else{
    V.yaw+=dx*0.006;
    V.pitch=Math.max(-1.48,Math.min(1.48,V.pitch+dy*0.006));
    syncYaw(); syncPitch();
  }
  invalidate();
});
cv.addEventListener('wheel',e=>{
  e.preventDefault();
  V.zoom=Math.max(0.3,Math.min(8,V.zoom*(e.deltaY>0?0.92:1.086)));
  $('zoom').value=V.zoom; $('zoomV').textContent=V.zoom.toFixed(2)+'x';
  invalidate();
},{passive:false});
cv.addEventListener('keydown',e=>{
  const big=e.shiftKey, k=e.key, stp=0.06, pan=24;
  if(k==='ArrowLeft'){ if(big){V.panX-=pan;}else{V.yaw-=stp;syncYaw();} }
  else if(k==='ArrowRight'){ if(big){V.panX+=pan;}else{V.yaw+=stp;syncYaw();} }
  else if(k==='ArrowUp'){ if(big){V.panY-=pan;}else{V.pitch=Math.min(1.48,V.pitch+stp);syncPitch();} }
  else if(k==='ArrowDown'){ if(big){V.panY+=pan;}else{V.pitch=Math.max(-1.48,V.pitch-stp);syncPitch();} }
  else if(k==='+'||k==='='){ V.zoom=Math.min(8,V.zoom*1.1); $('zoom').value=V.zoom; }
  else if(k==='-'||k==='_'){ V.zoom=Math.max(0.3,V.zoom/1.1); $('zoom').value=V.zoom; }
  else if(k==='q'||k==='Q'){ V.roll-=0.05; $('rol').value=V.roll*180/Math.PI; }
  else if(k==='e'||k==='E'){ V.roll+=0.05; $('rol').value=V.roll*180/Math.PI; }
  else if(k==='r'||k==='R'){ resetView(); }
  else if(k===' '){ V.autoRot=!V.autoRot; syncPlay(); }
  else return;
  invalidate(); e.preventDefault();
});

function syncYaw(){ const d=((V.yaw*180/Math.PI)%360+360)%360;
  $('yaw').value=d; $('yawV').textContent=d.toFixed(0)+'°'; }
function syncPitch(){ const d=V.pitch*180/Math.PI;
  $('pit').value=d; $('pitV').textContent=d.toFixed(0)+'°'; }
function syncPlay(){ $('play').classList.toggle('on',V.autoRot);
  $('play').textContent = V.autoRot?'Rotazione':'Ferma'; invalidate(); }
function resetView(){
  V.yaw=34*Math.PI/180; V.pitch=26*Math.PI/180; V.roll=0;
  V.zoom=1; V.panX=0; V.panY=0;
  syncYaw(); syncPitch();
  $('rol').value=0; $('rolV').textContent='0°';
  $('zoom').value=1; $('zoomV').textContent='1.00x';
  invalidate();
}
function bindRange(id,key,fmt,transform,onSet){
  const el=$(id), out=$(id+'V');
  const set=function(){ const raw=+el.value; V[key]=transform?transform(raw):raw;
                        out.textContent=fmt(raw); if(onSet)onSet(); invalidate(); };
  el.addEventListener('input',set); set();
}
bindRange('yaw','yaw',v=>v.toFixed(0)+'°',v=>v*Math.PI/180);
bindRange('pit','pitch',v=>v.toFixed(0)+'°',v=>v*Math.PI/180);
bindRange('rol','roll',v=>v.toFixed(0)+'°',v=>v*Math.PI/180);
bindRange('zoom','zoom',v=>v.toFixed(2)+'x');
bindRange('ex','zEx',v=>v.toFixed(1)+'x',null,invalidateShade);
bindRange('op','opacity',v=>(v*100).toFixed(0)+'%');
bindRange('dot','dot',v=>v.toFixed(1)+' px',null,invalidateSprite);
bindRange('mmt','mmThr',v=>v.toFixed(2));
bindRange('saz','sunAz',v=>v.toFixed(0)+'°',null,invalidateShade);
bindRange('sel','sunEl',v=>v.toFixed(0)+'°',null,invalidateShade);
$('thr').value=V.thr; bindRange('thr','thr',v=>v.toFixed(1)+' dB');
bindRange('sz','vsize',v=>v.toFixed(1));
bindRange('zlo','zlo',v=>(D.z_min+(D.z_max-D.z_min)*v/100).toFixed(0)+' m',v=>v/100);
bindRange('zhi','zhi',v=>(D.z_min+(D.z_max-D.z_min)*v/100).toFixed(0)+' m',v=>v/100);
bindRange('sec','sec',v=>v>=100?'nessuna':(v+'%'),v=>v/100);

$('play').onclick=()=>{V.autoRot=!V.autoRot;syncPlay();};
$('reset').onclick=resetView;
$('proj').onclick=e=>{V.persp=!V.persp;
  e.target.textContent=V.persp?'Prospettiva':'Ortografica';
  e.target.classList.toggle('on',!V.persp); invalidate();};
$('secAx').onclick=e=>{V.secAxis=V.secAxis==='E'?'N':'E';
  e.target.textContent='Sezione: '+(V.secAxis==='E'?'Est':'Nord'); invalidate();};
const toggles={lSurf:'surf',lVox:'vox',lPyr:'pyr',lGrid:'grid',
               lAxes:'axes',lMask:'onlyGood',lOnlyPyr:'onlyPyr',
               lMM:'mmp',lStr:'str',lDem:'dem'};
for(const id in toggles){
  (function(tid,key){
    $(tid).onclick=e=>{V[key]=!V[key];
      e.target.classList.toggle('on',V[key]); invalidate();};
  })(id,toggles[id]);
}
/* F38/F39: senza layer del suolo il pulsante resta disattivato invece che
   finto-attivo. Con F39 il layer c'e' anche senza rete (terreno = bilineare
   .xml + piramidi), quindi qui si arriva solo con cfg.suolo_dem=False. */
if(!DEM){ const b=$('lDem'); b.classList.remove('on'); b.disabled=true; b.title='layer del suolo non presente in questa uscita'; }
document.querySelectorAll('.attr').forEach(b=>{
  b.onclick=()=>{
    V.attr=b.dataset.a;
    document.querySelectorAll('.attr').forEach(x=>x.classList.toggle('on',x===b));
    $('attrV').textContent=b.textContent.toLowerCase();
    invalidateColor(); legend();
  };
});
function legend(){
  const lg=$('lg'); lg.innerHTML='';
  const ramp=RAMP[V.attr];
  for(let i=0;i<44;i++){
    const s=document.createElement('i');
    const c=ramp(i/43); s.style.background='rgb('+c[0]+','+c[1]+','+c[2]+')';
    lg.appendChild(s);
  }
  const lo=RANGE[V.attr][0], hi=RANGE[V.attr][1], lab=ATTR_LABEL[V.attr];
  $('lgA').textContent = V.attr==='h' ? lo.toFixed(0)+' m' : lab[0];
  $('lgB').textContent = V.attr==='h' ? hi.toFixed(0)+' m' : lab[1];
}
/* F33: rampa pieno/vuoto dei voxel, fissa (non dipende da V.attr, che sceglie
   solo il colore della nuvola di superficie), costruita una volta sola. */
function legendVox(){
  const lg=$('lgVox'); if(!lg) return;
  for(let i=0;i<44;i++){
    const s=document.createElement('i');
    const c=RAMP.sol(i/43); s.style.background='rgb('+c[0]+','+c[1]+','+c[2]+')';
    lg.appendChild(s);
  }
}

/* ---------- pannelli ---------- */
const MMD=D.mm||{};
$('mt').innerHTML =
  '<tr><td>quote .xml sul ritaglio</td><td>'+
    (D.href_range?D.href_range[0]+' .. '+D.href_range[1]+' m':'-')+'</td></tr>'+
  '<tr><td>banda osservabile</td><td>'+
    (MMD.banda_hz?MMD.banda_hz[0]+' .. '+MMD.banda_hz[1]+' Hz':'-')+'</td></tr>'+
  '<tr><td>righe del banco</td><td>'+
    ((MMD.righe||[]).filter(v=>v>0).map(v=>{
       const ok = MMD.banda_hz && v>=MMD.banda_hz[0] && v<=MMD.banda_hz[1];
       return '<span style="color:'+(ok?'var(--accent-2)':'var(--muted)')+'">'+
              v.toFixed(1)+'</span>';
     }).join(' ')||'-')+' Hz</td></tr>'+
  '<tr><td>concentrazione mediana</td><td>'+(MMD.conc_mediana!=null?MMD.conc_mediana:'-')+'</td></tr>'+
  '<tr><td>rumore bianco</td><td>'+(MMD.rumore!=null?MMD.rumore:'-')+'</td></tr>'+
  '<tr><td colspan="2" style="padding-top:6px;color:var(--muted);text-align:left">'+
    'Il banco ha N<sub>D</sub> = 12 marce, quindi sei righe distinte in tutto e '+
    'una sola dentro la banda osservabile (in azzurro). I punti poggiano sulla '+
    'superficie .xml: il micro-moto e\' un attributo di superficie, il banco di '+
    'sub-aperture non gli assegna una profondita\'.'+
  '</td></tr>';

const ST=D.structures||[];
$('st').innerHTML = ST.length
  ? ST.map(s=>'<tr><td><span style="color:'+s.colore+'">&#9632;</span> '+
      s.pyr.split(' ')[0]+' '+s.num+'</td><td style="text-align:left;color:var(--muted)">'+
      s.nome+'</td></tr>').join('')+
    '<tr><td colspan="2" style="padding-top:6px;color:var(--muted);text-align:left">'+
      'Geometria da piramide_cheope_3d.py e piramide_kefren_3d.py (Petrie 1883; '+
      'Legon, Vyse &amp; Perring, Lehner). E\' un RIFERIMENTO archeologico: con '+
      '&delta;z = '+D.budget.delta_z_vertical+' m e camere di pochi metri, nulla '+
      'di tutto questo e\' rilevabile da questi dati.</td></tr>'
  : '<tr><td colspan="2">moduli di geometria interna non caricati</td></tr>';

const b=D.budget;
$('bt').innerHTML =
  '<tr><td>area processata</td><td>'+(D.area_m?D.area_m[0].toFixed(0)+' x '+D.area_m[1].toFixed(0)+' m':'-')+'</td></tr>'+
  '<tr><td>celle</td><td>'+(D.n_cells||NN)+'</td></tr>'+
  '<tr><td>date impilate</td><td>'+b.n_baselines+'</td></tr>'+
  '<tr><td>escursione baseline</td><td>'+b.b_spread+' m</td></tr>'+
  '<tr><td>risoluzione &delta;z</td><td>'+b.delta_z_vertical+' m</td></tr>'+
  '<tr><td>precisione &sigma;h</td><td>'+b.sigma_h+' m</td></tr>'+
  '<tr><td>altezza di ambiguita\'</td><td>'+b.ambiguity+' m</td></tr>'+
  '<tr><td>altezza Cheope</td><td>'+b.target_height+' m</td></tr>'+
  '<tr><td>celle sul target</td><td>'+b.cells+'</td></tr>'+
  '<tr><td>segno di k<sub>z</sub></td><td>'+(D.sign>0?'+1':'-1')+'</td></tr>';

/* ---------- profilo dei lobi verticali (F41) -----------------------------
   Perche' la nuvola tomografica ha punti anche SOPRA il suolo. Due ragioni,
   e questo grafico mostra la seconda. (1) L'asse z e' un intervallo di
   RICERCA simmetrico attorno al datum (+-elev_max): nulla nell'inversione
   obbliga la soluzione a stare sottoterra, meta' dell'asse sta in aria per
   costruzione. (2) Dove finisce l'energia lo decide la PSF verticale
   dell'array di baseline, che e' la curva qui sotto: dB in orizzontale (0 =
   picco), quota relativa al datum in verticale, orientata come nella scena.
   Se la curva non attraversa mai la riga dei -3 dB non esiste un lobo
   principale separabile, e ogni struttura verticale nella nuvola e' un lobo
   dell'array, non stratigrafia. La nuvola NON viene ritagliata a z<=0: sarebbe
   un'immagine piu' credibile e piu' falsa.                                */
const LB = D.lobi;
if(LB && LB.z && LB.z.length){
  const cvp=$('psf'), g=cvp.getContext('2d');
  const dp=Math.min(2, window.devicePixelRatio||1);
  const wCss=Math.max(220, cvp.clientWidth||304), hCss=250;
  cvp.width=Math.round(wCss*dp); cvp.height=Math.round(hCss*dp);
  cvp.style.height=hCss+'px';
  g.setTransform(dp,0,0,dp,0,0);
  const F='9px ui-monospace,Consolas,monospace';
  const ML=42, MR=12, MT=14, MB=26;
  const W=wCss-ML-MR, H=hCss-MT-MB;
  const NZ=LB.z.length, zmin=LB.z[0], zmax=LB.z[NZ-1];
  let lo=0;
  for(let k=0;k<NZ;k++){
    if(LB.pyr[k]<lo) lo=LB.pyr[k];
    if(LB.tutto[k]<lo) lo=LB.tutto[k];
  }
  lo=Math.min(-4, Math.floor(lo)-1);
  const X=d=>ML+W*(1-Math.max(lo,Math.min(0,d))/lo);
  const Y=z=>MT+H*(1-(z-zmin)/(zmax-zmin));

  /* la meta' dell'asse che sta in aria, tinta appena: e' quella la domanda */
  g.fillStyle='rgba(245,158,11,.07)';
  g.fillRect(ML, Y(zmax), W, Y(0)-Y(zmax));

  g.font=F; g.strokeStyle='#1E2740'; g.lineWidth=1;
  g.fillStyle='#94A3B8'; g.textAlign='right'; g.textBaseline='middle';
  for(const zv of [zmin, zmin/2, 0, zmax/2, zmax]){
    const y=Y(zv);
    g.beginPath(); g.moveTo(ML,y); g.lineTo(ML+W,y); g.stroke();
    g.fillText(zv.toFixed(0), ML-5, y);
  }
  g.textAlign='center'; g.textBaseline='top';
  const passoDb = lo<=-10 ? 4 : 2;
  for(let d=0; d>=lo; d-=passoDb){
    const x=X(d);
    g.beginPath(); g.moveTo(x,MT); g.lineTo(x,MT+H); g.stroke();
    g.fillText(d.toFixed(0), x, MT+H+5);
  }
  g.textAlign='left'; g.fillText('dB', ML+2, MT+H+15);
  g.textAlign='right'; g.textBaseline='bottom'; g.fillText('m', ML-5, MT-1);

  /* la riga di meta' potenza: il criterio, disegnato invece che raccontato */
  if(lo < -3){
    g.save(); g.setLineDash([4,3]); g.strokeStyle='#DC2626'; g.lineWidth=1;
    const x3=X(-3);
    g.beginPath(); g.moveTo(x3,MT); g.lineTo(x3,MT+H); g.stroke(); g.restore();
    g.fillStyle='#F87171'; g.textAlign='left'; g.textBaseline='top';
    g.fillText('-3 dB', x3+3, MT+2);
  }

  const curva=(arr,col,lw)=>{
    g.beginPath();
    for(let k=0;k<NZ;k++){
      const x=X(arr[k]), y=Y(LB.z[k]);
      if(k) g.lineTo(x,y); else g.moveTo(x,y);
    }
    g.strokeStyle=col; g.lineWidth=lw; g.lineJoin='round'; g.stroke();
  };
  curva(LB.tutto,'#64748B',1);
  curva(LB.pyr,'#38BDF8',1.8);

  /* il datum: i punti sopra questa riga sono in aria */
  g.save(); g.setLineDash([3,3]); g.strokeStyle='#F59E0B'; g.lineWidth=1.1;
  g.beginPath(); g.moveTo(ML,Y(0)); g.lineTo(ML+W,Y(0)); g.stroke(); g.restore();
  g.fillStyle='#F59E0B'; g.textAlign='right'; g.textBaseline='bottom';
  g.fillText('superficie', ML+W-3, Y(0)-2);

  /* la cella di Rayleigh, in scala sullo stesso asse: si vede subito che e'
     piu' larga della piramide che dovrebbe contenere */
  if(LB.delta_z>0){
    const zc=LB.z_picco;
    const z1=Math.max(zmin, zc-LB.delta_z/2), z2=Math.min(zmax, zc+LB.delta_z/2);
    const xb=ML+10;
    g.strokeStyle='#F8FAFC'; g.lineWidth=1.2; g.beginPath();
    g.moveTo(xb,Y(z1)); g.lineTo(xb,Y(z2));
    g.moveTo(xb-4,Y(z1)); g.lineTo(xb+4,Y(z1));
    g.moveTo(xb-4,Y(z2)); g.lineTo(xb+4,Y(z2));
    g.stroke();
    g.fillStyle='#F8FAFC'; g.textAlign='left'; g.textBaseline='bottom';
    g.fillText('δz '+LB.delta_z+' m', xb-4, Y(z2)-3);
  }

  g.fillStyle='#38BDF8';
  g.beginPath(); g.arc(X(0),Y(LB.z_picco),2.8,0,6.2832); g.fill();
  g.strokeStyle='#334155'; g.lineWidth=1; g.strokeRect(ML+.5,MT+.5,W,H);

  const pc=v=>(100*v).toFixed(1)+' %';
  $('pt').innerHTML =
    '<tr><td>picco del profilo</td><td>'+LB.z_picco+' m</td></tr>'+
    '<tr><td>contrasto picco/minimo</td><td>'+LB.contrasto_db+' dB</td></tr>'+
    '<tr><td>livello ai bordi dell\'asse</td><td>'+LB.bordi_db[0]+' / '+
      LB.bordi_db[1]+' dB</td></tr>'+
    '<tr><td>larghezza a -3 dB</td><td>'+
      (LB.lobo_troncato?'&gt; ':'')+LB.lobo_larghezza+' m</td></tr>'+
    '<tr><td>lobo laterale peggiore</td><td>'+
      (LB.lobo_laterale_db==null?'nessun fuori-lobo':LB.lobo_laterale_db+' dB')+
      '</td></tr>'+
    '<tr><td>&delta;z di Rayleigh</td><td>'+LB.delta_z+' m</td></tr>'+
    '<tr><td>passo dell\'asse z</td><td>'+LB.passo_z+' m'+
      (LB.sovracamp?' ('+LB.sovracamp+'&times;)':'')+'</td></tr>'+
    '<tr><td>voxel disegnati</td><td>'+LB.n_voxel+'</td></tr>'+
    '<tr><td>di cui sopra la superficie</td><td>'+pc(LB.frac_sopra)+'</td></tr>'+
    '<tr><td>di cui sopra l\'apice</td><td>'+pc(LB.frac_sopra_apice)+'</td></tr>';

  $('psfnote').innerHTML =
    (LB.lobo_troncato
      ? 'Il fascio <b>non scende mai a -3 dB</b> dentro l\'asse: non esiste '+
        'un lobo principale separabile. '
      : 'Il lobo principale a -3 dB e\' largo '+LB.lobo_larghezza+' m (da '+
        LB.lobo_lo+' a '+LB.lobo_hi+' m), ma il <b>lobo laterale peggiore '+
        'sta a '+LB.lobo_laterale_db+' dB</b>, appena '+
        Math.abs(LB.lobo_laterale_db).toFixed(1)+' dB sotto il picco: un '+
        'picco che i suoi lobi quasi eguagliano non localizza niente in '+
        'profondita\'. ')+
    'Su tutti gli '+(LB.z[LB.z.length-1]-LB.z[0]).toFixed(0)+' m dell\'asse '+
    'il contrasto e\' <b>'+LB.contrasto_db+' dB</b>. E\' per questo che il '+
    pc(LB.frac_sopra)+' dei voxel disegnati finisce sopra il suolo e il '+
    pc(LB.frac_sopra_apice)+' sopra l\'apice di Cheope: sono lobi dell\'array '+
    'di baseline, non struttura. La cella di Rayleigh vale &delta;z = '+
    LB.delta_z+' m contro i '+D.budget.target_height+' m di Cheope, e l\'asse '+
    'e\' campionato a '+LB.passo_z+' m: la nuvola sembra dettagliata perche\' '+
    'e\' sovracampionata '+(LB.sovracamp||'?')+' volte, non perche\' contenga '+
    'quel dettaglio. Non viene ritagliata a z&le;0 di proposito: nasconderne '+
    'meta\' darebbe un\'immagine piu\' credibile e piu\' falsa. La correzione '+
    'vera e\' piu\' baseline.';
}


const R=D.regression||{};
$('rt').innerHTML = ('pendenza' in R)
  ? '<tr><td>celle di qualita\'</td><td>'+R.pixel+'</td></tr>'+
    '<tr><td>pendenza (attesa 1.0)</td><td>'+R.pendenza+' &plusmn; '+R.errore_standard_pendenza+'</td></tr>'+
    (R.pendenza_theil_sen!=null?'<tr><td>pendenza robusta</td><td>'+R.pendenza_theil_sen+'</td></tr>':'')+
    '<tr><td>correlazione r</td><td>'+R.correlazione_r+'</td></tr>'+
    '<tr><td>rmse</td><td>'+R.rmse_m+' m</td></tr>'+
    '<tr><td>soglia &gamma; (nullo)</td><td>'+D.gamma_thr+' / '+D.null_median+'</td></tr>'+
    '<tr><td>celle sopra soglia</td><td>'+(D.good_frac*100).toFixed(1)+' %</td></tr>'+
    '<tr><td colspan="2" style="padding-top:6px;color:var(--fg)">'+R.esito+'</td></tr>'
  : '<tr><td colspan="2">'+(R.esito||'regressione non disponibile')+'</td></tr>';

$('vt').innerHTML = D.validation.length
  ? '<tr><td colspan="2" style="color:var(--muted);padding-bottom:4px">sommita\': misurato vs simulato</td></tr>'+
    D.validation.map(v=>
      '<tr><td>'+v.nome.split(' ')[0]+'<br><span style="font-size:10px">simulato '+v.atteso+' m</span></td>'+
      '<td>'+v.misurato+' m<br><span class="pill '+(v.ok?'g':'r')+'">'+
      (v.errore>0?'+':'')+v.errore+' m</span></td></tr>').join('')
  : '<tr><td colspan="2">nessuna cella di sommita\' sopra la soglia di qualita\'</td></tr>';

/* F39: pannello del suolo. Mostra, per ciascuna piramide, la quota .xml letta
   dai file dello stack, la base ricavata dal terreno e l'apice che ne segue,
   piu' lo scarto dal valore di letteratura -- dichiarato, non corretto. */
{
  const SU = D.suolo_dem;
  const el = $('suot');
  if(!SU || !SU.info){
    el.innerHTML = '<tr><td colspan="2" style="color:var(--muted)">'+
      'layer del suolo non presente in questa uscita</td></tr>';
  } else {
    const I = SU.info, P = SU.plateau || [];
    const f = (v,n) => (v==null ? '-' : Number(v).toFixed(n==null?1:n));
    el.innerHTML =
      '<tr><td colspan="2" style="color:var(--muted);padding-bottom:4px">'+
        'datum .xml sul ritaglio: '+f(I.datum_xml_m&&I.datum_xml_m[0])+' .. '+
        f(I.datum_xml_m&&I.datum_xml_m[1])+' m, da '+(SU.n_date_lette||0)+
        ' annotation.xml'+(I.rilievo_esterno
          ? ' &nbsp; rilievo esterno riportato di '+f(-I.bias_datum_m)+' m'
          : ' &nbsp; senza rilievo esterno')+'</td></tr>'+
      P.map((d,k) =>
        '<tr><td>'+d.nome+'</td><td>base '+f(I.base_alt_m&&I.base_alt_m[k])+
        ' m &nbsp; apice '+f(I.apici_m&&I.apici_m[k])+' m &nbsp; <span style="color:var(--muted)">('+
        '.xml '+f(d.h_xml_master)+' m, letteratura '+f(d.base_alt_letteratura,0)+
        ' m, scarto '+(I.scarto_letteratura_m&&I.scarto_letteratura_m[k]>=0?'+':'')+
        f(I.scarto_letteratura_m&&I.scarto_letteratura_m[k])+' m)</span></td></tr>'
      ).join('')+
      '<tr><td colspan="2" style="padding-top:6px;color:var(--muted);text-align:left">'+
        'suolo disegnato: '+f(I.suolo_m&&I.suolo_m[0],0)+' .. '+f(I.suolo_m&&I.suolo_m[1],0)+
        ' m (terreno da solo '+f(I.terreno_m&&I.terreno_m[0],0)+' .. '+
        f(I.terreno_m&&I.terreno_m[1],0)+' m).</td></tr>';
  }
}

/* F35: nodi VERI (non spline-interpolati) della geolocation grid dello .xml,
   ordinati per distanza dal ritaglio delle piramidi. Provano che il quasi
   piano di "xml_ref" dentro il ritaglio e' un limite di risoluzione del
   reticolo (231 nodi su tutta la scena), non un errore di calcolo: il
   rilievo vero, altrove nella scena, arriva fino a GRAW.h_range[1] m. */
const GRAW = D.gcp_raw || {h:[], dist_km:[], line:[], pixel:[]};
{
  const idxs = GRAW.h.map((_, i) => i)
    .sort((a, b) => GRAW.dist_km[a] - GRAW.dist_km[b])
    .slice(0, 8);
  $('gcpt').innerHTML =
    '<tr><td colspan="2" style="color:var(--muted);padding-bottom:4px">'+
      'reticolo: '+(D.gcp_raw_spacing_px?D.gcp_raw_spacing_px.join(' x '):'-')+
      ' px &nbsp; rilievo vero sulla scena: '+
      (D.gcp_raw_h_range?D.gcp_raw_h_range[0]+' .. '+D.gcp_raw_h_range[1]:'-')+
      ' m (&sigma; '+(D.gcp_raw_h_std!=null?D.gcp_raw_h_std:'-')+' m)</td></tr>'+
    idxs.map(i =>
      '<tr><td>linea '+GRAW.line[i]+', px '+GRAW.pixel[i]+'</td><td>'+
      GRAW.dist_km[i].toFixed(1)+' km &nbsp; '+GRAW.h[i].toFixed(0)+' m</td></tr>'
    ).join('') +
    '<tr><td colspan="2" style="padding-top:6px;color:var(--muted);text-align:left">'+
      'Nodo piu\' vicino al ritaglio: '+(D.gcp_raw_near_km!=null?D.gcp_raw_near_km:'-')+
      ' km, quota '+(D.gcp_raw_near_h!=null?D.gcp_raw_near_h:'-')+' m.</td></tr>';
}

$('note').className = b.resolves ? 'note ok' : 'note';
$('note').innerHTML = b.resolves
  ? 'Risoluzione sufficiente per separare strutture interne.'
  : '<b>Risoluzione e precisione non sono la stessa cosa.</b> '+
    '&delta;z = '+b.delta_z_vertical+' m e\' la separazione minima fra <i>due</i> '+
    'diffusori nella stessa cella: nessuna struttura interna e\' separabile. '+
    '&sigma;h = '+b.sigma_h+' m e\' l\'errore sulla quota di <i>un</i> diffusore '+
    'dominante, ed e\' quello che rende misurabile la superficie qui disegnata. '+
    'L\'attributo &laquo;micro-moto&raquo; e\' un discriminante, non una rilevazione '+
    'di cavita\'. I voxel sono colorati per indice di pieno/vuoto (energia del '+
    'periodogramma &times; coerenza &times; (1 &minus; micro-moto)), che incorpora '+
    'il micro-moto ma resta lo stesso tipo di discriminante, non una rilevazione '+
    'risolta in profondita\'. '+D.n_fixes+' correzioni applicate rispetto alla versione precedente.';

/* ---------- avvio ---------- */
resize(); resetView(); legend(); legendVox();
if(reduce){ V.autoRot=false; }
syncPlay(); loop();
</script>
</body>
</html>
"""


# ==========================================================================
# 13.  Autotest  (F21)
# ==========================================================================

def _esito(superato: bool) -> bool:
    """Stampa l'esito di un controllo dell'autotest e lo restituisce."""
    print("ok" if superato else "FALLITO")
    return superato


def selftest() -> int:
    """Blocca le convenzioni che, se invertite, ribaltano silenziosamente i
    risultati: segno dello shift, segno del tracker, interpolatore orbitale."""
    ok = True
    rng = np.random.default_rng(7)

    # 1) apply_shift: apply_shift(f, d)(x) == f(x - d)
    n = 64
    x = np.arange(n)
    f = np.exp(-((x - 20.0) ** 2) / 8.0).astype(np.complex64)
    img = np.outer(f, np.ones(n)).astype(np.complex64)
    sh = apply_shift(img, complex(5.0, 0.0))
    peak = int(np.argmax(np.abs(sh[:, 0])))
    print(f"  apply_shift(+5) sposta il picco 20 -> {peak}", end="  ")
    ok &= _esito(peak == 25)

    # 2) tracker: slave(x) = master(x - d)  =>  ritorna +d
    m = rng.standard_normal((1, 32, 32)) + 1j * rng.standard_normal((1, 32, 32))
    m = m.astype(np.complex64)
    d = complex(3.0, -2.0)
    s = apply_shift(m[0], d)[None].astype(np.complex64)
    est = complex(_batch_subpixel_shift(m, s, refine="dft")[0])
    print(f"  tracker: atteso {d}, stimato ({est.real:+.2f},{est.imag:+.2f})", end="  ")
    ok &= _esito(abs(est - d) <= 0.25)

    # 3) Lagrange contro un'orbita polinomiale nota
    t = np.arange(0.0, 100.0, 10.0)
    coef = rng.standard_normal((4, 3)) * 1000.0
    pos = np.stack([np.polyval(coef[:, c], t) for c in range(3)], axis=1)
    vel = np.stack([np.polyval(np.polyder(coef[:, c]), t) for c in range(3)], axis=1)
    orb = Orbit(t=t, pos=pos, vel=vel, t0=datetime(2026, 1, 1), order=8)
    tq = np.array([37.3])
    p_hat, _ = orb.state(tq)
    p_true = np.array([np.polyval(coef[:, c], tq[0]) for c in range(3)])
    err = float(np.linalg.norm(p_hat[0] - p_true))
    lin = np.array([np.interp(tq[0], t, pos[:, c]) for c in range(3)])
    err_lin = float(np.linalg.norm(lin - p_true))
    print(f"  orbita: Lagrange err={err:.3e} m  vs  lineare err={err_lin:.1f} m", end="  ")
    ok &= _esito(err <= 1e-6 and err_lin > err)

    # 4) banco di sub-aperture dentro lo spettro, con Doppler rate realistico
    from piramidi_v01 import Burst as _B

    class _A:
        azimuth_bandwidth = 313.0
        prf = 1451.6
        lines_per_burst = 1508
        azimuth_time_interval = 2.0555e-3
        slant_range_time = 5.6419e-3
        range_sampling_rate = 6.4345e7
        azimuth_steering_rate_deg = 0.97986
        orbit_velocity = 7596.9
        wavelength = 0.0555
        bursts = [_B(0, 0, 1507, "", 0.0, (-2188.0,), 0.0, (0.0,), 0.0)]
    # _A e' uno stub duck-typed con i soli campi usati da plan_subapertures
    plan = plan_subapertures(_A(), 1024, Config(),  # type: ignore[arg-type]
                             burst_idx=0, pixel=2100.0)
    hi = -plan.b_cd / 2 + (plan.n_d - 1) * plan.step + plan.b_sub + plan.b_shift
    print(f"  sub-aperture: bordo slave {hi:+.2f} Hz  <= {plan.b_cd / 2:+.2f} Hz", end="  ")
    ok &= _esito(hi <= plan.b_cd / 2 + 1e-6)

    # 5) layover: il confronto e' con l'INCIDENZA, non con 90 - incidenza
    lay = layover_report(39.4)
    near = {r["piramide"]: r["faccia_vicina"] for r in lay["piramidi"]}
    look = lay["angolo_di_vista_off_nadir_deg"]
    print(f"  layover: soglia {lay['soglia_layover_deg']} gradi, "
          f"vista off-nadir {look} gradi", end="  ")
    ok &= _esito(all(v == "layover" for v in near.values())
                 and 30.0 < look < 40.0)

    # 6) F39: il profilo del suolo sulle piramidi. Due invarianti che, se
    #    rotte, fanno galleggiare o affondare il layer senza dare errore:
    #    l'apice sta a base + altezza, e la rotazione applicata qui e'
    #    l'inversa ESATTA di quella con cui pyramid_mesh() pone gli spigoli.
    lat0, lon0 = PYRAMIDS[0].lat, PYRAMIDS[0].lon
    base = [40.0, 41.0, 42.0]
    ep, np_ = _enu_offset(np.array([p.lat for p in PYRAMIDS]),
                          np.array([p.lon for p in PYRAMIDS]), lat0, lon0)
    z_apex = pyramid_profile_enu(ep, np_, lat0, lon0, base)
    att = np.array([b + p.height_m for b, p in zip(base, PYRAMIDS)])
    e_apex = float(np.max(np.abs(z_apex - att)))
    # spigoli della base da pyramid_mesh(): il profilo deve valere la base
    corner_err = 0.0
    for k, p in enumerate(PYRAMIDS):
        v = np.array(pyramid_mesh(p, lat0, lon0)["vertices"][:4], dtype=float)
        zc = pyramid_profile_enu(v[:, 0], v[:, 1], lat0, lon0, base)
        corner_err = max(corner_err, float(np.max(np.abs(zc - base[k]))))
    print(f"  suolo (F39): apice err={e_apex:.3e} m, spigoli di base "
          f"err={corner_err:.3e} m", end="  ")
    ok &= _esito(e_apex <= 1e-6 and corner_err <= 1e-6)

    # 7) F42: la finestra dello slave va PRESA all'offset, non traslata dopo.
    #    apply_shift() e' circolare: chiedergli l'inquadramento del prodotto
    #    (migliaia di linee) non riallinea nulla, ripiega il ritaglio su se
    #    stesso. Qui il "prodotto" e' sintetico e l'offset e' noto: la lettura
    #    all'offset deve ricostruire il master, la traslazione circolare no.
    n_r, n_c, off = 400, 48, -240
    y_g = np.arange(n_r)[:, None]
    x_g = np.arange(n_c)[None, :]
    prod_m = np.exp(-((y_g - 300.0) ** 2) / 50.0
                    - ((x_g - 24.0) ** 2) / 50.0).astype(np.complex64)
    prod_s = np.roll(prod_m, off, axis=0)          # stesso terreno, altro frame
    w0, w1 = 290, 310
    chip_m = prod_m[w0:w1 + 1]

    def _corr(a: np.ndarray) -> float:
        u = (a - a.mean()).ravel()
        v = (chip_m - chip_m.mean()).ravel()
        d = np.linalg.norm(u) * np.linalg.norm(v)
        return float(abs(np.vdot(v, u)) / d) if d > 0 else 0.0

    c_vecchio = _corr(apply_shift(prod_s[w0:w1 + 1], complex(-off, 0.0)))
    c_nuovo = _corr(prod_s[w0 + off:w1 + off + 1])
    print(f"  F42: finestra all'offset r={c_nuovo:.3f}  vs  "
          f"traslazione circolare r={c_vecchio:.3f}", end="  ")
    ok &= _esito(c_nuovo > 0.999 and c_vecchio < 0.1)

    print("\n  " + ("TUTTI GLI AUTOTEST SUPERATI" if ok else "AUTOTEST FALLITI"))
    return 0 if ok else 1


# ==========================================================================
# 14.  CLI
# ==========================================================================

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Tomografia 3D multi-baseline VH della piana di Giza con "
                    "ricostruzione della superficie reale.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--stack-dir", default=DEFAULT_STACK)
    ap.add_argument("--swath", default="iw2", choices=["iw1", "iw2", "iw3"])
    ap.add_argument("--platform", default="s1c")
    ap.add_argument("--pol", dest="polarisation", default="vh", choices=["vh", "vv"],
                    help="canale usato per i calcoli (vedere F26 sul costo di VH)")
    ap.add_argument("--dates", dest="n_dates", type=int, default=11)
    ap.add_argument("--n-elev", type=int, default=257)
    ap.add_argument("--elev-max", dest="elev_max_m", type=float, default=400.0)
    ap.add_argument("--nd", dest="n_d", type=int, default=12)
    ap.add_argument("--b-shift", dest="b_shift_fraction", type=float, default=0.09,
                    help="B_shift / B_cD: selettore della frequenza meccanica")
    ap.add_argument("--mm-lines", type=int, default=1024)
    ap.add_argument("--look-range", type=int, default=4)
    ap.add_argument("--look-azimuth", type=int, default=1)
    ap.add_argument("--gamma-min", type=float, default=0.0,
                    help="F44: pavimento AGGIUNTIVO sotto la soglia calibrata "
                         "sul nullo (F22). Va lasciato a 0: la soglia la "
                         "decide la distribuzione nulla, non una costante")
    ap.add_argument("--margin", dest="area_margin_m", type=float, default=150.0,
                    help="margine [m] attorno alle piramidi in geometria radar (F28)")
    ap.add_argument("--full-scene", action="store_true",
                    help="ripristina la finestra larga su tutta la piana")
    ap.add_argument("--out", dest="out_dir", default="out_piramidi_v02")
    ap.add_argument("--html", dest="html_name", default="tomografia_piramidi_3d.html")
    ap.add_argument("--report-only", action="store_true",
                    help="solo baseline e budget, senza leggere i .tiff")
    ap.add_argument("--no-calibrazione", dest="calibrazione",
                    action="store_false",
                    help="F40: non applicare calibration-*.xml e noise-*.xml "
                         "(ampiezze in conteggi DN grezzi)")
    ap.add_argument("--no-suolo", dest="suolo_dem", action="store_false",
                    help="F39: non produrre il layer 'suolo (DEM)'")
    ap.add_argument("--no-dem-esterno", dest="fetch_dem", action="store_false",
                    help="F39: niente rete, il terreno resta la bilineare .xml")
    ap.add_argument("--no-dem-piramidi", dest="dem_pyramids", action="store_false",
                    help="F39: il suolo non segue il profilo delle piramidi")
    ap.add_argument("--dem-grid", dest="dem_grid", type=int, default=24,
                    help="F39: nodi per lato del reticolo di query del DEM esterno")
    ap.add_argument("--fixes", action="store_true", help="elenca le correzioni applicate")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    print("=" * 76)
    print("piramidi_v02 - tomografia multi-baseline VH + superficie reale")
    print("=" * 76)

    if args.fixes:
        print()
        for i, g, d in FIXES:
            print(f"  {i}  [{g}]")
            for line in _wrap(d, 70):
                print(f"        {line}")
            print()
        return 0

    if args.selftest:
        print("\nautotest delle convenzioni (F21):")
        return selftest()

    cfg = Config(
        stack_dir=args.stack_dir, swath=args.swath, platform=args.platform,
        polarisation=args.polarisation,
        calibrazione=args.calibrazione,
        n_dates=args.n_dates, n_elev=args.n_elev, elev_max_m=args.elev_max_m,
        n_d=args.n_d, b_shift_fraction=args.b_shift_fraction,
        mm_lines=args.mm_lines, look_range=args.look_range,
        look_azimuth=args.look_azimuth, gamma_min=args.gamma_min,
        area_margin_m=args.area_margin_m, full_scene=args.full_scene,
        out_dir=args.out_dir, html_name=args.html_name,
        suolo_dem=args.suolo_dem, fetch_dem=args.fetch_dem,
        dem_pyramids=args.dem_pyramids, dem_grid=args.dem_grid,
    )

    if args.report_only:
        entries = discover_stack(cfg)[: cfg.n_dates]
        anns = [parse_annotation(e.annotation) for e in entries]
        orbits = [read_orbit(e.annotation) for e in entries]
        tgt = ecef_from_llh(PYRAMIDS[0].lat, PYRAMIDS[0].lon, PYRAMIDS[0].base_alt_m)
        probe = compute_baselines(list(zip([e.date for e in entries], orbits)), tgt, 0)
        mi = pick_supermaster(probe)
        bl = compute_baselines(list(zip([e.date for e in entries], orbits)), tgt, mi)
        print(f"\nstack VH: {len(entries)} acquisizioni  "
              f"supermaster {entries[mi].date} (F13)")
        for b in bl:
            print(f"  {b.date}: B_perp = {b.b_perp:+8.2f} m   "
                  f"B_par = {b.b_par:+9.2f} m   dt = {b.b_temp:+6.0f} d")
        budget = compute_tomo_budget(bl, anns[mi], anns[mi].incidence_mid,
                                     PYRAMIDS[0].height_m)
        print(budget.as_text())
        geo = Geocoder(anns[mi])
        win, burst = target_window(anns[mi], geo, cfg, n_lines=cfg.mm_lines)
        print(plan_subapertures(anns[mi], win.n_l, cfg, burst_idx=burst,
                                pixel=0.5 * (win.p0 + win.p1)).as_text())
        lay = layover_report(anns[mi].incidence_mid)
        print(f"\nGEOMETRIA (F16): incidenza {lay['angolo_di_incidenza_deg']} gradi, "
              f"vista off-nadir {lay['angolo_di_vista_off_nadir_deg']} gradi")
        print(f"  soglia layover {lay['soglia_layover_deg']} gradi, "
              f"soglia ombra {lay['soglia_ombra_deg']} gradi")
        for r in lay["piramidi"]:
            print(f"  {r['piramide']:24s} faccia {r['pendenza_faccia_deg']:5.2f} gradi"
                  f" -> vicina: {r['faccia_vicina']:14s} lontana: {r['faccia_lontana']}"
                  f"  (compressione {r['compressione']}x)")
        return 0

    # la cartella di uscita deve esistere PRIMA di run(): ground_dem_suolo()
    # ci scrive dentro la cache del DEM mentre run() e' ancora in corso
    os.makedirs(cfg.out_dir, exist_ok=True)
    res = run(cfg, verbose=not args.quiet)

    print("\n  [10] validazione a tre livelli")
    valid = validate(res, cfg)
    lvl1 = valid["livello_1_geometria"]
    print(f"      quota mediana della piana: {lvl1['_quota_mediana_piana_m']} m")
    reg = lvl1["regressione_misurato_vs_simulato"]
    if "pendenza" in reg:
        print(f"      regressione misurato vs simulato (F25): pendenza "
              f"{reg['pendenza']:+.3f} +- {reg['errore_standard_pendenza']:.3f} "
              f"(attesa 1.0, t = {reg['t_pendenza']}), r = {reg['correlazione_r']}, "
              f"rmse = {reg['rmse_m']} m su {reg['pixel']} celle")
        if "pendenza_theil_sen" in reg:
            print(f"      pendenza robusta (Theil-Sen) = "
                  f"{reg['pendenza_theil_sen']:+.3f}  IC95 {reg['ic95']}")
        for f in reg.get("mediane_per_fascia", []):
            hi = f["simulato_m"][1]
            print(f"        quota simulata {f['simulato_m'][0]:3d}-"
                  f"{(str(hi) + ' m') if hi else 'oltre':>6s}: "
                  f"{f['celle']:4d} celle, misurato mediano "
                  f"{f['misurato_mediano_m']:+6.1f} m")
        print(f"      -> {reg['esito']}")
    else:
        print(f"      regressione: {reg['esito']}")
    for p in PYRAMIDS:
        d = lvl1.get(p.name, {})
        if "rilievo_misurato_m" in d:
            flag = "entro 2 sigma" if d["entro_2_sigma_h"] else "FUORI da 2 sigma"
            print(f"      {p.name:24s} sommita': misurato "
                  f"{d['rilievo_misurato_m']:+7.1f} m  simulato "
                  f"{d['rilievo_simulato_m']:+7.1f} m  ({flag}, "
                  f"{d['celle_di_sommita']} celle)")
        else:
            print(f"      {p.name:24s} {d.get('esito', 'celle di qualita insufficienti')}"
                  f"  ({d.get('celle_di_qualita', 0)}/{d.get('celle_simulate', 0)} celle)")
    l3 = valid["livello_3_struttura"]["separabilita_piramide_vs_piana"]
    if l3:
        print(f"      separabilita' piramidi/piana: "
              f"{l3['delta_quota_mediana_m']:+.1f} m  (z = {l3['z_score']})")

    print("\n  [11] profili verticali e mappe")
    profiles = profile_analysis(res, cfg)
    surface_plot(res, cfg, valid)

    print("\n  [12] uscite")
    save_outputs(res, cfg, valid, profiles)
    path = build_html(res, cfg, valid)

    print(f"\n  HTML 3D interattivo : {os.path.abspath(path)}")
    print(f"  altri output in     : {os.path.abspath(cfg.out_dir)}")
    print(f"  completato in {float(res['elapsed_s']):.1f} s")
    return 0


def _wrap(text: str, width: int) -> List[str]:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return lines


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramide_kefren_3d.py
=====================
PIRAMIDE DI KEFREN (KHAFRE) - MODELLO 3D e **UNICA SORGENTE** delle quote
delle strutture interne note usate dal resto della pipeline.

Sistema di riferimento (lo stesso del modello di Cheope):
    X = Est (+) / Ovest (-)      [m]
    Y = Nord (+) / Sud (-)       [m]
    Z = quota rispetto al piano di base della piramide [m]

Il file ha DUE usi:

  * eseguito da solo (`python piramide_kefren_3d.py`) stampa gli ingombri e
    disegna il modello 3D;
  * importato come modulo espone `STRUTTURE` (la geometria, dato puro) e
    `riferimenti()` (gli INGOMBRI nel formato usato da piramide_unificato.py /
    piramide_acustica_vh.py).

Da qui, e SOLO da qui, piramide_acustica_vh.py prende le strutture interne di
Chefren: le vecchie tabelle copiate a mano dentro gli altri script sono state
cancellate, cosi' i numeri stanno in un posto solo e il disegno 3D e l'overlay
dei grafici SAR non possono piu' divergere.

Fonti delle quote: Legon (GM 110, Table I) per gli sviluppi dei passaggi;
Vyse & Perring 1837, Belzoni 1818, Lehner "The Complete Pyramids" per le
camere. Gli elementi la cui geometria non e' certa lo dichiarano nel campo
`certezza`.
"""
import math
from piramide_3d_comune import (riferimenti as _riferimenti,
                                stampa_riferimenti, guscio_piramide,
                                disegna_box, disegna_corridoio,
                                legenda_senza_doppioni)

# ============================================================
# 1. STRUTTURA ESTERNA
# ============================================================
BASE = 214.5       # m, lato di base (valore NOVA)
H = 143.5          # m, altezza ORIGINARIA ricostruita. L'altezza ATTUALE
                   # (~136 m) e' quella che vede il DEM: per la superficie
                   # "DEM + piramide" la pipeline usa il valore del preset,
                   # non questo.
L = BASE / 2.0

# Pendenza teorica ricavata dalle quote esterne
FACE_ANGLE = math.degrees(math.atan(H / L))

PIRAMIDE = dict(nome="Chefren / Khafre", base_m=BASE, h_m=H,
                face_angle_deg=FACE_ANGLE)

# ============================================================
# 2. ASSE INTERNO
# ============================================================
# Le due entrate sono sul lato nord e circa 12 m a est della mezzeria N-S.
# Questo offset e' ben attestato e vale per TUTTI i tratti dei due corridoi.
X_OFFSET = 12.0

# ============================================================
# 3. GEOMETRIA DEI PASSAGGI (calcolata, non trascritta)
# ============================================================
# Dati geometrici principali da Legon (GM 110, Table I):
#   Lower entrance:        34.94 m @ 21 gradi 40'
#   Lower horizontal:       7.88 m
#   Lower chamber passage:  6.71 m @ 21 gradi 19'
#   Upper horizontal:      39.37 m
#   Upper entrance:        36.95 m @ 26 gradi 28'
def _delta(length, angle_deg):
    a = math.radians(angle_deg)
    return length * math.cos(a), length * math.sin(a)


# Entrata superiore: quota 11.54 m sopra il piano di base, sulla FACCIA nord
# (non sullo spigolo di base: a quella quota la faccia rientra di
# 11.54/tan(FACE_ANGLE) m).
Z_UPPER = 11.54
Y_UPPER = L - Z_UPPER / math.tan(math.radians(FACE_ANGLE))

_dy, _dz = _delta(36.95, 26 + 28 / 60)
Y_JOIN = Y_UPPER - _dy
Z_JOIN = Z_UPPER - _dz

# Percorso orizzontale principale verso la camera funeraria (39.37 m, Legon)
Y_HORIZ_END = Y_JOIN - 39.37

# Entrata inferiore: nel lastricato, 9.33 m oltre la base (Legon), a quota 0
Y_LOWER = L + 9.33
Z_LOWER = 0.0
_dy, _dz = _delta(34.94, 21 + 40 / 60)
Y_LOW_END = Y_LOWER - _dy
Z_LOW_END = Z_LOWER - _dz
Y_LOW_H1 = Y_LOW_END - 7.88
_dy, _dz = _delta(6.71, 21 + 19 / 60)
Y_LOW_RISE = Y_LOW_H1 - _dy
Z_LOW_RISE = Z_LOW_END + _dz

# ============================================================
# 4. CAMERE
# ============================================================
# Camera funeraria: 14.15 (E-O) x 5.0 (N-S) m, scavata nel substrato roccioso,
# tetto a doppio spiovente con apice a 6.83 m dal pavimento (pareti 5.24 m).
CHAMBER_Z = -3.73
CHAMBER_Y = -5.0
CHAMBER_DX = 14.15
CHAMBER_DY = 5.00
CHAMBER_H_PARETI = 5.24
CHAMBER_H_COLMO = 6.83

# Camera secondaria (funzione incerta: serdab / deposito / corredo)
SUB_Y = 17.0
SUB_Z = -3.0

SEZIONE_W = 1.05     # larghezza standard dei corridoi
SEZIONE_H = 1.20     # altezza standard dei corridoi

# ============================================================
# 5. LE STRUTTURE - DATO UNICO, usato sia dal disegno sia da riferimenti()
# ============================================================
# tipo "box"      : (x, y, z=quota del pavimento, dx, dy, dz)
# tipo "cuspide"  : box + colmo (dz_colmo) per il tetto a doppio spiovente
# tipo "corridoio": (x, y0, z0, y1, z1, w, h), asse lungo Y,
#                   z0/z1 = quota del PAVIMENTO agli estremi
STRUTTURE = [
    dict(num=1, nome="Ingresso inferiore (lastricato, a nord)", tipo="corridoio",
         x=X_OFFSET, y0=Y_LOWER, z0=Z_LOWER, y1=Y_LOWER - 2.0, z1=Z_LOWER,
         w=SEZIONE_W, h=SEZIONE_H, colore="darkgreen",
         certezza="posizione attestata (Legon: 9.33 m oltre la base)"),
    dict(num=2, nome="Passaggio inferiore (34.94 m @ 21 gradi 40')", tipo="corridoio",
         x=X_OFFSET, y0=Y_LOWER, z0=Z_LOWER, y1=Y_LOW_END, z1=Z_LOW_END,
         w=SEZIONE_W, h=SEZIONE_H, colore="forestgreen",
         certezza="sviluppo e pendenza da Legon, Table I"),
    dict(num=3, nome="Passaggio orizzontale inferiore (7.88 m)", tipo="corridoio",
         x=X_OFFSET, y0=Y_LOW_END, z0=Z_LOW_END, y1=Y_LOW_H1, z1=Z_LOW_END,
         w=SEZIONE_W, h=SEZIONE_H, colore="limegreen",
         certezza="sviluppo da Legon, Table I"),
    dict(num=4, nome="Raccordo inferiore in salita (6.71 m @ 21 gradi 19')",
         tipo="corridoio",
         x=X_OFFSET, y0=Y_LOW_H1, z0=Z_LOW_END, y1=Y_LOW_RISE, z1=Z_LOW_RISE,
         w=SEZIONE_W, h=SEZIONE_H, colore="seagreen",
         certezza="sviluppo e pendenza da Legon, Table I"),
    dict(num=5, nome="Camera secondaria (funzione incerta)", tipo="box",
         x=X_OFFSET - 6.0, y=SUB_Y, z=SUB_Z, dx=10.4, dy=3.0, dz=2.8,
         colore="mediumpurple",
         certezza="INCERTA: dimensioni indicative, funzione non accertata"),
    dict(num=6, nome="Ingresso superiore (facciata nord, +11.54 m)", tipo="corridoio",
         x=X_OFFSET, y0=Y_UPPER, z0=Z_UPPER, y1=Y_UPPER - 2.0,
         z1=Z_UPPER - 2.0 * math.tan(math.radians(26 + 28 / 60)),
         w=SEZIONE_W, h=SEZIONE_H, colore="orangered",
         certezza="quota 11.54 m attestata; imbocco sulla FACCIA, non sullo spigolo"),
    dict(num=7, nome="Passaggio superiore (36.95 m @ 26 gradi 28')", tipo="corridoio",
         x=X_OFFSET, y0=Y_UPPER, z0=Z_UPPER, y1=Y_JOIN, z1=Z_JOIN,
         w=SEZIONE_W, h=SEZIONE_H, colore="darkorange",
         certezza="sviluppo e pendenza da Legon, Table I"),
    dict(num=8, nome="Galleria orizzontale superiore (39.37 m)", tipo="corridoio",
         x=X_OFFSET, y0=Y_JOIN, z0=Z_JOIN, y1=Y_HORIZ_END, z1=Z_JOIN,
         w=SEZIONE_W, h=SEZIONE_H, colore="gold",
         certezza="sviluppo da Legon, Table I"),
    dict(num=9, nome="Collegamento alla camera funeraria", tipo="corridoio",
         x=X_OFFSET, y0=Y_HORIZ_END, z0=Z_JOIN,
         y1=CHAMBER_Y + 2.5, z1=CHAMBER_Z + 1.0,
         w=SEZIONE_W, h=SEZIONE_H, colore="royalblue",
         certezza="SEMPLIFICATO: il rilievo reale ha altri segmenti e cambi di quota"),
    dict(num=10, nome="Camera funeraria (tetto a doppio spiovente)", tipo="cuspide",
         x=X_OFFSET, y=CHAMBER_Y, z=CHAMBER_Z,
         dx=CHAMBER_DX, dy=CHAMBER_DY, dz=CHAMBER_H_PARETI,
         dz_colmo=CHAMBER_H_COLMO, colore="crimson",
         certezza="14.15 x 5.0 m misurati; quota del pavimento semplificata"),
    dict(num=11, nome="Sarcofago in granito", tipo="box",
         x=X_OFFSET + 2.8, y=CHAMBER_Y - 0.5, z=CHAMBER_Z + 0.05,
         dx=2.15, dy=0.95, dz=0.70, colore="dimgray",
         certezza="volume semplificato"),
    dict(num=12, nome="Vano/cassone canopico (interpretazione)", tipo="box",
         x=X_OFFSET - 3.4, y=CHAMBER_Y + 0.9, z=CHAMBER_Z + 0.02,
         dx=1.00, dy=0.70, dz=0.35, colore="saddlebrown",
         certezza="INTERPRETAZIONE, non un rilievo"),
]


# ============================================================
# 6. INGOMBRI (bounding box) PER IL RESTO DELLA PIPELINE
# ============================================================
def riferimenti():
    """Gli ingombri di questa piramide nel formato di piramide_unificato.py /
    piramide_acustica_vh.py.

    Il calcolo sta in piramide_3d_comune.riferimenti(): era identico parola
    per parola nei due script 3D, insieme a _bbox(). Qui resta solo il legame
    con la tabella STRUTTURE di questo file, che e' l'unica cosa che cambia
    fra le due piramidi."""
    return _riferimenti(STRUTTURE)


# ============================================================
# 7. DISEGNO (solo se eseguito direttamente)
# ============================================================
#: trasparenze delle due primitive, unica cosa che distingueva questo disegno
#: da quello di Cheope nelle righe che ora stanno in piramide_3d_comune
ALPHA_CAMERA = 0.75
ALPHA_CORRIDOIO = 0.7


def _disegna():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("Planimetria 3D della Piramide di Kefren (quote principali)")

    guscio_piramide(ax, L, H)

    for s in STRUTTURE:
        if s["tipo"] == "corridoio":
            disegna_corridoio(ax, s["x"], s["y0"], s["z0"], s["y1"], s["z1"],
                              s["w"], s["h"], s["colore"], s["nome"],
                              ALPHA_CORRIDOIO)
        else:
            disegna_box(ax, s["x"], s["y"], s["z"],
                        (s["dx"], s["dy"], s["dz"]), s["colore"], s["nome"],
                        ALPHA_CAMERA)
            if s["tipo"] == "cuspide":
                # tetto a doppio spiovente, disegnato a parte
                zb = s["z"] + s["dz"]; za = s["z"] + s["dz_colmo"]
                hy = s["dy"] / 2.0
                x0, x1 = s["x"] - s["dx"]/2, s["x"] + s["dx"]/2
                roof = [[[x0, s["y"]-hy, zb], [x1, s["y"]-hy, zb],
                         [x1, s["y"], za], [x0, s["y"], za]],
                        [[x0, s["y"], za], [x1, s["y"], za],
                         [x1, s["y"]+hy, zb], [x0, s["y"]+hy, zb]]]
                ax.add_collection3d(Poly3DCollection(
                    roof, facecolors="lightgray", linewidths=0.4,
                    edgecolors="black", alpha=0.55,
                    label="Tetto a doppio spiovente"))

    ax.set_xlabel("Asse Est - Ovest (X) [m]")
    ax.set_ylabel("Asse Nord - Sud (Y) [m]")
    ax.set_zlabel("Quota rispetto al piano di base [m]")
    ax.set_xlim([-120, 120]); ax.set_ylim([-130, 130]); ax.set_zlim([-40, 150])
    ax.view_init(elev=22, azim=-45)
    legenda_senza_doppioni(ax)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.show()


if __name__ == "__main__":
    print(f"Angolo della faccia calcolato: {FACE_ANGLE:.2f} gradi")
    print(f"Ingresso superiore sulla faccia nord: y = {Y_UPPER:.2f} m, z = {Z_UPPER:.2f} m")
    print(f"Giunzione del passaggio superiore:    y = {Y_JOIN:.2f} m, z = {Z_JOIN:.2f} m")
    stampa_riferimenti(STRUTTURE)
    _disegna()

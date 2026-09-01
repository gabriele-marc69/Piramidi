#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramide_cheope_3d.py
=====================
PIRAMIDE DI CHEOPE / KHUFU - RICOSTRUZIONE 3D e **UNICA SORGENTE** delle quote
delle strutture interne note usate dal resto della pipeline.

Sistema di riferimento:
    X = Est (+) / Ovest (-)      [m]
    Y = Nord (+) / Sud (-)       [m]
    Z = quota rispetto al piano di base della piramide [m]

Il file ha DUE usi:

  * eseguito da solo (`python piramide_cheope_3d.py`) stampa i controlli
    geometrici e disegna il modello 3D;
  * importato come modulo espone `STRUTTURE` (la geometria, dato puro) e
    `riferimenti()` (gli INGOMBRI nel formato usato da piramide_unificato.py /
    piramide_acustica_vh.py).

Da qui, e SOLO da qui, piramide_acustica_vh.py prende le strutture interne di
Cheope: le vecchie tabelle copiate a mano dentro gli altri script sono state
cancellate, cosi' i numeri stanno in un posto solo e il disegno 3D e l'overlay
dei grafici SAR non possono piu' divergere.

IMPORTANTE: alcuni ambienti reali hanno geometrie piu' complesse di un semplice
parallelepipedo (Grande Galleria, tetto della Camera della Regina, camere di
scarico). Qui sono rappresentati gli INGOMBRI geometrici principali, non ogni
singolo blocco costruttivo. Le cavita' individuate via muografia dal progetto
ScanPyramids (Big Void 2017, North Face Corridor 2023) NON sono modellate: non
hanno una geometria pubblicata con la stessa precisione del rilievo Petrie.
"""
import math

# ================================================================
# 1. PIRAMIDE ESTERNA
# ================================================================
BASE = 230.30
H = 146.60         # altezza ORIGINARIA. L'altezza ATTUALE (~138.7 m) e' quella
                   # che vede il DEM: per la superficie "DEM + piramide" la
                   # pipeline usa il valore del preset, non questo.
L = BASE / 2.0

# Angolo della faccia esterna derivato dalle quote del modello
FACE_ANGLE = math.degrees(math.atan(H / L))

PIRAMIDE = dict(nome="Cheope / Khnum-Khufu", base_m=BASE, h_m=H,
                face_angle_deg=FACE_ANGLE)

# ================================================================
# 2. ASSE INTERNO
# ================================================================
# Offset E-O documentato: 287 pollici = 7.29 m verso EST della mezzeria N-S
# (Petrie 1883). NON 7.9 m: "15 cubiti reali" e' una lettura sbagliata,
# 287" = 13.9 cubiti.
X_OFFSET = 7.29

SEZIONE_W = 1.05     # larghezza standard dei corridoi
SEZIONE_H = 1.20     # altezza standard dei corridoi

# ================================================================
# 3. CAMERE (misure Petrie 1883)
# ================================================================
# CAMERA DEL RE: 10.47 m E-O x 5.24 m N-S x 5.85 h; pavimento +42.95 m.
# Il lato LUNGO corre EST-OVEST. Centro N-S = -11.02 m: Petrie misura la
# parete nord a 330.6 +/- 0.8 pollici (8.40 m) a SUD del centro della
# piramide -> centro = -8.40 - 5.24/2.
RE_Y, RE_Z = -11.02, 42.95
RE_L, RE_W, RE_H = 10.47, 5.24, 5.85

# CAMERA DELLA REGINA: 5.23 x 5.75 m, quota +21.19 m, ingombro semplificato
# fino al colmo del tetto a doppio spiovente.
REGINA_Y, REGINA_Z = 0.0, 21.19
REGINA_L, REGINA_W, REGINA_H = 5.23, 5.75, 6.26

# CAMERA SOTTERRANEA: pavimento -30.01 m, 14.06 x 8.28 m, incompiuta.
SOTTER_Y, SOTTER_Z = -5.27, -30.01
SOTTER_L, SOTTER_W, SOTTER_H = 14.06, 8.28, 3.50

# ANTICAMERA / PORTCULLIS: Petrie 116.30 (N-S) x 65.00 (E-O) x 149.35 (h)
# POLLICI = 2.95 x 1.65 x 3.79 m. NON 5.6 m: quella e' la stessa misura in
# CUBITI reali (116.30" = 5.64 cubiti) scambiata per metri.
ANTICAMERA_Y, ANTICAMERA_Z = -4.36, 42.95
ANTICAMERA_L, ANTICAMERA_W, ANTICAMERA_H = 1.65, 2.95, 3.79

# ================================================================
# 4. CORRIDOI
# ================================================================
# CORRIDOIO DISCENDENTE: ingresso sulla FACCIA nord a +16.97 m (non sullo
# spigolo di base: a quella quota la faccia, pendenza atan(146.6/115.15) =
# 51.85 gradi, rientra di 13.34 m -> y = 115.15 - 13.34 = 101.83 m);
# 105.2 m @ 26 gradi 31'23" -> fine inclinazione a (Y=7.71, Z=-30.01).
DISC_START_Y, DISC_START_Z = 101.83, 16.97
DISC_END_Y, DISC_END_Z = 7.71, -30.01

# Breve tratto orizzontale finale (~8.84 m) verso la parete nord della
# Camera Sotterranea.
SOTTER_NORTH_WALL_Y = SOTTER_Y + SOTTER_W / 2.0

# CORRIDOIO ASCENDENTE: 39.28 m @ 26 gradi 02'30", da (77.46, +3.95) a
# (42.16, +21.19).
ASC_START_Y, ASC_START_Z = 77.46, 3.95
ASC_END_Y, ASC_END_Z = 42.16, 21.19

# GRANDE GALLERIA: 46.68 m di SVILUPPO @ 26 gradi 02'30" -> 41.94 m
# orizzontali e 20.49 m di dislivello, da (42.16, +21.19) a (0.22, +41.68).
# Ingombro semplificato: 2.06 m di larghezza, 8.60 m di altezza.
GALL_START_Y, GALL_START_Z = 42.16, 21.19
GALL_END_Y, GALL_END_Z = 0.22, 41.68
GALL_W, GALL_H = 2.06, 8.60

# CORRIDOIO ORIZZONTALE DELLA REGINA: 39.29 m a quota costante +21.19 m,
# dallo svincolo (Y=42.16) alla parete nord della Camera della Regina (Y=2.88).
REG_START_Y, REG_END_Y, REG_Z = 42.16, 2.88, 21.19
REG_W, REG_H = 1.05, 1.15

# ================================================================
# 5. LE STRUTTURE - DATO UNICO, usato sia dal disegno sia da riferimenti()
# ================================================================
# tipo "box"      : (x, y, z=quota del pavimento, dx, dy, dz)
# tipo "corridoio": (x, y0, z0, y1, z1, w, h), asse lungo Y,
#                   z0/z1 = quota del PAVIMENTO agli estremi
STRUTTURE = [
    dict(num=1, nome="Entrata originale (faccia nord, +16.97 m)", tipo="corridoio",
         x=X_OFFSET, y0=DISC_START_Y, z0=DISC_START_Z,
         y1=DISC_START_Y - 2.0,
         z1=DISC_START_Z - 2.0 * math.tan(math.radians(26 + 31/60 + 23/3600)),
         w=SEZIONE_W, h=SEZIONE_H, colore="darkgreen",
         certezza="quota 668.3 pollici = 16.97 m (Petrie); imbocco sulla FACCIA"),
    dict(num=2, nome="Corridoio discendente (105.2 m @ 26 gradi 31'23\")",
         tipo="corridoio",
         x=X_OFFSET, y0=DISC_START_Y, z0=DISC_START_Z,
         y1=DISC_END_Y, z1=DISC_END_Z, w=SEZIONE_W, h=SEZIONE_H,
         colore="forestgreen", certezza="rilievo Petrie 1883"),
    dict(num=3, nome="Passaggio orizzontale finale (~8.84 m)", tipo="corridoio",
         x=X_OFFSET, y0=DISC_END_Y, z0=DISC_END_Z,
         y1=SOTTER_NORTH_WALL_Y, z1=DISC_END_Z, w=SEZIONE_W, h=SEZIONE_H,
         colore="darkseagreen", certezza="rilievo Petrie 1883"),
    dict(num=4, nome="Camera sotterranea (-30.01 m, incompiuta)", tipo="box",
         x=X_OFFSET, y=SOTTER_Y, z=SOTTER_Z,
         dx=SOTTER_L, dy=SOTTER_W, dz=SOTTER_H, colore="saddlebrown",
         certezza="14.06 x 8.28 m (Petrie 553.5x325.9\"); altezza modellata ~3.5 m"),
    dict(num=5, nome="Corridoio ascendente (39.28 m @ 26 gradi 02'30\")",
         tipo="corridoio",
         x=X_OFFSET, y0=ASC_START_Y, z0=ASC_START_Z,
         y1=ASC_END_Y, z1=ASC_END_Z, w=SEZIONE_W, h=SEZIONE_H,
         colore="darkorange", certezza="rilievo Petrie 1883"),
    dict(num=6, nome="Corridoio orizzontale della Regina (39.29 m)",
         tipo="corridoio",
         x=X_OFFSET, y0=REG_START_Y, z0=REG_Z, y1=REG_END_Y, z1=REG_Z,
         w=REG_W, h=REG_H, colore="gold", certezza="rilievo Petrie 1883"),
    dict(num=7, nome="Camera della Regina (+21.19 m)", tipo="box",
         x=X_OFFSET, y=REGINA_Y, z=REGINA_Z,
         dx=REGINA_L, dy=REGINA_W, dz=REGINA_H, colore="darkorchid",
         certezza="5.23 x 5.75 m (Petrie 205.9x226.5\"), colmo ~245\""),
    dict(num=8, nome="Grande Galleria (46.68 m @ 26 gradi 02'30\")",
         tipo="corridoio",
         x=X_OFFSET, y0=GALL_START_Y, z0=GALL_START_Z,
         y1=GALL_END_Y, z1=GALL_END_Z, w=GALL_W, h=GALL_H,
         colore="royalblue",
         certezza="46.68 m di SVILUPPO = 41.94 m orizzontali; ingombro semplificato"),
    dict(num=9, nome="Anticamera / Portcullis (+42.95 m)", tipo="box",
         x=X_OFFSET, y=ANTICAMERA_Y, z=ANTICAMERA_Z,
         dx=ANTICAMERA_L, dy=ANTICAMERA_W, dz=ANTICAMERA_H, colore="rosybrown",
         certezza="Petrie 116.30 x 65.00 x 149.35 pollici"),
    dict(num=10, nome="Camera del Re (+42.95 m)", tipo="box",
         x=X_OFFSET, y=RE_Y, z=RE_Z, dx=RE_L, dy=RE_W, dz=RE_H,
         colore="crimson",
         certezza="10.47 (E-O) x 5.24 (N-S) x 5.85 m (Petrie 412.66\"); "
                  "parete nord a 330.6\" a sud del centro"),
]

# Il Great Step (Grande Gradino) non e' un ingombro ma un PUNTO di controllo:
# Petrie lo colloca sulla mezzeria N-S della piramide, ed e' esattamente dove
# arriva la Grande Galleria (y = +0.22 m). Disegnato, non esportato.
GREAT_STEP = (X_OFFSET, GALL_END_Y, GALL_END_Z)


# ================================================================
# 6. INGOMBRI (bounding box) PER IL RESTO DELLA PIPELINE
# ================================================================
def _bbox(s):
    """Ingombro di una struttura -> (x_c, y_c, z_c, dx, dy, dz) nel sistema
    X=Est, Y=Nord, Z=quota sopra la base. Per i corridoi inclinati e' il
    bounding box del tratto, comprensivo dello spessore della sezione."""
    if s["tipo"] in ("box", "cuspide"):
        dz = s.get("dz_colmo", s["dz"])
        return (s["x"], s["y"], s["z"] + dz / 2.0, s["dx"], s["dy"], dz)
    if s["tipo"] == "corridoio":
        y_lo, y_hi = min(s["y0"], s["y1"]), max(s["y0"], s["y1"])
        z_lo, z_hi = min(s["z0"], s["z1"]), max(s["z0"], s["z1"]) + s["h"]
        return (s["x"], (y_lo + y_hi) / 2.0, (z_lo + z_hi) / 2.0,
                s["w"], y_hi - y_lo, z_hi - z_lo)
    raise ValueError(f"tipo di struttura sconosciuto: {s['tipo']}")


def riferimenti():
    """Strutture interne note nel formato atteso da piramide_unificato.py /
    piramide_acustica_vh.py:

        (num, nome, s_m, h_m, e_m, ns_m, ew_m, alt_m)

    con s_m = posizione NORD-SUD del centro dell'ingombro, h_m = quota del
    CENTRO sopra il piano di base, e_m = offset EST del centro, poi le
    dimensioni N-S x E-O x altezza. Gli offset sono riferiti al NORD VERO
    (la griglia SAR e' ruotata: chi li usa deve ruotarli, non sommarli)."""
    out = []
    for s in STRUTTURE:
        x_c, y_c, z_c, dx, dy, dz = _bbox(s)
        out.append((s["num"], s["nome"], y_c, z_c, x_c, dy, dx, dz))
    return out


def descrizione():
    """Righe di log leggibili con provenienza e grado di certezza."""
    return [f"{s['num']:2d}. {s['nome']} - {s['certezza']}" for s in STRUTTURE]


# ================================================================
# 7. CONTROLLI GEOMETRICI (la catena si chiude da sola)
# ================================================================
def controlli():
    """Verifiche di chiusura: i numeri non sono aggiustati a occhio, escono
    dal calcolo e si incastrano fra loro."""
    def _len_ang(y0, z0, y1, z1):
        return (math.hypot(y1 - y0, z1 - z0),
                math.degrees(math.atan2(abs(z1 - z0), abs(y1 - y0))))
    d_len, d_ang = _len_ang(DISC_START_Y, DISC_START_Z, DISC_END_Y, DISC_END_Z)
    a_len, a_ang = _len_ang(ASC_START_Y, ASC_START_Z, ASC_END_Y, ASC_END_Z)
    g_len, g_ang = _len_ang(GALL_START_Y, GALL_START_Z, GALL_END_Y, GALL_END_Z)
    return [
        ("Corridoio discendente", 105.20, d_len, 26 + 31/60 + 23/3600, d_ang),
        ("Corridoio ascendente", 39.28, a_len, 26 + 2/60 + 30/3600, a_ang),
        ("Grande Galleria (sviluppo)", 46.68, g_len, 26 + 2/60 + 30/3600, g_ang),
    ]


# ================================================================
# 8. DISEGNO (solo se eseguito direttamente)
# ================================================================
def _disegna():
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title("Piramide di Cheope - ricostruzione 3D")

    v = np.array([[-L, -L, 0.0], [L, -L, 0.0], [L, L, 0.0], [-L, L, 0.0],
                  [0.0, 0.0, H]])
    facce = [[v[0], v[1], v[4]], [v[1], v[2], v[4]], [v[2], v[3], v[4]],
             [v[3], v[0], v[4]], [v[0], v[1], v[2], v[3]]]
    ax.add_collection3d(Poly3DCollection(facce, facecolors="gold", linewidths=0.5,
                                         edgecolors="goldenrod", alpha=0.08))

    def box(x_c, y_c, z_b, dim, colore, etichetta):
        dx, dy, dz = dim
        x = np.array([x_c - dx/2, x_c + dx/2, x_c + dx/2, x_c - dx/2] * 2)
        y = np.array([y_c - dy/2, y_c - dy/2, y_c + dy/2, y_c + dy/2] * 2)
        z = np.array([z_b] * 4 + [z_b + dz] * 4)
        idx = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
               [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]]
        verts = [[[x[i], y[i], z[i]] for i in f] for f in idx]
        ax.add_collection3d(Poly3DCollection(verts, facecolors=colore, linewidths=0.5,
                                             edgecolors="black", alpha=0.72,
                                             label=etichetta))

    def corridoio(x_f, y0, z0, y1, z1, width, height, colore, etichetta):
        w = width / 2.0
        a = np.array([[x_f-w, y0, z0], [x_f+w, y0, z0],
                      [x_f+w, y0, z0+height], [x_f-w, y0, z0+height]])
        b = np.array([[x_f-w, y1, z1], [x_f+w, y1, z1],
                      [x_f+w, y1, z1+height], [x_f-w, y1, z1+height]])
        fs = [[a[0], a[1], a[2], a[3]], [b[0], b[1], b[2], b[3]],
              [a[0], a[1], b[1], b[0]], [a[3], a[2], b[2], b[3]],
              [a[0], a[3], b[3], b[0]], [a[1], a[2], b[2], b[1]]]
        ax.add_collection3d(Poly3DCollection(fs, facecolors=colore, linewidths=0.3,
                                             edgecolors="black", alpha=0.68,
                                             label=etichetta))

    for s in STRUTTURE:
        if s["tipo"] == "corridoio":
            corridoio(s["x"], s["y0"], s["z0"], s["y1"], s["z1"],
                      s["w"], s["h"], s["colore"], s["nome"])
        else:
            box(s["x"], s["y"], s["z"], (s["dx"], s["dy"], s["dz"]),
                s["colore"], s["nome"])

    ax.scatter([GREAT_STEP[0]], [GREAT_STEP[1]], [GREAT_STEP[2]],
               s=22, label="Great Step")

    ax.set_xlabel("Asse Est - Ovest (X) [m]")
    ax.set_ylabel("Asse Nord - Sud (Y) [m]")
    ax.set_zlabel("Quota Altimetrica (Z) [m]")
    ax.set_xlim([-120, 120]); ax.set_ylim([-120, 120]); ax.set_zlim([-40, 150])
    ax.view_init(elev=22, azim=-45)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="upper left", fontsize="small")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    print("\n=== PIRAMIDE DI CHEOPE - CONTROLLO GEOMETRICO ===")
    print(f"Base: {BASE:.2f} x {BASE:.2f} m")
    print(f"Altezza originaria: {H:.2f} m")
    print(f"Angolo faccia derivato: {FACE_ANGLE:.4f} gradi")
    print(f"Offset interno E: {X_OFFSET:.2f} m")
    for nome, len_dich, len_calc, ang_dich, ang_calc in controlli():
        print(f"\n--- {nome} ---")
        print(f"Lunghezza dichiarata: {len_dich:.2f} m | da coordinate: {len_calc:.4f} m")
        print(f"Angolo dichiarato: {ang_dich:.4f} gradi | da coordinate: {ang_calc:.4f} gradi")
    print("\n=== STRUTTURE INTERNE (ingombri esportati a piramide_acustica_vh.py) ===")
    for (num, nome, s_m, h_m, e_m, ns_m, ew_m, alt_m) in riferimenti():
        print(f"{num:2d}. {nome:52s} N-S {s_m:+8.2f} m | quota {h_m:+7.2f} m | "
              f"E {e_m:+6.2f} m | ingombro {ns_m:6.2f} x {ew_m:6.2f} x {alt_m:6.2f} m")
    _disegna()

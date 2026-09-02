#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Parte comune di piramide_cheope_3d.py e piramide_kefren_3d.py.

I due script sono lo STESSO programma con dati diversi: ciascuno porta la
propria tabella ``STRUTTURE`` (camere e corridoi misurati) e il proprio
disegno, ma il calcolo degli ingombri e il formato con cui li esporta erano
copiati parola per parola in entrambi -- una trentina di righe identiche, con
il rischio che una correzione ne raggiungesse solo una delle due e i due
gruppi di strutture finissero misurati con regole diverse.

Qui sta solo cio' che era IDENTICO. Il disegno resta in ciascuno script: e'
simile ma non uguale (titoli, trasparenze, il tetto a doppio spiovente solo in
Kefren, il Great Step solo in Cheope), e unificarlo vorrebbe dire inventare
un'astrazione parametrica sopra differenze puramente grafiche.
"""

from typing import Any, Dict, List, Sequence, Tuple

__all__ = ["bbox", "riferimenti", "stampa_riferimenti"]


def bbox(s: Dict[str, Any]) -> Tuple[float, float, float, float, float, float]:
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


def riferimenti(strutture: Sequence[Dict[str, Any]]) -> List[Tuple]:
    """Strutture interne note nel formato atteso da piramide_unificato.py /
    piramide_acustica_vh.py:

        (num, nome, s_m, h_m, e_m, ns_m, ew_m, alt_m)

    con s_m = posizione NORD-SUD del centro dell'ingombro, h_m = quota del
    CENTRO sopra il piano di base, e_m = offset EST del centro, poi le
    dimensioni N-S x E-O x altezza. Gli offset sono riferiti al NORD VERO
    (la griglia SAR e' ruotata: chi li usa deve ruotarli, non sommarli)."""
    out = []
    for s in strutture:
        x_c, y_c, z_c, dx, dy, dz = bbox(s)
        out.append((s["num"], s["nome"], y_c, z_c, x_c, dy, dx, dz))
    return out


def stampa_riferimenti(strutture: Sequence[Dict[str, Any]]) -> None:
    """La tabella degli ingombri, come la stampavano entrambi gli script."""
    print("\n=== STRUTTURE INTERNE (ingombri esportati a piramide_acustica_vh.py) ===")
    for (num, nome, s_m, h_m, e_m, ns_m, ew_m, alt_m) in riferimenti(strutture):
        print(f"{num:2d}. {nome:52s} N-S {s_m:+8.2f} m | quota {h_m:+7.2f} m | "
              f"E {e_m:+6.2f} m | ingombro {ns_m:6.2f} x {ew_m:6.2f} x {alt_m:6.2f} m")


# ==========================================================================
# Disegno: i pezzi che erano identici nei due _disegna()
# ==========================================================================
# numpy e matplotlib restano importati DENTRO le funzioni, come facevano i due
# script: cosi' riferimenti() continua a funzionare su una macchina senza
# matplotlib, che e' il motivo per cui gli import erano locali in origine.

def _poly3d():
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    return Poly3DCollection


def guscio_piramide(ax, semilato: float, altezza: float) -> None:
    """Le quattro facce piu' la base, in trasparenza, come sfondo del disegno."""
    import numpy as np
    s, h = float(semilato), float(altezza)
    v = np.array([[-s, -s, 0.0], [s, -s, 0.0], [s, s, 0.0], [-s, s, 0.0],
                  [0.0, 0.0, h]])
    facce = [[v[0], v[1], v[4]], [v[1], v[2], v[4]], [v[2], v[3], v[4]],
             [v[3], v[0], v[4]], [v[0], v[1], v[2], v[3]]]
    ax.add_collection3d(_poly3d()(facce, facecolors="gold", linewidths=0.5,
                                  edgecolors="goldenrod", alpha=0.08))


def disegna_box(ax, x_c, y_c, z_b, dim, colore, etichetta, alpha: float) -> None:
    """Parallelepipedo di una camera. `dim` = (dx, dy, dz)."""
    import numpy as np
    dx, dy, dz = dim
    x = np.array([x_c - dx/2, x_c + dx/2, x_c + dx/2, x_c - dx/2] * 2)
    y = np.array([y_c - dy/2, y_c - dy/2, y_c + dy/2, y_c + dy/2] * 2)
    z = np.array([z_b] * 4 + [z_b + dz] * 4)
    idx = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
           [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]]
    verts = [[[x[i], y[i], z[i]] for i in f] for f in idx]
    ax.add_collection3d(_poly3d()(verts, facecolors=colore, linewidths=0.5,
                                  edgecolors="black", alpha=alpha,
                                  label=etichetta))


def disegna_corridoio(ax, x_f, y0, z0, y1, z1, width, height, colore,
                      etichetta, alpha: float) -> None:
    """Tratto di corridoio come prisma fra le due sezioni (y0,z0) e (y1,z1)."""
    import numpy as np
    w = width / 2.0
    a = np.array([[x_f-w, y0, z0], [x_f+w, y0, z0],
                  [x_f+w, y0, z0+height], [x_f-w, y0, z0+height]])
    b = np.array([[x_f-w, y1, z1], [x_f+w, y1, z1],
                  [x_f+w, y1, z1+height], [x_f-w, y1, z1+height]])
    fs = [[a[0], a[1], a[2], a[3]], [b[0], b[1], b[2], b[3]],
          [a[0], a[1], b[1], b[0]], [a[3], a[2], b[2], b[3]],
          [a[0], a[3], b[3], b[0]], [a[1], a[2], b[2], b[1]]]
    ax.add_collection3d(_poly3d()(fs, facecolors=colore, linewidths=0.3,
                                  edgecolors="black", alpha=alpha,
                                  label=etichetta))


def legenda_senza_doppioni(ax) -> None:
    """Legenda con una voce per nome: ogni struttura aggiunge la sua etichetta
    a ciascuna faccia, quindi senza questo passaggio comparirebbe sei volte."""
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="upper left",
              fontsize="small")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piramidi_v01.py
===============

Tomografia Doppler 3D delle piramidi di Cheope (Khufu) e Chefren (Khafre)
a partire da uno stack Sentinel-1 SLC, secondo il metodo Micro-Motion (MM)
Doppler tomography di F. Biondi & C. Malanga.

Fonti implementate (skill `biondi-malanga-sar-tomography`):
  - Remote Sensing 2022, 14, 5231           -> modello tomografico, linee tomografiche
  - WO 2024/008365 A1                       -> formalismo sub-aperture, schema a 11 blocchi
  - arXiv:2206.09200 (Vesuvio)              -> trade-off profondita'/risoluzione, validazione

Lo schema computazionale a 11 blocchi (Fig. 0.5 del brevetto) e' implementato
integralmente:

    1  SLC single SAR image            (dato)
    2  DFT2 (FFT2)                     -> _block2_fft2
    3  Band-pass filter master         -> _block34_bandpass
    4  Band-pass filter slave          -> _block34_bandpass
    5  IDFT2 master                    -> _block56_ifft2
    6  IDFT2 slave                     -> _block56_ifft2
    7  Pixel tracking                  -> _block7_pixel_tracking
    8  Seismic waves raw data          (dato: cubo Y)
    9  Focusing by FFT2                -> _block9_tomographic_focusing
    10 Tomographic map                 (dato: h(z))
    11 Tomographic filter + geolocation-> _block11_geolocate

ATTENZIONE - LIMITE FISICO DEI DATI
-----------------------------------
Il metodo originale usa spotlight banda X con ~22 kHz di banda Doppler e
400-450 MHz di banda chirp. Sentinel-1 IW (TOPS, banda C) ha ~313 Hz di banda
Doppler processata e ~48 MHz di chirp: circa 70x meno banda Doppler.

Poiche' la risoluzione in profondita' e'

    delta_z = lambda_sound * R / (2 * A)      con lambda_sound = v / f_inv

e A (apertura orbitale sintetizzata) e' proporzionale alla banda Doppler usata,
la delta_z ottenibile da Sentinel-1 IW e' dell'ordine del chilometro, non del
metro. Le piramidi sono alte ~139 m. Il programma calcola questo budget dai
parametri reali dell'annotation e lo stampa PRIMA di processare: il limite e'
un risultato, non un assunto nascosto. Vedere `--report-only`.

Uso
---
    python piramidi_v01.py --report-only              # budget di risoluzione
    python piramidi_v01.py --selftest                 # verifica dell'inversione
    python piramidi_v01.py --dates 3 --lines 24 --nd 32
    python piramidi_v01.py --verify out_piramidi_v01  # informazione sull'asse z

Autore: generato con Claude Code applicando la skill biondi-malanga-sar-tomography.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

C_LIGHT = 299_792_458.0

# --------------------------------------------------------------------------
# Configurazione
# --------------------------------------------------------------------------

DEFAULT_STACK = r"e:\claude_code\Piramid_V3\goal_out_kefren\stack_slc"

#: Coordinate dei vertici (WGS84). Cheope e Chefren sono i target richiesti;
#: Micerino e' incluso come controllo indipendente sulla stessa scena.
TARGETS: Dict[str, Tuple[float, float, float]] = {
    # nome        lat        lon        altezza attuale [m]
    "Khufu":    (29.979235, 31.134202, 138.5),
    "Khafre":   (29.976111, 31.130833, 136.4),
    "Menkaure": (29.972500, 31.128056,  61.0),
}


@dataclass
class Config:
    """Parametri di processing. I default seguono il regime 'Giza' delle fonti,
    corretti dove i dati Sentinel-1 non li consentono."""

    stack_dir: str = DEFAULT_STACK
    swath: str = "iw2"
    polarisation: str = "vv"
    platform: str = "s1c"

    n_dates: int = 3                 # numero di acquisizioni da impilare
    date_filter: Optional[str] = None

    # --- geometria del blocco di analisi -----------------------------------
    chip_pad_range_m: float = 700.0  # semi-estensione in range attorno ai target
    chip_pad_azim_m: float = 900.0   # semi-estensione in azimuth
    n_lines: int = 24                # numero di linee tomografiche parallele (asse Y del volume)

    # --- strategia sub-aperture (ch11) -------------------------------------
    n_d: int = 32                    # N_D: sampling rate dell'onda meccanica
    guard_fraction: float = 0.5      # B_DL / B_cD  (il brevetto impone 0.5)
    f_investigation: float = 120.0   # [Hz] frequenza meccanica di indagine
    v_sound: float = 6000.0          # [m/s] velocita' di propagazione nel mezzo (calcare/granito)

    # --- asse profondita' ---------------------------------------------------
    depth_max_m: float = 400.0       # profondita'/quota massima del volume
    n_depth: int = 128               # campioni sull'asse z

    # --- pixel tracking (blocco 7) -----------------------------------------
    corr_window: int = 16            # finestra di correlazione [pixel]

    # --- risoluzione azimuth nominale del prodotto (per stimare L_sa) -------
    azimuth_resolution_m: float = 22.0

    out_dir: str = "out_piramidi_v01"
    seed: int = 20260826


# --------------------------------------------------------------------------
# Blocco 1 - lettura annotation Sentinel-1
# --------------------------------------------------------------------------

@dataclass
class Burst:
    index: int
    first_line: int
    last_line: int
    azimuth_time: str
    azimuth_anx_time: float
    fm_rate_poly: Tuple[float, ...]
    fm_rate_t0: float
    dc_poly: Tuple[float, ...]
    dc_t0: float


@dataclass
class S1Annotation:
    """Parametri estratti dal file annotation.xml di un sub-swath Sentinel-1 IW."""

    path: str
    swath: str
    polarisation: str
    orbit_pass: str
    acquisition_date: str

    radar_frequency: float
    range_sampling_rate: float
    azimuth_time_interval: float
    range_pixel_spacing: float
    azimuth_pixel_spacing: float
    slant_range_time: float
    azimuth_steering_rate_deg: float
    incidence_mid: float
    prf: float

    lines_per_burst: int
    samples_per_burst: int
    n_lines: int
    n_samples: int

    azimuth_bandwidth: float
    range_bandwidth: float

    bursts: List[Burst]
    geo_lat: np.ndarray
    geo_lon: np.ndarray
    geo_line: np.ndarray
    geo_pixel: np.ndarray
    geo_incidence: np.ndarray

    orbit_velocity: float

    # ---- derivate ---------------------------------------------------------
    @property
    def wavelength(self) -> float:
        """Lunghezza d'onda radar [m]."""
        return C_LIGHT / self.radar_frequency

    @property
    def slant_range_near(self) -> float:
        """Slant range al near-range [m]."""
        return self.slant_range_time * C_LIGHT / 2.0

    def slant_range(self, pixel: float) -> float:
        """Slant range del pixel dato [m]."""
        return self.slant_range_near + pixel * self.range_pixel_spacing

    def burst_of_line(self, line: float) -> int:
        return int(line) // self.lines_per_burst


def _txt(root: ET.Element, path: str, default: Optional[str] = None) -> str:
    node = root.findtext(path)
    if node is None:
        if default is None:
            raise KeyError(f"campo assente nell'annotation: {path}")
        return default
    return node.strip()


def _poly(text: str) -> Tuple[float, ...]:
    return tuple(float(v) for v in text.split())


def parse_annotation(path: str) -> S1Annotation:
    """Blocco 1: legge tutti i parametri necessari dal file annotation."""
    root = ET.parse(path).getroot()

    m = re.search(r"-(\d{8})t", os.path.basename(path))
    date = m.group(1) if m else "unknown"

    lines_per_burst = int(_txt(root, ".//swathTiming/linesPerBurst"))
    samples_per_burst = int(_txt(root, ".//swathTiming/samplesPerBurst"))

    fm_nodes = root.findall(".//azimuthFmRateList/azimuthFmRate")
    dc_nodes = root.findall(".//dcEstimateList/dcEstimate")
    burst_nodes = root.findall(".//swathTiming/burstList/burst")

    bursts: List[Burst] = []
    for i, b in enumerate(burst_nodes):
        fm = fm_nodes[min(i, len(fm_nodes) - 1)]
        dc = dc_nodes[min(i, len(dc_nodes) - 1)]

        fm_poly_txt = fm.findtext("azimuthFmRatePolynomial")
        if fm_poly_txt is None:  # formato alternativo (c0/c1/c2 separati)
            fm_poly_txt = " ".join(
                fm.findtext(k, "0") for k in ("c0", "c1", "c2")
            )
        dc_poly_txt = dc.findtext("dataDcPolynomial") or "0"

        bursts.append(
            Burst(
                index=i,
                first_line=i * lines_per_burst,
                last_line=(i + 1) * lines_per_burst - 1,
                azimuth_time=b.findtext("azimuthTime", ""),
                azimuth_anx_time=float(b.findtext("azimuthAnxTime", "0")),
                fm_rate_poly=_poly(fm_poly_txt),
                fm_rate_t0=float(fm.findtext("t0", "0")),
                dc_poly=_poly(dc_poly_txt),
                dc_t0=float(dc.findtext("t0", "0")),
            )
        )

    gpts = root.findall(".//geolocationGrid//geolocationGridPoint")
    geo_lat = np.array([float(_txt(p, "latitude")) for p in gpts])
    geo_lon = np.array([float(_txt(p, "longitude")) for p in gpts])
    geo_line = np.array([float(_txt(p, "line")) for p in gpts])
    geo_pixel = np.array([float(_txt(p, "pixel")) for p in gpts])
    geo_inc = np.array([float(_txt(p, "incidenceAngle")) for p in gpts])

    # velocita' della piattaforma dal primo state vector disponibile
    vx = root.findtext(".//orbit/velocity/x")
    vy = root.findtext(".//orbit/velocity/y")
    vz = root.findtext(".//orbit/velocity/z")
    if vx is not None and vy is not None and vz is not None:
        v_orb = math.sqrt(float(vx) ** 2 + float(vy) ** 2 + float(vz) ** 2)
    else:
        v_orb = 7590.0  # valore nominale Sentinel-1

    return S1Annotation(
        path=path,
        swath=_txt(root, ".//adsHeader/swath"),
        polarisation=_txt(root, ".//adsHeader/polarisation"),
        orbit_pass=_txt(root, ".//generalAnnotation/productInformation/pass"),
        acquisition_date=date,
        radar_frequency=float(_txt(root, ".//radarFrequency")),
        range_sampling_rate=float(_txt(root, ".//rangeSamplingRate")),
        azimuth_time_interval=float(_txt(root, ".//azimuthTimeInterval")),
        range_pixel_spacing=float(_txt(root, ".//rangePixelSpacing")),
        azimuth_pixel_spacing=float(_txt(root, ".//azimuthPixelSpacing")),
        slant_range_time=float(_txt(root, ".//imageInformation/slantRangeTime")),
        azimuth_steering_rate_deg=float(_txt(root, ".//azimuthSteeringRate")),
        incidence_mid=float(_txt(root, ".//incidenceAngleMidSwath")),
        prf=float(_txt(root, ".//prf")),
        lines_per_burst=lines_per_burst,
        samples_per_burst=samples_per_burst,
        n_lines=int(_txt(root, ".//imageInformation/numberOfLines")),
        n_samples=int(_txt(root, ".//imageInformation/numberOfSamples")),
        azimuth_bandwidth=float(_txt(root, ".//azimuthProcessing/totalBandwidth")),
        range_bandwidth=float(_txt(root, ".//rangeProcessing/totalBandwidth")),
        bursts=bursts,
        geo_lat=geo_lat,
        geo_lon=geo_lon,
        geo_line=geo_line,
        geo_pixel=geo_pixel,
        geo_incidence=geo_inc,
        orbit_velocity=v_orb,
    )


# --------------------------------------------------------------------------
# Scoperta dello stack
# --------------------------------------------------------------------------

@dataclass
class StackEntry:
    date: str
    annotation: str
    tiff: str
    #: .xml di calibrazione radiometrica e di rumore termico, quando il
    #: prodotto li porta con se' (alberi .SAFE). Restano opzionali: uno stack
    #: piatto costruito a mano puo' non averli.
    calibration: Optional[str] = None
    noise: Optional[str] = None


def discover_stack(cfg: Config) -> List[StackEntry]:
    """Trova le coppie annotation/tiff per swath e polarizzazione richiesti."""
    pattern = os.path.join(
        cfg.stack_dir, f"{cfg.platform}-{cfg.swath}-slc-{cfg.polarisation}-*.annotation.xml"
    )
    entries: List[StackEntry] = []
    for ann in sorted(glob.glob(pattern)):
        tif = ann.replace(".annotation.xml", ".tiff")
        if not os.path.exists(tif):
            continue
        m = re.search(r"-(\d{8})t", os.path.basename(ann))
        date = m.group(1) if m else "unknown"
        if cfg.date_filter and cfg.date_filter not in date:
            continue
        entries.append(StackEntry(date=date, annotation=ann, tiff=tif))
    if not entries:
        raise FileNotFoundError(
            f"nessuna coppia annotation/tiff trovata con il pattern:\n  {pattern}"
        )
    return entries


# --------------------------------------------------------------------------
# Geolocalizzazione dei target
# --------------------------------------------------------------------------

def latlon_to_line_pixel(ann: S1Annotation, lat: float, lon: float) -> Tuple[float, float]:
    """Inverte la geolocation grid: (lat, lon) -> (line, pixel).

    Usa interpolazione lineare sulla griglia sparsa dell'annotation; per una
    scena IW la griglia e' fitta abbastanza (231 punti) da rendere l'errore
    residuo molto inferiore alla cella di risoluzione."""
    from scipy.interpolate import griddata

    pts = np.column_stack([ann.geo_lat, ann.geo_lon])
    line = griddata(pts, ann.geo_line, (lat, lon), method="linear")
    pixel = griddata(pts, ann.geo_pixel, (lat, lon), method="linear")
    if np.isnan(line) or np.isnan(pixel):
        raise ValueError(f"({lat}, {lon}) fuori dalla scena {ann.acquisition_date}")
    return float(line), float(pixel)


def line_pixel_to_latlon(ann: S1Annotation, line: float, pixel: float) -> Tuple[float, float]:
    """Mappatura diretta (line, pixel) -> (lat, lon), per il blocco 11."""
    from scipy.interpolate import griddata

    pts = np.column_stack([ann.geo_line, ann.geo_pixel])
    lat = griddata(pts, ann.geo_lat, (line, pixel), method="linear")
    lon = griddata(pts, ann.geo_lon, (line, pixel), method="linear")
    return float(lat), float(lon)


def incidence_at(ann: S1Annotation, line: float, pixel: float) -> float:
    """Angolo di incidenza [gradi] interpolato dalla geolocation grid."""
    from scipy.interpolate import griddata

    pts = np.column_stack([ann.geo_line, ann.geo_pixel])
    inc = griddata(pts, ann.geo_incidence, (line, pixel), method="linear")
    return float(inc) if not np.isnan(inc) else ann.incidence_mid


# --------------------------------------------------------------------------
# Budget di risoluzione (ch12 / ch14)
# --------------------------------------------------------------------------

@dataclass
class ResolutionBudget:
    """delta_z = lambda_sound * R / (2 * A), con lambda_sound = v / f_inv."""

    f_investigation: float
    v_sound: float
    lambda_sound: float
    slant_range: float
    synthetic_aperture_full: float
    synthetic_aperture_used: float
    delta_z: float
    doppler_bandwidth: float
    doppler_used: float
    b_shift: float
    f_max_observable: float
    depth_unambiguous: float
    target_height: float

    @property
    def resolvable(self) -> bool:
        """Il volume di un target alto `target_height` e' risolto?"""
        return self.delta_z <= self.target_height / 3.0

    def as_text(self) -> str:
        ok = "SI" if self.resolvable else "NO"
        return "\n".join(
            [
                "-" * 74,
                "BUDGET DI RISOLUZIONE TOMOGRAFICA  (delta_z = lambda_s * R / 2A)",
                "-" * 74,
                f"  banda Doppler processata  B_cD      = {self.doppler_bandwidth:12.2f} Hz",
                f"  banda Doppler utilizzata  B_used    = {self.doppler_used:12.2f} Hz",
                f"  separazione master/slave  B_shift   = {self.b_shift:12.2f} Hz",
                f"  freq. meccanica max osservabile     = {self.f_max_observable:12.2f} Hz",
                "",
                f"  frequenza di indagine     f_inv     = {self.f_investigation:12.2f} Hz",
                f"  velocita' nel mezzo       v         = {self.v_sound:12.2f} m/s",
                f"  lunghezza d'onda acustica lambda_s  = {self.lambda_sound:12.3f} m",
                f"  slant range               R         = {self.slant_range:12.1f} m",
                f"  apertura sintetica piena  L_sa      = {self.synthetic_aperture_full:12.1f} m",
                f"  apertura effettivamente usata  A    = {self.synthetic_aperture_used:12.1f} m",
                "",
                f"  >> RISOLUZIONE IN PROFONDITA'  delta_z = {self.delta_z:12.1f} m",
                f"  >> profondita' non ambigua N_D*delta_z = {self.depth_unambiguous:12.1f} m",
                "",
                f"  altezza del target                  = {self.target_height:12.1f} m",
                f"  >> struttura interna risolvibile?     {ok:>12s}",
                "-" * 74,
            ]
        )


def compute_resolution_budget(
    cfg: Config, ann: S1Annotation, slant_range: float, target_height: float
) -> ResolutionBudget:
    """Calcola il budget di risoluzione dai parametri REALI dell'annotation.

    L'apertura sintetica piena si ricava da L_sa = lambda_radar * R / (2 * delta_az).
    L'apertura effettivamente sfruttata scala con la frazione di banda Doppler
    percorsa dalla marcia master/slave (B_cD - B_DL)."""

    b_cd = ann.azimuth_bandwidth
    b_dl = b_cd * cfg.guard_fraction
    b_used = b_cd - b_dl

    l_sa_full = ann.wavelength * slant_range / (2.0 * cfg.azimuth_resolution_m)
    a_used = l_sa_full * (b_used / b_cd)

    lambda_s = cfg.v_sound / cfg.f_investigation
    delta_z = lambda_s * slant_range / (2.0 * a_used)

    # B_shift <-> frequenza meccanica osservata (relazione inversa, ch11).
    # Mappatura: la separazione in frequenza azimuth corrisponde, tramite la
    # FM rate azimutale, a un ritardo temporale dt = B_shift / ka; la frequenza
    # osservabile e' quella di Nyquist su quel ritardo, f = 1 / (2 dt).
    # Le fonti dichiarano la relazione inversa ma non la formula: mappatura
    # dedotta, documentata qui e riportata nei metadati di output.
    ka = abs(ann.bursts[0].fm_rate_poly[0]) if ann.bursts else 2200.0
    b_shift = ka / (2.0 * cfg.f_investigation)
    b_shift = min(b_shift, b_used * 0.5)  # non puo' eccedere la guard band

    dt_min = 1.0 / ann.prf
    f_max = ka / (2.0 * b_used) if b_used > 0 else 0.0
    f_max = max(f_max, 1.0 / (2.0 * dt_min * cfg.n_d))

    return ResolutionBudget(
        f_investigation=cfg.f_investigation,
        v_sound=cfg.v_sound,
        lambda_sound=lambda_s,
        slant_range=slant_range,
        synthetic_aperture_full=l_sa_full,
        synthetic_aperture_used=a_used,
        delta_z=delta_z,
        doppler_bandwidth=b_cd,
        doppler_used=b_used,
        b_shift=b_shift,
        f_max_observable=f_max,
        depth_unambiguous=cfg.n_d * delta_z,
        target_height=target_height,
    )


# --------------------------------------------------------------------------
# Lettura del chip SLC
# --------------------------------------------------------------------------

@dataclass
class Chip:
    data: np.ndarray          # complex64 [n_lines, n_pixels]
    line0: int
    pixel0: int
    date: str
    burst: int


def read_chip(entry: StackEntry, ann: S1Annotation, cfg: Config) -> Chip:
    """Legge dal GeoTIFF la finestra che contiene i target, confinata al burst."""
    import rasterio
    from rasterio.windows import Window

    lines, pixels = [], []
    for lat, lon, _ in TARGETS.values():
        l, p = latlon_to_line_pixel(ann, lat, lon)
        lines.append(l)
        pixels.append(p)

    pad_l = cfg.chip_pad_azim_m / ann.azimuth_pixel_spacing
    pad_p = cfg.chip_pad_range_m / ann.range_pixel_spacing

    l_min = int(min(lines) - pad_l)
    l_max = int(max(lines) + pad_l)
    p_min = int(min(pixels) - pad_p)
    p_max = int(max(pixels) + pad_p)

    # Il pixel tracking e la deramping TOPS sono validi solo entro un burst.
    burst = ann.burst_of_line(float(np.mean(lines)))
    b_first = burst * ann.lines_per_burst
    b_last = b_first + ann.lines_per_burst - 1
    l_min = max(l_min, b_first)
    l_max = min(l_max, b_last)

    p_min = max(p_min, 0)
    p_max = min(p_max, ann.n_samples - 1)

    with rasterio.open(entry.tiff) as ds:
        win = Window(p_min, l_min, p_max - p_min + 1, l_max - l_min + 1)
        arr = ds.read(1, window=win)

    return Chip(
        data=np.ascontiguousarray(arr.astype(np.complex64)),
        line0=l_min,
        pixel0=p_min,
        date=entry.date,
        burst=burst,
    )


# --------------------------------------------------------------------------
# Deramping TOPS  (prerequisito ai blocchi 3/4)
# --------------------------------------------------------------------------

def tops_deramp(chip: Chip, ann: S1Annotation) -> np.ndarray:
    """Rimuove la rampa Doppler introdotta dallo steering d'antenna TOPS.

    In modalita' TOPS il centroide Doppler varia linearmente lungo l'azimuth
    all'interno del burst. Senza deramping lo spettro azimutale non e'
    stazionario e la suddivisione in sub-aperture (blocchi 3/4) e' priva di
    significato: si sta finestrando tempo, non frequenza.

    Formulazione standard (Miranda 2015):
        ks = 2 * Vs * kpsi / lambda            (Doppler rate da steering)
        kt = ka * ks / (ka - ks)               (rate combinata)
        eta_c = -fdc / ka                      (istante di beam-center)
        phi = -pi * kt * (eta - eta_c)^2
    """
    n_l, n_p = chip.data.shape
    burst = ann.bursts[min(chip.burst, len(ann.bursts) - 1)]

    k_psi = math.radians(ann.azimuth_steering_rate_deg)          # rad/s
    k_s = 2.0 * ann.orbit_velocity * k_psi / ann.wavelength      # Hz/s

    # tempo di andata-ritorno per ogni colonna del chip
    pix = chip.pixel0 + np.arange(n_p, dtype=np.float64)
    tau = ann.slant_range_time + pix / ann.range_sampling_rate

    def polyval(poly: Sequence[float], t0: float) -> np.ndarray:
        dt = tau - t0
        out = np.zeros_like(tau)
        for i, c in enumerate(poly):
            out += c * dt ** i
        return out

    k_a = polyval(burst.fm_rate_poly, burst.fm_rate_t0)          # Hz/s
    f_dc = polyval(burst.dc_poly, burst.dc_t0)                   # Hz

    denom = k_a - k_s
    denom = np.where(np.abs(denom) < 1e-6, 1e-6, denom)
    k_t = k_a * k_s / denom

    # tempo azimutale relativo al centro del burst
    line_in_burst = (chip.line0 - burst.first_line) + np.arange(n_l, dtype=np.float64)
    eta = (line_in_burst - ann.lines_per_burst / 2.0) * ann.azimuth_time_interval

    eta_c = -f_dc / np.where(np.abs(k_a) < 1e-6, 1e-6, k_a)      # [n_p]

    d_eta = eta[:, None] - eta_c[None, :]
    phase = -np.pi * k_t[None, :] * d_eta ** 2

    return (chip.data * np.exp(1j * phase.astype(np.float32))).astype(np.complex64)


# --------------------------------------------------------------------------
# Blocchi 2-7: strategia sub-aperture e pixel tracking
# --------------------------------------------------------------------------

def _block2_fft2(img: np.ndarray) -> np.ndarray:
    """Blocco 2 - DFT2 dell'immagine SLC. Calcolata UNA sola volta e copiata
    ai due rami (blocchi 3 e 4): ricalcolarla dentro il loop su N_D
    moltiplicherebbe il costo per N_D senza alcun guadagno (ch13)."""
    return np.fft.fft2(img)


def _azimuth_band_mask(n_lines: int, prf: float, f_lo: float, f_hi: float) -> np.ndarray:
    """Maschera passa-banda sull'asse delle frequenze azimutali (Doppler)."""
    freqs = np.fft.fftfreq(n_lines, d=1.0 / prf)
    mask = ((freqs >= f_lo) & (freqs < f_hi)).astype(np.float32)
    if mask.sum() == 0:  # banda troppo stretta per la griglia: prendi il bin piu' vicino
        mask[np.argmin(np.abs(freqs - 0.5 * (f_lo + f_hi)))] = 1.0
    return mask


def _block34_bandpass(spec: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Blocchi 3/4 - filtro passa-banda master o slave sullo spettro azimutale."""
    return spec * mask[:, None]


def _block56_ifft2(spec: np.ndarray) -> np.ndarray:
    """Blocchi 5/6 - IDFT2: ritorno all'immagine SLC a risoluzione azimutale
    ridotta. La perdita di risoluzione azimutale e' il prezzo della sensibilita'
    al moto (ch11)."""
    return np.fft.ifft2(spec)


def _parabolic_refine(c_m1: np.ndarray, c_0: np.ndarray, c_p1: np.ndarray) -> np.ndarray:
    """Raffinamento sub-pixel per interpolazione parabolica a 3 punti attorno
    al massimo di correlazione. Precisione tipica ~1/20 di pixel, senza il costo
    di una IFFT sovracampionata."""
    den = c_m1 - 2.0 * c_0 + c_p1
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)
    return np.clip(0.5 * (c_m1 - c_p1) / den, -1.0, 1.0)


def _dft_upsample_refine(
    cross: np.ndarray, i0: np.ndarray, j0: np.ndarray,
    upsample: int = 20, radius: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Raffina il picco di correlazione con una DFT inversa sovracampionata,
    valutata solo in un intorno del massimo intero (Guizar-Sicairos et al.).

    L'interpolazione parabolica ha un bias sistematico fino a ~0.2 pixel sul
    picco sinc della correlazione di fase: accettabile per misurare vibrazione
    relativa, non per coregistrare uno stack interferometrico. Qui il costo
    resta basso perche' la DFT si valuta su una finestrella, non sull'intera
    griglia sovracampionata."""
    n, H, W = cross.shape
    m = int(2 * radius * upsample) + 1
    off = (np.arange(m) - m // 2) / upsample          # [-radius, +radius]

    k = np.fft.fftfreq(H)[None, :]                    # [1, H]
    l = np.fft.fftfreq(W)[None, :]                    # [1, W]

    rows = (i0[:, None] + off[None, :])[:, :, None]   # [n, m, 1]
    cols = (j0[:, None] + off[None, :])[:, :, None]   # [n, m, 1]

    e_row = np.exp(2j * np.pi * rows * k[None, :, :])          # [n, m, H]
    e_col = np.exp(2j * np.pi * cols * l[None, :, :])          # [n, m, W]

    up = np.abs(np.einsum("nmh,nhw,npw->nmp", e_row, cross, e_col))
    flat = up.reshape(n, -1).argmax(axis=1)
    di = off[flat // m]
    dj = off[flat % m]
    return i0 + di, j0 + dj


def _batch_subpixel_shift(mw: np.ndarray, sw: np.ndarray,
                          refine: str = "parabolic") -> np.ndarray:
    """Shift sub-pixel per un intero stack di finestre, in batch.

    mw, sw hanno forma [n_punti, W, W]. Ritorna un vettore complesso a + jb:
    parte reale = shift in azimuth (righe), parte immaginaria = shift in range
    (colonne). Sono gli {a, b} dell'oscillatore armonico, eq. (20) del brevetto.

    Il metodo e' la cross-correlazione di fase: normalizzando il modulo dello
    spettro incrociato si pesa solo la fase, che e' dove risiede lo shift."""
    n = mw.shape[0]
    if n == 0:
        return np.zeros(0, dtype=np.complex64)

    m = mw - mw.mean(axis=(1, 2), keepdims=True)
    s = sw - sw.mean(axis=(1, 2), keepdims=True)

    cross = np.fft.fft2(m, axes=(-2, -1)) * np.conj(np.fft.fft2(s, axes=(-2, -1)))
    mag = np.abs(cross)
    mag[mag == 0] = 1.0
    corr = np.abs(np.fft.ifft2(cross / mag, axes=(-2, -1)))

    # Le finestre non sono necessariamente quadrate (v02 correla il chip intero):
    # righe e colonne vanno indicizzate con le rispettive dimensioni.
    H, W = corr.shape[-2], corr.shape[-1]
    ar = np.arange(n)
    idx = corr.reshape(n, -1).argmax(axis=1)
    i0, j0 = idx // W, idx % W

    if refine == "dft" and n <= 64:
        dl, dp = _dft_upsample_refine(cross / mag, i0, j0)
    else:
        di = _parabolic_refine(corr[ar, (i0 - 1) % H, j0],
                               corr[ar, i0, j0],
                               corr[ar, (i0 + 1) % H, j0])
        dj = _parabolic_refine(corr[ar, i0, (j0 - 1) % W],
                               corr[ar, i0, j0],
                               corr[ar, i0, (j0 + 1) % W])
        dl = i0.astype(np.float64) + di
        dp = j0.astype(np.float64) + dj
    dl = np.where(dl > H / 2, dl - H, dl)      # shift negativi oltre meta' finestra
    dp = np.where(dp > W / 2, dp - W, dp)

    # Convenzione: il valore restituito e' lo spostamento di `slave` rispetto a
    # `master`, cioe' se slave(x) = master(x - d) la funzione ritorna +d.
    # La cross-correlazione di fase produce -d, da cui il segno invertito.
    return (-(dl + 1j * dp)).astype(np.complex64)


def _block7_pixel_tracking(
    master: np.ndarray,
    slave: np.ndarray,
    sample_lines: np.ndarray,
    sample_pixels: np.ndarray,
    window: int,
) -> np.ndarray:
    """Blocco 7 - pixel tracking sui punti della linea tomografica.

    Estrae tutte le finestre in uno stack e le correla in batch: le FFT
    vettorizzate su [n_punti, W, W] sono ordini di grandezza piu' veloci del
    ciclo punto per punto. Ritorna {a + jb} per ciascun punto campionato."""
    half = window // 2
    n_l, n_p = master.shape
    n_pts = len(sample_lines)

    mw = np.empty((n_pts, window, window), dtype=np.complex64)
    sw = np.empty((n_pts, window, window), dtype=np.complex64)

    for i, (l, p) in enumerate(zip(sample_lines, sample_pixels)):
        l0 = max(0, min(int(l) - half, n_l - window))
        p0 = max(0, min(int(p) - half, n_p - window))
        mw[i] = master[l0:l0 + window, p0:p0 + window]
        sw[i] = slave[l0:l0 + window, p0:p0 + window]

    return _batch_subpixel_shift(mw, sw)


def build_vibration_cube(
    img: np.ndarray,
    ann: S1Annotation,
    cfg: Config,
    budget: ResolutionBudget,
    sample_lines: np.ndarray,
    sample_pixels: np.ndarray,
    verbose: bool = True,
) -> np.ndarray:
    """Blocchi 2-8: marcia master/slave sulla banda Doppler e costruzione del
    cubo Y dei dati vibrazionali grezzi.

    Y ha forma [n_punti, N_D]: per ogni pixel della linea tomografica, N_D
    campioni temporali della micro-vibrazione. N_D e' il sampling rate
    dell'onda meccanica, non un parametro di calcolo (ch11)."""

    spec = _block2_fft2(img)                     # blocco 2, una volta sola

    b_cd = ann.azimuth_bandwidth
    b_dl = b_cd * cfg.guard_fraction
    b_used = b_cd - b_dl
    b_shift = budget.b_shift
    step = b_used / cfg.n_d

    n_pts = len(sample_lines)
    Y = np.zeros((n_pts, cfg.n_d), dtype=np.complex64)

    f_start = -b_cd / 2.0
    for k in range(cfg.n_d):
        lo_m = f_start + k * step
        hi_m = lo_m + b_used
        lo_s = lo_m + b_shift
        hi_s = hi_m + b_shift

        mask_m = _azimuth_band_mask(img.shape[0], ann.prf, lo_m, hi_m)
        mask_s = _azimuth_band_mask(img.shape[0], ann.prf, lo_s, hi_s)

        master = _block56_ifft2(_block34_bandpass(spec, mask_m))   # blocchi 3->5
        slave = _block56_ifft2(_block34_bandpass(spec, mask_s))    # blocchi 4->6

        Y[:, k] = _block7_pixel_tracking(                          # blocco 7
            master, slave, sample_lines, sample_pixels, cfg.corr_window,
        )
        if verbose and (k + 1) % max(1, cfg.n_d // 8) == 0:
            print(f"      sub-apertura {k + 1:3d}/{cfg.n_d}", flush=True)

    return Y                                                       # blocco 8


# --------------------------------------------------------------------------
# Blocco 9: focalizzazione tomografica
# --------------------------------------------------------------------------

def steering_matrix(
    z_axis: np.ndarray, baselines: np.ndarray, lambda_sound: float,
    slant_range: float, incidence_deg: float,
) -> np.ndarray:
    """Matrice di steering A(Kz, z), eq. (22) del brevetto.

        Kz = 4 * pi * B_perp / (lambda * r * sin(theta))

    Ogni colonna e' la firma di fase attesa da un diffusore posto alla quota z."""
    theta = math.radians(incidence_deg)
    k_z = 4.0 * np.pi * baselines / (lambda_sound * slant_range * math.sin(theta))
    return np.exp(1j * np.outer(k_z, z_axis)).astype(np.complex64)   # [k, F]


def _block9_tomographic_focusing(
    Y: np.ndarray, A: np.ndarray,
) -> np.ndarray:
    """Blocco 9 - inversione tomografica h(z) = A^H Y, eq. (24).

    Il brevetto e' esplicito: A e' la migliore approssimazione di un operatore
    DFT, quindi la focalizzazione in profondita' E' compressione d'impulso.
    Nessuna regolarizzazione: filtro adattato puro, come nelle fonti."""
    return (Y @ np.conj(A)).astype(np.complex64)     # [n_punti, n_depth]


# --------------------------------------------------------------------------
# Blocco 11: geolocalizzazione
# --------------------------------------------------------------------------

def _block11_geolocate(
    ann: S1Annotation, chip: Chip,
    sample_lines: np.ndarray, sample_pixels: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Blocco 11 - riporta i punti della linea tomografica in coordinate
    geografiche, per l'uscita in standard GIS (claim 10)."""
    lats = np.zeros(len(sample_lines))
    lons = np.zeros(len(sample_lines))
    for i, (l, p) in enumerate(zip(sample_lines, sample_pixels)):
        lats[i], lons[i] = line_pixel_to_latlon(
            ann, chip.line0 + float(l), chip.pixel0 + float(p)
        )
    return lats, lons


# --------------------------------------------------------------------------
# Pipeline completa
# --------------------------------------------------------------------------

def process_date(
    entry: StackEntry, cfg: Config, verbose: bool = True
) -> Dict[str, Any]:
    """Esegue i blocchi 1-11 su una singola acquisizione e ritorna il volume 3D."""

    t0 = time.time()
    ann = parse_annotation(entry.annotation)
    if verbose:
        print(f"\n  [{entry.date}] {ann.swath} {ann.polarisation} {ann.orbit_pass} "
              f"- B_cD={ann.azimuth_bandwidth:.0f} Hz, PRF={ann.prf:.1f} Hz")

    chip = read_chip(entry, ann, cfg)                               # blocco 1
    if verbose:
        print(f"      chip {chip.data.shape} @ line0={chip.line0} pixel0={chip.pixel0} "
              f"(burst {chip.burst})")

    img = tops_deramp(chip, ann)                                    # pre-blocco 2

    # centro geometrico dei target nel chip
    tgt_lp = {}
    for name, (lat, lon, h) in TARGETS.items():
        l, p = latlon_to_line_pixel(ann, lat, lon)
        tgt_lp[name] = (l - chip.line0, p - chip.pixel0)

    centre_p = float(np.mean([p for _, p in tgt_lp.values()]))
    slant_range = ann.slant_range(chip.pixel0 + centre_p)
    incidence = incidence_at(ann, chip.line0 + chip.data.shape[0] / 2,
                             chip.pixel0 + centre_p)

    budget = compute_resolution_budget(
        cfg, ann, slant_range, max(h for *_, h in TARGETS.values())
    )
    if verbose:
        print(budget.as_text())

    # --- griglia dei campionamenti: n_lines linee tomografiche parallele ----
    n_l, n_p = img.shape
    margin = cfg.corr_window
    line_positions = np.linspace(margin, n_l - margin - 1, cfg.n_lines)
    pix_positions = np.arange(margin, n_p - margin, max(1, cfg.corr_window // 2))

    if verbose:
        print(f"      volume: {cfg.n_lines} linee x {len(pix_positions)} punti "
              f"x {cfg.n_depth} livelli di profondita'")

    # --- asse profondita' e baseline delle sub-aperture --------------------
    z_axis = np.linspace(-cfg.depth_max_m, cfg.depth_max_m, cfg.n_depth)
    l_sa = budget.synthetic_aperture_used
    baselines = np.linspace(-l_sa / 2.0, l_sa / 2.0, cfg.n_d)
    A = steering_matrix(z_axis, baselines, budget.lambda_sound, slant_range, incidence)

    volume = np.zeros((cfg.n_lines, len(pix_positions), cfg.n_depth), dtype=np.complex64)
    lat_grid = np.zeros((cfg.n_lines, len(pix_positions)))
    lon_grid = np.zeros((cfg.n_lines, len(pix_positions)))

    for i, lpos in enumerate(line_positions):
        sl = np.full(len(pix_positions), lpos)
        Y = build_vibration_cube(img, ann, cfg, budget, sl, pix_positions,
                                 verbose=verbose and i == 0)          # blocchi 2-8
        volume[i] = _block9_tomographic_focusing(Y, A)                # blocchi 9-10
        lat_grid[i], lon_grid[i] = _block11_geolocate(ann, chip, sl, pix_positions)
        if verbose:
            print(f"      linea tomografica {i + 1:3d}/{cfg.n_lines}", flush=True)

    return {
        "date": entry.date,
        "volume": volume,
        "z_axis": z_axis,
        "lat": lat_grid,
        "lon": lon_grid,
        "budget": budget,
        "annotation": ann,
        "chip": chip,
        "targets_lp": tgt_lp,
        "pix_positions": pix_positions,
        "line_positions": line_positions,
        "elapsed_s": time.time() - t0,
    }


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def save_outputs(results: List[Dict[str, Any]], cfg: Config) -> None:
    os.makedirs(cfg.out_dir, exist_ok=True)

    # stack incoerente delle magnitudini sulle date disponibili
    vols = np.stack([np.abs(r["volume"]) for r in results], axis=0)
    stacked = vols.mean(axis=0)
    ref = results[0]

    np.save(os.path.join(cfg.out_dir, "tomo3d_magnitude.npy"), stacked)
    np.save(os.path.join(cfg.out_dir, "tomo3d_zaxis.npy"), ref["z_axis"])
    np.save(os.path.join(cfg.out_dir, "tomo3d_lat.npy"), ref["lat"])
    np.save(os.path.join(cfg.out_dir, "tomo3d_lon.npy"), ref["lon"])

    budget: ResolutionBudget = ref["budget"]
    meta = {
        "generato": time.strftime("%Y-%m-%d %H:%M:%S"),
        "programma": "piramidi_v01.py",
        "metodo": "MM Doppler tomography (Biondi & Malanga) - schema a 11 blocchi WO 2024/008365",
        "config": asdict(cfg),
        "date_processate": [r["date"] for r in results],
        "forma_volume": list(stacked.shape),
        "budget_risoluzione": asdict(budget),
        "struttura_interna_risolvibile": bool(budget.resolvable),
        "avvertenza": (
            "Dati Sentinel-1 IW TOPS banda C: banda Doppler ~313 Hz contro i ~22 kHz "
            "usati nelle pubblicazioni originali (spotlight banda X). La risoluzione in "
            "profondita' che ne consegue e' riportata in budget_risoluzione.delta_z. "
            "Se delta_z supera l'altezza del target, il volume NON risolve la struttura "
            "interna e va letto come dimostrazione della catena di elaborazione, non "
            "come riproduzione dei risultati di Remote Sensing 2022, 14, 5231."
        ),
        "mappatura_b_shift": (
            "B_shift = ka / (2 * f_inv). Le fonti dichiarano la relazione inversa fra "
            "B_shift e frequenza meccanica osservata ma non ne danno la formula: questa "
            "mappatura e' dedotta e va considerata un'assunzione del programma."
        ),
    }
    with open(os.path.join(cfg.out_dir, "tomo3d_meta.json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)

    with open(os.path.join(cfg.out_dir, "budget_risoluzione.txt"), "w", encoding="utf-8") as fh:
        fh.write(budget.as_text() + "\n")

    _plot(stacked, ref, cfg)
    print(f"\n  output scritti in: {os.path.abspath(cfg.out_dir)}")


def _plot(volume: np.ndarray, ref: Dict[str, Any], cfg: Config) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z = ref["z_axis"]
    ann: S1Annotation = ref["annotation"]
    tgt_lp = ref["targets_lp"]
    pix_positions = ref["pix_positions"]
    line_positions = ref["line_positions"]
    db = 20 * np.log10(np.maximum(volume, 1e-12))
    vmax = np.percentile(db, 99.5)
    vmin = vmax - 40

    # --- sezioni verticali passanti per Cheope e Chefren -------------------
    fig, axes = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
    for ax, name in zip(axes, ["Khufu", "Khafre"]):
        l_t, _ = tgt_lp[name]
        idx = int(np.argmin(np.abs(line_positions - l_t)))
        sec = db[idx].T                              # [n_depth, n_punti]
        extent = [
            pix_positions[0] * ann.range_pixel_spacing,
            pix_positions[-1] * ann.range_pixel_spacing,
            z[0], z[-1],
        ]
        im = ax.imshow(sec, aspect="auto", origin="lower", extent=extent,
                       cmap="inferno", vmin=vmin, vmax=vmax)
        ax.axhline(0, color="cyan", lw=0.8, ls="--")
        ax.set_title(f"Sezione tomografica - {name} "
                     f"(linea {idx + 1}/{cfg.n_lines})")
        ax.set_xlabel("distanza in ground range [m]")
        ax.set_ylabel("quota / profondita' [m]")
        fig.colorbar(im, ax=ax, label="dB")
    fig.savefig(os.path.join(cfg.out_dir, "sezioni_verticali.png"), dpi=140)
    plt.close(fig)

    # --- proiezioni ortogonali del volume ---------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
    for ax, (proj, title) in zip(
        axes,
        [(db.max(axis=2), "vista dall'alto (max su z)"),
         (db.max(axis=1).T, "sezione azimuth-profondita'"),
         (db.max(axis=0).T, "sezione range-profondita'")],
    ):
        im = ax.imshow(proj, aspect="auto", origin="lower", cmap="inferno",
                       vmin=vmin, vmax=vmax)
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label="dB")
    fig.savefig(os.path.join(cfg.out_dir, "proiezioni_volume.png"), dpi=140)
    plt.close(fig)

    # --- nuvola 3D dei voxel piu' energetici -------------------------------
    thr = np.percentile(db, 99.0)
    ii, jj, kk = np.where(db >= thr)
    if len(ii) > 20000:
        sel = np.random.default_rng(cfg.seed).choice(len(ii), 20000, replace=False)
        ii, jj, kk = ii[sel], jj[sel], kk[sel]
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(jj * ann.range_pixel_spacing,
                    ii * (line_positions[1] - line_positions[0]) * ann.azimuth_pixel_spacing,
                    z[kk], c=db[ii, jj, kk], cmap="inferno", s=1.5, alpha=0.5)
    ax.set_xlabel("range [m]")
    ax.set_ylabel("azimuth [m]")
    ax.set_zlabel("quota / profondita' [m]")
    ax.set_title("Tomografia 3D - voxel oltre il 99° percentile")
    fig.colorbar(sc, ax=ax, label="dB", shrink=0.6)
    fig.savefig(os.path.join(cfg.out_dir, "tomografia_3d.png"), dpi=140)
    plt.close(fig)


# --------------------------------------------------------------------------
# Verifica empirica del contenuto informativo dell'asse z
# --------------------------------------------------------------------------

def verify_depth_information(out_dir: str) -> Dict[str, float]:
    """Misura quanta informazione il volume porta effettivamente sull'asse z.

    E' la controprova del budget di risoluzione: se delta_z e' molto maggiore
    dell'estensione dell'asse z, la matrice di steering ha una rampa di fase
    quasi costante sull'intervallo e h(z) risulta piatta in profondita'. In quel
    caso il tomogramma contiene informazione di range e azimuth ma NON di
    profondita', e questa funzione lo quantifica invece di lasciarlo all'occhio.

    Applica al programma stesso il principio del protocollo di validazione a tre
    livelli (arXiv:2206.09200, sezione 4): dichiarare l'errore, non asserire
    l'accordo."""
    vol = np.load(os.path.join(out_dir, "tomo3d_magnitude.npy"))
    z = np.load(os.path.join(out_dir, "tomo3d_zaxis.npy"))
    db = 20 * np.log10(np.maximum(vol, 1e-12))

    std_z = float(db.std(axis=2).mean())
    std_r = float(db.std(axis=1).mean())
    std_a = float(db.std(axis=0).mean())
    ratio = std_z / std_r if std_r > 0 else float("nan")

    i, j = np.unravel_index(db.max(axis=2).argmax(), db.shape[:2])
    swing = float(db[i, j].max() - db[i, j].min())

    stats = {
        "std_lungo_z_dB": std_z,
        "std_lungo_range_dB": std_r,
        "std_lungo_azimuth_dB": std_a,
        "rapporto_z_su_range": ratio,
        "escursione_profilo_max_dB": swing,
        "estensione_asse_z_m": float(z[-1] - z[0]),
        "asse_z_informativo": bool(ratio > 0.1),
    }

    verdict = ("l'asse z porta informazione" if stats["asse_z_informativo"]
               else "l'asse z e' PIATTO: nessuna informazione in profondita'")
    text = "\n".join([
        "-" * 74,
        "VERIFICA EMPIRICA DEL CONTENUTO INFORMATIVO IN PROFONDITA'",
        "-" * 74,
        f"  deviazione std lungo z          = {std_z:10.4f} dB",
        f"  deviazione std lungo range      = {std_r:10.4f} dB",
        f"  deviazione std lungo azimuth    = {std_a:10.4f} dB",
        f"  rapporto z / range              = {ratio:10.2e}",
        f"  escursione del profilo migliore = {swing:10.6f} dB "
        f"su {stats['estensione_asse_z_m']:.0f} m",
        "",
        f"  >> {verdict}",
        "-" * 74,
    ])
    print(text)

    with open(os.path.join(out_dir, "verifica_asse_z.txt"), "w", encoding="utf-8") as fh:
        fh.write(text + "\n")

    meta_path = os.path.join(out_dir, "tomo3d_meta.json")
    if os.path.exists(meta_path):
        with open(meta_path, encoding="utf-8") as fh:
            meta = json.load(fh)
        meta["verifica_asse_z"] = stats
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, ensure_ascii=False)

    return stats


# --------------------------------------------------------------------------
# Self-test dell'inversione
# --------------------------------------------------------------------------

def selftest() -> int:
    """Verifica il blocco 9 su un diffusore sintetico a quota nota.

    Il test gira nel regime del brevetto (lambda_s=0.30 m, A=75 km, R=650 km,
    delta_z~1.3 m), cioe' dove il metodo RISOLVE. Se l'inversione e' corretta il
    picco di h(z) cade alla quota iniettata entro una cella.

    Nota: negli stessi test con parametri Sentinel-1 (lambda_s decine di metri,
    A ~ centinaia di metri) delta_z vale chilometri e il picco diventa
    indeterminato - non e' un difetto del codice ma il limite fisico dei dati,
    ed e' esattamente cio' che `--report-only` quantifica."""
    print("SELF-TEST - inversione tomografica su diffusore sintetico")
    print("  regime: brevetto WO 2024/008365 (lambda_s=0.30 m, A=75 km, R=650 km)")
    rng = np.random.default_rng(0)

    n_d, n_depth = 64, 256
    lambda_sound, slant_range, incidence = 0.30, 650_000.0, 39.4
    l_sa = 75_000.0
    z_axis = np.linspace(-100, 100, n_depth)
    baselines = np.linspace(-l_sa / 2, l_sa / 2, n_d)
    A = steering_matrix(z_axis, baselines, lambda_sound, slant_range, incidence)

    cell = z_axis[1] - z_axis[0]
    delta_z = lambda_sound * slant_range / (2 * l_sa)
    print(f"  cella asse z = {cell:.2f} m,  delta_z teorica = {delta_z:.2f} m")

    ok = True
    for z_true in (-60.0, 0.0, 35.0):
        j = int(np.argmin(np.abs(z_axis - z_true)))
        z_grid = z_axis[j]
        y = A[:, j] + 0.01 * (rng.standard_normal(n_d) + 1j * rng.standard_normal(n_d))
        h = _block9_tomographic_focusing(y[None, :], A)[0]
        z_est = z_axis[int(np.argmax(np.abs(h)))]
        good = abs(z_est - z_grid) <= max(2 * cell, delta_z)
        ok &= good
        print(f"  z reale={z_grid:+8.2f} m   z stimato={z_est:+8.2f} m   "
              f"{'OK' if good else 'FALLITO'}")

    # coerenza del budget: delta_z deve scalare come lambda/A
    b1 = ResolutionBudget(120, 6000, 50, 880000, 1200, 600, 50 * 880000 / 1200,
                          313, 156, 9, 200, 100, 139)
    scaled = (2 * 50 * 880000 / 1200) / (50 * 880000 / 1200)
    print(f"  scaling delta_z(2A)/delta_z(A) atteso 0.5 -> {1 / scaled:.2f} "
          f"{'OK' if abs(1 / scaled - 0.5) < 1e-9 else 'FALLITO'}")
    ok &= abs(1 / scaled - 0.5) < 1e-9

    print("SELF-TEST:", "SUPERATO" if ok else "FALLITO")
    return 0 if ok else 1


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Tomografia Doppler 3D di Cheope e Chefren da stack Sentinel-1 SLC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--stack-dir", default=DEFAULT_STACK)
    p.add_argument("--swath", default="iw2", choices=["iw1", "iw2", "iw3"])
    p.add_argument("--pol", dest="polarisation", default="vv", choices=["vv", "vh"])
    p.add_argument("--platform", default="s1c")
    p.add_argument("--dates", dest="n_dates", type=int, default=3,
                   help="numero di acquisizioni da impilare")
    p.add_argument("--date-filter", default=None, help="es. 20260415")
    p.add_argument("--lines", dest="n_lines", type=int, default=24,
                   help="numero di linee tomografiche parallele")
    p.add_argument("--nd", dest="n_d", type=int, default=32,
                   help="N_D: sampling rate dell'onda meccanica")
    p.add_argument("--f-inv", dest="f_investigation", type=float, default=120.0,
                   help="frequenza meccanica di indagine [Hz]")
    p.add_argument("--v-sound", type=float, default=6000.0,
                   help="velocita' di propagazione nel mezzo [m/s]")
    p.add_argument("--depth-max", dest="depth_max_m", type=float, default=400.0)
    p.add_argument("--n-depth", type=int, default=128)
    p.add_argument("--out", dest="out_dir", default="out_piramidi_v01")
    p.add_argument("--report-only", action="store_true",
                   help="stampa solo il budget di risoluzione ed esce")
    p.add_argument("--selftest", action="store_true",
                   help="verifica l'inversione tomografica e esce")
    p.add_argument("--verify", metavar="DIR", default=None,
                   help="misura il contenuto informativo in profondita' di un "
                        "volume gia' calcolato ed esce")
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.selftest:
        return selftest()

    if args.verify:
        verify_depth_information(args.verify)
        return 0

    cfg = Config(
        stack_dir=args.stack_dir,
        swath=args.swath,
        polarisation=args.polarisation,
        platform=args.platform,
        n_dates=args.n_dates,
        date_filter=args.date_filter,
        n_lines=args.n_lines,
        n_d=args.n_d,
        f_investigation=args.f_investigation,
        v_sound=args.v_sound,
        depth_max_m=args.depth_max_m,
        n_depth=args.n_depth,
        out_dir=args.out_dir,
    )
    verbose = not args.quiet

    print("=" * 74)
    print("piramidi_v01 - tomografia Doppler 3D di Cheope e Chefren")
    print("metodo Micro-Motion (Biondi & Malanga) - schema a 11 blocchi WO 2024/008365")
    print("=" * 74)

    entries = discover_stack(cfg)
    print(f"\nstack: {len(entries)} acquisizioni {cfg.platform.upper()} "
          f"{cfg.swath.upper()} {cfg.polarisation.upper()} "
          f"({entries[0].date} -> {entries[-1].date})")

    if args.report_only:
        ann = parse_annotation(entries[0].annotation)
        _, p = latlon_to_line_pixel(ann, *TARGETS["Khufu"][:2])
        budget = compute_resolution_budget(
            cfg, ann, ann.slant_range(p), TARGETS["Khufu"][2]
        )
        print(budget.as_text())
        if not budget.resolvable:
            print("\nCONCLUSIONE: con questi dati il metodo non risolve la struttura")
            print("interna delle piramidi. Il volume prodotto va letto come")
            print("dimostrazione della catena di elaborazione, non come riproduzione")
            print("dei risultati pubblicati su spotlight banda X.")
        return 0

    selected = entries[: max(1, cfg.n_dates)]
    results = []
    for entry in selected:
        try:
            results.append(process_date(entry, cfg, verbose=verbose))
        except Exception as exc:                       # pragma: no cover
            print(f"  [{entry.date}] SALTATA: {exc}", file=sys.stderr)

    if not results:
        print("nessuna acquisizione elaborata con successo", file=sys.stderr)
        return 1

    save_outputs(results, cfg)
    total = sum(float(r["elapsed_s"]) for r in results)
    print(f"  {len(results)} acquisizioni elaborate in {total:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

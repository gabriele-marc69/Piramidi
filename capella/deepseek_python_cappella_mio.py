#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
 SAR Doppler Tomography — Estrazione Immagine 3D Tomografica (Magnitude)
=============================================================================
 Implementazione della pipeline descritta in:
   Biondi, F.; Malanga, C. "Synthetic Aperture Radar Doppler Tomography
   Reveals Details of Undiscovered High-Resolution Internal Structure of
   the Great Pyramid of Giza". Remote Sens. 2022, 14, 5231.
   https://doi.org/10.3390/rs14205231

 Pipeline:
   1. Lettura file CPHD (sarpy per Capella / h5py per HDF5)
   2. Focusing SLC                     [Eq. 1-2]  (backprojection CPHD FX su GPU)
   3. DFT 2D                           [Eq. 3]
   4. Doppler Sub-Aperture Strategy    [Eq. 4-5, Fig. 7]
   5. Micro-Motion (Pixel Tracking)    [Eq. 6-12]
   6. Tomographic Focusing             [Eq. 22-24]
   7. Magnitude Image
   8. Visualizzazione e salvataggio    [Fig. 29-38b]

 Uso:
   python deepseek_python_cappella_mio.py --cphd input.cphd --out output/
   python deepseek_python_cappella_mio.py --demo          # genera dati sintetici
=============================================================================
"""

import os
import sys
import time
import argparse
import logging
import base64
import dataclasses
import json
import re
from dataclasses import dataclass
from typing import Optional, Tuple, List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator

# --- torch/CUDA per i calcoli pesanti (FFT, sub-aperture, matmul) ---
try:
    import torch
    USE_GPU = torch.cuda.is_available()
except ImportError:
    torch = None
    USE_GPU = False


def _gpu_fallback(logger: logging.Logger, what: str, e: Exception) -> None:
    """Log + svuota la cache GPU; il chiamante prosegue su CPU."""
    logger.warning(f"  GPU fallita in {what} ({type(e).__name__}: "
                   f"{str(e)[:80]}) -> uso la CPU")
    torch.cuda.empty_cache()


# --- sarpy (lettura CPHD Capella / NGA standard) ---
try:
    from sarpy.io.phase_history.cphd import CPHDReader
    SARPY_AVAILABLE = True
except ImportError:
    SARPY_AVAILABLE = False

# --- h5py (fallback per CPHD in formato HDF5) ---
try:
    import h5py
    H5PY_AVAILABLE = True
except ImportError:
    H5PY_AVAILABLE = False


# =============================================================================
#  CONFIGURAZIONE
# =============================================================================

@dataclass
class Config:
    """Parametri di configurazione della pipeline."""

    # --- Parametri di acquisizione SAR (Tabella 1 del paper) ---
    fc: float = 9.6e9              # frequenza centrale [Hz]
    bw: float = 400e6              # chirp bandwidth [Hz]
    B_cD: float = 22e3             # Doppler bandwidth [Hz]
    prf: float = 2.0e3             # PRF [Hz]
    prt: float = 0.23e-3           # PRT [s]
    antenna_length: float = 6.0    # lunghezza antenna [m]
    v_platform: float = 7.0e3      # velocità piattaforma [m/s]
    r_ref: float = 650e3           # slant range di riferimento [m]

    # --- Dominio dati CPHD ---
    domain_type: str = "FX"        # "FX" (freq range, time az) o "TOA"
    # Se FX, i dati sono già in frequenza in range → serve FFT in range.

    # --- Parametri tomografici ---
    N_D: int = 8                   # numero di sub-aperture Doppler
    B_DL_ratio: float = 0.5        # B_DL = B_cD * ratio
    z_min: float = -3000.0         # profondità minima [m]
    z_MAX: float = 0.0             # superficie [m]
    dz: float = 0.92               # risoluzione tomografica [m]
    A_orbit: float = 42_000.0      # apertura orbitale [m]
    f_investigation: float = 12_500.0  # frequenza onde sismiche [Hz]
    v_seismic: float = 6000.0      # velocità onde sismiche [m/s]
    steer_2pi: bool = False        # True: exp(j*2*pi*Kz*z) (lettura letterale del paper);
                                   # False: exp(j*Kz*z), con Kz in rad/m (coerente con dz=lambda*R/2A)

    # --- Linea tomografica e pixel tracking (Sez. 4, blocchi 4 e 7) ---
    line_orient: str = "range"     # "range": riga a azimuth fisso; "azimuth": colonna
    line_pos: Optional[int] = None # indice riga/colonna; None = la piu' energetica
    line_len: int = 512            # max pixel "puri" sulla linea
    pure_ratio: float = 0.3        # |p|^2 / energia 3x3 minima per un pixel "puro"
    patch: int = 16                # lato [pixel] della patch di tracking

    # --- Parametri di processing ---
    polarization: str = "HH"       # i file Capella sono HH
    window: str = "rect"           # "rect" (Fig. 7 del paper) o "hann"
    B_shift_ratio: float = 1 / 3   # B_shift = B_DL * ratio (usato se chain_len = 1)
    chain_len: int = 0             # 0 = strategia del paper: B_DL diviso in N_D passi, una
                                   # catena slave -> master successivo di N_D coppie.
                                   # >=1: coppie per posizione (N_D x chain_len misure)
    k_frac: Optional[np.ndarray] = None   # posizioni normalizzate delle misure
                                          # (riempito da estimate_micromotion)
    lpf_sigma: float = 1.0         # passa-basso gaussiano su Y(k) dopo il tracking [campioni]; 0 = off
    shift_px: int = 0              # >0: shift master/slave fisso in px (es. 1) invece di B_DL/N_D
    osc: bool = False              # True: Y(k) = traiettoria di una massa unitaria (oscillatore
                                   # armonico smorzato) spinta dal vettore di tracking (a + jb)
    osc_cycles: float = 4.0        # cicli dell'oscillatore sull'intera catena (w0 = 2*pi*cicli/N)
    osc_damp: float = 0.1          # smorzamento lambda in unita' di w0
    vec_waves: bool = False        # onde da catena di delta di N_D vettori (Re, Im) per pixel
    pyr_height: float = 138.5      # altezza su cui stirare le onde [m] (Khufu)
    z_bottom: float = 0.0          # quota di arrivo dell'onda [m] (< 0 = sotto il livello zero)
    wave_lpf: float = np.pi / 2    # taglio del passa-basso sulle onde [rad per passo] (meta' Nyquist)
    wave_half: float = 100.0       # discesa [m] in cui l'onda di un tiro si dimezza (smorzamento)
    upsample_factor: int = 100

    # --- Output ---
    out_dir: str = "output"
    dpi: int = 300
    cmap: str = "jet"
    log_scale: bool = False
    save_npy: bool = True

    def __post_init__(self):
        self.lambda_sound = self.v_seismic / self.f_investigation
        self.delta_T = (self.lambda_sound * self.r_ref) / (2 * self.A_orbit)
        self.z = np.arange(self.z_min, self.z_MAX + self.dz, self.dz)

    def summary(self) -> str:
        return (
            f"=== Configurazione ===\n"
            f"  fc                : {self.fc/1e9:.2f} GHz\n"
            f"  B_cr              : {self.bw/1e6:.1f} MHz\n"
            f"  B_cD              : {self.B_cD/1e3:.1f} kHz\n"
            f"  PRF               : {self.prf/1e3:.1f} kHz\n"
            f"  V piattaforma     : {self.v_platform/1e3:.1f} km/s\n"
            f"  R rif.            : {self.r_ref/1e3:.0f} km\n"
            f"  Dominio dati      : {self.domain_type}\n"
            f"  N_D sub-aperture  : {self.N_D}\n"
            f"  λ (vibrazionale)  : {self.lambda_sound:.3f} m\n"
            f"  δ_T (tomografica) : {self.delta_T:.3f} m\n"
            f"  z range           : [{self.z_min:.0f}, {self.z_MAX:.0f}] m "
            f"({len(self.z)} campioni)\n"
            f"  Polarizzazione    : {self.polarization}\n"
        )


# =============================================================================
#  LOGGING
# =============================================================================

def setup_logger(out_dir: str, verbose: bool = False) -> logging.Logger:
    os.makedirs(out_dir, exist_ok=True)
    logger = logging.getLogger("SAR-Tomo")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s",
                            datefmt="%H:%M:%S")

    # output rediretto su file: la codifica di default (cp1252) non ha "λ", "δ", "→"
    # e ogni messaggio che li contiene darebbe un "Logging error"
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    fh = logging.FileHandler(os.path.join(out_dir, "processing.log"),
                             encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    return logger


# =============================================================================
#  STEP 1 — LETTURA CPHD (sarpy + fallback h5py)
# =============================================================================

def _try_read_with_sarpy(filepath: str, logger: logging.Logger,
                         pulse_start: int = 0, n_pulses: Optional[int] = None):
    """Lettura CPHD via sarpy (Capella, NGA standard), a blocco di impulsi."""
    reader = CPHDReader(filepath)
    xml_dict = reader.cphd_meta.to_dict()

    # --- Global ---
    g = xml_dict.get("Global", {})

    # --- Channel (primo canale) ---
    ch_list = xml_dict["Data"]["Channels"]
    ch = ch_list[0]

    # --- PVP disponibili ---
    pvp_names = list(xml_dict.get("PVP", {}).keys())
    logger.debug(f"PVP disponibili: {pvp_names}")

    # --- Estrazione metadati principali ---
    # CPHD 1.0: la banda in frequenza e' Global/FxBand {FxMin, FxMax}
    fc = bw = None
    fxb = g.get("FxBand")
    if fxb:
        fc = 0.5 * (float(fxb["FxMin"]) + float(fxb["FxMax"]))
        bw = float(fxb["FxMax"]) - float(fxb["FxMin"])

    # PRF: stimata piu' sotto dai tempi di trasmissione (TxTime)
    prf = None

    # Numero campioni/impulsi
    num_samples = int(ch.get("NumSamples", 0))
    num_vectors = int(ch.get("NumVectors", 0))

    # Polarizzazione: CPHD 1.0 la mette in Channel/Parameters[i]/Polarization {TxPol, RcvPol}
    # (Data/Channels contiene solo dimensioni e offset)
    pol = (xml_dict.get("Channel", {}).get("Parameters") or [{}])[0].get("Polarization") or {}
    tx_pol, rcv_pol = pol.get("TxPol"), pol.get("RcvPol")
    polarization = f"{tx_pol}{rcv_pol}" if tx_pol and rcv_pol else "UNK"

    # --- Lettura dati grezzi (phase history) ---
    # Il file intero e' enorme (es. 311079 x 27000 c64 = 62 GiB): leggo solo
    # un blocco di impulsi [pulse_start, pulse_start + n_pulses).
    p0 = max(0, int(pulse_start))
    p1 = num_vectors if n_pulses is None else min(num_vectors, p0 + int(n_pulses))
    pulse_slice = slice(p0, p1)
    est_gb = (p1 - p0) * num_samples * 8 / 2**30
    logger.info(f"  blocco impulsi : {p0}..{p1} di {num_vectors} "
                f"(~{est_gb:.2f} GiB raw)")
    raw_data = reader.read_chip(pulse_slice, slice(0, num_samples))
    num_vectors = p1 - p0

    # --- Lettura PVP per geometria ---
    pvp = {}
    for name in ("TxPos", "RcvPos", "SRPPos", "TxTime", "RcvTime",
                 "SRPTime", "TxPulseNumber"):
        if name in pvp_names:
            try:
                pvp[name] = np.ascontiguousarray(
                    reader.read_pvp_variable(name, 0, pulse_slice)
                )
            except Exception as e:
                logger.debug(f"PVP {name} non leggibile: {e}")

    # --- Calcolo r_ref e v_platform ---
    r_ref = None
    v_platform = None

    # PRF reale dai tempi di trasmissione
    if "TxTime" in pvp and pvp["TxTime"].shape[0] > 1:
        dt = float(np.median(np.diff(pvp["TxTime"].astype(np.float64))))
        if dt > 0:
            prf = 1.0 / dt

    # Slant range = distanza piattaforma -> SRP (NON |SRPPos|, che e' la
    # distanza dal centro della Terra, ~6370 km)
    if "SRPPos" in pvp and "TxPos" in pvp and pvp["SRPPos"].shape[0] > 1:
        srp = pvp["SRPPos"].astype(np.float64)  # [N, 3] ECF in metri
        tx = pvp["TxPos"].astype(np.float64)
        r_ref = float(np.mean(np.linalg.norm(tx - srp, axis=1)))
    elif "TxPos" in pvp and "RcvPos" in pvp:
        tx = pvp["TxPos"].astype(np.float64)
        rx = pvp["RcvPos"].astype(np.float64)
        r_ref = float(np.mean(np.linalg.norm(tx - rx, axis=1)))

    if "TxPos" in pvp and pvp["TxPos"].shape[0] > 1 and prf:
        tx = pvp["TxPos"].astype(np.float64)
        dr = np.diff(tx, axis=0)
        v_platform = float(np.mean(np.linalg.norm(dr, axis=1)) * prf)

    # --- Dominio dati ---
    domain_type = g.get("DomainType", ch.get("DomainType", "FX"))

    # --- Banda Doppler effettiva (stima) ---
    # B_cD = 2 * V * La / (λ * R)  — usiamo il valore di default se manca
    # (Bw e' la banda chirp in range, NON la banda Doppler: non usarla qui)
    B_cD = None  # None -> main() mantiene il default del paper (22 kHz)

    meta = {
        "fc":           fc,
        "bw":           bw,
        "prf":          prf,
        "B_cD":         B_cD,
        "num_pulses":   num_vectors,
        "num_samples":  num_samples,
        "polarization": polarization,
        "v_platform":   v_platform,
        "r_ref":        r_ref,
        "domain_type":  domain_type,
        "pvp":          pvp,
    }

    logger.info(f"  raw_data shape : {raw_data.shape}")
    logger.info(f"  dtype          : {raw_data.dtype}")
    logger.info(f"  polarizzazione : {polarization}")
    logger.info(f"  dominio        : {domain_type}")
    if fc:  logger.info(f"  fc             : {fc/1e9:.3f} GHz")
    if bw:  logger.info(f"  bw             : {bw/1e6:.1f} MHz")
    if prf: logger.info(f"  PRF            : {prf/1e3:.3f} kHz")
    if r_ref:      logger.info(f"  r_ref          : {r_ref/1e3:.1f} km")
    if v_platform: logger.info(f"  v_platform     : {v_platform:.1f} m/s")

    return raw_data, meta


def _try_read_with_h5py(filepath: str, logger: logging.Logger):
    """Fallback: CPHD in formato HDF5 (alcuni vendor)."""
    with h5py.File(filepath, "r") as f:
        data = f["/CPHD/0/Data"][:]
        meta = {
            "fc":           float(f["/CPHD/0/Global/fc"][()])
                            if "/CPHD/0/Global/fc" in f else None,
            "bw":           float(f["/CPHD/0/Global/Bw"][()])
                            if "/CPHD/0/Global/Bw" in f else None,
            "prf":          float(f["/CPHD/0/Global/PRF"][()])
                            if "/CPHD/0/Global/PRF" in f else None,
            "B_cD":         float(f["/CPHD/0/Global/B_cD"][()])
                            if "/CPHD/0/Global/B_cD" in f else 22e3,
            "num_pulses":   data.shape[0],
            "num_samples":  data.shape[1],
            "polarization": f["/CPHD/0/Global/Polarization"][()].decode()
                            if "/CPHD/0/Global/Polarization" in f else "UNK",
            "v_platform":   float(f["/CPHD/0/Global/PlatformVelocity"][()])
                            if "/CPHD/0/Global/PlatformVelocity" in f else None,
            "r_ref":        float(f["/CPHD/0/Global/ReferenceRange"][()])
                            if "/CPHD/0/Global/ReferenceRange" in f else None,
            "domain_type":  "FX",
            "pvp":          {},
        }
    return data, meta


def read_cphd(filepath: str, logger: logging.Logger,
              pulse_start: int = 0,
              n_pulses: Optional[int] = None) -> Tuple[np.ndarray, dict]:
    """
    Legge un file CPHD provando prima sarpy (Capella/NGA) e poi h5py.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"File CPHD non trovato: {filepath}")

    # --- Prova sarpy ---
    if SARPY_AVAILABLE:
        try:
            logger.info(f"Lettura CPHD (sarpy): {filepath}")
            return _try_read_with_sarpy(filepath, logger, pulse_start, n_pulses)
        except Exception as e:
            logger.warning(f"sarpy ha fallito: {e}")
            logger.warning("Provo fallback h5py...")

    # --- Fallback h5py (solo se il file e' davvero HDF5) ---
    with open(filepath, "rb") as fh:
        is_hdf5 = fh.read(8) == b"\x89HDF\r\n\x1a\n"
    if H5PY_AVAILABLE and is_hdf5:
        try:
            logger.info(f"Lettura CPHD (h5py): {filepath}")
            return _try_read_with_h5py(filepath, logger)
        except Exception as e:
            logger.warning(f"h5py ha fallito: {e}")

    raise RuntimeError(
        "Impossibile leggere il file CPHD. Installare 'sarpy' "
        "(pip install sarpy) oppure verificare il formato del file."
    )


# =============================================================================
#  UTILITA' DI MEMORIA (immagini enormi: 8 GB di RAM)
# =============================================================================

def _decimate(a: np.ndarray, max_px: int = 1024) -> np.ndarray:
    """|a| ridotto a block-max, al massimo ~max_px per lato (float32).
    Lavora a blocchi di righe: non crea mai |a| intero."""
    H, W = a.shape
    by = max(1, -(-H // max_px))
    bx = max(1, -(-W // max_px))
    h2, w2 = (H // by) * by, (W // bx) * bx
    out = np.empty((h2 // by, w2 // bx), dtype=np.float32)
    rows = 64 * by
    for r in range(0, h2, rows):
        blk = np.abs(a[r:min(r + rows, h2), :w2])
        n = blk.shape[0] // by
        out[r // by:r // by + n] = blk.reshape(n, by, w2 // bx, bx).max(axis=(1, 3))
    return out


# =============================================================================
#  STEP 2 — FOCUSING SLC
# =============================================================================

def focus_slc(raw_data: np.ndarray, cfg: Config,
              logger: logging.Logger) -> np.ndarray:
    """
    Focusing SAR:
      - Se dominio "FX": dati già in frequenza in range → FFT range
      - Se dominio "TOA": segnale gia' campionato in tempo d'arrivo, cioe' gia'
        compresso in range (CPHD 1.0): nessuna operazione in range
      - Azimuth: matched filter Doppler
    Output: S_SLC(k, x) [Eq. 1-2].
    """
    logger.info(f"Step 2: Focusing SLC (dominio={cfg.domain_type})...")
    raw_data = raw_data.astype(np.complex64, copy=False)
    n_pulses = raw_data.shape[0]

    lambda_em = 299_792_458.0 / cfg.fc
    ta = np.arange(n_pulses) / cfg.prf
    ka = 2 * cfg.v_platform**2 / (lambda_em * cfg.r_ref)
    az_ref = np.exp(1j * np.pi * ka * (ta - ta.mean())**2)
    phase_corr = np.complex64(np.exp(-1j * 4 * np.pi * cfg.r_ref / lambda_em))

    # --- Percorso GPU (stesse operazioni, FFT su CUDA) ---
    if USE_GPU:
        try:
            dev = "cuda"
            x = torch.from_numpy(np.ascontiguousarray(raw_data)).to(dev)
            if cfg.domain_type.upper() == "FX":
                x = torch.fft.ifft(x, dim=1)
            az_t = torch.from_numpy(np.fft.fft(az_ref).astype(np.complex64)).to(dev)
            x = torch.fft.fft(x, dim=0)
            x *= torch.conj(az_t)[:, None]
            x = torch.fft.ifft(x, dim=0)
            x *= complex(phase_corr)
            max_val = float(x.abs().max())
            if max_val > 0:
                x /= max_val
            slc = x.cpu().numpy()
            del x
            torch.cuda.empty_cache()
            logger.info(f"  SLC shape: {slc.shape}  [GPU]")
            return slc
        except RuntimeError as e:
            _gpu_fallback(logger, "focus_slc", e)

    # --- 2a. Range ---
    if cfg.domain_type.upper() == "FX":
        # Dati in frequenza → FFT per tornare nel tempo (range compression)
        range_compressed = np.fft.ifft(raw_data, axis=1)
    else:
        # Dominio TOA: i vettori sono gia' in tempo d'arrivo (compressi in range)
        range_compressed = raw_data.copy()

    # --- 2b. Azimuth compression (az_ref: lambda EM = c/fc, vedi sopra) ---
    # complex64 e operazioni in-place: niente copie complex128 da GB
    az_ref_fft = np.fft.fft(az_ref).astype(np.complex64)
    slc = np.fft.fft(range_compressed, axis=0)
    del range_compressed
    slc *= np.conj(az_ref_fft)[:, None]
    slc = np.fft.ifft(slc, axis=0)

    # --- 2c. Correzione fase residua ---
    slc *= phase_corr

    # Normalizzazione (per blocchi: evita |slc| intero)
    max_val = max(float(np.abs(slc[i:i + 256]).max())
                  for i in range(0, n_pulses, 256))
    if max_val > 0:
        slc /= np.float32(max_val)

    logger.info(f"  SLC shape: {slc.shape}")
    return slc


# =============================================================================
#  STEP 2 (alternativo) — FOCUSING PER BACKPROJECTION (CPHD FX, spotlight)
# =============================================================================

_C = 299_792_458.0


def _bp_geometry(pv: dict, n: int, N: int, d_sc: float, bw: float,
                 scene_az: float, scene_rg: float, max_px: int,
                 logger: logging.Logger, cx: float = 0.0, cy: float = 0.0,
                 az_os: int = 1):
    """Piano immagine al suolo attraverso l'SRP e geometria per impulso.

    Coordinate RELATIVE all'SRP (float64 qui): in float32 le coordinate ECF
    (~6.4e6 m) non danno la precisione di fase necessaria (lambda ~ 3 cm).
    az_os > 1: passo dei pixel in azimuth = risoluzione / az_os (sovracampionato)."""
    srp = pv["SRPPos"].mean(axis=0)
    if np.ptp(pv["SRPPos"], axis=0).max() > 1.0:
        logger.warning("  SRPPos non costante nel blocco: uso la media")

    a_e, b_e = 6378137.0, 6356752.314245                 # WGS84
    nrm = np.array([srp[0] / a_e**2, srp[1] / a_e**2, srp[2] / b_e**2])
    up = nrm / np.linalg.norm(nrm)                        # normale ellissoide

    mid = n // 2
    los = srp - pv["TxPos"][mid]
    los /= np.linalg.norm(los)
    los_h = los - (los @ up) * up
    u_rg = los_h / np.linalg.norm(los_h)                  # range crescente
    u_az = np.cross(up, u_rg)
    if u_az @ pv["TxVel"][mid] < 0:
        u_az = -u_az                                      # azimuth lungo il moto
    inc = float(np.arccos(max(-1.0, min(1.0, -los @ up))))  # angolo da verticale

    fc = pv["fc"]
    lam = _C / fc
    e0 = pv["TxPos"][0] - srp
    e1 = pv["TxPos"][-1] - srp
    theta = float(np.arccos(np.clip(e0 @ e1 / np.linalg.norm(e0) /
                                    np.linalg.norm(e1), -1, 1)))
    da = lam / (2 * max(theta, 1e-9))                     # risoluzione azimuth
    dr = _C / (2 * bw) / np.sin(inc)                      # risoluzione ground-range
    res_az = da
    da = da / az_os                                       # passo dei pixel in azimuth
    n_az = max(8, int(scene_az / da) | 1)
    n_rg = max(8, int(scene_rg / dr) | 1)
    if n_az * n_rg > max_px:
        k = np.sqrt(max_px / (n_az * n_rg))
        n_az, n_rg = max(8, int(n_az * k) | 1), max(8, int(n_rg * k) | 1)
        logger.warning(f"  scena troppo grande: ridotta a {n_az}x{n_rg} px")
    logger.info(f"  geometria: inc={np.degrees(inc):.1f} deg, apertura "
                f"{theta * 1e3:.2f} mrad -> res. azimuth {res_az:.2f} m "
                f"(pixel {da:.2f} m), "
                f"ground-range {dr:.2f} m; immagine {n_az}x{n_rg} px "
                f"({n_az * da:.0f} x {n_rg * dr:.0f} m)")

    def path(pos):
        t = pos - srp
        tn = np.linalg.norm(t, axis=1)
        u = t / tn[:, None]
        return tn, u @ u_rg, u @ u_az

    # cx, cy: centro dell'immagine rispetto all'SRP (solo la griglia dei pixel si
    # sposta; la fase resta riferita all'SRP, come nei dati)
    return dict(path_tx=path(pv["TxPos"]), path_rcv=path(pv["RcvPos"]),
                cx=cx, cy=cy, da=da, dr=dr, n_az=n_az, n_rg=n_rg, inc=inc, lam=lam,
                res_az=res_az,
                srp=srp, u_rg=u_rg, u_az=u_az, up=up)


def _bp_run(rd, p0: int, n: int, N: int, geo: dict, pv: dict, f_c: float,
            d_sc: float, sgn: int, os_: int, chunk: int, dev: str,
            logger: logging.Logger, tag: str = "") -> np.ndarray:
    """Backprojection a blocchi di impulsi; il segnale e' letto in streaming."""
    n_az, n_rg = geo["n_az"], geo["n_rg"]
    M = N * os_
    h = N // 2
    Y = ((torch.arange(n_az, device=dev, dtype=torch.float32) - n_az // 2)
         * geo["da"] + geo["cy"])[:, None]
    X = ((torch.arange(n_rg, device=dev, dtype=torch.float32) - n_rg // 2)
         * geo["dr"] + geo["cx"])[None, :]
    P2 = X * X + Y * Y
    win = torch.from_numpy(np.hamming(N).astype(np.float32)).to(dev)
    img = torch.zeros((n_az, n_rg), dtype=torch.complex64, device=dev)
    ttx, atx, btx = geo["path_tx"]
    trx, arx, brx = geo["path_rcv"]
    bins_per_m = d_sc * M / _C            # bin dell'asse di ritardo per metro di cammino
    cyc_per_m = f_c / _C
    ang_k = -sgn * 2.0 * np.pi
    amp = pv["AmpSF"]

    def d_path(tn, a, b):
        pu = Y * float(b) + X * float(a)               # u . P'
        num = P2 - (2.0 * float(tn)) * pu
        r = torch.sqrt(num + float(tn) * float(tn))
        return num / (r + float(tn))                   # |T-P| - |T-SRP|

    t0 = time.time()
    for c0 in range(0, n, chunk):
        c1 = min(n, c0 + chunk)
        S = rd.read_chip(slice(p0 + c0, p0 + c1), slice(0, N))
        St = torch.from_numpy(np.ascontiguousarray(S)).to(dev)
        St = St * win * torch.from_numpy(amp[c0:c1].astype(np.float32)).to(dev)[:, None]
        A = torch.zeros((c1 - c0, M), dtype=torch.complex64, device=dev)
        A[:, :N - h] = St[:, h:]
        A[:, M - h:] = St[:, :h]
        q = torch.fft.fft(A, dim=1) if sgn > 0 else torch.fft.ifft(A, dim=1) * M
        del A, St
        for j in range(c1 - c0):
            k = c0 + j
            dR = d_path(ttx[k], atx[k], btx[k]) + d_path(trx[k], arx[k], brx[k])
            bf = torch.remainder(dR * bins_per_m, float(M))
            i0f = torch.floor(bf)
            fr = bf - i0f
            i0 = torch.remainder(i0f.long(), M)   # remainder float32 puo' dare M
            i1 = torch.remainder(i0 + 1, M)
            qj = q[j]
            v = qj[i0] * (1.0 - fr) + qj[i1] * fr
            ang = torch.remainder(dR * cyc_per_m, 1.0) * ang_k
            img += v * torch.complex(torch.cos(ang), torch.sin(ang))
        del q
        if tag == "" and (c0 // chunk) % 40 == 0:
            logger.info(f"    BP {c1}/{n} impulsi  ({time.time() - t0:.0f} s)")
    return img


def backproject_cphd(filepath: str, cfg: Config, logger: logging.Logger,
                     pulse_start: int, n_pulses: int,
                     scene_az: float = 600.0, scene_rg: float = 600.0,
                     os_: int = 8, max_px: int = 4_000_000,
                     center: Tuple[float, float] = (0.0, 0.0),
                     sign_mult: int = -1, center_on: Optional[str] = None,
                     az_os: int = 1):
    """Focalizza un blocco di impulsi di un CPHD FX (Capella spotlight) con
    backprojection al suolo, centrata sull'SRP (o sulla piramide `center_on`, con
    la scena spostata a meta' del layover dell'apice). Output: SLC (az x range)
    complex64, normalizzato a max=1, con spettro azimutale rettangolare."""
    logger.info("Step 2: Focusing SLC per BACKPROJECTION (CPHD FX)...")
    if torch is None:
        raise RuntimeError("torch non disponibile: la backprojection lo richiede")
    dev = "cuda" if USE_GPU else "cpu"
    if dev == "cpu":
        logger.warning("  GPU non disponibile: la backprojection su CPU e' molto lenta")

    rd = CPHDReader(filepath)
    d = rd.cphd_meta.to_dict()
    ch = d["Data"]["Channels"][0]
    N, total = int(ch["NumSamples"]), int(ch["NumVectors"])
    p0 = max(0, int(pulse_start))
    p1 = min(total, p0 + int(n_pulses))
    n = p1 - p0
    sl = slice(p0, p1)

    pv = {k: np.asarray(rd.read_pvp_variable(k, 0, sl), dtype=np.float64)
          for k in ("TxPos", "RcvPos", "TxVel", "SRPPos", "AmpSF", "SC0", "SCSS", "TxTime")}
    if np.ptp(pv["SC0"]) > 1.0 or np.ptp(pv["SCSS"]) > 1e-3:
        logger.warning("  SC0/SCSS non costanti nel blocco: uso il primo impulso")
    sc0, d_sc = float(pv["SC0"][0]), float(pv["SCSS"][0])
    bw = N * d_sc
    f_c = sc0 + (N // 2) * d_sc                    # frequenza del campione centrale
    pv["fc"] = sc0 + 0.5 * bw
    sgn_file = int(d["Global"].get("SGN", 1))
    logger.info(f"  blocco impulsi {p0}..{p1}  N={N}  banda={bw / 1e6:.0f} MHz  "
                f"fc={pv['fc'] / 1e9:.3f} GHz  SGN(file)={sgn_file}")

    geo = _bp_geometry(pv, n, N, d_sc, bw, scene_az, scene_rg, max_px, logger,
                       cx=float(center[0]), cy=float(center[1]), az_os=az_os)
    if center_on:
        px, _, H = _pyr_px_fn(geo, center_on)
        r0, c0 = px(0.0, 0.0)
        geo["cy"] += (r0 - geo["n_az"] // 2) * geo["da"]
        geo["cx"] += (c0 - geo["n_rg"] // 2) * geo["dr"] - H / np.tan(geo["inc"]) / 2

    # Segno di fase = SGN del file: lo standard CPHD definisce
    # Phase(fx) = SGN * fx * dTOA. Il segno opposto dara' l'immagine riflessa
    # attorno all'SRP (non distinguibile dalla sola nitidezza). Validato con un
    # CPHD simulato (bersagli puntiformi in posizioni note, geometria reale).
    sgn = sgn_file * int(sign_mult)
    logger.info(f"  segno di fase usato: {sgn:+d} (SGN file {sgn_file:+d} x {sign_mult:+d})")
    # orbita per impulso (TxTime: secondi dall'inizio della raccolta) per --bp19
    geo.update(sgn=sgn, f_c=f_c, t_pulse=pv["TxTime"], tx=pv["TxPos"], rcv=pv["RcvPos"],
               vel=pv["TxVel"])
    img = _bp_run(rd, p0, n, N, geo, pv, f_c, d_sc, sgn, os_, 64, dev, logger)
    slc = img.cpu().numpy()
    del img
    if dev == "cuda":
        torch.cuda.empty_cache()
    mx = float(np.abs(slc).max())
    if mx > 0:
        slc /= np.float32(mx)
    logger.info(f"  SLC (backprojection) shape: {slc.shape}")
    return slc, geo


# =============================================================================
#  STEP 3 — DFT 2D
# =============================================================================

def dft2(slc: np.ndarray, logger: logging.Logger) -> np.ndarray:
    """DFT 2D con shift. Output: spettro rettangolare [Eq. 3]."""
    logger.info("Step 3: DFT 2D...")
    if USE_GPU:
        try:
            x = torch.fft.fftshift(torch.fft.fft2(torch.from_numpy(slc).cuda()))
            out = x.cpu().numpy()
            del x
            torch.cuda.empty_cache()
            return out
        except RuntimeError as e:
            _gpu_fallback(logger, "dft2", e)
    return np.fft.fftshift(np.fft.fft2(slc))


def select_tomographic_line(slc: np.ndarray, cfg: Config,
                            logger: logging.Logger):
    """Blocco 4 (filtro spaziale): sceglie la linea tomografica sull'SLC e
    tiene solo i pixel "puri" (Sez. 5.2: la cella deve venire da un solo
    scatteratore, per attenuare multipath/layover).

    Linea = riga a azimuth fisso (line_orient="range") o colonna a range fisso
    ("azimuth"). Pixel "puro": massimo locale 3x3 con potenza >= pure_ratio
    dell'energia 3x3. Ritorna (rows, cols) interi, ordinati lungo la linea."""
    n_az, n_rg = slc.shape
    h = cfg.patch // 2
    horiz = cfg.line_orient.lower() == "range"
    n_pos, n_line = (n_az, n_rg) if horiz else (n_rg, n_az)

    pos = cfg.line_pos
    if pos is None:
        # riga/colonna piu' energetica (a blocchi: niente |slc| intero)
        if horiz:
            e = np.concatenate([np.abs(slc[r:r + 256]).mean(axis=1)
                                for r in range(0, n_az, 256)])
        else:
            e = np.zeros(n_rg)
            for r in range(0, n_az, 256):
                e += np.abs(slc[r:r + 256]).sum(axis=0)
        pos = h + int(np.argmax(e[h:n_pos - h + 1]))
    pos = int(min(max(pos, h + 1), n_pos - h - 1))

    band = (np.abs(slc[pos - 1:pos + 2, :]) if horiz
            else np.abs(slc[:, pos - 1:pos + 2]).T).astype(np.float64)
    p = band ** 2
    c = p[1, 1:-1]
    e3 = (p[:, :-2] + p[:, 1:-1] + p[:, 2:]).sum(axis=0)
    ratio = c / (e3 + 1e-30)
    is_max = ((c >= p[0, 1:-1]) & (c >= p[2, 1:-1]) &
              (c >= p[1, :-2]) & (c >= p[1, 2:]))
    t = np.arange(1, n_line - 1)
    ok = is_max & (ratio >= cfg.pure_ratio) & (t >= h) & (t <= n_line - h)
    n_pure = int(ok.sum())

    if n_pure < cfg.N_D:
        logger.warning(f"  solo {n_pure} pixel 'puri' sulla linea "
                       f"(pos={pos}): uso i pixel piu' forti")
        ok = (t >= h) & (t <= n_line - h)
    cand = t[ok][np.argsort(c[ok])[::-1]]       # dal piu' forte
    cand = np.sort(cand[:cfg.line_len])

    if horiz:
        rows, cols = np.full(cand.size, pos), cand
    else:
        rows, cols = cand, np.full(cand.size, pos)
    logger.info(f"  linea tomografica: {cfg.line_orient} pos={pos}, "
                f"{cand.size} pixel (candidati puri: {n_pure})")
    return rows.astype(int), cols.astype(int)


# =============================================================================
#  STEP 4 — DOPPLER SUB-APERTURE
# =============================================================================

def _gather_patches(sub, rows: np.ndarray, cols: np.ndarray, h: int) -> np.ndarray:
    """|patch| (2h x 2h) attorno a ogni pixel della linea -> (N, 2h, 2h) f32.
    Funziona su array numpy e su tensori torch (GPU)."""
    off = np.arange(-h, h)
    rr = rows[:, None, None] + off[None, :, None]
    cc = cols[:, None, None] + off[None, None, :]
    if torch is not None and isinstance(sub, torch.Tensor):
        rr = torch.from_numpy(rr).to(sub.device)
        cc = torch.from_numpy(cc).to(sub.device)
        return sub[rr, cc].abs().float().cpu().numpy()
    return np.abs(sub[rr, cc]).astype(np.float32)


_LAST_STRATEGY: dict = {}      # descrizione dell'ultima strategia (per il grafico)


def generate_doppler_subapartures(
    spectrum: np.ndarray,
    cfg: Config,
    logger: logging.Logger,
    rows: np.ndarray,
    cols: np.ndarray,
) -> Tuple[list, List[np.ndarray], List[np.ndarray]]:
    """Strategia Doppler a CATENA master -> slave [Eq. 4-5, Fig. 7, blocchi 3-6].

    Per ognuna delle N_D posizioni di partenza a_i, le finestre azimutali (larghezza
    B_cD - B_DL) sono posizionate a a_i, a_i+u, a_i+2u, ... (u = spostamento unitario
    = B_shift). Ogni immagine e' lo SLAVE della coppia precedente e il MASTER della
    successiva: con chain_len = c si ottengono c coppie (misure) per posizione e
    N_D x c misure in totale (8 x 7 = 56, contro 8 con la singola coppia), da c+1
    immagini per posizione invece di 2c.

    chain_len = 0 (default, strategia del paper): B_DL diviso in N_D passi, una sola
    catena di N_D coppie con u = B_DL/N_D (o cfg.shift_px px se > 0).
    chain_len = 1: una sola coppia per posizione, B_shift = B_DL * B_shift_ratio.
    chain_len > 1: u = (corsa disponibile) / (chain_len + 2), cosi' la catena resta
    dentro la banda; i passi tra posizioni, step = u / N_D, rendono distinte le
    N_D x c coppie.

    Blocchi 5-6 = IDFT2 (ifft2). Di ogni coppia si tengono solo le patch attorno ai
    pixel della linea e lo spostamento sub-pixel (dy, dx, qualita') calcolato subito,
    senza conservare le patch.
    Ritorna (misure, anteprime master, anteprime slave) della prima catena."""
    n_pulses, n_samples = spectrum.shape
    B_DL = cfg.B_cD * cfg.B_DL_ratio
    B_processed = cfg.B_cD - B_DL
    h = cfg.patch // 2
    paper = int(cfg.chain_len) <= 0                  # strategia del paper (Fig. 7)
    chain = cfg.N_D if paper else max(1, int(cfg.chain_len))
    bw_px = max(2, int(round(B_processed / cfg.B_cD * n_pulses)))
    T = max(0, n_pulses - bw_px)                     # corsa totale delle finestre

    if paper:
        # B_DL (lasciata fuori dal filtro) divisa in N_D passi uguali: la finestra
        # master/slave scorre di u = B_DL / N_D, una sola catena di N_D coppie
        u = int(cfg.shift_px) if cfg.shift_px > 0 else max(1, T // cfg.N_D)
        while u > 1 and chain * u > T:
            u -= 1
        starts0 = [0]
    elif chain == 1:
        u = max(1, int(round(B_DL * cfg.B_shift_ratio / cfg.B_cD * n_pulses)))
        max_shift = max(0, T - u)                    # lo slave resta dentro la banda
        starts0 = [int(round(i * max_shift / max(1, cfg.N_D - 1)))
                   for i in range(cfg.N_D)]
    else:
        u = max(1, T // (chain + 2))
        step = max(1, u // cfg.N_D)
        while u > 1 and (cfg.N_D - 1) * step + (chain + 1) * u > T:
            u -= 1
            step = max(1, u // cfg.N_D)
        starts0 = [i * step for i in range(cfg.N_D)]
    B_shift = u / n_pulses * cfg.B_cD
    n_pos = len(starts0)
    logger.info(f"Step 4: Doppler Sub-Aperture a catena"
                f"{' (paper: B_DL/N_D)' if paper else ''}: {n_pos} posizioni x "
                f"{chain} coppie = {n_pos * chain} misure ({n_pos * (chain + 1)} "
                f"immagini); B_shift={B_shift:.0f} Hz ({u} px)")

    def az_band(idx_start: int, idx_end: int) -> np.ndarray:
        idx_start = max(0, min(idx_start, n_pulses - 1))
        idx_end = max(idx_start + 1, min(idx_end, n_pulses))
        win = np.zeros(n_pulses)
        if cfg.window == "hann":
            win[idx_start:idx_end] = np.hanning(idx_end - idx_start)
        else:
            win[idx_start:idx_end] = 1.0             # finestra rettangolare (Fig. 7)
        return win

    state = {"gpu": USE_GPU, "spec": None}
    if state["gpu"]:
        try:
            state["spec"] = torch.fft.ifftshift(torch.from_numpy(spectrum).cuda())
        except RuntimeError as e:
            _gpu_fallback(logger, "sub-aperture (upload)", e)
            state["gpu"] = False
    if not state["gpu"]:
        state["spec"] = np.fft.ifftshift(spectrum)   # una sola copia

    def one_image(w: np.ndarray, want_prev: bool):
        """ifft2(spec * w) -> (patch attorno ai pixel, anteprima|None)."""
        if state["gpu"]:
            try:
                sub = torch.fft.ifft2(state["spec"] * torch.from_numpy(w).cuda()[:, None])
                pt = _gather_patches(sub, rows, cols, h)
                pv = _decimate(sub.cpu().numpy()) if want_prev else None
                del sub
                return pt, pv
            except RuntimeError as e:
                _gpu_fallback(logger, "sub-aperture", e)
                state["gpu"] = False
                state["spec"] = state["spec"].cpu().numpy()
        sub = np.fft.ifft2(state["spec"] * w[:, None])
        pt = _gather_patches(sub, rows, cols, h)
        pv = _decimate(sub) if want_prev else None
        del sub
        return pt, pv

    shifts, prev = [], []
    for i, a0 in enumerate(starts0):
        p_prev = None
        for j in range(chain + 1):
            a = a0 + j * u
            # ifftshift(spec*win) == ifftshift(spec) * ifftshift(win)
            w = np.fft.ifftshift(az_band(a, a + bw_px)).astype(np.float32)
            pt, pv = one_image(w, i == 0)
            if i == 0:
                prev.append(pv)
            if p_prev is not None:                    # coppia (master = p_prev, slave = pt)
                dy, dx, q = _subpixel_shifts(p_prev, pt)
                shifts.append((dy, dx, q, a - u))     # a - u = inizio del master
            p_prev = pt                               # lo slave diventa il master
        del p_prev

    del state["spec"]
    if USE_GPU:
        torch.cuda.empty_cache()
    _LAST_STRATEGY.clear()
    _LAST_STRATEGY.update(n_pulses=n_pulses, bw_px=bw_px, u=u, chain=chain,
                          starts0=list(starts0), B_cD=cfg.B_cD, n_pairs=len(shifts),
                          gpu=state["gpu"])
    logger.info(f"  generate {len(shifts)} coppie master/slave  "
                f"[{'GPU' if state['gpu'] else 'CPU'}]")
    return shifts, prev[:-1][:8], prev[1:][:8]


# =============================================================================
#  ONDE DA CATENA DI VETTORI DELTA (massa 1)
# =============================================================================

def _mass_response(F: np.ndarray, cfg: Config) -> np.ndarray:
    """Traiettoria di una massa 1 su molla smorzata spinta dalla forza F[k, pixel]
    (complessa: Re = asse x, Im = asse y): x'' = F - w0^2 x - lambda x', dt = 1 passo
    (Euler semi-implicito), w0 = 2*pi*osc_cycles/N, lambda = osc_damp*w0 [Eq. 18-20]."""
    n = F.shape[0]
    w0 = 2 * np.pi * cfg.osc_cycles / n
    lam = cfg.osc_damp * w0
    x = np.zeros(F.shape[1], dtype=np.complex128)
    v = np.zeros_like(x)
    X = np.zeros(F.shape, dtype=np.complex128)
    for k in range(n):
        v += F[k] - w0 ** 2 * x - lam * v
        x += v
        X[k] = x
    return X


def _spring_waves(D: np.ndarray, t: np.ndarray, omega: float = 1.0,
                  gamma: float = 0.0) -> np.ndarray:
    """Onde dei vettori delta su una massa 1 legata a una molla (omega rad per passo),
    smorzata di gamma per passo. Il vettore k, (a, b) = (Re, Im), tira la massa al passo
    k e la lascia: onda a*cos(omega (t - k))*exp(-gamma (t - k)) sull'asse Est-Ovest e
    lo stesso con b sul Nord-Sud, per t >= k. Le onde di tutti i vettori sono SOMMATE,
    una per asse. D: (N, pixel) complesso; t: istanti in passi (anche frazionari) ->
    (len(t), pixel) complesso, Re = Est-Ovest, Im = Nord-Sud."""
    lag = t[:, None] - np.arange(D.shape[0])[None, :]          # t - k
    W = np.where(lag >= 0, np.cos(omega * lag) * np.exp(-gamma * np.clip(lag, 0, None)), 0.0)
    return W @ D.real + 1j * (W @ D.imag)


def _lowpass(X: np.ndarray, cut: float, os_: int) -> np.ndarray:
    """Passa-basso Butterworth a 4 poli applicato avanti e indietro (filtfilt: fase zero,
    ordine effettivo 8, -6 dB al taglio), lungo l'onda (asse 0), separato su Est-Ovest
    (Re) e Nord-Sud (Im). cut = taglio [rad per passo]; os_ = campioni per passo."""
    from scipy.signal import butter, filtfilt
    b, a = butter(4, cut / (np.pi * os_))
    return filtfilt(b, a, X.real, axis=0) + 1j * filtfilt(b, a, X.imag, axis=0)


def vector_chain_waveforms(spectrum: np.ndarray, cfg: Config, logger: logging.Logger,
                           rows: np.ndarray, cols: np.ndarray, out_dir: str,
                           samples: Optional[Tuple[np.ndarray, np.ndarray]] = None) -> dict:
    """N_D vettori (Re, Im) per pixel -> vettori di spostamento -> onde di una massa 1.

    1) V[j, p] = valore complesso del pixel p nell'immagine j (finestra azimutale
       shiftata di 1 px per volta, j = 0..N_D-1). Il vettore 0 e' il master e il 1 lo
       slave: delta_0 = V[1] - V[0]; poi lo slave diventa il master: delta_1 = V[2] - V[1],
       ... (N_D - 1 vettori di spostamento, salvati; normalizzati per l'rms del pixel).
    2) ogni delta (a, b) tira una massa 1 legata a una molla, smorzata (l'onda di un
       tiro si dimezza ogni cfg.wave_half m di discesa), sugli assi Est-Ovest (Re) e
       Nord-Sud (Im); 3) le onde sono sommate: una per asse (_spring_waves);
       4) passa-basso sulle 2 onde (_lowpass, taglio cfg.wave_lpf).
    L'onda e' la proiezione in VERTICALE del moto del pixel e parte dalla superficie
    (t = 0): gli N delta sono distribuiti su tutta l'altezza pyr_height (un passo =
    pyr_height / N, il tiro k agisce al passo k) e la molla fa UN ciclo ogni N passi, cioe'
    sull'altezza; sotto, l'onda continua smorzata fino a pyr_height - z_bottom di
    profondita'. Calcolo a blocchi di pixel, 4 campioni per passo. samples = (nj, depth):
    per il pixel p servono nj[p] profondita' [m] sotto la superficie (concatenate pixel per
    pixel) a cui si restituisce l'ampiezza filtrata sqrt(E-O^2 + N-S^2) (centri dei voxel)."""
    n_pulses = spectrum.shape[0]
    n_vec = cfg.N_D
    bw_px = max(2, int(round((1 - cfg.B_DL_ratio) * n_pulses)))
    n_vec = min(n_vec, n_pulses - bw_px + 1)
    spec = np.fft.ifftshift(spectrum)
    V = np.zeros((n_vec, len(rows)), dtype=np.complex64)
    for j in range(n_vec):
        win = np.zeros(n_pulses, dtype=np.float32)
        win[j:j + bw_px] = np.hanning(bw_px) if cfg.window == "hann" else 1.0
        sub = np.fft.ifft2(spec * np.fft.ifftshift(win)[:, None])
        V[j] = sub[rows, cols]
        del sub
    D = V[1:] - V[:-1]                                  # delta master -> slave
    scale = np.sqrt((np.abs(V) ** 2).mean(axis=0)) + 1e-30
    Dn = D / scale[None, :]
    n, npx = Dn.shape
    os_ = 4                                             # campioni dell'onda per passo
    dz = cfg.pyr_height / n                             # n tiri distribuiti sull'altezza
    omega = 2 * np.pi / n                               # un ciclo della molla = altezza
    n_t = int(np.ceil((cfg.pyr_height - cfg.z_bottom) / dz))
    t = np.arange(n_t * os_ + 1) / os_                  # passi 0 .. n_t (t = 0 alla superficie)
    gamma = np.log(2) * dz / cfg.wave_half              # smorzamento per passo
    sub = np.unique(np.linspace(0, npx - 1, min(npx, 1500)).astype(int))   # pixel del grafico
    Xw = np.zeros((n_t + 1, npx), dtype=np.complex64)   # onde sommate, 1 campione per passo
    Xs = np.zeros((n_t + 1, sub.size), dtype=np.complex64)   # filtrate, pixel del grafico
    xf_sum = np.zeros(t.size, dtype=np.complex128)
    amp_mean = np.zeros(npx)
    ss_f = 0.0
    if samples is not None:
        nj, ts = samples[0], samples[1] / dz            # profondita' -> passi
        off = np.concatenate([[0], np.cumsum(nj)])
        amp_s = np.zeros(ts.size, dtype=np.float32)
    for p0 in range(0, npx, 2048):                      # blocchi di pixel: picco ~100 MB
        p1 = min(npx, p0 + 2048)
        X = _spring_waves(Dn[:, p0:p1], t, omega, gamma)
        Xf = _lowpass(X, cfg.wave_lpf, os_)
        A = np.abs(Xf)
        Xw[:, p0:p1] = X[::os_]
        k = (sub >= p0) & (sub < p1)
        Xs[:, k] = Xf[::os_][:, sub[k] - p0]
        xf_sum += Xf.sum(axis=1)
        amp_mean[p0:p1] = A.mean(axis=0)
        ss_f += float((A ** 2).sum())
        if samples is not None:                         # ampiezza ai centri dei voxel
            lo, hi = off[p0], off[p1]
            lc = np.repeat(np.arange(p1 - p0), nj[p0:p1])
            u = ts[lo:hi] * os_
            i0 = np.clip(u.astype(int), 0, t.size - 2)
            f = np.clip(u - i0, 0.0, 1.0)
            amp_s[lo:hi] = A[i0, lc] * (1 - f) + A[i0 + 1, lc] * f
    _, first = np.unique(amp_mean, return_index=True)   # con il layover piu' colonne hanno lo
    top = first[np.argsort(amp_mean[first])[::-1][:3]]  # stesso pixel: 3 onde distinte
    X_top = _spring_waves(Dn[:, top], t, omega, gamma)
    np.save(os.path.join(out_dir, "vec_V.npy"), V)
    np.save(os.path.join(out_dir, "vec_delta.npy"), D)
    np.save(os.path.join(out_dir, "vec_waves.npy"), Xw)
    logger.info(f"Onde da vettori: {n_vec} vettori/pixel, {n} delta, {npx} pixel; "
                f"{n} tiri sull'altezza {cfg.pyr_height:g} m (passo {dz:.2f} m), molla: un ciclo "
                f"= altezza, dimezza ogni {cfg.wave_half:g} m; passa-basso {cfg.wave_lpf:.3f} "
                f"rad/passo; |delta| rms {np.sqrt((np.abs(Dn) ** 2).mean()):.3f}, |onda| rms "
                f"{np.sqrt((np.abs(Xw) ** 2).mean()):.3f} -> filtrata "
                f"{np.sqrt(ss_f / (t.size * npx)):.3f}; dalla superficie fino a "
                f"{cfg.pyr_height - cfg.z_bottom:g} m di profondita'")
    return {"V": V, "D": D, "dz": dz, "d": t * dz, "d1": t[::os_] * dz, "Xs": Xs,
            "top": top, "X_top": X_top, "Xf_top": _lowpass(X_top, cfg.wave_lpf, os_),
            "Xf_mean": xf_sum / npx, "amp_s": amp_s if samples is not None else None}


def plot_vector_waveforms(res: dict, cfg: Config, out_path: str) -> None:
    """Onde est-ovest (Re) e nord-sud (Im) filtrate passa-basso in funzione della
    profondita' sotto la superficie del pixel (proiezione verticale, passo res["dz"]);
    sotto, per i 3 pixel piu' energetici, la somma delle onde (sottile) e la filtrata."""
    d, d1, Xs, dz = res["d"], res["d1"], res["Xs"], res["dz"]
    d_last = (res["D"].shape[0] - 1) * dz               # profondita' dell'ultimo tiro
    fig, ax = plt.subplots(2, 2, figsize=(14, 9), gridspec_kw={"height_ratios": [1.2, 1]})
    for c, (name, part) in enumerate([("Est-Ovest (Re)", np.real), ("Nord-Sud (Im)", np.imag)]):
        a = ax[0, c]
        img = part(Xs)
        vmax = np.percentile(np.abs(img), 99) + 1e-30
        im = a.imshow(img.T, aspect="auto", cmap="seismic", vmin=-vmax, vmax=vmax,
                      extent=[d1[0], d1[-1], img.shape[1], 0], interpolation="nearest")
        a.axvline(d_last, color="k", lw=0.8, ls="--")
        a.set_title(f"Onda {name} filtrata: {img.shape[1]} pixel x profondita'")
        a.set_ylabel("pixel (campione)")
        a.set_xlabel("profondita' sotto la superficie del pixel [m] (passo %.2f m)" % dz)
        fig.colorbar(im, ax=a, fraction=0.046)
        b = ax[1, c]
        for i, p in enumerate(res["top"]):
            b.plot(d, part(res["X_top"][:, i]), lw=0.7, alpha=0.5, color=f"C{i}")
            b.plot(d, part(res["Xf_top"][:, i]), lw=1.6, color=f"C{i}", label=f"pixel {p}")
        b.plot(d, part(res["Xf_mean"]), "k", lw=2, label="media filtrata")
        b.axvline(d_last, color="gray", lw=0.8, ls="--", label="ultimo tiro")
        b.grid(alpha=0.3); b.legend(fontsize=7)
        b.set_ylabel("ampiezza (norm.)"); b.set_xlabel("profondita' sotto la superficie del pixel [m]")
        b.set_title(f"{name}: 3 pixel piu' energetici (sottile = somma, spessa = filtrata)")
    plt.suptitle(f"Onde verticali di una massa 1 su molla tirata dai delta master/slave "
                 f"({res['D'].shape[0]} delta), dimezza ogni {cfg.wave_half:g} m\n"
                 f"tiri distribuiti sull'altezza {cfg.pyr_height:g} m (passo {dz:.2f} m), un ciclo "
                 f"della molla = altezza, passa-basso {cfg.wave_lpf:.2f} rad/passo", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


# =============================================================================
#  STEP 5 — MICRO-MOTION (pixel tracking, blocco 7)
# =============================================================================

def _subpixel_shifts(pm: np.ndarray, ps: np.ndarray):
    """Spostamento sub-pixel (dy, dx) di `pm` rispetto a `ps`, per ogni patch.
    Cross-correlazione FFT (media tolta, finestra di Hann) + picco parabolico.
    pm, ps: (N, P, P) ampiezze. dy = azimuth (axis 0), dx = range (axis 1)."""
    N, P, _ = pm.shape
    w = np.outer(np.hanning(P), np.hanning(P)).astype(np.float32)
    a = (pm - pm.mean(axis=(1, 2), keepdims=True)) * w
    b = (ps - ps.mean(axis=(1, 2), keepdims=True)) * w
    r = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real
    r = np.fft.fftshift(r, axes=(1, 2))

    flat = r.reshape(N, -1).argmax(axis=1)
    iy, ix = np.unravel_index(flat, (P, P))
    n = np.arange(N)
    # qualita' del tracking: picco della correlazione / norma delle due patch
    q = r[n, iy, ix] / (np.sqrt((a ** 2).sum(axis=(1, 2)) * (b ** 2).sum(axis=(1, 2)))
                        + 1e-30)

    def refine(m1, c0, p1):
        den = m1 - 2 * c0 + p1
        d = np.where(np.abs(den) > 1e-20, 0.5 * (m1 - p1) / np.where(den == 0, 1, den), 0.0)
        return np.clip(d, -0.5, 0.5)

    inner_y = (iy > 0) & (iy < P - 1)
    inner_x = (ix > 0) & (ix < P - 1)
    iy_c, ix_c = np.clip(iy, 1, P - 2), np.clip(ix, 1, P - 2)
    dy = refine(r[n, iy_c - 1, ix], r[n, iy_c, ix], r[n, iy_c + 1, ix]) * inner_y
    dx = refine(r[n, iy, ix_c - 1], r[n, iy, ix_c], r[n, iy, ix_c + 1]) * inner_x
    return iy - P // 2 + dy, ix - P // 2 + dx, q


def estimate_micromotion(
    shifts: list,
    cfg: Config,
    logger: logging.Logger,
    return_quality: bool = False,
):
    """Pixel tracking master/slave -> Y [k x N_pixels] [Eq. 20-21].
    y(k, p) = a + j*b: spostamento sub-pixel complesso (a = range, b = azimuth)
    del pixel p alla coppia k. Le coppie sono ordinate per posizione della finestra
    master (l'asse 'tempo' del metodo); le posizioni relative, normalizzate a [0,1],
    restano in cfg.k_frac e fissano le basi B_perp della matrice di steering
    (campionamento non uniforme quando la catena e' attiva)."""
    logger.info("Step 5: Micro-Motion (Pixel Tracking)...")

    n_sub = len(shifts)
    n_pix = shifts[0][0].shape[0]
    order = np.argsort([s_[3] for s_ in shifts], kind="stable")
    Y = np.zeros((n_sub, n_pix), dtype=np.complex64)
    Q = np.zeros(n_pix, dtype=np.float64)
    starts = np.zeros(n_sub, dtype=np.float64)
    for r_, k_ in enumerate(order):
        dy, dx, q, st = shifts[k_]
        Y[r_, :] = dx + 1j * dy
        Q += q / n_sub
        starts[r_] = st
    if cfg.osc:
        # pixel = massa 1 su molla (Eq. 18-20): x'' = F - w0^2 x - lambda x', con la forza
        # F_k = dx + j*dy (parte reale = range, immaginaria = azimuth) e dt = 1 passo
        # (Euler semi-implicito). Y(k) = x(k): l'onda nasce dall'oscillatore, non dalla IFFT.
        Yo = _mass_response(Y, cfg)
        logger.info(f"  Oscillatore (massa 1): {cfg.osc_cycles:g} cicli su {n_sub}, "
                    f"lambda={cfg.osc_damp:g}*w0; "
                    f"rms F {np.sqrt((np.abs(Y) ** 2).mean()):.4f} -> "
                    f"rms x {np.sqrt((np.abs(Yo) ** 2).mean()):.4f}")
        Y = Yo.astype(np.complex64)
    if cfg.lpf_sigma > 0 and n_sub > 2:
        # passa-basso (gaussiano, fase zero) lungo k sulle parti reale e immaginaria
        # di Y: toglie il rumore del tracking tra coppie vicine
        from scipy.ndimage import gaussian_filter1d
        rms0 = np.sqrt((np.abs(Y) ** 2).mean())
        Y = (gaussian_filter1d(Y.real, cfg.lpf_sigma, axis=0, mode="nearest")
             + 1j * gaussian_filter1d(Y.imag, cfg.lpf_sigma, axis=0, mode="nearest")
             ).astype(np.complex64)
        logger.info(f"  LPF su k: sigma={cfg.lpf_sigma:g} campioni, rms "
                    f"{rms0:.4f} -> {np.sqrt((np.abs(Y) ** 2).mean()):.4f} px")
    span = starts.max() - starts.min()
    cfg.k_frac = ((starts - starts.min()) / span if span > 0
                  else np.linspace(0, 1, n_sub))

    logger.info(f"  Y shape: {Y.shape}  |Y| medio = {np.abs(Y).mean():.3f} px")
    return (Y, Q) if return_quality else Y


# =============================================================================
#  STEP 6 — TOMOGRAPHIC FOCUSING
# =============================================================================

def build_steering_matrix(K_z: np.ndarray, z: np.ndarray,
                          two_pi: bool = False) -> np.ndarray:
    """A[i, z] = exp(j * c * Kz_i * z). Kz = 4*pi*B_perp/(lambda*r*sin(theta)) e'
    GIA' in rad/m: c = 1. Con c = 2*pi (lettura letterale dell'Eq. 22) la fase
    per metro e' 2*pi volte piu' alta e la profondita' risulta sottocampionata."""
    c = 2 * np.pi if two_pi else 1.0
    A = np.zeros((len(K_z), len(z)), dtype=np.complex128)
    for i, k in enumerate(K_z):
        A[i, :] = np.exp(1j * c * k * z)
    return A


def _tomo_pinv(n_sub: int, cfg: Config) -> np.ndarray:
    """Pseudo-inversa A^+ (F x n_sub) della matrice di steering, complex64."""
    frac = getattr(cfg, "k_frac", None)
    if frac is not None and len(frac) == n_sub:
        B_perp = np.asarray(frac, dtype=float) * cfg.A_orbit      # posizioni reali
    else:
        B_perp = np.linspace(0, cfg.A_orbit, n_sub)
    theta = np.deg2rad(45.0)
    K_z = (4 * np.pi * B_perp) / (cfg.lambda_sound * cfg.r_ref * np.sin(theta))
    A = build_steering_matrix(K_z, cfg.z, cfg.steer_2pi)
    # rcond: scarta i valori singolari < 1e-3 del massimo (con molti passi la
    # matrice e' mal condizionata e la pinv piena amplificherebbe il rumore ~1e9)
    return np.linalg.pinv(A, rcond=1e-3).astype(np.complex64)


def tomographic_focusing(
    Y: np.ndarray,
    cfg: Config,
    logger: logging.Logger,
) -> np.ndarray:
    """h(z) = A† · Y  [Eq. 23-24]."""
    logger.info("Step 6: Tomographic Focusing...")

    n_sub, n_pix = Y.shape
    z = cfg.z
    F = len(z)

    A_pinv = _tomo_pinv(n_sub, cfg)
    Y64 = Y.astype(np.complex64)             # complex64: meta' memoria
    h = None
    if USE_GPU:
        try:
            h_t = torch.from_numpy(A_pinv).cuda() @ torch.from_numpy(Y64).cuda()
            h = h_t.cpu().numpy()
            del h_t
            torch.cuda.empty_cache()
        except RuntimeError as e:
            _gpu_fallback(logger, "tomographic_focusing", e)
    if h is None:
        h = A_pinv @ Y64

    logger.info(f"  h shape: {h.shape}  (F={F}, N_pix={n_pix})")
    return h


# =============================================================================
#  STEP 7 — MAGNITUDE
# =============================================================================

def generate_tomographic_magnitude(
    h: np.ndarray,
    cfg: Config,
    logger: logging.Logger,
) -> Tuple[np.ndarray, np.ndarray]:
    logger.info("Step 7: Magnitude Image...")
    magnitude = np.abs(h)                      # float32
    mn, mx = float(magnitude.min()), float(magnitude.max())
    magnitude /= np.float32(mx if mx > 0 else 1.0)    # in-place -> norm
    magnitude_db = np.maximum(magnitude, np.float32(1e-12))
    np.log10(magnitude_db, out=magnitude_db)
    magnitude_db *= np.float32(20)
    logger.info(f"  Magnitude range: [{mn:.3e}, {mx:.3e}]")
    return magnitude, magnitude_db


# =============================================================================
#  STEP 8 — VISUALIZZAZIONE
# =============================================================================

def plot_tomographic_result(
    tomogram: np.ndarray,
    cfg: Config,
    out_path: str,
    title: str = "3D Tomographic Imaging (magnitude)",
    log_scale: Optional[bool] = None,
) -> None:
    """
    Tomogramma in stile Figure 29-38b:
      X: Azimuth-depth (pixel)
      Y: Range-depth (pixel)
      Colorbar: Magnitude (jet)
    """
    if log_scale is None:
        log_scale = cfg.log_scale

    fig, ax = plt.subplots(figsize=(10, 6), dpi=cfg.dpi)

    full_h, full_w = tomogram.shape
    small = _decimate(tomogram, 2000)          # block-max per il plot
    data = small if not log_scale else 20 * np.log10(small + 1e-12)
    vmin, vmax = (np.nanmin(data), np.nanmax(data))

    im = ax.imshow(
        data, cmap=cfg.cmap, aspect="auto", origin="upper",
        vmin=vmin, vmax=vmax, interpolation="nearest",
        extent=(0, full_w, full_h, 0),         # assi in pixel reali
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Magnitude" if not log_scale else "Magnitude (dB)",
                   fontsize=11)

    ax.set_xlabel("Tomographic line (pixel index)", fontsize=12)
    ax.set_ylabel("Range-depth (pixel)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.xaxis.set_major_locator(MultipleLocator(max(50, tomogram.shape[1] // 10)))
    ax.yaxis.set_major_locator(MultipleLocator(max(50, tomogram.shape[0] // 10)))
    ax.grid(True, color="white", alpha=0.15, linewidth=0.5)

    info = (
        f"N_D = {cfg.N_D} | δ_T = {cfg.delta_T:.2f} m | "
        f"λ = {cfg.lambda_sound:.3f} m | Pol = {cfg.polarization}"
    )
    ax.text(0.01, -0.12, info, transform=ax.transAxes,
            fontsize=9, color="gray")

    plt.tight_layout()
    plt.savefig(out_path, dpi=cfg.dpi, bbox_inches="tight")
    plt.close(fig)


def plot_processing_overview(
    slc_prev: np.ndarray,
    spec_prev: np.ndarray,
    mag: np.ndarray,
    cfg: Config,
    out_path: str,
) -> None:
    """slc_prev / spec_prev: |SLC| e |spettro| gia' ridotti con _decimate;
    mag: magnitudine tomografica normalizzata."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=cfg.dpi)

    axes[0, 0].imshow(slc_prev, cmap="gray", aspect="auto")
    axes[0, 0].set_title("(a) SLC SAR image (magnitude)")
    axes[0, 0].set_xlabel("Range (pixel)")
    axes[0, 0].set_ylabel("Azimuth (pixel)")

    spec_db = 20 * np.log10(spec_prev + 1e-12)
    axes[0, 1].imshow(spec_db, cmap="viridis", aspect="auto")
    axes[0, 1].set_title("(b) Spettro 2D (dB)")
    axes[0, 1].set_xlabel("Range frequency (pixel)")
    axes[0, 1].set_ylabel("Azimuth frequency (pixel)")

    im = axes[1, 0].imshow(_decimate(mag, 2000), cmap=cfg.cmap, aspect="auto")
    axes[1, 0].set_title("(c) 3D Tomographic Imaging (magnitude)")
    axes[1, 0].set_xlabel("Tomographic line (pixel index)")
    axes[1, 0].set_ylabel("Range-depth (pixel)")
    fig.colorbar(im, ax=axes[1, 0], fraction=0.046, pad=0.04)

    prof = np.mean(mag, axis=1, dtype=np.float64)
    z_axis = np.linspace(cfg.z_min, cfg.z_MAX, len(prof))
    axes[1, 1].plot(prof, z_axis, "b-", linewidth=1.5)
    axes[1, 1].fill_betweenx(z_axis, 0, prof, alpha=0.3, color="blue")
    axes[1, 1].set_title("(d) Profilo verticale medio (magnitude)")
    axes[1, 1].set_xlabel("Magnitude (norm.)")
    axes[1, 1].set_ylabel("Quota z (m)")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=cfg.dpi, bbox_inches="tight")
    plt.close(fig)


def plot_subaparture_strategy(
    master_subaps: List[np.ndarray],
    slave_subaps: List[np.ndarray],
    cfg: Config,
    out_path: str,
) -> None:
    """Strategia Doppler a catena (Fig. 7 modificata).
    In alto: finestre azimutali della PRIMA catena (ogni immagine e' lo slave della
    coppia precedente e il master della successiva) e posizioni di partenza delle
    altre catene. In basso: anteprime |immagine| della prima catena
    (master_j -> slave_j = master_{j+1})."""
    info = _LAST_STRATEGY
    n = min(len(master_subaps), 8)
    fig = plt.figure(figsize=(max(14, 2.6 * n), 9), dpi=min(cfg.dpi, 150))
    gs = fig.add_gridspec(3, max(n, 1), height_ratios=[1.5, 1, 1])

    ax = fig.add_subplot(gs[0, :])
    if info:
        npx, bw, u, chain = info["n_pulses"], info["bw_px"], info["u"], info["chain"]
        kHz = lambda px: px / npx * info["B_cD"] / 1e3
        a0s = info["starts0"]
        cmap = plt.cm.viridis
        for j in range(chain + 1):                    # prima catena
            x0 = a0s[0] + j * u
            ax.barh(j, kHz(bw), left=kHz(x0), height=0.7,
                    color=cmap(j / max(1, chain)), alpha=0.85, edgecolor="k")
            lab = ("M1" if j == 0 else (f"S{j} = M{j + 1}" if j < chain else f"S{chain}"))
            ax.text(kHz(x0) + kHz(bw) / 2, j, lab, ha="center", va="center", fontsize=8,
                    color="white" if j < chain / 2 else "black")
        for i in range(1, len(a0s)):                  # altre catene: solo l'origine
            ax.plot([kHz(a0s[i])] * 2, [-0.5, chain + 0.5], color="gray", lw=0.6, alpha=0.6)
        ax.set_xlim(0, kHz(npx)); ax.set_ylim(-0.7, chain + 0.7)
        ax.set_xlabel("frequenza Doppler [kHz]  (banda B_cD; larghezza finestra = B_cD - B_DL)")
        ax.set_ylabel("immagine nella catena")
        if len(a0s) == 1:
            ax.set_title(f"B_DL diviso in N_D = {chain} passi | {info['n_pairs']} coppie "
                         f"master/slave | B_shift = u = B_DL/N_D = {kHz(u) * 1e3:.0f} Hz",
                         fontsize=10)
        else:
            ax.set_title(f"N_D = {len(a0s)} posizioni x {chain} coppie = {info['n_pairs']} misure | "
                         f"B_shift = u = {kHz(u) * 1e3:.0f} Hz | linee grigie: origine delle altre "
                         f"{len(a0s) - 1} catene", fontsize=10)
        ax.invert_yaxis()
    else:
        ax.axis("off")

    for i in range(n):
        a1 = fig.add_subplot(gs[1, i]); a2 = fig.add_subplot(gs[2, i])
        a1.imshow(master_subaps[i], cmap="gray", aspect="auto")
        a1.set_title(f"Master {i + 1}", fontsize=9); a1.axis("off")
        a2.imshow(slave_subaps[i], cmap="gray", aspect="auto")
        a2.set_title(f"Slave {i + 1} (= Master {i + 2})", fontsize=9); a2.axis("off")
    plt.suptitle("Doppler Sub-Aperture Strategy a catena (slave -> master successivo)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=min(cfg.dpi, 150), bbox_inches="tight")
    plt.close(fig)


def plot_profile_validation(mag: np.ndarray, cfg: Config, out_path: str) -> None:
    """mag: magnitudine tomografica (gia' |h|, float32)."""
    try:
        from scipy.ndimage import uniform_filter1d
    except ImportError:
        uniform_filter1d = None

    prof = np.mean(mag, axis=1, dtype=np.float64)
    prof_n = prof / (prof.max() + 1e-12)
    z_axis = np.linspace(cfg.z_min, cfg.z_MAX, len(prof))

    if uniform_filter1d is not None:
        prof_filt = uniform_filter1d(prof_n, size=5)
    else:
        prof_filt = prof_n

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), dpi=cfg.dpi)

    axes[0].plot(z_axis, prof_n, color="steelblue", alpha=0.6,
                 label="Unfiltered")
    axes[0].plot(z_axis, prof_filt, color="black", linewidth=2,
                 label="Filtered")
    axes[0].set_xlabel("Quota z (m)")
    axes[0].set_ylabel("Magnitude (norm.)")
    axes[0].set_title("(a) Profilo tomografico verticale")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    err = prof_n - prof_filt
    axes[1].plot(z_axis, err, color="black", linewidth=1.5)
    axes[1].fill_between(z_axis, 0, err, alpha=0.3, color="gray")
    axes[1].set_xlabel("Quota z (m)")
    axes[1].set_ylabel("Errore (norm.)")
    axes[1].set_title("(b) Errore unfiltered vs filtered")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=cfg.dpi, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
#  NULL TEST — la stessa catena su finestre senza piramide
# =============================================================================
#  Il paper non pubblica ne' un null ne' una soglia: le strutture sono lette a
#  occhio. Qui la catena (pixel puri -> sub-aperture -> tracking -> tomografia)
#  e' eseguita su TUTTA la scena; poi la regione della piramide (ROI) e' confrontata
#  con molte finestre di pari dimensione fuori dalla ROI. Se la ROI non si
#  distingue dal null, le mappe non contengono segnale oltre il rumore.

def pure_pixel_map(slc: np.ndarray, cfg: Config, max_pix: int = 80000):
    """Tutti i pixel 'puri' (massimo locale 3x3 con >= pure_ratio dell'energia
    3x3) lontani dai bordi almeno mezza patch. Ritorna (rows, cols, amp)."""
    h = cfg.patch // 2
    p = np.abs(slc).astype(np.float64) ** 2
    c = p[1:-1, 1:-1]
    e3 = np.zeros_like(c)
    ismax = np.ones(c.shape, dtype=bool)
    n_az, n_rg = p.shape
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            nb = p[1 + dy:n_az - 1 + dy, 1 + dx:n_rg - 1 + dx]
            e3 += nb
            if dy or dx:
                ismax &= c >= nb
    mask = ismax & (c / (e3 + 1e-30) >= cfg.pure_ratio)
    r, k = np.nonzero(mask)
    r, k = r + 1, k + 1
    ok = (r >= h) & (r <= n_az - h) & (k >= h) & (k <= n_rg - h)
    r, k = r[ok], k[ok]
    amp = np.sqrt(p[r, k])
    if r.size > max_pix:
        keep = np.argsort(amp)[::-1][:max_pix]
        r, k, amp = r[keep], k[keep], amp[keep]
    return r.astype(int), k.astype(int), amp


def _tomo_peak_ratio(Y: np.ndarray, cfg: Config, chunk: int = 4000) -> np.ndarray:
    """Per ogni pixel: log10( picco / mediana ) di |h(z)| lungo z."""
    A_pinv = _tomo_pinv(Y.shape[0], cfg)
    out = np.empty(Y.shape[1], dtype=np.float64)
    dev = "cuda" if USE_GPU else "cpu"
    Ap = torch.from_numpy(A_pinv).to(dev)
    for i in range(0, Y.shape[1], chunk):
        Yc = torch.from_numpy(np.ascontiguousarray(Y[:, i:i + chunk])
                              .astype(np.complex64)).to(dev)
        m = (Ap @ Yc).abs()
        pk = m.max(dim=0).values
        md = m.median(dim=0).values
        out[i:i + chunk] = torch.log10(pk / (md + 1e-30) + 1e-30).double().cpu().numpy()
    return out


def _neighbor_coherence(Y: np.ndarray, rows: np.ndarray, cols: np.ndarray,
                        dmin: int, dmax: int) -> np.ndarray:
    """Coerenza |<Yi, Yj>|/(|Yi||Yj|) con il pixel successivo sulla stessa riga a
    distanza [dmin, dmax] (oltre la patch: niente sovrapposizione). NaN se manca."""
    order = np.lexsort((cols, rows))
    key = rows[order].astype(np.int64) * 1_000_000 + cols[order]
    j = np.searchsorted(key, key + dmin)
    valid = j < key.size
    jj = np.minimum(j, key.size - 1)
    valid &= (rows[order][jj] == rows[order]) & (cols[order][jj] - cols[order] <= dmax)
    Ys = Y[:, order].astype(np.complex128)
    num = np.abs(np.sum(Ys * np.conj(Ys[:, jj]), axis=0))
    den = np.linalg.norm(Ys, axis=0) * np.linalg.norm(Ys[:, jj], axis=0) + 1e-30
    coh_sorted = np.where(valid, num / den, np.nan)
    out = np.empty_like(coh_sorted)
    out[order] = coh_sorted
    return out


def _stratified_z(stat: np.ndarray, covs, ref_mask: np.ndarray,
                  n_bins: int = 6, min_ref: int = 50):
    """z-score di `stat` dentro strati definiti da una o piu' covariate (bin a
    quantili), con media/sd ricavate SOLO dai pixel di riferimento (fuori ROI).
    Toglie l'effetto di SNR e qualita' del tracking, che cambiano la varianza
    dello stimatore indipendentemente dalla vibrazione."""
    if isinstance(covs, np.ndarray):
        covs = [covs]
    idx = np.zeros(stat.shape, dtype=np.int64)
    for k, c in enumerate(covs):
        edges = np.quantile(c, np.linspace(0, 1, n_bins + 1))
        bk = np.clip(np.searchsorted(edges, c, side="right") - 1, 0, n_bins - 1)
        idx += bk * (n_bins ** k)
    z = np.full(stat.shape, np.nan)
    n_dropped = 0
    for k in np.unique(idx):
        inb = idx == k
        ref = inb & ref_mask & np.isfinite(stat)
        if ref.sum() < min_ref:
            n_dropped += int((inb & np.isfinite(stat)).sum())
            continue
        mu, sd = stat[ref].mean(), stat[ref].std() + 1e-12
        z[inb] = (stat[inb] - mu) / sd
    return z, n_dropped


PYRAMIDS = {          # (lat, lon, lato base [m], altezza attuale [m])
    "khufu": (29.979167, 31.134167, 230.3, 138.5),
    "khafre": (29.976128, 31.130589, 215.3, 136.4),     # fit spigoli SLC 2 date (2026-10-08)
    "menkaure": (29.972402, 31.128073, 102.2, 61.0),    # fit spigoli SLC 2 date (2026-10-08)
}
def load_known_structures(logger: Optional[logging.Logger] = None) -> dict:
    """Parti note per piramide, da known_structures.json (accanto al programma): estratte
    da output_vox_20241004/khufu_3d.html e khafre_3d.html. Ogni struttura: nome, categoria,
    colore, fonte, incertezza, schematica, nota, c = centro, seg = segmenti del wireframe
    (x1, y1, z1, x2, y2, z2, ...; x Est, y Nord, z sopra la base, origine = centro base).
    Riferimenti da letteratura, NON derivati dai dati SAR."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "known_structures.json")
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError) as e:
        if logger:
            logger.warning(f"  parti note non caricate ({path}): {e}")
        return {}


_WGS_A, _WGS_F = 6378137.0, 1 / 298.257223563
_WGS_E2 = _WGS_F * (2 - _WGS_F)


def _ecf2geo(p):
    x, y, z = p
    lon = np.arctan2(y, x)
    r = np.hypot(x, y)
    lat = np.arctan2(z, r * (1 - _WGS_E2))
    for _ in range(8):
        n = _WGS_A / np.sqrt(1 - _WGS_E2 * np.sin(lat) ** 2)
        h = r / np.cos(lat) - n
        lat = np.arctan2(z, r * (1 - _WGS_E2 * n / (n + h)))
    n = _WGS_A / np.sqrt(1 - _WGS_E2 * np.sin(lat) ** 2)
    return lat, lon, r / np.cos(lat) - n


def _geo2ecf(lat, lon, h):
    n = _WGS_A / np.sqrt(1 - _WGS_E2 * np.sin(lat) ** 2)
    return np.array([(n + h) * np.cos(lat) * np.cos(lon),
                     (n + h) * np.cos(lat) * np.sin(lon),
                     (n * (1 - _WGS_E2) + h) * np.sin(lat)])


def pyramid_faces_px(geo: dict, name: str) -> dict:
    """Poligoni (riga, colonna) della base e delle 4 facce di una piramide nota,
    nell'immagine di backprojection (piano al suolo, con LAYOVER dell'apice:
    un punto a quota z compare spostato di z/tan(inc) verso il sensore).
    Lati allineati ai punti cardinali (come Giza). Le facce N,E,S,O sono i
    triangoli (spigolo di base, apice) e tassellano la base."""
    la, lo, side, H = PYRAMIDS[name]
    srp, u_rg, u_az, inc = geo["srp"], geo["u_rg"], geo["u_az"], geo["inc"]
    lat0, lon0, h0 = _ecf2geo(srp)
    east = np.array([-np.sin(lon0), np.cos(lon0), 0.0])
    north = np.array([-np.sin(lat0) * np.cos(lon0), -np.sin(lat0) * np.sin(lon0),
                      np.cos(lat0)])
    d = _geo2ecf(np.radians(la), np.radians(lo), h0) - srp   # stessa quota dell'SRP
    xc, yc = d @ u_rg, d @ u_az
    half = side / 2

    def px(e_off, n_off, z=0.0):
        v = e_off * east + n_off * north
        x = xc + v @ u_rg - z / np.tan(inc)
        y = yc + v @ u_az
        return (geo["n_az"] // 2 + (y - geo["cy"]) / geo["da"],
                geo["n_rg"] // 2 + (x - geo["cx"]) / geo["dr"])

    c = {"SW": px(-half, -half), "SE": px(half, -half),
         "NE": px(half, half), "NW": px(-half, half)}
    ap = px(0.0, 0.0, H)
    ctr = px(0.0, 0.0, 0.0)                            # centro della base, a terra
    return {"base": [c["SW"], c["SE"], c["NE"], c["NW"]],
            # layover: triangoli (spigolo di base, apice spostato); con l'apice vicino
            # a uno spigolo i grandi triangoli E/S ricoprono N/O
            "N": [c["NW"], c["NE"], ap], "E": [c["NE"], c["SE"], ap],
            "S": [c["SE"], c["SW"], ap], "W": [c["SW"], c["NW"], ap],
            # diagonali: le due diagonali della base dividono il quadrato in 4
            # triangoli (centro, spigolo, spigolo), uno per lato, senza sovrapposizioni
            "N_d": [c["NW"], c["NE"], ctr], "E_d": [c["NE"], c["SE"], ctr],
            "S_d": [c["SE"], c["SW"], ctr], "W_d": [c["SW"], c["NW"], ctr],
            "apex": ap, "center": ctr, "corners": c}


def polygon_mask(shape: Tuple[int, int], poly) -> np.ndarray:
    """Maschera booleana (righe, colonne) dei pixel dentro il poligono."""
    from matplotlib.path import Path
    pts = np.array([(p[1], p[0]) for p in poly])             # (x=col, y=riga)
    r0 = int(max(0, np.floor(pts[:, 1].min()))); r1 = int(min(shape[0], np.ceil(pts[:, 1].max()) + 1))
    c0 = int(max(0, np.floor(pts[:, 0].min()))); c1 = int(min(shape[1], np.ceil(pts[:, 0].max()) + 1))
    m = np.zeros(shape, dtype=bool)
    if r1 <= r0 or c1 <= c0:
        return m
    rr, cc = np.mgrid[r0:r1, c0:c1]
    inside = Path(pts).contains_points(np.column_stack([cc.ravel(), rr.ravel()]))
    m[r0:r1, c0:c1] = inside.reshape(rr.shape)
    return m


def _dilate(mask: np.ndarray, k: int) -> np.ndarray:
    """Dilatazione con quadrato (2k+1)x(2k+1), senza scipy."""
    out = mask.copy()
    for d in range(1, k + 1):
        out[d:, :] |= mask[:-d, :]
        out[:-d, :] |= mask[d:, :]
    tmp = out.copy()
    for d in range(1, k + 1):
        out[:, d:] |= tmp[:, :-d]
        out[:, :-d] |= tmp[:, d:]
    return out


def null_test_masks(slc: np.ndarray, cfg: Config, logger: logging.Logger,
                    masks: dict, exclude: np.ndarray, n_null: int = 300,
                    seed: int = 0, inject_px: float = 0.0,
                    tag: str = "faces", aspect: float = 1.0) -> dict:
    """Null test per regioni di forma qualsiasi.

    masks: {nome: maschera booleana}. exclude: zona da non usare ne' come
    riferimento ne' per le finestre null (le piramidi, dilatate). Per ogni
    regione: statistica = media degli z-score (stratificati per ampiezza, con
    riferimento = pixel fuori da `exclude`) dei pixel puri della regione; il null
    sono N traslazioni della STESSA forma, interamente fuori da `exclude`.
    aspect = da / dr (lato del pixel in azimuth / in range) per la mappa delle regioni."""
    logger.info(f"NULL TEST ({tag}): catena su tutta la scena, {len(masks)} regioni...")
    n_az, n_rg = slc.shape
    rows, cols, amp = pure_pixel_map(slc, cfg, max_pix=100000)
    logger.info(f"  pixel puri in scena: {rows.size}")

    spectrum = dft2(slc, logger)
    patches, _, _ = generate_doppler_subapartures(spectrum, cfg, logger, rows, cols)
    del spectrum
    Y, Q = estimate_micromotion(patches, cfg, logger, return_quality=True)
    del patches

    if inject_px > 0:
        ph = np.exp(1j * 2 * np.pi * 0.2 * np.arange(Y.shape[0]))
        sel = np.zeros(rows.size, dtype=bool)
        for m_ in masks.values():
            sel |= m_[rows, cols]
        Y[:, sel] += (inject_px * ph[:, None]).astype(Y.dtype)
        logger.info(f"  CONTROLLO POSITIVO: iniettati {inject_px:.3f} px coerenti "
                    f"in tutte le regioni")

    stats = {
        "|Y| rms [log10 px]": np.log10(np.sqrt(np.mean(np.abs(Y) ** 2, axis=0)) + 1e-12),
        "picco/mediana |h(z)| [log10]": _tomo_peak_ratio(Y, cfg),
        "coerenza tra pixel vicini": _neighbor_coherence(Y, rows, cols,
                                                         cfg.patch + 1, 3 * cfg.patch),
    }
    ref_mask = ~exclude[rows, cols]
    logger.info(f"  pixel di riferimento (fuori dalle piramidi): {int(ref_mask.sum())}")
    covs = [np.log10(amp + 1e-30), Q]            # ampiezza e qualita' del tracking
    zs = {name: _stratified_z(st, covs, ref_mask)[0] for name, st in stats.items()}

    rng = np.random.default_rng(seed)
    res = {"regions": {}}
    lines = []
    for rname, M in masks.items():
        ii, jj = np.nonzero(M)
        a0, a1, b0, b1 = ii.min(), ii.max() + 1, jj.min(), jj.max() + 1
        Mc = M[a0:a1, b0:b1]
        hh, ww = Mc.shape
        in_roi = M[rows, cols]
        offs = []
        for _ in range(30000):
            oa = int(rng.integers(0, n_az - hh + 1)); ob = int(rng.integers(0, n_rg - ww + 1))
            if not exclude[oa:oa + hh, ob:ob + ww][Mc].any():
                offs.append((oa, ob))
            if len(offs) >= n_null:
                break
        n_eff = max(1.0, float((~exclude).sum()) / max(1, int(M.sum())))
        entry = {"n_pix": int(in_roi.sum()), "n_null_win": len(offs), "n_eff": n_eff,
                 "stats": {}}
        for sname, z in zs.items():
            w_roi = float(np.nanmean(z[in_roi])) if np.isfinite(z[in_roi]).sum() >= 10 else np.nan
            w_null = []
            for oa, ob in offs:
                m = (rows >= oa) & (rows < oa + hh) & (cols >= ob) & (cols < ob + ww)
                m &= Mc[np.clip(rows - oa, 0, hh - 1), np.clip(cols - ob, 0, ww - 1)]
                if np.isfinite(z[m]).sum() >= 10:
                    w_null.append(float(np.nanmean(z[m])))
            w_null = np.array(w_null)
            if len(w_null) < 10 or not np.isfinite(w_roi):
                entry["stats"][sname] = dict(w_roi=w_roi, z=np.nan, p_emp=np.nan,
                                             w_null=w_null, w_null_mean=np.nan, w_null_sd=np.nan)
                lines.append(f"  {rname:8s} {sname:30s} dati insufficienti "
                             f"(pixel ROI={entry['n_pix']}, finestre null={len(w_null)})")
                continue
            zz = (w_roi - w_null.mean()) / (w_null.std() + 1e-12)
            p = (1 + np.sum(np.abs(w_null - w_null.mean()) >= abs(w_roi - w_null.mean()))) \
                / (len(w_null) + 1)
            entry["stats"][sname] = dict(w_roi=w_roi, w_null_mean=float(w_null.mean()),
                                         w_null_sd=float(w_null.std()), z=float(zz),
                                         p_emp=float(p), w_null=w_null)
            lines.append(f"  {rname:8s} {sname:30s} ROI={w_roi:+.3f} null={w_null.mean():+.3f}"
                         f"+-{w_null.std():.3f} z={zz:+.2f} p_emp={p:.3f}"
                         f"  [pix={entry['n_pix']}, finestre={len(w_null)}, aree~{n_eff:.1f}]")
        res["regions"][rname] = entry
    for ln in lines:
        logger.info(ln)

    # --- figure ---
    names = list(masks)
    fig, axes = plt.subplots(len(names), 3, figsize=(15, 3.2 * len(names)),
                             dpi=min(cfg.dpi, 110), squeeze=False)
    for r_, rname in enumerate(names):
        for c_, (sname, st) in enumerate(res["regions"][rname]["stats"].items()):
            ax = axes[r_, c_]
            if len(st["w_null"]):
                ax.hist(st["w_null"], bins=25, color="steelblue", alpha=0.8)
            if np.isfinite(st["w_roi"]):
                ax.axvline(st["w_roi"], color="red", lw=2,
                           label=f"{rname}  z={st['z']:+.1f}")
            ax.set_title(sname, fontsize=9)
            ax.legend(fontsize=8)
    plt.tight_layout()
    out_png = os.path.join(cfg.out_dir, f"null_{tag}.png")
    plt.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

    # mappa delle regioni sull'SLC
    a = np.abs(slc)
    fig, ax = plt.subplots(figsize=(14, 9), dpi=100)
    ax.imshow(20 * np.log10(a / np.percentile(a, 99.9) + 1e-6), cmap="gray",
              vmin=-45, vmax=0, aspect=aspect, extent=(0, n_rg, n_az, 0))
    cols_ = plt.cm.tab10(np.linspace(0, 1, 10))
    for k, (rname, M) in enumerate(masks.items()):
        ov = np.zeros(M.shape + (4,))
        ov[M] = (*cols_[k % 10][:3], 0.35)
        ax.imshow(ov, aspect=aspect, extent=(0, n_rg, n_az, 0))
        ii, jj = np.nonzero(M)
        ax.text(jj.mean(), ii.mean(), rname, color="white", fontsize=12, ha="center")
    plt.tight_layout()
    plt.savefig(os.path.join(cfg.out_dir, f"null_{tag}_regioni.png"), bbox_inches="tight")
    plt.close(fig)

    with open(os.path.join(cfg.out_dir, f"null_{tag}_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    logger.info(f"  → {out_png}")
    return res


def null_test(slc: np.ndarray, cfg: Config, logger: logging.Logger,
              roi: Tuple[int, int, int, int], n_null: int = 300,
              seed: int = 0, inject_px: float = 0.0, aspect: float = 1.0) -> dict:
    """Null test su una ROI rettangolare (az0, az1, rg0, rg1) in pixel."""
    az0, az1, rg0, rg1 = [int(v) for v in roi]
    M = np.zeros(slc.shape, dtype=bool)
    M[az0:az1, rg0:rg1] = True
    return null_test_masks(slc, cfg, logger, {"ROI": M}, M, n_null, seed,
                           inject_px, tag="roi", aspect=aspect)


def faces_test(slc: np.ndarray, cfg: Config, logger: logging.Logger, geo: dict,
               pyramids: List[str], n_null: int = 300, inject_px: float = 0.0,
               mode: str = "diag") -> dict:
    """Null test sulle 4 facce (N,E,S,O) delle piramidi indicate.
    mode "diag": 4 triangoli delimitati dalle diagonali della base (nessuna
    sovrapposizione); mode "layover": triangoli base-apice con l'apice spostato."""
    shape = slc.shape
    masks, excl = {}, np.zeros(shape, dtype=bool)
    for name in PYRAMIDS:                              # tutte le piramidi note sono escluse
        pf = pyramid_faces_px(geo, name)
        excl |= polygon_mask(shape, pf["base"])
    excl = _dilate(excl, 15)
    for name in pyramids:
        pf = pyramid_faces_px(geo, name)
        for f in ("N", "E", "S", "W"):
            m = polygon_mask(shape, pf[f + ("_d" if mode == "diag" else "")])
            if m.sum() >= 200:
                masks[f"{name[:4]}-{f}"] = m
            else:
                logger.warning(f"  faccia {name}-{f} fuori scena o troppo piccola")
    if not masks:
        raise RuntimeError("nessuna faccia dentro la scena (controlla --center_x/--center_y)")
    return null_test_masks(slc, cfg, logger, masks, excl, n_null, 0, inject_px,
                           tag="faces_" + "_".join(pyramids) + ("_diag" if mode == "diag" else ""),
                           aspect=geo["da"] / geo["dr"])


# =============================================================================
#  TOMOGRAMMA LUNGO UNA LINEA SULLA FACCIA DI UNA PIRAMIDE (come in Sez. 5)
# =============================================================================

def _pyr_px_fn(geo: dict, name: str):
    """Funzione (e, n, z) -> (riga, colonna) per la piramide nota `name`: offset
    in metri Est/Nord dal centro e quota z (layover: z sposta di z/tan(inc))."""
    la, lo, side, H = PYRAMIDS[name]
    srp, u_rg, u_az, inc = geo["srp"], geo["u_rg"], geo["u_az"], geo["inc"]
    lat0, lon0, h0 = _ecf2geo(srp)
    east = np.array([-np.sin(lon0), np.cos(lon0), 0.0])
    north = np.array([-np.sin(lat0) * np.cos(lon0), -np.sin(lat0) * np.sin(lon0),
                      np.cos(lat0)])
    d = _geo2ecf(np.radians(la), np.radians(lo), h0) - srp
    xc, yc = d @ u_rg, d @ u_az

    def px(e_off, n_off, z=0.0):
        v = e_off * east + n_off * north
        x = xc + v @ u_rg - z / np.tan(inc)
        y = yc + v @ u_az
        return (geo["n_az"] // 2 + (y - geo["cy"]) / geo["da"],
                geo["n_rg"] // 2 + (x - geo["cx"]) / geo["dr"])

    return px, side / 2, H


def face_line_points(geo: dict, name: str, line: str, h_frac: float = 0.3):
    """Pixel (righe, colonne) lungo una linea NORD->SUD su una faccia illuminata.

    "W":   linea orizzontale sulla faccia OVEST (illuminata in tutte e due le date), a
           quota h_frac*H, da Nord a Sud.
    "Nup": dal centro del bordo NORD alla cima, risalendo la faccia NORD (illuminata il
           16/09, quasi radente il 04/10): procede da Nord verso Sud salendo."""
    px, half, H = _pyr_px_fn(geo, name)
    if line == "W":
        z = h_frac * H
        e = -half * (1 - h_frac)
        span = 0.92 * half * (1 - h_frac)
        p0, p1 = px(e, +span, z), px(e, -span, z)
    elif line == "Nup":
        p0, p1 = px(0.0, half * 0.97, 0.03 * H), px(0.0, half * 0.12, 0.88 * H)
    else:
        raise ValueError(f"linea sconosciuta: {line}")
    n = int(2 * np.hypot(p1[0] - p0[0], p1[1] - p0[1])) + 2
    t = np.linspace(0, 1, n)
    r = np.round(p0[0] + t * (p1[0] - p0[0])).astype(int)
    c = np.round(p0[1] + t * (p1[1] - p0[1])).astype(int)
    keep = np.r_[True, (np.diff(r) != 0) | (np.diff(c) != 0)]
    return r[keep], c[keep]


def snap_to_ridge(slc: np.ndarray, rows: np.ndarray, cols: np.ndarray, h: int,
                  r_az: int = 12, r_rg: int = 45, smooth: int = 15):
    """Aggancia ogni campione della linea al massimo locale di ampiezza (media 3x3)
    nell'intorno (+-r_az righe, +-r_rg colonne): come nel paper, la linea e'
    tracciata sul tratto visibile dell'SLC. Gli spostamenti sono poi regolarizzati
    con un filtro mediano lungo la linea. Ritorna (righe, colonne) agganciate,
    senza duplicati consecutivi."""
    n_az, n_rg = slc.shape
    a = np.abs(slc).astype(np.float32)
    k = np.zeros_like(a)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            k += np.roll(np.roll(a, dy, axis=0), dx, axis=1)
    k /= 9.0
    dr_ = np.zeros(rows.size, dtype=int)
    dc_ = np.zeros(rows.size, dtype=int)
    for i, (r, c) in enumerate(zip(rows, cols)):
        a0, a1 = max(h, r - r_az), min(n_az - h, r + r_az + 1)
        b0, b1 = max(h, c - r_rg), min(n_rg - h, c + r_rg + 1)
        if a1 <= a0 or b1 <= b0:
            continue
        win = k[a0:a1, b0:b1]
        iy, ix = np.unravel_index(win.argmax(), win.shape)
        dr_[i], dc_[i] = a0 + iy - r, b0 + ix - c
    if smooth > 1 and rows.size >= smooth:
        pad = smooth // 2
        for arr in (dr_, dc_):
            ap = np.pad(arr, pad, mode="edge")
            arr[:] = np.array([np.median(ap[i:i + smooth]) for i in range(arr.size)]).astype(int)
    rs = np.clip(rows + dr_, h, n_az - h)
    cs = np.clip(cols + dc_, h, n_rg - h)
    keep = np.r_[True, (np.diff(rs) != 0) | (np.diff(cs) != 0)]
    return rs[keep], cs[keep]


def _detrend_steps(Y: np.ndarray) -> np.ndarray:
    """Toglie da Y (passi x pixel) la media e il trend lineare lungo i passi."""
    i = np.arange(Y.shape[0], dtype=np.float64) - (Y.shape[0] - 1) / 2.0
    Yd = Y.astype(np.complex128)
    Yd = Yd - Yd.mean(axis=0, keepdims=True)
    c1 = (Yd * i[:, None]).sum(axis=0) / (i ** 2).sum()
    return (Yd - c1[None, :] * i[:, None]).astype(Y.dtype)


def line_tomogram(slc: np.ndarray, cfg: Config, logger: logging.Logger, geo: dict,
                  pyramid: str = "khufu", line: str = "W", h_frac: float = 0.3,
                  tomo_nd: int = 192, z_range: Tuple[float, float] = (-10.0, 150.0),
                  n_ctrl: int = 3, a_pct: float = 75.0, seed: int = 0,
                  tomo_dz: float = 0.4, snap: bool = True,
                  tomo_chain: int = 1) -> dict:
    """Tomogramma |h(z)| normalizzato (jet 0-1) lungo una linea di pixel validi su
    una faccia, accanto a n_ctrl linee di CONTROLLO identiche per forma, traslate
    sul terreno fuori dalle piramidi (con lo stesso aggancio al tratto luminoso).
    Pixel validi: ampiezza >= percentile a_pct dell'immagine e qualita' di
    tracking >= mediana delle linee. Seconda versione con media e trend lineare
    sui passi rimossi (la componente costante finisce a z = 0)."""
    import dataclasses
    cfg2 = dataclasses.replace(cfg, N_D=int(tomo_nd), z_min=float(z_range[0]),
                               z_MAX=float(z_range[1]), dz=float(tomo_dz),
                               chain_len=int(tomo_chain))
    n_az, n_rg = slc.shape
    h = cfg.patch // 2
    logger.info(f"TOMOGRAMMA: {pyramid} linea {line} (h={h_frac}), N_D={tomo_nd}, "
                f"z in [{z_range[0]:.0f}, {z_range[1]:.0f}] m, aggancio={snap}")

    g0, g1 = face_line_points(geo, pyramid, line, h_frac)           # geometrica
    ok = (g0 >= h) & (g0 <= n_az - h) & (g1 >= h) & (g1 <= n_rg - h)
    g0, g1 = g0[ok], g1[ok]
    if g0.size < 20:
        raise RuntimeError("la linea cade fuori scena (controlla --center_x/--center_y)")
    r0, c0 = snap_to_ridge(slc, g0, g1, h) if snap else (g0, g1)

    shape = slc.shape
    excl = np.zeros(shape, dtype=bool)
    for nm in PYRAMIDS:
        excl |= polygon_mask(shape, pyramid_faces_px(geo, nm)["base"])
    excl = _dilate(excl, 25)
    rng = np.random.default_rng(seed)
    ctrl = []
    for _ in range(40000):
        if len(ctrl) >= n_ctrl:
            break
        dr_, dc_ = int(rng.integers(-n_az, n_az)), int(rng.integers(-n_rg, n_rg))
        rr, cc = g0 + dr_, g1 + dc_
        if rr.min() < h + 14 or rr.max() > n_az - h - 14 or cc.min() < h + 50 or cc.max() > n_rg - h - 50:
            continue
        if excl[rr, cc].any():
            continue
        if any(abs(dr_ - a) + abs(dc_ - b) < 60 for a, b in [x[2] for x in ctrl]):
            continue
        rs, cs = snap_to_ridge(slc, rr, cc, h) if snap else (rr, cc)
        if excl[rs, cs].any():
            continue
        ctrl.append((rs, cs, (dr_, dc_)))
    logger.info(f"  linea principale: {r0.size} campioni; controlli: {len(ctrl)}")

    lines = [(r0, c0)] + [(a, b) for a, b, _ in ctrl]
    sizes = [len(a) for a, _ in lines]
    rows = np.concatenate([a for a, _ in lines]); cols = np.concatenate([b for _, b in lines])

    spectrum = dft2(slc, logger)
    patches, _, _ = generate_doppler_subapartures(spectrum, cfg2, logger, rows, cols)
    del spectrum
    Y, Q = estimate_micromotion(patches, cfg2, logger, return_quality=True)
    del patches
    amp = np.abs(slc)[rows, cols]
    a_thr = float(np.percentile(np.abs(slc), a_pct))
    q_thr = float(np.median(Q))
    valid = (amp >= a_thr) & (Q >= q_thr)
    A_pinv = _tomo_pinv(Y.shape[0], cfg2)
    zc = cfg2.z

    dK = (4 * np.pi * (cfg2.A_orbit / (tomo_nd - 1))) / (
        cfg2.lambda_sound * cfg2.r_ref * np.sin(np.deg2rad(45.0)))
    z_amb = 2 * np.pi / dK / (2 * np.pi if cfg2.steer_2pi else 1.0)
    z_res = cfg2.delta_T
    logger.info(f"  risoluzione nominale lambda*R/2A = {z_res:.2f} m ; periodo di "
                f"ambiguita' in z = {z_amb:.0f} m ; soglie: ampiezza >= {a_thr:.3g}, "
                f"qualita' >= {q_thr:.2f}")

    secs, off = [], 0
    for k, n_ in enumerate(sizes):
        sel = valid[off:off + n_]
        Yk = Y[:, off:off + n_][:, sel]
        if sel.sum() < 3:
            Yk = np.zeros((Y.shape[0], 1), dtype=Y.dtype)
            logger.warning(f"  linea {k}: meno di 3 pixel validi, sezione vuota")
        sct = dict(Y=Yk, n=int(sel.sum()), rows=rows[off:off + n_][sel],
                   cols=cols[off:off + n_][sel])
        for tag, Yv in (("", Yk), ("d", _detrend_steps(Yk))):
            Hv = np.abs(A_pinv @ Yv.astype(np.complex64))
            prof = Hv.mean(axis=1)
            sct["H" + tag] = Hv
            sct["Hn" + tag] = Hv / (Hv.max() + 1e-30)
            sct["prof" + tag] = prof
            sct["score" + tag] = float(prof.max() / (np.median(prof) + 1e-30))
            sct["zpk" + tag] = float(zc[int(np.argmax(prof))])
        secs.append(sct)
        off += n_
    for k, sct in enumerate(secs):
        nm = "PRINCIPALE" if k == 0 else f"controllo {k}"
        logger.info(f"  {nm:12s}: {sct['n']:4d} validi su {sizes[k]:3d} | picco/mediana "
                    f"grezzo={sct['score']:.1f} (z={sct['zpk']:.0f}) | senza media/trend="
                    f"{sct['scored']:.1f} (z={sct['zpkd']:.0f})")

    # --- figura ---
    ncol = len(secs)
    fig = plt.figure(figsize=(5.0 * ncol, 17), dpi=min(cfg.dpi, 100))
    gs = fig.add_gridspec(5, ncol, height_ratios=[1.1, 1.3, 1.3, 0.9, 0.8])
    a_img = np.abs(slc)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.imshow(20 * np.log10(a_img / np.percentile(a_img, 99.9) + 1e-6), cmap="gray",
               vmin=-45, vmax=0, aspect=geo["da"] / geo["dr"], extent=(0, n_rg, n_az, 0))
    for k, sct in enumerate(secs):
        ax0.plot(sct["cols"], sct["rows"], ".", ms=3, color="yellow" if k == 0 else "cyan")
    pf = pyramid_faces_px(geo, pyramid)
    cs_ = ["SW", "SE", "NE", "NW", "SW"]
    ax0.plot([pf["corners"][c][1] for c in cs_], [pf["corners"][c][0] for c in cs_], "r-", lw=0.8)
    ax0.set_title(f"{pyramid}: linea {line} Nord->Sud (giallo, pixel validi) e controlli (ciano); "
                  f"base prevista in rosso")
    for k, sct in enumerate(secs):
        for row_, tag, ttl in ((1, "", "|h| normalizzato (come nel paper)"),
                               (2, "d", "|h| senza media e trend sui passi")):
            ax = fig.add_subplot(gs[row_, k])
            imx = ax.imshow(sct["Hn" + tag], cmap="jet", aspect="auto", origin="lower",
                            vmin=0, vmax=1, extent=(0, max(1, sct["n"]), zc[0], zc[-1]),
                            interpolation="nearest")
            ax.set_title(("PRINCIPALE" if k == 0 else f"controllo {k}") +
                         f" - {ttl} ({sct['n']} px)", fontsize=8)
            ax.set_xlabel("pixel lungo la linea"); ax.set_ylabel("z [m]")
            if k == ncol - 1:
                fig.colorbar(imx, ax=ax, fraction=0.05)
        axy = fig.add_subplot(gs[3, k])
        axy.imshow(np.abs(sct["Y"]), cmap="viridis", aspect="auto", origin="lower",
                   extent=(0, max(1, sct["n"]), 0, tomo_nd))
        axy.set_title("|Y|: spostamento (px) per passo k", fontsize=9)
        axy.set_xlabel("pixel lungo la linea"); axy.set_ylabel("passo k")
    half = max(1, ncol // 2)
    for idx, (tag, ttl) in enumerate((("", "profilo medio |h| / mediana (grezzo)"),
                                      ("d", "senza media e trend"))):
        axp = fig.add_subplot(gs[4, :half] if idx == 0 else gs[4, half:])
        for k, sct in enumerate(secs):
            axp.plot(zc, sct["prof" + tag] / (np.median(sct["prof" + tag]) + 1e-30),
                     lw=2 if k == 0 else 1, color="red" if k == 0 else None,
                     label="principale" if k == 0 else f"controllo {k}")
        axp.set_xlabel("z [m]"); axp.set_title(f"{ttl} (nominale {z_res:.1f} m, ambiguita' {z_amb:.0f} m)", fontsize=9)
        axp.legend(fontsize=7); axp.grid(alpha=0.3)
    plt.tight_layout()
    out_png = os.path.join(cfg.out_dir, f"tomo_{pyramid}_{line}.png")
    plt.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    np.savez(os.path.join(cfg.out_dir, f"tomo_{pyramid}_{line}.npz"), z=zc,
             **{f"H{k}": s_["H"] for k, s_ in enumerate(secs)},
             **{f"Hd{k}": s_["Hd"] for k, s_ in enumerate(secs)},
             **{f"Y{k}": s_["Y"] for k, s_ in enumerate(secs)})
    logger.info(f"  → {out_png}")
    return {"sections": secs, "z": zc, "z_amb": z_amb, "z_res": z_res}


# =============================================================================
#  VOLUME TOMOGRAFICO DI UNA PIRAMIDE (tomografia del paper su tutti i pixel)
# =============================================================================

def _zbandpass(V: np.ndarray, dz: float, p_lo: float, p_hi: float,
               att_db: float = 3.1, order: int = 4) -> np.ndarray:
    """Passa-banda lungo z (asse 1) di profili reali V (pixel x z), fase zero.
    Banda in periodi verticali [p_lo, p_hi] m; risposta tipo Butterworth di ordine
    `order` con attenuazione att_db ai due bordi. Il profilo e' specchiato (niente
    salto ai bordi) e la componente costante e' tolta."""
    F = V.shape[1]
    X = np.concatenate([V, V[:, ::-1]], axis=1)
    f = np.fft.rfftfreq(2 * F, dz)
    e2 = 10 ** (att_db / 10) - 1                   # |H|^2 = 1/(1+e2) ai bordi
    with np.errstate(divide="ignore"):
        Hf = 1 / np.sqrt((1 + e2 * (f * p_lo) ** (2 * order))
                         * (1 + e2 * (1 / (f * p_hi)) ** (2 * order)))
    Hf[0] = 0.0
    return np.fft.irfft(np.fft.rfft(X, axis=1) * Hf, 2 * F, axis=1)[:, :F].astype(np.float32)


def tomo_volume(slc: np.ndarray, geo: dict, cfg: Config, logger: logging.Logger,
                name: str = "khufu", band: Tuple[float, float] = (3.1, 30.0),
                block: int = 40000) -> dict:
    """Tomografia del paper (sub-aperture a catena -> tracking -> h = A+ Y) su TUTTI
    i pixel dell'impronta della piramide `name` (base + facce con layover), a blocchi
    di `block` pixel (il picco di memoria resta ~1 GB). Per ogni pixel: |h(z)| sulla
    griglia cfg.z e il suo passa-banda in z (periodi `band` m, -3.1 dB ai bordi):
    evidenzia le variazioni di riflettivita' lungo z (cambi di densita')."""
    n_az, n_rg = slc.shape
    h = cfg.patch // 2
    pf = pyramid_faces_px(geo, name)
    mask = np.zeros(slc.shape, dtype=bool)
    for k in ("base", "N", "E", "S", "W"):
        mask |= polygon_mask(slc.shape, pf[k])
    mask[:h], mask[n_az - h:], mask[:, :h], mask[:, n_rg - h:] = False, False, False, False
    rows, cols = np.nonzero(mask)
    npx, F = rows.size, len(cfg.z)
    dK = (4 * np.pi * (cfg.A_orbit / (cfg.N_D - 1))) / (
        cfg.lambda_sound * cfg.r_ref * np.sin(np.deg2rad(45.0)))
    z_amb = 2 * np.pi / dK / (2 * np.pi if cfg.steer_2pi else 1.0)
    logger.info(f"VOLUME {name}: {npx} pixel, N_D={cfg.N_D}, z in [{cfg.z[0]:.0f}, "
                f"{cfg.z[-1]:.0f}] m passo {cfg.dz:g} m ({F} campioni), ambiguita' "
                f"{z_amb:.0f} m; passa-banda {band[0]:g}-{band[1]:g} m (-3.1 dB), "
                f"blocchi di {block} pixel")
    if z_amb < cfg.z[-1] - cfg.z[0]:
        logger.warning("  l'intervallo in z supera il periodo di ambiguita': aumenta N_D")

    spectrum = dft2(slc, logger)
    Hraw = np.zeros((npx, F), dtype=np.float16)
    Hbp = np.zeros((npx, F), dtype=np.float16)
    Q = np.zeros(npx, dtype=np.float32)
    A_pinv = None
    t0 = time.time()
    for b0 in range(0, npx, block):
        b1 = min(npx, b0 + block)
        shifts, _, _ = generate_doppler_subapartures(spectrum, cfg, logger,
                                                     rows[b0:b1], cols[b0:b1])
        Y, Q[b0:b1] = estimate_micromotion(shifts, cfg, logger, return_quality=True)
        del shifts
        if A_pinv is None:
            A_pinv = _tomo_pinv(Y.shape[0], cfg)
        Hb = np.abs(A_pinv @ Y)                        # (F, pixel) float32
        Hraw[b0:b1] = Hb.T
        Hbp[b0:b1] = _zbandpass(Hb.T, cfg.dz, band[0], band[1])
        logger.info(f"  blocco {b1}/{npx} pixel ({time.time() - t0:.0f} s)")
    del spectrum
    return {"name": name, "rows": rows, "cols": cols, "z": cfg.z, "H": Hraw,
            "Hbp": Hbp, "Q": Q, "z_amb": z_amb, "band": band, "faces": pf}


def plot_tomo_volume(vol: dict, slc: np.ndarray, cfg: Config, out_path: str) -> None:
    """Fette orizzontali del volume passa-banda a 6 profondita' e due sezioni
    verticali per il centro dell'impronta (lungo il range e lungo l'azimuth), con l'SLC."""
    rows, cols, z, V = vol["rows"], vol["cols"], vol["z"], vol["Hbp"]
    r0, r1, c0, c1 = rows.min(), rows.max() + 1, cols.min(), cols.max() + 1
    lim = float(np.percentile(np.abs(V.astype(np.float32)), 99.5)) + 1e-30
    ar, ac = int(np.median(rows)), int(np.median(cols))   # centro dell'impronta
    fig, axs = plt.subplots(3, 3, figsize=(18, 17))
    a = np.abs(slc[r0:r1, c0:c1])
    axs[0, 0].imshow(20 * np.log10(a / np.percentile(a, 99.9) + 1e-6), cmap="gray",
                     vmin=-40, vmax=0, extent=(c0, c1, r1, r0), aspect="auto")
    axs[0, 0].axhline(ar, color="y", lw=0.8); axs[0, 0].axvline(ac, color="c", lw=0.8)
    axs[0, 0].set_title(f"SLC {vol['name']} (dB); sezioni: giallo = range, ciano = azimuth")
    axs[0, 0].set_xlabel("range (pixel)"); axs[0, 0].set_ylabel("azimuth (pixel)")
    for ax, zk in zip(axs.ravel()[1:7], np.linspace(z[-1] - 10, z[0] + 10, 6)):
        k = int(np.argmin(np.abs(z - zk)))
        img = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        img[rows - r0, cols - c0] = V[:, k]
        im = ax.imshow(img, cmap="seismic", vmin=-lim, vmax=lim,
                       extent=(c0, c1, r1, r0), aspect="auto", interpolation="nearest")
        ax.set_title(f"passa-banda, z = {z[k]:.1f} m")
        ax.set_xlabel("range (pixel)"); ax.set_ylabel("azimuth (pixel)")
    for ax, sel, x, lab in ((axs[2, 1], rows == ar, cols, "range (pixel)"),
                            (axs[2, 2], cols == ac, rows, "azimuth (pixel)")):
        o = np.argsort(x[sel])
        ax.imshow(V[sel][o].astype(np.float32).T, cmap="seismic", vmin=-lim, vmax=lim,
                  origin="lower", aspect="auto", interpolation="nearest",
                  extent=(x[sel].min(), x[sel].max() + 1, z[0], z[-1]))
        ax.set_title(f"sezione verticale per il centro lungo {lab.split()[0]}")
        ax.set_xlabel(lab); ax.set_ylabel("z [m]")
    fig.colorbar(im, ax=axs[:2, :], fraction=0.02, label="|h| passa-banda (unita' arbitrarie)")
    fig.suptitle(f"Volume tomografico {vol['name']}: {rows.size} pixel, passo z {cfg.dz:g} m, "
                 f"ambiguita' {vol['z_amb']:.0f} m, passa-banda {vol['band'][0]:g}-"
                 f"{vol['band'][1]:g} m (-3.1 dB)", fontsize=14)
    plt.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
#  VOXEL 3D DI UNA PIRAMIDE (onde della massa 1, passa-basso lungo z)
# =============================================================================

def plot_orientation(slc: np.ndarray, geo: dict, name: str, out_path: str,
                     map_path: Optional[str] = None) -> None:
    """Controllo N/S/E/O: SLC con la base prevista, gli spigoli (NO, NE, SE, SO),
    l'apice (con layover) e le direzioni Nord ed Est nell'immagine. map_path: json con
    la trasformazione (x = colonna, y = riga) -> pixel del png (x_png = ax x + bx, ecc.)."""
    px, half, H = _pyr_px_fn(geo, name)
    a = np.abs(slc)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.imshow(20 * np.log10(a / np.percentile(a, 99.9) + 1e-6), cmap="gray", vmin=-40, vmax=0,
              aspect=geo["da"] / geo["dr"])
    cor = {"SO": (-half, -half), "SE": (half, -half), "NE": (half, half), "NO": (-half, half)}
    pts = [px(e, n) for e, n in cor.values()]
    ax.plot([p[1] for p in pts + pts[:1]], [p[0] for p in pts + pts[:1]], "r-", lw=1.2)
    for lab, p in zip(cor, pts):
        ax.annotate(lab, (p[1], p[0]), color="yellow", fontsize=12, weight="bold")
    ap = px(0.0, 0.0, H)
    ax.plot(ap[1], ap[0], "c^", ms=10)
    ax.annotate("apice (layover)", (ap[1], ap[0]), color="cyan", fontsize=10)
    c0 = px(0.0, 0.0)
    for lab, (e, n), col in (("Nord 50 m", (0.0, 50.0), "lime"), ("Est 50 m", (50.0, 0.0), "orange")):
        p1 = px(e, n)
        ax.annotate("", (p1[1], p1[0]), (c0[1], c0[0]),
                    arrowprops=dict(arrowstyle="->", color=col, lw=2))
        ax.annotate(lab, (p1[1], p1[0]), color=col, fontsize=11, weight="bold")
    ax.set_title(f"{name}: orientamento (base prevista in rosso, frecce dal centro della base)")
    ax.set_xlabel("range (pixel, verso destra = lontano dal sensore)")
    ax.set_ylabel("azimuth (pixel)")
    if map_path:
        # trasformazione assi (colonna, riga dell'SLC) -> pixel del png salvato con
        # bbox_inches="tight" (riquadro stretto + pad_inches), per la pagina dell'ecografia
        fig.set_dpi(110)                                 # stessa risoluzione del salvataggio
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer()).padded(plt.rcParams["savefig.pad_inches"])
        (x0, y0), (x1, y1) = ax.transData.transform([(0.0, 0.0), (1.0, 1.0)]) / fig.dpi
        with open(map_path, "w") as f:
            json.dump({"ax": (x1 - x0) * 110, "bx": (x0 - bb.x0) * 110,
                       "ay": -(y1 - y0) * 110, "by": (bb.y1 - y0) * 110,
                       "n_az": int(slc.shape[0]), "n_rg": int(slc.shape[1])}, f)
    plt.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)


def known_occupancy(known: list, e: np.ndarray, z_b: float, cs: float):
    """Griglia delle parti note sul reticolo dei voxel: colonne (indice c = i_e*n_col + i_n)
    dentro il riquadro delle strutture e celle m (centro z = z_b + (m + 0.5) cs).
    Una cella e' occupata se il centro sta nell'involucro convesso dei vertici del wireframe
    di una struttura, allargato della sua incertezza. Ritorna (colonne, m0, occ[col, m])."""
    from scipy.spatial import ConvexHull
    hulls, P = [], []
    for s in known:
        p = np.asarray(s["seg"], dtype=float).reshape(-1, 3)
        hulls.append((ConvexHull(p, qhull_options="QJ").equations, float(s["incertezza"])))
        P.append(p)
    P = np.concatenate(P)
    r = max(i for _, i in hulls) + cs
    lo, hi = P.min(axis=0) - r, P.max(axis=0) + r
    ie = np.nonzero((e >= lo[0]) & (e <= hi[0]))[0]
    inn = np.nonzero((e >= lo[1]) & (e <= hi[1]))[0]
    m0 = int(np.floor((lo[2] - z_b) / cs))
    m1 = int(np.ceil((hi[2] - z_b) / cs)) + 1
    EI, NI, M = np.meshgrid(ie, inn, np.arange(m0, m1), indexing="ij")
    X = np.column_stack([e[EI.ravel()], e[NI.ravel()], z_b + (M.ravel() + 0.5) * cs])
    occ = np.zeros(X.shape[0], dtype=bool)
    for eq, inc in hulls:
        occ |= (X @ eq[:, :3].T + eq[:, 3]).max(axis=1) <= inc
    cols = (EI[:, :, 0] * e.size + NI[:, :, 0]).ravel()
    return cols, m0, occ.reshape(cols.size, m1 - m0)


def known_match_scan(vox: np.ndarray, nj: np.ndarray, n_col: int, cols: np.ndarray, m0: int,
                     occ: np.ndarray, s_max: int, k_max: int = 30, n_null: int = 500,
                     min_shift: float = 30.0, cs: float = 1.0, seed: int = 0) -> dict:
    """Corrispondenza tra i cubi scelti (i K con valore piu' basso per colonna, K = 1..k_max,
    come nella pagina con il cursore del valore al massimo) e le parti note, per ogni
    spostamento verticale s = -s_max..s_max cubi del campo dei voxel.
    hits(K, s) = cubi scelti dentro le parti note; atteso(K, s) = stessi cubi per colonna
    messi a caso nella colonna; z = (hits - atteso) / sqrt(atteso).
    Significativita' corretta per tutta la ricerca (K x s): lo stesso massimo di z e'
    calcolato su n_null controlli in cui lo schema dei cubi scelti (in profondita' sotto
    la superficie) e' traslato in orizzontale di almeno min_shift m sulla griglia;
    p = (1 + controlli con max >= osservato) / (1 + n_null)."""
    off = np.concatenate([[0], np.cumsum(nj)])
    col = np.repeat(np.arange(nj.size), nj)
    j = np.arange(off[-1]) - off[col]
    order = np.lexsort((j, vox, col))                   # come la pagina: valore, poi dal fondo
    rank = np.empty(off[-1], dtype=np.int32)
    rank[order] = np.arange(off[-1]) - off[col[order]]
    keep = rank < k_max
    d_all = (nj[col] - 1 - j)[keep]                      # schema: profondita' sotto la superficie
    c_all, r_all = col[keep], rank[keep]
    o2 = np.concatenate([[0], np.cumsum(np.bincount(c_all, minlength=nj.size))])
    srt = np.argsort(c_all, kind="stable")
    d_all, r_all = d_all[srt], r_all[srt]
    s = np.arange(-s_max, s_max + 1)
    nm = occ.shape[1]
    reg_nj = nj[cols].astype(float)
    P = np.concatenate([np.zeros((cols.size, 1)), np.cumsum(occ, axis=1)], axis=1)

    def cnt_s():                                         # celle occupate nella colonna, per s
        a = np.clip(s[None, :] - m0, 0, nm)
        b = np.clip(nj[cols][:, None] + s[None, :] - m0, 0, nm)
        return np.take_along_axis(P, b, 1) - np.take_along_axis(P, a, 1)
    CNT = cnt_s()                                        # (colonne, s)

    def stat(donor: np.ndarray):
        """donor[i] = colonna da cui prendere lo schema per la colonna cols[i]."""
        n_d = o2[donor + 1] - o2[donor]                  # cubi scelti nella colonna donatrice
        ci = np.repeat(np.arange(donor.size), n_d)
        src = np.repeat(o2[donor] - np.concatenate([[0], np.cumsum(n_d)[:-1]]), n_d) \
            + np.arange(n_d.sum())
        dd, rr = d_all[src], r_all[src]
        jj = nj[cols][ci] - 1 - dd
        v = jj >= 0
        ci, jj, rr = ci[v], jj[v], rr[v]
        m = jj[:, None] + s[None, :] - m0
        ok = (m >= 0) & (m < nm)
        hv = np.where(ok, occ[ci[:, None], np.clip(m, 0, nm - 1)], False).astype(float)
        H = np.zeros((k_max, s.size)); np.add.at(H, rr, hv)
        H = np.cumsum(H, axis=0)                          # hits per K = 1..k_max
        kc = np.zeros((cols.size, k_max)); np.add.at(kc, (ci, rr), 1.0)
        kc = np.cumsum(kc, axis=1)                        # cubi scelti per colonna, per K
        E = (kc / np.maximum(reg_nj, 1)[:, None]).T @ CNT
        Z = np.where(E >= 5, (H - E) / np.sqrt(np.maximum(E, 1e-9)), -np.inf)
        return Z, H, E

    Z, H, E = stat(cols.copy())
    kb, sb = np.unravel_index(np.argmax(Z), Z.shape)
    rng = np.random.default_rng(seed)
    ie, inn = cols // n_col, cols % n_col
    zmax = []
    while len(zmax) < n_null:
        de, dn = rng.integers(-n_col, n_col, 2)
        if np.hypot(de, dn) * cs < min_shift:
            continue
        dc = ((ie + de) % n_col) * n_col + (inn + dn) % n_col
        zmax.append(float(np.max(stat(dc)[0])))
    zmax = np.array(zmax)
    p = (1 + int((zmax >= Z[kb, sb]).sum())) / (1 + n_null)
    return {"Z": Z, "H": H, "E": E, "s": s, "K": int(kb) + 1, "shift": int(s[sb]),
            "z": float(Z[kb, sb]), "hits": float(H[kb, sb]), "exp": float(E[kb, sb]),
            "p": p, "zmax_null": zmax}


def plot_match_scan(res: dict, name: str, cs: float, out_path: str) -> None:
    """Mappa z(K, spostamento) e distribuzione del massimo nei controlli."""
    fig, ax = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={"width_ratios": [2.2, 1]})
    Z = np.where(np.isfinite(res["Z"]), res["Z"], np.nan)
    lim = np.nanmax(np.abs(Z))
    im = ax[0].imshow(Z, aspect="auto", origin="lower", cmap="seismic", vmin=-lim, vmax=lim,
                      extent=(res["s"][0] * cs, res["s"][-1] * cs, 0.5, Z.shape[0] + 0.5))
    ax[0].plot(res["shift"] * cs, res["K"], "k*", ms=14)
    ax[0].set_xlabel("spostamento verticale dei voxel [m] (+ = in alto)")
    ax[0].set_ylabel("massimo per colonna K")
    ax[0].set_title(f"{name}: z = (dentro - atteso) / sqrt(atteso); stella = massimo")
    fig.colorbar(im, ax=ax[0])
    ax[1].hist(res["zmax_null"], bins=40, color="gray")
    ax[1].axvline(res["z"], color="r", lw=2)
    ax[1].set_title(f"massimo nei {res['zmax_null'].size} controlli; osservato z = {res['z']:.2f}, "
                    f"p = {res['p']:.3f}", fontsize=10)
    ax[1].set_xlabel("max z su tutta la ricerca")
    plt.tight_layout()
    plt.savefig(out_path, dpi=110)
    plt.close(fig)


def voxel_pyramid(slc: np.ndarray, geo: dict, cfg: Config, logger: logging.Logger,
                  name: str, n_col: int, out_dir: str) -> dict:
    """Dati per la pagina 3D dei voxel della piramide `name`.

    Griglia n_col x n_col di colonne (Est, Nord) sulla base (n_col <= 0: massima
    risoluzione, una colonna per cella azimuth dell'SLC). La colonna prende i valori del
    pixel dove il radar vede la superficie della piramide in quella colonna,
    px(e, n, cap) con cap = H * (1 - max(|e|, |n|) / (side/2)) (layover: spostato verso
    il sensore di cap / tan(inc)). L'onda del pixel scende in verticale dalla superficie
    (t = 0) fino a z_bottom: i tiri sono distribuiti sull'altezza H e un ciclo della molla
    copre H. Voxel = cubi di lato = passo della griglia, tra z_bottom e la superficie;
    per ogni cubo si esporta
    l'ampiezza filtrata sqrt(E-O^2 + N-S^2) al suo centro (uint8 compresso a radice,
    riferimento = 99.5 percentile) e, per distanza dalla superficie (in cubi), il
    guadagno della normalizzazione per profondita'. Nel browser il cursore mostra i cubi
    con ampiezza SOTTO il valore scelto."""
    px, half, H = _pyr_px_fn(geo, name)
    if n_col <= 0:
        n_col = int(round(2 * half / geo["da"])) + 1
    cs = 2 * half / (n_col - 1)                         # lato del cubo = passo della griglia
    e = np.linspace(-half, half, n_col)
    EE, NN = np.meshgrid(e, e, indexing="ij")
    cap = (H * (1 - np.maximum(np.abs(EE), np.abs(NN)) / half)).ravel()   # superficie
    rc = np.array([px(a, b, c) for a, b, c in zip(EE.ravel(), NN.ravel(), cap)])
    rows, cols = np.round(rc[:, 0]).astype(int), np.round(rc[:, 1]).astype(int)
    ok = (rows >= 0) & (rows < slc.shape[0]) & (cols >= 0) & (cols < slc.shape[1])
    if ok.sum() < 0.5 * ok.size:
        raise RuntimeError(f"{name}: meno di meta' della base cade nella scena")
    logger.info(f"Voxel {name}: base {2 * half:.0f} m, altezza {H:g} m, griglia "
                f"{n_col}x{n_col} (passo {cs:.2f} m), {int(ok.sum())}/{ok.size} colonne "
                f"dentro la scena (pixel della superficie, con layover)")
    os.makedirs(out_dir, exist_ok=True)
    # cubi j = 0.. dal fondo: centro zc = z_bottom + (j + 0.5) cs, sotto la superficie
    nj = np.where(ok, np.floor((cap - cfg.z_bottom) / cs), 0).astype(int)
    j = np.concatenate([np.arange(k) for k in nj[ok]])
    depth = np.repeat(cap[ok], nj[ok]) - (cfg.z_bottom + (j + 0.5) * cs)   # sotto la superficie
    n_d = cfg.N_D if cfg.N_D > 0 else int(round(H / cs)) + 1   # 0: un tiro per cubo (passo = cs)
    cfg2 = dataclasses.replace(cfg, pyr_height=H, N_D=n_d)
    res = vector_chain_waveforms(dft2(slc, logger), cfg2, logger, rows[ok], cols[ok], out_dir,
                                 samples=(nj[ok], depth))
    plot_vector_waveforms(res, cfg2, os.path.join(out_dir, "vector_waveforms.png"))
    A = res["amp_s"]                                    # ampiezza filtrata al centro di ogni cubo
    ref = float(np.percentile(A, 99.5)) + 1e-30         # riferimento del cursore
    # normalizzazione per profondita' (interruttore nella pagina): A / rms dei cubi alla stessa
    # distanza dalla superficie (i = nj - 1 - j), riportata al suo 99.5 percentile
    i_top = np.repeat(nj[ok] - 1, nj[ok]) - j
    rms = np.sqrt(np.bincount(i_top, A.astype(np.float64) ** 2) / np.bincount(i_top)) + 1e-30
    gain = ref / (rms * np.percentile(A / rms[i_top], 99.5))
    vox = np.round(255 * np.sqrt(np.clip(A / ref / 1.5, 0, 1))).astype(np.uint8)   # A/ref = 1.5 (q/255)^2

    def b64(a, dt):
        return base64.b64encode(np.ascontiguousarray(a.astype(dt)).tobytes()).decode("ascii")

    plot_orientation(slc, geo, name, os.path.join(out_dir, "orientamento.png"))
    ks = load_known_structures(logger).get(name, [])
    known = [{"n": i + 1, "name": s["nome"], "cat": s["categoria"], "color": s["colore"],
              "fonte": s["fonte"], "inc": s["incertezza"], "schem": s["schematica"],
              "nota": s["nota"], "c": s["c"], "seg": s["seg"]}
             for i, s in enumerate(ks)]
    s_max = int(round(100.0 / cs))                      # cursore verticale: +-100 m a passi di un cubo
    match = None
    if ks:
        occ_cols, m0, occ = known_occupancy(ks, e, cfg.z_bottom, cs)
        sc = known_match_scan(vox, nj, n_col, occ_cols, m0, occ, s_max, cs=cs)
        plot_match_scan(sc, name, cs, os.path.join(out_dir, "corrispondenza.png"))
        logger.info(f"  corrispondenza con le parti note ({len(ks)} strutture, {int(occ.sum())} "
                    f"celle): massimo z = {sc['z']:.2f} con K = {sc['K']} e spostamento "
                    f"{sc['shift'] * cs:+.1f} m ({sc['hits']:.0f} cubi dentro, attesi "
                    f"{sc['exp']:.1f}); p corretto per la ricerca = {sc['p']:.3f} "
                    f"({sc['zmax_null'].size} controlli)")
        match = {"cols_b64": b64(occ_cols, "<u4"), "m0": m0, "nm": int(occ.shape[1]),
                 "occ_b64": b64(occ, "u1"), "K": sc["K"], "shift": sc["shift"],
                 "z": round(sc["z"], 2), "p": round(sc["p"], 4), "hits": sc["hits"],
                 "exp": round(sc["exp"], 1), "n_null": int(sc["zmax_null"].size)}
    # nel file solo i k_file cubi con valore piu' basso di ogni colonna (a parita', i primi
    # dal fondo): gli altri non servono al cursore "massimo per colonna" (1..k_file)
    k_file = 30
    col = np.repeat(np.arange(nj.size), nj)
    order = np.lexsort((j, A, col))                     # per colonna, A crescente, poi j
    rank = np.empty(A.size, dtype=np.int64)
    rank[order] = np.arange(A.size) - np.concatenate([[0], np.cumsum(nj)])[col[order]]
    keep = rank < k_file
    logger.info(f"  nel file {int(keep.sum())} cubi di {A.size} (max {k_file} per colonna)")
    return {"name": name, "label": name.capitalize(), "side": round(2 * half, 1), "H": H,
            "zmin": cfg.z_bottom, "n_col": n_col, "dstep": cs,
            "n": res["D"].shape[0], "dz": res["dz"], "half": cfg.wave_half,
            "lpf": cfg.wave_lpf, "ref": ref, "gain": [round(float(g), 4) for g in gain],
            "nj_b64": b64(nj, "<u2"), "nk_b64": b64(np.minimum(nj, k_file), "<u2"),
            "jj_b64": b64(j[keep], "<u2"), "vox_b64": b64(vox[keep], "u1"), "known": known,
            "s_max": s_max, "match": match}


def write_voxel_html(pyramids: list, title: str, out_path: str,
                     template: str = "voxel_template.html", data_ref: Optional[str] = None) -> None:
    """Pagina HTML 3D (Plotly + Tailwind) dal modello `template` (accanto al programma).
    Con `data_ref` la pagina non contiene i dati ma il nome del file .json da cui leggerli."""
    tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), template)
    with open(tpl, encoding="utf-8") as f:
        html = f.read()
    data = json.dumps(data_ref if data_ref else {"title": title, "pyramids": pyramids},
                      separators=(",", ":"))
    try:                      # plotly.js incorporato (pagina autonoma, anche offline)
        import plotly
        with open(os.path.join(os.path.dirname(plotly.__file__), "package_data",
                               "plotly.min.js"), encoding="utf-8") as f:
            lib = "<script>" + f.read() + "</script>"
    except Exception:
        lib = ('<script src="https://cdn.plot.ly/plotly-3.0.1.min.js" '
               'charset="utf-8"></script>')
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(html.replace("__DATA__", data).replace("__PLOTLY__", lib))


BP19_BIG = ("vox_idx_b64", "vox_val_b64", "vox_dep_b64", "png_b64")


def write_bp19_pages(pyramids: list, title: str, out_dir: str, max_mb: float = 20.0) -> list:
    """Pagine --bp19 senza dati (3D e spaccato) + dati in .json di al massimo max_mb MB:
    pyramids_bp19_dati.json (indice: titolo, piramidi senza i campi grandi, elenco dei file)
    e pyramids_bp19_dati_N.json con i campi grandi a pezzi {p, f, k, s}. Ogni .json e'
    racchiuso in bp19_dati(...): la pagina lo carica come script anche da file://."""
    lim = int(max_mb * 1e6)
    pieces = [(i, f, k, p[f][o:o + lim]) for i, p in enumerate(pyramids) for f in BP19_BIG
              for k, o in enumerate(range(0, max(1, len(p[f])), lim))]
    files, cur, size = [], [], 0
    for pc in pieces:                                    # pezzi raggruppati fino a lim
        if cur and size + len(pc[3]) > lim:
            files.append(cur); cur, size = [], 0
        cur.append(pc); size += len(pc[3])
    if cur:
        files.append(cur)
    names = [f"pyramids_bp19_dati_{n + 1}.json" for n in range(len(files))]
    out = []
    def dump(name, obj):
        out.append(os.path.join(out_dir, name))
        with open(out[-1], "w", encoding="utf-8", newline="\n") as f:
            f.write("bp19_dati(" + json.dumps(obj, separators=(",", ":"), ensure_ascii=False) + ");\n")

    for name, grp in zip(names, files):
        dump(name, {"parts": [{"p": i, "f": fl, "k": k, "s": s} for i, fl, k, s in grp]})
    dump("pyramids_bp19_dati.json", {"title": title, "files": names, "pyramids":
         [{k: v for k, v in p.items() if k not in BP19_BIG} for p in pyramids]})
    for tpl, name in (("bp19_template.html", "pyramids_bp19.html"),
                      ("bp19_spaccato_template.html", "pyramids_bp19_spaccato.html")):
        out.append(os.path.join(out_dir, name))
        write_voxel_html(pyramids, title, out[-1], template=tpl, data_ref="pyramids_bp19_dati.json")
    return out


# =============================================================================
#  BIONDI PROTOCOL v1.9 (--bp19): banco R/O, registrazione, fit a ellisse (GPU)
# =============================================================================
# Porting del pacchetto "Biondi Protocol v1.9.0" (github.com/BiondiProtocol/
# Replication-and-Verification-Biondi-Protocol, 10/10/2026): stesse formule di
# primitives.py, measurement.py, focus.py e statistics.py; registrazione e fit in torch.
# Adattamento all'SLC da backprojection (piano al suolo, sovracampionato in azimuth):
# - il bin Doppler dell'SLC del pacchetto qui e' la frequenza spaziale lungo le righe;
#   ogni bin corrisponde all'impulso il cui LOS ha quella frequenza (stessa relazione di
#   f_D = v . k), e lo zero-Doppler di riferimento e' il bin 0 (centro dello spettro);
# - il banco ha le stesse frazioni del supporto Doppler del preset di riferimento
#   (ipotesi: maschere del preset da bordo a bordo del supporto, 3655 bin);
# - fianchi: due archi espliciti a terra, a Nord e a Sud della base (le linee sono E-O).
BP19 = {
    # modello di profondita' di riferimento di Giza (calcare del Mokattam), dal pacchetto
    "layers": [(2.4, 0.036858797794585174), (None, 0.08074816571554078)],
    "depth": (0.1, 0.05),                        # min, passo [m] sotto la superficie
    "below_base_m": 300.0,                       # profondita' fino a questa quota sotto la base
    "ref_support_bins": 3655, "ref_mask": 1241, "ref_sep": 62, "ref_step": 8,
    "ref_K": 295, "ref_W": 175,
    "pixel_support": 8, "axis_ratio_min": 0.1,
    "support_az_rg": (128, 128), "core_az_rg": (24, 32), "upsample": 100,
    "planted_row_col": (0.73, -0.41), "control_limit_px": 0.05,
    "flank_count": 10, "flank_margin_m": 20.0,  # archi a +-(base/2 + margine) a Nord e Sud
}
BP19_SCALE = [(0, "#000000"), (0.15, "#4d4000"), (0.4, "#a68a00"), (0.7, "#e6c700"), (1, "#ffff00")]  # scala R^2: nero -> giallo


def bp19_layer_integral(z: np.ndarray, layers=BP19["layers"]) -> np.ndarray:
    """focus.layer_integral: I(z) = integrale da 0 a z di ds / lunghezza di fase."""
    out, top = np.zeros_like(z, dtype=np.float64), 0.0
    for bottom, length in layers:
        bottom = np.inf if bottom is None else bottom
        out += np.clip(z - top, 0, bottom - top) / length
        top = bottom
    return out


def bp19_rate(b, t, r, inc_deg) -> np.ndarray:
    """focus.time_variant_rate: G = 4 pi dB / (r sin(theta)) / dt * t, array [K, N]."""
    db = np.empty_like(b); dt = np.empty_like(t)
    db[1:] = np.diff(b, axis=0); db[0] = db[1]
    dt[1:] = np.diff(t, axis=0); dt[0] = dt[1]
    sine = np.sin(np.deg2rad(inc_deg))
    if np.any(dt == 0) or np.any(r <= 0) or np.any(sine <= 0):
        raise ValueError("bp19: incrementi di tempo, range o incidenza non validi")
    return 4 * np.pi * db / (r * sine) / dt * t


def bp19_pixel_support(best: np.ndarray, support: int) -> np.ndarray:
    """focus.pixel_support: tiene i membri di sequenze positive lunghe >= support."""
    n = best.shape[0]
    pos = np.isfinite(best) & (best > 0)
    keep = np.zeros_like(pos)
    for s in range(n - support + 1):
        keep[s:s + support] |= pos[s:s + support].all(axis=0)[None, :]
    return keep


def bp19_geometric_median(points: np.ndarray, tol: float = 1e-9, it: int = 200) -> np.ndarray:
    """primitives.geometric_median (Weiszfeld)."""
    v = np.asarray(points, dtype=np.float64)
    est = np.median(v, axis=0)
    for _ in range(it):
        d = np.linalg.norm(v - est, axis=1)
        if np.any(d < tol):
            return v[np.flatnonzero(d < tol)[0]].copy()
        new = np.sum(v / d[:, None], axis=0) / np.sum(1.0 / d)
        if float(np.linalg.norm(new - est)) < tol:
            return new
        est = new
    return est


def _bp19_polar(angle):
    """exp(j angle) complex64, con l'angolo ridotto in float64 prima del cast."""
    a = torch.remainder(angle, 2 * np.pi).float()
    return torch.polar(torch.ones_like(a), a)


def _bp19_pcc(ref, mov, up: int):
    """primitives.phase_cross_correlation_none a lotti (B, h, w): shift [riga, colonna]
    che allinea mov a ref, con il raffinamento DFT a matrice della versione verificata."""
    B, h, w = ref.shape
    dev = ref.device
    prod = torch.fft.fft2(ref) * torch.fft.fft2(mov).conj()
    i = torch.fft.ifft2(prod).abs().reshape(B, -1).argmax(1)
    sh = torch.stack([i // w, i % w], 1).double()
    size = torch.tensor([h, w], dtype=torch.float64, device=dev)
    sh = torch.where(sh > torch.trunc(size / 2), sh - size, sh)
    sh = torch.round(sh * up) / up
    R = int(np.ceil(up * 1.5))
    d0 = float(np.trunc(R / 2))
    off = d0 - sh * up
    r = torch.arange(R, dtype=torch.float64, device=dev)
    fw = torch.fft.ifftshift(torch.arange(w, dtype=torch.float64, device=dev))
    fh = torch.fft.ifftshift(torch.arange(h, dtype=torch.float64, device=dev))
    kc = _bp19_polar((r[None, :, None] - off[:, 1, None, None]) * fw * (-2 * np.pi / (w * up)))
    tmp = torch.einsum("brw,bhw->brh", kc, prod.conj())
    del kc
    kr = _bp19_polar((r[None, :, None] - off[:, 0, None, None]) * fh * (-2 * np.pi / (h * up)))
    c = torch.einsum("bsh,brh->bsr", kr, tmp).conj()
    j = c.abs().reshape(B, -1).argmax(1)
    return sh + (torch.stack([j // R, j % R], 1).double() - d0) / up


def _bp19_fshift(x, shift):
    """scipy.ndimage.fourier_shift + ifft: x (B, h, w) traslata di shift (B, 2) [riga, col]."""
    B, h, w = x.shape
    fr = torch.fft.fftfreq(h, dtype=torch.float64, device=x.device)
    fc = torch.fft.fftfreq(w, dtype=torch.float64, device=x.device)
    ph = -2 * np.pi * (shift[:, 0, None, None] * fr[:, None] + shift[:, 1, None, None] * fc)
    return torch.fft.ifft2(torch.fft.fft2(x) * _bp19_polar(ph))


def bp19_register(a, b, up: int, core: Tuple[int, int]):
    """measurement.register a lotti: shift dal riquadro intero, poi residuo sul core
    centrale dopo aver traslato b (shift-then-crop). Ritorna [riga, colonna] (B, 2)."""
    coarse = _bp19_pcc(a, b, up)
    al = _bp19_fshift(b, coarse)
    h, w = a.shape[1:]
    r0, c0 = (h - core[0]) // 2, (w - core[1]) // 2
    sl = (slice(None), slice(r0, r0 + core[0]), slice(c0, c0 + core[1]))
    return coarse + _bp19_pcc(a[sl].contiguous(), al[sl].contiguous(), up)


def bp19_fit(q: np.ndarray, spacing, rate: np.ndarray, integral: np.ndarray, W: int,
             P: int, ba_min: Optional[float], chunk: int = 100) -> dict:
    """focus.focus su GPU (float64): q [N, K, 2] pixel (range, azimuth), rate [K, N].
    Per ogni posizione e profondita' la finestra di W vettori con l'R^2 corretto massimo
    (fit a ellisse senza offset dei due componenti in metri), poi supporto di P pixel."""
    dev = "cuda" if USE_GPU else "cpu"
    n, k, _ = q.shape
    f64 = torch.float64
    y = torch.as_tensor(q * np.asarray(spacing, dtype=np.float64)[None, None, :],
                        dtype=f64, device=dev)
    yT = y.transpose(1, 2)                                  # N, 2, K
    Rt = torch.as_tensor(rate.T, dtype=f64, device=dev)     # N, K
    It = torch.as_tensor(integral, dtype=f64, device=dev)

    def win(a):          # somme su tutte le finestre di W vettori, finestra come primo asse
        c = torch.nn.functional.pad(torch.cumsum(a, -1), (1, 0))
        return (c[..., W:] - c[..., :-W]).movedim(-1, 0)

    yy, sm = win((y * y).sum(-1)), win(yT)                  # S, N  /  S, N, 2
    sst = yy - (sm * sm).sum(-1) / W
    fac = (2 * W - 1) / (2 * W - 5)
    D = len(integral)
    best = np.empty((n, D)); starts = np.empty((n, D), np.int32)
    for lo in range(0, D, chunk):
        hi = min(D, lo + chunk)
        ph = Rt[:, None, :] * It[None, lo:hi, None]         # N, d, K
        c, s = torch.cos(ph), torch.sin(ph)
        cc, ss, cs = win(c * c), win(s * s), win(c * s)     # S, N, d
        yc = win(c[:, :, None, :] * yT[:, None])            # S, N, d, 2
        ys = win(s[:, :, None, :] * yT[:, None])
        det = cc * ss - cs * cs
        valid = det > 1e-10 * (cc + ss) ** 2
        inv = torch.where(valid, 1 / det, torch.zeros_like(det))[..., None]
        bc = (ss[..., None] * yc - cs[..., None] * ys) * inv
        bs = (cc[..., None] * ys - cs[..., None] * yc) * inv
        sse = torch.clamp(yy[:, :, None] - (bc * yc + bs * ys).sum(-1), min=0)
        st = sst[:, :, None]
        score = torch.where(st > 0, 1 - sse / st * fac, torch.full_like(sse, -np.inf))
        score = torch.where(valid, score, torch.full_like(score, -np.inf))
        if ba_min is not None:
            aa, bb, ab = (bc * bc).sum(-1), (bs * bs).sum(-1), (bc * bs).sum(-1)
            tr = aa + bb
            disc = torch.sqrt(torch.clamp((aa - bb) ** 2 + 4 * ab ** 2, min=0))
            major, minor = (tr + disc) / 2, torch.clamp((tr - disc) / 2, min=0)
            ba = torch.sqrt(torch.where(major > 0, minor / major, torch.zeros_like(major)))
            score = torch.where(ba > ba_min, score, torch.full_like(score, -np.inf))
        mx, am = score.max(0)          # a parita' il primo (come numpy argmax)
        best[:, lo:hi] = mx.cpu().numpy(); starts[:, lo:hi] = am.cpu().numpy()
    sup = bp19_pixel_support(best, P)
    return {"score": np.where(sup, best, 0).astype(np.float32), "best": best,
            "start": np.where(sup, starts, -1), "support": sup}


def bp19_known_match(vol: np.ndarray, e: np.ndarray, H: float, half: float, known: list,
                     cs: float, n_null: int = 1000, min_shift: float = 30.0,
                     seed: int = 0) -> list:
    """Corrispondenza tra l'R^2 corretto dei cubi (vol [iE, iN, i dalla superficie]) e le
    parti note. Celle di una struttura: cubi il cui centro sta nell'involucro convesso del
    wireframe allargato della sua incertezza (come known_occupancy). Statistica: R^2 medio
    nelle celle. Controlli: le stesse celle (stessa profondita' sotto la superficie)
    traslate in orizzontale di almeno min_shift m, con giro sul bordo della griglia;
    p = (1 + controlli con media >= osservata) / (1 + n_null). Ultima riga: tutte insieme."""
    from scipy.spatial import ConvexHull
    n_col, nz = vol.shape[0], vol.shape[2]
    rng = np.random.default_rng(seed)
    sh = rng.integers(-n_col + 1, n_col, (8 * n_null, 2))
    sh = sh[np.hypot(sh[:, 0], sh[:, 1]) * cs >= min_shift][:n_null]
    rows, cells = [], []
    for s in known:
        p = np.asarray(s["seg"], dtype=float).reshape(-1, 3)
        eq = ConvexHull(p, qhull_options="QJ").equations
        r = float(s["incertezza"]) + cs
        ie = np.flatnonzero((e >= p[:, 0].min() - r) & (e <= p[:, 0].max() + r))
        inn = np.flatnonzero((e >= p[:, 1].min() - r) & (e <= p[:, 1].max() + r))
        EI, NI = (a.ravel() for a in np.meshgrid(ie, inn, indexing="ij"))
        cap = H * (1 - np.maximum(np.abs(e[EI]), np.abs(e[NI])) / half)
        i0 = np.clip(np.floor((cap - p[:, 2].max() - r) / cs), 0, nz).astype(int)
        i1 = np.clip(np.ceil((cap - p[:, 2].min() + r) / cs) + 1, 0, nz).astype(int)
        n = np.maximum(i1 - i0, 0)
        c = np.repeat(np.arange(EI.size), n)
        I = np.repeat(i0, n) + np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
        X = np.column_stack([e[EI[c]], e[NI[c]], cap[c] - (I + 0.5) * cs])
        ok = (X @ eq[:, :3].T + eq[:, 3]).max(axis=1) <= float(s["incertezza"])
        cells.append(np.column_stack([EI[c][ok], NI[c][ok], I[ok]]))
    for k, cl in enumerate(cells + [np.unique(np.concatenate(cells), axis=0)]):
        if not len(cl):
            continue
        obs = float(vol[cl[:, 0], cl[:, 1], cl[:, 2]].mean())
        ctrl = np.array([vol[(cl[:, 0] + a) % n_col, (cl[:, 1] + b) % n_col, cl[:, 2]].mean()
                         for a, b in sh])
        rows.append({"n": k + 1 if k < len(known) else 0,
                     "nome": known[k]["nome"] if k < len(known) else "tutte le parti note",
                     "celle": int(len(cl)), "R2_medio": round(obs, 4),
                     "frazione_positivi": round(float((vol[cl[:, 0], cl[:, 1], cl[:, 2]] > 0).mean()), 3),
                     "R2_medio_controlli": round(float(ctrl.mean()), 4),
                     "percentile": round(float((ctrl < obs).mean() * 100), 1),
                     "p": round((1 + int((ctrl >= obs).sum())) / (1 + len(sh)), 4)})
    return rows


def _bp19_ky(geo: dict, P: np.ndarray, chunk: int = 1024, step: int = 8):
    """Frequenza spaziale lungo le righe [cicli/m] che l'impulso porta nel punto P (relativo
    all'SRP): dalla fase della backprojection, ky = sgn f_c/c (u_az . u_tx + u_az . u_rcv),
    con u = versore dal punto al sensore. Un impulso ogni `step`; a blocchi di punti."""
    n = len(geo["t_pulse"])
    idx = np.unique(np.r_[np.arange(0, n, step), n - 1])
    tx, rc = geo["tx"][idx] - geo["srp"], geo["rcv"][idx] - geo["srp"]
    k = geo["sgn"] * geo["f_c"] / _C
    for i0 in range(0, len(P), chunk):
        p = P[i0:i0 + chunk, None, :]
        ut, ur = tx[None] - p, rc[None] - p
        ky = k * ((ut @ geo["u_az"]) / np.linalg.norm(ut, axis=2) +
                  (ur @ geo["u_az"]) / np.linalg.norm(ur, axis=2))
        yield i0, idx, ky


def bp19_geometry(geo: dict, P: np.ndarray, ky_looks: np.ndarray) -> dict:
    """geometry.from_orbit per l'SLC da backprojection: per ogni punto P (M, 3) e look
    (frequenza ky del bin centrale della coppia) l'istante dell'impulso con quella ky
    (la ky e' monotona lungo l'apertura), poi baseline B sulla tangente perpendicolare
    al LOS di riferimento (ky = 0), range, incidenza locale e G. Array [K, M]."""
    n = len(geo["t_pulse"])
    pos = 0.5 * (geo["tx"] + geo["rcv"]) - geo["srp"]
    K, M = len(ky_looks), len(P)
    out = {a: np.empty((K, M)) for a in ("b", "t", "r", "inc")}
    want = np.r_[0.0, ky_looks]
    for i0, idx, ky in _bp19_ky(geo, P):
        for i in range(ky.shape[0]):
            x, f = (ky[i], idx) if ky[i, -1] > ky[i, 0] else (ky[i, ::-1], idx[::-1])
            if not (np.all(np.diff(x) > 0) and x[0] <= want.min() and want.max() <= x[-1]):
                raise RuntimeError("bp19: look fuori dall'apertura o ky non monotona")
            fi = np.interp(want, x, f)                     # indice d'impulso frazionario
            j = np.minimum(fi.astype(int), n - 2)
            w = (fi - j)[:, None]
            s = pos[j] * (1 - w) + pos[j + 1] * w
            t = geo["t_pulse"][j] * (1 - w[:, 0]) + geo["t_pulse"][j + 1] * w[:, 0]
            v0 = geo["vel"][j[0]] * (1 - w[0]) + geo["vel"][j[0] + 1] * w[0]
            p = P[i0 + i]
            los = p - s[0]; los /= np.linalg.norm(los)
            tan = v0 - (v0 @ los) * los; tan /= np.linalg.norm(tan)
            d = s[1:] - p
            r = np.linalg.norm(d, axis=1)
            out["b"][:, i0 + i] = (s[1:] - s[0]) @ tan
            out["t"][:, i0 + i] = t[1:]
            out["r"][:, i0 + i] = r
            out["inc"][:, i0 + i] = np.degrees(np.arccos(np.clip(d @ geo["up"] / r, -1, 1)))
    out["rate"] = bp19_rate(out["b"], out["t"], out["r"], out["inc"])
    return out


def bp19_bank(geo: dict, P: np.ndarray) -> dict:
    """Banco R/O con le frazioni del supporto Doppler del preset di riferimento. Il
    supporto (in bin dello spettro lungo le righe) e' quello comune a tutti i punti."""
    lo, hi = -np.inf, np.inf
    for _, _, ky in _bp19_ky(geo, P):
        lo = max(lo, float(ky.min(axis=1).max())); hi = min(hi, float(ky.max(axis=1).min()))
    bin_ky = 1.0 / (geo["n_az"] * geo["da"])               # cicli/m per bin
    b_lo, b_hi = int(np.ceil(lo / bin_ky)), int(np.floor(hi / bin_ky))
    S = b_hi - b_lo + 1
    f = S / BP19["ref_support_bins"]
    mask = int(round(BP19["ref_mask"] * f)) | 1
    sep = max(2, 2 * int(round(BP19["ref_sep"] * f / 2)))
    step = max(1, int(round(BP19["ref_step"] * f)))
    K = (S - mask - sep) // step + 1
    span = (K - 1) * step + sep + mask
    mid = b_lo + sep // 2 + mask // 2 + (S - span) // 2 + step * np.arange(K)
    W = int(round(BP19["ref_W"] / BP19["ref_K"] * K))
    return {"support_bins": [b_lo, b_hi], "S": S, "mask": mask, "sep": sep, "step": step,
            "K": K, "W": W, "mid": mid, "ky": mid * bin_ky}


def bp19_measure(slc: np.ndarray, rows: np.ndarray, cols: np.ndarray, bank: dict,
                 logger: logging.Logger, batch: int = 768):
    """measurement.measure su GPU: per ogni coppia (maschere rettangolari di `mask` bin
    centrate a mid -+ sep/2, ricostruzione con IFFT intera) registra la patch di supporto
    attorno a ogni punto (righe, colonne). Ritorna shift [K, M, 2] in pixel nell'ordine
    (range, azimuth), coerenza complessa [K, M] e il controllo con shift iniettato."""
    dev = "cuda" if USE_GPU else "cpu"
    h = np.array(BP19["support_az_rg"]) // 2
    core, up = BP19["core_az_rg"], BP19["upsample"]
    N, n_rg = slc.shape
    if (rows.min() < h[0] or rows.max() + h[0] > N or cols.min() < h[1] or
            cols.max() + h[1] > n_rg):
        raise RuntimeError("bp19: le patch di supporto escono dall'immagine")
    x = torch.from_numpy(np.ascontiguousarray(slc)).to(dev)
    spec = torch.fft.fftshift(torch.fft.fft(x, dim=0), dim=0)
    del x
    mask, sep, mid = bank["mask"], bank["sep"], bank["mid"]

    def image(c):
        lo = N // 2 + int(c) - mask // 2
        m = torch.zeros_like(spec)
        m[lo:lo + mask] = spec[lo:lo + mask]
        return torch.fft.ifft(torch.fft.ifftshift(m, dim=0), dim=0)

    Rt = torch.from_numpy(rows.astype(np.int64)).to(dev)[:, None, None]
    Ct = torch.from_numpy(cols.astype(np.int64)).to(dev)[:, None, None]
    dR = torch.arange(-h[0], h[0], device=dev)[None, :, None]
    dC = torch.arange(-h[1], h[1], device=dev)[None, None, :]

    def frames(img, i0, i1):
        return img[Rt[i0:i1] + dR, Ct[i0:i1] + dC]

    K, M = len(mid), len(rows)
    a = frames(image(mid[K // 2] - sep // 2), 0, 1)
    inj = torch.tensor([BP19["planted_row_col"]], dtype=torch.float64, device=dev)
    rec = bp19_register(a, _bp19_fshift(a, inj), up, core)
    err = float((rec + inj).abs().max())
    control = {"injected_shift_row_col": list(BP19["planted_row_col"]),
               "recovered_alignment_shift_row_col": rec[0].tolist(),
               "max_error_pixels": err, "passed": err <= BP19["control_limit_px"]}
    logger.info(f"  controllo con shift iniettato {BP19['planted_row_col']}: errore "
                f"{err:.4f} px ({'superato' if control['passed'] else 'FALLITO'})")
    if not control["passed"]:
        raise RuntimeError("bp19: controllo con shift iniettato fallito")
    Y = np.empty((K, M, 2), np.float32)
    coh = np.empty((K, M), np.float32)
    t0 = time.time()
    for k, m in enumerate(mid):
        ref, mov = image(m - sep // 2), image(m + sep // 2)
        for i0 in range(0, M, batch):
            i1 = min(M, i0 + batch)
            a, b = frames(ref, i0, i1), frames(mov, i0, i1)
            Y[k, i0:i1] = bp19_register(a, b, up, core).flip(1).cpu().numpy()
            num = (a.conj() * b).sum((1, 2)).abs()
            den = torch.sqrt((a.abs() ** 2).sum((1, 2)) * (b.abs() ** 2).sum((1, 2)))
            coh[k, i0:i1] = (num / den).cpu().numpy()
        del ref, mov
        if (k + 1) % max(1, K // 10) == 0 or k + 1 == K:
            el = time.time() - t0
            logger.info(f"    coppia {k + 1}/{K}  ({el:.0f} s, restano ~{el / (k + 1) * (K - k - 1):.0f} s)")
    del spec
    if dev == "cuda":
        torch.cuda.empty_cache()
    return Y, coh, control


def plot_bp19_line(score: np.ndarray, x: np.ndarray, z: np.ndarray, label: str,
                   out_path: str) -> None:
    """cli.plot_result: R^2 corretto lungo la linea (nero = zero / non supportato)."""
    xe = np.r_[x[0] - (x[1] - x[0]) / 2, (x[:-1] + x[1:]) / 2, x[-1] + (x[-1] - x[-2]) / 2]
    ze = np.r_[z[0] - (z[1] - z[0]) / 2, (z[:-1] + z[1:]) / 2, z[-1] + (z[-1] - z[-2]) / 2]
    cmap = LinearSegmentedColormap.from_list("nero_giallo", BP19_SCALE)
    cmap.set_bad("black")
    fig, ax = plt.subplots(figsize=(16, 6), layout="constrained")
    im = ax.pcolormesh(xe, ze, np.ma.masked_less_equal(score.T, 0), shading="flat", cmap=cmap,
                       vmin=0, vmax=float(score.max()) or 1, rasterized=True)
    ax.set_ylim(float(z[-1]), 0)
    ax.set_xlabel("Est lungo la linea [m]")
    ax.set_ylabel("Profondita' nominale sotto la superficie [m]")
    ax.set_title(label)
    fig.colorbar(im, ax=ax).set_label("R$^2$ corretto (nero = zero / non supportato)")
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def bp19_pyramid(args, cfg: Config, logger: logging.Logger, name: str, out_dir: str) -> dict:
    """Protocollo v1.9 su una piramide: linee E-O di colonne (passo = risoluzione azimuth)
    sulla base; per ogni colonna il pixel dove il radar vede la superficie (layover),
    profondita' misurata dalla superficie. Dati per la pagina 3D (cubi di R^2 corretto).
    Con --bp19_refit rifa solo il fit dalle misure gia' salvate in out_dir."""
    if args.bp19_refit:
        with open(os.path.join(out_dir, "report.json"), encoding="utf-8") as f:
            old = json.load(f)
        meta = {k: old[k] for k in BP19_MEASURE_KEYS}
        meta["durata_misura_s"] = old.get("durata_misura_s", old.get("durata_s"))
        with np.load(os.path.join(out_dir, "result.npz")) as f:
            e, lines = f["east_m"], f["lines"]
        with np.load(os.path.join(out_dir, "bank.npz")) as f:
            q = f["q"]
        with np.load(os.path.join(out_dir, "geometry.npz")) as f:
            rate = f["rate"]
        logger.info(f"bp19 {name}: nuovo fit dalle misure salvate in {out_dir} ({len(lines)} "
                    f"linee, K={meta['banco']['K']}, W={meta['banco']['W']})")
        return bp19_volume(logger, name, out_dir, meta, e, lines, q, rate)
    side, H = PYRAMIDS[name][2], PYRAMIDS[name][3]
    t_start = time.time()
    # scena: deve contenere i fianchi (a +-(lato/2 + margine) N/S, in diagonale nell'immagine)
    # con le patch (128 px ~ 17 m az, 27 m rg) e, in range, il layover dell'apice
    rad = (side / 2 + BP19["flank_margin_m"]) * np.sqrt(2)
    slc, geo = backproject_cphd(args.cphd, cfg, logger, args.pulse_start, args.n_pulses,
                                2 * (rad + 20), 2 * (rad + 30) + 1.3 * H, sign_mult=args.bp_sign,
                                center_on=name, az_os=args.bp19_os)
    px, half, H = _pyr_px_fn(geo, name)
    lat0, lon0, _ = _ecf2geo(geo["srp"])
    east = np.array([-np.sin(lon0), np.cos(lon0), 0.0])
    north = np.array([-np.sin(lat0) * np.cos(lon0), -np.sin(lat0) * np.sin(lon0), np.cos(lat0)])
    la, lo = PYRAMIDS[name][:2]
    ctr = _geo2ecf(np.radians(la), np.radians(lo), _ecf2geo(geo["srp"])[2]) - geo["srp"]
    n_col = int(round(2 * half / geo["res_az"])) + 1
    cs = 2 * half / (n_col - 1)
    e = np.linspace(-half, half, n_col)
    lines = np.arange(n_col)
    if args.bp19_lines > 0:                              # prova: solo le linee centrali
        lines = lines[max(0, n_col // 2 - args.bp19_lines // 2):][:args.bp19_lines]
    EE = np.tile(e, len(lines))                          # target: linea per linea, da Ovest a Est
    NN = np.repeat(e[lines], n_col)
    cap = H * (1 - np.maximum(np.abs(EE), np.abs(NN)) / half)
    rc = np.array([px(a, b, c) for a, b, c in zip(EE, NN, cap)])
    J, fn = BP19["flank_count"], half + BP19["flank_margin_m"]
    fe = np.linspace(-half, half, J)
    frc = np.array([px(a, s * fn, 0.0) for s in (1, -1) for a in fe])   # Nord, poi Sud
    rows = np.floor(np.r_[rc[:, 0], frc[:, 0]] + 0.5).astype(int)
    cols = np.floor(np.r_[rc[:, 1], frc[:, 1]] + 0.5).astype(int)
    P = ctr + EE[:, None] * east + NN[:, None] * north + cap[:, None] * geo["up"]
    bank = bp19_bank(geo, P)
    logger.info(f"bp19 {name}: {len(lines)} linee E-O x {n_col} colonne (passo {cs:.2f} m), "
                f"fianchi {J}+{J} a N/S a {fn:.0f} m dal centro; supporto Doppler {bank['S']} bin, "
                f"K={bank['K']} coppie, maschera {bank['mask']}, separazione {bank['sep']}, "
                f"passo {bank['step']} bin, W={bank['W']}, P={BP19['pixel_support']}")
    g = bp19_geometry(geo, P, bank["ky"])
    logger.info(f"  geometria: t {g['t'].min():.3f}..{g['t'].max():.3f} s dall'inizio della raccolta, "
                f"B {g['b'].min():.0f}..{g['b'].max():.0f} m, r {g['r'].mean() / 1e3:.1f} km, "
                f"incidenza {g['inc'].min():.1f}..{g['inc'].max():.1f} deg")
    Y, coh, control = bp19_measure(slc, rows, cols, bank, logger)
    n_t = len(EE)
    med = np.array([[bp19_geometric_median(Y[k, n_t + s * J:n_t + (s + 1) * J]) for s in (0, 1)]
                    for k in range(bank["K"])])
    corr = 0.5 * (med[:, 0] + med[:, 1])                 # G_k = media delle mediane dei fianchi
    q = (Y[:, :n_t] - corr[:, None, :]).transpose(1, 0, 2)
    spacing = np.array([geo["dr"], geo["da"]])           # metri per pixel (range, azimuth)
    os.makedirs(out_dir, exist_ok=True)
    np.savez_compressed(os.path.join(out_dir, "bank.npz"), q=q, raw_shift_px=Y,
                        coherence=coh.astype(np.float16), flank_geometric_median_px=med,
                        arc_correction_px=corr, midpoint_bins=bank["mid"], ky_looks=bank["ky"],
                        rows=rows, cols=cols, east_m=EE, north_m=NN, surface_m=cap)
    np.savez_compressed(os.path.join(out_dir, "geometry.npz"), rate=g["rate"],
                        **{k: g[k].astype(np.float32) for k in ("b", "t", "r", "inc")})
    meta = {"protocollo": "Biondi Protocol v1.9.0 (porting GPU)", "piramide": name,
            "cphd": os.path.basename(args.cphd), "impulsi": [args.pulse_start, args.n_pulses],
            "sovracampionamento_azimuth": args.bp19_os, "linee": int(len(lines)),
            "colonne_per_linea": int(n_col), "passo_colonne_m": round(cs, 3),
            "spacing_range_azimuth_m": spacing.tolist(),
            "banco": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                      for k, v in bank.items() if k not in ("mid", "ky")},
            "midpoint_bins": [int(bank["mid"][0]), int(bank["mid"][-1])],
            "controllo": control, "registrazioni": int(bank["K"] * len(rows)),
            "coerenza_mediana": round(float(np.median(coh[:, :n_t])), 4),
            "durata_misura_s": round(time.time() - t_start, 1)}
    del slc, Y
    return bp19_volume(logger, name, out_dir, meta, e, lines, q, g["rate"])


BP19_MEASURE_KEYS = ("protocollo", "piramide", "cphd", "impulsi", "sovracampionamento_azimuth",
                     "linee", "colonne_per_linea", "passo_colonne_m", "spacing_range_azimuth_m",
                     "banco", "midpoint_bins", "controllo", "registrazioni", "coerenza_mediana")


def bp19_volume(logger: logging.Logger, name: str, out_dir: str, meta: dict, e: np.ndarray,
                lines: np.ndarray, q: np.ndarray, rate: np.ndarray) -> dict:
    """Fit linea per linea dalla superficie fino a below_base_m sotto la base (ogni linea
    fino alla sua quota massima + below_base_m), cubi di R^2 massimo, corrispondenza con le
    parti note, file e dati della pagina."""
    t_start = time.time()
    side, H = PYRAMIDS[name][2], PYRAMIDS[name][3]
    half, n_col = side / 2, len(e)
    cs = 2 * half / (n_col - 1)                          # come in bp19_pyramid
    bank = meta["banco"]
    spacing = np.array(meta["spacing_range_azimuth_m"])
    z0, dz = BP19["depth"]
    zb = BP19["below_base_m"]
    depths = np.round(z0 + np.arange(int(np.floor((H + zb - z0) / dz + 1e-9)) + 1) * dz, 12)
    integral = bp19_layer_integral(depths)
    nz = int(np.floor((H + zb) / cs)) + 1
    cid = np.floor(depths / cs).astype(int)
    vol = np.zeros((n_col, n_col, nz), np.float32)       # [iE, iN, cubo dalla superficie]
    kmax = np.zeros((n_col, n_col, nz), np.uint16)       # indice in depths del massimo del cubo
    mid_line = len(lines) // 2
    logger.info(f"  fit: profondita' {z0}..{depths[-1]:.2f} m (passo {dz}), fino a {zb:.0f} m "
                f"sotto la base; {nz} cubi per colonna")
    for li, iN in enumerate(lines):
        nd = int(np.searchsorted(depths, H * (1 - abs(e[iN]) / half) + zb, side="right"))
        sl = slice(li * n_col, (li + 1) * n_col)
        r = bp19_fit(q[sl], spacing, rate[:, sl], integral[:nd], bank["W"],
                     BP19["pixel_support"], BP19["axis_ratio_min"])
        st = np.flatnonzero(np.r_[True, np.diff(cid[:nd]) > 0])
        vol[:, iN, cid[st]] = np.maximum.reduceat(r["score"], st, axis=1)
        for a, b in zip(st, np.r_[st[1:], nd]):
            kmax[:, iN, cid[a]] = a + r["score"][:, a:b].argmax(1)
        if li == mid_line:
            central = dict(r, iN=int(iN), nd=nd)
        if (li + 1) % max(1, len(lines) // 10) == 0 or li + 1 == len(lines):
            el = time.time() - t_start
            logger.info(f"    linea {li + 1}/{len(lines)}  ({el:.0f} s, restano "
                        f"~{el / (li + 1) * (len(lines) - li - 1):.0f} s)")
    cap = H * (1 - np.maximum(np.abs(e)[:, None], np.abs(e)[None, :]) / half)   # [iE, iN]
    valid = np.arange(nz)[None, None, :] * cs < cap[..., None] + zb   # fino a zb sotto la base
    vol[~valid] = 0
    lab = (f"{name.capitalize()} - linea E-O a N {e[central['iN']]:+.1f} m - v1.9 - "
           f"K={bank['K']}, W={bank['W']}, P={BP19['pixel_support']}")
    plot_bp19_line(central["score"], e, depths[:central["nd"]], lab,
                   os.path.join(out_dir, "tomogramma_linea.png"))
    np.savez_compressed(os.path.join(out_dir, "result.npz"), cube_score=vol, cube_m=cs,
                        cube_depth_max_m=np.where(vol > 0, depths[kmax], 0).astype(np.float32),
                        depth_m=depths, layer_integral=integral, east_m=e, lines=lines,
                        line_iN=central["iN"], line_score=central["score"],
                        line_best=central["best"], line_start=central["start"],
                        line_support=central["support"])
    done = vol[:, lines]
    pk = np.unravel_index(int(done.argmax()), done.shape)
    peak = None if done.max() <= 0 else {
        "est_m": round(float(e[pk[0]]), 2), "nord_m": round(float(e[lines[pk[1]]]), 2),
        "profondita_cubo_m": [round(float(pk[2] * cs), 2), round(float((pk[2] + 1) * cs), 2)],
        "profondita_massimo_m": round(float(depths[kmax[pk[0], lines[pk[1]], pk[2]]]), 2),
        "quota_z_m": round(float(cap[pk[0], lines[pk[1]]] - (pk[2] + 0.5) * cs), 2),
        "R2_corretto": round(float(done.max()), 4)}
    ks = load_known_structures(logger).get(name, [])
    match = bp19_known_match(vol, e, H, half, ks, cs) if ks and len(lines) == n_col else []
    for m in match:
        logger.info(f"  parte nota {m['n']:2d} {m['nome'][:48]:48s}: {m['celle']:5d} celle, "
                    f"R2 medio {m['R2_medio']:.3f} (controlli {m['R2_medio_controlli']:.3f}, "
                    f"percentile {m['percentile']:.1f}), p = {m['p']:.4f}")
    report = dict(meta, config={k: v for k, v in BP19.items()},
                  copertura={"posizioni": int(n_col * len(lines)), "coppie": int(bank["K"]),
                             "finestre": int(bank["K"] - bank["W"] + 1),
                             "profondita": int(len(depths))},
                  cubi_positivi_frazione=round(float((done[valid[:, lines]] > 0).mean()), 4),
                  picco=peak,
                  corrispondenza_parti_note=match,
                  durata_fit_s=round(time.time() - t_start, 1),
                  nota="i punteggi di profondita' non provano vuoti, penetrazione o profondita' "
                       "assolute (README del pacchetto); l'R^2 corretto non e' una densita'")
    with open(os.path.join(out_dir, "report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    logger.info(f"  {name}: cubi positivi {report['cubi_positivi_frazione'] * 100:.1f} %, "
                f"picco {peak}; fit {report['durata_fit_s']:.0f} s")
    with open(os.path.join(out_dir, "tomogramma_linea.png"), "rb") as f:
        png = base64.b64encode(f.read()).decode("ascii")
    known = [{"n": i + 1, "name": s["nome"], "cat": s["categoria"], "color": s["colore"],
              "fonte": s["fonte"], "inc": s["incertezza"], "schem": s["schematica"],
              "nota": s["nota"], "c": s["c"], "seg": s["seg"]} for i, s in enumerate(ks)]
    vox = np.round(65535 * np.clip(vol, 0, 1)).astype(np.uint16).ravel()   # R^2 a 16 bit
    idx = np.flatnonzero(vox).astype(np.uint32)          # solo i cubi con R^2 > 0
    b64 = lambda a: base64.b64encode(np.ascontiguousarray(a, a.dtype.newbyteorder("<")).tobytes()).decode("ascii")
    return {"name": name, "label": name.capitalize(), "side": round(2 * half, 1), "H": H,
            "n_col": n_col, "cs": cs, "nz": nz, "zb": zb, "lines": len(lines),
            "z0": z0, "dz": dz, "vox_idx_b64": b64(idx), "vox_val_b64": b64(vox[idx]),
            "vox_dep_b64": b64(kmax.ravel()[idx]),
            "png_b64": png, "known": known, "match": match, "K": bank["K"], "W": bank["W"],
            "P": BP19["pixel_support"], "mask": bank["mask"], "sep": bank["sep"],
            "step": bank["step"], "S": bank["S"], "os": meta["sovracampionamento_azimuth"],
            "supp": list(BP19["support_az_rg"]), "core": list(BP19["core_az_rg"]),
            "ctrl": round(meta["controllo"]["max_error_pixels"], 4),
            "pos": report["cubi_positivi_frazione"], "peak": peak,
            "line_n": round(float(e[central["iN"]]), 1)}


# =============================================================================
#  ECOGRAFIA (--eco): test dell'ombra e sezione verticale con rifrazione
# =============================================================================

ECO_GRAZ = 75.0      # incidenza locale [gradi] oltre la quale una faccia conta come zona scura


def _eco_setup(filepath: str, name: str, pulse_start: int, n_pulses: int,
               logger: logging.Logger, sign_mult: int = -1) -> dict:
    """Blocco CPHD (PVP) e geometria dell'SLC centrata sulla piramide, come per --voxels:
    i pixel sono quelli di orientamento.png."""
    rd = CPHDReader(filepath)
    d = rd.cphd_meta.to_dict()
    ch = d["Data"]["Channels"][0]
    N, total = int(ch["NumSamples"]), int(ch["NumVectors"])
    p0 = max(0, int(pulse_start))
    n = min(total, p0 + int(n_pulses)) - p0
    pv = {k: np.asarray(rd.read_pvp_variable(k, 0, slice(p0, p0 + n)), dtype=np.float64)
          for k in ("TxPos", "RcvPos", "TxVel", "SRPPos", "AmpSF", "SC0", "SCSS")}
    sc0, d_sc = float(pv["SC0"][0]), float(pv["SCSS"][0])
    bw = N * d_sc
    pv["fc"] = sc0 + 0.5 * bw
    side, H = PYRAMIDS[name][2], PYRAMIDS[name][3]
    geo = _bp_geometry(pv, n, N, d_sc, bw, side * 1.6, side * 1.6 + H * 1.2, 4_000_000, logger)
    px, _, H = _pyr_px_fn(geo, name)
    r0, c0 = px(0.0, 0.0)
    geo["cy"] += (r0 - geo["n_az"] // 2) * geo["da"]
    geo["cx"] += (c0 - geo["n_rg"] // 2) * geo["dr"] - H / np.tan(geo["inc"]) / 2
    # versori Est, Nord, Su e centro della base nel riferimento (u_rg, u_az, up) dell'SLC
    lat0, lon0, h0 = _ecf2geo(geo["srp"])
    east = np.array([-np.sin(lon0), np.cos(lon0), 0.0])
    north = np.array([-np.sin(lat0) * np.cos(lon0), -np.sin(lat0) * np.sin(lon0), np.cos(lat0)])
    B = np.stack([geo["u_rg"], geo["u_az"], geo["up"]])
    c = _geo2ecf(np.radians(PYRAMIDS[name][0]), np.radians(PYRAMIDS[name][1]), h0) - geo["srp"]
    inc = geo["inc"]
    return dict(rd=rd, pv=pv, N=N, n=n, p0=p0, d_sc=d_sc, bw=bw, sc0=sc0, geo=geo, B=B,
                sgn=int(d["Global"].get("SGN", 1)) * int(sign_mult),   # segno di fase usato
                ue=B @ east, un=B @ north,
                uu=B @ np.cross(east, north), cen=np.array([c @ geo["u_rg"], c @ geo["u_az"], 0.0]),
                k=np.array([np.sin(inc), 0.0, -np.cos(inc)]))   # sensore -> scena (locale)


def _eco_px(S: dict, e, n, z):
    """(Est, Nord, z) dal centro della base -> (riga, colonna) dell'SLC (vettoriale, layover)."""
    geo, ue, un, cen = S["geo"], S["ue"], S["un"], S["cen"]
    x = cen[0] + e * ue[0] + n * un[0] - z / np.tan(geo["inc"])
    y = cen[1] + e * ue[1] + n * un[1]
    return (geo["n_az"] // 2 + (y - geo["cy"]) / geo["da"],
            geo["n_rg"] // 2 + (x - geo["cx"]) / geo["dr"])


def _eco_faces(S: dict, half: float, H: float) -> dict:
    """Normali esterne delle 4 facce nel riferimento locale."""
    b = np.arctan2(H, half)
    out = {"Est": (1, 0), "Ovest": (-1, 0), "Nord": (0, 1), "Sud": (0, -1)}
    return {f: np.sin(b) * (h[0] * S["ue"] + h[1] * S["un"]) + np.cos(b) * S["uu"]
            for f, h in out.items()}


def _eco_height(e, n, half: float, H: float):
    return np.clip(H * (1 - np.maximum(np.abs(e), np.abs(n)) / half), 0, None)


def _eco_refract(k: np.ndarray, nrm: np.ndarray, nr: float):
    """Direzione rifratta (Snell vettoriale) e angoli di incidenza / rifrazione [gradi]."""
    ci = -k @ nrm
    eta = 1.0 / nr
    ct = np.sqrt(1 - eta ** 2 * (1 - ci ** 2))
    return eta * k + (eta * ci - ct) * nrm, np.degrees(np.arccos(ci)), np.degrees(np.arccos(ct))


def _eco_bp(S: dict, geo: dict, logger: logging.Logger) -> np.ndarray:
    """Backprojection SENZA normalizzazione (scala assoluta comune a scena e rumore)."""
    dev = "cuda" if USE_GPU else "cpu"
    f_c = S["sc0"] + (S["N"] // 2) * S["d_sc"]
    img = _bp_run(S["rd"], S["p0"], S["n"], S["N"], geo, S["pv"], f_c, S["d_sc"],
                  S["sgn"], 8, 64, dev, logger)
    out = img.cpu().numpy()
    del img
    if dev == "cuda":
        torch.cuda.empty_cache()
    return out


def eco_shadow_test(S: dict, name: str, nr: float, out_dir: str, logger: logging.Logger):
    """Test dell'ombra (A). Le facce con incidenza locale < ECO_GRAZ danno eco di superficie
    forte; i pixel dove ricadono solo facce quasi radenti sono la zona scura. Per ogni pixel
    scuro L_min = cammino minimo nella pietra che un'eco interna (ingresso da una faccia
    illuminata, raggio rifratto) deve fare per finirci. Un'eco interna vera cala di decine
    di dB per metro; il rumore si misura con una patch fuori scena (stessa backprojection).
    Ritorna (slc, potenza media del rumore)."""
    geo, k = S["geo"], S["k"]
    _, half, H = _pyr_px_fn(geo, name)
    n_az, n_rg = geo["n_az"], geo["n_rg"]
    F = _eco_faces(S, half, H)
    cg = np.cos(np.radians(ECO_GRAZ))
    # superfici campionate a 0.25 m: suolo + facce
    ext = 1.6 * half + H
    g = np.arange(-ext, ext, 0.25)
    E, Nn = (a.ravel() for a in np.meshgrid(g, g, indexing="ij"))
    Z = _eco_height(E, Nn, half, H)
    on_pyr = Z > 0
    fe = np.where(np.abs(E) >= np.abs(Nn), np.where(E > 0, "Est", "Ovest"),
                  np.where(Nn > 0, "Nord", "Sud"))
    lit = np.zeros(E.size, dtype=bool)
    for f, nrm in F.items():
        lit[on_pyr & (fe == f)] = (-k @ nrm) > cg
    # suolo: illuminato se il raggio verso il sensore non incontra la piramide
    kE = np.array([k @ S["ue"], k @ S["un"], k @ S["uu"]])
    gnd = ~on_pyr
    vis = np.ones(int(gnd.sum()), dtype=bool)
    Eg, Ng = E[gnd], Nn[gnd]
    for s in np.arange(0.5, H / (-kE[2]) + 1, 0.5):
        vis &= ~(-s * kE[2] < _eco_height(Eg - s * kE[0], Ng - s * kE[1], half, H))
    lit[gnd] = vis
    rr, cc = _eco_px(S, E, Nn, Z)
    ri, ci = np.round(rr).astype(int), np.round(cc).astype(int)
    ok = (ri >= 0) & (ri < n_az) & (ci >= 0) & (ci < n_rg)
    n_lit = np.bincount(ri[ok & lit] * n_rg + ci[ok & lit], minlength=n_az * n_rg).reshape(n_az, n_rg)
    n_all = np.bincount(ri[ok] * n_rg + ci[ok], minlength=n_az * n_rg).reshape(n_az, n_rg)
    dark = (n_all > 0) & (n_lit == 0)
    # L_min: ingressi dalle facce illuminate, raggio rifratto, eco spostata in range di n L
    Lmin = np.full((n_az, n_rg), np.inf)
    ent = np.nonzero(on_pyr & lit)[0][::4]
    rep = {"pyr": name, "inc_deg": round(float(np.degrees(geo["inc"])), 2), "eps": round(nr ** 2, 2),
           "vista_da_azimut_deg": round(float(np.degrees(np.arctan2(-kE[0], -kE[1])) % 360), 1),
           "facce": {}}
    for f, nrm in F.items():
        ci_f = float(-k @ nrm)
        rep["facce"][f] = {"incidenza_locale_deg": round(float(np.degrees(np.arccos(np.clip(ci_f, -1, 1)))), 1),
                           "stato": "illuminata" if ci_f > cg else ("radente" if ci_f > 0 else "in ombra")}
        if ci_f <= cg:
            continue
        t, _, tt = _eco_refract(k, nrm, nr)
        rep["facce"][f]["rifrazione_deg"] = round(float(tt), 1)
        sel = ent[fe[ent] == f]
        r0, c0 = _eco_px(S, E[sel], Nn[sel], Z[sel])
        tE = np.array([t @ S["ue"], t @ S["un"], t @ S["uu"]])
        for L in np.arange(0.0, 2.2 * half, 0.5):
            inside = Z[sel] + L * tE[2] < _eco_height(E[sel] + L * tE[0], Nn[sel] + L * tE[1], half, H) + 1e-6
            ri = np.round(r0 + L * t[1] / geo["da"]).astype(int)
            ci = np.round(c0 + nr * L / np.sin(geo["inc"]) / geo["dr"]).astype(int)
            m = inside & (ri >= 0) & (ri < n_az) & (ci >= 0) & (ci < n_rg)
            np.minimum.at(Lmin, (ri[m], ci[m]), L)
    rep["eco_interna"] = [{"L_m": L, "cammino_equiv_m": round(nr * L, 2),
                           "ritardo_extra_ns": round(2 * nr * L / _C * 1e9, 1),
                           "spost_ground_range_m": round(nr * L / np.sin(geo["inc"]), 1),
                           "spost_px": round(nr * L / np.sin(geo["inc"]) / geo["dr"], 1)}
                          for L in (5.0, 10.0, 20.0)]
    for tand in (0.01, 0.05):
        rep[f"attenuazione_2vie_dB_m_tand{tand}"] = round(2 * 8.686 * np.pi * S["pv"]["fc"] / _C * nr * tand, 1)

    # SLC e patch di rumore: il ritardo piu' lontano dal centro con profilo in range a fondo
    N, n = S["N"], S["n"]
    Sx = S["rd"].read_chip(slice(S["p0"] + n // 2, S["p0"] + n // 2 + 256), slice(0, N))
    prof = (np.abs(np.fft.fftshift(np.fft.ifft(Sx * np.hamming(N)[None, :], axis=1), axes=1)) ** 2).mean(0)
    r_sl = (np.arange(N) - N // 2) / (N * S["d_sc"]) * _C / 2
    pdb = 10 * np.log10(np.convolve(prof, np.ones(64) / 64, "same") + 1e-30)
    cand = np.nonzero(pdb < np.percentile(pdb, 5) + 1.0)[0]
    cand = cand[np.abs(r_sl[cand]) < 0.9 * np.abs(r_sl).max()]
    i_n = cand[np.argmax(np.abs(r_sl[cand]))]
    slc = _eco_bp(S, geo, logger)
    noise = _eco_bp(S, dict(geo, cx=float(r_sl[i_n] / np.sin(geo["inc"])), cy=0.0, n_az=101, n_rg=301), logger)
    pn = float((np.abs(noise) ** 2).mean())
    rep["rumore_slant_m"] = round(float(r_sl[i_n]), 0)
    pw = np.abs(slc) ** 2 / pn
    rows = []
    for lo, hi in zip([0, 5, 10, 20, 40, 80], [5, 10, 20, 40, 80, 400]):
        m = dark & (Lmin >= lo) & (Lmin < hi)
        if m.sum() >= 20:
            rows.append((lo, hi, int(m.sum()), float(10 * np.log10(pw[m].mean()))))
    db_of = lambda m: round(float(10 * np.log10(pw[m].mean())), 1) if m.sum() > 20 else None
    rep["pixel_zona_scura"] = int(dark.sum())
    rep["illuminato_dB"] = db_of(n_lit > 0)
    rep["zona_scura_dB"] = db_of(dark)
    rep["zona_scura_irraggiungibile_dB"] = db_of(dark & ~np.isfinite(Lmin))
    rep["zona_scura_per_L_min"] = [dict(L=f"{lo}-{hi} m", pixel=c, dB_sul_rumore=round(v, 1))
                                   for lo, hi, c, v in rows]
    logger.info("  test dell'ombra: " + json.dumps(rep, ensure_ascii=True))
    with open(os.path.join(out_dir, "test_ombra.json"), "w") as f:
        json.dump(rep, f, indent=1)

    fig, ax = plt.subplots(1, 3, figsize=(20, 5.5), gridspec_kw={"width_ratios": [1.4, 1.4, 1]})
    asp = geo["da"] / geo["dr"]
    ax[0].imshow(10 * np.log10(pw + 1e-12), cmap="gray", vmin=-5, vmax=40, aspect=asp)
    ax[0].contour(dark, levels=[0.5], colors="r", linewidths=0.8)
    ax[0].set_title(f"{name}: SLC [dB sul rumore], zona scura in rosso")
    im = ax[1].imshow(np.where(dark & np.isfinite(Lmin), Lmin, np.nan), cmap="viridis", aspect=asp)
    plt.colorbar(im, ax=ax[1], label="L_min [m]")
    ax[1].set_title(f"zona scura: cammino minimo nella pietra (eps = {nr ** 2:g})")
    for x in ax[:2]:
        x.set_xlabel("range (pixel)"); x.set_ylabel("azimuth (pixel)")
    sel = dark & np.isfinite(Lmin)
    ax[2].plot(Lmin[sel], 10 * np.log10(pw[sel] + 1e-12), ",", alpha=0.3, color="0.5")
    if rows:
        ax[2].plot([(r[0] + min(r[1], Lmin[sel].max())) / 2 for r in rows], [r[3] for r in rows],
                   "o-", color="C1", label="media per intervallo")
    ax[2].axhline(0, color="k", ls="--", label="rumore (patch fuori scena)")
    if rep["illuminato_dB"] is not None:
        ax[2].axhline(rep["illuminato_dB"], color="C2", ls=":", label="superfici illuminate")
    ax[2].set_xlabel("L_min [m]"); ax[2].set_ylabel("potenza [dB sul rumore]")
    ax[2].set_title("un'eco interna calerebbe di decine di dB/m"); ax[2].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(out_dir, "test_ombra.png"), dpi=100); plt.close(fig)
    return slc, pn


def eco_pick_line(slc: np.ndarray, S: dict, name: str):
    """Finestra con l'immagine di orientamento.png: due clic = inizio e fine della retta."""
    plt.switch_backend("TkAgg")
    geo = S["geo"]
    px, half, _ = _pyr_px_fn(geo, name)
    a = np.abs(slc)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.imshow(20 * np.log10(a / np.percentile(a, 99.9) + 1e-6), cmap="gray", vmin=-40, vmax=0,
              aspect=geo["da"] / geo["dr"])
    cor = {"SO": (-half, -half), "SE": (half, -half), "NE": (half, half), "NO": (-half, half)}
    pts = [px(e, n) for e, n in cor.values()]
    ax.plot([p[1] for p in pts + pts[:1]], [p[0] for p in pts + pts[:1]], "r-", lw=1.2)
    for lab, p in zip(cor, pts):
        ax.annotate(lab, (p[1], p[0]), color="yellow", fontsize=12, weight="bold")
    ax.set_title("clicca i 2 punti della retta (inizio, fine)")
    ax.set_xlabel("range (pixel)"); ax.set_ylabel("azimuth (pixel)")
    p = plt.ginput(2, timeout=0)
    plt.close(fig)
    plt.switch_backend("Agg")
    if len(p) < 2:
        raise SystemExit("--eco_pick: servono 2 punti")
    return p[0][0], p[0][1], p[1][0], p[1][1]


def _eco_pix_to_en(S: dict, name: str, x_px: float, y_px: float):
    """(colonna range, riga azimuth) al suolo -> (Est, Nord) dal centro della base."""
    r0, c0 = _eco_px(S, 0.0, 0.0, 0.0)
    r1, c1 = _eco_px(S, 1.0, 0.0, 0.0)
    r2, c2 = _eco_px(S, 0.0, 1.0, 0.0)
    return np.linalg.solve(np.array([[c1 - c0, c2 - c0], [r1 - r0, r2 - r0]]),
                           np.array([x_px - c0, y_px - r0]))


def eco_section(S: dict, name: str, line, step: float, z_lo: float, nr: float,
                logger: logging.Logger, progress=None) -> dict:
    """Ecografia interna (B): backprojection (SAFT) sul piano verticale sotto la retta
    (x0, y0, x1, y1 in pixel di orientamento.png, al suolo). Punto in aria: ritardo in spazio
    libero. Punto nella pietra: per ogni faccia illuminata f il raggio rifratto t_f (Snell,
    indice nr, onda piana per impulso) e' ripercorso all'indietro; in un poliedro convesso esce
    dal piano piu' vicino, e se quel piano e' f l'ingresso E e' valido: ritardo = |T - E| + nr L.
    Punti non raggiungibili da nessuna faccia illuminata restano a zero.
    progress(impulsi fatti, impulsi totali): avanzamento (pagina dell'ecografia)."""
    dev = "cuda" if USE_GPU else "cpu"
    geo, ue, un, uu, cen = S["geo"], S["ue"], S["un"], S["uu"], S["cen"]
    _, half, H = _pyr_px_fn(geo, name)
    e0, n0 = _eco_pix_to_en(S, name, line[0], line[1])
    e1, n1 = _eco_pix_to_en(S, name, line[2], line[3])
    Ltot = float(np.hypot(e1 - e0, n1 - n0))
    logger.info(f"  retta: ({line[0]:.0f}, {line[1]:.0f}) -> ({line[2]:.0f}, {line[3]:.0f}) px = "
                f"(E {e0:+.1f}, N {n0:+.1f}) -> (E {e1:+.1f}, N {n1:+.1f}) m, lunga {Ltot:.1f} m")
    s = np.arange(0.0, Ltot + 1e-9, step)
    z = np.arange(z_lo, H + 5.0 + 1e-9, step)
    Sg, Zg = np.meshgrid(s, z, indexing="ij")
    Eg, Ng = e0 + (e1 - e0) * Sg / Ltot, n0 + (n1 - n0) * Sg / Ltot
    inside = (Zg < _eco_height(Eg, Ng, half, H)) & (Zg >= 0)
    fl = lambda a: torch.tensor(np.asarray(a, dtype=np.float32), device=dev)
    Pt = fl(cen[None, :] + Eg.ravel()[:, None] * ue + Ng.ravel()[:, None] * un
            + Zg.ravel()[:, None] * uu)
    P2 = (Pt * Pt).sum(1)
    ins = torch.tensor(inside.ravel(), device=dev)
    F = _eco_faces(S, half, H)
    names = list(F)
    dP = {g: (Pt - fl(cen + H * uu)) @ fl(v) for g, v in F.items()}   # < 0 dentro
    srp, B, pv = geo["srp"], S["B"], S["pv"]

    def unit_local(pos):
        t = pos - srp
        tn = np.linalg.norm(t, axis=1)
        return tn, (t / tn[:, None]) @ B.T                # versore SRP -> sensore (locale)

    def dpath(tn, u, Q, Q2):
        num = Q2 - 2.0 * tn * (Q @ u)
        return num / (torch.sqrt(num + tn * tn) + tn)     # |T - Q| - |T - SRP|

    def el_path(tn, u_np):
        k = -u_np
        u = fl(u_np)
        best_L = torch.full((Pt.shape[0],), float("inf"), device=dev)
        best_E = torch.zeros_like(Pt)
        for f, nv in F.items():
            if -k @ nv <= 0:
                continue                                  # faccia non illuminata
            tv, _, _ = _eco_refract(k, nv, nr)
            Lx = torch.full_like(best_L, float("inf"))
            arg = torch.full(best_L.shape, -1, dtype=torch.long, device=dev)
            for gi, g in enumerate(names):
                ng_t = float(F[g] @ tv)
                if ng_t >= 0:
                    continue
                Lg = dP[g] / ng_t
                closer = Lg < Lx
                Lx = torch.where(closer, Lg, Lx)
                arg = torch.where(closer, torch.full_like(arg, gi), arg)
            better = (arg == names.index(f)) & (Lx >= 0) & (Lx < best_L)
            best_L = torch.where(better, Lx, best_L)
            best_E = torch.where(better[:, None], Pt - Lx[:, None] * fl(tv), best_E)
        reach = torch.isfinite(best_L)
        el = dpath(tn, u, best_E, (best_E * best_E).sum(1)) + nr * torch.where(reach, best_L, 0.0)
        return torch.where(ins & reach, el, dpath(tn, u, Pt, P2)), ins & ~reach, best_L

    tn_tx, u_tx = unit_local(pv["TxPos"])
    tn_rx, u_rx = unit_local(pv["RcvPos"])
    N, n, M = S["N"], S["n"], S["N"] * 8
    h = N // 2
    win = torch.from_numpy(np.hamming(N).astype(np.float32)).to(dev)
    img = torch.zeros(Pt.shape[0], dtype=torch.complex64, device=dev)
    f_c = S["sc0"] + (N // 2) * S["d_sc"]
    bins_per_m, cyc_per_m, ang_k = S["d_sc"] * M / _C, f_c / _C, -S["sgn"] * 2.0 * np.pi
    amp = pv["AmpSF"]
    Lm = None
    t0 = time.time()
    for c0 in range(0, n, 64):
        c1 = min(n, c0 + 64)
        Sx = S["rd"].read_chip(slice(S["p0"] + c0, S["p0"] + c1), slice(0, N))
        St = torch.from_numpy(np.ascontiguousarray(Sx)).to(dev)
        St = St * win * torch.from_numpy(amp[c0:c1].astype(np.float32)).to(dev)[:, None]
        A = torch.zeros((c1 - c0, M), dtype=torch.complex64, device=dev)
        A[:, :N - h] = St[:, h:]
        A[:, M - h:] = St[:, :h]
        q = torch.fft.fft(A, dim=1) if S["sgn"] > 0 else torch.fft.ifft(A, dim=1) * M
        del A, St
        for j in range(c1 - c0):
            kk = c0 + j
            d_tx, dead, Lb = el_path(float(tn_tx[kk]), u_tx[kk])
            d_rx, _, _ = el_path(float(tn_rx[kk]), u_rx[kk])
            dR = d_tx + d_rx
            bf = torch.remainder(dR * bins_per_m, float(M))
            i0f = torch.floor(bf)
            fr = bf - i0f
            i0 = torch.remainder(i0f.long(), M)
            i1 = torch.remainder(i0 + 1, M)
            v = q[j][i0] * (1.0 - fr) + q[j][i1] * fr
            v = torch.where(dead, torch.zeros_like(v), v)
            ang = torch.remainder(dR * cyc_per_m, 1.0) * ang_k
            img += v * torch.complex(torch.cos(ang), torch.sin(ang))
            if kk == n // 2:
                Lm = Lb.cpu().numpy().reshape(Sg.shape)
        del q
        if progress:
            progress(c1, n)
        if (c0 // 64) % 40 == 0:
            logger.info(f"    sezione {c1}/{n} impulsi ({time.time() - t0:.0f} s)")
    out = img.cpu().numpy().reshape(Sg.shape)
    del img, Pt
    if dev == "cuda":
        torch.cuda.empty_cache()
    return dict(s=s, z=z, img=out, inside=inside, Lm=Lm, line=np.array(line, dtype=float),
                en=np.array([e0, n0, e1, n1]),
                hs=_eco_height(e0 + (e1 - e0) * s / Ltot, n0 + (n1 - n0) * s / Ltot, half, H))


def plot_eco_section(sec: dict, pn: float, name: str, eps: float, out_dir: str,
                     logger: logging.Logger) -> None:
    """Sezione in dB sul rumore, la stessa normalizzata per profondita' nella pietra (come il
    TGC di un ecografo), parti note tagliate dal piano (+-3 m) e potenza contro profondita'."""
    s, z, inside, Lm, hs = sec["s"], sec["z"], sec["inside"], sec["Lm"], sec["hs"]
    e0, n0, e1, n1 = sec["en"]
    pw = np.abs(sec["img"]) ** 2 / pn
    db = 10 * np.log10(pw + 1e-12)
    ok = inside & np.isfinite(Lm)
    tgc = np.full_like(db, np.nan)
    Lb = np.floor(np.where(ok, Lm, -1)).astype(int)
    for L in np.unique(Lb[Lb >= 0]):
        m = Lb == L
        tgc[m] = db[m] - np.median(db[m])
    st = {"retta_px": [round(float(v), 1) for v in sec["line"]],
          "retta_EN_m": [round(float(v), 1) for v in sec["en"]], "eps": eps,
          "aria_sopra_dB": round(float(np.median(db[~inside & (z[None, :] > hs[:, None] + 2)])), 1),
          "pietra_dB_per_profondita": {}}
    for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 20), (20, 40), (40, 200)):
        m = ok & (Lm >= lo) & (Lm < hi)
        if m.sum() > 10:
            st["pietra_dB_per_profondita"][f"{lo}-{hi} m"] = round(float(10 * np.log10(pw[m].mean())), 1)
    logger.info("  sezione: " + json.dumps(st, ensure_ascii=True))
    with open(os.path.join(out_dir, "sezione.json"), "w") as f:
        json.dump(st, f, indent=1)
    # parti note tagliate dal piano
    Ll = np.hypot(e1 - e0, n1 - n0)
    ud = np.array([e1 - e0, n1 - n0]) / Ll
    kp = []
    for k_ in load_known_structures(logger).get(name, []):
        for pa, pb in np.asarray(k_["seg"], dtype=float).reshape(-1, 2, 3):
            q = pa[None, :] + np.linspace(0, 1, max(2, int(np.linalg.norm(pb - pa) / 0.5)))[:, None] * (pb - pa)
            rel = q[:, :2] - np.array([e0, n0])
            sq, dq = rel @ ud, rel @ np.array([-ud[1], ud[0]])
            m = (np.abs(dq) < 3.0) & (sq >= 0) & (sq <= Ll)
            kp += [(a, b, k_["colore"] or "#4ade80") for a, b in zip(sq[m], q[m, 2])]
    ext = [s[0], s[-1], z[0], z[-1]]
    # larghezza della figura dalle proporzioni della sezione (assi in scala 1:1): niente vuoti
    w_fig = float(np.clip(2.2 + 4.3 * (ext[1] - ext[0]) / (ext[3] - ext[2]), 8, 16))
    fig, ax = plt.subplots(3, 1, figsize=(w_fig, 13), gridspec_kw={"height_ratios": [1, 1, 0.8]})
    im = ax[0].imshow(db.T, origin="lower", extent=ext, cmap="gray", vmin=-5, vmax=35, aspect="equal")
    plt.colorbar(im, ax=ax[0], label="dB sul rumore")
    ax[0].set_title(f"{name}: sezione verticale (SAFT con rifrazione, eps = {eps:g})")
    im = ax[1].imshow(tgc.T, origin="lower", extent=ext, cmap="RdBu_r", vmin=-10, vmax=10, aspect="equal")
    plt.colorbar(im, ax=ax[1], label="dB sulla mediana alla stessa profondita'")
    ax[1].set_title("dentro la pietra, normalizzata per profondita' (TGC)")
    for x in ax[:2]:
        x.plot(s, hs, "y-", lw=1)
        x.axhline(0, color="c", lw=0.6)
        if kp:
            x.scatter([p[0] for p in kp], [p[1] for p in kp], s=2, c=[p[2] for p in kp],
                      label="parti note (entro 3 m dal piano)")
            x.legend(loc="upper right", fontsize=8)
        x.set_xlim(ext[0], ext[1]); x.set_ylim(ext[2], ext[3])
        x.set_xlabel("distanza lungo la retta [m]"); x.set_ylabel("quota z [m]")
    Lv, Pv = Lm[ok], pw[ok]
    edges = np.arange(0, max(2.0, float(Lv.max()) if Lv.size else 2.0) + 1, 1.0)
    idx = np.digitize(Lv, edges) - 1
    mean = [Pv[idx == i].mean() if (idx == i).sum() > 5 else np.nan for i in range(len(edges) - 1)]
    ax[2].plot(edges[:-1] + 0.5, 10 * np.log10(np.array(mean)), "o-", ms=3)
    ax[2].axhline(0, color="k", ls="--", label="rumore")
    ax[2].set_xlabel("cammino nella pietra L [m]"); ax[2].set_ylabel("potenza media [dB sul rumore]")
    ax[2].set_title("decadimento con la profondita' (un'eco interna vera cala di decine di dB/m)")
    ax[2].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "sezione.png"), dpi=100)
    plt.close(fig)


def eco_run(args, logger: logging.Logger, out_root: str) -> None:
    """--eco: test dell'ombra sulla piramide e, se c'e' una retta (--eco_line o --eco_pick),
    la sezione verticale. SLC e rumore sono salvati e riusati se ci sono gia'."""
    nm, nr = args.eco_pyr, float(np.sqrt(args.eco_eps))
    out_dir = os.path.join(out_root, f"eco_{nm}")
    os.makedirs(out_dir, exist_ok=True)
    S = _eco_setup(args.cphd, nm, args.pulse_start, args.n_pulses, logger, args.bp_sign)
    f_slc = os.path.join(out_dir, "slc_eco.npz")
    key = np.array([S["p0"], S["n"], S["sgn"]])         # blocco di impulsi e segno di fase
    r = np.load(f_slc) if os.path.exists(f_slc) else None
    if r is not None and "key" in r and np.array_equal(r["key"], key):
        slc, pn = r["slc"], float(r["pn"])
        logger.info(f"  SLC e rumore riusati da {f_slc} (test dell'ombra gia' fatto)")
    else:
        slc, pn = eco_shadow_test(S, nm, nr, out_dir, logger)
        np.savez(f_slc, slc=slc, pn=pn, key=key)
    plot_orientation(slc, S["geo"], nm, os.path.join(out_dir, "orientamento.png"),
                     os.path.join(out_dir, "orientamento.json"))
    if args.eco_web:
        eco_web(S, nm, nr, pn, args, out_dir, logger)
        return
    if args.eco_pick:
        line = eco_pick_line(slc, S, nm)
    elif args.eco_line:
        line = tuple(float(v) for v in args.eco_line.split(","))
    else:
        logger.info("  nessuna retta (--eco_line, --eco_pick o --eco_web): solo test dell'ombra")
        return
    _eco_do_section(S, nm, line, args, nr, pn, out_dir, logger)


def _eco_do_section(S: dict, nm: str, line, args, nr: float, pn: float, out_dir: str,
                    logger: logging.Logger, progress=None) -> str:
    """Sezione lungo `line`, salvata in <out_dir>/sezione_<punti>/; ritorna il nome della cartella."""
    sec = eco_section(S, nm, line, args.eco_step, args.eco_zlo, nr, logger, progress)
    tag = "sezione_" + "_".join(f"{v:.0f}" for v in line)
    sub = os.path.join(out_dir, tag)
    os.makedirs(sub, exist_ok=True)
    np.savez(os.path.join(sub, "sezione.npz"), **{k: v for k, v in sec.items()})
    plot_eco_section(sec, pn, nm, args.eco_eps, sub, logger)
    logger.info(f"  → {sub}")
    return tag


def eco_web(S: dict, nm: str, nr: float, pn: float, args, out_dir: str,
            logger: logging.Logger) -> None:
    """Pagina locale (solo 127.0.0.1) dell'ecografia: orientamento.png, due clic = la retta,
    "Esegui ecografia" = sezione calcolata qui (GPU, una alla volta) e mostrata nella
    pagina. Il programma resta in ascolto finche' non si preme Ctrl+C."""
    import threading
    import webbrowser
    from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
    tpl = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eco_template.html")
    root = os.path.realpath(out_dir)
    with open(os.path.join(out_dir, "orientamento.json")) as f:
        omap = json.load(f)
    date = re.search(r"_(\d{8})\d{6}_", os.path.basename(args.cphd))
    job = {"busy": False, "done": 0, "tot": 1, "tag": None, "error": None, "lock": threading.Lock()}

    def sections():
        out = []
        for d in sorted(os.listdir(out_dir)):
            fj = os.path.join(out_dir, d, "sezione.json")
            if d.startswith("sezione_") and os.path.exists(fj):
                with open(fj) as f:
                    out.append({"tag": d, **json.load(f)})
        return out

    def worker(line):
        try:
            tag = _eco_do_section(S, nm, line, args, nr, pn, out_dir, logger,
                                  lambda a, b: job.update(done=a, tot=b))
            job.update(tag=tag)
        except Exception as e:                           # l'errore va mostrato nella pagina
            logger.exception(f"  sezione fallita: {e}")
            job.update(error=f"{type(e).__name__}: {e}")
        finally:
            job.update(busy=False)

    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):                       # niente log per ogni richiesta
            pass

        def _send(self, code, body, ctype):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _json(self, obj, code=200):
            self._send(code, json.dumps(obj).encode(), "application/json")

        def do_GET(self):
            path = self.path.split("?")[0]
            if path == "/":
                with open(tpl, encoding="utf-8") as f:
                    data = {"title": f"Ecografia {nm}" + (f" - {date.group(1)}" if date else ""),
                            "map": omap, "eps": args.eco_eps, "step": args.eco_step,
                            "sections": sections()}
                    page = f.read().replace("__DATA__", json.dumps(data))
                self._send(200, page.encode("utf-8"), "text/html; charset=utf-8")
            elif path == "/status":
                self._json({k: job[k] for k in ("busy", "done", "tot", "tag", "error")}
                           | {"sections": sections()})
            elif path.startswith("/files/"):
                f = os.path.realpath(os.path.join(root, path[len("/files/"):]))
                ok = f.startswith(root + os.sep) and os.path.isfile(f) and \
                    f.endswith((".png", ".json"))
                if not ok:
                    return self._send(404, b"non trovato", "text/plain")
                with open(f, "rb") as fh:
                    self._send(200, fh.read(), "image/png" if f.endswith(".png") else "application/json")
            else:
                self._send(404, b"non trovato", "text/plain")

        def do_POST(self):
            if self.path != "/run":
                return self._send(404, b"non trovato", "text/plain")
            try:
                q = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))))
                line = tuple(float(q[k]) for k in ("x0", "y0", "x1", "y1"))
            except (ValueError, KeyError, TypeError):
                return self._json({"error": "retta non valida"}, 400)
            if np.hypot(line[2] - line[0], line[3] - line[1]) < 5:
                return self._json({"error": "retta troppo corta"}, 400)
            with job["lock"]:
                if job["busy"]:
                    return self._json({"error": "c'e' gia' una sezione in calcolo"}, 409)
                job.update(busy=True, done=0, tot=1, tag=None, error=None)
            threading.Thread(target=worker, args=(line,), daemon=True).start()
            self._json({"ok": True})

    srv = None
    for port in range(8770, 8790):
        try:
            srv = ThreadingHTTPServer(("127.0.0.1", port), H)
            break
        except OSError:
            continue
    if srv is None:
        raise RuntimeError("nessuna porta libera tra 8770 e 8789")
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    logger.info(f"  pagina dell'ecografia: {url}  (Ctrl+C per chiudere)")
    if not args.eco_no_browser:
        webbrowser.open(url)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        logger.info("  pagina chiusa")
    finally:
        srv.server_close()


# =============================================================================
#  DEMO
# =============================================================================

def generate_synthetic_data(cfg: Config, logger: logging.Logger) -> np.ndarray:
    logger.info("Generazione dati SAR sintetici (demo)...")
    n_pulses, n_samples = 512, 1024
    raw = (np.random.randn(n_pulses, n_samples) +
           1j * np.random.randn(n_pulses, n_samples)).astype(np.complex64) * 0.1

    k0, x0 = n_samples // 2, n_pulses // 2
    t = np.arange(n_pulses)
    vib = 0.8 * np.exp(1j * 2 * np.pi * 0.05 * t)
    chirp_range = np.sinc((np.arange(n_samples) - k0) / 8)
    chirp_az = np.sinc((t - x0) / 16)
    target = np.outer(chirp_az * vib, chirp_range) * 5.0

    for (kx, xx, amp) in [(k0-200, x0-100, 2.0),
                          (k0+150, x0+80, 1.5),
                          (k0-50, x0+200, 1.2)]:
        cr = np.sinc((np.arange(n_samples) - kx) / 10)
        ca = np.sinc((t - xx) / 20)
        target += amp * np.outer(ca, cr)

    raw += target.astype(np.complex64)
    logger.info(f"  dati sintetici shape: {raw.shape}")
    return raw


# =============================================================================
#  PIPELINE
# =============================================================================

def run_pipeline(cfg: Config, logger: logging.Logger,
                 raw_data: Optional[np.ndarray] = None,
                 slc_in: Optional[np.ndarray] = None) -> dict:
    t0 = time.time()
    os.makedirs(cfg.out_dir, exist_ok=True)

    logger.info(cfg.summary())

    slc = slc_in if slc_in is not None else focus_slc(raw_data, cfg, logger)
    slc_prev = _decimate(slc)
    spectrum = dft2(slc, logger)
    spec_prev = _decimate(spectrum)
    rows, cols = select_tomographic_line(slc, cfg, logger)
    if cfg.save_npy:
        np.save(os.path.join(cfg.out_dir, "slc.npy"), slc)
    del slc
    if cfg.vec_waves:
        res_v = vector_chain_waveforms(spectrum, cfg, logger, rows, cols, cfg.out_dir)
        plot_vector_waveforms(res_v, cfg, os.path.join(cfg.out_dir, "vector_waveforms.png"))
        logger.info(f"  → {os.path.join(cfg.out_dir, 'vector_waveforms.png')}")
    patches, master, slave = generate_doppler_subapartures(
        spectrum, cfg, logger, rows, cols)
    del spectrum
    Y = estimate_micromotion(patches, cfg, logger)
    h = tomographic_focusing(Y, cfg, logger)
    mag, mag_db = generate_tomographic_magnitude(h, cfg, logger)

    logger.info("Step 8: Visualizzazione...")

    out_mag = os.path.join(cfg.out_dir, "tomographic_magnitude.png")
    plot_tomographic_result(mag, cfg, out_mag,
                            title="3D Tomographic Imaging (magnitude)",
                            log_scale=cfg.log_scale)
    logger.info(f"  → {out_mag}")

    out_over = os.path.join(cfg.out_dir, "processing_overview.png")
    plot_processing_overview(slc_prev, spec_prev, mag, cfg, out_over)
    logger.info(f"  → {out_over}")

    out_sub = os.path.join(cfg.out_dir, "doppler_subaparture_strategy.png")
    plot_subaparture_strategy(master, slave, cfg, out_sub)
    logger.info(f"  → {out_sub}")

    out_prof = os.path.join(cfg.out_dir, "profile_validation.png")
    plot_profile_validation(mag, cfg, out_prof)
    logger.info(f"  → {out_prof}")

    if cfg.save_npy:
        np.save(os.path.join(cfg.out_dir, "tomographic_magnitude.npy"), mag)
        np.save(os.path.join(cfg.out_dir, "tomographic_magnitude_db.npy"), mag_db)
        np.save(os.path.join(cfg.out_dir, "h_complex.npy"), h)
        np.save(os.path.join(cfg.out_dir, "line_pixels.npy"),
                np.stack([rows, cols]))
        np.save(os.path.join(cfg.out_dir, "Y_tracking.npy"), Y)
        logger.info(f"  → file .npy salvati in {cfg.out_dir}/")

    elapsed = time.time() - t0
    logger.info(f"Pipeline completata in {elapsed:.2f} s")

    return {
        "slc_prev": slc_prev, "spec_prev": spec_prev,
        "master_prev": master, "slave_prev": slave,
        "Y": Y, "h": h,
        "magnitude": mag, "magnitude_db": mag_db,
        "elapsed_s": elapsed,
    }


# =============================================================================
#  CLI
# =============================================================================

DEFAULT_CPHD = (
    r"e:\claude_code\Biondi_Malanga_1\Capella\capella_data"
    r"\CAPELLA_C13_SP_CPHD_HH_20240916215711_20240916215743.cphd"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="SAR Doppler Tomography — Estrazione immagine 3D "
                    "tomografica (magnitude) da CPHD.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--cphd", type=str, default=DEFAULT_CPHD,
                   help="Percorso file CPHD di input.")
    p.add_argument("--pulse_start", type=int, default=100000,
                   help="Primo impulso del blocco da leggere (centro scena).")
    p.add_argument("--n_pulses", type=int, default=16384,
                   help="Numero di impulsi del blocco (backprojection: letti "
                        "in streaming; 'legacy': caricati in RAM, 2048 ~ 0.4 GiB).")
    p.add_argument("--focus", choices=["bp", "legacy"], default="bp",
                   help="bp = backprojection (corretta per spotlight); "
                        "legacy = vecchia catena FFT approssimata.")
    p.add_argument("--scene_az", type=float, default=600.0,
                   help="Estensione scena in azimuth [m] (backprojection).")
    p.add_argument("--scene_rg", type=float, default=600.0,
                   help="Estensione scena in ground-range [m] (backprojection).")
    p.add_argument("--out", type=str, default="output",
                   help="Directory di output.")
    p.add_argument("--demo", action="store_true",
                   help="Esegue con dati SAR sintetici.")
    p.add_argument("--N_D", type=int, default=8,
                   help="Numero di posizioni di partenza delle sub-aperture Doppler.")
    p.add_argument("--chain", type=int, default=0,
                   help="0 = strategia del paper (B_DL/N_D, N_D coppie). >=1: coppie "
                        "concatenate per posizione, N_D x chain misure.")
    p.add_argument("--tomo_chain", type=int, default=1,
                   help="Come --chain ma per il tomogramma di linea (con N_D grande "
                        "il costo cresce come N_D x chain).")
    p.add_argument("--z_min", type=float, default=-3000.0)
    p.add_argument("--z_MAX", type=float, default=0.0)
    p.add_argument("--dz", type=float, default=0.92)
    p.add_argument("--pol", type=str, default="HH",
                   choices=["HH", "VV"],
                   help="Polarizzazione (default HH per Capella).")
    p.add_argument("--log", action="store_true",
                   help="Usa scala logaritmica (dB) per il tomogramma.")
    p.add_argument("--line_orient", choices=["range", "azimuth"], default="range",
                   help="Linea tomografica: riga a azimuth fisso (range) o "
                        "colonna a range fisso (azimuth).")
    p.add_argument("--line_pos", type=int, default=None,
                   help="Indice di riga/colonna della linea (default: la piu' "
                        "energetica).")
    p.add_argument("--line_len", type=int, default=512,
                   help="Max pixel 'puri' sulla linea.")
    p.add_argument("--B_shift_ratio", type=float, default=1 / 3,
                   help="B_shift = B_DL * ratio (frequenza osservata).")
    p.add_argument("--lpf_sigma", type=float, default=1.0,
                   help="Passa-basso gaussiano su Y(k) dopo il tracking [campioni]; 0 = off.")
    p.add_argument("--shift_px", type=int, default=0,
                   help="Shift master/slave in px (es. 1 con --N_D 64: 64 coppie shiftate di 1 px).")
    p.add_argument("--osc", action="store_true",
                   help="Y(k) dall'oscillatore: pixel = massa 1 spinta dal vettore (a + jb).")
    p.add_argument("--osc_cycles", type=float, default=4.0,
                   help="Cicli dell'oscillatore sulla catena (w0 = 2*pi*cicli/N_D).")
    p.add_argument("--osc_damp", type=float, default=0.1,
                   help="Smorzamento lambda in unita' di w0.")
    p.add_argument("--vec_waves", action="store_true",
                   help="Onde da N_D vettori (Re, Im) per pixel: delta master/slave a catena "
                        "applicati a una massa 1, assi E-O (Re) e N-S (Im), stirate sull'altezza.")
    p.add_argument("--pyr_height", type=float, default=138.5,
                   help="Altezza [m] su cui stirare le onde (Khufu 138.5).")
    p.add_argument("--voxels", action="store_true",
                   help="Pagina HTML 3D con i voxel di tutte le piramidi (cursore: si vedono i "
                        "voxel con ampiezza filtrata sotto il valore). Solo dati reali (bp).")
    p.add_argument("--vox_pyr", choices=["all"] + list(PYRAMIDS), default="all")
    p.add_argument("--vox_cols", type=int, default=0,
                   help="Colonne (Est e Nord) della griglia sulla base di ogni piramide; "
                        "0 = massima risoluzione (una colonna per cella azimuth dell'SLC).")
    p.add_argument("--vox_nd", type=int, default=0,
                   help="Vettori per pixel (shift di 1 px) per le onde dei voxel; 0 = massima "
                        "risoluzione (un tiro per cubo: passo dei tiri = lato del cubo).")
    p.add_argument("--vox_depth", type=float, default=100.0,
                   help="Profondita' [m] sotto il livello zero fino a cui scende l'onda.")
    p.add_argument("--vox_lpf", type=float, default=np.pi / 2,
                   help="Taglio del passa-basso sulle onde [rad per passo]; default pi/2 = "
                        "meta' della Nyquist (la molla fa un ciclo sull'altezza).")
    p.add_argument("--vox_half", type=float, default=100.0,
                   help="Smorzamento in verticale: discesa [m] in cui l'onda di un tiro si "
                        "dimezza (passo verticale = altezza della piramide / numero di delta).")
    p.add_argument("--bp19", action="store_true",
                   help="Biondi Protocol v1.9 (porting GPU) sulle piramidi: banco di coppie "
                        "Doppler R/O, registrazione, correzione con i fianchi, fit a ellisse "
                        "con R^2 corretto; pagina 3D pyramids_bp19.html e spaccato "
                        "pyramids_bp19_spaccato.html, dati in pyramids_bp19_dati*.json.")
    p.add_argument("--bp19_pyr", choices=["all"] + list(PYRAMIDS), default="all")
    p.add_argument("--bp19_os", type=int, default=4,
                   help="Sovracampionamento azimuth della backprojection per --bp19.")
    p.add_argument("--bp19_lines", type=int, default=0,
                   help="Solo le N linee E-O centrali (prova rapida); 0 = tutte.")
    p.add_argument("--bp19_refit", action="store_true",
                   help="--bp19: rifa solo il fit (e la pagina) dalle misure gia' salvate in "
                        "<out>/bp19_<piramide> (bank.npz, geometry.npz, report.json).")
    p.add_argument("--eco", action="store_true",
                   help="Ecografia dal CPHD: test dell'ombra sulla piramide --eco_pyr e, con "
                        "--eco_line o --eco_pick, sezione verticale lungo la retta (SAFT con "
                        "rifrazione nella pietra). Uscita in <out>/eco_<piramide>/.")
    p.add_argument("--eco_pyr", choices=list(PYRAMIDS), default="khufu")
    p.add_argument("--eco_line", type=str, default=None,
                   help="Retta x0,y0,x1,y1 negli assi di orientamento.png (x = range, y = "
                        "azimuth, in pixel).")
    p.add_argument("--eco_pick", action="store_true",
                   help="Apre l'immagine di orientamento.png: due clic = inizio e fine della retta.")
    p.add_argument("--eco_web", action="store_true",
                   help="Pagina locale nel browser: orientamento.png, disegni la retta con due clic e "
                        "lanci l'ecografia (calcolata qui, una sezione alla volta). Ctrl+C per chiudere.")
    p.add_argument("--eco_no_browser", action="store_true",
                   help="Con --eco_web: non aprire il browser (solo l'indirizzo nel log).")
    p.add_argument("--eco_eps", type=float, default=7.0,
                   help="Permittivita' relativa della pietra (indice n = sqrt(eps)).")
    p.add_argument("--eco_step", type=float, default=0.5, help="Passo della sezione [m].")
    p.add_argument("--eco_zlo", type=float, default=-20.0, help="Quota minima della sezione [m].")
    p.add_argument("--tomo3d", action="store_true",
                   help="Volume tomografico (paper: --N_D, --z_min/--z_MAX/--dz) su tutti i "
                        "pixel di --tomo_pyr (base + facce con layover), a blocchi di pixel. "
                        "Solo dati reali (bp).")
    p.add_argument("--zband", type=float, nargs=2, default=[3.1, 30.0], metavar=("LO", "HI"),
                   help="Passa-banda in z del volume: periodi verticali [m], -3.1 dB ai bordi.")
    p.add_argument("--tomo3d_block", type=int, default=40000,
                   help="Pixel per blocco del volume tomografico.")
    p.add_argument("--window", choices=["rect", "hann"], default="rect")
    p.add_argument("--center_x", type=float, default=0.0,
                   help="Centro immagine: offset in ground-range [m] dall'SRP.")
    p.add_argument("--center_y", type=float, default=0.0,
                   help="Centro immagine: offset in azimuth [m] dall'SRP.")
    p.add_argument("--bp_sign", type=int, choices=[-1, 1], default=-1,
                   help="Moltiplicatore di SGN per la backprojection. -1 (default): "
                        "verificato con la posizione nota di Chefren, che con +1 "
                        "compare riflessa attorno all'SRP.")
    p.add_argument("--steer_2pi", action="store_true",
                   help="Matrice di steering letterale exp(j*2*pi*Kz*z) (default: "
                        "exp(j*Kz*z), Kz in rad/m).")
    p.add_argument("--f_inv", type=float, default=None,
                   help="Frequenza di investigazione [Hz] (paper: 12500 nel testo, "
                        "25000 nella formula).")
    p.add_argument("--A_orbit", type=float, default=None,
                   help="Apertura orbitale A [m] (paper: 42000 nel testo, 84000 nella "
                        "formula).")
    p.add_argument("--tomo_line", nargs="+", choices=["W", "Nup"], default=None,
                   help="Tomogramma lungo una linea Nord->Sud su una faccia "
                        "illuminata: W = orizzontale sulla faccia ovest; Nup = "
                        "risale la faccia nord.")
    p.add_argument("--tomo_pyr", choices=list(PYRAMIDS), default="khufu")
    p.add_argument("--tomo_h", type=float, default=0.3,
                   help="Quota relativa (0-1) della linea W sulla faccia.")
    p.add_argument("--tomo_nd", type=int, default=192,
                   help="Passi Doppler N_D per il tomogramma di linea: servono "
                        "molti passi perche' il periodo di ambiguita' in z "
                        "(~N_D * risoluzione) copra l'altezza della piramide.")
    p.add_argument("--tomo_z0", type=float, default=-10.0)
    p.add_argument("--tomo_z1", type=float, default=150.0)
    p.add_argument("--tomo_dz", type=float, default=0.4,
                   help="Passo in z [m]; deve stare sotto pi/Kz_max (~0.4 m) o i "
                        "picchi tra due campioni vanno persi.")
    p.add_argument("--null", action="store_true",
                   help="Esegue il null test (ROI vs finestre fuori ROI).")
    p.add_argument("--roi", type=int, nargs=4, metavar=("AZ0", "AZ1", "RG0", "RG1"),
                   default=None,
                   help="ROI della piramide in pixel dell'SLC (per --null).")
    p.add_argument("--n_null", type=int, default=300,
                   help="Numero di finestre null.")
    p.add_argument("--faces", choices=["khufu", "khafre", "both"], default=None,
                   help="Null test sulle 4 facce (N,E,S,O) della/e piramide/i: i "
                        "poligoni sono calcolati dalla geometria (coordinate note, "
                        "layover dell'apice). Richiede che la scena le contenga "
                        "(usa --center_x/--center_y/--scene_*).")
    p.add_argument("--face_mode", choices=["diag", "layover"], default="diag",
                   help="Definizione delle facce: diag = 4 triangoli delimitati dalle "
                        "diagonali della base (default); layover = triangoli base-apice "
                        "con l'apice spostato dal layover.")
    p.add_argument("--null_inject", type=float, default=0.0,
                   help="Controllo positivo: px di vibrazione coerente iniettati "
                        "nelle regioni testate.")
    p.add_argument("--cpu", action="store_true",
                   help="Disabilita la GPU (default: usa CUDA se disponibile).")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def main() -> int:
    global USE_GPU
    args = parse_args()
    if args.cpu:
        USE_GPU = False

    cfg = Config(
        N_D=args.N_D,
        z_min=args.z_min,
        z_MAX=args.z_MAX,
        dz=args.dz,
        polarization=args.pol,
        out_dir=args.out,
        log_scale=args.log,
        line_orient=args.line_orient,
        line_pos=args.line_pos,
        line_len=args.line_len,
        B_shift_ratio=args.B_shift_ratio,
        window=args.window,
        steer_2pi=args.steer_2pi,
        chain_len=args.chain,
        lpf_sigma=args.lpf_sigma,
        shift_px=args.shift_px,
        osc=args.osc,
        osc_cycles=args.osc_cycles,
        osc_damp=args.osc_damp,
        vec_waves=args.vec_waves,
        pyr_height=args.pyr_height,
    )
    if args.f_inv is not None or args.A_orbit is not None:
        if args.f_inv is not None:
            cfg.f_investigation = args.f_inv
        if args.A_orbit is not None:
            cfg.A_orbit = args.A_orbit
        cfg.__post_init__()
    logger = setup_logger(cfg.out_dir, verbose=args.verbose)

    try:
        if args.demo:
            raw = generate_synthetic_data(cfg, logger)
        else:
            # bp: per i metadati basta un blocco minuscolo (il segnale viene
            # poi letto in streaming); legacy: carica il blocco intero in RAM
            raw, meta = read_cphd(args.cphd, logger, args.pulse_start,
                                  64 if args.focus == "bp" else args.n_pulses)

            # Override parametri dai metadati se disponibili
            if meta.get("fc"):         cfg.fc = float(meta["fc"])
            if meta.get("bw"):         cfg.bw = float(meta["bw"])
            if meta.get("prf"):        cfg.prf = float(meta["prf"])
            if meta.get("B_cD"):       cfg.B_cD = float(meta["B_cD"])
            if meta.get("v_platform"): cfg.v_platform = float(meta["v_platform"])
            if meta.get("r_ref"):      cfg.r_ref = float(meta["r_ref"])
            if meta.get("polarization") not in (None, "", "UNK"):
                cfg.polarization = str(meta["polarization"])
            if meta.get("domain_type"):
                cfg.domain_type = str(meta["domain_type"])

            # Ricalcola parametri derivati
            cfg.__post_init__()

        if args.eco:
            if args.demo:
                raise SystemExit("--eco richiede dati reali (file CPHD)")
            eco_run(args, logger, cfg.out_dir)
            return 0

        if args.voxels:
            if args.demo or args.focus != "bp":
                raise SystemExit("--voxels richiede dati reali con --focus bp")
            cfg_v = dataclasses.replace(cfg, N_D=args.vox_nd, z_bottom=-args.vox_depth,
                                        wave_lpf=args.vox_lpf, wave_half=args.vox_half)
            names = list(PYRAMIDS) if args.vox_pyr == "all" else [args.vox_pyr]
            pyrs = []
            for nm in names:
                side, H = PYRAMIDS[nm][2], PYRAMIDS[nm][3]
                slc, geo = backproject_cphd(args.cphd, cfg_v, logger, args.pulse_start,
                                            args.n_pulses, side * 1.6, side * 1.6 + H * 1.2,
                                            sign_mult=args.bp_sign, center_on=nm)
                pyrs.append(voxel_pyramid(slc, geo, cfg_v, logger, nm, args.vox_cols,
                                          os.path.join(cfg.out_dir, nm)))
                del slc
            date = re.search(r"_(\d{8})\d{6}_", os.path.basename(args.cphd))
            out_html = os.path.join(cfg.out_dir, "pyramids_3d.html")
            write_voxel_html(pyrs, "Voxel 3D delle piramidi" +
                             (f" - {date.group(1)}" if date else ""), out_html)
            logger.info(f"  → {out_html}")
            return 0

        if args.bp19:
            if args.demo or args.focus != "bp":
                raise SystemExit("--bp19 richiede dati reali con --focus bp")
            names = list(PYRAMIDS) if args.bp19_pyr == "all" else [args.bp19_pyr]
            pyrs = [bp19_pyramid(args, cfg, logger, nm, os.path.join(cfg.out_dir, f"bp19_{nm}"))
                    for nm in names]
            date = re.search(r"_(\d{8})\d{6}_", os.path.basename(args.cphd))
            title = "Biondi Protocol v1.9 sulle piramidi" + (f" - {date.group(1)}" if date else "")
            for out in write_bp19_pages(pyrs, title, cfg.out_dir):
                logger.info(f"  → {out} ({os.path.getsize(out) / 1e6:.1f} MB)")
            return 0

        if args.tomo3d:
            if args.demo or args.focus != "bp":
                raise SystemExit("--tomo3d richiede dati reali con --focus bp")
            nm = args.tomo_pyr
            side, H = PYRAMIDS[nm][2], PYRAMIDS[nm][3]
            slc, geo = backproject_cphd(args.cphd, cfg, logger, args.pulse_start,
                                        args.n_pulses, side * 1.6, side * 1.6 + H * 1.2,
                                        sign_mult=args.bp_sign, center_on=nm)
            vol = tomo_volume(slc, geo, cfg, logger, nm, tuple(args.zband), args.tomo3d_block)
            out_png = os.path.join(cfg.out_dir, f"tomo3d_{nm}.png")
            plot_tomo_volume(vol, slc, cfg, out_png)
            np.savez(os.path.join(cfg.out_dir, f"tomo3d_{nm}.npz"), rows=vol["rows"],
                     cols=vol["cols"], z=vol["z"], H=vol["H"], Hbp=vol["Hbp"], Q=vol["Q"],
                     slc_shape=np.array(slc.shape), da=geo["da"], dr=geo["dr"])
            logger.info(f"  → {out_png}")
            return 0

        if not args.demo and args.focus == "bp":
            slc, geo = backproject_cphd(args.cphd, cfg, logger, args.pulse_start,
                                   args.n_pulses, args.scene_az, args.scene_rg,
                                   center=(args.center_x, args.center_y),
                                   sign_mult=args.bp_sign)
            if args.null:
                if args.roi is None:
                    raise SystemExit("--null richiede --roi AZ0 AZ1 RG0 RG1")
                null_test(slc, cfg, logger, tuple(args.roi), args.n_null,
                          inject_px=args.null_inject, aspect=geo["da"] / geo["dr"])
            if args.tomo_line:
                for ln in args.tomo_line:
                    line_tomogram(slc, cfg, logger, geo, args.tomo_pyr, ln,
                                  args.tomo_h, args.tomo_nd,
                                  (args.tomo_z0, args.tomo_z1),
                                  tomo_dz=args.tomo_dz, tomo_chain=args.tomo_chain)
            if args.faces:
                pyr = ["khufu", "khafre"] if args.faces == "both" else [args.faces]
                faces_test(slc, cfg, logger, geo, pyr, args.n_null,
                           inject_px=args.null_inject, mode=args.face_mode)
            run_pipeline(cfg, logger, slc_in=slc)
        else:
            run_pipeline(cfg, logger, raw_data=raw)
        return 0

    except Exception as e:
        logger.exception(f"Errore durante l'esecuzione: {e}")
        return 2


if __name__ == "__main__":
    sys.exit(main())

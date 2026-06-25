#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sar_tomo.py — micro-Doppler SAR tomography library (Biondi & Malanga, RS 2022, 14, 5231).

Self-contained core used by the `sar-doppler-tomography` skill:
  - GeometriaSAR              : acquisition geometry + tomographic resolution / Kz
  - steering_matrix           : A(Kz,z) = exp(j Kz z)            (Eq. 22)
  - sotto_aperture_doppler    : split an SLC azimuth line into Doppler sub-apertures
  - pixel_tracking            : sub-pixel master/slave coherence+phase  -> Y sample
  - tomografia_beamforming    : h(z) = Aᴴ Y  (Eq. 24); tomografia_capon = MVDR
  - genera_Y_sintetico        : forward model to VERIFY the chain recovers known depths
  - leggi_geometria           : read geometry dict from a Sentinel-1 annotation.xml

Dependencies: numpy, scipy (Capon/freq optional). matplotlib only in the image script.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass

C = 299_792_458.0


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
@dataclass
class GeometriaSAR:
    f_em: float = 5.405e9          # carrier [Hz] (Sentinel-1 C-band)
    c: float = C
    R0: float = 850_000.0          # reference slant range [m]
    theta_deg: float = 39.5        # incidence angle [deg]
    V: float = 7_600.0             # platform speed [m/s]
    B_doppler: float = 1_500.0     # synthesized Doppler bandwidth [Hz]
    apertura_orbitale: float = 84_000.0   # A: orbital aperture for the sonic synthesis [m]
    # sonic/seismic wave used by the tomography (NOT the EM wave)
    v_sonic: float = 6_000.0
    f_sonic: float = 12_500.0

    @property
    def lambda_em(self) -> float:
        return self.c / self.f_em

    @property
    def lambda_sonic(self) -> float:
        # v / (2 f): round trip, ~0.24 m with v=6000, f=12500
        return self.v_sonic / (2.0 * self.f_sonic)

    @property
    def theta(self) -> float:
        return np.deg2rad(self.theta_deg)

    def risoluzione_tomografica(self) -> float:
        """delta_z = lambda_sonic * R / (2 A)."""
        return self.lambda_sonic * self.R0 / (2.0 * self.apertura_orbitale)


def baseline_ortogonali(k: int, B_perp_max: float = 42_000.0) -> np.ndarray:
    """k orthogonal baselines B_perp, symmetric over the orbital aperture (Fig. 8d)."""
    return np.linspace(-B_perp_max, B_perp_max, k)


def wavenumber_verticale(geo: GeometriaSAR, B_perp: np.ndarray,
                         r_i: np.ndarray | None = None) -> np.ndarray:
    """Kz = 4 pi B_perp / (lambda_sonic r_i sin theta)   (Eq. 22)."""
    if r_i is None:
        r_i = np.full_like(B_perp, geo.R0)
    return 4.0 * np.pi * B_perp / (geo.lambda_sonic * r_i * np.sin(geo.theta))


def steering_matrix(geo: GeometriaSAR, z: np.ndarray, k: int) -> np.ndarray:
    """A(Kz,z) in C^{k x F}: A[i,f] = exp(j Kz_i z_f)."""
    Kz = wavenumber_verticale(geo, baseline_ortogonali(k))
    return np.exp(1j * np.outer(Kz, z))


def altezza_ambiguita(geo: GeometriaSAR, k: int) -> float:
    Kz = wavenumber_verticale(geo, baseline_ortogonali(k))
    return float(2 * np.pi / abs(Kz[1] - Kz[0]))


# ---------------------------------------------------------------------------
# Inversion
# ---------------------------------------------------------------------------
def tomografia_beamforming(A: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """h(z) = Aᴴ Y  (matched filter / DFT, Eq. 24)."""
    return A.conj().T @ Y


def tomografia_pinv(A: np.ndarray, Y: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(A) @ Y


def tomografia_capon(A: np.ndarray, Y: np.ndarray, look: np.ndarray | None = None,
                     diag_load: float = 1e-3) -> np.ndarray:
    """Capon / MVDR: P(z) = 1 / (aᴴ R⁻¹ a). Pass snapshots in `look` (k x L) for a
    full-rank covariance; otherwise R = Y Yᴴ + diagonal loading."""
    y = Y.reshape(-1, 1)
    R = (look @ look.conj().T / look.shape[1]) if look is not None else (y @ y.conj().T)
    R += diag_load * np.trace(R) / R.shape[0] * np.eye(R.shape[0])
    Rinv = np.linalg.inv(R)
    out = np.empty(A.shape[1])
    for f in range(A.shape[1]):
        a = A[:, f:f+1]
        out[f] = 1.0 / np.real(a.conj().T @ Rinv @ a).item()
    return out.astype(complex)


# ---------------------------------------------------------------------------
# SLC -> data vector Y
# ---------------------------------------------------------------------------
def sotto_aperture_doppler(slc_line: np.ndarray, n_sub: int,
                           frazione_band: float = 0.5) -> list[np.ndarray]:
    """Split an SLC azimuth line into n_sub Doppler sub-apertures: FFT in azimuth,
    n_sub equispaced band-pass over the central frazione_band of B_CD, IFFT back."""
    npix = slc_line.shape[-1]
    Spec = np.fft.fftshift(np.fft.fft(slc_line, axis=-1), axes=-1)
    band = int(npix * frazione_band)
    start0 = (npix - band) // 2
    win = max(1, band // n_sub)
    sub = []
    for i in range(n_sub):
        mask = np.zeros(npix)
        a = start0 + i * win
        b = min(a + win, npix)
        mask[a:b] = np.hanning(b - a) if (b - a) > 1 else 1.0
        si = np.fft.ifft(np.fft.ifftshift(Spec * mask, axes=-1), axis=-1)
        sub.append(si)
    return sub


def pixel_tracking(master: np.ndarray, slave: np.ndarray, pixel_idx: int,
                   search: int = 4, ovs: int = 16) -> complex:
    """One tomographic observation (Eq. 21) for a pixel: sub-pixel shift magnitude
    weights the coherence, the inter-look phase carries Kz·z.
    Returns (1+|shift|)·coh·exp(j·phase) — never zero, preserving displacement
    sensitivity even when co-registered sub-apertures have integer shift 0."""
    i0 = max(0, pixel_idx - search)
    i1 = min(master.shape[-1], pixel_idx + search + 1)
    m, s = master[..., i0:i1], slave[..., i0:i1]
    M = np.fft.fft(m, n=(i1 - i0) * ovs, axis=-1)
    S = np.fft.fft(s, n=(i1 - i0) * ovs, axis=-1)
    xc = np.abs(np.fft.ifft(M * np.conj(S), axis=-1)).ravel()
    peak = int(np.argmax(xc))
    if peak > xc.size // 2:
        peak -= xc.size
    shift = peak / ovs
    g = np.vdot(s.ravel(), m.ravel())                       # phase = phi_m - phi_s
    coh = np.abs(g) / (np.linalg.norm(m) * np.linalg.norm(s) + 1e-12)
    return (1.0 + abs(shift)) * coh * np.exp(1j * np.angle(g))


def costruisci_Y(slc_line: np.ndarray, pixel_idx: int, k: int) -> np.ndarray:
    """k-look data vector Y for a pixel: tracking between consecutive sub-aperture pairs."""
    sub = sotto_aperture_doppler(slc_line, n_sub=k + 1)
    return np.array([pixel_tracking(sub[i], sub[i + 1], pixel_idx) for i in range(k)],
                    dtype=complex)


# ---------------------------------------------------------------------------
# Synthetic forward model — verify the chain recovers known depths
# ---------------------------------------------------------------------------
def genera_Y_sintetico(geo: GeometriaSAR, k: int,
                       profondita=(12.0, 41.0, 78.0), ampiezze=(1.0, 0.7, 0.9),
                       snr_db: float = 20.0, n_snapshot: int = 1, seed: int = 0):
    """Y = A(z_reflectors) @ amps + noise. Inversion must recover the injected depths."""
    rng = np.random.default_rng(seed)
    A_r = steering_matrix(geo, np.asarray(profondita, float), k)
    amp = np.asarray(ampiezze, float)
    cols = []
    for _ in range(n_snapshot):
        Y = A_r @ (amp * np.exp(1j * rng.uniform(0, 2 * np.pi, amp.size)))
        p = np.mean(np.abs(Y) ** 2) / (10 ** (snr_db / 10))
        cols.append(Y + np.sqrt(p / 2) * (rng.standard_normal(k) + 1j * rng.standard_normal(k)))
    M = np.stack(cols, 1)
    return M[:, 0] if n_snapshot == 1 else M


# ---------------------------------------------------------------------------
# Annotation reader
# ---------------------------------------------------------------------------
def leggi_geometria(ann_path: str) -> dict:
    """Geometry dict from a Sentinel-1 annotation.xml."""
    import xml.etree.ElementTree as ET
    r = ET.parse(ann_path).getroot()
    f = lambda p: float(r.find(p).text)
    g = {"f_em": f(".//radarFrequency"), "dt_az": f(".//azimuthTimeInterval"),
         "slantRangeTime": f(".//imageInformation/slantRangeTime"),
         "incid_mid": f(".//incidenceAngleMidSwath"),
         "dR": f(".//rangePixelSpacing"), "dA": f(".//azimuthPixelSpacing")}
    g["R_near"] = C * g["slantRangeTime"] / 2.0
    vs = []
    for o in r.findall(".//orbit"):
        v = o.find("velocity")
        if v is not None:
            vs.append(sum(float(v.find(c).text) ** 2 for c in "xyz") ** 0.5)
    g["V"] = float(np.mean(vs)) if vs else 7600.0
    return g


if __name__ == "__main__":
    # self-test: recover injected reflectors at 12 / 41 / 78 m
    geo = GeometriaSAR()
    k = 48
    Y = genera_Y_sintetico(geo, k, snr_db=20.0)
    z = np.linspace(0, 100, 400)
    h = np.abs(tomografia_beamforming(steering_matrix(geo, z, k), Y))
    from scipy.signal import find_peaks
    pk, _ = find_peaks(h / h.max(), height=0.3, distance=8)
    print("z_amb =", round(altezza_ambiguita(geo, k), 1), "m")
    print("recovered peaks at z =", np.round(z[pk], 1), "m  (injected 12/41/78)")

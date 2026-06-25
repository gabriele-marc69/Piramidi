---
name: sar-doppler-tomography
description: Reconstruct the vertical (depth) reflectivity profile of a target from a single or multi-temporal Sentinel-1 SLC stack using the Biondi & Malanga micro-motion Doppler-tomography method (Remote Sens. 2022, 14, 5231; arXiv:2206.09200). Use to build a tomographic reflectivity volume V[azimuth,range,z], B-scan depth sections, horizontal slices, and a tomographic-height map — e.g. to reproduce the Great Pyramid / Khafre internal-structure figures from a real SLC chip stack.
---

# SAR micro-Doppler tomography (Biondi–Malanga)

**Principle.** EM radar does not penetrate solids — the SLC image is surface-only. But the
structure *vibrates* (seismic/anthropic/wind micro-motion, MM). MM imprints a *micro-Doppler*
on the SLC phase. Estimating that vibration per pixel across Doppler sub-apertures and
refocusing it along the vertical wavenumber yields a *sonic* tomography of the interior —
the SAR is used as a receiver of seismic/acoustic waves, not as a penetrating illuminator.

## Computational chain (paper Fig. 9)
```
SLC chip  ──FFT azimuth──▶ N Doppler sub-apertures (master/slave pairs)
          ──pixel tracking──▶ per-pixel observation vector Y (k looks, complex)
          ──A(Kz,z)──▶ h(z) = Aᴴ·Y  (beamforming) or Capon/MVDR  ──▶ |h(z)|
stack the per-pixel profiles ──▶ reflectivity volume V[az, rng, z]
```

Key equations (in `scripts/sar_tomo.py`):
- Vertical wavenumber **Kz = 4π·B⊥ / (λ_sonic · r · sinθ)** (Eq. 22).
- Steering matrix **A[i,f] = exp(j·Kz_i·z_f)** (k looks × F depths).
- Inversion **h(z) = Aᴴ·Y** (matched filter / DFT, Eq. 24); `tomografia_capon` for MVDR.
- Tomographic resolution **δz = λ_sonic·R / (2A)**; ambiguity height **z_amb = 2π/|ΔKz|**.

## Workflow
1. Get the complex chip stack `(n_scene, H, W)` and geometry — use the
   `sentinel1-slc-reader` skill.
2. Build `GeometriaSAR` from the annotation (f_em, V, θ, R0 = R_near + col·dR).
3. Choose `K_LOOK` (number of Doppler sub-apertures = tomographic looks), `Z_MAX`
   (stay within `z_amb`!), and depth oversampling.
4. For each pixel/scene: `sotto_aperture_doppler` → `pixel_tracking` pairs → `Y`,
   then accumulate `|Aᴴ·Y|` into `V[az,rng,:]`. Average over scenes.
5. Render the **publication-style images** with `scripts/tomographic_images.py`:
   - `sezione_range_quota.png` / `sezione_azimuth_quota.png` — vertical B-scans
     (reflectivity vs depth), the analogue of the paper's tomographic-line sections;
   - `slice_orizzontali.png` — horizontal slices of V at increasing depth;
   - `mappa_quota.png` — dominant-scatterer height (argmax|h|) + peak-intensity map;
   - `volume_tomografico.npy` — the full V for 3-D scatterer clouds.

```bash
python scripts/tomographic_images.py --stack /path/to/stack_slc --outdir tomo_img \
  --zmax 8.5 --klook 12 --ovs 4
```

## Critical limit — read before interpreting results
The ambiguity height with few azimuth samples is small. For the Giza IW2 box (~26 azimuth
samples, 6 scenes) `z_amb ≈ 8.5 m`: depths beyond that **alias**. So:
- the B-scans are valid vertical reflectograms **only within 0…z_amb**;
- absolute pyramid-scale height (~140 m) is **not** recoverable from this stack — values
  read as height *modulo* `z_amb`. Reaching pyramid scale needs hundreds of looks
  (far more acquisitions / baseline diversity) and ideally a Capon/MVDR estimator.

Be honest about this in any write-up: the method is mathematically validated (it recovers
reflectors injected at known depths in the synthetic forward model `genera_Y_sintetico`),
but a single open Sentinel-1 stack reproduces the *image types* of the paper, not its
deep, high-resolution interior — that used many calibrated narrow Doppler sub-apertures.

## Files
- `scripts/sar_tomo.py` — library: geometry, steering matrix, sub-apertures, pixel
  tracking, beamforming/Capon, synthetic forward model for verification.
- `scripts/tomographic_images.py` — end-to-end volume + the four publication figures.

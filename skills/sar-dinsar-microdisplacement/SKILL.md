---
name: sar-dinsar-microdisplacement
description: Estimate per-pixel line-of-sight micro-displacements (mm) of a target from the interferometric phase of a multi-temporal Sentinel-1 SLC chip stack, via differential interferometry (DInSAR). Use to produce the micro-motion waveforms, displacement maps, and 3-D displacement clouds that underlie the Biondi & Malanga InSAR figures (Figs. 10–17), or to assess temporal coherence/stability of a structure before tomography.
---

# DInSAR micro-displacement from an SLC chip stack

The micro-motion (MM) that feeds the Doppler tomography also shows up directly as
interferometric phase between repeat passes. This skill turns a complex chip stack
`(n_scene, az, rng)` (from the `sentinel1-slc-reader` skill) into LOS displacements.

## Method
- Reference = earliest scene `z_ref`.
- Per scene `s`, per pixel: interferometric phase `phi = angle(z_s · conj(z_ref))`.
  The conjugate product cancels the common scattering phase; subtracting raw phases
  would **not**.
- LOS displacement `d = -(λ / 4π) · phi`  [mm], with `λ = c/f_em ≈ 55.5 mm` for
  Sentinel-1 C-band → `λ/4π ≈ 4.41 mm/rad`, one fringe `2π = 27.7 mm`.
- Temporal coherence per pixel `γ = |Σ z_s conj(z_ref)| / (‖z_s‖‖z_ref‖)`; linear
  regression of `d` vs time → velocity [mm/yr] and a modeled waveform per pixel.

## Workflow
1. Build/Load the chip stack `.npz` (`vv`, `dR`, `dA`, `incid_mid`, `vv_scene`).
2. Run `scripts/dinsar.py --npz box.npz` → per-pixel/date CSV, per-pixel velocity CSV,
   the waveform figure and a 3-D displacement cloud.
3. **Prefer a single-satellite sub-stack** (e.g. only S1D scenes) for cleaner phase:
   `--sat s1d`. Mixed S1A/S1C/S1D have uncompensated baselines → topo + atmosphere leak
   into the phase and coherence drops.

## Honest limits (state them in any write-up)
- Phase is **wrapped** (no unwrapping): values are relative, ambiguous modulo `λ/2 = 27.7 mm`.
- Mixed-sensor stack: coherence ~0.36 (noise-dominated). S1D-only: ~0.47, sub-cm stability
  over ~24 days on the coherent pixels — illustrative, not geodetic-grade.
- These are **micro**-displacements of the surface scatterers, not absolute deformation.

## File
- `scripts/dinsar.py` — DInSAR estimation + figures from the chip-stack `.npz`.

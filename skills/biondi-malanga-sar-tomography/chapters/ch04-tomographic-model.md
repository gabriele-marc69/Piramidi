# Section 4: Tomographic Model

## Core Idea
Each pixel on a chosen "tomographic line" is modeled as a damped, driven harmonic oscillator (mass on a spring) whose vibration is excited by ambient seismic background, nearby human activity (Cairo), and wind; inverting the multi-baseline vibration measurements along that line (analogous to classical multi-baseline SAR tomography, but using vibration/phase instead of direct backscatter height) produces a depth-resolved tomogram.

## Frameworks Introduced
- **Harmonic-oscillator vibrational model of the Earth/structure**: every pixel along the tomographic line is treated as a mass-spring system; nonlinear restoring force is derived and reduced (under a low-nonlinearity assumption) to a two-degree-of-freedom linear harmonic oscillator `r(t) = (a·cos(ω₀t), b·sin(ω₀t))·exp(−γt/2)`, whose amplitude/phase parameters `{a, b}` come directly from the sub-aperture coregistration of Section 3.
  - When to use: whenever you want to convert a *single* SAR image's Doppler-domain vibration measurements into a physically-motivated, invertible tomographic model.
  - How: (1) build the multi-frequency data vector `Y` from `k` sub-aperture time samples; (2) build a steering matrix `A(Kz, z)` encoding phase vs. elevation `z` for each orthogonal baseline; (3) invert `Y = A(Kz,z)·h(z)` via matched-filter approximation `h(z) ≈ Aᴴ(Kz,z)·Y` to get the tomographic reflectivity profile.
- **Tomographic resolution formula**: `δz = λ·R / (2·A)`, where `λ` is the sound (vibration) wavelength, `R` is slant range, and `A` is the orbital aperture used for tomographic synthesis — directly analogous to classical SAR-tomography vertical resolution formulas, but applied to acoustic/vibrational wavelength instead of electromagnetic wavelength.

## Key Concepts
- **Tomographic line**: a chosen line of contiguous pixels (from surface to depth) along which the vibration-based inversion is performed; its orientation (vertical/horizontal, N/S/E/W) determines what part of the pyramid's interior becomes visible.
- **Orthogonal baseline (`B_i`)**: per-sub-aperture geometric baseline used to build the steering matrix, analogous to multi-pass baseline diversity in classical SAR tomography — but here synthesized from Doppler sub-apertures of a *single* pass rather than from multiple satellite passes.
- **Computational block chain**: SLC image → range/Doppler bandpass filters (per Section 3 sub-aperture strategy) → spatial filter isolating the tomographic line → vibration-metric algorithm → focusing of estimated vibrations into a depth tomogram (paper's Figure 9 flowchart, described narratively here).

## Reference Tables
Worked numeric example from the paper (Section 4 resolution calculation):

| Quantity | Value | Meaning |
|---|---|---|
| Seismic wave propagation speed `v` | ≈ 6000 m/s | Assumed medium propagation speed |
| Investigation frequency `f` | 12,500 Hz | Chosen vibration frequency of interest |
| Derived wavelength `λ = v/f` | ≈ 0.24 m | Vibration wavelength used in resolution calc |
| Slant range `R` | 650,000 m | Satellite-to-target distance |
| Orbital aperture `A` (half total synthesized orbit) | ≈ 84,000 m | Aperture used for tomographic synthesis |
| Resulting resolution `δz = λR/2A` | ≈ 0.92 m | Vertical/depth tomographic resolution used for all Section 5 results |

## Key Takeaways
1. The tomographic depth resolution actually used throughout the paper's results (Section 5) is ≈0.92 m per pixel, derived from a specific assumed vibration wavelength and orbital aperture — this number is the practical bound on how finely internal features can be distinguished.
2. The oscillator model is deliberately reduced from a nonlinear restoring-force expansion to a linear 2-DOF damped oscillator to keep the inversion tractable.
3. The tomographic inversion `h(z) ≈ Aᴴ·Y` is a matched-filter (Fourier-like) approximation, not an exact least-squares or compressive-sensing solve — a methodological simplification worth flagging when evaluating result reliability.

## Connects To
- **Section 3**: supplies the per-pixel time-domain vibration vector `Y` that this section's steering matrix inverts.
- **Section 5.19 (Metric Determination)**: applies the ≈0.92 m resolution derived here to report physical dimensions of discovered structures.

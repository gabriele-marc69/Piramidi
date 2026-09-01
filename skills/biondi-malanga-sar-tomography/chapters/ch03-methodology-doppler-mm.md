# Section 3: Methodology — Doppler Sub-Aperture MM Estimation

## Core Idea
A single SAR single-look-complex (SLC) image is split into multiple partially/non-overlapped Doppler sub-apertures; tracking the sub-pixel displacement of a target across these sub-apertures (via MCA in the Doppler direction) yields a time-domain vibration signal for that pixel, which becomes the raw input to tomography.

## Frameworks Introduced
- **Doppler Sub-Aperture Decomposition**: split the full synthesized Doppler bandwidth `B_cD` into `N_D` narrower sub-bands (leaving out a guard band `B_DL = B_cD/2` around the matched-filter boundary for sensitivity), each acting as an independent lower-resolution "look" of the same scene at a different point in the SAR observation interval.
  - When to use: whenever you need a time series of a target's position/phase from a *single* SAR acquisition (no repeat-pass needed) to estimate motion.
  - How: (1) focus the SLC image; (2) partition Doppler bandwidth into master/slave sub-aperture matrices; (3) sub-pixel coregister each sub-aperture pair; (4) track pixel displacement across sub-apertures over time.
- **Target motion model (range/azimuth anomalies)**: a moving/vibrating point target produces three characteristic SAR artifacts — azimuth displacement (constant range velocity), azimuth smearing (azimuth velocity/range acceleration), and range-walk/defocusing (range speed). These are formalized via a Taylor-expanded range-history equation relating target velocity/acceleration components `{v_rs, v_a, a_r, a_a}` to focusing-domain parameters `{r1 (range velocity), r2 (range acceleration), c1 (azimuth velocity)}`.

## Key Concepts
- **SLC (Single-Look-Complex)**: the focused, complex-valued SAR image format used as the sole raw input (only one acquisition needed — no interferometric pair required for the MM step).
- **Chirp bandwidth vs. Doppler bandwidth**: range-direction (`B_cr`) vs. azimuth-direction (`B_cD`) bandwidths of the focused SAR signal; both appear as sinc functions in the 2D spectrum.
- **Pixel/sub-pixel tracking / coregistration**: the core measurement operation — estimating the displacement of a pixel of interest between "master" and "slave" sub-aperture realizations.
- **Guard band (`B_DL`)**: bandwidth deliberately excluded from the matched filter to preserve sensitivity to target motion.

## Reference Tables
Sub-aperture strategy quantities (paraphrased):

| Symbol | Meaning |
|---|---|
| `B_cr`, `B_cD` | Total chirp (range) / Doppler (azimuth) bandwidth |
| `B_DL = B_cD/2` | Bandwidth withheld from matched filter for motion sensitivity |
| `N_D` | Number of Doppler sub-apertures the guard band is split into |
| `N_c` | Number of rigid master–slave shifts along azimuth bandwidth |
| `L_sa = 2Nd` | Total synthetic aperture length |

## Worked Example
The paper demonstrates the method on one pixel (circled in yellow, "pixel 1") in the SAR image of the pyramid: the raw time-domain magnitude displacement trend is plotted in range and azimuth separately, each shown as (a) an unfiltered blue trace and (b) a smoothed positive-envelope red trace. This single-pixel vibration extraction is the atomic operation repeated across the whole tomographic line to build one tomogram (see Section 4).

## Key Takeaways
1. Only **one** SAR SLC image is required per tomogram — motion/vibration is extracted from Doppler sub-apertures of that single acquisition, not from a repeat-pass pair.
2. Three canonical SAR motion artifacts (azimuth displacement, azimuth smearing, range-walk) are the observable signatures the method deliberately measures rather than corrects away.
3. The displacement model uses a second-order Taylor expansion of the range history, dropping a cubic term, to linearize velocity/acceleration recovery.

## Connects To
- **Section 4 (Tomographic Model)**: takes the per-pixel vibration signal produced here as input to the harmonic-oscillator inversion.
- **Section 1.3**: cites the ship MM-estimation and bridge/dam structural-health-monitoring papers this sub-aperture strategy was originally developed for.

# Section 5.1: External Results — InSAR Fringe Analysis

## Core Idea
Classical repeat-pass InSAR interferometry (not the MM tomography method) applied to all three Giza pyramids' exterior faces reveals that each face is not flat but bows slightly inward toward the ground and is split into two indented half-faces — meaning each pyramid has 8 facets rather than the conventional 4.

## Frameworks Introduced
- **Interferometric fringe inclination analysis**: measure the tilt/inclination of interferometric phase fringes across a pyramid face; a change in fringe inclination partway up the face indicates a change in surface slope, revealing hidden facet boundaries invisible to the naked eye.
  - When to use: detecting subtle large-scale surface-shape deviations (bowing, faceting, settling) on any large static structure via repeat-pass SAR.
  - How: acquire two SAR images separated by a full orbital repeat cycle with a suitable spatial baseline; compute interferometric phase and coherence maps; read fringe spacing/inclination changes as slope discontinuities.

## Key Concepts
- **Repeat-pass interferometry**: two acquisitions separated by a fixed orbital cycle (here, 16 days for COSMO-SkyMed Second Generation) with a spatial baseline, used to generate interferometric fringes.
- **Coherence map**: quality measure of the interferometric pair; high coherence (reported near 1 outside radar-shadow areas) supports low-noise, reliable fringe measurements.
- **Radar shadow**: geometric occlusion (e.g., the south face of Menkaure) preventing measurement from a given acquisition geometry.

## Reference Tables
Acquisitions used for external results (Table 2 of the paper, external rows only, paraphrased):

| Pair | Orbit type | Beam | Polarization |
|---|---|---|---|
| 28 Oct / 13 Nov 2021 | Right-descending | 06 | HH |
| 27 Oct / 12 Nov 2021 | Right-descending | 08 | HH |
| 24 Jul / 9 Aug 2021 | Right-ascending | 39 | HH |

## Worked Example
For validation, the authors compared the SAR-derived displacement trend of one interferometric fringe on Khufu's west face against an independent LiDAR-measured profile along the same trajectory: both show the same trend, with a quantified maximum error of ~0.35 m (typical error ~0.1 m) between the two independent measurement methods.

## Key Takeaways
1. All three pyramids (Khufu, Kefren, Menkaure) were found to have 8 facets, not 4 — each face divided by a subtle inward bow near the base.
2. This finding used **interferometry** (needs an image pair), distinct from the **MM tomography** method (needs only a single image) used for the interior results in Sections 5.2+.
3. Cross-validation against LiDAR gives the paper's main independent accuracy check, with sub-half-meter agreement.
4. This 8-facet finding becomes a load-bearing premise in the Section 6.2 speculative water-basin interpretation (used to argue for a purpose-built water-shedding geometry).

## Connects To
- **Section 6.2**: cites the 8-facet finding as supporting evidence for the hydraulic/acoustic interpretation of pyramid function.
- **Section 5.2**: switches from interferometric (multi-image) to MM tomographic (single-image) technique for the interior.

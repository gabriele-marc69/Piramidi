# ESA Earthnet / Third Party Missions — Project Proposal (DRAFT)

> **Status: ready for submission — public repository published at
> https://github.com/gabriele-marc69/Piramidi .**
> Target submission channel: ESA EO Sign In → earth.esa.int/eogateway → Third Party
> Missions data access (project proposal). Parallel submissions worth considering:
> ASI COSMO-SkyMed Open Call (same science case), DLR TerraSAR-X Science Service.

---

## 1. Title

**Independent, Reproducible Validation of Micro-Motion SAR Doppler Tomography on
Large Static Structures: an Open-Source Test Case on the Giza Pyramids**

## 2. Principal Investigator

- Name: Gabriele Marchini
- Affiliation: Independent researcher, Italy
- E-mail: gabrielemarchini69@gmail.com
- Open-source project repository: https://github.com/gabriele-marc69/Piramidi
  (the full processing chain described in Sections 6–7, publicly released)

## 3. Executive summary (≤ 200 words)

Micro-Motion (MM) SAR Doppler tomography (Biondi & Malanga, *Remote Sensing* 2022,
14, 5231) claims depth-resolved imaging of the interior of large static structures
from a single SAR acquisition, by re-reading micro-Doppler defocusing artifacts as a
vibration signal and inverting it with a steering-matrix (beamforming) model. The
method is potentially transformative for non-invasive structural assessment
(monuments, dams, bridges), but no independent, fully reproducible verification
exists. We have built an open-source end-to-end pipeline and applied it to free
Sentinel-1 IW SLC stacks over the Giza pyramids. This established a quantitative
information budget showing that IW TOPS data are **intrinsically insufficient**: the
per-target Doppler dwell supports only ~12 independent looks, giving an ambiguity
height of ~8.5 m — far short of the ~140 m scale of the monuments. We request
high-dwell Spotlight (preferred) or Stripmap SLC acquisitions over the Giza plateau
to (a) validate the MM tomographic chain against the *known* internal structures of
the Great Pyramid, (b) publish the first independent, fully open assessment of the
method's real capabilities and limits. All code, intermediate products, and results
will be released openly.

## 4. Scientific background and objectives

### 4.1 Background

Biondi & Malanga (2022) split the Doppler bandwidth of a single SLC image into
sub-apertures, track each pixel sub-pixel-wise across them to extract a vibration
time series, model each pixel on a tomographic line as a damped harmonic
oscillator, and invert the multi-look data vector `Y` with a matched-filter
steering matrix, `h(z) ≈ Aᴴ(Kz, z)·Y`, at a nominal depth resolution
`δz = λ_sonic·R/(2A) ≈ 0.92 m`. Their catalog of internal structures of the Great
Pyramid is presented by the authors themselves as evidence pending field
confirmation. The method has not been independently reproduced; its published
results rest on high-dwell X-band data not openly available.

### 4.2 Why open Sentinel-1 IW data cannot settle the question (information budget)

Our pipeline quantifies the limit precisely. For the Khafre-pyramid ground box
(~26 azimuth × 81 range pixels, IW2 VV/VH, 6–12 scenes, one track):

| Quantity                                | Value       | Driver                                                                                                                            |
| --------------------------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Independent Doppler looks per target, k | ~12         | TOPS burst dwell (~0.5 s) caps per-target azimuth bandwidth; box supports ~13 usable spectral bins at the paper's 50 % guard band |
| Nominal depth resolution δz             | ~1.2 m      | δz = λ_sonic·R/(2A), λ_sonic = 0.24 m, R ≈ 850 km, A ≈ 84 km                                                                      |
| Ambiguity height z_amb                  | **≈ 8.5 m** | z_amb ∝ (k−1); depths alias modulo z_amb                                                                                          |
| Resolvable depth cells ≈ z_amb/δz       | ~7          | ≈ k−1, invariant under any re-parameterisation                                                                                    |
| Structure scale to resolve              | ~140 m      | Khafre height ≈ 136 m, base 215 m                                                                                                 |

The number of resolvable depth cells is bounded by the number of independent looks;
no re-grouping of parameters, multi-modulus fusion, or post-processing can exceed
it (we verified this on the pipeline's synthetic forward model: multi-k "puzzle"
disambiguation fails in 20/20 noise realisations even under idealised conditions).
Reaching monument scale at metric resolution requires **k of order 10²**, i.e.
Spotlight-class dwell, or a large Stripmap time series.

### 4.3 Objectives

1. **O1 — Blind validation on known structures.** Apply the open MM chain to
   high-dwell SLC data of the Great Pyramid (Khufu) and test whether it recovers the
   *documented* interior features (Grand Gallery, King's/Queen's chambers, relieving
   chambers) at their known positions, before any claim about unknown ones.
2. **O2 — Quantified capability envelope.** Measure, on real high-dwell data, the
   achieved z_amb, δz, cross-scene stability (our pipeline computes per-voxel
   cross-scene coefficient of variation), and false-alarm behaviour vs. the
   synthetic forward model.
3. **O3 — Cross-method consistency.** Compare against published independent
   evidence (2017 muon spectroscopy, thermal and microgravimetric surveys),
   explicitly reporting non-matches as well as matches.
4. **O4 — Open release.** Publish code, processing parameters, and derived products
   so any group can reproduce every figure.

We deliberately do **not** frame this project as a search for undiscovered chambers.
If the method fails validation on known structures, that negative result is equally
publishable and equally valuable to the EO community.

## 5. Data requirements

| Item                                      | Request                                                                                                                   |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Missions (any of, in order of preference) | COSMO-SkyMed / COSMO-SkyMed SG (Spotlight-2), TerraSAR-X / TanDEM-X (Staring/High-Res SpotLight), PAZ (SpotLight)         |
| Fallback                                  | Stripmap SLC time series (TerraSAR-X SM, CSK Himage); Sentinel-1 SM on-demand from L0 where archive acquisitions exist    |
| Product level                             | SLC (complex, single look) — mandatory; no detected/geocoded products                                                     |
| Polarisation                              | Single co-pol sufficient (HH or VV); dual-pol welcome                                                                     |
| Area of interest                          | Giza plateau, Egypt — centre ≈ 29.976° N, 31.132° E; footprint ≥ 2 × 2 km covering Khufu and Khafre                       |
| Number of acquisitions                    | ≥ 6 same-geometry passes (new tasking or archive), to enable cross-scene stability statistics; archive-only is acceptable |
| Time span                                 | Any; ≥ 3 months between first and last pass preferred                                                                     |
| Estimated volume                          | ≤ 20 scenes, ≈ 40 GB total (Spotlight/Stripmap SLC ≈ 1–2 GB per scene)                                                    |

**Justification of mode choice:** the tomographic information content scales with
per-target Doppler dwell (Section 4.2). Spotlight dwell (seconds) supports the
k ≈ 10² independent sub-apertures needed for monument-scale unambiguous depth;
TOPS IW (dwell ~0.5 s) demonstrably does not. This is the single parameter that
open Copernicus data cannot provide, and the only reason third-party data are
requested.

## 6. Methodology and work plan

1. **WP1 — Ingestion** (month 1–2): extend the existing open pipeline (Sentinel-1
   reader → co-registered complex box) to the granted mission's SLC format.
2. **WP2 — Synthetic calibration** (month 2–3): re-derive the information budget
   for the actual acquisition parameters; verify the inversion chain end-to-end on
   the forward model (`genera_Y_sintetico`) with reflectors at known depths.
3. **WP3 — Blind validation** (month 3–8): O1 on Khufu; pre-registered feature
   positions from the literature; cross-scene CV maps as the stability metric.
4. **WP4 — Capability envelope & cross-method comparison** (month 8–11): O2, O3.
5. **WP5 — Open release & reporting** (month 11–12): code, products, final report,
   peer-reviewed submission.

Processing is fully implemented and public: Doppler sub-aperture decomposition,
sub-pixel tracking, steering-matrix beamforming (plus Capon/MVDR), synthetic
forward-model verification, incoherent multi-scene averaging with per-voxel
cross-scene stability maps.

## 7. Feasibility (previous work)

- Working end-to-end open-source pipeline on Sentinel-1 IW SLC (search → download →
  co-registered complex box → DInSAR → MM tomographic volume and B-scans).
- Synthetic forward model validating the inversion (reflectors at known depths
  recovered within z_amb; alias behaviour characterised).
- Quantified negative result on IW data (z_amb ≈ 8.5 m) — the analysis that
  motivates this request and demonstrates the team understands the method's
  physics rather than over-claiming from insufficient data.

## 8. Expected results and deliverables

- D1: Validation report — does MM tomography recover Giza's known interior on
  high-dwell data? (positive or negative, with uncertainty quantification)
- D2: Open-source processing chain for the granted mission's SLC format.
- D3: Derived products (tomographic volumes, stability maps) under an open licence
  compatible with the mission's data policy.
- D4: Peer-reviewed publication acknowledging ESA/TPM data provision.

## 9. Data handling and publication commitment

Original third-party products will be used and stored per the mission's licence and
not redistributed; only derived products and code are released. Results will be
reported faithfully, including negative or inconclusive outcomes, and all published
depth claims will carry the applicable ambiguity/resolution bounds.

## 10. References

1. Biondi, F.; Malanga, C. *Synthetic Aperture Radar Doppler Tomography Reveals
   Details of Undiscovered High-Resolution Internal Structure of the Great Pyramid
   of Giza.* Remote Sens. 2022, 14, 5231. (arXiv:2206.09200)
2. Morishima, K. et al. *Discovery of a big void in Khufu's Pyramid by observation
   of cosmic-ray muons.* Nature 2017, 552, 386–390.
3. Procureur, S.; Morishima, K. et al. *Precise characterization of a
   corridor-shaped structure in Khufu's Pyramid by observation of cosmic-ray
   muons.* Nat. Commun. 2023, 14, 1144.
4. Bui, H.D. *Imaging the Cheops Pyramid.* Solid Mechanics and Its Applications,
   Vol. 182; Springer, 2012. (EDF/CPGF microgravimetry campaigns, 1986–1987)
5. ScanPyramids mission (HIP Institute / Ministry of Antiquities, Egypt).
   Infrared thermography campaign reports, 2015–2016 (thermal anomalies at the
   NE ground-level edge of Khufu).
6. Project repository: https://github.com/gabriele-marc69/Piramidi

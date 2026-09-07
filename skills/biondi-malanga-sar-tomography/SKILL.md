---
name: biondi-malanga-sar-tomography
description: "Knowledge base from three primary sources by Filippo Biondi (with Corrado Malanga on the Giza paper): the paper \"Synthetic Aperture Radar Doppler Tomography Reveals Details of Undiscovered High-Resolution Internal Structure of the Great Pyramid of Giza\" by Filippo Biondi & Corrado Malanga (Remote Sensing, MDPI, 2022), the underlying patent WO 2024/008365 A1, and the originating preprint arXiv:2206.09200 \"Scanning Volcanoes by Synthetic Aperture Radar\" (Vesuvius case study + three-layer validation). Use when applying the Doppler sub-aperture / Micro-Motion (MM) SAR tomography method, tuning its parameters (B_shift, N_D, guard band, depth resolution), implementing its 11-block processing pipeline, interpreting the Giza pyramid interior-structure findings, or referencing SAR interferometry vs. tomography technique choices for large static structures."
---

<!-- argument-hint: [topic, framework name, or section number] -->

# SAR Doppler Tomography — Giza Paper, WO 2024/008365 Patent & Vesuvius Preprint
**Authors**: Filippo Biondi, Corrado Malanga (Giza paper) | Filippo Biondi (patent, sole inventor; preprint, sole author)
**Sources**: Remote Sensing (MDPI) 2022, 14, 5231 (44 pp.) + WIPO PCT WO 2024/008365 A1, publ. 11 Jan 2024 (30 pp.) + arXiv:2206.09200v2 [eess.SP], 18 Jul 2022 (22 pp.)
**Sections**: 16 | **Generated**: 2026-08-26 | **Verified against the source PDFs**: 2026-09-07 (see *Errata in the Sources*)

## How to Use This Skill

- **Without arguments** — load the core frameworks below for reference
- **With a topic** — ask about `B_shift`, `doppler sub-apertures`, `tomographic resolution`, `computational scheme`, `patent claims`, `discovered structures`, `insar fringes`, or another indexed topic; I find and read the relevant section file
- **With a section** — ask for `ch11` or "the claims section"; I load that specific file
- **Browse** — ask "what sections do you have?" to see the full index

Ch 1–9 come from the Giza paper; Ch 10–13 from the patent; Ch 14–16 from the Vesuvius preprint. The three sources use different parameter regimes rather than contradicting each other — see the three-regime table in [cheatsheet.md](cheatsheet.md).

---

## Core Frameworks & Mental Models

**Micro-Motion (MM) Doppler Tomography** — the central method. Reframes SAR micro-Doppler artifacts (normally treated as defocusing noise) as an information carrier: split a single SAR SLC image's Doppler bandwidth into sub-apertures, sub-pixel-track a target pixel across them to recover a time-domain vibration signal, then invert that vibration along a chosen depth line to produce a tomogram. Use this whenever you need to infer internal structure of a large static solid (pyramid, dam, bridge, volcano) non-invasively from a *single* SAR pass.

**Photon-versus-phonon framing** — EM waves cannot penetrate solid matter for kilometres, but vibrations do. The surface is not the target; it is the read-out membrane of an underground acoustic field. The illumination source is the Earth's own "seismic perpetual ripple" plus man-made vibration, not the radar. The patent calls the result **Space-Sonar** (see ch10).

**Doppler Sub-Aperture Decomposition with a guard band** — split total Doppler bandwidth `B_cD` into `N_D` sub-bands, withholding `B_DL = B_cD/2` around the matched-filter boundary to preserve motion sensitivity. Focusing the whole band restores azimuth resolution and destroys the measurement (see ch03, ch11).

**`B_shift` as the frequency selector** — master and slave sub-bands are held rigidly separated by `B_shift` along the azimuth-frequency axis; `B_shift` *is* the vibrational frequency you observe, in inverse relation (higher `B_shift` gives lower mechanical frequency). `N_D` is not a compute parameter but the **digital sampling rate of the mechanical wave**. `B_shift` appears **only in the patent**; `N_D` and the guard band `B_DL` *are* defined in the paper (§3.1, eqs 4–5), but the paper never says what `N_D` physically is — that reading is the patent's (see ch11).

**Controlled linearization** — the true model is a nonlinear Duffing oscillator with cubic restoring force; it reduces to the tractable damped 2-DOF linear oscillator `r(t) = (a·cos ω0t, b·sin ω0t)·exp(−λt/2)` **only when nonlinearity is low**, i.e. when the spring is meaningfully tensioned (`L` well away from `L0`). The `{a, b}` are exactly the shifts the coregistrator measures (see ch12).

**Tomographic inversion as pulse compression** — `Y = A(K_z,z)·h(z)`, solved as `h(z) = A^H·Y`. The steering matrix is explicitly the best approximation of a DFT operator, so depth focusing *is* Fourier focusing. Resolution `δ_z = λR/2A` where **λ is the sound wavelength in the medium, not the radar wavelength** — about 36 m (Vesuvius), 0.92 m (Giza paper) or 1.30 m (patent) depending on the parameter regime (see ch04, ch12, ch14). The Giza paper's 0.92 m does **not** follow from the parameters that paper itself states — see *Errata in the sources* below.

**The 11-block pipeline** — `SLC → FFT2 → {BPF master → IFFT2, BPF slave → IFFT2} → pixel tracking → raw seismic data → FFT2 focusing → tomographic map → filter + geolocation`. Blocks 8 and 10 are data states, not operations. Compute the forward FFT2 once; only the two inverse transforms and the tracker run inside the `N_D` loop (see ch13).

**Depth-versus-resolution tuning** — `δ_z = λR/2A` with `λ = v/f`, so the investigation frequency is the knob you control for trading cell size against penetration: 200 Hz gives ~36 m cells and ~3 km reach at Vesuvius; 12.5 kHz gives metre-scale cells and pyramid-interior reach at Giza. It is not the *only* term that moves: the assumed propagation speed `v` differs by a factor ~6 across the corpus (972 m/s at Vesuvius vs. 6000–6600 m/s at Giza), and the synthesized aperture `A` by ~2. Choose `f` from the depth requirement before any processing (see ch12, ch14).

**Three-layer validation** — validate a novel imaging chain at three independent levels: geometry against a DEM, the raw measurement against a co-located ground sensor (with a plotted error, not an assertion), and the imaged structure against a different modality. Each can fail without the others; only all three constrain the whole chain (see ch15).

**Confirm features by changing a processing parameter** — a candidate structure that survives both a change of tomographic-line orientation and a change of investigation frequency is unlikely to be an artifact. Line orientation is a processing parameter, so this costs no new data (see ch14).

**Interferometric Fringe-Inclination Mapping** — separate from MM tomography: uses a *repeat-pass pair* to detect exterior shape deviations via fringe-inclination changes. Used to find that all three Giza pyramids have 8 facets, not 4. Choose interferometry when you have two well-baselined passes and want exterior shape, not interior structure (see ch05).

**Unique-Tag Cataloging** — assign a persistent ID to every candidate feature at first detection and carry it through all figures, tables and discussion. The paper catalogs 20 tagged interior structures this way, enabling traceable cross-referencing and self-flagging of at least one false alarm (see ch07).

**Cross-Method Consistency Checking (including non-matches)** — every finding is checked against prior independent methods (thermal anomaly studies, microgravimetry, 1998 electrogravitic measurements, 2017 muon spectroscopy) and disagreements are reported explicitly rather than omitted (see ch08).

**Epistemic framing** — the structural catalog is presented as evidence-based; the functional/hydraulic interpretation (Quincke's-tube water circuit, Helmholtz-resonator acoustics, therapeutic use) is labelled a starting hypothesis requiring excavation. Preserve that distinction. The same discipline applies to the patent: WO 2024/008365 is a *published application carrying a category-X search opinion against all 10 claims* — cite it as disclosure of a method, never as granted validity (see ch08–ch10).

---

## Section Index

| # | Title | Source | Key Frameworks |
|---|-------|--------|----------------|
| [ch01](chapters/ch01-intro-motivation.md) | Introduction & Motivation | Paper | MM Doppler tomography origin, micro-Doppler reframing |
| [ch02](chapters/ch02-giza-plateau-known-structures.md) | Giza Plateau & Known Structures | Paper | Baseline known-structure diagram |
| [ch03](chapters/ch03-methodology-doppler-mm.md) | Methodology — Doppler Sub-Aperture MM | Paper | Sub-aperture decomposition, motion-artifact model |
| [ch04](chapters/ch04-tomographic-model.md) | Tomographic Model | Paper | Harmonic-oscillator inversion, resolution formula |
| [ch05](chapters/ch05-external-results-insar.md) | External Results — InSAR | Paper | Fringe-inclination mapping, 8-facet finding |
| [ch06](chapters/ch06-internal-imaging-known-structures.md) | Internal Imaging — Known Structures | Paper | Overlap/non-overlap validation |
| [ch07](chapters/ch07-discovered-structures-catalog.md) | Discovered Structures Catalog | Paper | Unique-tag cataloging, metric determination |
| [ch08](chapters/ch08-discussion-interpretation.md) | Discussion — Analysis & Interpretation | Paper | Cross-method consistency, speculative hydraulic model |
| [ch09](chapters/ch09-conclusions-future-work.md) | Conclusions & Future Work | Paper | Validation path, epistemic caveats |
| [ch10](chapters/ch10-patent-scope-claims.md) | Patent Scope, Claims & Legal Frame | Patent | 10-claim set, Space-Sonar, ISR category-X finding |
| [ch11](chapters/ch11-patent-subaperture-formalism.md) | Doppler Sub-Aperture Formalism | Patent | B_shift, N_D, guard band, eqs (1)–(12) |
| [ch12](chapters/ch12-patent-oscillator-inversion.md) | Nonlinear Spring Model & Inversion | Patent | Duffing oscillator, steering matrix, eqs (13)–(24) |
| [ch13](chapters/ch13-patent-pipeline-applications.md) | 11-Block Scheme & Application Domains | Patent | Reproducible pipeline, volcano/undersea/pyramid geometries |
| [ch14](chapters/ch14-vesuvius-case-study.md) | Vesuvius Case Study — Parameters & Findings | Preprint | Depth/resolution trade, line orientation, layover test |
| [ch15](chapters/ch15-three-layer-validation.md) | The Three-Layer Validation Protocol | Preprint | DEM / in-situ seismograph / cross-modality validation |
| [ch16](chapters/ch16-provenance-tooling-caveats.md) | Provenance, Tooling & Preprint Caveats | Preprint | Priority chain, Malanga attribution, template residue |

## Topic Index

- **B_shift (frequency selection)** → ch11, cheatsheet
- **Claims / claim set** → ch10
- **Coherence map** → ch05
- **Complex structure (tags 9, 10)** → ch07
- **Computational scheme (11 blocks)** → ch13, cheatsheet
- **Doppler centroid anomaly** → ch01, ch03
- **Errata / source arithmetic** -> SKILL.md (Errata in the Sources), ch04, ch12, ch14, cheatsheet
- **Doppler sub-apertures** → ch03, ch04, ch11
- **Duffing / nonlinear oscillator** → ch12, ch14
- **False alarm (tag 19 area)** → ch07
- **Foreshortening / layover** → ch14
- **Figure inventory (0.1–0.8)** → ch13
- **Geocoding / GIS output** → ch10, ch13
- **Giza 8-facet finding** → ch05, ch08
- **Grand Gallery / King's / Queen's Chamber** → ch02, ch06
- **Guard band (B_DL)** → ch03, ch11
- **Harmonic oscillator model** → ch04, ch12
- **International Search Report / novelty** → ch10
- **K_z / steering matrix** → ch04, ch12
- **LiDAR cross-validation** → ch05
- **Master/slave pixel tracking** → ch11, ch13
- **Metric determination** → ch07
- **Micro-Motion (MM) tomography** → ch01, ch03, ch04, ch10
- **Motion artifacts (smearing, range-walking)** → ch11
- **Investigation frequency (depth/resolution trade)** → ch12, ch14, cheatsheet
- **Lava conduits / vent apertures** → ch14
- **Magnetotelluric cross-validation** → ch15
- **Muon spectroscopy discrepancy** → ch08
- **N_D (mechanical sampling rate)** → ch11
- **Patent / funding disclosure** → ch09, ch10
- **Priority chain / arXiv versions** → ch16
- **SRTM DEM validation** → ch15
- **SARPROZ / tooling** → ch16
- **Phonon / Space-Sonar** → ch10
- **Pulse compression** → ch12, ch13
- **SAR acquisition parameters** → ch01, ch11, cheatsheet
- **Structure catalog (20 tags)** → ch07
- **Tomographic resolution formula** → ch04, ch12, cheatsheet
- **Validation protocol (three layers)** → ch15
- **Vesuvius** → ch14, ch15
- **Volcano / undersea / under-ice domains** → ch13, ch14
- **Zed (monument)** → ch02, ch06, ch08

## Supporting Files

- [glossary.md](glossary.md) — all key terms with definitions
- [patterns.md](patterns.md) — reusable SAR/imaging and disclosure techniques with trade-offs
- [cheatsheet.md](cheatsheet.md) — decision rules, paper-vs-patent parameter tables, method-selection matrix

---

## Errata in the Sources

These are arithmetic and bookkeeping defects **in the primary documents**, verified against the PDFs. The skill reports the published numbers *and* the discrepancy; do not quietly propagate either one alone.

| Where | What the source says | The problem |
|---|---|---|
| Giza paper, §4 (δ_z derivation) | `λ = v/f ≈ 6000/25,000 ≈ 0.24 m` | the same sentence sets the investigation frequency to **12,500 Hz**; dividing by 25,000 is dividing by `2f`. With the stated `f`, `λ = 0.48 m`. |
| Giza paper, §4 (δ_z derivation) | "maximum orbital aperture equal to half the total length of the orbit, therefore about **42,000 m**" then `δ_z = 0.24·650,000 / (2·84,000)` | the denominator uses `A = 84,000 m`, twice the aperture the same sentence states. The Vesuvius preprint and the patent both write `2A` with their own stated `A` (42,000 and 75,000 m). |
| Giza paper, §4 vs. §5.2 | `δ_z ≈ 0.92 m`, but "each pixel corresponds to 1 m of spatial resolution" | two different resolution figures for the same results. |
| Giza paper, §4 vs. Table 1 | tomograms "synthesized at 24 kHz"; Table 1 lists Doppler bandwidth **22 kHz** | unreconciled. |
| Vesuvius preprint, §2.3 (δ_z derivation) | `λ = v/f ≈ 3500/200 ≈ 4.86 m` | `3500` is km/h; 4.86 m is `972/200`, i.e. the metres-per-second value. The printed fraction mixes units, the result is right. |
| Vesuvius preprint, §2.3 (δ_z derivation) | `δ_z = 4.86·650,000/(2·42,000) ≈ 36 m` | the fraction evaluates to 37.6 m; 36 m is a loose rounding. |
| Patent, [0004] vs. [0005] | withheld band called `B_C_L` in [0004], `B_D_L` in [0005] | same quantity, two symbols. Both are `B_cD/2`. |
| Paper §3.1 vs. patent [0005] | paper: "`N_c` rigid shifts of the master–slave system"; patent: "`N_D` rigid shifts" | the paper's `N_c` is never used again; the patent collapses it onto `N_D`. |

**Net effect on the headline number**: recomputing `δ_z = λR/2A` from the Giza paper's own stated `v = 6000 m/s`, `f = 12,500 Hz`, `A = 42,000 m`, `R = 650,000 m` gives **≈3.7 m**, four times coarser than the published 0.92 m. Every metric claim in the paper's §5 rests on the 0.92 m figure. Cite 0.92 m as *the value the paper used*, never as a value that follows from its parameters.

---

## Scope & Limits

This skill covers three primary sources only — one peer-reviewed paper, one patent application, and one arXiv preprint, all by the same lead author. It is a primary-source summary, not an independent verification of the archaeological claims (which the authors themselves describe as hypotheses pending field confirmation) nor of the patent's validity (whose international search report is adverse on all claims). The Vesuvius preprint is not peer-reviewed and carries unremoved LaTeX template filler in section 2.1 (see ch16). For broader SAR/InSAR theory or other Giza-pyramid research, consult additional sources.

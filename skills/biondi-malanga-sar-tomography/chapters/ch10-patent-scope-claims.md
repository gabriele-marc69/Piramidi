# Chapter 10: Patent WO 2024/008365 A1 — Scope, Claims & Legal Frame

*Source: WIPO PCT publication WO 2024/008365 A1, "Synthetic Aperture Radar Underground, Undersea, Underice, and Inside Distributed Targets Tomographic Doppler Imaging". Application PCT/EP2023/064345, filed 29 May 2023, published 11 Jan 2024. Priority PT4451, 04 July 2022 (MT). Sole inventor/applicant: **BIONDI, Filippo** (MT/IT). IPC: G01S 13/88, G01S 13/90, G01V 1/00.*

## Core Idea
The patent is the legal and methodological parent of the Giza paper: it claims the *general* method — extracting phonon/vibrational information from the Doppler content of a **single** SAR image to tomographically image the inside of any body — with the pyramid appearing only as one illustrative geometry among volcanoes, seabed, ice and man-made structures.

## Frameworks Introduced

- **"New observation domain" framing**: the invention is positioned not as a better SAR product but as an *added observation domain* bolted onto any existing SAR satellite — underground, undersea, under-ice, and inside solid bodies.
  - When to use: framing a processing-only invention that requires no new hardware.
  - How: claim the method (claim 1), then the platform-agnosticism (claim 2), then the acquisition-agnosticism (claim 3), then system/software embodiments (claims 4–6), then each processing stage separately (claims 7–10).

- **Photon-versus-phonon interaction**: the stated physical premise. Electromagnetic waves cannot penetrate solid matter for kilometres, but *vibrations* do; the EM wave interacts with the vibrating surface and thereby carries out sound information from depth.
  - When to use: justifying why depth information can exist in a surface-only EM measurement.
  - How: treat the surface as the read-out membrane of an underground acoustic field, not as the target itself.

- **Ambient-seismic-field exploitation**: the illumination source is not the radar. It is the Earth's own "seismic perpetual ripple" plus man-made vibration, always present and exploitable even at very low magnitude.
  - Consequence: the method is passive-acoustic / active-EM — a hybrid the authors call **Space-Sonar**.

## Key Concepts
- **Claim 1 (independent)** — a method to obtain complex tomographic imaging below/inside the Earth, ice, sea, and inside any subsurface or surface man-made object (buildings, bridges, dams) or any other object made of matter, *starting from a single or multiple SAR images*.
- **Claim 2** — SAR sensor carried by air/space/satellite platform.
- **Claim 3** — any frequency, any mode (ScanSAR, Stripmap, Spotlight), any polarization, any geometry (incidence angle, squinted).
- **Claim 4/5/6** — system, software product, and compiled program embodiments.
- **Claim 7** — chirp-Doppler sub-aperture refocusing of SAR data at any processing level.
- **Claim 8** — complex vibrational extrapolation by pixel-displacement measurement.
- **Claim 9** — pulse compression of estimated complex vibration data.
- **Claim 10** — complex tomographic extrapolation output at any GIS standard.
- **Space-Sonar** — the authors' name for the recast concept: sonar-like sounding performed from orbit via EM-matter interaction.
- **Processing level agnosticism** — accepts SLC *and/or raw* SAR data; the minimum is one image.

## Reference Tables

| Claim | Type | What it locks down |
|---|---|---|
| 1 | Method (independent) | Tomographic imaging inside/below matter from ≥1 SAR image |
| 2 | Dependent | Air/space/satellite platform |
| 3 | Dependent | Any frequency / mode / polarization / geometry |
| 4 | System | Apparatus implementing claims 1–3 |
| 5 | Software product | Loadable code portions implementing claims 1–4 |
| 6 | System program | Any programming/compiling language, claims 1–5 |
| 7 | Processing stage | Chirp-Doppler sub-aperture refocusing |
| 8 | Processing stage | Vibrational extrapolation by pixel displacement |
| 9 | Processing stage | Pulse compression of vibration data |
| 10 | Output | GIS-standard tomographic extrapolation |

## Anti-patterns

- **Publishing your own method before the priority date.** The International Search Report (EPO, completed 31 Aug 2023, mailed 08 Sep 2023, officer Rudolf, Hans) cites exactly one document — Biondi's own arXiv preprint *"Scanning Inside Volcanoes by Synthetic Aperture Radar Echography Tomographic Doppler Imaging"*, 18 June 2022, XP091253966 — as category **X against claims 1–10**, i.e. the whole claim set lacks novelty or inventive step over it taken alone. The arXiv posting predates the 04 July 2022 priority date by about two weeks. **Why it fails**: in the EPO/PCT system the applicant's own prior publication is prior art against them; there is no general grace period. Sequence the disclosure after the priority filing, not before.
- **Claiming "any" everywhere.** The claim set leans heavily on "any frequency, any mode, any polarization, any geometry, any programming language". Broad drafting is normal, but paired with an X-cited self-disclosure it leaves little fallback ground to narrow toward during prosecution.

## Worked Example — reading the ISR verdict

The ISR box C contains a single row:

```
Category: X
Citation:  FILIPPO BIONDI: "Scanning Inside Volcanoes by Synthetic Aperture
           Radar Echography Tomographic Doppler Imaging", ARXIV.ORG,
           18 June 2022 (2022-06-18), XP091253966, the whole document
Relevant to claim No.: 1-10
```

How to read it: **X** = "of particular relevance; the claimed invention cannot be considered novel or cannot be considered to involve an inventive step when the document is taken alone." Applied to *1-10* it covers every claim, independent and dependent. Only one document was found, and it is the inventor's own. So the examiner's objection is not "someone else got there first" — it is "you published it yourself, first." The practical consequence for anyone citing this patent: **WO 2024/008365 A1 is a published application with an adverse search opinion, not a granted patent.** Cite it as disclosure of method, never as evidence of granted, examined validity.

## Key Takeaways
1. The patent covers the *general* MM Doppler tomography method; Giza is an application, not the subject.
2. Filippo Biondi is the sole inventor/applicant — Corrado Malanga is not named on the patent, only on the Remote Sensing paper.
3. The minimum input is one SAR image at any processing level (SLC or raw) — this single-image property is the invention's core distinguishing feature versus in-situ sounding systems and versus InSAR.
4. Claims 7–10 decompose the pipeline into separately-claimed stages: sub-aperture refocusing → vibrational extrapolation → pulse compression → GIS output.
5. The EPO search report is category X against all 10 claims over the inventor's own arXiv preprint — treat published-application status accordingly.
6. Priority 04.07.2022 sits *after* the 18.06.2022 arXiv posting; the two-week gap is the whole novelty problem.

## Connects To
- **Ch 9**: the paper's patent/funding disclosure — this is the patent it points at.
- **Ch 11–13**: the technical body of this same document.
- **Ch 3, Ch 4**: the paper's methodology and tomographic model, here given in fuller mathematical form.

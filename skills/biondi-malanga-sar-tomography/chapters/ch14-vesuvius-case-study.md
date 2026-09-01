# Chapter 14: Vesuvius Case Study — Parameters & Findings

*Source: arXiv:2206.09200v2 [eess.SP], "Scanning Volcanoes by Synthetic Aperture Radar", Filippo Biondi (University of Strathclyde), 18 July 2022. This is v2; the version cited as category-X prior art against the patent is the 18 June 2022 posting — see ch16.*

## Core Idea
The Vesuvius study is the method's first published application and its most permissive regime: by dropping the investigation frequency to 200 Hz the tomogram reaches ~3 km depth at ~36 m resolution — the mirror image of the Giza parameter choice, and the clearest demonstration that depth and resolution are traded against each other through a single knob.

## Frameworks Introduced

- **Depth-versus-resolution trade via investigation frequency**: `δ_z = λR/2A` with `λ = v/f`. Lowering `f` lengthens the acoustic wavelength, which coarsens `δ_z` but buys penetration.
  - When to use: at the very start of an acquisition plan, before any processing.
  - How: decide the depth you must reach, pick `f` accordingly, accept the resulting `δ_z`. Vesuvius: 200 Hz, 3 km depth, 36 m cells. Giza: 12.5 kHz, pyramid interior, 0.92 m cells.
  - Why it works: the sound wavelength in the medium sets the resolution floor; you cannot have metre-scale cells and kilometre-scale reach from the same tomogram.

- **Tomographic line orientation as an experimental variable**: the study runs the same scene through range-oriented, azimuth-oriented, oblique, and "inclined to span the crater diameter" lines, and reports what each reveals.
  - When to use: whenever a single line gives an ambiguous or partial answer.
  - How: re-run the tomogram along a differently-oriented line on the same SLC image. No new acquisition is needed — line choice is a processing parameter, not an acquisition one.

- **Layover/foreshortening robustness test**: one line is deliberately oriented pure-range to stress the reconstruction under heavy layover and foreshortening. Result: the layover side (eastern crater) images *better* than the foreshortening side (western).
  - Use as a siting rule: prefer geometries that put your target on the layover side.

- **Non-resonating material as a feature class**: vent conduits are identified as "material that does not resonate" — absence of vibrational energy is the detection signature, not presence.
  - Complementary tell: dense material shows as *accumulated* vibrational energy (the crater cap).

## Reference Tables

### Acquisition and processing parameters (Vesuvius)

| Parameter | Value |
|---|---|
| Satellite | COSMO-SkyMed Second Generation (CSG) |
| Mode / polarization | Spotlight-2A, HH |
| Doppler band | ~22.5 kHz |
| Chirp band | ~450 MHz |
| SLC Doppler synthesis for tomography | 24 kHz |
| SAR acquisition duration | ~14 seconds |
| Seismic propagation speed `v` | ~3500 km/h ≈ 972 m/s |
| Investigation frequency `f` | 200 Hz |
| Sound wavelength `λ = v/f` | ≈ 4.86 m |
| Orbit aperture `A` | ~42 000 m (half total orbit length) |
| Slant range `R` | 650 000 m |
| **Tomographic resolution `δ_z`** | **≈ 36 m** |
| Depth of investigation | ~3 km below peak topographic height |
| Acquisition month | February 2022 |

### The three parameter regimes across the corpus

| Quantity | Vesuvius (ch14) | Giza paper (ch04) | Patent (ch12) |
|---|---|---|---|
| `v` | ~972 m/s | ~6000 m/s | ~6600 m/s |
| `f` | 200 Hz | 12 500 Hz | ~22 000 Hz |
| `λ` | 4.86 m | not stated | 0.30 m |
| `A` | 42 000 m | not stated | 75 000 m |
| `δ_z` | ~36 m | ~0.92 m | ~1.30 m |
| Reach | ~3 km | pyramid interior | "several km" |

Note the assumed propagation speed also differs by a factor of ~6 between the volcanic edifice and the limestone/granite pyramid — the medium, not just the frequency, sets `λ`.

## Key Concepts
- **Cap / plug** — a vibrational-energy singularity below the main crater, read as denser material blocking the main conduit; also seen as "cooled compact lava" forming a blockage mass.
- **Magma chamber** — visible at the bottom of the tomograms of the eastern slope.
- **Lava tube splitting** — the deepest part of the imaged conduit shows the tube dividing.
- **Vent apertures** — 11 candidate lava apertures editable from the western/external side; 5 of them inside the main crater; apertures 9, 10, 11 possibly extending outside it.
- **Soil layering** — layered material on the inner western crater slope, resolved down to ~2 km relative to the top surface.
- **Secondary conduit** — a second lava tube alongside the main (currently plugged) one.

## Worked Example — reading a Vesuvius tomogram

Take Figure 13 (a)/(b), the pure-range line chosen as the robustness stress test:

1. **Line placement** — a full range line spanning near-range to far-range across the crater, drawn on the SLC magnitude image (Fig. 13 a).
2. **The tomogram is split into three labelled portions** (red boxes 2, 3, 4), each traced back to its segment on the tomographic line:
   - *Portion 2* — the eastern crater side, the one under **layover**. Best imaged of the three.
   - *Portion 3* — the interior of the main crater.
   - *Portion 4* — the descending side of the volcanic cone.
3. **Reading portion 3** — a substantial blockage mass appears deep inside the main volcanic conduit (box 5): high vibrational energy relative to background, interpreted as rock-like cooled compact lava.
4. **Reading the western side** — under **foreshortening**, the reconstruction is visibly poorer than the layover side. The conclusion drawn is geometric, not physical: possible lava conduits remain visible on the eastern side only.
5. **Cross-check the same feature on a second line** — the oblique line of Figure 10 (b) shows the same cap inside yellow circle 2, and the main conduits rising from east and west (arrows 3, 4). Agreement across two independent line orientations is what upgrades a bright patch to a reported feature.
6. **Frequency cross-check** — the same result is reproduced by recomputing the sonic tomogram at a *lower* vibrational frequency (Fig. 11). A feature that survives a change of investigation frequency is unlikely to be a processing artifact.

## Anti-patterns
- **Reading vibrational brightness as material identity.** Energy accumulation says "this vibrates more than background", which is *consistent with* denser material — the paper's own hedged phrasing ("could be associated with", "probably", "possible"). Keep the hedge.
- **Comparing a 14-second SAR snapshot against month-long in-situ records as if they were the same measurement.** The paper flags this explicitly when correlating with the INGV earthquake catalog (see ch15).
- **Trusting the foreshortened side of a range-oriented tomogram.** The study demonstrates it is measurably worse; treat features there as unconfirmed.

## Key Takeaways
1. `δ_z` at Vesuvius is ~36 m — about 40x coarser than Giza's ~0.92 m, by design, to reach 3 km depth.
2. Investigation frequency is the depth/resolution knob; the medium's propagation speed sets the rest.
3. Line orientation is a *processing* parameter — multiple orientations from one SLC image, no re-acquisition.
4. Layover images better than foreshortening for this method; site the target accordingly.
5. Detection works in both directions: energy accumulation (dense cap) and energy absence (non-resonating vent conduits).
6. Features are confirmed by surviving both a change of line orientation and a change of investigation frequency.
7. Atmospheric phase delay is argued away as time-invariant across the ~14 s acquisition, and the Doppler-domain scan within a single image is claimed robust to it regardless.

## Connects To
- **Ch 12, Ch 04**: the same `δ_z = λR/2A` formula at two other parameter settings.
- **Ch 13**: this is the volcano geometry of the patent's Figure 0.6, worked out on real data.
- **Ch 15**: the validation layers applied to these tomograms.
- **Ch 16**: provenance of this preprint and its relationship to the patent's search report.
- **Ch 07**: the Giza catalog, the same "tag and report every candidate" discipline applied to lava apertures.

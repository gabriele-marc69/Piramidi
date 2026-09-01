# Section 1: Introduction & Motivation

## Core Idea
No SAR/microwave technique had non-invasively produced a *complete, metrically-measured 3D map* of the Great Pyramid of Khufu's interior; prior methods (muon tomography, microgravimetry, georadar) gave hints of voids but not full geometry — this paper's Micro-Motion (MM) Doppler tomography aims to close that gap.

## Frameworks Introduced
- **Micro-Motion (MM) Doppler Tomography**: reuse of a ship/bridge/volcano vibration-monitoring SAR technique, repurposed to image *inside* a solid megalithic structure instead of measuring surface vibration only.
  - When to use: any large, static, non-metallic solid body where you want subsurface structure without physical access (bridges, dams, volcanoes, monuments).
  - How: treat the SAR-observed micro-Doppler vibration signature of a target as an information-bearing artifact (not noise), then invert it tomographically along a chosen line of sight.

## Key Concepts
- **Micro-Doppler effect**: extra Doppler shift caused by target micro-movement (vibration/rotation) superimposed on the bulk Doppler shift from platform motion; classically a "defocusing artifact" in SAR imaging.
- **Doppler centroid anomaly**: deviation in the Doppler centroid frequency used during azimuth focusing, exploited here as the vibration signal carrier.
- **Multi-Chromatic Analysis (MCA)**: technique (coined in prior work) performed in the range direction to retrieve unambiguous height/phase information from a single or paired SAR image(s).
- **Stop-and-go approximation breakdown**: when a target moves/vibrates during the SAR observation, the standard "stop-and-go" imaging assumption fails, producing exploitable artifacts (azimuth smearing, range-walk).
- **Persistent Scatterer Interferometry (PSI)**: prior-generation SAR structural-health-monitoring technique the authors build on for bridges/dams before extending it to tomography.

## Mental Models
- Think of SAR imaging artifacts (defocusing, smearing) as *encoded information* rather than noise to be filtered out — the paper's central reframing.
- Use "poor EM penetration into solid rock is a known SAR limitation" as the baseline problem this paper claims to sidestep by measuring *vibration* (acoustic/phonon-like) rather than direct electromagnetic backscatter through the material.

## Anti-patterns
- **Relying on destructive confirmation**: earlier "definitive" claims about the pyramid's interior (e.g. muon-detected voids) could not be drilled to confirm; the authors flag independent non-destructive cross-checks as necessary before treating any single method's output as ground truth.

## Key Takeaways
1. The MM/Doppler tomography approach is adapted from ship micro-motion estimation and infrastructure (bridge/dam) health monitoring — the same core math is reused across very different targets.
2. The claimed advantage over georadar/muon methods is deeper penetration via vibration, not electromagnetic backscatter, at the cost of very heavy compute.
3. Two research fronts are combined: (a) external InSAR fringe analysis of all 3 Giza pyramids, (b) internal MM tomography of Khufu's interior alone.
4. SAR acquisition parameters used: 9.6 GHz center frequency, 400 MHz chirp bandwidth, 22 kHz Doppler bandwidth, spotlight mode, HH/VV polarization, COSMO-SkyMed satellite (see cheatsheet for full table).

## Connects To
- **Section 3 (Methodology)**: formalizes the Doppler sub-aperture strategy introduced conceptually here.
- **Section 6 (Discussion)**: revisits the muon/microgravimetry literature cited here to argue consistency with the new results.

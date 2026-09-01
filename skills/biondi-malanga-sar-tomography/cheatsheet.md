# Cheatsheet — Decision Rules & Quick Reference

## Decision Rules

- **When you have only one SAR acquisition** -> use Doppler sub-aperture MM tomography (Ch 3-4, Ch 11-12), not InSAR (needs a pair).
- **When you have two repeat-pass acquisitions with a good baseline** -> use interferometric fringe analysis (Ch 5) for exterior shape/faceting questions.
- **When you must choose which mechanical frequency to observe** -> set B_shift. Higher B_shift = lower observed frequency (Ch 11).
- **When choosing N_D** -> treat it as the sampling rate of the mechanical wave and apply Nyquist to the vibration you want, not to your compute budget (Ch 11).
- **When focusing master and slave** -> always withhold the guard band B_DL = B_cD/2. Focusing the full Doppler band maximizes azimuth resolution and destroys motion sensitivity (Ch 11).
- **When computing depth resolution** -> lambda in delta_z = lambda*R/(2*A) is the SOUND wavelength in the medium, never the radar wavelength (Ch 12).
- **When applying the linear 2-DOF oscillator** -> first check that L is not close to L0; that is the regime where the cubic term dominates (Ch 12).
- **When implementing the pipeline** -> compute the forward FFT2 once (block 2) and copy to both branches; only the two IFFT2s and the tracker run inside the N_D loop (Ch 13).
- **When one tomographic line misses a feature** -> add an acquisition geometry, do not add processing. Line orientation gates visibility (Ch 13).
- **When a tomographic feature looks unusually bright/large near a pyramid edge or corner** -> check for a geometric false-alarm cause (e.g. ascending-angle artifact) before cataloging it (Ch 7).
- **When a new finding overlaps prior-literature claims (thermal, gravimetric, muon, etc.)** -> explicitly check for both agreement and disagreement; report non-matches (Ch 8).
- **When presenting a new imaging method's output** -> pair an overlapped-on-ground-truth view with a blind/non-overlapped view (Ch 6).
- **When claiming a functional/interpretive narrative beyond raw measurements** -> label it explicitly as hypothesis pending physical validation (Ch 8-9).
- **When citing WO 2024/008365** -> cite it as a published application disclosing a method, never as a granted or validated patent (Ch 10).
- **When you intend to both publish and patent** -> file priority first, post the preprint second (Ch 10).
- **When planning an acquisition** -> pick the investigation frequency from the depth you must reach; accept the resulting cell size. 200 Hz buys 3 km at 36 m; 12.5 kHz buys metre cells and no depth (Ch 14).
- **When you doubt a feature is real** -> recompute along a different line orientation AND at a different investigation frequency; survival of both is the confirmation test (Ch 14).
- **When siting a target in the SAR geometry** -> prefer the layover side over the foreshortening side; the Vesuvius stress test showed layover reconstructs measurably better (Ch 14).
- **When validating a novel imaging chain** -> validate geometry, measurement and structure separately, and say which layers are quantitative and which are only visual (Ch 15).
- **When citing the volcano preprint** -> name the version: v1 (18 Jun 2022) is the patent's prior art, v2 (18 Jul 2022) is the retitled document (Ch 16).

## Thresholds & Defaults

### SAR acquisition (Giza paper, CSG satellite)

| Parameter | Value |
|---|---|
| SAR center frequency | 9.6 GHz |
| Chirp bandwidth | 400 MHz |
| Doppler bandwidth | 22 kHz |
| PRF | 2.0 kHz |
| Antenna length | 6 m |
| Acquisition mode | Spotlight, HH/VV |
| Platform velocity | 7 km/s |
| Observation height | 650 000 m |
| Repeat-pass cycle (InSAR) | 16 days |
| Full-resolution tomogram compute time | ~6 days (i7, 32 GB RAM) |
| SAR/LiDAR external validation error | ~0.1 m typical, ~0.35 m max |

### Tomographic parameters — three regimes

| Quantity | Vesuvius preprint | Giza paper (RS 2022) | Patent (WO 2024/008365) |
|---|---|---|---|
| Target | volcano edifice | pyramid | generic |
| Vibration propagation speed v | ~972 m/s (3500 km/h) | ~6000 m/s | ~6600 m/s |
| Investigation frequency f | 200 Hz | 12 500 Hz | ~22 000 Hz |
| Sound wavelength lambda = v/f | ~4.86 m | not stated | ~0.30 m |
| Slant range R | 650 000 m | 650 000 m | 650 000 m |
| Orbit aperture A | 42 000 m | not stated | 75 000 m (half orbit) |
| **Depth resolution delta_z** | **~36 m** | **~0.92 m** | **~1.30 m** |
| Reach | ~3 km | pyramid interior | several km (axis to -3 km) |
| Satellite / mode | CSG spotlight-2A, HH | CSG spotlight, HH/VV | any (claim 3) |
| Doppler band | ~22.5 kHz (synth. 24 kHz) | 22 kHz | 22 kHz |
| Chirp band | ~450 MHz | 400 MHz | any |

Rule of thumb: delta_z scales as v/(f) x R/(2A). Depth reach and cell size move in opposite directions through f.

### Sub-aperture processing constants

| Symbol | Rule |
|---|---|
| B_DL | = B_cD / 2 (always withheld) |
| Master focus band | B_cr x (B_cD - B_DL) |
| Doppler step per shift | (B_cD - B_DL) / N_D |
| B_shift | = selected vibrational frequency (inverse relation) |
| N_D | = mechanical sampling rate |
| delta_D | ~ 1/B_cD = lambda*R / (2*L_sa) |
| K_z | = 4*pi*B_perp / (lambda * r_i * sin(theta)) |

## Method Selection Matrix

| Question | Method | Data needed |
|---|---|---|
| Exterior shape/faceting | Interferometric fringe analysis | 2 SAR passes, good baseline |
| Interior structure/voids | MM Doppler tomography | 1 SAR SLC (or raw) image per line |
| Which mechanical frequency | B_shift selection | Same single image |
| Cross-check exterior displacement | LiDAR comparison | Independent LiDAR survey |
| Validate known interior features | Overlap/non-overlap tomogram rendering | Known architectural schematic |
| Reach a feature one line misses | Second acquisition geometry | Additional satellite/pass |

## Pipeline at a Glance (Fig. 0.5)

`1 SLC -> 2 FFT2 -> {3 BPF master -> 5 IFFT2, 4 BPF slave -> 6 IFFT2} -> 7 pixel tracking -> 8 raw seismic data* -> 9 FFT2 focusing -> 10 tomographic map* -> 11 filter + geolocation`
(*8 and 10 are data states, not computations.)

## Tells & Smells

- A "big void" claim from one modality (e.g. muon) not reproduced by another (e.g. SAR) at the same location -> flag as an open discrepancy, don't silently pick a winner.
- A structure described in prior literature only via "intuitions not confirmed by objective evidence" -> corroboration from direct measurement is meaningfully stronger, still not definitive without field confirmation.
- Heavy compute time (multi-day) per result -> expect the survey to be built incrementally across many acquisitions, not from one pass.
- Raw vibrational magnitude that looks like structureless oscillation (Fig. 0.8 a, c) -> that is the expected pre-compression state, not a failed measurement.
- Azimuth displacement / azimuth smearing / range-walking in a focused image -> tells you which motion component (range velocity / azimuth velocity or range acceleration / range speed) dominates.
- A depth resolution that looks ~10x too good -> someone used the radar wavelength instead of the sound wavelength in delta_z = lambda*R/(2*A).
- A tomogram feature that disappears when you rotate the tomographic line -> processing artifact, not structure.
- A validation section with overlays but no error plot -> visual agreement only; weigh it below a quantitative check.

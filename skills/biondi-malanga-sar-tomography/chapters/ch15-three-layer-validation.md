# Chapter 15: The Three-Layer Validation Protocol

*Source: arXiv:2206.09200v2, section 4 "Validation of the Sonic Tomographic SAR Results" (4.1 Topographic, 4.2 Vibration, 4.3 Tomographic).*

## Core Idea
This is the strongest validation evidence anywhere in the corpus and the part the Giza paper has no equivalent of: the method is checked at three independent levels — the surface it should reproduce, the vibration signal it claims to measure, and the subsurface structure it claims to image — each against an outside instrument.

## Frameworks Introduced

- **Layered validation of a novel imaging modality**: validate the *geometry*, the *measurement*, and the *interpretation* separately, because each can fail without the others.
  - When to use: any time you publish results from an instrument or processing chain nobody else has run.
  - How: (1) show the output reproduces something already known and independently measured; (2) show the raw quantity you claim to measure agrees with a ground-truth sensor measuring the same quantity; (3) show the derived structure agrees with an independent modality imaging the same volume.
  - Why it works: layer 1 alone proves only that your surface geocoding is right. Layer 2 alone proves you measured vibration, not that you inverted it correctly. Only all three together constrain the whole chain.

- **Topographic self-consistency check (layer 1)**: the *top* of a subsurface tomogram must trace the actual ground surface. Extract the elevation profile along the exact tomographic line from an independent DEM and overlay it.
  - Instrument: SRTM (Shuttle Radar Topography Mission) global DEM, converted to slant-range coordinates.
  - Result: the surface component of the tomogram overlaps the DEM line "almost perfectly" across three separate line orientations.
  - Why it is the cheapest and most necessary check: if the tomogram's surface is misplaced, everything below it is misplaced by the same amount.

- **Direct measurement validation against in-situ sensors (layer 2)**: compare the radar-estimated vibrational streaming against seismographs measuring the same ground at the same time.
  - Instruments: two INGV seismograph stations, the only two both radar-visible and inside the imaged scene.
  - Comparison is made in three forms: full-bandwidth spectrum, low-pass filtered spectrum, and a ~1 kHz narrow-band detail; plus time-domain streaming over ~1 second of synchronized data, plus an explicit error plot.
  - Reported outcome: the error oscillates around zero and stays confined to very low values.

- **Cross-modality structural validation (layer 3)**: overlay the SAR tomograms on tomography from a physically unrelated technique, and separately correlate against an independent event catalog.
  - Instruments: magnetotelluric tomography of the same volume (50% overlay comparisons at three separate detail locations), plus the INGV earthquake catalog for Jan/Feb 2022 plotted in 3-D with magnitude-scaled radii.
  - The catalog check is deliberately weakened by the authors: a ~14-second SAR acquisition is not temporally comparable to month-long seismic monitoring; they claim only spatial and depth correspondence between February events and the estimated magmatic consistency.

## Reference Tables

### The three validation layers

| Layer | What is validated | Independent reference | Form of the check |
|---|---|---|---|
| 4.1 Topographic | Geometry / geocoding | SRTM DEM | Elevation profile overlaid on tomogram top, 3 line orientations |
| 4.2 Vibration | The raw measurement | INGV seismographs IV-VCRE, IV-VBKN | Spectrum (native / filtered / 1 kHz), time-domain streaming, error function |
| 4.3 Tomographic | The imaged structure | Magnetotelluric tomography + INGV earthquake catalog | 50% overlays at 3 detail sites; 3-D event plot vs. tomogram |

### In-situ stations used

| Station | Network | Channel | Location (Lat, Lon) |
|---|---|---|---|
| Vesuvius — East Crater | IV-VCRE | HH | 40.818999, 14.431419 |
| Vesuvius — North Bunker | IV-VBKN | HH | 40.829959, 14.429881 |

All other INGV stations were excluded — either outside the SAR scene or not in radar line-of-sight. Data source: `https://eida.ingv.it/it/`.

### Earthquake catalog overlay conventions (Figure 25)

| Element | Meaning |
|---|---|
| Peg position | Actual 3-D location of the event (lat, lon, depth) |
| Circle radius | Richter magnitude |
| Scale 0 / 1 / 2 / 3 | Magnitude 0–1 / 1–2 / 2–3 / >3 |
| Blue circles | Events throughout January 2022 |
| Red circles | Events in February 2022 — the SAR acquisition month |

## Worked Example — the vibration validation, end to end

This is the layer worth reproducing, because it is the one that tests whether the method measures anything real at all:

1. **Select comparable sensors.** Of the INGV network at Vesuvius, keep only stations that are both inside the SAR scene *and* in radar line-of-sight. Two survive: IV-VCRE (east crater) and IV-VBKN (north side).
2. **Locate them in both domains.** Show each station on the optical image and on the SAR magnitude image (arrows 1 in Figs. 17 and 18) — proof that the pixel whose vibration you extracted is the pixel the seismograph sits on.
3. **Synchronize in space and time.** Align the radar-estimated vibrational stream with the in-situ stream over roughly 1 second of overlapping data.
4. **Compare in the frequency domain, three ways** (Fig. 19 for IV-VCRE, Fig. 21 for IV-VBKN): native full-bandwidth spectrum; low-pass filtered spectrum; a ~1 kHz narrow-band detail. Blue = in-situ, brown = radar.
5. **Compare in the time domain** (Fig. 20 a/b): native and filtered synchronized streaming.
6. **Plot the error explicitly** (Fig. 20 c): unfiltered error in blue, filtered in brown. Report that it oscillates around zero and stays small — do not merely assert agreement.
7. **State the residual limitation.** The authors concede that perfectly synchronizing a ground instrument with a space instrument is very difficult, and characterize the achieved synchronization as "very acceptable" rather than exact.

What makes this convincing: both instruments are measuring the *seismic background* of the same area — a signal neither of them controls. The agreement is not between a model and its own output.

## Anti-patterns
- **Validating only the surface.** A perfect DEM overlay proves geocoding, nothing about depth. It is necessary and nowhere near sufficient.
- **Overlaying without an error metric.** Layer 3's magnetotelluric comparison is a *visual* 50% overlay with no quantitative agreement figure — weaker evidence than layer 2's explicit error function. Note the asymmetry when citing it.
- **Correlating instruments with incompatible time supports.** A 14-second acquisition versus a month of seismicity supports a spatial claim only; the paper says so, and repeating the claim without that caveat overstates it.
- **Dropping the stations that disagree.** Here the exclusion is defensible and stated (visibility and scene coverage), but the principle matters: state the selection rule, not just the survivors.

## Key Takeaways
1. Three layers — geometry (SRTM DEM), measurement (INGV seismographs), structure (magnetotelluric + earthquake catalog) — each against an outside instrument.
2. Layer 2 is the strongest: radar-versus-seismograph spectra and time series with a plotted error oscillating around zero.
3. Only two of the network's stations qualified; the selection rule (radar line-of-sight + inside scene) is stated.
4. Layer 3's overlays are visual, not quantitative — weigh them accordingly.
5. The temporal-support mismatch against the earthquake catalog is acknowledged by the authors and limits that check to spatial/depth correspondence.
6. Validating against a signal neither instrument controls (ambient seismic background) is what makes layer 2 hard to fake.

## Connects To
- **Ch 05, Ch 06**: the Giza paper's validation approach — LiDAR for exterior, overlap/non-overlap rendering for known interior structures — weaker than this three-layer protocol, and worth contrasting.
- **Ch 08**: cross-method consistency checking, of which layer 3 is the concrete instance.
- **Ch 14**: the tomograms these layers validate.
- **Ch 09**: the drilling target proposed at Giza is the physical validation this protocol substitutes for.

# Chapter 13: The 11-Block Computational Scheme & Application Domains

*Source: WO 2024/008365 A1, paragraph [0007]-[0008], Figures 0.5, 0.6, 0.7, 0.8.*

## Core Idea
The patent publishes the complete processing pipeline as 11 numbered, interconnected blocks — explicitly "to ensure the repeatability of the experiments" — and then demonstrates it on two geometries the Giza paper never shows: a volcano interior and a below-ground void beside a pyramid superstructure.

## Frameworks Introduced

- **The 11-block pipeline**: a reproducible, implementation-ready decomposition. Two of the blocks (8 and 10) are explicitly *non-computational* — they are data states, not operations. Naming data states as blocks alongside operations is what makes the diagram runnable.
  - When to use: as the reference implementation checklist for any MM Doppler tomography build.
  - How: see the Worked Example below.

- **Symmetric master/slave branch topology**: blocks 3-5 and 4-6 are mirror-image chains (band-pass filter -> IDFT2) that rejoin at the pixel tracker. The single DFT2 at block 2 is computed once and copied to both branches — blocks 3 and 4 "represent the copy of the DFT2 of 1, and therefore contain the same data, having a common source."
  - Efficiency consequence: one forward FFT2 per image, then two inverse FFT2s per sub-aperture step. Do not recompute the forward transform inside the `N_D` loop.

- **Multi-geometry tomographic lines**: Figure 0.7 shows two satellites (SAR 1, SAR 2) at ranges R1 and R2 producing two independent tomographic lines and two tomographic maps of the same scene — one into a below-ground void, one into an above-ground pyramid.
  - When to use: when one line's orientation cannot reach the feature of interest.
  - Principle: the tomographic line's orientation determines which internal features become visible; multiple acquisition geometries give multiple, differently-oriented lines.

- **Above-ground and below-ground are the same problem**: the method images "below the Earth surface *and* a pyramid-shaped superstructure extending above the Earth's surface" with identical processing. There is no separate above-ground mode.

## Reference Tables — Figure 0.5, the 11 blocks

| # | Block | Computational? | Role |
|---|---|---|---|
| 1 | SLC single SAR image | data | Sole required input |
| 2 | DFT2 (in practice FFT2) | yes | Forward 2-D transform, computed once |
| 3 | Band-pass filter, master | yes | Black-squared Doppler sub-band (Fig. 0.2) |
| 4 | Band-pass filter, slave | yes | Blue-squared sub-band, offset by `B_shift` |
| 5 | IDFT2 (IFFT2) of master | yes | Back to lower-azimuth-resolution SLC |
| 6 | IDFT2 (IFFT2) of slave | yes | Back to lower-azimuth-resolution SLC |
| 7 | Pixel tracking algorithm | yes | Sub-pixel coregistration -> complex vectors |
| 8 | Seismic waves raw data collection | **no** | Raw tomographic complex data (time domain) |
| 9 | Focusing processing by FFT2 | yes | Depth/elevation compression, orthogonal to slant range |
| 10 | Tomographic map | **no** | Focused tomographic image |
| 11 | Tomographic filter and in-depth geolocation | yes | Geocoding into 3-D geographic reference system |

Flow: `1 -> 2 -> {3 -> 5, 4 -> 6} -> 7 -> 8 -> 9 -> 10 -> 11`

## Figure inventory (what each drawing sheet actually shows)

| Figure | Sheet | Content |
|---|---|---|
| 0.1 | 1/7 | SAR acquisition geometry: height, range, azimuth, moving target T1 with `v_r`/`v_a`, swath, orbit, `G_sa`, `d_a`, `d`, `theta` |
| 0.2 | 1/7 | Frequency allocation strategy in 3 positions — the `B_shift` band march |
| 0.3 | 2/7 | Master-slave pixel tracking: box 1 (master), box 2 (slave), displacement `d` at angle `theta` |
| 0.4 | 3/7 (also on cover) | Vibrational estimation geometries (a)-(d): spring-mass surface, vibration traces, single spring `L`/`L0`, orbital sub-aperture scheme |
| 0.5 | 4/7 | The 11-block computational scheme |
| 0.6 | 5/7 | Experimental geometry 1: satellite imaging a **volcano**, acoustic-wave ripple rising through it, tomographic line and map |
| 0.7 | 6/7 | Experimental geometry 2: SAR 1 + SAR 2, one tomographic map into a below-ground void, one into a **pyramid** |
| 0.8 | 7/7 | Experimental results: (a),(c) raw complex vibrational magnitude in time domain; (b),(d) compressed/focused tomograms revealing targets T1 (depth, to -1000) and T2 (height, to ~500) |

## Worked Example — running the pipeline end to end

Given one SLC image and a chosen tomographic line:

1. **Block 1** — load the single SLC SAR image. No pair, no repeat pass.
2. **Block 2** — take FFT2 once. Result is the rectangular spectrum of ch11 equation (3).
3. **Blocks 3 & 4** — copy that spectrum to two branches. Band-pass the master over `B_cr x (B_cD - B_DL)`; band-pass the slave over the same width offset by `B_shift`. `B_shift` encodes the mechanical frequency you want to observe.
4. **Blocks 5 & 6** — IFFT2 each branch back to image domain. Both are now lower-azimuth-resolution SLC images of the same scene at two instants of the aperture.
5. **Block 7** — sub-pixel coregister master against slave for every pixel on the tomographic line. Output: complex displacement vectors `{a, b}`.
6. **Repeat 3-7** `N_D` times, sliding the rigid master-slave pair by `(B_cD - B_DL)/N_D` each step. Do **not** re-run block 2.
7. **Block 8** — assemble the `N_D` samples per pixel into `Y`, the raw tomographic complex data. This is the "seismic waves raw data collection" — sound in the time domain, not yet an image. Figure 0.8 (a) and (c) show what it looks like: dense oscillatory traces with no visible structure.
8. **Block 9** — FFT2 focusing along elevation/depth, orthogonal to the slant range. This is `h(z) = A^dagger * Y` from ch12.
9. **Block 10** — the tomogram. Figure 0.8 (b) and (d) show the same data after compression: a clear basin down to -1000 with target T1 marked, and a peaked above-ground profile up to ~500 with T2 marked. The structure was invisible in (a)/(c) and obvious in (b)/(d) — pulse compression is what turns noise-looking traces into a tomogram.
10. **Block 11** — apply the tomographic filter and geolocate in depth into a 3-D geographic coordinate system; release as GIS-standard data (claim 10).

## Application domains claimed and illustrated

- **Volcanoes** — Figure 0.6; the acoustic-wave ripple through the edifice is the illumination. This is the domain of the arXiv preprint the ISR cites against the claims (see ch10).
- **Underground voids** — Figure 0.7, tomographic map 1.
- **Pyramids / above-ground superstructures** — Figure 0.7, tomographic map 2. This is the Giza application.
- **Undersea and under-ice** — claimed in claim 1 and the abstract, not separately illustrated.
- **Man-made structures** — buildings, bridges, dams; the continuity with Persistent Scatterer Interferometry structural-health monitoring (ch01).
- Depth reach claimed: "several kilometres from the Earth's surface"; Figure 0.4 (a)/(b) annotate the depth axis to -3 km.

## Anti-patterns

- **Recomputing the forward FFT2 inside the sub-aperture loop.** Blocks 3 and 4 are copies of block 2's output by design; treating them as independent transforms multiplies cost by `N_D` for no gain.
- **Reading Figure 0.8 (a)/(c) as a failed result.** Raw vibrational magnitude is *supposed* to look like structureless oscillation; judging the method before block 9 is judging it before the measurement exists.
- **Expecting one tomographic line to reveal everything.** Line orientation gates visibility — this is why Figure 0.7 uses two satellites and two lines for one scene.

## Key Takeaways

1. The pipeline is 11 blocks, published for reproducibility; blocks 8 and 10 are data states, not operations.
2. One forward FFT2 per image; two IFFT2 per sub-aperture step; `N_D` steps total.
3. Block 7 (pixel tracking) is the measurement; block 9 (FFT2 focusing) is what makes it legible.
4. Block 11 geocodes into 3-D and emits GIS-standard output — the method is designed to end in a GIS, not in a figure.
5. Volcano imaging, not pyramid imaging, is the patent's lead illustrated application.
6. Above-ground and below-ground targets are processed identically.
7. Multiple satellites/geometries give multiple tomographic lines and therefore multiple views of one structure.

## Connects To

- **Ch 3, Ch 4**: the paper's method sections, which describe this pipeline without publishing the block diagram.
- **Ch 6, Ch 7**: the Giza internal-imaging results this pipeline produced.
- **Ch 11, Ch 12**: blocks 3-7 implement ch11; blocks 8-10 implement ch12.
- **Ch 10**: claims 7-10 map one-to-one onto pipeline stages.

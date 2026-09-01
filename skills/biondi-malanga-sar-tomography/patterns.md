# Patterns & Techniques

## Doppler Sub-Aperture Vibration Extraction
**When to use**: You need a time series of a target's micro-motion from a *single* SAR SLC acquisition (no repeat pass available or needed).
**How**: Split the focused SAR signal's Doppler bandwidth into N_D sub-apertures (withholding a guard band B_DL for sensitivity); sub-pixel coregister master/slave sub-aperture realizations; track pixel displacement across sub-apertures to build a time-domain vibration trace.
**Trade-offs**: Avoids needing a second acquisition (unlike InSAR), but each sub-aperture has lower azimuth resolution than the full image, and full-resolution tomogram processing is computationally heavy (~days per tomogram in the paper's example).

## Harmonic-Oscillator Tomographic Inversion
**When to use**: Converting per-pixel vibration measurements along a chosen depth line into a physically-interpretable depth profile ("tomogram").
**How**: Model each pixel as a damped 2-DOF linear harmonic oscillator; assemble the multi-baseline data vector Y and steering matrix A(Kz,z); approximate the inverse via matched filtering, h(z) ≈ Aᴴ(Kz,z)·Y.
**Trade-offs**: Matched-filter inversion is fast but is an approximation, not an exact/regularized solve — sensitive to model assumptions (linear oscillator, single dominant scatterer per resolution cell).

## Interferometric Fringe-Inclination Surface Mapping
**When to use**: Detecting subtle shape deviations (bowing, hidden facets) on a large static structure's exterior using two SAR passes.
**How**: Acquire a repeat-pass pair with a suitable spatial baseline; compute interferometric phase and coherence; read changes in fringe spacing/inclination along a face as slope-change boundaries.
**Trade-offs**: Needs two acquisitions (unlike the single-image MM method) and radar-shadow-free geometry; validated here against independent LiDAR to ~0.1–0.35 m accuracy.

## Overlap / Non-Overlap Dual Rendering for Validation
**When to use**: Presenting a new imaging modality's output against a known ground-truth diagram.
**How**: Render the new measurement twice — once overlaid on the accepted schematic, once standalone ("blind") — so agreement can be judged without the schematic biasing the reading of the raw result.
**Trade-offs**: Purely a presentation/validation technique, not a processing method; doesn't itself increase accuracy, only reviewer confidence.

## Unique-Tag Cataloging for Exploratory Imaging
**When to use**: Any project generating many candidate features from a new imaging technique (structures, defects, anomalies).
**How**: Assign a persistent unique ID to every candidate feature at first detection; carry that ID consistently through all figures, tables, and discussion so any feature can be traced end-to-end.
**Trade-offs**: Requires discipline to keep tags stable across drafts/revisions; pays off heavily once cross-referencing prior literature or flagging false alarms (as the paper does for the tag-19 area).

## Cross-Method Consistency Check (Including Non-Matches)
**When to use**: Interpreting new results in light of prior independent methods addressing the same object.
**How**: For each new finding, explicitly search prior literature (different modality: thermal, gravimetric, muon, etc.) for corroboration or conflict; report both matches and non-matches rather than curating only agreements.
**Trade-offs**: Slower and can surface inconvenient disagreements (as with the muon "big void" non-detection here), but substantially strengthens the credibility of the overall interpretation.

## B_shift Frequency Selection
**When to use**: You need the tomogram to be sensitive to a specific mechanical frequency rather than whatever the data happens to contain.
**How**: Focus master and slave from the same SLC image, holding them rigidly separated by B_shift along the azimuth-frequency axis. B_shift *is* the vibrational frequency selector: raise it to observe a lower mechanical frequency, lower it to observe a higher one. Slide the rigid pair across the withheld band in N_D steps to build the time series.
**Trade-offs**: Gives explicit control the Giza paper never exposes, but couples frequency selection to the amount of guard band available; the wider the separation you need, the less room remains for the N_D march.

## Controlled Linearization of a Nonlinear Oscillator
**When to use**: Converting a physically nonlinear vibration model into something a matched-filter inversion can actually solve.
**How**: Start from the true cubic-restoring-force (Duffing) equation r_ddot + lambda*r_dot + omega^2*(1 + xi*r^2)*r = f(omega*t). Check the nonlinearity regime: it dominates when L is close to L0. Only where nonlinearity is sufficiently low, reduce to the 2-DOF damped linear oscillator r(t) = (a*cos(omega0*t), b*sin(omega0*t))*exp(-lambda*t/2), whose {a,b} are the coregistrator's measured shifts.
**Trade-offs**: The linear form is what makes the inversion tractable, but the validity condition is a real boundary, not a formality — near L ~= L0 the cubic term is not negligible and the error propagates silently into the tomogram.

## Publish the Block Diagram for Reproducibility
**When to use**: Disclosing a multi-stage signal-processing method you want others to be able to re-run.
**How**: Decompose the pipeline into numbered, interconnected blocks and state explicitly which blocks are computational operations and which are merely data states. WO 2024/008365 uses 11 blocks and flags 8 (raw seismic data) and 10 (tomographic map) as non-computational.
**Trade-offs**: Costs you disclosure detail a competitor can implement from, but converts an unreproducible narrative into an implementable spec; also maps cleanly onto stage-by-stage patent claims (claims 7-10 here).

## Multi-Geometry Tomographic Lines
**When to use**: One tomographic line's orientation cannot reach the internal feature you care about.
**How**: Acquire the same scene from more than one satellite/geometry (different slant ranges R1, R2, different squint) and run an independent tomographic line per acquisition. Each line yields its own tomographic map of the same object.
**Trade-offs**: Multiplies acquisition and compute cost linearly, but line orientation gates visibility — no amount of processing recovers a feature the chosen line never crosses.

## Sequence Priority Filing Before Self-Disclosure
**When to use**: Any time you intend to both publish research and patent it.
**How**: File the priority application first, then post the preprint. The EPO/PCT system treats the applicant's own earlier publication as prior art with no general grace period.
**Trade-offs**: Delays public disclosure by the filing lead time. Skipping it is not recoverable: WO 2024/008365 (priority 04.07.2022) drew a category-X search opinion against all 10 claims over the inventor's own arXiv preprint of 18.06.2022 — a two-week gap that covers the entire claim set.

## Depth-vs-Resolution Tuning via Investigation Frequency
**When to use**: At acquisition-planning time, before any processing, once you know the depth you must reach.
**How**: delta_z = lambda*R/(2*A) with lambda = v/f. Pick the investigation frequency f from the depth requirement: low f gives a long acoustic wavelength, coarse cells, deep reach; high f gives fine cells, shallow reach. Worked settings across the corpus: 200 Hz -> ~36 m cells, ~3 km depth (Vesuvius); 12.5 kHz -> ~0.92 m cells, pyramid interior (Giza); ~22 kHz -> ~1.30 m cells (patent).
**Trade-offs**: The medium's propagation speed v is not yours to choose and shifts lambda by up to a factor of 6 between a volcanic edifice (~972 m/s assumed) and pyramid stone (~6000-6600 m/s assumed). Getting v wrong scales the entire depth axis.

## Three-Layer Validation of a Novel Imaging Modality
**When to use**: Publishing results from an instrument or processing chain nobody else has run.
**How**: Validate geometry, measurement, and structure separately against three outside references. (1) Geometry: overlay an independent DEM profile along the exact tomographic line onto the tomogram's surface. (2) Measurement: compare the raw quantity you claim to measure against a ground sensor measuring the same quantity, in both spectrum and time domain, and plot the error. (3) Structure: overlay an independent modality imaging the same volume, plus an independent event catalog.
**Trade-offs**: Layer 1 is cheap and only proves geocoding; layer 2 is the expensive, convincing one and requires a co-located ground instrument; layer 3 is often only visual and therefore weakest. Report which layers are quantitative and which are visual rather than presenting them as equal.

## Confirm a Feature by Changing a Processing Parameter
**When to use**: Deciding whether a bright or dark patch in a tomogram is a structure or an artifact.
**How**: Recompute the same scene along a differently-oriented tomographic line, and separately at a different investigation frequency. A feature that survives both changes is unlikely to be a processing artifact. Line orientation is a processing parameter, not an acquisition one, so this costs no new data.
**Trade-offs**: Multiplies compute time per confirmed feature; cannot rescue a feature that lies outside every reachable line orientation.

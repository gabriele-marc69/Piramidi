# Chapter 11: Doppler Sub-Aperture Formalism (Patent Mathematics)

*Source: WO 2024/008365 A1, paragraphs [0004]–[0005], equations (1)–(12), Figures 0.1–0.3.*

## Core Idea
The patent gives the full signal-model derivation the Giza paper compresses into a paragraph: the focused SLC spectrum is rectangular, so it can be sliced into master/slave sub-bands held a fixed distance `B_shift` apart, and marching that pair across the Doppler axis turns one static image into a time series of the target's micro-motion.

## Frameworks Introduced

- **Master/slave band pair with fixed separation `B_shift`**: the master is focused on the full chirp band `B_cT` and half the Doppler band; the slave is the same image focused on a band offset in azimuth frequency by `B_shift`.
  - When to use: whenever you want to select *which* mechanical frequency the tomogram is sensitive to.
  - How: set `B_shift` to the vibrational frequency you wish to observe. **The higher `B_shift`, the lower the observed mechanical frequency** — this is the tuning knob, and it is stated nowhere in the paper.
  - Why it works: the two bands sample the same scene at two different points of the synthetic aperture; their sub-pixel disparity is the displacement accumulated between those two instants.

- **Guard band `B_DL = B_cD / 2`**: half the Doppler bandwidth is deliberately *not* processed by the matched filter.
  - When to use: always, in this method.
  - How: focus master and slave using range-azimuth bandwidth `B_cr, B_cD − B_DL`; the withheld `B_DL` is the room in which the pair can slide.
  - Failure mode if skipped: focusing on the entire Doppler band restores maximum azimuth resolution but destroys motion sensitivity — there is no unused band left to shift through.

- **`N_D` as the mechanical sampling rate**: the withheld bandwidth is divided into `N_D` equally-distributed steps; `N_D` rigid shifts of the master-slave system populate the data matrices.
  - Key reframing: `N_D` is not a processing convenience — it *is* the digital sampling rate at which the Earth's mechanical wave is being sampled. Choose it by Nyquist against the vibration you want, not by compute budget.
  - Each Doppler frequency step equals `(B_cD − B_DL) / N_D`.

- **Motion-artifact taxonomy**: a moving target betrays itself in the focused image through three distinct signatures — azimuth displacement (constant range velocity), azimuth smearing (azimuth velocity or range acceleration), and range-walking / range defocusing (range speed spreading energy across range cells.)
  - Use as a diagnostic checklist: which artifact you see tells you which motion component dominates.

## Reference Tables

| Symbol | Meaning |
|---|---|
| `r` | zero-Doppler distance (constant) |
| `R`, `R₀` | slant-range; reference range at t = 0 |
| `d_a` | physical antenna aperture length |
| `V` | platform velocity |
| `d` | distance between two range acquisitions |
| `G_sa`, `L_sa = 2Nd` | total synthetic aperture length |
| `t`, `T` | acquisition time variable; observation duration |
| `L = λr/d_a` | azimuth electromagnetic footprint width |
| `θ` | incidence angle of the radiation pattern |
| `B_cr` / `B_cT` | total chirp (range) bandwidth |
| `B_cD = 4Nd/λr` | total Doppler bandwidth |
| `B_DL = B_cD/2` | withheld (unprocessed) Doppler band |
| `B_shift` | master–slave azimuth-frequency separation = selected vibrational frequency |
| `N_D` | number of Doppler sub-aperture refocused images = mechanical sampling rate |
| `δ_D ≈ 1/B_cD = λR/2L_sa` | azimuth resolution |

## Code Examples — the equations in sequence

Focused SLC signal of a stationary point target (1):
```
s_SLC(k,x) = 2Nτ · exp[−j(4π/λ)r] · sinc[πB_cr(k − 2R/c)] · sinc[πB_cD·x]
             for x = kt, k ∈ {0…N−1}, x ∈ {0…M−1}
```

Same target off beam-centre, at slant-range/azimuth position `L_cg, L_Dh` (2):
```
s_SLC(k,x) = 2Nτ · exp[−j4πr/λ] · sinc[πB_cr(k − L_cg)] · sinc[πB_cD(x − L_Dh)]
```

Its 2-D DFT collapses to a **rectangular spectrum** with linear phase ramps (3):
```
S_SLC_F(n,q) = 2Nτ·exp[−j4πr/λ] · (1/πB_cr)·rect[n/πB_cr]
                                  · (1/πB_cD)·rect[q/πB_cD]
                                  · exp(−j2πn·L_cg) · exp(−j2πq·L_Dh)
```
- **What it demonstrates**: because the spectrum is rectangular, sub-banding is a clean windowing operation — this rectangularity is what licenses the whole sub-aperture strategy.

Master and slave sub-aperture matrices (4)–(5), one row of `N_D` entries each:
```
S_SLC(k,x)_M = [ S_M{1,1}  S_M{1,2}  S_M{1,3} … S_M{1,N_D} ]
S_SLC(k,x)_S = [ S_S{1,1}  S_S{1,2}  S_S{1,3} … S_S{1,N_D} ]
```

Range history of a target moving with velocity/acceleration `{v_r, v_a}, {a_r, a_a}` (7), Taylor-expanded (9)–(11) to:
```
|R(Vt)| = R₀ − v_r t + (t²/2R₀)·[ (V − v_a)² − R₀ a_r ]
```
and recast in `x = Vt` (12):
```
|R(x)| = R₀ − ε_r1·x + [ (1 − ε_c1)² − ε_r2 ]·x²/2R₀
   ε_r1 = v_r/V      (range velocity)
   ε_r2 = a_r R₀/V²  (range acceleration)
   ε_c1 = v_c/V      (azimuth velocity)
```
- **What it demonstrates**: the three ε terms are the complete first-order description of how motion perturbs the focused signal — they are what the coregistrator actually estimates, and they feed the oscillator model of ch12.

## Worked Example — the three-position band march (Figure 0.2)

Figure 0.2 shows the same rectangular spectrum three times, labelled (1), (2), (3). In each panel a black square (master) and a blue/green square (slave) sit side by side, rigidly separated by `B_shift`, inside the azimuth-frequency axis that runs 0 → `B_cD`.

1. **Position (1)** — the pair starts at the low-frequency edge. Master is focused over `B_cr` × `(B_cD − B_DL)`; slave is the same scene focused `B_shift` further along the azimuth-frequency axis. Coregistering the two gives one sub-pixel displacement sample, `{a, b}` — one instant of the vibration.
2. **Position (2)** — both squares slide by `(B_cD − B_DL)/N_D`. Coregister again → second sample.
3. **Position (3)** — repeat until `N_D` shifts have been made and the withheld band is exhausted.

The ordered set of `N_D` displacement estimates *is* the time-domain vibration trace of that pixel (Figure 0.3 shows the master–slave pixel tracking geometry: box 1 = master pixel, box 2 = slave pixel, displaced by `d` at angle `θ`). Repeat over every pixel of the tomographic line, and you have the raw tomographic data cube that ch12 focuses in depth.

Note what makes this work: the orbital shift along the aperture provides the *time* axis. The satellite's own motion is the clock.

## Key Takeaways
1. `B_shift` selects the observed mechanical frequency — higher `B_shift` → lower observed frequency. This is the single most actionable parameter the paper omits.
2. `B_DL = B_cD/2` is not a lossy compromise but the precondition for motion sensitivity; half the azimuth resolution is the price of the depth dimension.
3. `N_D` is the mechanical sampling rate — pick it against the vibration bandwidth you need, then accept the compute cost.
4. The rectangular spectrum of a focused SAR image is what makes sub-banding mathematically clean.
5. `{ε_r1, ε_r2, ε_c1}` — range velocity, range acceleration, azimuth velocity — are the estimated quantities that carry the vibration.
6. Three motion artifacts (azimuth displacement, azimuth smearing, range-walking) diagnose which motion component dominates a given target.

## Connects To
- **Ch 3**: the paper's shorter account of the same sub-aperture decomposition.
- **Ch 12**: consumes `{ε_r1, ε_r2, ε_c1}` as input to the oscillator and inversion.
- **Ch 13**: block 3/4 (band-pass filters) and block 7 (pixel tracking) of the computational scheme implement this chapter.

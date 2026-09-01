# Chapter 12: Nonlinear Spring Model & Tomographic Inversion (Patent Mathematics)

*Source: WO 2024/008365 A1, paragraphs [0006]-[0007], equations (13)-(24), Figure 0.4.*

## Core Idea
The patent derives the harmonic oscillator that the Giza paper simply asserts: it starts from a genuinely **nonlinear** spring (cubic restoring force, a Duffing-type equation), then states the explicit condition under which it collapses to the tractable 2-DOF linear oscillator whose `{a, b}` shifts the coregistrator measures.

## Frameworks Introduced

- **Spring-lattice surface model**: the Earth's surface is modelled as a field of springs held by masses, oscillating perpetually, each contributing its own oscillatory component to one image pixel. A tomographic line is a strip of adjacent such pixels.
  - When to use: reasoning about what a single pixel's vibration trace physically represents.
  - Key detail: `L` is the spring length at maximum tension, `L0` its length with no mass, `xi` its elastic constant.

- **Standing-wave energy premise**: an impulse perturbs the spring; the wave travels, reflects off a constraint, returns, and interferes with itself. Two counter-propagating sine waves of equal amplitude and frequency produce an **ideal perpetual standing wave** — which is why the surface can be treated as continuously "ringing" rather than needing an active acoustic source.
  - Corollary the patent states explicitly: longitudinal (tension-domain) oscillations propagate at roughly **twice** the frequency of the transverse ones, and transverse-longitudinal coupling is essentially nonlinear.

- **Controlled linearization**: nonlinearity dominates when `L ~= L0`. The reduction to a linear oscillator is licensed only when nonlinearity is "sufficiently low" — i.e. when the spring is meaningfully tensioned. State this condition when you apply the model; it is the model's real boundary.

- **Steering-matrix DFT equivalence**: the inversion is not a generic regularized solve. `A(Kz, z)` is explicitly "the best approximation of a matrix operator performing the DFT of Y", so the tomogram is obtained *by pulse compression*. Depth focusing and Fourier focusing are the same operation here.

## Code Examples — the derivation chain

Vibrational force on mass `m1` (13):

```
F = -4*xi*r*( 1 - L0 / sqrt(L^2 + 4*r^2) )
```

For `r << L`, series-expand (14) and keep the cubic restoring term (15):

```
F = m*r_ddot ~= -4*xi*r*(L - L0)*(r/L)*[ 1 + (2*L0/(L - L0))*(r/L)^2 ]
```

Define natural frequency (16) and nonlinearity coefficient (17):

```
omega0 = (4*xi/m)*[ (L - L0)/L ]
xi     = (2*L0/L^2)*(L - L0)
```

Undamped nonlinear equation (18), then with damping `lambda` and forcing `f(omega*t)` (19) — a Duffing-form oscillator:

```
r_ddot + omega^2*r*(1 + xi*r^2) = 0
r_ddot + lambda*r_dot + omega^2*(1 + xi*r^2)*r = f(omega*t)
```

If nonlinearity is sufficiently low, reduce to the 2-DOF damped linear oscillator (20):

```
r(t) = ( a*cos(omega0*t), b*sin(omega0*t) ) * exp(-lambda*t/2)
```

- **What it demonstrates**: `{a, b}` are exactly *the instantaneous shifts estimated by the coregistrator* — the bridge from ch11's pixel tracking to physics. They approximate the `{eps_r1, eps_r2, eps_c1}` displacement parameters of equation (12).

Multi-frequency data vector from `k` time samples (21), steering matrix (22):

```
Y = [ y(1), ..., y(k) ]  in C^(k x 1)

A(Kz, z) = [ 1, exp(j*2*pi*kz2*t*z0),     ..., exp(j*2*pi*kz(k-1)*t*z0)     ]
           [ 1, exp(j*2*pi*kz2*t*z1),     ..., exp(j*2*pi*kz(k-1)*t*z1)     ]
           [                          ...                                    ]
           [ 1, exp(j*2*pi*kz2*t*z(F-1)), ..., exp(j*2*pi*kz(k-1)*t*z(F-1)) ]

with  Kz = 4*pi*B_perp / (lambda * r_i * sin(theta)),   i = 1...k
```

where `B_perp` is the i-th orthogonal baseline and `r_i` the i-th slant-range distance.

Sonic tomographic model (23) and its solution (24):

```
Y    = A(Kz, z) * h(z)
h(z) = A(Kz, z)^dagger * Y        (dagger = matched-filter / pseudo-inverse)
```

Tomographic resolution:

```
delta_T = lambda*R / (2*A)
```

- **Critical clarification the paper leaves implicit**: here `lambda` is the **sound wavelength over the Earth**, not the radar wavelength. `R` is slant range and `A` is the orbit aperture used in the tomographic synthesis — and `A` is proportional to the Doppler bandwidth used to synthesize the sub-apertures.

## Reference Tables — patent worked parameters vs. the Giza paper

| Quantity | Patent (WO 2024/008365) | Giza paper (RS 2022, 14, 5231) |
|---|---|---|
| Seismic propagation speed `v` | ~6600 m/s | ~6000 m/s |
| Max investigation frequency `f` | ~22 000 Hz | 12 500 Hz |
| Sound wavelength `lambda = v/f` | ~0.30 m | not stated |
| Slant range `R` | 650 000 m | 650 000 m |
| Orbit aperture `A` | 75 000 m (half total orbit length) | not stated |
| SLC Doppler synthesis | 22 kHz | 22 kHz |
| **Tomographic resolution `delta_z`** | **~1.30 m** | **~0.92 m** |

Read the difference as a range, not a contradiction: `delta_z = lambda*R/(2*A)` is fully determined by the assumed sound speed, the chosen investigation frequency, and how much orbit you are willing to synthesize. The patent takes the aggressive-aperture / high-frequency case; the paper takes a lower investigation frequency.

## Worked Example — computing `delta_z` the patent's way

The patent walks the arithmetic explicitly, and it is worth reproducing because it shows which knobs actually move the resolution:

1. Assume average seismic propagation speed `v ~= 6600 m/s`.
2. Take the maximum observable investigation frequency `f ~= 22 000 Hz` (matching the 22 kHz Doppler synthesis of the SLC data).
3. Sound wavelength: `lambda = v/f = 6600 / 22000 ~= 0.30 m`.
4. Extend the tomography to the maximum orbital aperture — half the total orbit length — giving `A ~= 75 000 m`.
5. With `R = 650 000 m`:
   `delta_z = lambda*R / (2*A) = (0.30 * 650 000) / (2 * 75 000) = 195 000 / 150 000 ~= 1.30 m`

"This is the tomographic resolution set to calculate all the results estimated by the invention."

What the arithmetic tells you: resolution improves linearly with **higher investigation frequency** (shorter sound wavelength) and with **longer synthesized aperture**. It degrades linearly with slant range. Sound speed is the one term you do not control — it is a property of the medium, and getting it wrong scales your depth axis directly.

## Anti-patterns

- **Treating `lambda` in `delta_z = lambda*R/(2*A)` as the radar wavelength.** It is the acoustic wavelength in the medium. Using the 9.6 GHz radar wavelength (~3 cm) gives a resolution roughly 10x too optimistic.
- **Applying the linear oscillator (20) when `L ~= L0`.** That is precisely the regime where the patent says nonlinearity dominates; the cubic term is not negligible there and the inversion inherits the error silently.
- **Assuming a single dominant scatterer per resolution cell.** Equation (23) is a single-source-per-elevation model; the matched-filter inverse (24) is fast but unregularized.

## Key Takeaways

1. The oscillator is nonlinear (Duffing-type, cubic restoring force) and only *reduces* to the linear 2-DOF form when tension keeps `L` well away from `L0`.
2. `{a, b}` in equation (20) are the coregistrator's measured shifts — the measurement-to-physics bridge.
3. `Kz = 4*pi*B_perp/(lambda*r_i*sin(theta))` carries the orthogonal baseline into the depth phase ramp.
4. `h(z) = A^dagger * Y` is pulse compression — depth focusing *is* Fourier focusing here.
5. `delta_z = lambda*R/(2*A)` with lambda = **sound** wavelength; patent case ~1.30 m, paper case ~0.92 m.
6. Longitudinal oscillations run at ~2x the transverse frequency; their coupling is nonlinear.

## Connects To

- **Ch 4**: the paper's version of the harmonic-oscillator tomographic model and resolution formula.
- **Ch 11**: supplies `{eps_r1, eps_r2, eps_c1}` and the `N_D` samples that populate `Y`.
- **Ch 13**: blocks 8-10 (raw tomographic data -> FFT2 focusing -> tomographic map) implement equations (21)-(24).

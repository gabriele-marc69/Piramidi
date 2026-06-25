# Skills — SAR micro-Doppler tomography of the Giza pyramids

Skills (anthropics/skills format: `SKILL.md` with `name`/`description` frontmatter + bundled
`scripts/`) distilled from the validated pipeline that reproduces the Biondi & Malanga
(Remote Sens. 2022, 14, 5231; arXiv:2206.09200) figures from a real Sentinel-1 SLC stack.

Pipeline order:

1. **sentinel1-slc-reader** — read SLC TIFFs + annotations, geolocate a lat/lon box via
   GCPs, extract a co-registered complex chip stack `(n_scene, az, rng)` → `.npz`.
2. **sar-doppler-tomography** — Doppler sub-apertures + micro-Doppler tomographic inversion
   `h(z)=Aᴴ·Y` → reflectivity volume + B-scan sections, slices, height map (the
   publication tomographic images). Includes a synthetic forward model that verifies the
   chain recovers reflectors injected at known depths.
3. **sar-dinsar-microdisplacement** — interferometric-phase LOS micro-displacements (mm),
   waveforms and 3-D cloud (the InSAR figures).

All three honestly document the fundamental limit: an open Sentinel-1 IW2 stack over a
small box gives few tomographic looks → low ambiguity height (`z_amb ≈ 8.5 m` for the
Khafre box), so absolute pyramid-scale depth is not recoverable — the image *types* are
reproduced, not the paper's deep high-resolution interior.

Target data: `Piramid/stack_slc/` (6 VV scenes, IW2, May–Jun 2026, Giza plateau).

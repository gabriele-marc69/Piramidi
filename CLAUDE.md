# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Python pipeline for SAR micro-Doppler tomography and DInSAR micro-displacement of the Giza pyramids (Cheope/Khufu and Kefren/Khafre) from open Sentinel-1 IW SLC stacks, inspired by Biondi & Malanga (Remote Sens. 2022, 14, 5231; arXiv:2206.09200). Code comments, docstrings, CLI help, and console output are in **Italian** — keep new code consistent with that.

## Commands

No test suite, linter, or build system is configured. Development = running the pipeline scripts directly.

```bash
pip install -r requirements.txt   # numpy, scipy, matplotlib, rasterio, plotly (Python >= 3.10)

# Full pipeline (search + download + processing). --pyramid picks the default
# coordinate box and output dir (goal_out_Cheope/ or goal_out_kefren/):
python piramide_unificato.py --pyramid cheope --download
python piramide_unificato.py --pyramid kefren --download

# Processing only, on data already in the output dir's stack_slc/ (no credentials needed):
python piramide_unificato.py --pyramid kefren --steps 3-6

# Paper-faithful Doppler tomography (steering matrix, h(z)=A^H*Y) in addition to chosen steps:
python piramide_unificato.py --pyramid kefren --steps 3 --vera-tomografia

# Windows wrapper (defaults: kefren, --steps 4-6, prompts for/caches CDSE credentials in .cdse.env):
piramide_unificato.bat [cheope|kefren] [extra args passed through]
```

Standalone skill scripts (each runnable on its own):

```bash
python skills/sentinel1-slc-reader/scripts/extract_box.py --stack <stack_slc> --nw D M S H D M S H --se D M S H D M S H --pol vv --out box.npz
python skills/sar-doppler-tomography/scripts/tomographic_images.py --stack <stack_slc> --outdir tomo_img --zmax 200 --klook 12 --ovs 4
python skills/sar-dinsar-microdisplacement/scripts/dinsar.py --npz box.npz
```

Download credentials (step 2 only) come from env vars `CDSE_USER` / `CDSE_PASS` (free account at dataspace.copernicus.eu) or `--cdse-user`/`--cdse-pass`. Never hardcode or commit them; `.cdse.env` is gitignored.

## Architecture

- **`piramide_unificato.py`** — the single entry point (~2900 lines), a 6-step pipeline selected with `--steps` (e.g. `1-6`, `3,5`):
  1. OData search on Copernicus Data Space (CDSE) for Sentinel-1 IW SLC products containing the coordinate box
  2. Resumable download; extracts measurement TIFFs + annotation.xml from `.SAFE` zips (also reuses zips already downloaded by `scarica_alta_risoluzione.py` in `alta_risoluzione_out/`)
  3. Box extraction via product GCPs → co-registered complex stack `(n_scene, az, rng)` → `box.npz`
  4. DInSAR: interferometric phase → vertical displacement layers `(n_strati, ny, nx)`
  5. Layers treated as Fourier coefficients of a depth-stretched waveform per pixel; 3-D Plotly "springs" colored by instantaneous (Hilbert) frequency
  6. Per-pixel PIENO/VUOTO (solid/void) classification at real DEM heights
  Pyramid presets (`PRESETS` dict: box DMS coordinates, output dir, real base/height geometry) live at the top of the file; the known internal structures of Khufu used as visual overlay are hardcoded there too.
- **`skills/`** — three reusable skills (anthropics/skills format, `SKILL.md` + `scripts/`) mirroring the pipeline stages: `sentinel1-slc-reader` (step 3), `sar-dinsar-microdisplacement` (step 4), `sar-doppler-tomography` (the paper's actual method, invoked by `--vera-tomografia`). The `sar_tomo.py` library includes a synthetic forward model (`genera_Y_sintetico`) that validates the tomographic inversion against reflectors at known depths — use it to verify changes to the tomography chain.
- **`scarica_alta_risoluzione.py`** — standalone search/download for the highest-resolution SLC products on CDSE (Stripmap SLC → on-demand production from RAW L0 → IW SLC fallback). `scarica_dati.bat` wraps it.
- **`cerca_cosmo_skymed.py`** — searches public catalogs (Zenodo, OpenAIRE, Brave/Bing web APIs) for pointers to COSMO-SkyMed material. It deliberately does NOT download restricted scenes or bypass access controls/CAPTCHAs — keep it that way.
- **Output dirs** `goal_out_Cheope/` and `goal_out_kefren/` hold per-pyramid data: `box.npz` and `unificato/` (results). The SLC TIFF stack lives in `goal_out_kefren/stack_slc/` and is **shared by both pyramids** (the same Sentinel-1 scenes cover all of Giza): the `cheope` preset sets `stack_outdir_name="goal_out_kefren"` so it searches TIFFs there instead of downloading its own copy (`--stack` still overrides). Raw data, `.npz`, `.png`, `.html` are all gitignored (SAFE products are ~8 GB each).

## Critical methodological distinction

Steps 5–6 of `piramide_unificato.py` are an **exploratory** analysis (DInSAR layers reinterpreted as Fourier coefficients over an arbitrary depth range + Hilbert instantaneous frequency) inspired by the paper's vocabulary but **not** its method. The paper-coherent Doppler tomography (sub-apertures, steering matrix, `h(z)=Aᴴ·Y`) is only in `skills/sar-doppler-tomography/` / `--vera-tomografia`. Physical limit stated throughout the repo: `sar_tomo.GeometriaSAR.apertura_orbitale` (the synthesis aperture `A` in both `δz=λR/2A` and `Kz`/`z_amb`) must be the REAL aperture of the source data, not the paper's own COSMO-SkyMed Spotlight value (84 km) — Sentinel-1 IW only synthesizes `A≈1.1 km` per burst (`sar_tomo.apertura_sintetica_em`, same formula as `_profondita_massima_fourier`). With that real aperture: ambiguity height `z_amb ≈ 660 m` (klook=12) but per-cell resolution `δz ≈ 24 m` (ovs=4) — coarse enough that meter-scale internal chambers are still far below what this stack can resolve, even though depths up to `z_amb` no longer alias. (Prior to a 2026-07 fix, `baseline_ortogonali` hardcoded the paper's aperture regardless of platform, silently giving an unrealistic `z_amb≈8.5 m` / `δz≈0.3 m` for Sentinel-1 data — if you see that stated anywhere, it's stale.) Any docs, figures, or write-ups you produce must preserve these honest disclaimers — do not present the exploratory output as validated depth measurement or as evidence of undiscovered structures.

## Conventions

- Steps 1–2 must run on a machine without scipy/matplotlib/plotly: heavy imports are lazy, inside the steps that use them. Keep it that way when editing.
- Coordinates are passed as DMS quadruples (`--nw 29 58 38.0 N --nw-lon 31 7 45.4 E ...`); geolocation goes through rasterio `GCPTransformer` and reads only the box window, never the whole ~1.5 GB scene.
- The `biondi-malanga-sar-tomography` skill (available via the Skill tool) is the knowledge base for the source paper — use it when touching the tomography math or interpreting results.

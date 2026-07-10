# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Python pipeline for SAR micro-Doppler tomography and DInSAR micro-displacement of the Giza pyramids from Sentinel-1 IW SLC stacks, inspired by Biondi & Malanga (Remote Sens. 2022, 14, 5231). Code comments, CLI help, and output messages are in **Italian** — keep that convention.

There is no build system, test suite, or linter. Verification is done by running the pipeline steps against local data and a synthetic model (see below).

## Commands

```bash
pip install -r requirements.txt   # numpy, scipy, matplotlib, rasterio, plotly (Python >= 3.10)

# Full pipeline (search + download + steps 1-6); needs CDSE_USER/CDSE_PASS env vars
python piramide_unificato.py --pyramid cheope --download
python piramide_unificato.py --pyramid kefren --download

# Processing only, on data already in <outdir>/stack_slc/ (no network)
python piramide_unificato.py --pyramid kefren --steps 3-6

# Paper-faithful Doppler tomography (steering matrix, h(z)=A^H*Y) in addition to chosen steps
python piramide_unificato.py --pyramid kefren --steps 3 --vera-tomografia

# Higher-resolution product search/download (Stripmap SLC > on-demand L0 > IW SLC fallback)
python scarica_alta_risoluzione.py --pyramid cheope --dry-run
scarica_dati.bat [cheope|kefren]   # Windows wrapper; caches credentials in .cdse.env (gitignored)

# Standalone box extraction from already-downloaded TIFFs
python skills/sentinel1-slc-reader/scripts/extract_box.py --stack <stack_slc> --nw ... --se ... --pol vv --out box.npz

# Synthetic self-test of the tomographic inversion (reflectors at known depths)
python skills/sar-doppler-tomography/scripts/sar_tomo.py
```

Custom AOI: pass `--nw/--nw-lon/--se/--se-lon` (DMS, 4 tokens each: D M S H) instead of `--pyramid`. Other key flags: `--pol vv,vh` (first is primary), `--layers`, `--start/--end`, `--dmin/--dmax` (step-6 density band), `--zmax/--klook/--ovs` (tomography), `--no-track-filter`, `--debug`.

## Architecture

Almost everything lives in `piramide_unificato.py` (~1300 lines), a single 6-step pipeline selected via `--steps`; the `PRESETS` dict at the top defines the per-pyramid AOI boxes and real geometry (`base_m`/`h_m` — used as the internal metric reference; wrong values falsify heights):

1. **step1_cerca** — OData search on Copernicus Data Space (CDSE) for SLC products containing the AOI.
2. **step2_scarica** — resumable download (retries, token refresh per attempt, waits for network up to 2h); extracts measurement TIFFs + annotation XMLs from the `.SAFE` zips into `<outdir>/stack_slc/`.
3. **step3_estrai_box** — geolocates the box via product GCPs, extracts a co-registered complex chip stack `(n_scene, az, rng)` per polarization → `box.npz`; filters scenes to the dominant track.
4. **step4_array_12** — classic DInSAR: interferograms between acquisition pairs → N-layer LOS→vertical displacement array `(n_strati, ny, nx)`.
5. **step5_grafico_onde** — treats the N layers as Fourier coefficients of a depth waveform stretched over an arbitrary range; 3D plotly surface + "springs". VH traces are interleaved between VV columns.
6. **step6_variazioni** — Hilbert instantaneous frequency → points of frequency (density) variation, plotted at real heights (Copernicus DEM + pyramid overlay).

**Critical methodological distinction** (documented in the module docstring and README): steps 5–6 are an *exploratory* extension, NOT the paper's Doppler sub-aperture tomography. The paper-faithful method (steering matrix, `h(z)=Aᴴ·Y`, Kz from synthesized orthogonal baselines) lives in `skills/sar-doppler-tomography/scripts/` and is invoked via `--vera-tomografia` (`step_vera_tomografia` calls `tomographic_images.py`). Preserve this distinction in any docs or output text; never present step 5–6 depths as physically validated. Physical limit to respect: few tomographic looks on a small box ⇒ ambiguity height `z_amb ≈ 8.5 m`, so keep `--zmax` below it.

`scarica_alta_risoluzione.py` is a separate search/download tool that prefers higher-resolution products (Stripmap SLC, then on-demand production from L0 RAW, then IW SLC; Academic Torrents as last resort). Its `PIRAMIDI` presets must stay consistent with `PRESETS` in `piramide_unificato.py`.

`skills/` contains three reusable skills (anthropics/skills format, `SKILL.md` + `scripts/`) mirroring pipeline stages: `sentinel1-slc-reader` (box extraction), `sar-doppler-tomography` (paper inversion + synthetic validation), `sar-dinsar-microdisplacement` (LOS micro-displacements). A fourth knowledge skill, `biondi-malanga-sar-tomography`, is available in the Claude Code session for the paper's method/notation.

## Conventions and constraints

- Only SLC products are usable (they keep phase); GRD is amplitude-only and useless here. Each `.SAFE` is ~8 GB.
- Raw data (`.SAFE`, `.tiff`, `.zip`), intermediate `.npz`, and renders (`.png`, `.html`) are gitignored; outputs go to `goal_out_Cheope/` or `goal_out_kefren/` (subfolder `unificato/`), and `alta_risoluzione_out/` for the high-res downloader.
- CDSE credentials come from `CDSE_USER`/`CDSE_PASS` env vars or `--cdse-user/--cdse-pass`; never hardcode them (a past commit removed hardcoded ones). `.cdse.env` is local-only.
- matplotlib/scipy/plotly are imported lazily inside the steps that need them so steps 1–2 run on a machine with only numpy.
- The reference PDFs in the repo root (arXiv 2206.09200, 2208.00811, 2007.05326) are the source papers for the method.

---
name: sentinel1-slc-reader
description: Read Sentinel-1 IW SLC measurement GeoTIFFs (complex data with phase), geolocate a ground box via the product GCPs, and extract a co-registered complex chip stack from a multi-temporal stack. Use when working with .tiff SLC files (e.g. s1a-iw2-slc-vv-*.tiff) plus their *.annotation.xml, when you need the complex (amplitude+phase) values inside a lat/lon box, or before any InSAR / Doppler-tomography processing.
---

# Sentinel-1 SLC reader & box extractor

Sentinel-1 IW SLC products are the only open Sentinel-1 data that keep the **phase**
(GRD is amplitude-detected and useless for tomography/InSAR). Each subswath/polarization
is one big complex GeoTIFF (`*-iw2-slc-vv-*.tiff`, ~1.5 GB) paired with an
`*.annotation.xml` carrying geometry and the geolocation grid.

## Workflow

1. **Pair each TIFF with its annotation.** Match on sensor (`s1[abcd]`), subswath
   (`iw[1-3]`) and date (`YYYYMMDD`) parsed from the file name.
2. **Read geometry from the XML**: `radarFrequency` (carrier), `rangePixelSpacing` (dR),
   `azimuthPixelSpacing` (dA), `azimuthTimeInterval` (dt_az), `incidenceAngleMidSwath`,
   `slantRangeTime` (→ near slant range `R_near = c·slantRangeTime/2`), and orbit
   `velocity` (platform speed V). These feed the SAR geometry of downstream steps.
3. **Geolocate the box.** Convert the box corners (DMS → decimal lon/lat) and map them
   to pixel (row=azimuth, col=range) with rasterio's `GCPTransformer` over `ds.gcps[0]`.
   Take the bounding rectangle, clamp to the raster.
4. **Windowed complex read.** Read only the box window — never the whole 1.5 GB scene.
   SLC TIFFs come as `complex_int16`/`complex64` or a 2-field structured dtype (I,Q);
   always normalize to `complex64` (see `leggi_complesso`).
5. **Stack & co-register.** Crop every scene to the common minimum `(H, W)` and stack to
   `(n_scene, azimuth, range)` `complex64`. For sub-pixel co-registration of the stack,
   use the Doppler sub-aperture cross-correlation in the `sar-doppler-tomography` skill.

## Script

`scripts/extract_box.py` does all of the above end-to-end and saves an `.npz`
(`vv` complex array + geometry scalars) plus a per-pixel CSV.

```bash
python scripts/extract_box.py \
  --stack /path/to/stack_slc \
  --nw 29 58 38.0 N 31 7 45.4 E \
  --se 29 58 29.0 N 31 7 55.4 E \
  --pol vv --out box.npz
```

Defaults target the Khafre (Chefren) pyramid box on the Giza plateau. Output array
shape is `(n_scene, H, W)`; `H` is azimuth samples (few → few tomographic looks),
`W` is range samples.

## Gotchas
- IW2 over Giza: the box maps to ~26 (azimuth) × ~81 (range) native pixels — slant range
  step 2.33 m (→ ~3.66 m ground at θ≈39°), azimuth step 13.95 m.
- A scene whose corners fall outside the raster does not cover the box → skip it.
- `slantRangeTime` and the column index give the per-pixel slant range `R0` used by the
  tomographic steering matrix.

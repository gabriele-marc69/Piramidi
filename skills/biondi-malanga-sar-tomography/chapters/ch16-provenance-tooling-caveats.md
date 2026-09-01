# Chapter 16: Provenance, Tooling & Preprint Caveats

*Source: arXiv:2206.09200v2, front matter, section 5 (Discussion), section 6 (Acknowledgements), plus cross-reference to the WO 2024/008365 search report (ch10).*

## Core Idea
This preprint is the origin document of the whole corpus: it is where the method first appears in public, where the Biondi–Malanga collaboration is explained, and — because of its posting date — it is the single document that the EPO cited to reject all ten claims of the patent.

## The priority chain

| Date | Event |
|---|---|
| 18 June 2022 | arXiv:2206.09200 **v1** posted — cited by the EPO as XP091253966, "Scanning Inside Volcanoes by Synthetic Aperture Radar Echography Tomographic Doppler Imaging" |
| 04 July 2022 | Patent priority date (PT4451, MT) — **16 days after the preprint** |
| 18 July 2022 | arXiv:2206.09200 **v2** posted — the version in hand, retitled "Scanning Volcanoes by Synthetic Aperture Radar" |
| 29 May 2023 | PCT application PCT/EP2023/064345 filed |
| 31 Aug 2023 | EPO international search completed — category X against claims 1–10 |
| 11 Jan 2024 | WO 2024/008365 A1 published |
| 2022 | Remote Sensing 14, 5231 (Giza paper) published |

**Reading the versions**: the document merged here is v2. The examiner cited the June v1, under a longer title. When citing "the volcano preprint" be explicit about which version — the titles differ, and only v1 is prior art against the patent.

## Frameworks Introduced

- **Attribution of the originating idea**: the acknowledgements state that Prof. Corrado Malanga (Industrial Chemistry, University of Pisa) guided the formulation of the idea of *extrapolating phononic information through photonic processing of SAR data*. Malanga is a co-author of the Giza paper and is absent from the patent — this acknowledgement is the link between the two.
  - Practical use: when attributing the conceptual origin of the method, the phonon-from-photon framing traces to Malanga; the SAR signal processing traces to Biondi.

- **Atmospheric invariance argument**: a single SAR acquisition takes ~14 seconds. Any electromagnetic phase delay from the atmosphere is therefore treated as constant over the acquisition and cancels in the differential measurement; additionally, scanning in the Doppler domain *within a single image* is claimed robust to atmospheric interaction regardless.
  - When to use: defending single-image differential measurements against the usual InSAR atmospheric objection.
  - Limitation to keep: this is an argument from time-invariance, not a measurement of the residual.

- **No commercial tooling exists**: the authors state plainly that no commercial software can extract the phonon information embedded in SAR data, and that they used purpose-built code. SARPROZ (Prof. Daniele Perissin) is credited for supporting calculations, but not as the tomographic engine.
  - Consequence for anyone reproducing this: expect to write the sub-aperture and tomographic focusing stages yourself. The 11-block scheme of ch13 is the spec.

## Data and tool provenance

| Resource | Provider |
|---|---|
| SAR data (CSG spotlight-2A) | Italian Space Agency (ASI) |
| Supporting processing software | SARPROZ, courtesy of Prof. Daniele Perissin |
| In-situ seismic streams | INGV, `https://eida.ingv.it/it/` |
| DEM for topographic validation | SRTM |
| Magnetotelluric tomography | prior published CSAMT work at Vesuvius |
| Tomographic focusing code | authors' own, unpublished |

## Anti-patterns

- **Publishing before filing.** The 16-day gap between v1 and the priority date is the entire novelty problem of the patent. See ch10.
- **Leaving template boilerplate in a released preprint.** Section 2.1, titled "Headings: second level", contains lorem-ipsum filler text ("Fusce mauris. Vestibulum luctus nibh at lectus…") and an unrelated equation (1) — a hidden-Markov-model forward-backward expression `ξij(t) = P(x_t = i, x_t+1 = j | y, v, w; θ)` — inherited from the PRIME AI paper LaTeX template and never removed. It has nothing to do with the method.
  - **Why it matters**: numbered equation (1) in this preprint is not part of the method. Anyone citing "equation (1) of Biondi 2022" from this document would be citing template residue. The real derivation starts at section 2.3.
  - It also signals the document did not receive a final proofing pass — reasonable grounds to prefer the peer-reviewed Giza paper or the patent for formal citation.
- **Citing this preprint as peer-reviewed.** It is an arXiv posting in eess.SP. The Giza paper (Remote Sensing, MDPI) is the peer-reviewed member of the corpus.

## Scope claims made in the discussion

Beyond volcanoes, the discussion proposes the technique could extend to: detecting crude-oil or natural-gas underground pockets; rapid search for metal and rare-earth veins; assessing the material consistency of large infrastructure; and building accurate subsurface models to strengthen volcanological and seismological predictive models. These are stated as prospects, not results — the same epistemic framing to preserve as with the Giza hydraulic hypothesis (ch08).

## Key Takeaways
1. v1 (18 June 2022) is the prior art; v2 (18 July 2022) is this document, retitled — always name the version.
2. The preprint predates the patent priority by 16 days, which is why the EPO issued category X on all claims.
3. Malanga is credited with the phonon-from-photon idea; Biondi with the SAR processing. This acknowledgement links the volcano preprint, the Giza paper, and the patent.
4. Atmospheric delay is argued away as time-invariant over a ~14 s acquisition — an argument, not a measured residual.
5. No commercial software implements this; reproduction requires building the pipeline from the 11-block spec.
6. Section 2.1 and equation (1) of this preprint are unremoved LaTeX template filler — do not cite them.
7. This is the first published claim of SAR-based soil-consistency estimation to 3 km depth from space.

## Connects To
- **Ch 10**: the search report that cites this document.
- **Ch 13**: the 11-block scheme that substitutes for the unavailable commercial tooling.
- **Ch 14, Ch 15**: the results and validation this document reports.
- **Ch 09**: the Giza paper's own patent and funding disclosure.

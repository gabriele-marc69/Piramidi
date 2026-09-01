# Section 5.2–5.3: Internal Tomography — Validation Against Known Structures

## Core Idea
Before presenting newly discovered structures, the authors validate the MM tomography method by re-imaging the *already known* interior features (King's Chamber, Zed, Queen's Chamber, Grotto, unfinished room) and checking that the tomograms line up with the accepted architectural diagram.

## Frameworks Introduced
- **Overlap validation strategy**: for each known feature, produce two views — a tomogram partially overlapped on the accepted schematic, and the same tomogram non-overlapped — so a reader/reviewer can directly judge visual agreement without the schematic biasing interpretation of the raw tomogram.
  - When to use: whenever presenting a new imaging modality's output, pair it against ground truth in both an overlapped and a "blind" rendering to demonstrate the method isn't just reproducing the expected shape by construction.

## Key Concepts
- **Tomographic line orientation**: the internal tomography in this paper uses 8+ different acquisitions (6 different geometries) so each investigates the pyramid from a different side/angle — no single line can reveal the whole interior.
- **Per-pixel metric resolution**: 1 pixel corresponds to 1 m of spatial resolution (paper's stated working assumption for interpreting all Section 5 tomograms).
- **Multipath / layover / foreshortening**: known SAR geometric distortions the authors say they mitigated by deliberately choosing tomographic lines through "pure" single-scattering pixels, rather than by explicit multipath modeling (explicitly deferred to future work).

## Reference Tables
Acquisitions used for internal results (Table 2, internal rows, paraphrased):

| Date | Orbit | Polarization |
|---|---|---|
| 25 Feb 2022 | Left-descending | VV |
| 16 Nov 2021 | Right-descending | HH |
| 22 Feb 2022 | Right-descending | VV |
| 16 Feb 2022 | Right-ascending | VV |
| 25 Mar 2022 | Right-descending | VV |
| 26 Apr 2022 | Right-descending | VV |

## Worked Example
The King's Chamber/Zed complex is imaged from a tomographic line running from ground to apex; the resulting tomogram (shown both overlapped and non-overlapped on the known schematic) also picks up the Queen's Chamber, the "Grotto" void, and — with a weaker signal — the underground "unfinished room," plus a previously-undescribed corridor (later formalized as structure #3 in the new catalog) linking the Grotto to the room below.

## Key Takeaways
1. Known structures (King's Chamber, Zed, Queen's Chamber, Grotto, unfinished room) were all recovered by the tomography, which the authors present as their core validation step.
2. Signal strength varies by structure — the unfinished room registered only weakly, flagged explicitly by the authors rather than presented as equally certain to the strong detections.
3. Six distinct acquisition geometries (E/N/W/S sides, vertical and horizontal line orientations) were needed to build a composite internal 3D picture — no single pass suffices.
4. A full-resolution tomogram computation took ~6 days on a workstation (Dell i7, 32 GB RAM) — a major practical cost of the method (see Section 6.1).

## Connects To
- **Section 5.4–5.19**: builds the new-structure catalog using the same tomographic-line technique validated here.
- **Section 2**: the known-structure baseline diagram this section's tomograms are checked against.

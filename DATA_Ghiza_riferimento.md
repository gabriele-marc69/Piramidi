# DATA_Ghiza — elenco dei file di riferimento

Questa cartella **non sta nel repository**: sono decine di GB di prodotti
Sentinel-1 e il `.gitignore` la esclude. Qui ci sono solo i NOMI, perche' chi
vuole rieseguire `piramidi_v02.py` sappia esattamente cosa deve avere sul disco
e con quale struttura.

Generato il 2026-09-02 interrogando il catalogo CDSE con `scarica_ghiza_cdse.py --list-only`.

## Come ottenerli

```
python scarica_ghiza_cdse.py --out DATA_Ghiza --mode swath \
                             --swath iw2 --polarisation vh --start 2026-01-01
```

Servono le credenziali CDSE in `.cdse.env` (`CDSE_USER` / `CDSE_PASS`), un
account gratuito su <https://dataspace.copernicus.eu>. La modalita' `swath`
scarica solo i file utili invece del prodotto intero: ~1.5 GB per data invece
di ~8 GB.

## Cosa contiene ogni prodotto

Di ciascun `.SAFE` si tengono **cinque file**, non l'archivio completo:

```
S1A_IW_SLC__1SDV_20260103T155650_20260103T155717_062605_07D8F1_C11C.SAFE/
  annotation/calibration/calibration-s1a-iw2-slc-vh-20260103t155650-20260103t155715-062605-07d8f1-002.xml
  annotation/calibration/noise-s1a-iw2-slc-vh-20260103t155650-20260103t155715-062605-07d8f1-002.xml
  annotation/s1a-iw2-slc-vh-20260103t155650-20260103t155715-062605-07d8f1-002.xml
  manifest.safe
  measurement/s1a-iw2-slc-vh-20260103t155650-20260103t155715-062605-07d8f1-002.tiff
```

Cioe', in forma generica:

```
<PRODOTTO>.SAFE/
  manifest.safe
  annotation/<sensore>-iw2-slc-vh-<tempi>-<orbita>-<id>-002.xml
  annotation/calibration/calibration-<stesso nome>.xml
  annotation/calibration/noise-<stesso nome>.xml
  measurement/<stesso nome>.tiff
```

Il `.tiff` e' l'immagine SLC complessa; i due `.xml` di calibrazione servono
alla radiometria (F40); l'`annotation` porta orbita, geolocation grid e
parametri Doppler, cioe' tutta la geometria del calcolo.

## La pila coerente — orbita relativa 58, ascendente

Una tomografia multi-baseline sta su UNA sola orbita relativa. Queste sono le
46 acquisizioni utili nella finestra dal 2026-01-01:

| # | data | sat | dim. prodotto | nome |
|---|------|-----|---------------|------|
| 1 | 2026-01-03 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260103T155650_20260103T155717_062605_07D8F1_C11C.SAFE` |
| 2 | 2026-01-09 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260109T155541_20260109T155609_005829_00BAA7_7A9A.SAFE` |
| 3 | 2026-01-15 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260115T155649_20260115T155716_062780_07DFB2_208B.SAFE` |
| 4 | 2026-01-21 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260121T155540_20260121T155608_006004_00C0A4_92CF.SAFE` |
| 5 | 2026-01-27 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260127T155648_20260127T155715_062955_07E5FA_44AD.SAFE` |
| 6 | 2026-02-02 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260202T155540_20260202T155608_006179_00C66E_D261.SAFE` |
| 7 | 2026-02-08 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260208T155647_20260208T155714_063130_07EC93_D775.SAFE` |
| 8 | 2026-02-14 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260214T155540_20260214T155607_006354_00CC71_5010.SAFE` |
| 9 | 2026-02-20 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260220T155647_20260220T155714_063305_07F31D_EE0F.SAFE` |
| 10 | 2026-02-26 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260226T155540_20260226T155607_006529_00D28A_4CCB.SAFE` |
| 11 | 2026-03-04 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260304T155646_20260304T155713_063480_07F9CD_E7A2.SAFE` |
| 12 | 2026-03-10 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260310T155540_20260310T155608_006704_00D88B_CACB.SAFE` |
| 13 | 2026-03-16 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260316T155647_20260316T155714_063655_080065_7550.SAFE` |
| 14 | 2026-03-22 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260322T155540_20260322T155608_006879_00DE8A_FEA1.SAFE` |
| 15 | 2026-03-28 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260328T155647_20260328T155714_063830_0806F1_F0BB.SAFE` |
| 16 | 2026-04-03 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260403T155540_20260403T155608_007054_00E47B_9B40.SAFE` |
| 17 | 2026-04-09 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260409T155647_20260409T155714_064005_080D78_65C9.SAFE` |
| 18 | 2026-04-15 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260415T155541_20260415T155609_007229_00EA62_166D.SAFE` |
| 19 | 2026-04-21 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260421T155647_20260421T155714_064180_0813FF_BC4E.SAFE` |
| 20 | 2026-04-27 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260427T155541_20260427T155609_007404_00F04E_5B9C.SAFE` |
| 21 | 2026-04-28 | S1D | 8.0 GB | `S1D_IW_SLC__1SDV_20260428T155551_20260428T155619_002549_00439C_DAD7.SAFE` |
| 22 | 2026-05-03 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260503T155647_20260503T155714_064355_081A81_901B.SAFE` |
| 23 | 2026-05-09 | S1C | 7.7 GB | `S1C_IW_SLC__1SDV_20260509T155542_20260509T155610_007579_00F639_1302.SAFE` |
| 24 | 2026-05-10 | S1D | 8.0 GB | `S1D_IW_SLC__1SDV_20260510T155552_20260510T155620_002724_004960_6E00.SAFE` |
| 25 | 2026-05-15 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260515T155646_20260515T155713_064530_0820BD_6DAA.SAFE` |
| 26 | 2026-05-21 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260521T155543_20260521T155611_007754_00FC1D_286A.SAFE` |
| 27 | 2026-05-22 | S1D | 8.0 GB | `S1D_IW_SLC__1SDV_20260522T155553_20260522T155621_002899_004F38_3F0E.SAFE` |
| 28 | 2026-05-27 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260527T155646_20260527T155713_064705_0826DA_37E3.SAFE` |
| 29 | 2026-06-02 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260602T155544_20260602T155611_007929_0101F2_60F1.SAFE` |
| 30 | 2026-06-03 | S1D | 8.0 GB | `S1D_IW_SLC__1SDV_20260603T155554_20260603T155621_003074_005502_0EDC.SAFE` |
| 31 | 2026-06-08 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260608T155645_20260608T155712_064880_082CFB_D308.SAFE` |
| 32 | 2026-06-15 | S1D | 8.0 GB | `S1D_IW_SLC__1SDV_20260615T155554_20260615T155622_003249_005AD1_514C.SAFE` |
| 33 | 2026-06-20 | S1A | 7.7 GB | `S1A_IW_SLC__1SDV_20260620T155644_20260620T155711_065055_083320_7470.SAFE` |
| 34 | 2026-06-27 | S1D | 8.8 GB | `S1D_IW_SLC__1SDV_20260627T155543_20260627T155614_003424_006088_24A7.SAFE` |
| 35 | 2026-06-27 | S1D | 7.7 GB | `S1D_IW_SLC__1SDV_20260627T155612_20260627T155639_003424_006088_BC07.SAFE` |
| 36 | 2026-07-03 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260703T155542_20260703T155610_008381_010953_C590.SAFE` |
| 37 | 2026-07-09 | S1D | 7.7 GB | `S1D_IW_SLC__1SDV_20260709T155556_20260709T155623_003599_00666F_A304.SAFE` |
| 38 | 2026-07-15 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260715T155543_20260715T155611_008556_010F01_1EE7.SAFE` |
| 39 | 2026-07-21 | S1D | 7.7 GB | `S1D_IW_SLC__1SDV_20260721T155557_20260721T155624_003774_006C65_EBDF.SAFE` |
| 40 | 2026-07-27 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260727T155544_20260727T155612_008731_0114D2_C22F.SAFE` |
| 41 | 2026-08-02 | S1D | 7.7 GB | `S1D_IW_SLC__1SDV_20260802T155557_20260802T155624_003949_007286_ADFD.SAFE` |
| 42 | 2026-08-08 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260808T155544_20260808T155612_008906_011AA5_22A3.SAFE` |
| 43 | 2026-08-14 | S1D | 7.7 GB | `S1D_IW_SLC__1SDV_20260814T155558_20260814T155625_004124_007891_88F7.SAFE` |
| 44 | 2026-08-20 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260820T155545_20260820T155613_009081_012067_77F2.SAFE` |
| 45 | 2026-08-26 | S1D | 7.7 GB | `S1D_IW_SLC__1SDV_20260826T155559_20260826T155626_004299_007EB0_F9B2.SAFE` |
| 46 | 2026-09-01 | S1C | 8.0 GB | `S1C_IW_SLC__1SDV_20260901T155546_20260901T155614_009256_012648_FDE1.SAFE` |

## Acquisizioni discendenti (orbita 167) — da NON impilare con le precedenti

Il catalogo le restituisce perche' coprono la stessa area, ma la baseline
rispetto alla pila ascendente e' oltre quella critica: coerenza nulla.
`piramidi_v02.py` le scarta da solo e lo segnala.

| data | sat | dim. prodotto | nome |
|------|-----|---------------|------|
| 2026-01-29 | S1C | 10.6 GB | `S1C_IW_SLC__1SDV_20260129T035132_20260129T035210_006113_00C443_C579.SAFE` |
| 2026-02-10 | S1C | 10.6 GB | `S1C_IW_SLC__1SDV_20260210T035131_20260210T035209_006288_00CA38_FD9B.SAFE` |
| 2026-02-22 | S1C | 10.6 GB | `S1C_IW_SLC__1SDV_20260222T035131_20260222T035209_006463_00D052_503B.SAFE` |
| 2026-03-06 | S1C | 10.6 GB | `S1C_IW_SLC__1SDV_20260306T035131_20260306T035209_006638_00D664_3CEE.SAFE` |
| 2026-05-05 | S1C | 10.6 GB | `S1C_IW_SLC__1SDV_20260505T035134_20260505T035212_007513_00F415_11CD.SAFE` |

## Le 8 date della corsa pubblicata in `out_ghiza_8date/`

Per riprodurre esattamente quei risultati servono queste, e il master e'
il 2026-02-14 (scelto da `pick_supermaster`, non a mano):

| data | sat | ruolo | nome |
|------|-----|-------|------|
| 2026-01-03 | S1A | slave | `S1A_IW_SLC__1SDV_20260103T155650_20260103T155717_062605_07D8F1_C11C.SAFE` |
| 2026-01-09 | S1C | slave | `S1C_IW_SLC__1SDV_20260109T155541_20260109T155609_005829_00BAA7_7A9A.SAFE` |
| 2026-01-15 | S1A | slave | `S1A_IW_SLC__1SDV_20260115T155649_20260115T155716_062780_07DFB2_208B.SAFE` |
| 2026-01-21 | S1C | slave | `S1C_IW_SLC__1SDV_20260121T155540_20260121T155608_006004_00C0A4_92CF.SAFE` |
| 2026-02-14 | S1C | **master** | `S1C_IW_SLC__1SDV_20260214T155540_20260214T155607_006354_00CC71_5010.SAFE` |
| 2026-02-20 | S1A | slave | `S1A_IW_SLC__1SDV_20260220T155647_20260220T155714_063305_07F31D_EE0F.SAFE` |
| 2026-02-26 | S1C | slave | `S1C_IW_SLC__1SDV_20260226T155540_20260226T155607_006529_00D28A_4CCB.SAFE` |
| 2026-03-04 | S1A | slave | `S1A_IW_SLC__1SDV_20260304T155646_20260304T155713_063480_07F9CD_E7A2.SAFE` |

**Attenzione a come si selezionano.** `--dates N` non sceglie *quali* date:
prende le prime N in ordine cronologico fra quelle presenti nella cartella
(`discover_stack(cfg)[:cfg.n_dates]`). Se scarichi tutte le 46 acquisizioni e
lanci `--dates 8`, ottieni le prime otto del catalogo, che NON sono queste, e i
risultati non coincideranno con quelli pubblicati -- verificato.

Per riprodurre esattamente la corsa, metti in una cartella **solo** questi otto
`.SAFE` e punta li' lo stack:

```
python piramidi_v02.py --stack-dir DATA_Ghiza_8date --dates 8 --out out_ghiza_8date
```

Anche il master non si sceglie: `pick_supermaster` (F13) lo deriva dalle
baseline, quindi cambia da solo se cambia l'insieme delle date.

## Quante date servono

Non tutte, ma nemmeno poche. Con 4 baseline la soglia di significativita'
non lascia passare nulla; il budget stampato da `--report-only` dice quanto
vale la risoluzione verticale per il numero di date che hai, prima ancora di
leggere un `.tiff`:

```
python piramidi_v02.py --report-only --stack-dir DATA_Ghiza --dates 8
```

## La corsa a 29 date del 2026-09-02 -- i file effettivamente usati

Comando: `python piramidi_v02.py --stack-dir DATA_Ghiza --dates 99 --out out_piramidi_v02` (log completo in `run_29date.log`, indice dei prodotti in `out_piramidi_v02/meta.json`).

Delle cartelle `.SAFE` presenti in `DATA_Ghiza/`, `discover_stack()` ne ha usate **29**. Nessun filtro e' stato scritto a mano: `_stack_safe()` salta da solo le cartelle senza `.tiff` (vuote, o con lo scarico ancora in `.part`) e `_omogenea()` (F43) toglie le acquisizioni che non stanno sulla traccia della pila.

Il master **non si sceglie**: `pick_supermaster()` (F13) lo deriva dalle baseline, e per questo insieme e' il **2026-03-28** (S1A). La colonna *offset finestra* e' la correzione F42: la distanza, in pixel del prodotto, fra dove cade il bersaglio nel master e dove cade in quella data.

| # | data | sat | ruolo | B_perp | dt | offset finestra (linee, pixel) | burst | prodotto |
|---|------|-----|-------|--------|----|--------------------------------|-------|----------|
| 1 | 2026-01-03 | S1A | slave | +101.6 m | -84 d | -8, -10 | 2 | `S1A_IW_SLC__1SDV_20260103T155650_20260103T155717_062605_07D8F1_C11C.SAFE` |
| 2 | 2026-01-09 | S1C | slave | -23.2 m | -78 d | +4501, +63 | 5 | `S1C_IW_SLC__1SDV_20260109T155541_20260109T155609_005829_00BAA7_7A9A.SAFE` |
| 3 | 2026-01-15 | S1A | slave | +63.0 m | -72 d | -6, -1 | 2 | `S1A_IW_SLC__1SDV_20260115T155649_20260115T155716_062780_07DFB2_208B.SAFE` |
| 4 | 2026-01-21 | S1C | slave | +53.8 m | -66 d | +4502, +34 | 5 | `S1C_IW_SLC__1SDV_20260121T155540_20260121T155608_006004_00C0A4_92CF.SAFE` |
| 5 | 2026-01-27 | S1A | slave | -38.1 m | -60 d | -13, +27 | 2 | `S1A_IW_SLC__1SDV_20260127T155648_20260127T155715_062955_07E5FA_44AD.SAFE` |
| 6 | 2026-02-02 | S1C | slave | +77.8 m | -54 d | +4502, +15 | 5 | `S1C_IW_SLC__1SDV_20260202T155540_20260202T155608_006179_00C66E_D261.SAFE` |
| 7 | 2026-02-08 | S1A | slave | -51.1 m | -48 d | -9, +25 | 2 | `S1A_IW_SLC__1SDV_20260208T155647_20260208T155714_063130_07EC93_D775.SAFE` |
| 8 | 2026-02-14 | S1C | slave | +17.3 m | -42 d | +4502, +40 | 5 | `S1C_IW_SLC__1SDV_20260214T155540_20260214T155607_006354_00CC71_5010.SAFE` |
| 9 | 2026-02-20 | S1A | slave | +0.6 m | -36 d | -2, +4 | 2 | `S1A_IW_SLC__1SDV_20260220T155647_20260220T155714_063305_07F31D_EE0F.SAFE` |
| 10 | 2026-02-26 | S1C | slave | -37.2 m | -30 d | +4502, +50 | 5 | `S1C_IW_SLC__1SDV_20260226T155540_20260226T155607_006529_00D28A_4CCB.SAFE` |
| 11 | 2026-03-04 | S1A | slave | +40.6 m | -24 d | -2, -7 | 2 | `S1A_IW_SLC__1SDV_20260304T155646_20260304T155713_063480_07F9CD_E7A2.SAFE` |
| 12 | 2026-03-10 | S1C | slave | -84.7 m | -18 d | +4501, +62 | 5 | `S1C_IW_SLC__1SDV_20260310T155540_20260310T155608_006704_00D88B_CACB.SAFE` |
| 13 | 2026-03-16 | S1A | slave | -7.8 m | -12 d | +0, +2 | 2 | `S1A_IW_SLC__1SDV_20260316T155647_20260316T155714_063655_080065_7550.SAFE` |
| 14 | 2026-03-22 | S1C | slave | -44.9 m | -6 d | +4504, +53 | 5 | `S1C_IW_SLC__1SDV_20260322T155540_20260322T155608_006879_00DE8A_FEA1.SAFE` |
| 15 | 2026-03-28 | S1A | **master** | +0.0 m | +0 d | +0, +0 | 2 | `S1A_IW_SLC__1SDV_20260328T155647_20260328T155714_063830_0806F1_F0BB.SAFE` |
| 16 | 2026-04-03 | S1C | slave | +18.0 m | +6 d | +4504, +34 | 5 | `S1C_IW_SLC__1SDV_20260403T155540_20260403T155608_007054_00E47B_9B40.SAFE` |
| 17 | 2026-04-09 | S1A | slave | -13.7 m | +12 d | -3, +6 | 2 | `S1A_IW_SLC__1SDV_20260409T155647_20260409T155714_064005_080D78_65C9.SAFE` |
| 18 | 2026-04-15 | S1C | slave | -62.4 m | +18 d | +4504, +57 | 5 | `S1C_IW_SLC__1SDV_20260415T155541_20260415T155609_007229_00EA62_166D.SAFE` |
| 19 | 2026-04-21 | S1A | slave | -42.3 m | +24 d | +3, -2 | 2 | `S1A_IW_SLC__1SDV_20260421T155647_20260421T155714_064180_0813FF_BC4E.SAFE` |
| 20 | 2026-04-27 | S1C | slave | -11.8 m | +30 d | +4505, +24 | 5 | `S1C_IW_SLC__1SDV_20260427T155541_20260427T155609_007404_00F04E_5B9C.SAFE` |
| 21 | 2026-04-28 | S1D | slave | +16.4 m | +31 d | +4535, +15 | 5 | `S1D_IW_SLC__1SDV_20260428T155551_20260428T155619_002549_00439C_DAD7.SAFE` |
| 22 | 2026-05-03 | S1A | slave | +101.5 m | +36 d | -7, -43 | 2 | `S1A_IW_SLC__1SDV_20260503T155647_20260503T155714_064355_081A81_901B.SAFE` |
| 23 | 2026-05-09 | S1C | slave | -75.2 m | +42 d | +4505, +44 | 5 | `S1C_IW_SLC__1SDV_20260509T155542_20260509T155610_007579_00F639_1302.SAFE` |
| 24 | 2026-05-21 | S1C | slave | -28.7 m | +54 d | +4531, +39 | 5 | `S1C_IW_SLC__1SDV_20260521T155543_20260521T155611_007754_00FC1D_286A.SAFE` |
| 25 | 2026-05-22 | S1D | slave | -41.1 m | +55 d | +4534, +40 | 5 | `S1D_IW_SLC__1SDV_20260522T155553_20260522T155621_002899_004F38_3F0E.SAFE` |
| 26 | 2026-05-27 | S1A | slave | +188.2 m | +60 d | -5, -55 | 2 | `S1A_IW_SLC__1SDV_20260527T155646_20260527T155713_064705_0826DA_37E3.SAFE` |
| 27 | 2026-06-08 | S1A | slave | +130.3 m | +72 d | -4, -42 | 2 | `S1A_IW_SLC__1SDV_20260608T155645_20260608T155712_064880_082CFB_D308.SAFE` |
| 28 | 2026-06-15 | S1D | slave | -1.2 m | +79 d | +4535, +29 | 5 | `S1D_IW_SLC__1SDV_20260615T155554_20260615T155622_003249_005AD1_514C.SAFE` |
| 29 | 2026-06-20 | S1A | slave | +106.7 m | +84 d | -6, -42 | 2 | `S1A_IW_SLC__1SDV_20260620T155644_20260620T155711_065055_083320_7470.SAFE` |

Gli offset di ~+4500 linee sono il motivo per cui esiste F42: sono 63 km di
volo lungo l'orbita, e separano i prodotti S1C/S1D da quelli S1A (il segno
dipende solo da quale delle due famiglie fa da master: qui il master e' S1A,
quindi le 16 date S1C/S1D escono positive; con un master S1C uscivano negative
le 13 date S1A). Prima della correzione quelle date venivano lette nel punto
sbagliato del prodotto -- deserto a 63 km di distanza, ripiegato dalla
traslazione circolare e dirampato col burst sbagliato -- e contribuivano solo
rumore.

### Scartate perche' su un'altra traccia (5)

Il verso e il numero di orbita relativa li legge `_traccia()` dal `manifest.safe`. La pila sta sulla **58 ascendente**.

| data | traccia | prodotto |
|------|---------|----------|
| 2026-01-29 | 167 descending | `S1C_IW_SLC__1SDV_20260129T035132_20260129T035210_006113_00C443_C579.SAFE` |
| 2026-02-10 | 167 descending | `S1C_IW_SLC__1SDV_20260210T035131_20260210T035209_006288_00CA38_FD9B.SAFE` |
| 2026-02-22 | 167 descending | `S1C_IW_SLC__1SDV_20260222T035131_20260222T035209_006463_00D052_503B.SAFE` |
| 2026-03-06 | 167 descending | `S1C_IW_SLC__1SDV_20260306T035131_20260306T035209_006638_00D664_3CEE.SAFE` |
| 2026-05-05 | 167 descending | `S1C_IW_SLC__1SDV_20260505T035134_20260505T035212_007513_00F415_11CD.SAFE` |

### Stato alla corsa a 43 date del 2026-09-03 (`out_piramidi_v02/`, `run_43date.log`)

**La cartella e' un bersaglio mobile**: fra due corse a mezz'ora di distanza
la pila e' passata da 28 a 29 date e il master e' cambiato da 2026-03-22 a
2026-03-28; il 2026-09-03, con 44 prodotti completi, il master e' 2026-05-03.
Rilanciare lo stesso comando quando lo scarico e' finito include da solo le
date nuove.

Delle 46 ascendenti in catalogo, 44 erano complete sul disco e **43 sono
entrate in pila**. Fuori:

| data | sat | motivo | prodotto |
|------|-----|--------|----------|
| 2026-06-27 | S1D | seconda fetta della stessa passata, non contiene Giza (F47) | `S1D_IW_SLC__1SDV_20260627T155612_20260627T155639_003424_006088_BC07.SAFE` |
| 2026-05-10 | S1D | `.tiff` mai completato: il server chiude il socket (SSL) e non concede la ripresa, 10 tentativi falliti | `S1D_IW_SLC__1SDV_20260510T155552_20260510T155620_002724_004960_6E00.SAFE` |
| 2026-05-15 | S1A | idem | `S1A_IW_SLC__1SDV_20260515T155646_20260515T155713_064530_0820BD_6DAA.SAFE` |

Per completarle basta rilanciare lo scarico (riprende dai file mancanti) e
poi la corsa:

```
python scarica_ghiza_cdse.py --out DATA_Ghiza --mode swath --swath iw2 \
                             --polarisation vh --start 2026-01-01 --relative-orbit 58
python piramidi_v02.py --stack-dir DATA_Ghiza --out out_piramidi_v02
```

# Piramidi di Giza da Capella — Biondi Protocol v1.9 (`--bp19`)

Due prese spotlight del satellite **Capella-13**, scaricate come **CPHD** (dati di
fase prima della focalizzazione) da **Capella Space Open Data** (licenza CC BY 4.0),
elaborate con il **Biondi Protocol v1.9** (pacchetto pubblicato il 10/10/2026)
portato su GPU. Per Cheope, Chefren e Micerino il programma calcola l'R² corretto
del fit a ellisse dalla superficie della piramide fino a 300 m sotto la base, ogni
0,05 m, e lo raccoglie in cubi di circa 1 m di lato.

**L'R² corretto misura quanto il moto misurato delle coppie Doppler segue l'ellisse
del modello: non è una densità e non prova vuoti.** Il metodo di partenza è stato
ritirato il 10 agosto 2026; i risultati non sono verificati.

Riferimento: il pacchetto
[Replication and Verification — Biondi Protocol](https://github.com/BiondiProtocol/Replication-and-Verification-Biondi-Protocol)
v1.9.0 (BiondiProtocol, 10/10/2026), di cui il programma riprende formule e preset
(`primitives.py`, `measurement.py`, `focus.py`, `statistics.py`). Grazie agli autori
per i consigli.

## I dati

| | 2024-09-16 | 2024-10-04 |
|---|---|---|
| file | `CAPELLA_C13_SP_CPHD_HH_20240916215711_20240916215743.cphd` | `CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012.cphd` |
| satellite, modo | Capella-13, spotlight | Capella-13, spotlight |
| banda, polarizzazione | X, 9,650 GHz, 500 MHz, HH | X, 9,600 GHz, 600 MHz, HH |
| lato, incidenza | destro, 44,9° | sinistro, 38,6° |
| azimut di vista | 315,6° (satellite a NO) | 231,6° (satellite a SO) |
| dimensione | 33,7 GB | 39,2 GB |
| sorgente | `s3://capella-open-data/data/2024/9/16/CAPELLA_C13_SP_CPHD_HH_20240916215711_20240916215743/` | `s3://capella-open-data/data/2024/10/4/CAPELLA_C13_SP_CPHD_HH_20241004001939_20241004002012/` |

Sono gli unici CPHD del catalogo Capella Open Data che contengono Giza.

## I file

| file | cosa è |
|---|---|
| `deepseek_python_cappella_mio.py` | il programma (lettura CPHD, retroproiezione su GPU, `--bp19` e le altre modalità) |
| `bp19_template.html` | modello della pagina 3D: cubi colorati per R² (scala dal nero al giallo, minimo e massimo regolabili), soglia, massimo per colonna, cursore che toglie i cubi dall'alto, parti note, tabella di corrispondenza, tomogramma |
| `bp19_spaccato_template.html` | modello della pagina dello spaccato verticale lungo una retta A-B scelta sulla pianta |
| `voxel_template.html` | modello della pagina della modalità `--voxels` |
| `known_structures.json` | le parti note delle piramidi (wireframe da letteratura, con incertezza), per il confronto |

## Uso

```bash
# una piramide alla volta sulla GPU (Cheope ~1 h 20 min su una GTX 1050 Ti)
python deepseek_python_cappella_mio.py --cphd CAPELLA_C13_SP_CPHD_HH_20240916215711_20240916215743.cphd \
       --bp19 --bp19_pyr khufu --out output_bp19_20240916_khufu

# rifa solo il fit dalle misure salvate (bank.npz, geometry.npz) in output_.../bp19_khufu
python deepseek_python_cappella_mio.py --cphd ... --bp19 --bp19_refit --bp19_pyr khufu --out output_bp19_20240916_khufu
```

Il percorso CPHD predefinito nel programma è quello della macchina di sviluppo:
passare sempre `--cphd`. Dipendenze: `numpy`, `scipy`, `matplotlib`, `torch` (CUDA
consigliata), `sarpy` per leggere i CPHD Capella, `plotly` (per incorporare
plotly.js nelle pagine).

## Uscita di `--bp19`

In `--out`:

- `pyramids_bp19.html` (3D) e `pyramids_bp19_spaccato.html` (spaccato): le pagine
  non contengono i dati, così nessun file supera i 25 MB;
- `pyramids_bp19_dati.json`: indice (titolo, elenco dei file con i cubi, parametri
  della corsa, parti note, tabella di corrispondenza, picco);
- `pyramids_bp19_dati_N.json`: i cubi con R² > 0 in base64 — indice del cubo
  (uint32), R² × 65535 (uint16), indice della profondità del massimo (uint16,
  profondità = 0,1 + 0,05·k m) — e il PNG del tomogramma, in file di al massimo
  20 MB. Ogni `.json` è racchiuso in `bp19_dati(...)` perché la pagina possa
  caricarlo con `<script>` anche aperta dal disco;
- `bp19_<piramide>/`: `bank.npz`, `geometry.npz`, `result.npz` (cubi in float32 e
  profondità del massimo), `report.json`, `tomogramma_linea.png`.

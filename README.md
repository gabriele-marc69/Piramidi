# Piramidi — SAR micro-Doppler tomography del plateau di Giza

Pipeline Python per la **tomografia micro-Doppler** e la **micro-displacement DInSAR** a
partire da uno stack SLC Sentinel-1, ispirata a Biondi & Malanga
(*Remote Sens.* 2022, 14, 5231; [arXiv:2206.09200](https://arxiv.org/abs/2206.09200)).

Il progetto riproduce i *tipi* di immagine pubblicati nell'articolo (sezioni B-scan,
slice di riflettività, mappe di altezza, forme d'onda DInSAR) partendo da dati Sentinel-1
aperti, documentando onestamente il limite fisico: pochi look tomografici su un box piccolo
⇒ bassa altezza di ambiguità (`z_amb ≈ 8.5 m` per il box di Khafre), quindi la profondità
assoluta a scala di piramide non è recuperabile.

## Struttura

```
.
├── skills/                       # Skill riusabili (formato anthropics/skills)
│   ├── sentinel1-slc-reader/     # lettura SLC + geolocalizzazione box via GCP → .npz
│   ├── sar-doppler-tomography/   # sub-aperture Doppler + inversione tomografica
│   └── sar-dinsar-microdisplacement/  # micro-spostamenti LOS da fase interferometrica
├── goal_out/                     # script della pipeline e output (output non versionati)
│   ├── goal_pipeline.py          # pipeline a 6 step (vedi skills/goal.txt)
│   ├── piramide_unificato.py     # pipeline unificata
│   ├── trasformata_fourier.py    # parametrizzazione Fourier degli strati
│   ├── tomografia_verticale.py   # sezioni verticali
│   └── ...                       # grafici 3D, picchi, variazioni di frequenza
├── requirements.txt
└── .gitignore
```

## Pipeline (6 step)

1. Scarica i dati satellitari Sentinel-1 alle coordinate richieste.
2. Ritaglia l'area dai file `.tiff` (geolocalizzazione box via GCP).
3. Genera un array 3D `(x, y, strati)` degli spostamenti ricavati dai `.tiff`.
4. Parametrizza gli strati `y` come funzioni di Fourier stirate su un range di 1000 m.
5. Disegna un grafico 3D dinamico delle funzioni di Fourier generate.
6. Estrae i punti con variazione di frequenza (e quindi di densità).

## Requisiti

```bash
pip install -r requirements.txt
```

Python ≥ 3.10. Dipendenze principali: `numpy`, `scipy`, `matplotlib`, `rasterio`, `plotly`.

## Scaricare i dati Sentinel-1

I dati di input sono prodotti **Sentinel-1 IW SLC** (Single Look Complex): gli unici
dati Sentinel-1 aperti che conservano la **fase**, indispensabile per tomografia e InSAR
(i prodotti GRD sono solo ampiezza → inutili qui). Si scaricano dal
**Copernicus Data Space Ecosystem (CDSE)**, gratuito.

### 1. Crea un account CDSE (gratuito)

Registrati su <https://dataspace.copernicus.eu> → otterrai username (email) e password.

### 2. Imposta le credenziali come variabili d'ambiente

La pipeline legge le credenziali da `CDSE_USER` / `CDSE_PASS` (oppure dai flag
`--cdse-user` / `--cdse-pass`). **Non committare mai le credenziali nel repo.**

PowerShell (Windows):
```powershell
$env:CDSE_USER = "tua-email@example.com"
$env:CDSE_PASS = "la-tua-password"
```

Bash (Linux/macOS):
```bash
export CDSE_USER="tua-email@example.com"
export CDSE_PASS="la-tua-password"
```

### 3. Cerca e scarica i prodotti per la tua area (AOI)

`goal_out/piramide_unificato.py` esegue ricerca (OData CDSE) + download dei prodotti
SLC che **contengono** il box di coordinate richiesto. L'AOI si passa in DMS
(gradi-minuti-secondi) con gli angoli NW e SE:

```bash
python goal_out/piramide_unificato.py \
  --nw 29 58 38.0 N  --nw-lon 31 7 45.4 E \
  --se 29 58 29.0 N  --se-lon 31 7 55.4 E \
  --start 2026-01-01 --end 2026-06-25 \
  --pol vv --layers 12 \
  --download                      # scarica davvero (servono le credenziali CDSE)
```

- Senza `--download` la ricerca viene comunque eseguita, ma **non** scarica i file:
  riusa lo stack già presente in `goal_out/stack_slc/`.
- I valori di default puntano già al box della piramide di **Khafre (Chefren)** sul
  plateau di Giza.
- ⚠️ Ogni prodotto `.SAFE` pesa **~8 GB**: assicurati di avere spazio e banda.

### 4. (Alternativa) lettura/ritaglio di un singolo TIFF già scaricato

Se hai già i `.tiff` SLC + i rispettivi `*.annotation.xml`, puoi estrarre direttamente
il chip complesso co-registrato nel box con la skill `sentinel1-slc-reader`:

```bash
python skills/sentinel1-slc-reader/scripts/extract_box.py \
  --stack /percorso/stack_slc \
  --nw 29 58 38.0 N 31 7 45.4 E \
  --se 29 58 29.0 N 31 7 55.4 E \
  --pol vv --out box.npz
```

Lo script abbina ogni TIFF alla sua annotazione, geolocalizza il box via GCP e legge
**solo** la finestra del box (mai l'intera scena da ~1.5 GB).

> 💡 In alternativa puoi scaricare i prodotti a mano dal browser CDSE
> (<https://browser.dataspace.copernicus.eu>) filtrando per *Sentinel-1 → SLC → IW*,
> poi scompattare i `.SAFE` in `goal_out/stack_slc/` ed eseguire la pipeline con
> `--steps 3-6` (solo elaborazione, dati locali).

## Dati

> ⚠️ I dati grezzi (`.SAFE`, `.tiff`, `.zip`), gli array intermedi (`.npz`) e i render
> (`.png`, `.html`) **non sono versionati** (vedi `.gitignore`) per via delle dimensioni.

Dati target: stack di 6 scene VV Sentinel-1 IW2 (mag–giu 2026, plateau di Giza).

## Licenza

Distribuito sotto licenza MIT — vedi [LICENSE](LICENSE).

---

*Questo lavoro riproduce tipologie di immagine a scopo di studio e non costituisce
evidenza di strutture sotterranee a scala di piramide.*

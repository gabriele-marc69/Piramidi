# Piramidi — SAR micro-Doppler tomography del plateau di Giza

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Data: Sentinel-1 SLC](https://img.shields.io/badge/Data-Sentinel--1%20SLC-1f6feb.svg)](https://dataspace.copernicus.eu)
[![Made with NumPy](https://img.shields.io/badge/SciPy-NumPy%20%7C%20rasterio%20%7C%20plotly-013243.svg)](requirements.txt)

Pipeline Python per la **tomografia micro-Doppler** e la **micro-displacement DInSAR** a
partire da uno stack SLC Sentinel-1, ispirata a Biondi & Malanga
(*Remote Sens.* 2022, 14, 5231; [arXiv:2206.09200](https://arxiv.org/abs/2206.09200)).

Il progetto riproduce i *tipi* di immagine pubblicati nell'articolo (sezioni B-scan,
slice di riflettività, mappe di altezza, forme d'onda DInSAR) partendo da dati Sentinel-1
aperti, documentando onestamente il limite fisico: pochi look tomografici su un box piccolo
⇒ bassa altezza di ambiguità (`z_amb ≈ 8.5 m` per il box di Khafre), quindi la profondità
assoluta a scala di piramide non è recuperabile.

## Quickstart

```bash
# 1. Clona il repository
git clone https://github.com/gabriele-marc69/Piramidi.git
cd Piramidi

# 2. Crea un ambiente virtuale e installa le dipendenze
python -m venv .venv
# Windows:  .venv\Scripts\activate
# Linux/mac: source .venv/bin/activate
pip install -r requirements.txt

# 3. Imposta le credenziali Copernicus CDSE (account gratuito su dataspace.copernicus.eu)
#    PowerShell:  $env:CDSE_USER="email"; $env:CDSE_PASS="password"
export CDSE_USER="tua-email@example.com"
export CDSE_PASS="la-tua-password"

# 4. Esegui la pipeline completa (ricerca + download + 6 step).
#    UN SOLO script per entrambe le piramidi: --pyramid sceglie il box di default
#    (cheope|kefren) e la cartella di output (goal_out_Cheope/ o goal_out_kefren/).
python piramide_unificato.py --pyramid cheope --download   # Cheope (Khufu)
python piramide_unificato.py --pyramid kefren --download   # Kefren (Khafre)

# Solo elaborazione su dati già presenti nello stack_slc/ della cartella (niente download):
python piramide_unificato.py --pyramid kefren --steps 3-6

# Tomografia Doppler "vera" del paper (steering matrix + h(z)=A^H*Y), in aggiunta
# agli step scelti — vedi la nota metodologica in cima a piramide_unificato.py:
python piramide_unificato.py --pyramid kefren --steps 3 --vera-tomografia
```

Output (grafici 3D `.html`/`.png`, array `.npz`) finiscono nella sottocartella di output
della rispettiva piramide (`unificato/`, non versionata).
Per AOI, date e polarizzazione personalizzate (box custom) passa esplicitamente
`--nw/--nw-lon/--se/--se-lon` al posto di `--pyramid`; vedi
[Scaricare i dati Sentinel-1](#scaricare-i-dati-sentinel-1).

> ⚠️ **Nota sul metodo**: gli step 5-6 (`piramide_unificato.py`) sono un'analisi
> esplorativa (Fourier stirata + frequenza istantanea di Hilbert) ispirata al lessico
> di Biondi & Malanga ma **non** implementano la loro tomografia Doppler a
> sub-aperture. Per la tomografia coerente col paper (steering matrix, inversione
> `h(z)=A^H*Y`, validata su un modello sintetico) usa `--vera-tomografia` oppure
> `skills/sar-doppler-tomography/scripts/tomographic_images.py` direttamente.

## Visualizzazioni interattive

I grafici 3D sono salvati come HTML Plotly interattivi (ruota/zoom nel browser).
GitHub non li renderizza nella vista del repo, quindi aprili tramite un proxy che li
serve come pagina:

- 👁️ **Step 6 — PIENO/VUOTO per frequenza istantanea (3D):**
  [apri nel browser](https://raw.githack.com/gabriele-marc69/Piramidi/main/step6_variazioni_3d.html)
  · [pagina GitHub Pages](https://gabriele-marc69.github.io/Piramidi/step6_variazioni_3d.html)
  *(quest'ultima attiva dopo aver abilitato GitHub Pages su branch `main`/root)*

## Struttura

```
.
├── piramide_unificato.py         # UNICO script per entrambe le piramidi (--pyramid cheope|kefren)
├── skills/                       # Skill riusabili (formato anthropics/skills)
│   ├── sentinel1-slc-reader/     # lettura SLC + geolocalizzazione box via GCP → .npz
│   ├── sar-doppler-tomography/   # sub-aperture Doppler + inversione tomografica (metodo del paper)
│   └── sar-dinsar-microdisplacement/  # micro-spostamenti LOS da fase interferometrica
├── goal_out_Cheope/              # SOLO output per la piramide di Cheope (Khufu)
│   └── ...                       # box.npz, stack_slc/, unificato/ (non versionati)
├── goal_out_kefren/              # SOLO output per la piramide di Kefren (Khafre)
│   └── ...                       # (output .npz/.png/.html non versionati)
├── requirements.txt
└── .gitignore
```

`piramide_unificato.py` sostituisce le due copie divergenti che vivevano prima dentro
`goal_out_Cheope/` e `goal_out_kefren/`, oltre agli script standalone superati
(`goal_pipeline.py`, `forma_donda_rettangolo.py`, `tomografia_verticale.py`,
`trasformata_fourier.py`, `grafico_*.py`, `picchi_onde.py`, `variazioni_onde.py`,
`fetch_dem.py`) che erano duplicati identici in entrambe le cartelle.

## Pipeline (6 step)

1. Scarica i dati satellitari Sentinel-1 alle coordinate richieste.
2. Ritaglia l'area dai file `.tiff` (geolocalizzazione box via GCP).
3. Genera un array 3D `(x, y, strati)` degli spostamenti ricavati dai `.tiff`.
4. Parametrizza gli strati `y` come funzioni di Fourier stirate su un range di 1000 m.
5. Disegna un grafico 3D dinamico delle funzioni di Fourier generate (colore =
   frequenza istantanea: PIENO/alta frequenza, VUOTO/bassa frequenza).
6. Per ogni pixel dell'area selezionata sceglie il campione con inviluppo massimo
   e lo classifica PIENO/VUOTO in base alla sua frequenza istantanea (un punto
   per pixel, non un sottoinsieme filtrato).

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

`piramide_unificato.py` esegue ricerca (OData CDSE) + download dei prodotti
SLC che **contengono** il box di coordinate richiesto. L'AOI si passa in DMS
(gradi-minuti-secondi) con gli angoli NW e SE:

```bash
python goal_out_kefren/piramide_unificato.py \
  --nw 29 58 38.0 N  --nw-lon 31 7 45.4 E \
  --se 29 58 29.0 N  --se-lon 31 7 55.4 E \
  --start 2026-01-01 --end 2026-06-25 \
  --pol vv --layers 12 \
  --download                      # scarica davvero (servono le credenziali CDSE)
```

- Senza `--download` la ricerca viene comunque eseguita, ma **non** scarica i file:
  riusa lo stack già presente in `stack_slc/` della cartella.
- I box di default sono già impostati per piramide: **Cheope (Khufu)** in
  `goal_out_Cheope/` e **Khafre (Chefren)** in `goal_out_kefren/`, entrambe sul
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
> poi scompattare i `.SAFE` nello `stack_slc/` della cartella della piramide ed eseguire la pipeline con
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

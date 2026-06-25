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

## Dati

> ⚠️ I dati grezzi (`.SAFE`, `.tiff`, `.zip`), gli array intermedi (`.npz`) e i render
> (`.png`, `.html`) **non sono versionati** (vedi `.gitignore`) per via delle dimensioni.

Dati target: stack di 6 scene VV Sentinel-1 IW2 (mag–giu 2026, plateau di Giza).

## Licenza

Distribuito sotto licenza MIT — vedi [LICENSE](LICENSE).

---

*Questo lavoro riproduce tipologie di immagine a scopo di studio e non costituisce
evidenza di strutture sotterranee a scala di piramide.*

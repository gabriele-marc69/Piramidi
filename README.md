# Tomografia SAR della piana di Giza — Sentinel-1 IW SLC

Implementazione e verifica critica del metodo di **tomografia Doppler a
micro-moto** di Filippo Biondi e Corrado Malanga, applicato alle tre piramidi
di Giza con dati **Sentinel-1 IW SLC, canale VH**, scaricati dal Copernicus
Data Space Ecosystem.

Il programma non prova a confermare le fonti: prova a **misurare** che cosa i
dati Sentinel-1 permettono davvero di dire, e a dichiarare i limiti dove
esistono. La risposta breve, dopo 29 acquisizioni, è che la superficie della
piana si ricostruisce con precisione metrica ma **le piramidi non vengono
ricostruite**, e la ragione è geometrica, non algoritmica.

## I programmi

| file | cosa fa |
|------|---------|
| `piramidi_v02.py` | catena principale: stack interferometrico multi-baseline, periodogramma in quota, superficie reale, banco di sub-aperture Doppler, uscita 3D interattiva |
| `piramidi_v01.py` | prima versione, ancora usata come libreria (lettura annotation, deramping TOPS, blocchi 2-7 delle fonti) |
| `scarica_ghiza_cdse.py` | scarico selettivo da CDSE: prende i 5 file utili di ogni `.SAFE` invece del prodotto intero (~1.5 GB per data invece di ~8 GB) |
| `piramide_cheope_3d.py`, `piramide_kefren_3d.py`, `piramide_3d_comune.py` | rendering 3D delle singole piramidi |

### Uso

```bash
# budget e baseline senza leggere un solo .tiff -- sempre da qui
python piramidi_v02.py --stack-dir DATA_Ghiza --report-only

# corsa completa su tutte le date presenti (default: --dates 0 = tutte, F48)
python piramidi_v02.py --stack-dir DATA_Ghiza --out out_piramidi_v02

# elenco delle correzioni applicate, e autotest delle convenzioni di segno
python piramidi_v02.py --fixes
python piramidi_v02.py --selftest
```

I dati **non** stanno nel repository: sono decine di GB. `DATA_Ghiza_riferimento.md`
elenca i prodotti, la struttura attesa sul disco, quelli effettivamente usati
nell'ultima corsa e come riscaricarli.

## Che cosa dicono i dati (43 date, master 2026-05-03)

| grandezza | valore |
|---|---|
| acquisizioni impilate | 43 (S1A + S1C + S1D, traccia 58 ascendente, 2026-01-03 → 2026-09-01) |
| escursione baseline ortogonali | 297,1 m |
| risoluzione verticale `δ_z` (Rayleigh) | 132,5 m |
| coerenza mediana misurata | 0,419 |
| precisione sulla quota `σ_h` | 8,3 m (dispersione robusta misurata sulla piana 12,3 m) |
| soglia di qualità dalla distribuzione **nulla** | 0,432 |
| celle sopra soglia | 3901 / 8370 (46,6 %) |
| pendenza Theil-Sen misurato *vs* simulato | **+0,041**, IC95 [+0,005, +0,077] — attesa 1,0 |

**La catena funziona** — il 47 % delle celle supera una soglia calibrata per
Monte Carlo sulla distribuzione nulla del periodogramma (contro l'1 % atteso
per caso), e la piana esce piatta sul datum degli `annotation.xml`.

**Le piramidi no.** Le facce a ~52° superano l'angolo di incidenza di 37°: sono
in **layover pieno**, con fino a 700 punti di superficie ripiegati nella stessa
cella di risoluzione. Con `δ_z` = 132 m non esiste un diffusore dominante da
localizzare, e il periodogramma restituisce il centro di fase della miscela.
È geometria, non un difetto del processing. La pendenza Theil-Sen esclude lo
zero ma vale un venticinquesimo di quella attesa: le mediane per fascia
salgono da +3,4 a +8,5 m per quote simulate da 5 a oltre 90 m.

La coerenza a 0,42 è il frutto della correzione **F46**: il residuo sub-pixel
di coregistrazione non viene più dal picco della cross-correlazione di fase
(che sul VH di Giza sbagliava di un pixel intero e abbassava la coerenza) ma
dallo spostamento che massimizza la coerenza con il master. Su 41 di queste
date la versione precedente dava coerenza 0,29, `σ_h` 12,7 m e 28 % di celle
sopra soglia. Resta vero che aggiungere date allarga le baseline ma porta
decorrelazione temporale (±120 giorni, tre satelliti): circa 12 date su 42
stanno al pavimento dello stimatore anche dopo la correzione.

Il datum verticale letto dagli `annotation.xml` dipende dal prodotto: la
bilineare locale sulle piramidi vale 64 m nei prodotti S1A e 44 m in altri
(F49). Il layer del suolo usa quello del master, con questa incertezza.

## Metodo e fonti

Il programma segue le fonti alla lettera dove sono applicabili — banda di
guardia `B_DL = B_cD/2`, `B_shift` come selettore della frequenza meccanica,
`N_D` come frequenza di campionamento della vibrazione, schema a 11 blocchi con
la FFT2 diretta calcolata una sola volta fuori dal ciclo, `δ_z = λR/(2A)` con λ
**acustica**, protocollo di validazione a tre livelli — e **dichiara le
divergenze** invece di nasconderle. La principale: con il TOPS di Sentinel-1 la
banda Doppler è di 313 Hz contro i ~22 kHz dello spotlight delle fonti, quindi
la profondità qui viene dalle **baseline orbitali**, non dal micro-moto, che
resta un attributo di superficie.

* F. Biondi, C. Malanga, *Synthetic Aperture Radar Doppler Tomography Reveals
  Details of Undiscovered High-Resolution Internal Structure of the Great
  Pyramid of Giza*, Remote Sensing 2022, **14**, 5231
* WO 2024/008365 A1 — domanda **pubblicata**, con rapporto di ricerca di
  categoria X su tutte e 10 le rivendicazioni: va citata come divulgazione di
  un metodo, mai come brevetto concesso
* arXiv:2206.09200 — *Scanning Volcanoes by Synthetic Aperture Radar*

## Avvertenza epistemica

La superficie ricostruita e la separabilità piramidi/piana sono **misure**, con
il loro errore. L'indice di solidità è un discriminante multi-attributo
dichiarato come tale, **non** una rilevazione di cavità risolta in profondità:
le camere note misurano metri, due ordini di grandezza sotto la risoluzione di
questi dati.

## Dati e licenze

Dati Sentinel-1 © Copernicus / ESA, distribuiti sotto la licenza del
[Copernicus Data Space Ecosystem](https://dataspace.copernicus.eu). Gli
articoli e il brevetto citati non sono ridistribuiti qui: restano presso i
rispettivi editori.

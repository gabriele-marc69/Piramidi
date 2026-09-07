# Tomografia SAR della piana di Giza — Sentinel-1 IW SLC

Implementazione e verifica critica del metodo di **tomografia Doppler a
micro-moto** di Filippo Biondi e Corrado Malanga, applicato alle tre piramidi
di Giza con dati **Sentinel-1 IW SLC, canale VH**, scaricati dal Copernicus
Data Space Ecosystem.

Il programma non prova a confermare le fonti: prova a **misurare** che cosa i
dati Sentinel-1 permettono davvero di dire, e a dichiarare i limiti dove
esistono. La risposta breve, dopo 43 acquisizioni, è che la superficie della
piana si ricostruisce con precisione metrica ma **le piramidi non vengono
ricostruite**, e la ragione è geometrica, non algoritmica.

## I programmi

| file | cosa fa |
|------|---------|
| `piramidi_v03.py` | terza versione: da **modulo e fase radiante** di ogni pixel complesso a una **sinusoide verticale per data**, somma, **FFT** (NUFFT di tipo 1) per il profilo in quota, indice **pieno/vuoto** come residuo rispetto alla PSF di un solo diffusore, nuovo 3D interattivo e disegno delle onde sotto Cheope e Chefren |
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

# v03: sinusoidi da modulo e fase, FFT, pieno/vuoto, nuovo 3D
python piramidi_v03.py --stack-dir DATA_Ghiza --out out_piramidi_v03
python piramidi_v03.py --novita
python piramidi_v03.py --selftest
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
decorrelazione temporale (±120 giorni, tre satelliti): 11 date su 42
restano al pavimento dello stimatore (coerenza 0,18-0,19) anche dopo la
correzione.

Il datum verticale letto dagli `annotation.xml` dipende dal prodotto: la
bilineare locale sulle piramidi spazia fra 44,0 e 64,0 m sulle 43 date (F49).
Non è una differenza fra missioni — il nodo grezzo più vicino a Cheope vale
64 m in tutti e 43 i prodotti impilati, S1A, S1C e S1D indifferentemente — ma
di reticolo: i prodotti non condividono la stessa geolocation grid, quindi la
cella che contiene il ritaglio non è la stessa e i suoi quattro nodi cambiano
da una data all'altra. Il layer del suolo usa il datum del master, con questa
incertezza.

## v03 — sinusoidi, FFT e pieno/vuoto

`piramidi_v03.py` riscrive l'inversione nel dominio in cui la si puo'
disegnare. Da ogni pixel dell'interferogramma multilooked prende **modulo** e
**fase in radianti**; ogni data diventa una sinusoide lungo la verticale,
`s_i(z) = A_i cos(φ_i − κ_i z)` con `κ_i = ±k_z,i`; la somma delle sinusoidi
viene trasformata con una **FFT vera** — le righe stanno a `κ_i`, che non è un
reticolo uniforme perché le baseline non lo sono, quindi vengono spalmate su
una griglia uniforme con un nucleo gaussiano e deapodizzate (NUFFT di tipo 1).

Non è un secondo metodo: la somma delle sinusoidi **è** `Re[h(z)]` del
periodogramma di v02 e il suo inviluppo è `|h(z)|`. La verifica è nei numeri,
non nell'argomento: profilo identico al calcolo diretto entro `1,5·10⁻⁵`, e
quota finale identica a quella di v02 entro **1 cm sul 100 % delle celle**,
passando per lo stesso estimatore. Cambia il costo: 0,7 s contro i ~5 s del
doppio ciclo sullo stesso volume.

L'indice **pieno/vuoto** è il residuo fra il profilo misurato e la risposta
attesa da **un solo diffusore** alla quota del picco (stesse ampiezze, stesse
baseline, fase puramente geometrica), normalizzato dalla dispersione misurata
sulle celle di sola piana: un z-score con il suo nullo empirico, come F22 per
`γ`.

| grandezza | valore |
|---|---|
| soglia \|z-score\| (p99 del nullo, 3663 celle di piana) | 3,72 |
| celle anomale — piramidi *vs* piana | 28,4 % *vs* 19,9 % (z = **+3,15**) |
| solo **vuoto** — piramidi *vs* piana | 0,4 % *vs* 0,3 % (z = +0,33) |
| solo **pieno** — piramidi *vs* piana | 27,5 % *vs* 19,3 % (z = +3,09) |
| correlazione dell'eccesso con il layover / ampiezza / coerenza / quota simulata | +0,07 / −0,23 / −0,25 / −0,08 |

Le celle delle piramidi mostrano davvero più energia di quanta la PSF di un
solo diffusore ne spieghi, e **solo** in eccesso: sul difetto le due
popolazioni sono indistinguibili. Ma nessuna delle grandezze di controllo
spiega quell'eccesso — la correlazione più forte vale il 6 % della varianza.
La lettura corretta è che il modello a un diffusore è troppo semplice per
celle in layover pieno, non che sotto le piramidi si sia visto qualcosa: con
`δ_z` = 132 m sotto quella scala non c'è nulla da risolvere.

Le due colonne disegnate a parte (`sezione_onde_sotto_piramidi.png`) mettono
la sagoma della piramide in scala vera sopra il suolo e, sotto, le 43
sinusoidi, la loro somma, l'inviluppo e la barra della cella di Rayleigh. Il
diffusore dominante sta a +9,4 m sopra il riferimento sia per Cheope sia per
Chefren: appena sopra il deserto, non a metà della piramide.

### Il campo di onde, pixel per pixel

`onde_3d_superficie_piramidi.png` e `onde_3d_piramidi.html` portano la stessa
cosa su **tutte** le celle della superficie delle piramidi in geometria radar:
706 colonne (Cheope 335, Chefren 301, Micerino 70; 236 sopra la soglia di
qualità), ognuna la sinusoide risultante `Re[h(z)]` del suo pixel,
normalizzata al proprio massimo e disegnata alla propria posizione est/nord.
Lo spostamento orizzontale è un artificio — l'ampiezza di un interferogramma
non è una lunghezza — e nella pagina interattiva è un cursore, così si vede
che cambiando scala cambia il disegno e non il dato.

Guardato di taglio, il piano dei diffusori dominanti sta appena sopra il
deserto per tutte e tre le piramidi e non segue in alcun modo il profilo delle
facce. È l'osservazione centrale di v02 — la pendenza misurato *vs* simulato
compatibile con zero — disegnata invece che riassunta in un numero.

## Pagine pubblicate

Le due catene stanno su due progetti distinti, ognuno con il proprio indice.
I dati della corsa (riassunti JSON e budget) sono **dentro** la pagina, non in
file da scaricare a parte, e i numeri dell'indice sono letti da quei riassunti
invece che riscritti a mano.

| | |
|---|---|
| **v02** — superficie e tomogramma | <https://tomografia-giza-v02.gabrielemarchini69.workers.dev> |
| **v03** — sinusoidi, FFT e pieno/vuoto | <https://tomografia-giza-v03.gabrielemarchini69.workers.dev> |

Il progetto v03 porta anche **una pagina per piramide** (`onde_cheope.html`,
`onde_chefren.html`, `onde_micerino.html`): le sole colonne di quella piramide
in 3D e, sotto, due grafici SVG disegnati dentro la pagina — la somma delle
sinusoidi della colonna rappresentativa con il suo inviluppo, e il profilo
medio incoerente di tutte le sue celle.

```bash
python pubblica_sito.py --profilo v03 --out sito_v03
python pubblica_sito.py --profilo v02 --out sito_v02
```

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

### Un errore nelle fonti

Il calcolo della risoluzione tomografica dell'articolo di Giza non torna con i
suoi stessi parametri. La stessa frase fissa la frequenza di indagine a
12&nbsp;500 Hz e poi scrive `λ = v/f ≈ 6000/25 000 ≈ 0,24 m`: divide per `2f`.
Dichiara l'apertura orbitale "circa 42&nbsp;000 m" e poi sostituisce `2·84 000`
al denominatore: il doppio. I due scarti spingono nella stessa direzione.
Rifacendo `δ_z = λR/(2A)` con i valori dichiarati — v = 6000 m/s,
f = 12&nbsp;500 Hz, A = 42&nbsp;000 m, R = 650&nbsp;000 m — viene **3,7 m**,
quattro volte più grossolano degli 0,92 m pubblicati; e il paragrafo 5.2 dello
stesso articolo dichiara ancora un'altra cifra, 1 m per pixel.

Le catene del preprint del Vesuvio (36 m) e del brevetto (1,30 m) tornano
invece riga per riga. Gli 0,92 m vanno citati come *la costante con cui
l'articolo ha scalato le sue figure*, non come un risultato che i suoi
parametri sostengono: tutte le misure del catalogo di 20 strutture poggiano su
quella scala. L'elenco completo degli errata — compresi i 22 contro 24 kHz
della banda Doppler, le due sigle per la stessa banda di guardia nel brevetto e
`N_c` contro `N_D` per gli stessi scorrimenti — sta in
`skills/biondi-malanga-sar-tomography/SKILL.md`.

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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cerca_cosmo_skymed.py
======================
Cerca PUNTATORI PUBBLICI (pagine, dataset, metadati, DOI) a materiale che
cita/riguarda prodotti COSMO-SkyMed (CSK/CSG, costellazione SAR italiana
dell'Agenzia Spaziale Italiana - ASI, commercializzata da e-GEOS), su:

  1. cataloghi accademici/open-data con API pubbliche ufficiali (Zenodo,
     OpenAIRE) -- spesso contengono dataset/derivati COSMO-SkyMed che un
     motore di ricerca generico indicizza poco o mette in fondo ai
     risultati ("web poco indicizzato" nel senso di scarsa visibilita' nei
     motori generici, non la rete Tor/onion);
  2. il web generico (fonte 'web', attiva di default) via un'API di ricerca
     UFFICIALE: Brave Search API (chiave personale gratuita in
     BRAVE_SEARCH_KEY, consigliata) oppure Bing Web Search (BING_SEARCH_KEY,
     API ritirata da Microsoft nell'agosto 2025: solo accessi legacy). La
     fonte 'web' esegue ANCHE le query mirate site:drive.google.com e
     site:dropbox.com (piu' altri siti di condivisione con
     --siti-condivisione), sia con le keyword sia coi prefissi dei nomi file
     CSKS1..4/CSG_SSAR1-2. Senza chiave questa fonte viene saltata con un
     avviso: MAI richieste "a mano" che aggirano un anti-bot/CAPTCHA --
     DuckDuckGo e Google, ad es., bloccano le richieste automatizzate senza
     API a pagamento con una verifica anti-bot, e questo programma non tenta
     di eluderla.

*** IMPORTANTE -- COSMO-SkyMed NON e' un dato aperto come Sentinel-1 ***
COSMO-SkyMed (CSK/CSG) e' un prodotto ad accesso ristretto di proprieta'
dell'Agenzia Spaziale Italiana (ASI) e commercializzato da e-GEOS. L'accesso
alle scene vere e proprie richiede un account autorizzato (es. il programma
scientifico "COSMO-SkyMed Open Call" di ASI, o un contratto commerciale con
e-GEOS): vedi https://www.asi.it e https://www.e-geos.it.

Questo programma:
  - NON scarica scene SAR e NON bypassa alcun controllo di accesso/licenza;
  - NON cerca su reti anonime (Tor/onion) ne' tenta di aggirare CAPTCHA o
    protezioni anti-bot di alcun servizio;
  - raccoglie SOLO puntatori (URL/DOI) a pagine/dataset GIA' pubblicati
    apertamente da terzi (articoli, dataset accademici, repository
    istituzionali) che citano COSMO-SkyMed, tramite le API pubbliche
    ufficiali dei cataloghi interrogati.
Verifica sempre licenza e provenienza di quanto trovi prima di riusarlo.

Oltre alle query keyword ("COSMO-SkyMed <keywords>"), il programma interroga
di default le fonti anche con i NOMI FILE dei prodotti COSMO-SkyMed secondo
la convenzione ufficiale di denominazione (vedi il blocco "nomi file
COSMO-SkyMed" piu' sotto): prefissi satellite CSKS1..CSKS4 (1a generazione) e
CSG_SSAR1/CSG_SSAR2 (2a generazione) e TUTTE le combinazioni con i tipi
prodotto RAW_B/SCS_B/SCS_U/DGM_B/GEC_B/GTC_B (es. CSKS1_SCS_B) -- cosi'
vengono trovati anche dataset/pagine che elencano le scene per nome file
senza mai scrivere "COSMO-SkyMed" per esteso. I tipi prodotto DA SOLI non
vengono interrogati (troppo rumore nei full-text dei cataloghi). Si regola
con --nomi-file tutti|base|nessuno (default: tutti).

Uso:
  python cerca_cosmo_skymed.py --keywords "pyramid Giza"
  python cerca_cosmo_skymed.py --keywords "Khufu Cheope" --sources zenodo,openaire
  python cerca_cosmo_skymed.py --nomi-file base        # meno query sui nomi file
  python cerca_cosmo_skymed.py --nomi-file nessuno     # solo query keyword
  python cerca_cosmo_skymed.py --nw 29 58 48.0 N --nw-lon 31 7 58.4 E \
                                --se 29 58 41.0 N --se-lon 31 8 07.7 E
  BING_SEARCH_KEY=... python cerca_cosmo_skymed.py --sources web,zenodo,openaire

Le fonti disponibili (--sources, default "zenodo,openaire,web"):
  zenodo   - https://zenodo.org/api/records (nessuna chiave richiesta)
  openaire - https://api.openaire.eu/search/publications (nessuna chiave richiesta)
  web      - Brave Search API (BRAVE_SEARCH_KEY, gratuita) o Bing legacy
             (BING_SEARCH_KEY); chiavi MAI hardcoded qui. Include le query
             site:drive.google.com e site:dropbox.com; senza chiave la
             fonte viene saltata con un avviso
"""
import os, sys, json, time, argparse
import urllib.request, urllib.parse

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

USER_AGENT = "cerca-cosmo-skymed/1.0 (ricerca puntatori pubblici; no scraping anti-bot)"


def log(*a):
    print(*a, flush=True)


def dms2dec(d, m, s, h):
    v = float(d) + float(m) / 60 + float(s) / 3600
    return -v if str(h).upper() in ("S", "W") else v


def _http_get_json(url, headers=None, timeout=25):
    hh = {"User-Agent": USER_AGENT}
    if headers:
        hh.update(headers)
    req = urllib.request.Request(url, headers=hh)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


# ===================================================================== Zenodo
def cerca_zenodo(query, max_results=10):
    """Dataset/pubblicazioni open-access (API pubblica ufficiale, nessuna
    chiave richiesta): https://developers.zenodo.org/#list36. Se la rete o
    il formato risposta non sono disponibili ritorna [] (mai un'eccezione
    fatale: e' solo UNA delle fonti interrogate)."""
    try:
        url = "https://zenodo.org/api/records?" + urllib.parse.urlencode(
            {"q": query, "size": max_results, "sort": "mostrecent"})
        data = _http_get_json(url)
        out = []
        for hit in data.get("hits", {}).get("hits", []):
            md = hit.get("metadata", {})
            link = (hit.get("links", {}).get("self_html")
                    or hit.get("links", {}).get("html", ""))
            out.append({"fonte": "Zenodo", "titolo": md.get("title", "(senza titolo)"),
                        "url": link, "doi": md.get("doi", ""), "query": query})
        return out
    except Exception as ex:
        log(f"    ! ricerca Zenodo fallita per '{query}' ({ex})")
        return []


# =================================================================== OpenAIRE
def _openaire_titolo(meta):
    t = meta.get("title")
    if isinstance(t, list) and t:
        t = t[0]
    if isinstance(t, dict):
        return t.get("$", "(senza titolo)")
    return t or "(senza titolo)"


def _openaire_url(meta):
    """La struttura OpenAIRE cambia da record a record; provo prima il DOI
    (pid.$, il piu' affidabile e sempre risolvibile via doi.org), poi cerco
    ricorsivamente un URL nei metadati come fallback best-effort."""
    pid = meta.get("pid")
    if isinstance(pid, dict) and pid.get("@classid") == "doi" and pid.get("$"):
        return f"https://doi.org/{pid['$']}"
    if isinstance(pid, list):
        for p in pid:
            if isinstance(p, dict) and p.get("@classid") == "doi" and p.get("$"):
                return f"https://doi.org/{p['$']}"
    return _trova_primo_url(meta)


def _trova_primo_url(obj, profondita=0):
    if profondita > 6:
        return ""
    if isinstance(obj, str):
        return obj if obj.startswith("http") else ""
    if isinstance(obj, dict):
        for v in obj.values():
            u = _trova_primo_url(v, profondita + 1)
            if u:
                return u
    if isinstance(obj, list):
        for v in obj:
            u = _trova_primo_url(v, profondita + 1)
            if u:
                return u
    return ""


def cerca_openaire(query, max_results=10):
    """Pubblicazioni/dataset di ricerca europei (API pubblica ufficiale,
    nessuna chiave richiesta): https://api.openaire.eu. Copre molta
    letteratura/dataset che un motore generico indicizza poco o mette in
    fondo ai risultati."""
    try:
        url = "https://api.openaire.eu/search/publications?" + urllib.parse.urlencode(
            {"keywords": query, "size": max_results, "format": "json"})
        data = _http_get_json(url)
        risultati = (data.get("response", {}).get("results", {}) or {}).get("result", [])
        out = []
        for r in risultati:
            try:
                meta = r["metadata"]["oaf:entity"]["oaf:result"]
                titolo = _openaire_titolo(meta)
                url_pub = _openaire_url(meta)
            except Exception:
                titolo, url_pub = "(record OpenAIRE)", ""
            out.append({"fonte": "OpenAIRE", "titolo": titolo, "url": url_pub,
                        "doi": "", "query": query})
        return out
    except Exception as ex:
        log(f"    ! ricerca OpenAIRE fallita per '{query}' ({ex})")
        return []


# ==================================================== web search (opzionale)
def cerca_web_brave(query, api_key, max_results=10):
    """Ricerca sul web generico via l'API UFFICIALE Brave Search (chiave
    personale gratuita su https://api-dashboard.search.brave.com, va in
    BRAVE_SEARCH_KEY -- MAI hardcoded qui). E' il provider consigliato per la
    fonte 'web': l'API Bing Web Search e' stata RITIRATA da Microsoft
    (agosto 2025) e resta qui solo per chi ha ancora un accesso legacy.
    Supporta le query site:drive.google.com / site:dropbox.com / ecc."""
    if not api_key:
        return []
    try:
        url = "https://api.search.brave.com/res/v1/web/search?" + \
              urllib.parse.urlencode({"q": query, "count": max_results})
        data = _http_get_json(url, headers={"X-Subscription-Token": api_key,
                                            "Accept": "application/json"})
        out = []
        for w in (data.get("web", {}) or {}).get("results", []):
            out.append({"fonte": "web (Brave)", "titolo": w.get("title", ""),
                        "url": w.get("url", ""), "doi": "", "query": query})
        return out
    except Exception as ex:
        log(f"    ! ricerca web (Brave) fallita per '{query}' ({ex})")
        return []


def cerca_web_bing(query, api_key, max_results=10):
    """Ricerca sul web generico via l'API UFFICIALE Bing Web Search (serve
    una chiave personale in BING_SEARCH_KEY, MAI hardcoded qui). Usata anche
    per query mirate site:drive.google.com / site:github.com / ecc., dove
    capita che ricercatori condividano estratti di dataset. Nessun tentativo
    di scraping diretto dei motori di ricerca: quelli bloccano le richieste
    automatizzate senza API con una verifica anti-bot, che questo programma
    non tenta di eludere."""
    if not api_key:
        return []
    try:
        url = "https://api.bing.microsoft.com/v7.0/search?" + urllib.parse.urlencode(
            {"q": query, "count": max_results, "responseFilter": "Webpages"})
        data = _http_get_json(url, headers={"Ocp-Apim-Subscription-Key": api_key})
        out = []
        for w in data.get("webPages", {}).get("value", []):
            out.append({"fonte": "web (Bing)", "titolo": w.get("name", ""),
                        "url": w.get("url", ""), "doi": "", "query": query})
        return out
    except Exception as ex:
        log(f"    ! ricerca web (Bing) fallita per '{query}' ({ex})")
        return []


# ==================================================== nomi file COSMO-SkyMed
# Convenzione di denominazione dei PRODOTTI (nomi file) COSMO-SkyMed, dai
# Product Handbook ASI/e-GEOS:
#
# CSK -- 1a generazione (satelliti COSMO-SkyMed 1-4), es.:
#   CSKS1_SCS_B_HI_09_HH_RA_SF_20110105180325_20110105180332
#   \___/ \___/ \/ \/ \/ \/ \/ \____________/ \____________/
#   sat.  tipo  modo|  pol |  consegna inizio (UTC)  fine (UTC)
#                  swath  lato/orbita (R/L + A/D)
#   - satellite: CSKS1..CSKS4
#   - tipo prodotto: RAW_B, SCS_B, SCS_U (Single-look Complex Slant),
#     DGM_B (Detected Ground Multi-look), GEC_B (Geocoded Ellipsoid
#     Corrected), GTC_B (Geocoded Terrain Corrected)
#   - modo: HI (HIMAGE/Stripmap), PP (PingPong), WR (WideRegion ScanSAR),
#     HR (HugeRegion ScanSAR), S2 (Spotlight-2)
#   - polarizzazione: HH, VV, HV, VH (CO/CH/CV per PingPong)
#
# CSG -- 2a generazione (COSMO-SkyMed Second Generation), es.:
#   CSG_SSAR1_DGM_B_0301_STR_008_HH_RA_D_20201118050553_20201118050600
#   - satellite: CSG_SSAR1, CSG_SSAR2
#   - tipi prodotto analoghi (RAW_B/SCS_B/SCS_U/DGM_B/GEC_B/GTC_B)
#   - modi: SPOTLIGHT-2A/B/C, STR (Stripmap), PNG (PingPong), QPS (QuadPol),
#     SC1/SC2 (ScanSAR)
#
# Ogni nome file CSK/CSG INIZIA col token satellite: una query col solo
# prefisso satellite ha gia' il richiamo massimo (trova QUALUNQUE nome file di
# quel satellite); le combinazioni satellite_tipo affinano la ricerca su
# pagine/dataset che riportano i nomi troncati o tokenizzati. Modo/pol/date
# sono sottostringhe del nome completo: aggiungerle non troverebbe pagine in
# piu' (solo sottoinsiemi di quelle gia' trovate dai prefissi), quindi non
# vengono enumerate nelle query.
SATELLITI_CSK = ["CSKS1", "CSKS2", "CSKS3", "CSKS4"]
SATELLITI_CSG = ["CSG_SSAR1", "CSG_SSAR2"]
TIPI_PRODOTTO = ["RAW_B", "SCS_B", "SCS_U", "DGM_B", "GEC_B", "GTC_B"]
MODI_CSK = ["HI", "PP", "WR", "HR", "S2"]
POLARIZZAZIONI = ["HH", "VV", "HV", "VH"]


def query_nomi_file(livello="tutti"):
    """Query costruite sui NOMI FILE dei prodotti COSMO-SkyMed (convenzione
    sopra): un dataset/una pagina che elenca scene CSK/CSG contiene questi
    token anche quando non nomina mai 'COSMO-SkyMed' per esteso.
      - 'nessuno': nessuna query sui nomi file (solo le query keyword);
      - 'base'   : prefissi satellite (CSKS1..4, CSG_SSAR1/2);
      - 'tutti'  : base + TUTTE le combinazioni satellite_tipo (default).
    Ogni pattern e' racchiuso tra virgolette (ricerca a FRASE ESATTA).
    I tipi prodotto DA SOLI ('SCS_B', 'RAW_B', ...) non vengono interrogati:
    testati su Zenodo/OpenAIRE restituiscono quasi solo rumore (il full-text
    dei cataloghi li pesca dentro allegati/OCR estranei), mentre ancorati al
    satellite (CSKS1_SCS_B, ...) restano precisi."""
    if livello == "nessuno":
        return []
    q = list(SATELLITI_CSK) + list(SATELLITI_CSG)
    if livello == "tutti":
        q += [f"{s}_{t}" for s in SATELLITI_CSK + SATELLITI_CSG
              for t in TIPI_PRODOTTO]
    return [f'"{t}"' for t in q]


# ============================================================== query builder
# Drive pubblici in internet dove capita che dataset/estratti vengano condivisi
# apertamente: SEMPRE inclusi nelle query site: della fonte 'web' (Bing).
SITI_DRIVE = ["drive.google.com", "dropbox.com"]
# Altri siti di condivisione, aggiunti solo con --siti-condivisione.
SITI_CONDIVISIONE = ["onedrive.live.com", "github.com", "researchgate.net",
                     "figshare.com"]


def costruisci_query(keywords, aoi_hint=None, con_siti_condivisione=False,
                     nomi_file_web=()):
    """Costruisce le query: (base, web_extra). Le query base (+ variante con
    suggerimento AOI) vanno a TUTTE le fonti; le varianti mirate site: sui
    drive/siti di condivisione (drive.google.com e dropbox.com SEMPRE, gli
    altri con --siti-condivisione) hanno senso solo per un motore di ricerca
    e vanno SOLO alla fonte 'web', non alle API accademiche. Su ogni sito
    site: si cerca sia la query keyword sia i pattern in `nomi_file_web`
    (prefissi satellite CSKS1..4/CSG_SSAR1-2: ogni nome file CSK/CSG inizia
    col satellite, quindi hanno il richiamo massimo sui file condivisi)."""
    base = f"COSMO-SkyMed {keywords}".strip()
    q = [base]
    if aoi_hint:
        q.append(f"{base} {aoi_hint}")
    siti = list(SITI_DRIVE) + (SITI_CONDIVISIONE if con_siti_condivisione else [])
    q_web = [f"site:{sito} {qq}" for sito in siti
             for qq in [base, *nomi_file_web]]
    return q, q_web


def main():
    ap = argparse.ArgumentParser(
        description="Cerca puntatori pubblici a materiale che cita/riguarda "
                    "prodotti COSMO-SkyMed, su cataloghi accademici e (con "
                    "chiave opzionale) sul web generico. NON scarica scene "
                    "SAR e NON bypassa alcun controllo di accesso.")
    ap.add_argument("--keywords", default="pyramid Giza",
                    help="parole chiave aggiuntive (default: 'pyramid Giza')")
    ap.add_argument("--nw", nargs=4, metavar=("D", "M", "S", "H"), default=None,
                    help="angolo NW lat (DMS), opzionale: aggiunge un suggerimento AOI alla query")
    ap.add_argument("--nw-lon", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--se", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--se-lon", nargs=4, metavar=("D", "M", "S", "H"), default=None)
    ap.add_argument("--nomi-file", choices=["tutti", "base", "nessuno"],
                    default="tutti",
                    help="query aggiuntive sui NOMI FILE dei prodotti "
                         "COSMO-SkyMed (convenzione CSKSn_/CSG_SSARn_...): "
                         "'tutti' (default) = prefissi satellite + TUTTE le "
                         "combinazioni satellite_tipo (es. CSKS1_SCS_B); "
                         "'base' = solo i prefissi satellite; 'nessuno' = "
                         "solo le query keyword")
    ap.add_argument("--sources", default="zenodo,openaire,web",
                    help="fonti da interrogare separate da virgola tra "
                         "zenodo,openaire,web (default: tutte e tre; 'web' "
                         "esegue anche le query site:drive.google.com / "
                         "site:dropbox.com e richiede BRAVE_SEARCH_KEY o "
                         "BING_SEARCH_KEY, altrimenti viene saltata)")
    ap.add_argument("--max-results", type=int, default=10,
                    help="risultati massimi per query per fonte")
    ap.add_argument("--pause", type=float, default=1.0,
                    help="pausa in secondi tra una query/fonte e la successiva "
                         "(rispetto dei servizi pubblici interrogati)")
    ap.add_argument("--siti-condivisione", action="store_true",
                    help="con --sources web: oltre a site:drive.google.com e "
                         "site:dropbox.com (SEMPRE inclusi), aggiunge anche "
                         "site:onedrive.live.com / site:github.com / "
                         "site:researchgate.net / site:figshare.com")
    ap.add_argument("--brave-key", default=None,
                    help="chiave Brave Search API (altrimenti env BRAVE_SEARCH_KEY; "
                         "mai hardcoded nel codice). Provider consigliato per la "
                         "fonte 'web': chiave gratuita su "
                         "https://api-dashboard.search.brave.com")
    ap.add_argument("--bing-key", default=None,
                    help="chiave Bing Web Search (altrimenti env BING_SEARCH_KEY; "
                         "mai hardcoded nel codice). NB: API ritirata da Microsoft "
                         "nell'agosto 2025, resta per chi ha un accesso legacy")
    ap.add_argument("--outdir", default=".",
                    help="cartella dove salvare il report JSON (default: cartella corrente)")
    a = ap.parse_args()

    aoi_hint = None
    if a.nw and a.nw_lon and a.se and a.se_lon:
        lat_nw = dms2dec(*a.nw[:3], a.nw[3]); lon_nw = dms2dec(*a.nw_lon[:3], a.nw_lon[3])
        lat_se = dms2dec(*a.se[:3], a.se[3]); lon_se = dms2dec(*a.se_lon[:3], a.se_lon[3])
        aoi_hint = f"{(lat_nw + lat_se) / 2:.4f}N {(lon_nw + lon_se) / 2:.4f}E"

    log("=" * 74)
    log("cerca_cosmo_skymed.py -- ricerca puntatori pubblici (NON un downloader)")
    log("=" * 74)
    log("COSMO-SkyMed e' un prodotto ASI/e-GEOS ad accesso ristretto: questo")
    log("programma NON scarica dati e NON bypassa alcun controllo di accesso,")
    log("cerca solo pagine/dataset GIA' pubblici che lo citano, tramite API")
    log("ufficiali. Per l'accesso vero alle scene vedi https://www.asi.it")
    log("(programma 'Open Call') o https://www.e-geos.it.")
    log("=" * 74)

    bing_key = a.bing_key or os.environ.get("BING_SEARCH_KEY")
    brave_key = a.brave_key or os.environ.get("BRAVE_SEARCH_KEY")
    sorgenti = {s.strip().lower() for s in a.sources.split(",") if s.strip()}
    if "web" in sorgenti and not (brave_key or bing_key):
        log("[cerca] fonte 'web' (incluse le query site:drive.google.com / "
            "site:dropbox.com) saltata: nessuna chiave in BRAVE_SEARCH_KEY o "
            "BING_SEARCH_KEY (nessuno scraping di motori di ricerca senza API: "
            "verrebbe bloccato da una verifica anti-bot, e questo programma "
            "non tenta di eluderla). Chiave Brave gratuita su "
            "https://api-dashboard.search.brave.com")
        sorgenti.discard("web")
    web_nome = "Brave" if brave_key else "Bing"

    def cerca_web(q):
        if brave_key:
            return cerca_web_brave(q, brave_key, a.max_results)
        return cerca_web_bing(q, bing_key, a.max_results)

    # site: sui drive: keyword base + prefissi satellite dei nomi file
    nomi_web = ([f'"{s}"' for s in SATELLITI_CSK + SATELLITI_CSG]
                if a.nomi_file != "nessuno" else [])
    query_base, query_web = costruisci_query(
        a.keywords, aoi_hint, con_siti_condivisione=a.siti_condivisione,
        nomi_file_web=nomi_web)
    query_nomi = query_nomi_file(a.nomi_file)
    query_base += query_nomi
    n_query = len(query_base) + (len(query_web) if "web" in sorgenti else 0)
    log(f"[cerca] {n_query} query, fonti: {sorted(sorgenti) or ['(nessuna)']}")
    if query_nomi:
        log(f"[cerca] di cui {len(query_nomi)} sui nomi file dei prodotti "
            f"COSMO-SkyMed (--nomi-file {a.nomi_file}: CSKS1..4, CSG_SSAR1/2"
            f"{' e tutte le combinazioni satellite_tipo' if a.nomi_file == 'tutti' else ''})")
    if "web" in sorgenti:
        log(f"[cerca] query site: sui drive/siti di condivisione (solo fonte web): "
            f"{len(query_web)}")

    tutti = []
    for q in query_base:
        log(f"[cerca] '{q}'")
        if "zenodo" in sorgenti:
            r = cerca_zenodo(q, a.max_results)
            log(f"    zenodo: {len(r)} risultati"); tutti += r
            time.sleep(a.pause)
        if "openaire" in sorgenti:
            r = cerca_openaire(q, a.max_results)
            log(f"    openaire: {len(r)} risultati"); tutti += r
            time.sleep(a.pause)
        if "web" in sorgenti:
            r = cerca_web(q)
            log(f"    web ({web_nome}): {len(r)} risultati"); tutti += r
            time.sleep(a.pause)
    # query site:drive.google.com / site:dropbox.com / ... -> solo fonte web
    if "web" in sorgenti:
        for q in query_web:
            log(f"[cerca] '{q}'")
            r = cerca_web(q)
            log(f"    web ({web_nome}): {len(r)} risultati"); tutti += r
            time.sleep(a.pause)

    visti, dedup = set(), []
    for r in tutti:
        chiave = r.get("url") or (r.get("fonte"), r.get("titolo"))
        if chiave in visti:
            continue
        visti.add(chiave); dedup.append(r)

    os.makedirs(a.outdir, exist_ok=True)
    outpath = os.path.join(a.outdir, "cosmo_skymed_risultati.json")
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(dedup, f, indent=2, ensure_ascii=False)

    log("=" * 74)
    log(f"[cerca] {len(dedup)} risultati unici (su {len(tutti)} totali) -> {outpath}")
    log("=" * 74)
    for r in dedup:
        log(f"  [{r['fonte']}] {r['titolo']}")
        if r.get("url"):
            log(f"      {r['url']}")
    log("=" * 74)
    log("Ricontrolla sempre licenza/provenienza prima di riusare quanto trovato.")
    log("Per l'accesso autorizzato alle scene COSMO-SkyMed: https://www.asi.it "
        "(COSMO-SkyMed Open Call) / https://www.e-geos.it")


if __name__ == "__main__":
    main()

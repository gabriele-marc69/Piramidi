# -*- coding: utf-8 -*-
"""
scarica_alta_risoluzione.py — ricerca e download dei prodotti SAR a piu' alta
risoluzione disponibili su Copernicus Data Space (CDSE) per la tomografia
micro-Doppler delle piramidi.

PERCHE' QUESTO PROGRAMMA (requisiti dati dal paper arXiv 2208.00811v1,
Biondi & Malanga, "SAR Doppler Tomography ... Great Pyramid of Giza"):

  Il paper lavora su dati COSMO-SkyMed Second Generation (CSG), banda X,
  acquisizione Stripmap/Staring-Spotlight, in formato *SLC* (Single Look
  Complex).  I file necessari al metodo sono quindi:

    1. SLC (complessi, con fase) — indispensabili: le sub-aperture Doppler
       si sintetizzano dallo spettro complesso; i prodotti GRD (solo
       ampiezza, detected) NON vanno bene.
    2. Massima banda Doppler / risoluzione azimutale possibile — la
       risoluzione tomografica e' delta_z = lambda*R/(2A): piu' banda
       azimutale => tomogramma piu' fine.
    3. Singola polarizzazione co-polare (HH o VV) ad alta PRF.

  I dati CSG Spotlight del paper non sono liberi.  Sul catalogo CDSE
  (gratuito) l'equivalente migliore, in ordine di preferenza, e':

    a) Sentinel-1 *Stripmap* SLC  (S1_SLC__1S .. S6_SLC__1S):
       ~1.7-3.6 m (range) x ~4.3-4.9 m (azimut) — il piu' vicino al setup
       del paper; acquisito solo su siti selezionati.
    b) se lo Stripmap SLC non esiste a catalogo ma esiste il grezzo
       Stripmap RAW L0 (S?_RAW__0S): si *ordina la produzione on-demand*
       dell'SLC con la On-Demand Production API
       (https://documentation.dataspace.copernicus.eu/APIs/On-Demand%20Production%20API.html)
       e si scarica il risultato.
    c) Sentinel-1 IW SLC (~2.3 x 14.1 m in azimut): il fallback gia' usato
       dalla pipeline (piramide_unificato.py).

  Se il download HTTP dal CDSE fallisce del tutto, come ultima spiaggia il
  programma cerca il prodotto per nome su Academic Torrents (dataset
  scientifici distribuiti via BitTorrent) e lo scarica con libtorrent.

USO:
  python scarica_alta_risoluzione.py --pyramid cheope --dry-run
  python scarica_alta_risoluzione.py --coords "29.9792,31.1342" --download
  python scarica_alta_risoluzione.py --coords-file punti.txt --download
  python scarica_alta_risoluzione.py --torrent-search "S1A IW SLC 20250321"
  python scarica_alta_risoluzione.py --torrent-search "..." --download
  (credenziali: --cdse-user/--cdse-pass oppure env CDSE_USER / CDSE_PASS)

  punti.txt: una coppia "lat lon" (o "lat,lon") per riga; l'AOI e' il
  rettangolo che racchiude tutti i punti (+ margine --margin, default 200 m).
"""

import argparse
import json
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

CATALOGO = "https://catalogue.dataspace.copernicus.eu/odata/v1"
ODP = "https://odp.dataspace.copernicus.eu/odata/v1"
IDENTITY = ("https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
            "protocol/openid-connect/token")

# preset coerenti con piramide_unificato.py (DMS -> decimale)
PIRAMIDI = {
    "cheope": dict(nw=("29", "58", "48.0", "N"), nw_lon=("31", "7", "58.4", "E"),
                   se=("29", "58", "41.0", "N"), se_lon=("31", "8", "07.7", "E"),
                   label="Cheope (Khnum-Khufu)"),
    "kefren": dict(nw=("29", "58", "38.0", "N"), nw_lon=("31", "7", "45.4", "E"),
                   se=("29", "58", "29.0", "N"), se_lon=("31", "7", "55.4", "E"),
                   label="Chefren (Khafre)"),
}

# tipi di prodotto in ordine di preferenza (vedi intestazione):
#   modo SM = beam S1..S6 nel nome SAFE (es. S1A_S3_SLC__1SSH_...)
TIPI_SM_SLC = [f"S{i}_SLC__1S" for i in range(1, 7)]
TIPI_SM_RAW = [f"S{i}_RAW__0S" for i in range(1, 7)]
TIPO_IW_SLC = "IW_SLC__1S"


def log(*a):
    print(*a, flush=True)


def dms2dec(d, m, s, h):
    v = float(d) + float(m) / 60 + float(s) / 3600
    return -v if str(h).upper() in ("S", "W") else v


# ================================================================ COORDINATE
def leggi_coordinate(args):
    """Ritorna l'AOI (lonW, lonE, latS, latN) dalle opzioni --pyramid,
    --coords o --coords-file. Un punto singolo diventa un rettangolo con
    margine --margin metri."""
    punti = []
    if args.pyramid:
        p = PIRAMIDI[args.pyramid]
        latN = dms2dec(*p["nw"]); lonW = dms2dec(*p["nw_lon"])
        latS = dms2dec(*p["se"]); lonE = dms2dec(*p["se_lon"])
        log(f"[coord] preset {p['label']}: lon {lonW:.5f}..{lonE:.5f}, "
            f"lat {latS:.5f}..{latN:.5f}")
        return lonW, lonE, latS, latN
    if args.coords:
        lat, lon = [float(x) for x in re.split(r"[,\s]+", args.coords.strip())]
        punti.append((lat, lon))
    if args.coords_file:
        with open(args.coords_file, encoding="utf-8") as f:
            for riga in f:
                riga = riga.strip()
                if not riga or riga.startswith("#"):
                    continue
                lat, lon = [float(x) for x in re.split(r"[,\s]+", riga)[:2]]
                punti.append((lat, lon))
    if not punti:
        raise SystemExit("specifica --pyramid, --coords o --coords-file")
    lats = [p[0] for p in punti]; lons = [p[1] for p in punti]
    # margine in gradi (~1 deg lat = 111 km); sul lon correggo per la latitudine
    import math
    mlat = args.margin / 111_000.0
    mlon = args.margin / (111_000.0 * max(0.1, math.cos(math.radians(sum(lats) / len(lats)))))
    aoi = (min(lons) - mlon, max(lons) + mlon, min(lats) - mlat, max(lats) + mlat)
    log(f"[coord] {len(punti)} punto/i -> AOI lon {aoi[0]:.5f}..{aoi[1]:.5f}, "
        f"lat {aoi[2]:.5f}..{aoi[3]:.5f} (margine {args.margin} m)")
    return aoi


# ================================================================== CATALOGO
def _get_json(url, token=None, timeout=90):
    # User-Agent esplicito: alcuni servizi (es. Academic Torrents) rifiutano
    # con 403 lo UA di default di urllib
    headers = {"User-Agent": "piramid-v2/1.0 (+https://dataspace.copernicus.eu)"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def cerca_catalogo(aoi, tipi, start, end, top=100):
    """Cerca su CDSE i prodotti Sentinel-1 dei tipi dati che intersecano
    l'AOI. Ritorna [{Name, Id, ContentLength, Start}] ordinati per data."""
    lonW, lonE, latS, latN = aoi
    poly = (f"POLYGON(({lonW} {latS},{lonE} {latS},{lonE} {latN},"
            f"{lonW} {latN},{lonW} {latS}))")
    per_tipo = " or ".join(f"contains(Name,'_{t}')" for t in tipi)
    flt = ("Collection/Name eq 'SENTINEL-1' "
           f"and OData.CSC.Intersects(area=geography'SRID=4326;{poly}') "
           f"and ({per_tipo}) "
           f"and ContentDate/Start gt {start}T00:00:00.000Z "
           f"and ContentDate/Start lt {end}T23:59:59.000Z")
    url = CATALOGO + "/Products?" + urllib.parse.urlencode(
        {"$filter": flt, "$orderby": "ContentDate/Start desc", "$top": str(top)})
    data = _get_json(url)
    return [{"Name": f["Name"], "Id": f["Id"],
             "ContentLength": f.get("ContentLength", 0),
             "Start": f.get("ContentDate", {}).get("Start", "")}
            for f in data.get("value", [])]


# ============================================================ AUTENTICAZIONE
def cdse_token(user, pwd):
    body = urllib.parse.urlencode({"client_id": "cdse-public",
                                   "grant_type": "password",
                                   "username": user, "password": pwd}).encode()
    with urllib.request.urlopen(urllib.request.Request(IDENTITY, data=body),
                                timeout=60) as r:
        return json.load(r)["access_token"]


def _rete_ok():
    try:
        socket.getaddrinfo("catalogue.dataspace.copernicus.eu", 443)
        return True
    except OSError:
        return False


def _attendi_rete(max_wait=7200, check_every=30):
    if _rete_ok():
        return True
    log(f"    ! rete assente: attendo (controllo ogni {check_every}s, "
        f"max {max_wait // 60} min)...")
    t0 = time.time()
    while time.time() - t0 < max_wait:
        time.sleep(check_every)
        if _rete_ok():
            log(f"    rete tornata dopo {time.time() - t0:.0f}s")
            return True
    return False


class _AuthRedirect(urllib.request.HTTPRedirectHandler):
    """Conserva l'Authorization sui redirect cross-host del CDSE (altrimenti 401)."""
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None:
            auth = req.get_header("Authorization")
            if auth:
                new.add_unredirected_header("Authorization", auth)
        return new


_OPENER = urllib.request.build_opener(_AuthRedirect())


def scarica_ripristinabile(url, dest, get_token, expected=0, n_try=60):
    """Download con ripresa (header Range) e token fresco a ogni tentativo —
    stessa logica collaudata di piramide_unificato.py."""
    part = dest + ".part"
    for att in range(1, n_try + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if expected and have >= expected:
            break
        try:
            headers = {"Authorization": f"Bearer {get_token()}"}
            if have:
                headers["Range"] = f"bytes={have}-"
                log(f"    riprendo da {have / 1e9:.2f} GB (tentativo {att}/{n_try})")
            req = urllib.request.Request(url, headers=headers)
            with _OPENER.open(req, timeout=1800) as r:
                mode = "ab" if (have and r.status == 206) else "wb"
                if mode == "wb":
                    have = 0
                with open(part, mode) as fo:
                    while True:
                        chunk = r.read(1 << 20)
                        if not chunk:
                            break
                        fo.write(chunk)
                        have += len(chunk)
            if not expected or have >= expected:
                break
            log(f"    ! ricevuti {have / 1e9:.2f}/{expected / 1e9:.1f} GB; ritento")
        except Exception as ex:
            wait = min(60, 5 * att)
            log(f"    ! interrotto ({ex}); ritento tra {wait}s ({att}/{n_try})")
            time.sleep(wait)
            if not _attendi_rete():
                return False
    else:
        log(f"    ! download non completato dopo {n_try} tentativi")
        return False
    if expected and os.path.getsize(part) < expected:
        log("    ! download incompleto")
        return False
    os.replace(part, dest)
    return True


# ==================================================== ON-DEMAND PRODUCTION API
def odp_workflows(token):
    """Elenco dei workflow di produzione on-demand disponibili."""
    data = _get_json(ODP + "/Workflows", token=token)
    return data.get("value", [])


def odp_trova_workflow(token, input_tipo, output_tipo):
    """Cerca il workflow che trasforma `input_tipo` (es. S3_RAW__0S) in
    `output_tipo` (es. S3_SLC__1S). I campi esposti dall'API sono
    InputProductType / OutputProductType."""
    for wf in odp_workflows(token):
        win = (wf.get("InputProductType") or "").upper()
        wout = (wf.get("OutputProductType") or "").upper()
        if input_tipo.upper() in win and output_tipo.upper() in wout:
            return wf
    return None


def odp_ordina(token, workflow_name, prodotto_nome, priorita=1):
    """Sottomette un ordine di produzione on-demand. Ritorna l'Id ordine."""
    body = json.dumps({
        "WorkflowName": workflow_name,
        "InputProductReference": {"Reference": prodotto_nome},
        "Priority": priorita,
        "Name": f"odp_{prodotto_nome[:40]}",
    }).encode()
    req = urllib.request.Request(
        ODP + "/ProductionOrder/OData.CSC.Order", data=body,
        headers={"Authorization": f"Bearer {token}",
                 "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        risp = json.load(r)
    oid = risp.get("Id") or risp.get("value", {}).get("Id")
    log(f"    ordine sottomesso: Id={oid}")
    return oid


def odp_attendi_e_scarica(get_token, order_id, dest, poll_s=120, max_h=12):
    """Polla lo stato dell'ordine finche' completed, poi scarica il risultato.
    I risultati restano in staging 14 giorni: l'ordine non va ripetuto se il
    programma viene interrotto e rilanciato."""
    t0 = time.time()
    while time.time() - t0 < max_h * 3600:
        try:
            o = _get_json(f"{ODP}/ProductionOrders({order_id})",
                          token=get_token())
            stato = (o.get("Status") or "").lower()
            log(f"    stato ordine {order_id}: {stato}"
                + (f" (stima: {o['EstimatedDate']})" if o.get("EstimatedDate") else ""))
            if stato == "completed":
                url = f"{ODP}/ProductionOrder({order_id})/Product/$value"
                return scarica_ripristinabile(url, dest, get_token)
            if stato in ("failed", "cancelled"):
                log(f"    ! ordine terminato con stato '{stato}': "
                    f"{o.get('StatusMessage', '')}")
                return False
        except Exception as ex:
            log(f"    ! polling fallito ({ex}); riprovo")
            _attendi_rete()
        time.sleep(poll_s)
    log(f"    ! ordine non completato entro {max_h} h")
    return False


# ============================================================ FALLBACK TORRENT
# Ricerca nella rete torrent: il DHT di BitTorrent e' indicizzato per infohash,
# non per parola chiave, quindi la "ricerca nella rete" passa dai motori che
# crawlano il DHT / dagli indici pubblici. Si interrogano piu' provider e si
# uniscono i risultati; il download poi avviene davvero in rete (DHT + tracker
# aperti), senza dipendere dal sito che ha fornito l'infohash.

_TRACKER_APERTI = [
    "udp://tracker.opentrackr.org:1337/announce",
    "udp://open.demonii.com:1337/announce",
    "udp://tracker.torrent.eu.org:451/announce",
    "udp://exodus.desync.com:6969/announce",
]


def _magnet(infohash, nome):
    tr = "".join(f"&tr={urllib.parse.quote(t)}" for t in _TRACKER_APERTI)
    return f"magnet:?xt=urn:btih:{infohash}&dn={urllib.parse.quote(nome)}{tr}"


def _prov_academic(q):
    """Academic Torrents (dataset scientifici)."""
    data = _get_json("https://academictorrents.com/apiv2/entries?" +
                     urllib.parse.urlencode({"search": q, "limit": "20"}),
                     timeout=60)
    out = []
    for v in data.get("entries") or data.get("value") or []:
        ih = v.get("infohash") or v.get("infoHash")
        if ih:
            out.append((v.get("name", q), ih, int(v.get("seeders", 0) or 0)))
    return out


def _prov_torrents_csv(q):
    """Torrents-CSV (indice aperto, API JSON senza chiave)."""
    data = _get_json("https://torrents-csv.com/service/search?" +
                     urllib.parse.urlencode({"q": q, "size": "20"}), timeout=60)
    return [(t.get("name", q), t["infohash"], int(t.get("seeders", 0) or 0))
            for t in data.get("torrents", []) if t.get("infohash")]


def _prov_apibay(q):
    """apibay (indice generalista; utile perche' molti mirror di dati aperti
    finiscono li'). Scarta la riga-sentinella 'No results returned'."""
    data = _get_json("https://apibay.org/q.php?" +
                     urllib.parse.urlencode({"q": q}), timeout=60)
    out = []
    for t in data if isinstance(data, list) else []:
        ih = (t.get("info_hash") or "").strip()
        if ih and set(ih) != {"0"}:
            out.append((t.get("name", q), ih, int(t.get("seeders", 0) or 0)))
    return out


_PROVIDER = [("Academic Torrents", _prov_academic),
             ("Torrents-CSV", _prov_torrents_csv),
             ("apibay", _prov_apibay)]


def torrent_cerca(query, min_token=2):
    """Interroga tutti gli indici della rete torrent per `query` e ritorna la
    lista [(nome, magnet, seeders)] dei candidati pertinenti, ordinata per
    numero di seeder. La pertinenza richiede che almeno `min_token` token
    significativi della query compaiano nel nome del torrent (evita di
    scaricare un file qualsiasi che l'indice ritorna per match parziale)."""
    token = [t.lower() for t in re.split(r"[^0-9A-Za-z]+", query) if len(t) >= 3]
    soglia = min(min_token, len(token)) or 1
    candidati = {}
    for nome_prov, prov in _PROVIDER:
        try:
            voci = prov(query)
        except Exception as ex:
            log(f"    ! {nome_prov}: ricerca fallita ({ex})")
            continue
        log(f"    {nome_prov}: {len(voci)} risultati grezzi")
        for nome, ih, seed in voci:
            punteggio = sum(1 for t in token if t in nome.lower())
            if punteggio >= soglia:
                ih = ih.lower()
                if ih not in candidati or seed > candidati[ih][2]:
                    candidati[ih] = (nome, _magnet(ih, nome), seed)
    out = sorted(candidati.values(), key=lambda c: -c[2])
    for nome, _, seed in out[:10]:
        log(f"    candidato: {nome}  ({seed} seeder)")
    return out


def torrent_scarica(magnet, destdir, timeout_h=24):
    """Scarica un magnet URI con libtorrent. Ritorna True se completato."""
    try:
        import libtorrent as lt
    except ImportError:
        log("    ! libtorrent non installato: pip install libtorrent  "
            "(oppure python-libtorrent dal gestore pacchetti)")
        return False
    os.makedirs(destdir, exist_ok=True)
    # DHT attivo con i router di bootstrap: i metadati e i peer si trovano
    # nella rete anche se i tracker del magnet non rispondono
    ses = lt.session({
        "listen_interfaces": "0.0.0.0:6881",
        "enable_dht": True,
        "dht_bootstrap_nodes": ("router.bittorrent.com:6881,"
                                "dht.transmissionbt.com:6881,"
                                "router.utorrent.com:6881"),
    })
    # API libtorrent 2.x: parse_magnet_uri + add_torrent
    # (add_magnet_uri e' deprecata)
    params = lt.parse_magnet_uri(magnet)
    params.save_path = destdir
    h = ses.add_torrent(params)
    log("    recupero metadati torrent (DHT + tracker)...")
    t0 = time.time()
    while not h.status().has_metadata:
        if time.time() - t0 > 600:
            log("    ! nessun metadato in 10 min (nessun seeder?): rinuncio")
            return False
        time.sleep(2)
    log(f"    scarico '{h.status().name}' via BitTorrent...")
    while time.time() - t0 < timeout_h * 3600:
        s = h.status()
        if s.is_seeding:
            log("    torrent completato")
            return True
        log(f"    {s.progress * 100:.1f}%  ({s.download_rate / 1e6:.2f} MB/s, "
            f"{s.num_peers} peer)")
        time.sleep(10)
    log(f"    ! torrent non completato entro {timeout_h} h")
    return False


def _query_varianti(nome_safe):
    """Varianti di ricerca per un prodotto SAFE, dalla piu' specifica alla
    piu' generica: nome completo, poi missione+modo+data (i mirror raramente
    indicizzano il nome SAFE integrale)."""
    nome = nome_safe.replace(".SAFE", "")
    varianti = [nome]
    m = re.match(r"(S1[A-D])_(\w{2})_(\w{3})__\w+_(\d{8})T", nome)
    if m:
        mis, modo, tipo, data = m.groups()
        varianti += [f"{mis} {modo} {tipo} {data}", f"{mis} {tipo} {data}",
                     f"sentinel-1 {tipo} {data}"]
    return varianti


def fallback_torrent(prodotti, outdir, magnet=None):
    """Ultima spiaggia quando il CDSE non e' scaricabile: cerca i prodotti
    negli indici della rete torrent e scarica via DHT/tracker con libtorrent.
    Con --magnet si salta la ricerca e si scarica direttamente quel link."""
    log("[torrent] fallback BitTorrent (ricerca multi-indice + libtorrent/DHT)")
    if magnet:
        return torrent_scarica(magnet, outdir)
    ok = False
    for p in prodotti:
        trovato = False
        for q in _query_varianti(p["Name"]):
            log(f"[torrent] cerco '{q}'...")
            candidati = torrent_cerca(q)
            for nome, mg, seed in candidati:
                log(f"[torrent] provo '{nome}' ({seed} seeder)")
                if torrent_scarica(mg, outdir):
                    ok = trovato = True
                    break
            if trovato:
                break
        if not trovato:
            log(f"    nessun torrent utilizzabile per {p['Name']}")
    if not ok:
        log("[torrent] nessun prodotto recuperato via torrent. I dati "
            "Sentinel non hanno mirror torrent ufficiali: riprovare il "
            "download CDSE piu' tardi resta la via principale.")
    return ok


# ======================================================================= MAIN
def main():
    ap = argparse.ArgumentParser(
        description="Cerca e scarica i prodotti SAR a piu' alta risoluzione "
                    "(SLC Stripmap > on-demand da RAW L0 > SLC IW) per l'AOI dato.")
    ap.add_argument("--pyramid", choices=sorted(PIRAMIDI), help="preset AOI")
    ap.add_argument("--coords", help='punto "lat,lon" (decimale)')
    ap.add_argument("--coords-file", help='file con righe "lat lon"')
    ap.add_argument("--margin", type=float, default=200.0,
                    help="margine AOI in metri attorno ai punti (default 200)")
    ap.add_argument("--start", default="2014-10-01", help="data inizio ricerca")
    ap.add_argument("--end", default=time.strftime("%Y-%m-%d"), help="data fine")
    ap.add_argument("--max-products", type=int, default=12,
                    help="quanti prodotti scaricare al massimo (default 12)")
    ap.add_argument("--outdir", default="alta_risoluzione_out",
                    help="cartella di destinazione")
    ap.add_argument("--download", action="store_true",
                    help="scarica davvero (servono credenziali CDSE)")
    ap.add_argument("--dry-run", action="store_true",
                    help="solo ricerca a catalogo, nessun download/ordine")
    ap.add_argument("--cdse-user", default=None)
    ap.add_argument("--cdse-pass", default=None)
    ap.add_argument("--magnet", default=None,
                    help="magnet URI da scaricare direttamente col fallback torrent")
    ap.add_argument("--torrent-search", default=None, metavar="QUERY",
                    help="cerca QUERY negli indici della rete torrent, stampa i "
                         "magnet trovati ed esce (con --download scarica il migliore)")
    args = ap.parse_args()

    # ---- modalita' ricerca torrent pura (non serve il catalogo CDSE) -------
    if args.torrent_search:
        log(f"[torrent] ricerca '{args.torrent_search}' negli indici...")
        candidati = torrent_cerca(args.torrent_search)
        if not candidati:
            log("[torrent] nessun risultato pertinente")
            return
        for nome, mg, seed in candidati:
            log(f"\n  {nome}  ({seed} seeder)\n  {mg}")
        if args.download:
            os.makedirs(args.outdir, exist_ok=True)
            nome, mg, _ = candidati[0]
            log(f"\n[torrent] scarico il candidato migliore: {nome}")
            torrent_scarica(mg, args.outdir)
        return

    aoi = leggi_coordinate(args)
    os.makedirs(args.outdir, exist_ok=True)

    # ---- 1) ricerca a catalogo, in ordine di qualita' ----------------------
    log("[1] cerco Sentinel-1 STRIPMAP SLC (risoluzione massima su CDSE)...")
    sm_slc = cerca_catalogo(aoi, TIPI_SM_SLC, args.start, args.end)
    log(f"    trovati {len(sm_slc)} prodotti SM SLC")
    log("[1] cerco Sentinel-1 STRIPMAP RAW L0 (per produzione on-demand)...")
    sm_raw = cerca_catalogo(aoi, TIPI_SM_RAW, args.start, args.end)
    log(f"    trovati {len(sm_raw)} prodotti SM RAW L0")
    log("[1] cerco Sentinel-1 IW SLC (fallback, gia' usato dalla pipeline)...")
    iw_slc = cerca_catalogo(aoi, [TIPO_IW_SLC], args.start, args.end)
    log(f"    trovati {len(iw_slc)} prodotti IW SLC")

    piano = []                       # (etichetta, prodotto, via)
    for p in sm_slc[:args.max_products]:
        piano.append(("SM_SLC diretto", p, "diretto"))
    if not sm_slc:
        for p in sm_raw[:args.max_products]:
            piano.append(("SM_RAW -> SLC on-demand", p, "odp"))
    if not sm_slc and not sm_raw:
        for p in iw_slc[:args.max_products]:
            piano.append(("IW_SLC diretto (fallback)", p, "diretto"))

    with open(os.path.join(args.outdir, "piano_download.json"), "w",
              encoding="utf-8") as f:
        json.dump([{"via": v, "tipo": e, **p} for e, p, v in piano], f, indent=2)
    log(f"[1] piano: {len(piano)} prodotti "
        f"(salvato in {args.outdir}/piano_download.json)")
    for e, p, _ in piano:
        log(f"      [{e}] {p['Name']}  ({(p['ContentLength'] or 0) / 1e9:.1f} GB)")

    if args.dry_run or not args.download:
        if not args.dry_run:
            log("[fine] usa --download per scaricare (credenziali CDSE richieste)")
        return

    # ---- 2) credenziali ----------------------------------------------------
    user = args.cdse_user or os.environ.get("CDSE_USER")
    pwd = args.cdse_pass or os.environ.get("CDSE_PASS")
    if not (user and pwd):
        raise SystemExit("credenziali mancanti: --cdse-user/--cdse-pass o env "
                         "CDSE_USER/CDSE_PASS (account gratuito su "
                         "dataspace.copernicus.eu)")
    try:
        cdse_token(user, pwd)
    except Exception as ex:
        raise SystemExit(f"login CDSE fallito: {ex}")
    get_token = lambda: cdse_token(user, pwd)

    # ---- 3) esecuzione del piano -------------------------------------------
    falliti = []
    for etichetta, p, via in piano:
        dest = os.path.join(args.outdir, p["Name"].replace(".SAFE", "") + ".zip")
        if os.path.exists(dest):
            log(f"[3] {p['Name']}: gia' scaricato, salto")
            continue
        log(f"[3] [{etichetta}] {p['Name']}")
        ok = False
        try:
            if via == "diretto":
                url = f"{CATALOGO}/Products({p['Id']})/$value"
                ok = scarica_ripristinabile(url, dest, get_token,
                                            expected=int(p.get("ContentLength") or 0))
            else:                                    # produzione on-demand
                # il beam e' nel nome (es. _S3_RAW__0S) -> workflow S3_RAW->S3_SLC
                m = re.search(r"_(S\d)_RAW__0S", p["Name"])
                beam = m.group(1) if m else "S1"
                wf = odp_trova_workflow(get_token(), f"{beam}_RAW__0S",
                                        f"{beam}_SLC__1S")
                if wf is None:
                    # alcuni workflow dichiarano tipi generici: ripiega sul primo
                    # workflow S1 L0->L1 SLC disponibile
                    wf = next((w for w in odp_workflows(get_token())
                               if "RAW" in (w.get("InputProductType") or "")
                               and "SLC" in (w.get("OutputProductType") or "")),
                              None)
                if wf is None:
                    log("    ! nessun workflow on-demand RAW->SLC disponibile")
                else:
                    log(f"    workflow: {wf['Name']} "
                        f"({wf.get('InputProductType')} -> {wf.get('OutputProductType')})")
                    oid = odp_ordina(get_token(), wf["Name"], p["Name"])
                    ok = odp_attendi_e_scarica(get_token, oid, dest)
        except urllib.error.HTTPError as ex:
            log(f"    ! HTTP {ex.code}: {ex.reason}")
        except Exception as ex:
            log(f"    ! errore: {ex}")
        if not ok:
            falliti.append(p)

    # ---- 4) fallback torrent per cio' che non e' arrivato via HTTP ---------
    if falliti or args.magnet:
        fallback_torrent(falliti, args.outdir, magnet=args.magnet)

    n_ok = len(piano) - len(falliti)
    log(f"[fine] scaricati {n_ok}/{len(piano)} prodotti in {args.outdir}")


if __name__ == "__main__":
    main()

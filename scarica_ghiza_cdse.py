#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scarico dei prodotti Sentinel-1 sulla piana di Giza dal Copernicus Data
Space Ecosystem (CDSE).

Interroga il catalogo OData di CDSE per i prodotti che coprono le coordinate
delle tre piramidi, a partire da una data di inizio (default 2026-01-01), e li
salva in una cartella locale.

Due modalita' di scarico:

  --mode full   il prodotto .SAFE completo, come archivio .zip (~8.3 GB per
                data). E' cio' che si ottiene dal portale web.

  --mode swath  solo il sottoinsieme che serve alla tomografia di questo
                repository -- manifest.safe, l'annotation .xml e il
                measurement .tiff di UNA sola combinazione swath/polarizzazione
                (default iw2/vv) -- ricostruito come albero .SAFE sul disco
                (~1.5 GB per data, cioe' un sesto). Usa l'API Nodes di CDSE.

Autenticazione: password grant sul realm Keycloak "CDSE"; le credenziali sono
lette da .cdse.env (CDSE_USER / CDSE_PASS), sovrascrivibili dalle variabili
d'ambiente omonime. L'access token dura 1800 s e viene rinnovato da solo
durante gli scarichi lunghi.

Esempi
------
    # cosa c'e' da scaricare, senza scaricare nulla
    python scarica_ghiza_cdse.py --list-only

    # pila coerente (stessa orbita relativa), solo iw2/vv
    python scarica_ghiza_cdse.py --mode swath --relative-orbit 58

    # prodotti completi, primi 10
    python scarica_ghiza_cdse.py --mode full --max-products 10
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
import zipfile
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional, Sequence, Tuple

import requests

# ==========================================================================
# Costanti di servizio
# ==========================================================================

IDENTITY_URL = ("https://identity.dataspace.copernicus.eu/auth/realms/CDSE"
                "/protocol/openid-connect/token")
CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1"
DOWNLOAD_URL = "https://download.dataspace.copernicus.eu/odata/v1"
OIDC_CLIENT_ID = "cdse-public"

#: Vertici delle piramidi (WGS84). Sono le stesse coordinate usate da
#: piramidi_v01/v02, cosi' l'area interrogata e' esattamente quella analizzata.
PIRAMIDI: Tuple[Tuple[str, float, float], ...] = (
    ("Cheope (Khufu)",      29.979235, 31.134202),
    ("Chefren (Khafre)",    29.976111, 31.130833),
    ("Micerino (Menkaure)", 29.972500, 31.128056),
)

DEFAULT_OUT = r"e:\claude_code\Biondi_Malanga\DATA_Ghiza"
DEFAULT_ENV = ".cdse.env"

#: dimensione dei blocchi letti dal socket
CHUNK = 1 << 20                       # 1 MiB

#: le righe di avanzamento si riscrivono in posto solo su un terminale: in un
#: file di log lascerebbero una riga sola, illeggibile
TTY = sys.stdout.isatty()

#: un prodotto IW SLC contiene 3 swath x 2 polarizzazioni; scaricandone uno
#: solo si prende circa questa frazione del totale (misurata: 1.47 / 8.32 GB)
FRAZIONE_SWATH = 5.6


@dataclass
class Config:
    """Parametri dello scarico."""

    out_dir: str = DEFAULT_OUT
    env_file: str = DEFAULT_ENV

    # --- finestra temporale ------------------------------------------------
    start: str = "2026-01-01"
    end: Optional[str] = None          # None = adesso

    # --- selezione del prodotto -------------------------------------------
    collection: str = "SENTINEL-1"
    product_type: str = "IW_SLC__1S"
    platform: Optional[str] = None     # "A" | "C" | "D" ...
    relative_orbit: Optional[int] = None
    orbit_direction: Optional[str] = None
    margin_deg: float = 0.008          # ~900 m attorno all'inviluppo dei vertici

    # --- modalita' di scarico ---------------------------------------------
    mode: str = "swath"                # "swath" | "zip" | "full"
    swath: str = "iw2"
    polarisation: str = "vh"

    # --- limiti e sicurezza ------------------------------------------------
    max_products: int = 0              # 0 = nessun limite
    reserve_gb: float = 10.0           # spazio da lasciare libero sul disco
    verify_md5: bool = False
    retries: int = 10
    timeout_connect: float = 30.0
    timeout_read: float = 300.0


# ==========================================================================
# 1.  Credenziali e token
# ==========================================================================

def carica_credenziali(env_file: str) -> Tuple[str, str]:
    """Legge CDSE_USER / CDSE_PASS dal file .env, con l'ambiente che vince.

    Il file e' nel formato KEY=VALUE, una riga per chiave; le righe vuote e
    quelle che iniziano con '#' sono ignorate."""
    valori: Dict[str, str] = {}
    if os.path.exists(env_file):
        with open(env_file, "r", encoding="utf-8") as fh:
            for riga in fh:
                riga = riga.strip()
                if not riga or riga.startswith("#") or "=" not in riga:
                    continue
                k, v = riga.split("=", 1)
                valori[k.strip()] = v.strip().strip('"').strip("'")

    user = os.environ.get("CDSE_USER") or valori.get("CDSE_USER", "")
    pwd = os.environ.get("CDSE_PASS") or valori.get("CDSE_PASS", "")
    if not user or not pwd:
        raise SystemExit(
            f"credenziali assenti: servono CDSE_USER e CDSE_PASS in {env_file} "
            "oppure nelle variabili d'ambiente")
    return user, pwd


class TokenCDSE:
    """Token di accesso con rinnovo automatico.

    L'access token dura 1800 s e il refresh token 3600 s: uno scarico da 8 GB
    puo' superarli entrambi, quindi si rinnova con il refresh finche' e' valido
    e si rifa' il login da capo quando scade anche quello."""

    def __init__(self, user: str, password: str, session: requests.Session):
        self._user = user
        self._password = password
        self._s = session
        self._access = ""
        self._refresh = ""
        self._scade = 0.0               # epoch di scadenza dell'access token
        self._scade_refresh = 0.0
        self._login()

    def _richiedi(self, dati: Dict[str, str]) -> None:
        r = self._s.post(IDENTITY_URL, data=dati, timeout=60)
        if not r.ok:
            raise SystemExit(f"autenticazione CDSE fallita ({r.status_code}): "
                             f"{r.text[:300]}")
        j = r.json()
        ora = time.time()
        self._access = j["access_token"]
        self._refresh = j.get("refresh_token", "")
        # 60 s di margine per non usare un token che scade a meta' richiesta
        self._scade = ora + float(j.get("expires_in", 600)) - 60.0
        self._scade_refresh = ora + float(j.get("refresh_expires_in", 3600)) - 60.0

    def _login(self) -> None:
        self._richiedi({"client_id": OIDC_CLIENT_ID, "grant_type": "password",
                        "username": self._user, "password": self._password})

    def valore(self) -> str:
        """Il bearer token, rinnovato se necessario."""
        ora = time.time()
        if ora < self._scade:
            return self._access
        if self._refresh and ora < self._scade_refresh:
            try:
                self._richiedi({"client_id": OIDC_CLIENT_ID,
                                "grant_type": "refresh_token",
                                "refresh_token": self._refresh})
                return self._access
            except SystemExit:
                pass                    # refresh rifiutato: si rifa' il login
        self._login()
        return self._access

    def header(self) -> Dict[str, str]:
        return {"Authorization": f"Bearer {self.valore()}"}


# ==========================================================================
# 2.  Interrogazione del catalogo
# ==========================================================================

def _iso(giorno: str, fine: bool = False) -> str:
    """Normalizza una data 'YYYY-MM-DD' (o gia' ISO) nel formato OData."""
    giorno = giorno.strip().replace("/", "-")
    if "T" in giorno:
        return giorno if giorno.endswith("Z") else giorno + "Z"
    ora = "23:59:59.999" if fine else "00:00:00.000"
    return f"{giorno}T{ora}Z"


def poligono_aoi(margin_deg: float) -> str:
    """Rettangolo WGS84 che racchiude tutte e tre le piramidi, piu' un margine.

    Un solo POINT basterebbe per Cheope, ma i tre vertici distano ~800 m e la
    pipeline li analizza insieme: si chiede quindi che il prodotto intersechi
    l'inviluppo di tutti e tre."""
    lat = [p[1] for p in PIRAMIDI]
    lon = [p[2] for p in PIRAMIDI]
    s, n = min(lat) - margin_deg, max(lat) + margin_deg
    w, e = min(lon) - margin_deg, max(lon) + margin_deg
    return (f"POLYGON(({w:.6f} {s:.6f},{e:.6f} {s:.6f},{e:.6f} {n:.6f},"
            f"{w:.6f} {n:.6f},{w:.6f} {s:.6f}))")


def costruisci_filtro(cfg: Config) -> str:
    """Il $filter OData completo."""
    poly = poligono_aoi(cfg.margin_deg)
    parti = [
        f"Collection/Name eq '{cfg.collection}'",
        f"OData.CSC.Intersects(area=geography'SRID=4326;{poly}')",
        f"ContentDate/Start gt {_iso(cfg.start)}",
    ]
    if cfg.end:
        parti.append(f"ContentDate/Start lt {_iso(cfg.end, fine=True)}")
    if cfg.product_type:
        parti.append("Attributes/OData.CSC.StringAttribute/any("
                     f"a:a/Name eq 'productType' and a/Value eq '{cfg.product_type}')")
    if cfg.relative_orbit is not None:
        parti.append("Attributes/OData.CSC.IntegerAttribute/any("
                     "a:a/Name eq 'relativeOrbitNumber' and "
                     f"a/Value eq {cfg.relative_orbit})")
    if cfg.orbit_direction:
        parti.append("Attributes/OData.CSC.StringAttribute/any("
                     "a:a/Name eq 'orbitDirection' and "
                     f"a/Value eq '{cfg.orbit_direction.upper()}')")
    if cfg.platform:
        parti.append("Attributes/OData.CSC.StringAttribute/any("
                     "a:a/Name eq 'platformSerialIdentifier' and "
                     f"a/Value eq '{cfg.platform.upper()}')")
    return " and ".join(parti)


def cerca_prodotti(cfg: Config, session: requests.Session,
                   verbose: bool = True) -> List[Dict[str, Any]]:
    """Tutti i prodotti che soddisfano il filtro, seguendo la paginazione.

    Il catalogo e' pubblico: questa fase non richiede il token."""
    params: Dict[str, Any] = {
        "$filter": costruisci_filtro(cfg),
        "$orderby": "ContentDate/Start asc",
        "$expand": "Attributes",
        "$top": 200,
        "$count": "true",
    }
    url: Optional[str] = f"{CATALOGUE_URL}/Products"
    out: List[Dict[str, Any]] = []
    atteso: Optional[int] = None

    while url:
        r = _get_con_retry(session, url, cfg, params=params, timeout_read=180.0)
        with r:
            j = r.json()
        if atteso is None:
            atteso = j.get("@odata.count")
        out.extend(j.get("value", []))
        url = j.get("@odata.nextLink")
        params = {}                     # il nextLink porta gia' la query
        if verbose and atteso:
            _avanzamento(f"      catalogo: {len(out)}/{atteso} prodotti")
    if verbose:
        _fine_avanzamento()
    return out


def attributo(prod: Dict[str, Any], nome: str) -> Any:
    """Valore di un attributo esteso del prodotto, o None."""
    for a in prod.get("Attributes", []):
        if a.get("Name") == nome:
            return a.get("Value")
    return None


def checksum_md5(prod: Dict[str, Any]) -> Optional[str]:
    """L'MD5 dichiarato dal catalogo, se presente."""
    for c in prod.get("Checksum", []):
        if str(c.get("Algorithm", "")).upper() == "MD5":
            val = c.get("Value")
            return str(val) if val else None
    return None


# ==========================================================================
# 3.  HTTP con ritentativi
# ==========================================================================

def _get_con_retry(session: requests.Session, url: str, cfg: Config,
                   params: Optional[Dict[str, Any]] = None,
                   headers: Optional[Dict[str, str]] = None,
                   stream: bool = False,
                   timeout_read: Optional[float] = None) -> requests.Response:
    """GET con backoff esponenziale su errori di rete, 429 e 5xx.

    Il servizio chiude ogni tanto la connessione a meta' trasferimento (si
    osservano errori SSL sui prodotti da 8 GB): sul volume di dati in gioco il
    ritentativo e' la norma, non l'eccezione."""
    attesa = 5.0
    ultimo = ""
    for tentativo in range(1, cfg.retries + 1):
        try:
            r = session.get(url, params=params, headers=headers, stream=stream,
                            timeout=(cfg.timeout_connect,
                                     timeout_read or cfg.timeout_read),
                            allow_redirects=True)
            if r.status_code in (429, 500, 502, 503, 504):
                ultimo = f"HTTP {r.status_code}"
                r.close()
            else:
                r.raise_for_status()
                return r
        except requests.HTTPError as ex:
            raise RuntimeError(f"{url}: {ex}") from ex
        except requests.RequestException as ex:
            ultimo = f"{type(ex).__name__}: {ex}"
        if tentativo < cfg.retries:
            print(f"      ritento fra {attesa:.0f} s ({ultimo})")
            time.sleep(attesa)
            attesa = min(attesa * 2.0, 120.0)
    raise RuntimeError(f"{url}: fallito dopo {cfg.retries} tentativi ({ultimo})")


# ==========================================================================
# 4.  Scarico di un flusso su file, con ripresa
# ==========================================================================

def _avanzamento(riga: str) -> None:
    """Riga di avanzamento riscritta in posto sul terminale, saltata altrove."""
    if TTY:
        print(riga, end="\r")


def _fine_avanzamento() -> None:
    """Ripulisce l'ultima riga di avanzamento."""
    if TTY:
        print(" " * 110, end="\r")


def _fmt(byte: float) -> str:
    """Dimensione leggibile."""
    for u in ("B", "KB", "MB", "GB", "TB"):
        if abs(byte) < 1024.0:
            return f"{byte:.1f} {u}"
        byte /= 1024.0
    return f"{byte:.1f} PB"


def _apri_flusso(url: str, gia: int, cfg: Config, session: requests.Session,
                 token: TokenCDSE) -> requests.Response:
    """Apre il flusso, chiedendo la ripresa da `gia` se non e' zero."""
    headers = dict(token.header())
    if gia > 0:
        headers["Range"] = f"bytes={gia}-"
    r = session.get(url, headers=headers, stream=True,
                    timeout=(cfg.timeout_connect, cfg.timeout_read),
                    allow_redirects=True)
    r.raise_for_status()
    return r


def scarica_flusso(url: str, destinazione: str, atteso: int, cfg: Config,
                   session: requests.Session, token: TokenCDSE,
                   etichetta: str) -> int:
    """Scarica `url` in `destinazione`, riprendendo un .part gia' presente.

    Restituisce i byte scritti in totale. Il file viene costruito in un
    '<nome>.part' e rinominato solo a trasferimento completo, cosi' la presenza
    del file finale e' garanzia che sia integro.

    Il ritentativo avvolge TUTTO il trasferimento, non solo l'apertura della
    connessione: su un measurement da 1.5 GB il servizio chiude regolarmente il
    socket a meta' strada (SSLError, ChunkedEncodingError) e riaprire e basta
    non serve a nulla. A ogni tentativo si riparte dai byte gia' sul disco con
    un Range; se il server risponde 200 invece di 206 la ripresa non e'
    concessa e si ricomincia da zero -- meglio ritrasferire che incollare due
    tronconi disallineati."""
    parziale = destinazione + ".part"
    attesa = 5.0
    ultimo = ""

    for tentativo in range(1, cfg.retries + 1):
        gia = os.path.getsize(parziale) if os.path.exists(parziale) else 0
        if atteso and gia == atteso:
            os.replace(parziale, destinazione)
            return gia
        try:
            r = _apri_flusso(url, gia, cfg, session, token)
        except requests.RequestException as ex:
            ultimo = f"{type(ex).__name__}: {ex}"
            if tentativo < cfg.retries:
                print(f"      ritento fra {attesa:.0f} s ({ultimo})")
                time.sleep(attesa)
                attesa = min(attesa * 2.0, 120.0)
            continue

        with r:
            if gia > 0 and r.status_code != 206:
                print("      ripresa non concessa dal server, riparto da 0")
                gia = 0
            lunghezza = r.headers.get("Content-Length")
            totale = gia + int(lunghezza) if lunghezza else atteso
            if gia:
                print(f"      riprendo da {_fmt(gia)} / {_fmt(totale)}")

            scritti = gia
            t0 = time.time()
            ultimo_avviso = 0.0
            try:
                with open(parziale, "ab" if gia > 0 else "wb") as fh:
                    for blocco in r.iter_content(CHUNK):
                        if not blocco:
                            continue
                        fh.write(blocco)
                        scritti += len(blocco)
                        ora = time.time()
                        # su terminale ogni 5 s, su file di log ogni 60 s
                        if ora - ultimo_avviso >= (5.0 if TTY else 60.0):
                            dt = max(ora - t0, 1e-6)
                            vel = (scritti - gia) / dt
                            pct = 100.0 * scritti / totale if totale else 0.0
                            resta = (totale - scritti) / vel if vel > 0 else 0.0
                            riga = (f"      {etichetta}: {pct:5.1f} %  "
                                    f"{_fmt(scritti)} / {_fmt(totale)}  "
                                    f"{_fmt(vel)}/s  mancano "
                                    f"{timedelta(seconds=int(resta))}   ")
                            if TTY:
                                print(riga, end="\r")
                            else:
                                print(riga, flush=True)
                            ultimo_avviso = ora
            except requests.RequestException as ex:
                _fine_avanzamento()
                ultimo = f"{type(ex).__name__}: {ex}"
                if tentativo < cfg.retries:
                    print(f"      interrotto a {_fmt(scritti)} / {_fmt(totale)}"
                          f"; ritento fra {attesa:.0f} s ({ultimo})")
                    time.sleep(attesa)
                    attesa = min(attesa * 2.0, 120.0)
                continue

        _fine_avanzamento()
        if totale and scritti != totale:
            ultimo = f"trasferimento incompleto: {scritti} / {totale} byte"
            if tentativo < cfg.retries:
                print(f"      {ultimo}; ritento fra {attesa:.0f} s")
                time.sleep(attesa)
                attesa = min(attesa * 2.0, 120.0)
            continue
        os.replace(parziale, destinazione)
        return scritti

    raise RuntimeError(f"{os.path.basename(destinazione)}: fallito dopo "
                       f"{cfg.retries} tentativi ({ultimo})")


def md5_file(path: str) -> str:
    """MD5 di un file letto a blocchi."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for blocco in iter(lambda: fh.read(CHUNK), b""):
            h.update(blocco)
    return h.hexdigest()


# ==========================================================================
# 5.  Modalita' "full": il prodotto .SAFE completo come .zip
# ==========================================================================

def scarica_completo(prod: Dict[str, Any], cfg: Config,
                     session: requests.Session,
                     token: TokenCDSE) -> Dict[str, Any]:
    """Il prodotto intero, cosi' come lo serve /Products(id)/$value."""
    pid = prod["Id"]
    nome_zip = str(prod["Name"]).replace(".SAFE", "") + ".zip"
    dest = os.path.join(cfg.out_dir, nome_zip)
    atteso = int(prod.get("ContentLength", 0))

    if os.path.exists(dest) and os.path.getsize(dest) == atteso:
        print("      gia' presente e completo, salto")
        return {"stato": "presente", "file": dest, "byte": atteso}

    url = f"{DOWNLOAD_URL}/Products({pid})/$value"
    byte = scarica_flusso(url, dest, atteso, cfg, session, token, nome_zip[:34])

    rec: Dict[str, Any] = {"stato": "scaricato", "file": dest, "byte": byte}
    if cfg.verify_md5:
        atteso_md5 = checksum_md5(prod)
        if atteso_md5:
            _avanzamento("      verifica MD5...")
            ottenuto = md5_file(dest)
            rec["md5_atteso"] = atteso_md5
            rec["md5_ottenuto"] = ottenuto
            rec["md5_ok"] = (ottenuto == atteso_md5)
            esito = "ok" if rec["md5_ok"] else "NON CORRISPONDE"
            print(f"      MD5 {esito}                    ")
    return rec


# ==========================================================================
# 6.  Modalita' "swath": solo iw2/vv, via l'API Nodes
# ==========================================================================

def _nodi(session: requests.Session, cfg: Config, token: TokenCDSE,
          url: str) -> List[Dict[str, Any]]:
    """Elenco dei nodi figli di un percorso dentro il prodotto."""
    r = _get_con_retry(session, url, cfg, headers=token.header(),
                       timeout_read=120.0)
    with r:
        return list(r.json().get("result", []))


def scarica_swath(prod: Dict[str, Any], cfg: Config,
                  session: requests.Session,
                  token: TokenCDSE) -> Dict[str, Any]:
    """Ricostruisce in locale un .SAFE ridotto al solo swath/pol richiesto.

    Sul disco si ottiene

        <nome>.SAFE/manifest.safe
        <nome>.SAFE/annotation/<...>-<swath>-slc-<pol>-<...>.xml
        <nome>.SAFE/annotation/calibration/{calibration,noise}-<...>.xml
        <nome>.SAFE/measurement/<...>-<swath>-slc-<pol>-<...>.tiff

    che e' esattamente cio' che parse_annotation()/read_window() leggono. Un
    prodotto intero pesa ~8.3 GB, questo sottoinsieme ~1.5 GB."""
    pid = prod["Id"]
    nome = str(prod["Name"])
    radice = os.path.join(cfg.out_dir, nome)
    base = f"{DOWNLOAD_URL}/Products({pid})/Nodes({nome})"
    marca = f"-{cfg.swath}-slc-{cfg.polarisation}-".lower()

    os.makedirs(os.path.join(radice, "annotation", "calibration"), exist_ok=True)
    os.makedirs(os.path.join(radice, "measurement"), exist_ok=True)

    da_prendere: List[Tuple[str, str, int]] = []   # (url, destinazione, byte)

    # manifest.safe -- contiene la lista delle componenti del prodotto
    for n in _nodi(session, cfg, token, f"{base}/Nodes"):
        if n.get("Name") == "manifest.safe":
            da_prendere.append((f"{base}/Nodes(manifest.safe)/$value",
                                os.path.join(radice, "manifest.safe"),
                                int(n.get("ContentLength", 0))))

    # annotation: lo .xml dello swath scelto, piu' calibration e noise
    ann = f"{base}/Nodes(annotation)"
    for n in _nodi(session, cfg, token, f"{ann}/Nodes"):
        nm = str(n.get("Name", ""))
        if nm.lower().endswith(".xml") and marca in nm.lower():
            da_prendere.append((f"{ann}/Nodes({nm})/$value",
                                os.path.join(radice, "annotation", nm),
                                int(n.get("ContentLength", 0))))
        elif nm == "calibration":
            cal = f"{ann}/Nodes(calibration)"
            for c in _nodi(session, cfg, token, f"{cal}/Nodes"):
                cn = str(c.get("Name", ""))
                if marca in cn.lower():
                    da_prendere.append(
                        (f"{cal}/Nodes({cn})/$value",
                         os.path.join(radice, "annotation", "calibration", cn),
                         int(c.get("ContentLength", 0))))

    # measurement: il solo .tiff dello swath scelto
    mea = f"{base}/Nodes(measurement)"
    for n in _nodi(session, cfg, token, f"{mea}/Nodes"):
        nm = str(n.get("Name", ""))
        if marca in nm.lower():
            da_prendere.append((f"{mea}/Nodes({nm})/$value",
                                os.path.join(radice, "measurement", nm),
                                int(n.get("ContentLength", 0))))

    if not any(d.endswith(".tiff") for _, d, _ in da_prendere):
        raise RuntimeError(f"nessun measurement {cfg.swath}/{cfg.polarisation} "
                           f"in {nome}")

    totale = 0
    for url, dest, atteso in da_prendere:
        if os.path.exists(dest) and (not atteso or os.path.getsize(dest) == atteso):
            totale += os.path.getsize(dest)
            continue
        totale += scarica_flusso(url, dest, atteso, cfg, session, token,
                                 os.path.basename(dest)[:34])
    return {"stato": "scaricato", "file": radice, "byte": totale,
            "componenti": len(da_prendere)}


# ==========================================================================
# 6b. Modalita' "zip": scarica il prodotto intero, ne estrae il solo
#     swath/polarizzazione utile e cancella l'archivio
# ==========================================================================

def membri_utili(zf: zipfile.ZipFile, cfg: Config) -> List[str]:
    """I soli membri dell'archivio che servono alla pipeline.

    Sono il measurement .tiff dello swath/polarizzazione scelto, il suo
    annotation .xml, i due .xml di calibration e noise corrispondenti, e il
    manifest.safe del prodotto."""
    marca = f"-{cfg.swath}-slc-{cfg.polarisation}-".lower()
    out: List[str] = []
    for nome in zf.namelist():
        if nome.endswith("/"):
            continue
        basso = nome.lower()
        if basso.endswith("/manifest.safe"):
            out.append(nome)
        elif marca in basso and (basso.endswith(".tiff") or basso.endswith(".xml")):
            out.append(nome)
    return out


def scarica_ed_estrai(prod: Dict[str, Any], cfg: Config,
                      session: requests.Session,
                      token: TokenCDSE) -> Dict[str, Any]:
    """Scarica il .zip completo, ne estrae i soli file utili, poi lo cancella.

    Un prodotto alla volta: l'archivio da ~7.7 GB esiste sul disco solo per il
    tempo dell'estrazione, dopo restano i ~1.5 GB che servono davvero. Il .zip
    va comunque trasferito per intero -- l'indice di uno zip sta in fondo --
    quindi si paga tutta la banda; l'alternativa senza questo costo e'
    --mode swath, che chiede al servizio i singoli file."""
    pid = prod["Id"]
    nome = str(prod["Name"])
    radice = os.path.join(cfg.out_dir, nome)
    zip_path = os.path.join(cfg.out_dir, nome.replace(".SAFE", "") + ".zip")
    atteso = int(prod.get("ContentLength", 0))

    # gia' estratto in una passata precedente?
    mea = os.path.join(radice, "measurement")
    if os.path.isdir(mea) and any(f.endswith(".tiff") for f in os.listdir(mea)):
        print("      gia' estratto, salto")
        return {"stato": "presente", "file": radice,
                "byte": _dimensione_albero(radice)}

    url = f"{DOWNLOAD_URL}/Products({pid})/$value"
    scarica_flusso(url, zip_path, atteso, cfg, session, token,
                   os.path.basename(zip_path)[:34])

    try:
        with zipfile.ZipFile(zip_path) as zf:
            membri = membri_utili(zf, cfg)
            if not any(m.lower().endswith(".tiff") for m in membri):
                raise RuntimeError(
                    f"nessun measurement {cfg.swath}/{cfg.polarisation} "
                    f"nell'archivio di {nome}")
            print(f"      estraggo {len(membri)} file utili su "
                  f"{len(zf.namelist())}")
            for m in membri:
                zf.extract(m, cfg.out_dir)
    finally:
        # l'archivio se ne va comunque: se l'estrazione e' fallita non serve
        # tenersi 7.7 GB, il prodotto si riscarica
        if os.path.exists(zip_path):
            os.remove(zip_path)
            print(f"      cancellato {os.path.basename(zip_path)} "
                  f"({_fmt(atteso)})")

    tenuti = _dimensione_albero(radice)
    print(f"      tenuti {_fmt(tenuti)} su {_fmt(atteso)} scaricati")
    return {"stato": "scaricato", "file": radice, "byte": tenuti,
            "byte_trasferiti": atteso, "componenti": len(membri)}


def _dimensione_albero(radice: str) -> int:
    """Byte occupati da una cartella, ricorsivamente."""
    tot = 0
    for cartella, _dirs, files in os.walk(radice):
        for f in files:
            try:
                tot += os.path.getsize(os.path.join(cartella, f))
            except OSError:
                pass
    return tot


# ==========================================================================
# 7.  Programma principale
# ==========================================================================

def stampa_elenco(prodotti: Sequence[Dict[str, Any]]) -> None:
    """Tabella di cio' che il catalogo ha restituito."""
    print(f"\n  {'data':<12} {'sat':<4} {'orb.rel':>7} {'passo':<11} "
          f"{'dim':>9}  nome")
    print("  " + "-" * 104)
    for p in prodotti:
        nome = str(p["Name"])
        inizio = str(p.get("ContentDate", {}).get("Start", ""))[:10]
        dim = int(p.get("ContentLength", 0))
        print(f"  {inizio:<12} {nome[:3]:<4} "
              f"{str(attributo(p, 'relativeOrbitNumber')):>7} "
              f"{str(attributo(p, 'orbitDirection') or '')[:10]:<11} "
              f"{_fmt(dim):>9}  {nome}")


def riepilogo_orbite(prodotti: Sequence[Dict[str, Any]]) -> None:
    """Quante acquisizioni per orbita relativa: una pila tomografica coerente
    deve stare tutta sulla stessa, altrimenti le baseline non hanno senso."""
    conteggio: Dict[Tuple[Any, Any], int] = {}
    for p in prodotti:
        k = (attributo(p, "relativeOrbitNumber"), attributo(p, "orbitDirection"))
        conteggio[k] = conteggio.get(k, 0) + 1
    print("\n  orbite relative presenti (una pila coerente sta su UNA sola):")
    for (orb, direz), n in sorted(conteggio.items(), key=lambda x: -x[1]):
        print(f"    orbita relativa {orb} {direz}: {n} acquisizioni")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Scarica i prodotti Sentinel-1 su Giza da Copernicus (CDSE).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT, help="cartella di destinazione")
    ap.add_argument("--env-file", default=DEFAULT_ENV,
                    help="file con CDSE_USER / CDSE_PASS")
    ap.add_argument("--start", default="2026-01-01",
                    help="data di inizio (YYYY-MM-DD)")
    ap.add_argument("--end", default=None, help="data di fine (default: adesso)")
    ap.add_argument("--collection", default="SENTINEL-1")
    ap.add_argument("--product-type", default="IW_SLC__1S",
                    help="stringa vuota per non filtrare sul tipo di prodotto")
    ap.add_argument("--platform", default=None,
                    help="identificativo della piattaforma: A, C, D ...")
    ap.add_argument("--relative-orbit", type=int, default=None,
                    help="numero di orbita relativa (per una pila coerente)")
    ap.add_argument("--orbit-direction", default=None,
                    choices=["ASCENDING", "DESCENDING"])
    ap.add_argument("--margin-deg", type=float, default=0.008,
                    help="margine attorno all'inviluppo dei tre vertici")
    ap.add_argument("--mode", default="swath",
                    choices=["swath", "zip", "full"],
                    help="swath: chiede al servizio i soli file utili; "
                         "zip: scarica il prodotto intero, estrae i file utili "
                         "e cancella l'archivio; full: tiene il .zip intero")
    ap.add_argument("--swath", default="iw2", help="non usato con --mode full")
    ap.add_argument("--polarisation", default="vh",
                    help="non usato con --mode full")
    ap.add_argument("--max-products", type=int, default=0,
                    help="0 = tutti quelli trovati")
    ap.add_argument("--reserve-gb", type=float, default=10.0,
                    help="spazio da lasciare libero sul disco")
    ap.add_argument("--verify-md5", action="store_true",
                    help="ricalcola l'MD5 dello .zip e lo confronta col catalogo")
    ap.add_argument("--retries", type=int, default=10)
    ap.add_argument("--list-only", action="store_true",
                    help="mostra cosa c'e' da scaricare e termina")
    a = ap.parse_args(argv)

    cfg = Config(out_dir=a.out, env_file=a.env_file, start=a.start, end=a.end,
                 collection=a.collection, product_type=a.product_type,
                 platform=a.platform, relative_orbit=a.relative_orbit,
                 orbit_direction=a.orbit_direction, margin_deg=a.margin_deg,
                 mode=a.mode, swath=a.swath, polarisation=a.polarisation,
                 max_products=a.max_products, reserve_gb=a.reserve_gb,
                 verify_md5=a.verify_md5, retries=a.retries)

    print("=" * 78)
    print("  Sentinel-1 su Giza -- Copernicus Data Space Ecosystem")
    print("=" * 78)
    for nome, lat, lon in PIRAMIDI:
        print(f"    {nome:<22} {lat:.6f} N  {lon:.6f} E")
    print(f"\n    area interrogata : {poligono_aoi(cfg.margin_deg)}")
    print(f"    finestra         : da {cfg.start} a {cfg.end or 'adesso'}")
    print(f"    prodotto         : {cfg.collection} / "
          f"{cfg.product_type or 'qualsiasi'}")
    dettaglio = ("" if cfg.mode == "full"
                 else f"  (si tiene solo {cfg.swath}/{cfg.polarisation})")
    print(f"    modalita'        : {cfg.mode}{dettaglio}")
    print(f"    destinazione     : {cfg.out_dir}")

    session = requests.Session()
    session.headers["User-Agent"] = "giza-cdse-downloader/1.0"

    print("\n  [1] interrogazione del catalogo")
    prodotti = cerca_prodotti(cfg, session)
    if not prodotti:
        print("      nessun prodotto trovato per questi filtri.")
        return 1
    dim_totale = sum(int(p.get("ContentLength", 0)) for p in prodotti)
    print(f"      {len(prodotti)} prodotti, {_fmt(dim_totale)} complessivi")
    stampa_elenco(prodotti)
    riepilogo_orbite(prodotti)

    if cfg.max_products > 0:
        prodotti = list(prodotti[:cfg.max_products])
        print(f"\n  limitato ai primi {len(prodotti)} prodotti (--max-products)")

    # stima dell'occupazione a valle della modalita' scelta
    dim_trasferita = sum(int(p.get("ContentLength", 0)) for p in prodotti)
    dim_scelta = dim_trasferita
    if cfg.mode in ("swath", "zip"):
        # a regime resta solo lo swath/pol utile; in modalita' zip la banda
        # spesa resta pero' quella del prodotto intero
        dim_scelta = int(dim_trasferita / FRAZIONE_SWATH)
        if cfg.mode == "zip":
            # picco: l'archivio piu' grande convive con la sua estrazione
            piu_grande = max(int(p.get("ContentLength", 0)) for p in prodotti)
            dim_scelta += piu_grande

    os.makedirs(cfg.out_dir, exist_ok=True)
    unita = os.path.splitdrive(os.path.abspath(cfg.out_dir))[0] or "/"
    libero = shutil.disk_usage(cfg.out_dir).free
    print(f"\n  [2] spazio: servono ~{_fmt(dim_scelta)}, liberi {_fmt(libero)} "
          f"su {unita}")
    if cfg.mode == "zip":
        print(f"      (in rete passano pero' {_fmt(dim_trasferita)}: il .zip va "
              "trasferito intero prima di poterlo aprire)")
    if dim_scelta > libero - cfg.reserve_gb * 1e9:
        print("      ATTENZIONE: non ci stanno tutti. Lo scarico procede in "
              "ordine di data e si ferma")
        print("      prima di esaurire il disco, lasciando liberi "
              f"{cfg.reserve_gb:.0f} GB.")
        if cfg.mode == "full":
            print("      Con --mode swath ne servirebbe circa un sesto.")

    if a.list_only:
        print("\n  --list-only: nessuno scarico eseguito.")
        return 0

    print("\n  [3] autenticazione")
    user, pwd = carica_credenziali(cfg.env_file)
    token = TokenCDSE(user, pwd, session)
    print(f"      autenticato come {user}")

    print(f"\n  [4] scarico di {len(prodotti)} prodotti in {cfg.out_dir}")
    esiti: List[Dict[str, Any]] = []
    scaricati = 0
    for i, prod in enumerate(prodotti, 1):
        nome = str(prod["Name"])
        atteso = int(prod.get("ContentLength", 0))
        data = str(prod.get("ContentDate", {}).get("Start", ""))[:10]
        print(f"\n    [{i}/{len(prodotti)}] {data}  {nome}")

        libero = shutil.disk_usage(cfg.out_dir).free
        if cfg.mode == "full":
            richiesto = atteso
        elif cfg.mode == "zip":
            richiesto = atteso + int(atteso / FRAZIONE_SWATH)
        else:
            richiesto = int(atteso / FRAZIONE_SWATH)
        if libero - richiesto < cfg.reserve_gb * 1e9:
            print(f"      spazio insufficiente ({_fmt(libero)} liberi, servono "
                  f"~{_fmt(richiesto)}")
            print(f"      lasciandone {cfg.reserve_gb:.0f} GB): mi fermo qui.")
            esiti.append({"nome": nome, "id": prod["Id"], "data": data,
                          "stato": "saltato_spazio"})
            break

        rec: Dict[str, Any]
        try:
            if cfg.mode == "full":
                rec = scarica_completo(prod, cfg, session, token)
            elif cfg.mode == "zip":
                rec = scarica_ed_estrai(prod, cfg, session, token)
            else:
                rec = scarica_swath(prod, cfg, session, token)
            if rec["stato"] == "scaricato":
                scaricati += 1
                print(f"      completato: {_fmt(rec['byte'])}")
        except Exception as ex:                             # pragma: no cover
            print(f"      ERRORE: {type(ex).__name__}: {ex}")
            rec = {"stato": "errore", "errore": f"{type(ex).__name__}: {ex}"}
        rec.update({"nome": nome, "id": prod["Id"], "data": data,
                    "orbita_relativa": attributo(prod, "relativeOrbitNumber"),
                    "direzione": attributo(prod, "orbitDirection")})
        esiti.append(rec)

    manifest = os.path.join(cfg.out_dir, "manifest_scarico.json")
    with open(manifest, "w", encoding="utf-8") as fh:
        json.dump({"generato": datetime.now(timezone.utc).isoformat(),
                   "filtro": costruisci_filtro(cfg),
                   "modalita": cfg.mode,
                   "swath": cfg.swath if cfg.mode == "swath" else None,
                   "polarizzazione": (cfg.polarisation if cfg.mode == "swath"
                                      else None),
                   "prodotti": esiti}, fh, indent=2, ensure_ascii=False)

    errori = [e for e in esiti if e["stato"] == "errore"]
    presenti = sum(1 for e in esiti if e["stato"] == "presente")
    print(f"\n  [5] fatto: {scaricati} scaricati, {presenti} gia' presenti, "
          f"{len(errori)} in errore")
    print(f"      riepilogo in {manifest}")
    if errori:
        print("      rilanciare il comando riprende dai file incompleti (.part).")
    return 1 if errori else 0


if __name__ == "__main__":
    sys.exit(main())

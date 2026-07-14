@echo off
rem ============================================================================
rem piramide_unificato.bat - lancia piramide_unificato.py (pipeline unica,
rem step 1-6: ricerca -> download -> box -> DInSAR -> forme d'onda 3D VV/VH ->
rem variazioni di densita')
rem
rem Uso:
rem   piramide_unificato.bat                        -> Kefren, step 4-6, dati locali (default)
rem   piramide_unificato.bat cheope                  -> Cheope, step 4-6, dati locali
rem   piramide_unificato.bat kefren --steps 1-6       -> pipeline completa (ricerca compresa)
rem   piramide_unificato.bat kefren --steps 2-6       -> ri-estrae/ricostruisce il box (usa anche
rem                                                      gli zip SAFE gia' in alta_risoluzione_out\)
rem   piramide_unificato.bat kefren --download        -> scarica anche i prodotti SLC dal CDSE
rem   piramide_unificato.bat kefren --vera-tomografia  -> + tomografia Doppler vera del paper
rem   piramide_unificato.bat cheope --pol vv --steps 4-6
rem   piramide_unificato.bat kefren --steps 4-6 --nz 512 --step-x 2
rem
rem Il primo argomento (opzionale) e' la piramide (cheope|kefren, default kefren);
rem TUTTI gli argomenti successivi sono passati cosi' come sono a
rem piramide_unificato.py (--steps, --pol, --layers, --download,
rem --alta-risoluzione-dir, --nz, --step-x, --step-y, --soglia,
rem --vera-tomografia, --debug, ...).
rem Se non viene passato esplicitamente --steps, il default e' --steps 4-6
rem (solo elaborazione: DInSAR + forme d'onda 3D VV/VH + variazioni di densita',
rem sul box gia' estratto -- niente ricerca/download dal CDSE).
rem
rem Lo step 2 (incluso in --steps 1-6 o 2-6) legge in automatico i .zip SAFE
rem gia' scaricati da scarica_alta_risoluzione.py in alta_risoluzione_out\
rem (accanto a questo file) e ne estrae i measurement VV/VH mancanti PRIMA di
rem provare il download CDSE: nessuna credenziale necessaria per questa parte.
rem I file .zip.part (download ancora in corso in un'altra finestra) vengono
rem ignorati automaticamente. Per usare un'altra cartella: --alta-risoluzione-dir
rem "C:\percorso\cartella" (oppure --alta-risoluzione-dir "" per disattivare).
rem
rem Vedi "python piramide_unificato.py --help" per l'elenco completo.
rem
rem Le credenziali CDSE (necessarie SOLO con --download) vengono chieste una
rem volta sola e salvate in .cdse.env accanto a questo file (gia' in
rem .gitignore, non finisce mai nel repository). Per cambiarle cancella
rem .cdse.env. NB: password che contengono "=" non sono supportate.
rem ============================================================================
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"

set "PIRAMIDE=kefren"
set "ARGS="
set "FIRST=1"
for %%A in (%*) do (
    if "!FIRST!"=="1" (
        set "FIRST=0"
        if /i "%%~A"=="cheope" (
            set "PIRAMIDE=cheope"
        ) else if /i "%%~A"=="kefren" (
            set "PIRAMIDE=kefren"
        ) else (
            set "ARGS=!ARGS! %%A"
        )
    ) else (
        set "ARGS=!ARGS! %%A"
    )
)

rem ---- default --steps 4-6 se l'utente non specifica esplicitamente gli step ----
echo !ARGS! | findstr /c:"--steps" >nul
if errorlevel 1 set "ARGS=!ARGS! --steps 4-6"

rem ---- credenziali CDSE: servono solo se e' stato chiesto --download ----
echo !ARGS! | findstr /c:"--download" >nul
if not errorlevel 1 (
    if exist ".cdse.env" (
        for /f "usebackq tokens=1,* delims==" %%X in (".cdse.env") do set "%%X=%%Y"
    )
    if not defined CDSE_USER (
        echo Credenziali CDSE non trovate: le chiedo ora ^(verranno salvate in .cdse.env^).
        set /p CDSE_USER=Email account CDSE:
    )
    if not defined CDSE_PASS (
        set /p CDSE_PASS=Password CDSE:
    )
    if not defined CDSE_USER goto :err_cred
    if not defined CDSE_PASS goto :err_cred
    (
        echo CDSE_USER=!CDSE_USER!
        echo CDSE_PASS=!CDSE_PASS!
    ) > ".cdse.env"
)

set "PYTHONIOENCODING=utf-8"
echo.
echo === piramide_unificato.py --pyramid %PIRAMIDE%!ARGS! ===
python piramide_unificato.py --pyramid %PIRAMIDE% !ARGS!
if errorlevel 1 goto :err_run

echo.
echo === Completato: output in goal_out_Cheope\unificato\ (Cheope) o goal_out_kefren\unificato\ (Kefren) ===
pause
exit /b 0

:err_cred
echo ERRORE: credenziali CDSE mancanti.
pause
exit /b 1

:err_run
echo.
echo ERRORE: la pipeline non e' terminata correttamente (vedi i messaggi sopra).
echo Se le credenziali CDSE sono sbagliate cancella .cdse.env e rilancia.
pause
exit /b 1

@echo off
rem ============================================================================
rem scarica_dati.bat - download automatico dei prodotti Sentinel-1 SLC dal CDSE
rem
rem Uso:
rem   scarica_dati.bat            -> scarica per Cheope (default)
rem   scarica_dati.bat kefren     -> scarica per Kefren
rem
rem Le credenziali CDSE (account gratuito su dataspace.copernicus.eu) vengono
rem chieste UNA volta e salvate in .cdse.env accanto a questo file; il file e'
rem gia' in .gitignore e non finisce mai nel repository. Per cambiarle basta
rem cancellare .cdse.env. NB: password che contengono "=" non sono supportate.
rem ============================================================================
setlocal EnableExtensions
cd /d "%~dp0"

set "PIRAMIDE=%~1"
if "%PIRAMIDE%"=="" set "PIRAMIDE=cheope"

rem ---- credenziali: da .cdse.env se esiste, altrimenti chiedile e salvale ----
if exist ".cdse.env" (
    for /f "usebackq tokens=1,* delims==" %%A in (".cdse.env") do set "%%A=%%B"
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
    echo CDSE_USER=%CDSE_USER%
    echo CDSE_PASS=%CDSE_PASS%
) > ".cdse.env"

rem ---- download (ripristinabile: rilanciare questo .bat riprende da dove era) ----
set "PYTHONIOENCODING=utf-8"
echo.
echo === Download prodotti per %PIRAMIDE% (interrompibile: riprende al rilancio) ===
python scarica_alta_risoluzione.py --pyramid %PIRAMIDE% --download
if errorlevel 1 goto :err_run

echo.
echo === Completato: prodotti in alta_risoluzione_out\ ===
pause
exit /b 0

:err_cred
echo ERRORE: credenziali mancanti.
pause
exit /b 1

:err_run
echo.
echo ERRORE: il download non e' andato a buon fine (vedi i messaggi sopra).
echo Se le credenziali sono sbagliate cancella .cdse.env e rilancia.
pause
exit /b 1

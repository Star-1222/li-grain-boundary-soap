@echo off
setlocal
set "PYTHON=C:\Users\1\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
if not exist "%PYTHON%" set "PYTHON=python"
"%PYTHON%" "%~dp0tools\artifact_ledger.py" %*
set EXITCODE=%ERRORLEVEL%
if not %EXITCODE%==0 (
    echo.
    echo Ledger generation failed. See messages above.
    pause
)
endlocal & exit /b %EXITCODE%

@echo off
setlocal
rem Locate a Python interpreter without hardcoding one machine's layout.
set "PYTHON="
where python >nul 2>nul && set "PYTHON=python"
if not defined PYTHON (
    set "BUNDLED=%USERPROFILE%\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
    if exist "%BUNDLED%" set "PYTHON=%BUNDLED%"
)
if not defined PYTHON (
    echo Python not found. Install Python or set RESEARCH_ROOT/PYTHON.
    pause
    exit /b 1
)
"%PYTHON%" "%~dp0tools\artifact_ledger.py" %*
set EXITCODE=%ERRORLEVEL%
if not %EXITCODE%==0 (
    echo.
    echo Ledger generation failed. See messages above.
    pause
)
endlocal & exit /b %EXITCODE%

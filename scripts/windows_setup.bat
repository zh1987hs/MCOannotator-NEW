@echo off
setlocal

where py >nul 2>nul
if %ERRORLEVEL%==0 (
  set PYCMD=py -3
) else (
  where python >nul 2>nul
  if %ERRORLEVEL%==0 (
    set PYCMD=python
  ) else (
    echo Python not found. Install Python 3.10+ and add to PATH.
    exit /b 1
  )
)

echo [1/3] Creating venv (.venv)...
%PYCMD% -m venv .venv || exit /b 1

echo [2/3] Upgrading pip...
.venv\Scripts\python.exe -m pip install --upgrade pip || exit /b 1

echo [3/3] Installing dependencies...
.venv\Scripts\python.exe -m pip install -r requirements.txt || exit /b 1

echo Done. No activation needed:
echo .venv\Scripts\python.exe -m mco_mnox.cli --help
echo .venv\Scripts\python.exe -m streamlit run mco_mnox\gui.py

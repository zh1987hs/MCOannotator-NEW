@echo off
setlocal

set REPO_ROOT=%~dp0..
pushd "%REPO_ROOT%"

where py >nul 2>nul
if %ERRORLEVEL%==0 (
  set PYCMD=py -3
) else (
  where python >nul 2>nul
  if %ERRORLEVEL%==0 (
    set PYCMD=python
  ) else (
    echo Python not found. Install Python 3.10+ and add to PATH.
    popd
    exit /b 1
  )
)

echo [1/3] Creating venv (.venv)...
%PYCMD% -m venv .venv || (popd & exit /b 1)

if not exist .venv\Scripts\python.exe (
  echo venv created but .venv\Scripts\python.exe not found.
  popd
  exit /b 1
)

echo [2/3] Installing dependencies...
.venv\Scripts\python.exe -m pip install --upgrade pip || (popd & exit /b 1)
.venv\Scripts\python.exe -m pip install -r requirements.txt || (popd & exit /b 1)

echo [3/3] Done. Test commands:
echo .venv\Scripts\python.exe -m mco_mnox.cli --help
echo .venv\Scripts\python.exe -m streamlit run mco_mnox\gui.py

popd

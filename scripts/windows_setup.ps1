# Windows one-shot setup without requiring Activate.ps1
$ErrorActionPreference = 'Stop'

Write-Host '[1/3] Detecting Python launcher...'
$pyCmd = $null
if (Get-Command py -ErrorAction SilentlyContinue) {
  $pyCmd = 'py -3'
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
  $pyCmd = 'python'
} else {
  throw 'Python not found. Please install Python 3.10+ and ensure py or python is in PATH.'
}

Write-Host "Using: $pyCmd"

Write-Host '[2/3] Creating venv (.venv)...'
Invoke-Expression "$pyCmd -m venv .venv"

Write-Host '[3/3] Installing dependencies...'
& .\.venv\Scripts\python.exe -m pip install --upgrade pip
& .\.venv\Scripts\python.exe -m pip install -r requirements.txt

Write-Host 'Done. Use these commands (no activation needed):'
Write-Host '.\.venv\Scripts\python.exe -m mco_mnox.cli --help'
Write-Host '.\.venv\Scripts\python.exe -m streamlit run mco_mnox/gui.py'

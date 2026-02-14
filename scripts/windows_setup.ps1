# Windows one-shot setup without requiring Activate.ps1
$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $RepoRoot
try {
  Write-Host "Repository root: $RepoRoot"
  Write-Host '[1/4] Detecting Python launcher...'
  $pyCmd = $null
  if (Get-Command py -ErrorAction SilentlyContinue) {
    $pyCmd = 'py -3'
  } elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $pyCmd = 'python'
  } else {
    throw 'Python not found. Please install Python 3.10+ and ensure py or python is in PATH.'
  }

  Write-Host "Using: $pyCmd"

  Write-Host '[2/4] Creating venv (.venv)...'
  Invoke-Expression "$pyCmd -m venv .venv"

  if (-not (Test-Path '.\.venv\Scripts\python.exe')) {
    throw 'venv created but .venv\Scripts\python.exe not found. Please check Python installation.'
  }

  Write-Host '[3/4] Installing dependencies...'
  & .\.venv\Scripts\python.exe -m pip install --upgrade pip
  & .\.venv\Scripts\python.exe -m pip install -r requirements.txt

  Write-Host '[4/4] Done. Test commands:'
  Write-Host '& .\.venv\Scripts\python.exe -m mco_mnox.cli --help'
  Write-Host '& .\.venv\Scripts\python.exe -m streamlit run mco_mnox/gui.py'
}
finally {
  Pop-Location
}

# End-to-end Windows PowerShell helper:
# 1) ensure venv exists, 2) install deps if needed, 3) launch GUI
$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot
Push-Location $RepoRoot
try {
  if (-not (Test-Path '.\.venv\Scripts\python.exe')) {
    Write-Host 'No .venv found. Running setup...'
    & .\windows_setup.ps1
  }

  Write-Host 'Python:'
  & .\.venv\Scripts\python.exe --version

  Write-Host 'Checking streamlit...'
  $ok = $true
  try {
    & .\.venv\Scripts\python.exe -c "import streamlit" | Out-Null
  } catch {
    $ok = $false
  }
  if (-not $ok) {
    Write-Host 'Installing requirements...'
    & .\.venv\Scripts\python.exe -m pip install -r requirements.txt
  }

  Write-Host 'Launching GUI at http://localhost:8501 ...'
  & .\.venv\Scripts\python.exe -m streamlit run mco_mnox/gui.py
}
finally {
  Pop-Location
}

# Wrapper: run Windows setup from repository root
$ErrorActionPreference = 'Stop'
$script = Join-Path $PSScriptRoot 'scripts\windows_setup.ps1'
if (-not (Test-Path $script)) { throw "Cannot find $script" }
& $script

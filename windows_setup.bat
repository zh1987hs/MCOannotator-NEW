@echo off
REM Wrapper: run Windows setup from repository root
if exist scripts\windows_setup.bat (
  call scripts\windows_setup.bat
) else (
  echo Cannot find scripts\windows_setup.bat
  exit /b 1
)

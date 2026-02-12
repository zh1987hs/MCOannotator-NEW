@echo off
REM Launch Streamlit GUI on Windows (prefer venv python, no activation required)
if exist .venv\Scripts\python.exe (
  .venv\Scripts\python.exe -m streamlit run mco_mnox\gui.py
) else (
  python -m streamlit run mco_mnox\gui.py
)

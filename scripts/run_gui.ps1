# Launch Streamlit GUI on Windows PowerShell (prefer venv python, no activation required)
if (Test-Path '.\.venv\Scripts\python.exe') {
  & .\.venv\Scripts\python.exe -m streamlit run mco_mnox/gui.py
} else {
  python -m streamlit run mco_mnox/gui.py
}

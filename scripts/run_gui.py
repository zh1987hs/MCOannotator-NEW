"""Launch Streamlit GUI for MCO MnOx PU workflow."""

if __name__ == "__main__":
    import subprocess
    import sys

    cmd = [sys.executable, "-m", "streamlit", "run", "mco_mnox/gui.py"]
    raise SystemExit(subprocess.call(cmd))

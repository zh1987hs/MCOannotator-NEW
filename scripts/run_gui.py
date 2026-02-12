"""Launch Streamlit GUI for MCO MnOx PU workflow."""

from pathlib import Path
import subprocess
import sys


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[1]
    app_path = repo_root / "mco_mnox" / "gui.py"
    cmd = [sys.executable, "-m", "streamlit", "run", str(app_path)]
    raise SystemExit(subprocess.call(cmd, cwd=str(repo_root)))

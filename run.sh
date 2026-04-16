#!/bin/bash
set -e
APP_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$APP_DIR"
if [ ! -d ".venv" ]; then python3 -m venv .venv; fi
source .venv/bin/activate
python -m pip install --upgrade pip >/dev/null
python -m pip install -r requirements.txt
python startup_check.py
PORT="${PORT:-8501}"
exec python -m streamlit run app.py --server.port "$PORT"

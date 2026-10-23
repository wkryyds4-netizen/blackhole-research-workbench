@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Please create the Python environment described in README.md first.
  pause
  exit /b 1
)
echo Open http://127.0.0.1:8765 after the Local message appears.
".venv\Scripts\python.exe" server.py
pause

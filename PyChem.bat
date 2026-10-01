@echo off
setlocal enabledelayedexpansion

cd /d "%~dp0"

echo.
echo =============================
echo   PyChem Launcher
echo =============================

:: Step 1: Check Python
echo.
echo Checking for Python...

set PYTHON=

for %%P in (python python3) do (
    where %%P >nul 2>&1
    if !errorlevel! == 0 (
        %%P -c "import sys; exit(0 if sys.version_info >= (3,10) else 1)" >nul 2>&1
        if !errorlevel! == 0 (
            set PYTHON=%%P
            goto :found
        )
    )
)

:: Step 2: Python not found
echo.
echo Python 3.10+ not found.
echo Please install Python from https://www.python.org/downloads/
echo Make sure to check "Add Python to PATH"
pause
exit /b 1

:found
echo Found Python: %PYTHON%

:: Step 3: Create virtual environment
if not exist venv (
    echo.
    echo Creating virtual environment...
    %PYTHON% -m venv venv
)

:: Step 4: Install dependencies
echo.
echo Installing dependencies...
venv\Scripts\python -m pip install --upgrade pip
venv\Scripts\python -m pip install -r requirements.txt

:: Step 5: Run project
echo.
echo Launching PyChem...
venv\Scripts\python main.py

pause

@echo off
setlocal enabledelayedexpansion
chcp 65001 >nul 2>&1

:: PyChem launcher for Windows — double-click to run PyChem.
:: Automatically installs Python and dependencies on first run.

cd /d "%~dp0"

echo.
echo   =============================================
echo    PyChem Launcher
echo   =============================================

:: ── Step 1: Find Python 3.10+ ─────────────────────────────────────────────
echo.
echo   ^>^>  Checking for Python 3.10+...

set PYTHON=

:: Search PATH first
for %%P in (python3.13 python3.12 python3.11 python3.10 python3 python) do (
    where %%P >nul 2>&1
    if !errorlevel! == 0 (
        %%P -c "import sys; exit(0 if sys.version_info >= (3,10) else 1)" >nul 2>&1
        if !errorlevel! == 0 (
            set PYTHON=%%P
            goto :python_found
        )
    )
)

:: Check common install locations not yet on PATH (e.g. fresh winget install)
for %%V in (313 312 311 310) do (
    for %%D in (
        "%LOCALAPPDATA%\Programs\Python\Python%%V\python.exe"
        "%ProgramFiles%\Python%%V\python.exe"
        "%ProgramFiles(x86)%\Python%%V\python.exe"
    ) do (
        if exist %%D (
            %%D -c "import sys; exit(0 if sys.version_info >= (3,10) else 1)" >nul 2>&1
            if !errorlevel! == 0 (
                set PYTHON=%%D
                goto :python_found
            )
        )
    )
)

:: ── Step 2: Auto-install Python ───────────────────────────────────────────
echo   X   Python 3.10+ not found.
echo   ^>^>  Attempting to install Python automatically...
echo.

:: Try winget (built into Windows 10/11)
where winget >nul 2>&1
if !errorlevel! == 0 (
    echo   ^>^>  Installing Python 3.12 via winget...
    winget install --id Python.Python.3.12 -e --silent --accept-package-agreements --accept-source-agreements
    :: Refresh PATH from registry so the new Python is visible
    for /f "tokens=*" %%i in ('powershell -NoProfile -Command "[System.Environment]::GetEnvironmentVariable(\"PATH\",\"Machine\") + \";\" + [System.Environment]::GetEnvironmentVariable(\"PATH\",\"User\")"') do set PATH=%%i
    :: Re-check known locations
    for %%V in (312 313 311 310) do (
        for %%D in (
            "%LOCALAPPDATA%\Programs\Python\Python%%V\python.exe"
            "%ProgramFiles%\Python%%V\python.exe"
        ) do (
            if exist %%D (
                set PYTHON=%%D
                goto :python_found
            )
        )
    )
    where python >nul 2>&1
    if !errorlevel! == 0 (
        set PYTHON=python
        goto :python_found
    )
)

:: Fallback: download Python installer silently
echo   ^>^>  winget not available. Downloading Python 3.12 installer...
set INSTALLER=%TEMP%\python_installer.exe
powershell -NoProfile -Command "Invoke-WebRequest -Uri 'https://www.python.org/ftp/python/3.12.9/python-3.12.9-amd64.exe' -OutFile '%INSTALLER%'" >nul 2>&1
if exist "%INSTALLER%" (
    echo   ^>^>  Running Python installer (this may take a moment)...
    "%INSTALLER%" /quiet InstallAllUsers=0 PrependPath=1 Include_test=0
    del "%INSTALLER%" >nul 2>&1
    set PYTHON=%LOCALAPPDATA%\Programs\Python\Python312\python.exe
    if exist "!PYTHON!" goto :python_found
)

:: ── Step 3: Final check — all methods failed ──────────────────────────────
echo.
echo   X   Could not install Python automatically.
echo       Please install Python 3.10+ from https://www.python.org/downloads/
echo       Make sure to check "Add Python to PATH" during installation.
echo       Then double-click PyChem.bat again.
echo.
powershell -NoProfile -Command "Add-Type -AssemblyName PresentationFramework; $r = [System.Windows.MessageBox]::Show('PyChem could not install Python automatically.`n`nPlease install Python 3.10+ from python.org`nMake sure to check \"Add Python to PATH\" during installation,`nthen open PyChem again.', 'PyChem — Python Not Found', 'OK', 'Error'); if ($r -eq 'OK') { Start-Process 'https://www.python.org/downloads/' }" 2>nul
pause
exit /b 1

:python_found
for /f "tokens=*" %%V in ('"%PYTHON%" --version 2^>^&1') do echo   OK  Using %%V

:: ── Step 4: Create shortcut with icon (once) ──────────────────────────────
if not exist ".launcher_icon_set" (
    set "PS1=%TEMP%\pychem_lnk.ps1"
    echo $ws = New-Object -ComObject WScript.Shell                    > "%TEMP%\pychem_lnk.ps1"
    echo $s  = $ws.CreateShortcut("%~dp0PyChem - Launch.lnk")       >> "%TEMP%\pychem_lnk.ps1"
    echo $s.TargetPath      = "%~f0"                                  >> "%TEMP%\pychem_lnk.ps1"
    echo $s.IconLocation    = "%~dp0assets\icon.ico"                  >> "%TEMP%\pychem_lnk.ps1"
    echo $s.WorkingDirectory = "%~dp0"                                >> "%TEMP%\pychem_lnk.ps1"
    echo $s.Description     = "Launch PyChem"                         >> "%TEMP%\pychem_lnk.ps1"
    echo $s.Save()                                                     >> "%TEMP%\pychem_lnk.ps1"
    powershell -NoProfile -ExecutionPolicy Bypass -File "%TEMP%\pychem_lnk.ps1" >nul 2>&1
    del "%TEMP%\pychem_lnk.ps1" >nul 2>&1
    echo. > ".launcher_icon_set"
    echo   OK  Shortcut created: PyChem - Launch.lnk  (use this for the icon^)
)

:: ── Step 5: Create venv and install dependencies (first run only) ──────────
if not exist venv (
    echo.
    echo   ^>^>  First run -- installing dependencies, please wait ~2 minutes...
    "%PYTHON%" -m venv venv
    venv\Scripts\python -m pip install --upgrade pip -q
    venv\Scripts\python -m pip install -r requirements.txt
    echo   OK  Dependencies installed.
)

:: ── Step 6: Launch PyChem ──────────────────────────────────────────────────
echo.
echo   ^>^>  Launching PyChem...
venv\Scripts\python main.py

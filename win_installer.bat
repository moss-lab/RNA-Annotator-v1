@echo off
REM --- RNA-Annotator: Minimal Installer for Windows ---

echo --- Welcome to the RNA-Annotator Setup for Windows ---

:: --- Configuration ---
:: !!! IMPORTANT: Replace this URL with your actual Zenodo direct download link !!!
set "DATA_URL=https://zenodo.org/records/16953760/files/resources_data_sets.zip?download=1"
set "DATA_ARCHIVE=resources_data_sets.zip"

:: --- Pre-flight Checks ---
where conda >nul 2>nul
if %errorlevel% neq 0 (
    echo ERROR: Conda is not found. Please run this script from an Anaconda Prompt.
    echo See: https://docs.conda.io/projects/miniconda/en/latest/
    pause
    exit /b 1
)

where curl >nul 2>nul
if %errorlevel% neq 0 (
    echo ERROR: 'curl' was not found. Please ensure you are on a modern version of Windows or install curl.
    pause
    exit /b 1
)

echo.
echo Pre-flight checks passed.
echo This script will now download the required data package (~16 GB).
echo This may take a long time.
pause

:: --- Step 1: Download and Extract the Data Package ---
echo.
echo --^> Step 1 of 2: Downloading data from Zenodo...
curl -L -o %DATA_ARCHIVE% %DATA_URL%
if %errorlevel% neq 0 (
    echo ERROR: Data download failed. Please check your internet connection.
    pause
    exit /b 1
)

echo.
echo Download complete. Extracting files...
tar -xf %DATA_ARCHIVE%
if %errorlevel% neq 0 (
    echo ERROR: Failed to unzip the data archive.
    pause
    exit /b 1
)
del %DATA_ARCHIVE%
echo Data extracted successfully.

:: --- Step 2: Create the Conda Environment ---
echo.
echo --^> Step 2 of 2: Creating the 'rna-annotator' conda environment...
call conda env create -f environment.yml
if %errorlevel% neq 0 (
    echo ERROR: Failed to create the conda environment.
    pause
    exit /b 1
)
echo Conda environment created successfully.

:: --- Final Instructions ---
echo.
echo --- SETUP COMPLETE! ---
echo.
echo To use the tool, open a new Anaconda Prompt and run:
echo     conda activate rna-annotator
pause

#!/bin/bash

# --- RNA-Annotator: A One-Command Installer for macOS & Linux ---

echo "--- Welcome to the RNA-Annotator Setup ---"

# --- Configuration ---
# !!! IMPORTANT: Replace this URL with your actual Zenodo direct download link !!!
DATA_URL="https://zenodo.org/records/16953760/files/resources_data_sets.zip?download=1"
DATA_ARCHIVE="resources_data_sets.zip"

# --- Pre-installation Checks ---
# 1. Check if Conda is installed and available in the PATH.
if ! command -v conda >/dev/null 2>&1; then
    echo "❌ Error: Conda is not found. Please install Miniconda or Anaconda first."
    echo "   See: https://docs.conda.io/projects/miniconda/en/latest/"
    exit 1
fi

# 2. Check for a download tool (prefer curl, fallback to wget).
if command -v curl >/dev/null 2>&1; then
    DOWNLOAD_CMD="curl -L -o"
elif command -v wget >/dev/null 2>&1; then
    DOWNLOAD_CMD="wget -O"
else
    echo "❌ Error: Neither 'curl' nor 'wget' was found. Please install one to proceed."
    exit 1
fi

echo "✅ All requirements are met."
echo "The script will now download the required data package (~16 GB) and create the conda environment."
read -p "Press [Enter] to continue or Ctrl+C to cancel."
echo

# --- Step 1: Download and Extract the Data Package ---
echo "--> Step 1 of 3: Downloading data from Zenodo... (This will take a long time)"
$DOWNLOAD_CMD "$DATA_ARCHIVE" "$DATA_URL"
if [ $? -ne 0 ]; then
    echo "❌ Error: Data download failed. Please check your internet connection."
    exit 1
fi

echo "✅ Download complete. Extracting files..."
unzip -q "$DATA_ARCHIVE"  # -q for "quiet" mode
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to unzip the data archive. Please make sure 'unzip' is installed."
    exit 1
fi
rm "$DATA_ARCHIVE"
echo "✅ Data extracted successfully. The 'resources_data_sets' folder is now in place."

# --- Step 2: Create the Conda Environment ---
echo
echo "--> Step 2 of 3: Creating the 'rna-annotator' conda environment..."
conda env create -f environment.yml
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to create the conda environment."
    exit 1
fi
echo "✅ Conda environment created successfully."

# --- Step 3: Index the Data Files (if needed) ---
# This step is commented out because best practice is to include the .tbi files in your zip.
# If you did NOT include them, you would uncomment the following lines.
# echo
# echo "--> Step 3 of 3: Indexing data files..."
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate rna-annotator
# tabix -p vcf resources_data_sets/clinvar.vcf.gz
# tabix -p bed resources_data_sets/hg38.phastCons100way.bedGraph.gz
# conda deactivate
# echo "✅ Indexing complete."


# --- Final Instructions ---
echo
echo "--- ✅ SETUP COMPLETE! ---"
echo
echo "To use the tool, you must first activate the conda environment in any new terminal:"
echo "    conda activate rna-annotator"
echo 
echo "You can now run RNA-Annotator commands. For help, type:"
echo "    rna-annotator --help"

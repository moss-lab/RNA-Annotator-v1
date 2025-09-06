    
#!/bin/bash

# --- RNA-Annotator: A One-Command Installer for macOS & Linux ---

echo "--- Welcome to the RNA-Annotator Setup ---"

# --- Configuration ---
DATA_URL="https://zenodo.org/records/16953760/files/resources_data_sets.zip?download=1"
DATA_ARCHIVE="resources_data_sets.zip"
RETRY_DELAY=10 # Seconds to wait before retrying a failed download

# --- Pre-installation Checks ---
# 1. Check if Conda is installed and available in the PATH.
if ! command -v conda >/dev/null 2>&1; then
    echo "❌ Error: Conda is not found. Please install Miniconda or Anaconda first."
    echo "   See: https://docs.conda.io/projects/miniconda/en/latest/"
    exit 1
fi

# 2. Check for a download tool and set it up for resuming.
if command -v curl >/dev/null 2>&1; then
    # -L: follow redirects, -C -: continue/resume, -o: output to file
    DOWNLOAD_CMD="curl -L -C - -o"
elif command -v wget >/dev/null 2>&1; then
    # -c: continue/resume, -O: output to file
    DOWNLOAD_CMD="wget -c -O"
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

# This 'until' loop will automatically retry the download until it succeeds.
until $DOWNLOAD_CMD "$DATA_ARCHIVE" "$DATA_URL"
do
    echo "❌ Download interrupted. Retrying in $RETRY_DELAY seconds..."
    sleep $RETRY_DELAY
done

echo "✅ Download complete. Extracting files..."
unzip -q "$DATA_ARCHIVE"
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to unzip the data archive. Please make sure 'unzip' is installed."
    exit 1
fi
rm "$DATA_ARCHIVE"
echo "✅ Data extracted successfully. The 'resources_data_sets' folder is now in place."

# --- Step 2: Create the Conda Environment ---
echo
echo "--> Step 2 of 3: Creating the 'rna-annotator' conda environment..."
# ... (rest of the script is unchanged) ...
conda env create -f environment.yml
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to create the conda environment."
    exit 1
fi
echo "✅ Conda environment created successfully."

# --- Step 3: Index the Data Files (if needed) ---
# ... (rest of the script is unchanged) ...

# --- Final Instructions ---
echo
echo "--- ✅ SETUP COMPLETE! ---"
echo
echo "To use the tool, you must first activate the conda environment in any new terminal:"
echo "    conda activate rna-annotator"
echo
echo "You can now run RNA-Annotator commands. For help, type:"
echo "    rna-annotator --help"

  

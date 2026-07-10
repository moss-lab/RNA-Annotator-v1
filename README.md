
<img width="728" height="355" alt="Screenshot 2025-08-19 at 8 10 12 PM" src="https://github.com/user-attachments/assets/1c5306e9-81b5-4a41-9901-5090da04d00b" />


# RNA-Annotator: A Comprehensive RNA Annotation Pipeline

RNA-Annotator is a powerful, command-line tool designed to rapidly fetch, process, and visualize a rich set of genomic and functional annotations for any given human RNA transcript.

By providing either a genomic coordinate or an Ensembl Transcript ID, you can generate a full suite of publication-ready data tracks for analysis in the Integrative Genomics Viewer (IGV). This pipeline automates the complex task of aggregating data from multiple major bioinformatics databases and APIs—including specialized RNA structure analysis from the **ScanFold database**—into a single, cohesive view.

## Output Files

The pipeline generates a main `[query]_results/` directory containing up to three subfolders:

1.  **`bed_tracks/`**: Contains `.bed`, `.bedGraph`, `.wig`, and `.bp` files formatted for direct visualization as tracks in the IGV genome browser.
2.  **`detailed_results/`**: Contains `.tsv` and `.csv` files with the raw, detailed information for every feature found in the analyzed region.
3.  **`[transcript_id]_ScanFold/`**: This directory is created when using the `-download_scanfold` option. It contains all the raw output files unzipped from the [StructuromeDB](https://www.structurome.bb.iastate.edu/azt/) database for the specified transcript.

## Features

*   Accepts genomic coordinates (`chr:start-end`) or Ensembl Transcript IDs (`ENST...`).
*   Fetches data from major resources like Ensembl, ClinVar, Dfam, and RMBase.
*   **Integrates RNA structure analysis from ScanFold**, generating tracks for base-pairing (arcs), MFE, z-score, and ensemble diversity.
*   Integrates local and remote datasets for conservation, variants, modifications, and protein binding sites.
*   Generates standard `.bed`, `.bedGraph`, `.wig`, and interaction (`.bp`) track files.
*   Automates launching and loading all generated tracks into IGV (macOS, Linux & Windows).

## Installation

This pipeline is distributed as a Conda environment and includes installer scripts to automate the setup process.

> **Note:** The setup will download a large data package (~29 GB). Please ensure you have a stable internet connection and sufficient disk space before you begin.

### Step 1: Clone the Repository

First, clone this repository to your local machine using Git:
```bash
git clone https://github.com/moss-lab/RNA-Annotator-v1.git
cd RNA-Annotator-v1
```

### Step 2: Run the Installer

The installer will automatically download the required data, create the `rna-annotator` conda environment, and set up all necessary software. Please choose the command for your operating system.

#### For macOS and Linux Users:

In your terminal, run the following command. This will take a significant amount of time due to the large data download.
```bash
chmod +x mac_linux_installer.sh
./mac_linux_installer.sh
```

#### For Windows Users:

Open **Anaconda Prompt** from your Start Menu and run the following command:
```cmd
win_installer.bat
```

## Usage

### 1. Activate the Conda Environment

Before running the tool, you **must** activate the conda environment in every new terminal session.
```bash
conda activate rna-annotator
```

### 2. Run the Tool

The basic command structure is: `python rna-annotator.py [QUERY] [OPTIONS]`

#### Arguments Explained

**`[QUERY]` (Required)**
The genomic region or transcript you want to analyze. This is the only required argument. It can be in one of three formats:
1.  **Genomic Coordinates:** `chr:start-end` (e.g., `chr19:58345183-58353492`)
2.  **Ensembl Transcript ID (unversioned):** `ENST...` (e.g., `ENST00000263100`)
3.  **Ensembl Transcript ID (versioned):** `ENST...` (e.g., `ENST00000263100.8`)

**`[OPTIONS]` (Optional Flags)**
These flags tell the tool which analyses to perform. You can combine as many as you like.

##### ScanFold RNA Structure Analysis
> **Note:** All ScanFold options (`-scanfold_bp`, `-scanfold_mfe`, etc.) **require** the `-download_scanfold` flag to be used in the same command. The input query for these options must be an **Ensembl Transcript ID**.

*   `-download_scanfold`: Downloads precomputed ScanFold data for a full transcript.
*   `-scanfold_bp`: Generates a track to visualize ScanFold base-pairing as arcs in IGV.
*   `-scanfold_mfe`: Generates a wig track for ScanFold Minimum Free Energy (MFE).
*   `-scanfold_zscore`: Generates a wig track for ScanFold z-score.
*   `-scanfold_ed`: Generates a wig track for ScanFold Ensemble Diversity (ED).
##### Standard Genomic Annotations

* `-refseq_functional`: Extracts RefSeq functional element annotations.
* `-eclips`: Extracts eCLIP-seq peaks to identify RNA-binding protein (RBP) sites.
* `-SNP`: Fetches known SNPs for the region from the Ensembl REST API.
* `-miRNA`: Extracts known miRNA annotations from a local GFF3 file.
* `-chem_mod`: Fetches known RNA chemical modifications from the RMBase database.
* `-polyA`: Extracts polyadenylation sites from the PolyASite 2.0 database.
* `-repeated_element`: Fetches repetitive element annotations from the Dfam API.
* `-chemical_prop`: Extracts local chemical probing/RNA structural reactivity tracks from local WIG files.
* `-chemical_prop_list`: Lists all available chemical probing/RNA structural reactivity dataset names.
* `-chemical_prop_group_list`: Lists available chemical probing/RNA structural reactivity method groups.
* `-clinvar`: Extracts clinical variants from a local ClinVar VCF file.
* `-target_scan`: Extracts predicted miRNA binding sites from local TargetScan BED files.
* `-phastCons`: Fetches evolutionary conservation scores from the phastCons 100-way vertebrate track.
* `-CpG_islands`: Extracts CpG islands from the UCSC CpG Islands Track.
* `-SpliceVar`: Extracts splice variants from the UCSC SpliceVarDB Track.
* `-Alt_Events`: Extracts alternative splicing events from the UCSC Alt Events Track.
* `-TFs`: Extracts transcription factor binding sites from the UCSC TFBS Track.
* `-GTEX`: Extracts RNA expression coverage per tissue from UCSC GTEx tracks.
* `-gtex_list`: Lists all available GTEx tissue names.

##### Selecting GTEx and Chemical Probing Tracks

By default, `-GTEX` processes all available GTEx tissue tracks. To avoid generating too many tracks, users can first list available tissues:

```bash
python rna-annotator.py -gtex_list
```

Then, selected tissues can be provided after `-GTEX`:

```bash
python rna-annotator.py chr19:58345183-58353492 -GTEX Liver
```

Multiple tissues can be provided as a comma-separated list:

```bash
python rna-annotator.py chr19:58345183-58353492 -GTEX Liver,Ovary
```

Similarly, `-chemical_prop` alone processes all available chemical probing/RNA structural reactivity tracks. Because this can generate many tracks, users can first list all available full dataset names:

```bash
python rna-annotator.py -chemical_prop_list
```

Users can also list only the available chemical probing method groups:

```bash
python rna-annotator.py -chemical_prop_group_list
```

Then, selected groups or datasets can be provided after `-chemical_prop`. For example, to select all icSHAPE-related tracks:

```bash
python rna-annotator.py chr19:58345183-58353492 -chemical_prop icSHAPE
```

To select multiple method groups:

```bash
python rna-annotator.py chr19:58345183-58353492 -chemical_prop DMS-seq,PARS
```

To select one specific full dataset name:

```bash
python rna-annotator.py chr19:58345183-58353492 -chemical_prop icSHAPE_hg38_Cell_Research_2021_NAI-N3_in_vivo_plus_transcriptome-wide_HeLa-pos
```

Multiple full dataset names can also be provided as a comma-separated list:

```bash
python rna-annotator.py chr19:58345183-58353492 -chemical_prop dataset_name_1,dataset_name_2
```

##### IGV Integration
*   `-igv`: After analysis, automatically launch and load results into IGV.
*   `--igv-path PATH`: **Required if using `-igv`**. Provide the full path to your IGV application.

### Examples

**Example 1: Analyze a region by coordinates with selected analyses.**
```bash
python rna-annotator.py chr19:58345183-58353492 -clinvar -phastCons -SNP
```

**Example 2: Analyze a transcript by ID and run only the ScanFold analyses.**
```bash
python rna-annotator.py ENST00000263100.8 -download_scanfold -scanfold_bp -scanfold_mfe -scanfold_zscore -scanfold_ed
```

**Example 3: Run *every possible analysis* on a transcript and visualize in IGV.**
This command demonstrates using all available flags for a comprehensive annotation.
```bash
# Note: The path to IGV will be different on your system.
python rna-annotator.py ENST00000263100.8 \
-download_scanfold \
-scanfold_bp \
-scanfold_mfe \
-scanfold_zscore \
-scanfold_ed \
-refseq_functional \
-eclips \
-SNP \
-miRNA \
-chem_mod \
-polyA \
-repeated_element \
-chemical_prop \
-GTEX \
-clinvar \
-target_scan \
-phastCons \
-Alt_Events \
-CpG_islands \
-SpliceVar \
-TFs \
-igv \
--igv-path "/Applications/IGV_2.16.2.app"
```

---
## Visualizing Results with IGV

To view the graphical representation of the results, you will need the **Integrative Genomics Viewer (IGV)**.

If you do not have IGV installed, you can download it for free from the official Broad Institute website. The tool will guide you through the simple installation process.
*   **Download IGV here:** [https://software.broadinstitute.org/software/igv/download](https://software.broadinstitute.org/software/igv/download)

You can then use the `-igv` and `--igv-path` flags as shown in the examples above to automatically load your results.

---

## Citation

If you use RNA-Annotator in your research, please cite:
> Dapour, A. O., Manigrasso, J., Rouse, W. B., La Sala, G., Aronnax, M. V., De Maria, L., Hamilton, R. S., Martinez Cuesta, S., Hogner, A., & Moss, W. N. (2026). RNAStructuromeDB: a transcriptome-wide database of predicted RNA secondary structures with integrated APIs for functional annotation and RNA-targeted drug discovery. *NAR Genomics and Bioinformatics, 8*(2), lqag044. https://doi.org/10.1093/nargab/lqag044

## Contact
For questions, bug reports, or suggestions, please contact Abdelraouf at: **raouf@iastate.edu**.

## License

This project is released under the **MIT License**. See the [LICENSE](LICENSE) file for details.

```

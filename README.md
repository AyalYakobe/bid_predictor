# Project Structure Overview

This project is organized into three main folders, each serving a distinct purpose:

## 1. Jupyter Notebooks and Bash Scripts

This folder contains Jupyter notebooks and bash scripts for downloading accession data. The notebooks provide a step-by-step guide and interactive environment for data retrieval and preliminary analysis. Bash scripts automate the downloading of data from various sources.

**accession_scraping.ipynb**

This Jupyter notebook is designed to scrape accession numbers from the NCBI database. It specifically extracts SRR (Sequence Read Archive) numbers associated with each accession. Users only need to update the `start_url` variable to point to the desired study on the NCBI website from which they wish to scrape data. This allows for efficient and targeted data gathering tailored to specific research needs.

**sra.sh**

This script automates the download of FASTQ files from the NCBI database using accession numbers listed in an input file. It ensures streamlined retrieval of sequence data for subsequent bioinformatics analysis.

- **Automated Directory Creation**: Creates a directory `fastq_files` for storing downloaded FASTQ files.
- **Detailed Logging**: Maintains a log of the download process in `fastq_download.log`, capturing both successful operations and errors.

To efficiently use this script for downloading FASTQ files, follow these steps:

1. **Configure `fastq-dump`:** Ensure that `fastq-dump` is properly installed and configured on your system to handle downloads from NCBI by following the instructions here:
https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump 


2. **Prepare the Accession List:**
   - Navigate to the `bin` directory within the provided ZIP folder.
   - Create and open the `accession_list.txt` file in this filepath.
   - Copy and paste your list of SRA accession numbers (`SRXXXX`) into this file.

3. **Execute the Script:**
   - Run the script by entering the following command in your terminal:
     ```bash
     ./sra
     ```
     
## 2. Data
Located within this folder are the datasets necessary for running the classification models. This includes training data, test data, and any other datasets required by the scripts. Included in this folder is also the dataset of `bacteria.txt` which houses the bacteria which has been noticed to have the most abundance in our model from `classify_ibd.py`.

Furthermore, the following .csv files for the three primary datasets (Liu/Lloyd-Price, Liu/Lloyd-Price/Gevers1, and Gevers1) following DADA2 abundance calculation:

`abundance_data_with_taxonomy_wspecies_LiuLloyd.csv`- Output from AbundanceCalculation.R following DADA2 for data from the Liu and Lloyd-Price datasets.
`abundance_data_with_taxonomy_wspecies_LiuLloydGevers1.csv`- Output from AbundanceCalculation.R following DADA2 for data from the Liu, Lloyd-Price, and Gevers datasets.
`abundance_data_with_taxonomy_wspecies_Gevers1.csv`- Output from AbundanceCalculation.R following DADA2 for data from the Gevers dataset.
`LiuLloyd_AbundanceData_WSpecies_WMetadata_IBD.csv`- Output following combining output from Metadata.py when adapted to the Liu and Lloyd-Price datasets with the corresponding abundance data. Input for PreProcessing.py.
`LiuLloydGevers1_AbundanceData_WSpecies_WMetadata_IBD.csv`- Output following combining output from Metadata.py when adapted to the Liu, Lloyd-Price, and Gevers datasets with the corresponding abundance data. Input for PreProcessing.py.
`Gevers1_AbundanceData_WSpecies_WMetadata_IBD.csv`- Output following combining output from Metadata.py when adapted to the Gevers datasets with the corresponding abundance data. Input for PreProcessing.py.
`LiuLloyd_AbundanceData_WSpecies_WMetadata_IBD_Normalized_WStudy_WNA.csv`- Output following PreProcessing.py when adapted to respective input file for Liu and Lloyd-Price datasets. Rows are normalized to get relative abundance, taxa not in a study are set to NaN, and study metadata is added.
`LiuLloydGevers1_AbundanceData_WSpecies_WMetadata_IBD_Normalized_WStudy_WNA.csv`- Output following PreProcessing.py when adapted to respective input file for Liu, Lloyd-Price, and Gevers1 datasets. Rows are normalized to get relative abundance, taxa not in a study are set to NaN, and study metadata is added.
`LiuLloyd_AbundanceData_WSpecies_WMetadata_IBD_Normalized.csv`- Output following PreProcessing.py when adapted to respective input file for Gevers dataset. Rows are normalized to get relative abundance. Since only one study, taxa not in a study are not set to NaN and study metadata is not added.

Additionally, subfolder `accession_files` contains the text files set to accession_file in order for .sra bash script to run. Lastly, `metadata_files` contains metadata files used by Metadata.py and are labeled by their respective study. 
## 3. Classifier_Scripts

This folder contains all the scripts used for running classification models. These scripts are crucial for processing the data and applying the machine learning algorithms.

**classify_ibd.py**
The function `current_ibd` is designed to classify individuals as having Inflammatory Bowel Disease (IBD) or not and is housed in this file.

- **Running the Function:**
  The `current_ibd` function is executed through the `main` function in the same file.
- **Output:**
  When run, `main` not only classifies the cases but also provides relevant statistics about the model's performance, aiding in the evaluation and understanding of the results.

**future_microbiome_state.py**

This script processes time-series microbiome data to understand subject repetition and perform statistical modeling. It primarily:
1. Preprocesses the data by removing NaN values and calculating the frequency of subject entries.
2. Utilizes Principal Component Analysis (PCA) to reduce dimensionality before applying a Vector Autoregression (VAR) model to forecast future microbiome states based on historical data.
3. Everything is executed through the `var_attempt` function, indicating streamlined testing and deployment.

Note: This script is under development and does not work at the moment due challenges with insufficient time-series data per subject and issues with PCA due to limited observations.

## 3. Data_Processing

This folder contains all scripts used for DADA2 calculations, metadata extraction, and pre-processing abundance data.

**AbundanceCalculation.R**
This script is an R script that takes in .fastq.gz files downloaded via .sra script using the accession files in `accession_files` under `data`. Additionally, this script requires taxonomic train datasets downloaded from the SILVA website (`silva_nr99_v138.1_train_set.fa.gz` and `silva_species_assignment_v138.1.fa.gz`). For each study, this script undergoes the full DADA2 protocol. For the Gevers study, due to RAM limitations this study was separated into 6 batches. This can be edited if RAM limitation are no longer a concern. Following DADA2 abundacne calculation and taxonomic assignment, outputs were imported into phyloseq objects, merged, and written to .csv files as an abundance table with taxonomic names. All packages that require installation are installed within the script.

**Metadata.py**
This script is a Python script that extracts metadata from appropriate files collected from the NCBI database, from individual studies, and the Qita database. These files containing metadata are within the `metadata_files` folder within `data`. Following the proper extraction of metadata, the metadata is arranged by study and sorted chronologically. A secondary purpose of this script was to determine the accession list for the Gevers study as only a subset of the files within the BioProject fit our requirements.

**PreProcessing.py**
This script is a Python script that takes in .csv files containing metadata and abundance data. These files are held within `data`. The script is currently adapted for `LiuLloydGevers1_AbundanceData_WSpecies_WMetadata_IBD.csv` but can be edited to be adapted to other datasets. The script removes rows that have no abundance data (sum to 0), normalizes the row to 1 in order to calculate relative abundance for each sample, and then columns where a study did not have that taxa found in any of its samples, those elements for that study were changed to NaN's instead of 0's. 





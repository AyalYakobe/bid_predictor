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

## 3. Scripts

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







#!/bin/bash

# Directory where the FASTQ files will be stored
output_dir="fastq_files"
mkdir -p "$output_dir"

# File containing the list of accession numbers
accession_file="accession_list.txt"

# Log file for the processing
log_file="fastq_download.log"
echo "Start processing at $(date)" > "$log_file"

# Check if the accession file exists
if [ ! -f "$accession_file" ]; then
    echo "Accession file does not exist: $accession_file"
    echo "Accession file does not exist: $accession_file" >> "$log_file"
    exit 1
fi

# Loop through each accession number in the file
while IFS= read -r accession; do
    if [ -z "$accession" ]; then
        continue # Skip empty lines
    fi

    echo "Processing $accession..."
    echo "Processing $accession..." >> "$log_file"

    # Run fastq-dump to split the files by read (assuming paired-end data)
    fastq-dump --split-files --gzip --outdir "$output_dir" "$accession"
    
    # Check if fastq-dump was successful
    if [ $? -eq 0 ]; then
        echo "$accession completed successfully."
        echo "$accession completed successfully." >> "$log_file"
    else
        echo "Error processing $accession."
        echo "Error processing $accession." >> "$log_file"
    fi
done < "$accession_file"

echo "All accessions processed."
echo "All accessions processed at $(date)" >> "$log_file"


#!/bin/bash

#This script counts the length of the different reads of every fastq.gz file in a directory

directori="./"

# Read all fastq.gz files in the directory
for archivo in "$directori"/*.fastq.gz; do
    # Verify it is a file, not a directory
    if [ -f "$archivo" ]; then
	echo "Reading $(basename "$archivo"), it may take a while."    
        # Count the length of the reads and their frequencies
        frecuencias=$(zcat "$archivo" | awk 'NR % 4 == 2' | awk '{print length}' | sort -n | uniq -c | awk '{print $2, $1}')
        # Keep the result into a file
	echo "$frecuencias" >> "Leng_$(basename "$archivo" .fastq.gz).txt"
	echo "Frequencies saved in Leng_$(basename "$archivo" .fastq.gz).txt"
    fi
done

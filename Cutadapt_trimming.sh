#!/bin/bash

#Laod module bcbio, please make sure you have this tool installed in your system
module load bcbio

# Path of files to process
directorio="./"

# Output path for procesed files 
dir_salida="$(pwd)/Processed_files/"

# Make the output directory
mkdir -p "$dir_salida"

# Process all the files in the directory
for archivo in "$directorio"/*.fastq.gz; do
    # Verify if it is a file and not a directory
    if [ -f "$archivo" ]; then
        nombre_archivo=$(basename "$archivo")
        echo "Processing $nombre_archivo, this will take quite long."

        # Perform adaptor trimming and generata files
        cutadapt -a ADAPTER=TGGAATTCTCGGGTGCCAAGG --minimum-length=17 --output="$dir_salida/sRNA_$nombre_archivo" --too-short-output="$dir_salida/tooshort_$nombre_archivo" --untrimmed-output="$dir_salida/untrimm_$nombre_archivo" "$archivo"
    fi
done
module unload bcbio
echo "Processing finished!"

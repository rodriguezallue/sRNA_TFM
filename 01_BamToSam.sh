#!/bin/bash

#This Script converts the BAM files in all subdirectories of the /Results directory into SAM files

#Laod module bcbio, please make sure you have this tool installed in your system
module load bcbio

madre="./Results"

# Looking for bam files in all subdirectories
find "$madre" -type f -name "*.bam" | while read -r archivo_bam; do
    # Collect the name of the file
    nombre_sin_extension=$(basename "$archivo_bam" .bam)
    
    # Convert file to sam format
    echo "Converting $archivo_bam to SAM..."
    samtools view -h $archivo_bam > "$(dirname "$archivo_bam")/$nombre_sin_extension.sam"
    echo "Conversion done for $archivo_bam"
done

# Unload bcbio module
module unload bcbio
echo "All conversions done!"

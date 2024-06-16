#!/bin/bash

# Define searching function
buscar_secuencia_en_sam() {
    local secuencia="$1"
    local archivo_sam="$2"
    grep -w "$secuencia" "$archivo_sam" | cut -f1
}
# Paths and file names
madre_sam="./FilesSeqs_SAM"
madre_fastq="./Samples"
archivoCSV="FirstSeqs.txt"
resultado="LlistaSeqs_FS.txt"

# Reading seqs in CSV file
mapfile -t secuencias < <(awk 'NR > 1 {print $1}' "$archivoCSV")

# Iteration on sequences
for secuencia in "${secuencias[@]}"; do
    echo "Sequence: $secuencia"
    # Define lists
    lista_identificadores="" 
    lista_fastq=""
    
    #Explore SAM files
    for archivo_sam in "$madre_sam"/*.sam; do
        identificadores=$(buscar_secuencia_en_sam "$secuencia" "$archivo_sam")
        lista_identificadores+="$(printf "%s\n" "$identificadores")"
        if [ -n "$identificadores" ]; then
            lista_identificadores+=$'\n'  # add a \n only if there are IDs in SAM
        fi
        nombre_base_sam=$(basename "$archivo_sam")
        # Generate the name of FastQ file  based on SAM name
        if [[ $nombre_base_sam =~ _MGloc_mapped_sort\.sam$ ]]; then
           nombre_base_fastq=$(echo "$nombre_base_sam" | sed 's/_MGloc_mapped_sort\.sam$/.fastq/')
           elif [[ $nombre_base_sam =~ _PGloc_mapped_sort\.sam$ ]]; then
              nombre_base_fastq=$(echo "$nombre_base_sam" | sed 's/_PGloc_mapped_sort\.sam$/.fastq/')
           else
              echo "Nombre de archivo SAM no reconocido: $nombre_base_sam"
              continue
        fi
        archivo_fastq="$madre_fastq/$nombre_base_fastq"
        # Get the seqs in thefastq file only if IDs are not empty
        if [ -n "$identificadores" ]; then
            while IFS= read -r id; do
                secuencia=$(grep -F -A1 "$id" "$archivo_fastq" | grep -v '^--$' | sed -n '2p')
                lista_fastq+=$(echo "$secuencia" | tr '\n' ',')
            done <<< "$identificadores"
        fi
    done
    # Keep original seqs and lists of fasts seqs in a file
    echo -e "$secuencia\t${lista_fastq// /,}" >> "$resultado"
done
echo "Done!"

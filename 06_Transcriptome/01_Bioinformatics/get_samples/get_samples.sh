#!/bin/bash

# #raw data path
# raw_data="${RNASEQ_RAW_DIR}"  # see sample_renaming.R

#new folder to save files
mkdir 00_rawdat

#loop over sample list
while read line; do

  #file_name=$(awk '{print $1}') #file name
  file_path=$(echo ${line} | awk '{print $1}') #file location
  file_name=$(basename ${file_path})
  new_name=$(echo ${line} | awk '{print $2}') #new file name
  new_path=00_rawdat/${new_name} #new path and name to save file

  #copy file and rename
  cp ${file_path} ${new_path}
  echo oldname: ${file_name}   new name: ${new_name}

done < get_samples.tab #list containing sample paths and new names

#list with all sample names
ls 00_rawdat | sed -E 's/.L00(1|2)_R(1|2).fastq.gz//' | sort -u > ../samples.txt

#!/bin/bash
###############################################
## Discription: Replace the repetitive command line with a SCRIPT.
##                    -- Assist : Frame by Frame MD trajectory Analysis (mkvmd_xxx.sh)
## Author:
##          Shang Chun, Oct 2022
##
##
################################################

## Content of analysis
FEATURE="contactnum"
## Prefix of file name
PREFIX="TCR_MHCPepContact"
## Collection of structure-folders
FOLDERS=(01_7QPJ 02_7PDW 03_7PBC 04_7RM4 05_7RTR 06_7N6E 07_7N1F 08_7N1E 09_6RSY 10_4WUU 12_6RP9 13_6R2L 15_6VMA 16_6VMC 19_6VM9)
## Folder list 1
FOLD1="01_7QPJ"

for folder in ${FOLDERS[*]};
do 
	if [ "$folder" == "$FOLD1" ]; then
		cd "$FOLD1"/run/output/Analysis/
		bash ../../../../../Analysis/mkvmd_"$PREFIX".sh ../../../pdb2namd/vmd_solvate/ionized.psf ../md.dcd "$FEATURE".out
	else 
		cd ../../../../"$folder"/run/output/Analysis
		bash ../../../../../Analysis/mkvmd_"$PREFIX".sh ../../../pdb2namd/vmd_solvate/ionized.psf ../md.dcd "$FEATURE".out
	fi
done



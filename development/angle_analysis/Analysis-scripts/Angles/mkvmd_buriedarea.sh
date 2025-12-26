#!/bin/bash
#############################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author:
##          Shang Chun, May 2022
##
##
############################################# 

# Selections: pMHC/TCR
SEL1="segname PROA or segname PROB or segname PROC"
SEL2="segname PROD or segname PROE"

PSF="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PSF] [TRJ] [OUTPUT]\n      By default, the selections are:\n      Selection 1: $SEL1\n      Selection 2: $SEL2\n"; exit 1; }

if [ ! -f $PSF ]; then
	echo -e "$PSF \nStructure not found!"
	exit 0
fi

if [ ! -f $TRJ ]; then
	echo -e "$TRJ \nTrajectory not found!"
	exit 0
fi

rm ${OUTPUT}

cat > BuriedArea.tcl << EOF
proc measureBuried { nf } {
	set selpmhc [atomselect top "$SEL1" frame \$nf]
	set seltcr [atomselect top "$SEL2" frame \$nf]
	set selcpx [atomselect top "$SEL1 or $SEL2" frame \$nf]

	return [format "%.2f" [expr ([measure sasa 1.4 \$selpmhc] + [measure sasa 1.4 \$seltcr] - [measure sasa 1.4 \$selcpx]) / 2]]
}

# /---------------------/
# /      Main Body      /
# /---------------------/

mol new $PSF waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something ..."
for {set nf 0} {\$nf < \$total_frame} {incr nf} {

	set outf [open ${OUTPUT} "a"]


	set out_line [format "%d" \$nf]

	puts "Calculating buried surface areas ..."
	lappend out_line [measureBuried \$nf]

	puts "Writing to outfile ..."
	puts \$outf "\$out_line"

	close \$outf
}

quit
EOF

vmd -dispdev text -e BuriedArea.tcl

#!/bin/bash
#############################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author:
##          Shang Chun, May 2022
##
##
############################################# 

# Selections:
## MHC: the best-fit straight line through the Cα atoms from the two MHC helices.
##      Class I: Cα atoms A50–A86 and A140–A176
##      Class II: A46–78, B54–64, and B67–91
SEL1="segname PROA and resid 50 to 86 140 to 176"


## TCR: mainly TCR ALPHA CHAIN V REGION
##      Alpha:
##      Beta:
SEL2="segname PROD and resid 3 to 115"
SEL3="segname PROE and resid 3 to 115"

PSF="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PSF] [TRJ] [OUTPUT]\n      By default, the selections are:\n      Selection 1: $SEL1\n      Selection 2: $SEL2\n      Selection 3: $SEL3\n"; exit 1; }

if [ ! -f $PSF ]; then
	echo -e "$PSF \nStructure not found!"
	exit 0
fi

if [ ! -f $TRJ ]; then
	echo -e "$TRJ \nTrajectory not found!"
	exit 0
fi

rm ${OUTPUT}

cat > TiltAngle.tcl << EOF
# Calc angles between two vectors
proc angle { a b } {
	# get Pi
	global M_PI

	# Angle between two vectors
	set cosine [expr [vecdot \$a \$b] / ( [veclength \$a] * [veclength \$b])]
	return [expr acos(\$cosine)*(180.0/\$M_PI)]
}

# Calc the Tilt Angle between MHC&TCR, corresponds to frame $nf
proc measureTCRTiltAngle {nf} {
	set selmhc [atomselect top "$SEL1" frame \$nf]
	set seltcr [atomselect top "$SEL2 or $SEL3" frame \$nf]

	set mhcnorm [lindex [Orient::calc_principalaxes \$selmhc] 0]
	set tcrnorm [lindex [Orient::calc_principalaxes \$seltcr] 1]

	if {[angle \$mhcnorm \$tcrnorm] > 90} {
		set mhcnorm [vecinvert \$mhcnorm]
	}

	\$selmhc delete
	\$seltcr delete

	return [format "%.2f" [angle \$mhcnorm \$tcrnorm]]
}

# /---------------------/
# /      Main Body      /
# /---------------------/

# Load packages for principle-axis calculation
package require Orient
namespace import Orient::orient

mol new $PSF waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something ..."
for {set nf 0} {\$nf < \$total_frame} {incr nf} {

	set outf [open ${OUTPUT} "a"]

	set out_line [format "%d" \$nf]
	lappend out_line [measureTCRTiltAngle \$nf]

	puts \$outf "\$out_line"

	close \$outf
}

quit
EOF

vmd -dispdev text -e TiltAngle.tcl

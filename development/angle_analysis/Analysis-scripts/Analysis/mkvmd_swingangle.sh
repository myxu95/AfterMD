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


## TCR: a vector between the centroids(质心) of the conserved sulfuratoms in the light&heavy chains.
## measure center：使用给定的权重返回选择中原子的几何中心
##      Alpha/Light: L22-L90
##      Beta/Heavy: H23-H92
SEL2="segname PROD and resid 89 to 94 and resname CYS and name CA" # TCR alpha
SEL3="segname PROD and resid 20 to 25 and resname CYS and name CA" 
SEL4="segname PROE and resid 89 to 94 and resname CYS and name CA" # TCR beta
SEL5="segname PROE and resid 20 to 25 and resname CYS and name CA" 

## Pep: for referance
SEL6="segname PROC" #选不了CA？？？

PSF="$1"
TRJ="$2"
OUTPUT="$3"
[ $# -ne 3 ] && { echo -e "mkvmd> Usage: $0 [PSF] [TRJ] [OUTPUT]\n      By default, the selections are:\n      Selection 1: $SEL1\n      Selection 2: $SEL2 or $SEL3\n      Selection 3: $SEL4 or $SEL5\n     Selection 3: $SEL6"; exit 1; }

if [ ! -f $PSF ]; then
	echo -e "$PSF \nStructure not found!"
	exit 0
fi

if [ ! -f $TRJ ]; then
	echo -e "$TRJ \nTrajectory not found!"
	exit 0
fi

rm ${OUTPUT}

cat > SwingAngle.tcl << EOF
# Calc angles between two vectors
proc angle { a b } {
	# get Pi
	global M_PI

	# Angle between two vectors
	set cosine [expr [vecdot \$a \$b] / ( [veclength \$a] * [veclength \$b])]
	return [expr acos(\$cosine)*(180.0/\$M_PI)]
}

# Calc the Swing Angle between MHC&TCR, corresponds to frame $nf
proc measureTCRSwingAngle {nf} {
	set selmhc [atomselect top "$SEL1" frame \$nf]

	## TCR alpha / segname PROD
	set seltcra [atomselect top "$SEL2 or $SEL3" frame \$nf]
	## TCR beta / segname PROE
	set seltcrb [atomselect top "$SEL4 or $SEL5" frame \$nf]
	set seltcrab [atomselect top "$SEL2 or $SEL3 or $SEL4 or $SEL5" frame \$nf] 

	set selpep [atomselect top "$SEL6" frame \$nf]
	

	# 投影
	############################################################
	##  \  TCR上的直线  \ MHC上投影平面的法线 \ MHC上的角度另一边 \
	##  / TCR二硫键连线 /   MHC主轴方向"3"   /  MHC主轴方向"1"  /
	set tcr [vecsub [measure center \$seltcrb weight mass] [measure center \$seltcra weight mass]]  
	## vecsub $end $start ## segD指向segE ## 轻链指向重链 ## alpha指向beta
	set mhc [lindex [Orient::calc_principalaxes \$selmhc] 2]
	set agl [lindex [Orient::calc_principalaxes \$selmhc] 0]
	############################################################

	## 计算投影向量
	set p1 [vecnorm \$mhc]
	set p2 [vecdot \$tcr \$p1]
	set p3 [vecscale \$p1 \$p2]
	set swi [vecsub \$tcr \$p3]	

	## 定义参考“竖直”方向 ## 从Pep指向TCR
	set ref [vecsub [measure center \$seltcrab weight mass] [measure center \$selpep weight mass]]

	## 使agl指向为“上”
	if { [angle \$ref \$agl] > 90 } {
		set agl [vecinvert \$agl]
	}

	set swiag [angle \$swi \$agl]

	\$selmhc delete
	\$seltcra delete
	\$seltcrb delete
	\$seltcrab delete
	\$selpep delete

	return [format "%.2f" \$swiag]
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
	lappend out_line [measureTCRSwingAngle \$nf]

	puts \$outf "\$out_line"

	close \$outf
}

quit
EOF

vmd -dispdev text -e SwingAngle.tcl

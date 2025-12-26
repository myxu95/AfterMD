#!/bin/bash
#############################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author:
##          Shang Chun, Jun 2022
##
##
############################################# 
## SCRIPT_MD：TCR两链上与Pep的接触氨基酸对（接触原子）统计打印

# Selection for alpha-chain
SEL1="segname PROD"
# Selection for beta-chain
SEL2="segname PROE"
# Selection for Pep
SEL3="segname PROC"

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

rm ${OUTPUT}_a_pep ${OUTPUT}_b_pep
## pair=OUTPUT
cat > tcl << EOF
proc PepContactPair { nf sel_stats outf } {
    puts "\$sel_stats -- Pep contact pairs : frame \$nf"
    puts \$outf "\$sel_stats -- Pep contact pairs : frame \$nf"

    set slist []
    set sel [atomselect top "segname PROC"]
    foreach res_i [lsort -unique -integer [\$sel get resid]] {

        set sel_icontact [atomselect top "noh and \$sel_stats and within 4 of (segname PROC and resid \$res_i)" frame \$nf]
        set sel_i [atomselect top "name CA and segname PROC and resid \$res_i"]
        set list [concat \$res_i [\$sel_i get resname] [\$sel_icontact get resid] [\$sel_icontact get resname]]
        puts \$outf \$list
        #lappend slist \$list \$\n
        #lappend slist \$res_i [\$sel_i get resname] [\$sel_icontact get resid] [\$sel_icontact get resname] \$\n
    }
    #return \$slist
}

##### 测试用
#set sel [atomselect top "segname PROC"]
#foreach res_i [lsort -unique -integer [$sel get resid]] {
#    set sel_icontact [atomselect top "noh and segname PROD and within 4 of (segname PROC and resid $res_i)"]
#    set sel_i [atomselect top "segname PROC and resid $res_i and name CA"]
#    set list [concat $res_i [$sel_i get resname] [$sel_icontact get resid] [$sel_icontact get resname]]
#    puts $list
#}    
#puts " "

# /------------------/
# /     Main Body    /
# /------------------/

mol new $PSF waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nf 0} {\$nf < \$total_frame} {incr nf} {
    set outf1 [open ${OUTPUT}_a_pep "a"]
    set outf2 [open ${OUTPUT}_b_pep "a"]

    # Write TIME at the very first of a line
    #set out_line [format "%d" \$nf]
    puts \$outf1 "[PepContactPair \$nf "$SEL1" \$outf1]"
    puts \$outf2 "[PepContactPair \$nf "$SEL2" \$outf2]"

    close \$outf1
    close \$outf2

}
quit
EOF

vmd -dispdev text -e tcl
#rm tcl

#!/bin/bash
#############################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author:
##          Shang Chun, Jun 2022
##
##
############################################# 
## SCRIPT_MD：TCR两链上与MHC/Pep的接触数目 （“前置条件“）

# Selection for alpha-chain
SEL1="segname PROD"
# Selection for beta-chain
SEL2="segname PROE"
# Selection for TCR
SEL3="segname PROD PROE"
# Selection for MHC
SEL4="segname PROA"
# Selection for Pep
SEL5="segname PROC"

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

rm ${OUTPUT}_mhc ${OUTPUT}_pep
## contactnum.out=OUTPUT
cat > tcl << EOF
proc CountTCR_Contact { nf sel_ref } {
    ## 注意list要用{}
    set cal_list {"$SEL1" "$SEL2" "$SEL3"}
    set dstlist {3.5 4 4.5}

    set numlist []
    foreach sel_cal \$cal_list {
        set count 0
        foreach dst \$dstlist {
            ## 删掉same residue as, Kevin, Jun 2022
            set selheavy [atomselect top "noh and \$sel_cal and within \$dst of \$sel_ref" frame \$nf]
            incr count [\$selheavy num]
            \$selheavy delete
        }
        lappend numlist [expr \$count / [llength \$dstlist]]
    }
    set total [lindex \$numlist 2] 
    puts "Total number of \$sel_ref nearby-residue-atoms: \$total"
    return \$numlist

}
    


# /------------------/
# /     Main Body    /
# /------------------/

mol new $PSF waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

puts "mkvmd> Computing something..."
for {set nf 0} {\$nf < \$total_frame} {incr nf} {
        set outf1 [open ${OUTPUT}_mhc "a"]
    set outf2 [open ${OUTPUT}_pep "a"]

    # Write TIME at the very first of a line
    set out_line [format "%d" \$nf]

    puts "Writing to outfile ..."
    puts \$outf1 "\$out_line [CountTCR_Contact \$nf "$SEL4"]"
    puts \$outf2 "\$out_line [CountTCR_Contact \$nf "$SEL5"]"

    close \$outf1
    close \$outf2

}
quit
EOF

vmd -dispdev text -e tcl
rm tcl

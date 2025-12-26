#!/bin/bash
#############################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author:
##          Shang Chun, Jul 2022
##
##
############################################# 
## SCRIPT_MD：TCR alpha 链上与 Pep 的接触氨基酸对（接触原子）频次/频率统计

# Selection for alpha-chain
SEL1="segname PROD"
# Selection for beta-chain
SEL2="segname PROE"
# Selection for Pep
SEGNAME="PROC"
SEL3="segname $SEGNAME"

PSF="$1"
TRJ="$2"
OUTPUT="$3"
PDB="$4"
[ $# -ne 4 ] && { echo -e "mkvmd> Usage: $0 [PSF] [TRJ] [OUTPUT] [PDB]\n      By default, the selections are:\n      Selection 1: $SEL1\n      Selection 2: $SEL2\n"; exit 1; }

if [ ! -f $PSF ]; then
    echo -e "$PSF \nStructure not found!"
    exit 0
fi

if [ ! -f $TRJ ]; then
    echo -e "$TRJ \nTrajectory not found!"
    exit 0
fi

if [ ! -f $PDB ]; then
    echo -e "$PDB \nStructure not found!"
    exit 0
fi

length_of_peptide=`grep $SEGNAME $PDB | grep "CA" -c`

rm ${OUTPUT}_a_pep
## map_CDR=OUTPUT

####### 格式： pep_posi:\$res_i        posi_resname:[\$sel_i get resname]
####### contact_residue_id:[\$sel_icontact get resid]      contact_residue_name:[\$sel_icontact get resname]

## 创建字典
## 例：set colours [dict create colour1 "black" colour2 "white"]

# set contact_res_name [lindex [\$sel_icontactj get resname] 0]
# set addnum [llength [\$sel_icontactj get resname]] 名称一样会传出几个值？？

## global：
## 全局变量NumMap在过程proc中被访问；在过程中对NumMap的改变会直接反映到全局上.

cat > tcl << EOF
set resnameshort [dict create "ARG" "R" "HIS" "H" "LYS" "K" "ASP" "D" "GLU" "E" "SER" "S" "THR" "T" "ASN" "N" "GLN" "Q" "CYS" "C" "GLY" "G" "TYR" "Y" "ALA" "A" "VAL" "V" "ILE" "I" "LEU" "L" "MET" "M" "PHE" "F" "TRP" "W" "PRO" "P"]

proc PepContactMap { nf sel_stats outf } {
  global resnameshort
  global NumMap
    puts "\$sel_stats -- Pep contact Map : frame \$nf"
    puts \$outf "\$sel_stats -- Pep contact Map : frame \$nf"
    puts \$outf [format "%-3s%-12s%-6s%-28s" "-" CDR1 CDR2 CDR3]
    set sel [atomselect top "$SEL3"]
    set resid_list [lsort -unique -integer [\$sel get resid]] 
    foreach res_i \$resid_list {
        set sel_i [atomselect top "name CA and $SEL3 and resid \$res_i"]
        set sel_icontact [atomselect top "noh and \$sel_stats and within 4 of ($SEL3 and resid \$res_i)" frame \$nf]
        set contactresid_list [lsort -unique -integer [\$sel_icontact get resid]]
        set contactresname_list [] 
        foreach contact_res_id \$contactresid_list {
          set sel_j [atomselect top "name CA and \$sel_stats and resid \$contact_res_id"]
          set j_resname [\$sel_j get resname]
          lappend contactresname_list [dict get \$resnameshort "\$j_resname"]
        }
        set print_ilist [format "%-2d" \$res_i]
        set cdr1_list []
        set cdr2_list []
        set cdr3_list []
        foreach id \$contactresid_list name \$contactresname_list {
          if { (25<=\$id)&&(\$id<=35) } {
            lappend cdr1_list [format "%2d%1s" \$id \$name]
          } elseif { (45<=\$id)&&(\$id<=55) } {
            lappend cdr2_list [format "%2d%1s" \$id \$name]
          } else {
            lappend cdr3_list [format "%2d%1s" \$id \$name]
          }
        }
        lappend print_ilist [format "%-12s%-6s%-28s" \$cdr1_list \$cdr2_list \$cdr3_list]
        puts \$outf \$print_ilist
    }
}


# /------------------/
# /     Main Body    /
# /------------------/

mol new $PSF waitfor all
mol addfile $TRJ waitfor all
set total_frame [molinfo top get numframes]

## “数组”初始化
## set NumMap(\$i,\$j) 0
for {set i 1} {\$i <= $length_of_peptide} {incr i} {
    for {set j 1} {\$j <= 20} {incr j} {
      set NumMap(\$i,\$j) 0
    }
}

puts "mkvmd> Computing something..."
for {set nf 0} {\$nf < \$total_frame} {incr nf} {
    set outf1 [open ${OUTPUT}_a_pep "a"]

    puts \$outf1 "[PepContactMap \$nf "$SEL1" \$outf1]"
    close \$outf1

}
quit
EOF

vmd -dispdev text -e tcl
rm tcl

#!/bin/bash
#############################################
## Description: Main body of a "boss run" VMD/tcl measure script
## Author:
##          Shang Chun, Jul 2022
##
##
############################################# 
## SCRIPT_MD：TCR 上与 Pep 的"反应"氨基酸(种类)对(侧链间接触)频次/频率统计，Shang Chun, Aug 2022

# Selection for alpha-chain
SEL1="segname PROD"
# Selection for beta-chain
SEL2="segname PROE"
# Selection for Pep
SEGNAME="PROC"
SEL3="segname $SEGNAME"
###不区分alpha/beta链，SC, Jan 2023
SEL0="segname PROD PROE"

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

rm ${OUTPUT}_pep
## map_shannon=OUTPUT

####### 格式： pep_posi:\$res_i        posi_resname:[\$sel_i get resname]
####### contact_residue_id:[\$sel_icontact get resid]      contact_residue_name:[\$sel_icontact get resname]

## 创建字典
## 例：set colours [dict create colour1 "black" colour2 "white"]

# set contact_res_name [lindex [\$sel_icontactj get resname] 0]
# set addnum [llength [\$sel_icontactj get resname]] 名称一样会传出几个值？？

## global：
## 全局变量NumMap在过程proc中被访问；在过程中对NumMap的改变会直接反映到全局上.

#set formatStr {%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s}
#puts \$outf [format \$formatStr "\\" "R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "G" "Y" "A" "V" "I" "L" "M" "F" "W" "P"]
# set addnum [expr double(\$NumMap(\$i,\$j)) / (\$nf + 1)]

#puts \$outf [format \$formatStr2 \$ilist]
#set formatStr2 {%4s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f}

cat > tcl << EOF
set residnum [dict create "ARG" "1" "HSD" "2" "LYS" "3" "ASP" "4" "GLU" "5" "SER" "6" "THR" "7" "ASN" "8" "GLN" "9" "CYS" "10" "GLY" "11" "TYR" "12" "ALA" "13" "VAL" "14" "ILE" "15" "LEU" "16" "MET" "17" "PHE" "18" "TRP" "19" "PRO" "20"]
set reslist {"R" "H" "K" "D" "E" "S" "T" "N" "Q" "C" "G" "Y" "A" "V" "I" "L" "M" "F" "W" "P"}
set resshort [dict create "ARG" "R" "HSD" "H" "LYS" "K" "ASP" "D" "GLU" "E" "SER" "S" "THR" "T" "ASN" "N" "GLN" "Q" "CYS" "C" "GLY" "G" "TYR" "Y" "ALA" "A" "VAL" "V" "ILE" "I" "LEU" "L" "MET" "M" "PHE" "F" "TRP" "W" "PRO" "P"]

proc PepContactMap { nf sel_stats outf } {
  package require csv
  global reslist
  global residnum
  global resshort
  global NumMap
    puts "\$sel_stats -- Pep contact Map : frame \$nf"
    #puts \$outf "\$sel_stats -- Pep contact Map : frame \$nf"
    set sel [atomselect top "$SEL3"]
    foreach res_i [lsort -unique -integer [\$sel get resid]] {
        set seq_\$res_i [] 
        set sel_i [atomselect top "name CA and $SEL3 and resid \$res_i"]
        set i_res_name [\$sel_i get resname]
        append seq_\$res_i [dict get \$resshort "\$i_res_name"]
        set sel_icontact [atomselect top "noh and \$sel_stats and sidechain and within 4 of (noh and $SEL3 and resid \$res_i and sidechain)" frame \$nf]
        foreach contact_res_id [lsort -unique -integer [\$sel_icontact get resid]] {
          set sel_j [atomselect top "name CA and \$sel_stats and resid \$contact_res_id"]
          set contact_res_name [\$sel_j get resname]
          set contact_res_num [dict get \$residnum "\$contact_res_name"]
          set newcount [expr \$NumMap(\$res_i,\$contact_res_num) + 1]
          set NumMap(\$res_i,\$contact_res_num) \$newcount
        }
    }
    if { \$nf==99 } { 
    set formatStr1 {%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s}
    set formatStr2 {%4s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f}
    puts \$outf [format \$formatStr1 "--" 1 2 3 4 5 6 7 8 9]
    foreach s1 \$seq_1 s2 \$seq_2 s3 \$seq_3 s4 \$seq_4 s5 \$seq_5 s6 \$seq_6 s7 \$seq_7 s8 \$seq_8 s9 \$seq_9 {
      puts \$outf [format \$formatStr1 "--" \$s1 \$s2 \$s3 \$s4 \$s5 \$s6 \$s7 \$s8 \$s9]
    }
    for {set i 1} {\$i <= $length_of_peptide} {incr i} {
      set list_\$i []
      for {set j 1} {\$j <= 20} {incr j} {
        set tcl_precision 2
        set addnum [expr \$NumMap(\$i,\$j) / (\$nf.0 + 1.0)]
        lappend list_\$i \$addnum
      }
    }
    foreach resitem \$reslist t1 \$list_1 t2 \$list_2 t3 \$list_3 t4 \$list_4 t5 \$list_5 t6 \$list_6 t7 \$list_7 t8 \$list_8 t9 \$list_9 {
      puts \$outf [csv::join [list \$resitem \$t1 \$t2 \$t3 \$t4 \$t5 \$t6 \$t7 \$t8 \$t9]]
    }
    } else {
      ##no opration
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
    set outf [open ${OUTPUT}.csv w]

    puts \$outf "[PepContactMap \$nf "$SEL1" \$outf]"
    puts \$outf "[PepContactMap \$nf "$SEL2" \$outf]"
    close \$outf

}
quit
EOF

vmd -dispdev text -e tcl
#rm tcl

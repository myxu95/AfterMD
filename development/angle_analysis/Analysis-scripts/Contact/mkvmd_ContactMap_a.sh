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
## map=OUTPUT

####### 格式： pep_posi:\$res_i        posi_resname:[\$sel_i get resname]
####### contact_residue_id:[\$sel_icontact get resid]      contact_residue_name:[\$sel_icontact get resname]

## 创建字典
## 例：set colours [dict create colour1 "black" colour2 "white"]

# set contact_res_name [lindex [\$sel_icontactj get resname] 0]
# set addnum [llength [\$sel_icontactj get resname]] 名称一样会传出几个值？？

## global：
## 全局变量NumMap在过程proc中被访问；在过程中对NumMap的改变会直接反映到全局上.

cat > tcl << EOF
set residnum [dict create "ARG" "1" "HIS" "2" "LYS" "3" "ASP" "4" "GLU" "5" "SER" "6" "THR" "7" "ASN" "8" "GLN" "9" "CYS" "10" "GLY" "11" "TYR" "12" "ALA" "13" "VAL" "14" "ILE" "15" "LEU" "16" "MET" "17" "PHE" "18" "TRP" "19" "PRO" "20"]

proc parray1 {a outf {pattern *}} {
    upvar 1 \$a array
    if {![array exists array]} {
      error "\"\$a\" isn't an array"
    }
    set maxl 0
    foreach name [lsort [array names array \$pattern]] {
      if {[string length \$name] > \$maxl} {
        set maxl [string length \$name]
      }
    }
    set maxl [expr {\$maxl + [string length \$a] + 2}]
    foreach name [lsort [array names array \$pattern]] {
      set nameString [format %s(%s) \$a \$name]
      puts \$outf [format "%-*s = %s" \$maxl \$nameString \$array(\$name)]
    }
}

proc PepContactMap { nf sel_stats outf } {
  global residnum
  global NumMap
    puts "\$sel_stats -- Pep contact Map : frame \$nf"
    puts \$outf "\$sel_stats -- Pep contact Map : frame \$nf"
    set sel [atomselect top "$SEL3"]
    foreach res_i [lsort -unique -integer [\$sel get resid]] {
        set sel_i [atomselect top "name CA and $SEL3 and resid \$res_i"]
        set sel_icontact [atomselect top "noh and \$sel_stats and within 4 of ($SEL3 and resid \$res_i)" frame \$nf]
        foreach contact_res_id [lsort -unique -integer [\$sel_icontact get resid]] {
          set sel_icontactj [atomselect top "noh and \$sel_stats and resid \$contact_res_id and within 4 of ($SEL3 and resid \$res_i)" frame \$nf]
          set sel_j [atomselect top "name CA and \$sel_stats and resid \$contact_res_id"]
          set contact_res_name [\$sel_j get resname]
          set contact_res_num [dict get \$residnum "\$contact_res_name"]
          set addcount [llength [\$sel_icontactj get name]]
          set newcount [expr \$NumMap(\$res_i,\$contact_res_num) + \$addcount]
          set NumMap(\$res_i,\$contact_res_num) \$newcount
        }
    }
    parray1 NumMap \$outf
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

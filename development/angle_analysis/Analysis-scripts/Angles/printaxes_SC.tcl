## 用于在个人MAC上，画出加载VMD结构的主轴 
package require Orient 
namespace import Orient::orient

set SEL1 "chain A and resid 50 to 86 140 to 176"
set SEL2 "chain D and resid 1 to 119 or chain E and resid 1 to 107"

set selmhc [atomselect top "$SEL1"]
set seltcr [atomselect top "$SEL2"]

## 画主轴

draw principalaxes $seltcr  
#puts "TCR axes are $I2"
draw principalaxes $selmhc 
#puts "MHC axes are $I1"

#set mhc [Orient::calc_principalaxes $selmhc] ## list形式
#set tcr [Orient::calc_principalaxes $seltcr] ## list形式

set mhc [lindex [Orient::calc_principalaxes $selmhc] 1]
set agl [lindex [Orient::calc_principalaxes $selmhc] 0]
set tcr [lindex [Orient::calc_principalaxes $seltcr] 1]

set p1 [vecnorm $mhc]
set p2 [vecdot $tcr $p1]
set p3 [vecscale $p1 $p2]
set tst [vecsub $tcr $p3]
# puts "Where wrong"
# set tst [vecsub $tcr [vecscale [vecnorm $mhc] [vecdot $tcr [vecnorm $mhc]]]]


## 求夹角函数
proc angle { a b } {
	# get Pi
	global M_PI

	# Angle between two vectors
	set cosine [expr [vecdot $a $b] / ( [veclength $a] * [veclength $b])]
	return [expr acos($cosine)*(180.0/$M_PI)]
}


## Twist投影向量
set tstag [angle $tst $agl]
if { $tstag > 90} {
		set tstag [expr (180 - $tstag)]
	}
puts "Twist Angle is $tstag"



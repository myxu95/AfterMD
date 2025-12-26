package provide Orient 1.0
package require La

## import：输入；export：输出

## 名字空间是一个命令和变量的集合，通过其封装来保证它们不会影响其它名字空间的命令和变量
## namespace eval允许用户创建一个新的namespace

## 注："::"可以当做linux系统里的文件夹"/"来理解；即绝对/相对路径
namespace eval ::Orient:: {
    namespace export orient
}

# package require Orient
# namespace import Orient::orient
# ... load your molecules and make a selection ...
#
# set I [draw principalaxes $sel]           <--- show/calc the principal axes
# set A [orient $sel [lindex $I 2] {0 0 1}] <--- rotate axis 2 to match Z
# $sel move $A
# set I [draw principalaxes $sel]           <--- recalc principal axes to check
# set A [orient $sel [lindex $I 1] {0 1 0}] <--- rotate axis 1 to match Y
# $sel move $A
# set I [draw principalaxes $sel]           <--- recalc principal axes to check#

## 算出体系的质心坐标
proc Orient::sel_com { sel weights } {
    set x [ $sel get x ]
    set y [ $sel get y ]
    set z [ $sel get z ]
    set m $weights
    
    set comx 0
    set comy 0
    set comz 0
    set totalm 0
    foreach xx $x yy $y zz $z mm $m {
        # use the abs of the weights
        set mm [expr abs($mm)]
	set comx [ expr "$comx + $xx*$mm" ]
	set comy [ expr "$comy + $yy*$mm" ]
	set comz [ expr "$comz + $zz*$mm" ]
	set totalm [ expr "$totalm + $mm" ]
    }
    set comx [ expr "$comx / $totalm" ]    ## x方向的$sel体系质心
    set comy [ expr "$comy / $totalm" ]
    set comz [ expr "$comz / $totalm" ]
    puts "Total weight: $totalm"
    return [list $comx $comy $comz]
}

## 理论力学里，转动惯量为一二阶张量(3*3的矩阵)
## COM：Center Of Mass，质心 
proc Orient::sel_it { sel COM weights} {
    set x [ $sel get x ]
    set y [ $sel get y ]
    set z [ $sel get z ]
    set m $weights

    # compute I
    ## I=mr^2，I转动惯量，r为距离转动轴线的垂直距离
    set Ixx 0
    set Ixy 0
    set Ixz 0
    set Iyy 0
    set Iyz 0
    set Izz 0
    foreach xx $x yy $y zz $z mm $m {
        # use the abs of the weights
        set mm [expr abs($mm)]    ## absolute value：绝对值
        
        # subtract the COM   ## subtract：减去
        set xx [expr $xx - [lindex $COM 0]]
        set yy [expr $yy - [lindex $COM 1]]
        set zz [expr $zz - [lindex $COM 2]]

        set rr [expr $xx + $yy + $zz]

        set Ixx [expr $Ixx + $mm*($yy*$yy+$zz*$zz)]
        set Ixy [expr $Ixy - $mm*($xx*$yy)]
        set Ixz [expr $Ixz - $mm*($xx*$zz)]
        set Iyy [expr $Iyy + $mm*($xx*$xx+$zz*$zz)]
        set Iyz [expr $Iyz - $mm*($yy*$zz)]
        set Izz [expr $Izz + $mm*($xx*$xx+$yy*$yy)]

    }
    
    return [list 2 3 3 $Ixx $Ixy $Ixz $Ixy $Iyy $Iyz $Ixz $Iyz $Izz]
}

## veclength：返回向量的模
## vecadd 2 5 9 10：返回参数的加和
## vecadd {1 2 3} {4 5 6} {7 8 9}：向量元素一对一相加，返回一个新向量 12 15 18
## vecsub：向量元素一对一相减
## vecscale {-5 4 -3 2} -2：对每个向量元素乘以-2

## arrow：箭头
proc vmd_draw_arrow {mol start end} {
    set scaling [expr [veclength [vecsub $end $start]]/100]
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.8 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius [expr 2*$scaling]
    puts [list cone $middle $end radius [expr 5*$scaling]]
    graphics $mol cone $middle $end radius [expr 5*$scaling]
}

## vector：矢量
proc vmd_draw_vector { mol pos val } {
    set end   [ vecadd $pos [ vecscale +1 $val ] ]
    vmd_draw_arrow $mol $pos $end
}    

# find the max of some numbers
proc Orient::max { args } {
    set maxval [lindex $args 0]
    foreach arg $args {
        if { $arg > $maxval } {
            set maxval $arg
        }
    }
    return $maxval
}

# draws the three principal axes
proc vmd_draw_principalaxes { mol sel {weights domass} } {
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    set I [Orient::calc_principalaxes $sel $weights]    ## I即主轴list
    set a1 [lindex $I 0]
    set a2 [lindex $I 1]
    set a3 [lindex $I 2]

    # find the size of the system
    set minmax [measure minmax $sel]
    set ranges [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    set scale [expr .7*[Orient::max [lindex $ranges 0] \
                             [lindex $ranges 1] \
                             [lindex $ranges 2]]]
    set scale2 [expr 1.02 * $scale]

    # draw some nice vectors
    graphics $mol delete all
    graphics $mol color yellow
    set COM [Orient::sel_com $sel $weights]
    vmd_draw_vector $mol $COM [vecscale $scale $a1]
    vmd_draw_vector $mol $COM [vecscale $scale $a2]
    vmd_draw_vector $mol $COM [vecscale $scale $a3]

    graphics $mol color white
    graphics $mol text [vecadd $COM [vecscale $scale2 $a1]] "1"
    graphics $mol text [vecadd $COM [vecscale $scale2 $a2]] "2"
    graphics $mol text [vecadd $COM [vecscale $scale2 $a3]] "3"
    
    return [list $a1 $a2 $a3]  ##返回主轴的list
}


## 区别：惯性张量是绕着一点的，而转动惯量是绕着一个轴的
## 联系：绕着x、y、z三个轴转动惯量恰好就是惯性张量的三个对角元

###################### 语句解析：La::mevsvd_br I evals
# 特征向量作为 I 的列返回(以列形式返回)
# 特征值作为向量返回
# 一个特征向量对应一个特征值

#### 主轴定理：
#### 指在几何学和线性代数中，与椭圆和双曲线的长轴和短轴有关的线。这些轴能够将椭圆和双曲线准确的描述出来，它们是正交的。
####    设A是一个n*n实对称矩阵，那么存在一个正交变量代换x=Py，
####    它将二次型x'Ax变换为不含交叉项的二次型y'Dy（即D为对角矩阵）
####    将P的诸列称为二次型x'Ax的主轴
####    以主轴作为R^n的一组单位正交基，则y是向量x在这组正交基下的坐标。

# returns the three principal axes
proc Orient::calc_principalaxes { sel {weights domass} } {
    puts "Calculating principal axes."
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    puts "Getting the center-of-mass..."
    # get the COM
    set COM [Orient::sel_com $sel $weights]
    puts "Computing the inertia tensor..."           ## inertia tensor：惯性张量
    # get the I
    set I [Orient::sel_it $sel $COM $weights]        ## I：惯性张量

    puts "Drawing the principal components..."
    La::mevsvd_br I evals                           ## La：Linear Algebra and similar math，Hume Integration Software 的线性代数包
    # now $I holds in its columns the principal axes
    set a1 "[lindex $I 3] [lindex $I 6] [lindex $I 9]"
    set a2 "[lindex $I 4] [lindex $I 7] [lindex $I 10]"
    set a3 "[lindex $I 5] [lindex $I 8] [lindex $I 11]"

    return [list $a1 $a2 $a3]
}

# rotate a selection about its COM, taking <vector1> to <vector2>
# e.g.: orient $sel [lindex $I 2] {0 0 1}
# (this aligns the third principal axis with z)

## vecnorm：返回与向量方向相同的单位向量
## veccross：向量叉乘
##  vecdot：向量点乘
proc Orient::orient { sel vector1 vector2 {weights domass}} {
    if { $weights == "domass" } {
        set weights [ $sel get mass ]
    }

    set COM [Orient::sel_com $sel $weights]

    set vec1 [vecnorm $vector1]

    set vec2 [vecnorm $vector2]

    # compute the angle and axis of rotation
    set rotvec [veccross $vec1 $vec2]
    set sine   [veclength $rotvec]
    set cosine [vecdot $vec1 $vec2]
    set angle [expr atan2($sine,$cosine)] ##已知正切，求弧度
    
    # return the rotation matrix
    return [trans center $COM axis $rotvec $angle rad]
}

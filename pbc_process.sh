#!/bin/bash

# 比较精细的GROMACS轨迹PBC处理脚本
# 作者: xmy
# 用途: 处理MD轨迹的周期性边界条件

# 设置变量
TPR_FILE="md.tpr"           # TPR文件名
XTC_FILE="md.xtc"           # 原始轨迹文件
OUTPUT_XTC="md_pbc.xtc"     # 处理后的轨迹文件
GRO_FILE="md.gro"           # 结构文件

# 检查必要文件是否存在
echo "检查输入文件..."
if [ ! -f "$TPR_FILE" ]; then
    echo "错误: 找不到TPR文件 $TPR_FILE"
    exit 1
fi

if [ ! -f "$XTC_FILE" ]; then
    echo "错误: 找不到轨迹文件 $XTC_FILE"
    exit 1
fi

echo "开始PBC处理流程..."

# 第一步: 去除跳跃 (Remove jumps)
# 修正由于PBC导致的原子/分子跳跃
echo "步骤1: 去除PBC跳跃..."
echo "0" | gmx trjconv -s $TPR_FILE -f $XTC_FILE -o temp1.xtc -pbc nojump
if [ $? -ne 0 ]; then
    echo "错误: trjconv -pbc nojump 失败"
    exit 1
fi

# 第二步: 将分子完整化
# 确保分子不被盒子边界分割
echo "步骤2: 分子完整化..."
echo "0" | gmx trjconv -s $TPR_FILE -f temp1.xtc -o temp2.xtc -pbc mol
if [ $? -ne 0 ]; then
    echo "错误: trjconv -pbc mol 失败"
    exit 1
fi

# 第三步: 居中处理 (可选)
# 将感兴趣的分组放在盒子中心
echo "步骤3: 居中处理 (自动选择蛋白质)..."
echo "自动选择蛋白质分组 (分组1) 进行居中处理"

echo -e "1\n0" | gmx trjconv -s $TPR_FILE -f temp2.xtc -o temp3.xtc -center
if [ $? -ne 0 ]; then
    echo "错误: trjconv -center 失败"
    exit 1
fi

# 第四步: 装配处理
# 将分子装配到主盒子中
echo "步骤4: 分子装配..."
echo "0" | gmx trjconv -s $TPR_FILE -f temp3.xtc -o $OUTPUT_XTC -pbc atom
if [ $? -ne 0 ]; then
    echo "错误: trjconv -pbc atom 失败"
    exit 1
fi

# 清理临时文件
echo "清理临时文件..."
rm -f temp1.xtc temp2.xtc temp3.xtc

# 可选: 生成处理后的结构文件
echo "生成最终结构文件..."
echo "0" | gmx trjconv -s $TPR_FILE -f $OUTPUT_XTC -o processed_final.gro -dump 0
if [ $? -ne 0 ]; then
    echo "警告: 生成最终结构文件失败，但主要处理已完成"
fi

echo "PBC处理完成!"
echo "输出文件: $OUTPUT_XTC"
echo "最终结构: processed_final.gro"

# 可选: 生成处理报告
echo "生成处理报告..."
cat << EOF > pbc_processing_report.txt
GROMACS PBC处理报告
==================
处理时间: $(date)
输入文件: $XTC_FILE
输出文件: $OUTPUT_XTC
TPR文件: $TPR_FILE

处理步骤:
1. 去除PBC跳跃 (nojump)
2. 分子完整化 (mol)
3. 居中处理 (center, 分组: 1-蛋白质)
4. 原子装配 (atom)

处理状态: 成功完成
EOF

echo "报告已保存至: pbc_processing_report.txt"
echo "完成!"
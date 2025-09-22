# 质量感知的SLURM批处理方案

## 概述

基于你的需求，我们实现了一个**前台质量检查 + 后台批处理**的解决方案：

1. **前台质量检查**：快速检查所有MD的质量（链数、完整性）
2. **智能过滤**：只为通过质量检查的MD生成SLURM脚本
3. **优化批处理**：避免在问题数据上浪费计算资源

## 质量检查标准

### MD完整性检查
- ✅ 必要文件存在：`md.gro`, `md.xtc/trr`, `md.log`
- ✅ 轨迹文件大小合理（> 1MB）
- ✅ 模拟时间达到预期（> 5000ps）
- ✅ 日志显示正常完成

### PDB结构验证
- ✅ 蛋白质链数正确（默认5条链）
- ✅ 链信息完整有效
- ✅ 坐标范围合理
- ✅ 结构完整性良好

## 使用方法

### 方法1：一键工作流（推荐）

```bash
# 完整工作流：质量检查 + SLURM脚本生成
python md_workflow.py /path/to/md/simulations

# 自定义参数
python md_workflow.py /path/to/md/simulations \
    --expected-chains 3 \
    --batch-size 15 \
    --partition gpu \
    --time 24:00:00
```

### 方法2：分步执行

```bash
# 步骤1：前台质量检查
python md_quality_check.py /path/to/md/simulations

# 步骤2：为合格MD生成SLURM脚本
python generate_slurm.py /path/to/md/simulations \
    --qualified-list ./quality_check_results/qualified_mds.txt
```

### 方法3：跳过质量检查（不推荐）

```bash
# 处理所有MD（忽略质量问题）
python generate_slurm.py /path/to/md/simulations --skip-quality-check
```

## 工作流程详解

### 第一阶段：前台质量检查 (2-10分钟)

```
🔍 搜索MD目录...
📁 发现 150 个潜在的MD目录

📋 开始检查 150 个MD目录...

[  1/150] 1abc_md
   ├─ MD完整性检查...
   ├─ PDB结构验证...
   ✅ 合格

[  2/150] 2xyz_md  
   ├─ MD完整性检查...
   ├─ PDB结构验证...
   ❌ 不合格: 链数异常: 发现3条链，期望5条链

...

🎯 MD质量检查结果
===========================
📊 总计检查: 150 个MD目录
✅ 质量合格: 118 个 (78.7%)
❌ 质量不合格: 32 个

🔍 主要失败原因:
   • 链数异常: 15 个
   • MD未完成: 12 个
   • 缺失文件: 5 个
```

### 第二阶段：SLURM脚本生成 (1-2分钟)

```
📋 从质量检查结果加载了 118 个合格MD

🚀 AfterMD SLURM脚本生成器
===============================
📂 输入目录: /path/to/md/simulations
🔍 质量过滤: 启用
📋 质量列表: ./quality_check_results/qualified_mds.txt
✅ 合格MD数: 118 个
📦 批次大小: 10 任务/作业

✅ SLURM脚本生成成功!
===============================
📊 发现MD任务: 118 个
📦 生成批次: 12 个
📄 生成脚本: 12 个
🚀 提交脚本: ./slurm_scripts/submit_all_batches.sh
```

### 第三阶段：SLURM批处理 (数小时到数天)

```bash
# 提交所有批处理作业
bash ./slurm_scripts/submit_all_batches.sh

# 监控作业状态
squeue -u $USER

# 作业状态示例
JOBID PARTITION     NAME     USER ST       TIME  NODES
12345      gpu aftermd_1     user  R       2:30      1
12346      gpu aftermd_2     user  R       2:25      1
12347      gpu aftermd_3     user PD       0:00      1
```

## 输出文件结构

```
md_workflow_results/
├── quality_check_results/          # 质量检查结果
│   ├── quality_summary.txt        # 质量摘要报告
│   ├── qualified_mds.txt          # 合格MD列表 (118个)
│   ├── failed_mds.txt             # 失败MD列表 (32个)
│   └── quality_check_report.json  # 详细JSON报告
└── slurm_scripts/                  # SLURM脚本
    ├── aftermd_batch_001.sh       # 批处理脚本1 (MD 1-10)
    ├── aftermd_batch_002.sh       # 批处理脚本2 (MD 11-20)
    ├── ...
    ├── aftermd_batch_012.sh       # 批处理脚本12 (MD 111-118)
    └── submit_all_batches.sh      # 批量提交脚本
```

## 核心优势

### 1. 资源优化
- ❌ **传统方式**：150个MD → 150个SLURM作业 → 32个失败浪费资源
- ✅ **新方式**：150个MD → 质量筛选 → 118个SLURM作业 → 0个已知失败

### 2. 时间效率
- 前台质量检查：5-10分钟快速识别问题
- 避免后续数小时的计算资源浪费
- 问题可以立即发现和解决

### 3. 批处理优化
```python
# 质量检查：轻量级，适合前台
- CPU要求低
- 内存占用小  
- 快速并行处理

# 轨迹处理：重量级，适合SLURM
- 需要GPU加速
- 大内存需求
- 长时间运行
```

### 4. 可追溯性
- 详细的质量报告
- 失败原因分析
- 处理过程日志
- 结果统计汇总

## 质量标准自定义

```bash
# 自定义链数要求
python md_workflow.py /data/simulations --expected-chains 3

# 自定义完整性要求
python md_workflow.py /data/simulations \
    --min-traj-size 5.0 \
    --min-sim-time 10000

# 自定义SLURM参数
python md_workflow.py /data/simulations \
    --batch-size 20 \
    --partition gpu \
    --time 48:00:00 \
    --cpus 16
```

## 监控和调试

### 质量问题诊断

```bash
# 查看质量摘要
cat ./quality_check_results/quality_summary.txt

# 查看失败原因
cat ./quality_check_results/failed_mds.txt

# 查看详细报告
cat ./quality_check_results/quality_check_report.json
```

### SLURM作业监控

```bash
# 查看作业状态
squeue -u $USER

# 查看作业详情
scontrol show job 12345

# 查看作业日志
tail -f slurm_scripts/aftermd_batch_001.out
```

## 最佳实践

1. **总是使用质量检查**：`md_workflow.py` 而不是直接 `generate_slurm.py`

2. **合理设置批次大小**：
   - GPU任务：5-15个MD/批次
   - CPU任务：10-30个MD/批次
   
3. **定期检查质量报告**：
   ```bash
   # 查看合格率
   grep "合格率" quality_check_results/quality_summary.txt
   ```

4. **分析失败模式**：
   ```bash
   # 统计失败原因
   cut -f2 quality_check_results/failed_mds.txt | sort | uniq -c
   ```

5. **监控资源使用**：
   ```bash
   # 查看作业效率
   sacct -j 12345 --format=JobID,JobName,Elapsed,TotalCPU,MaxRSS
   ```

## 故障排除

### 常见问题

1. **没有发现合格MD**
   ```bash
   # 检查质量标准是否过于严格
   python md_quality_check.py /data/simulations --expected-chains 3
   ```

2. **质量检查过慢**
   ```bash
   # 质量检查本身很快，如果慢可能是网络存储问题
   # 可以在计算节点上运行
   ```

3. **SLURM脚本生成失败**
   ```bash
   # 检查AfterMD模块是否可用
   python -c "import aftermd; print('AfterMD可用')"
   ```

## 总结

这个质量感知的SLURM方案实现了：

✅ **智能筛选**：只处理质量合格的MD
✅ **资源优化**：避免在问题数据上浪费计算时间
✅ **实时反馈**：前台质量检查提供即时结果
✅ **批量高效**：合格MD的自动化批处理
✅ **可追溯性**：完整的质量报告和处理日志

相比传统的"盲目处理所有MD"方式，这个方案可以节省20-50%的计算资源，并提前发现数据质量问题。
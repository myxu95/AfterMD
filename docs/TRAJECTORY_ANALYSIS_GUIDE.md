# MD轨迹分析流水线使用指南

本指南介绍了如何使用AfterMD的轨迹分析流水线进行RMSD、回转半径、距离等后续分析。

## 分析流程概述

```
质量检测 → PBC预处理 → 轨迹分析 → 结果可视化
    ↓           ↓           ↓           ↓
合格MD筛选  → 标准轨迹   → 分析数据   → 图表报告
```

## 输入文件要求

轨迹分析需要以下标准化输入文件：

```
{MD_Product}_processed/
├── {MD_Product}_processed.xtc    # PBC处理后的轨迹
├── md.gro                        # 参考结构文件 (首选)
├── (或其他.gro/.tpr/.pdb文件)     # 备选参考结构
└── analysis/                     # 分析结果输出目录 (自动创建)
    ├── rmsd/
    ├── radius_gyration/
    ├── distances/
    ├── plots/
    ├── data/
    └── reports/
```

## 使用方法

### 1. 完整工作流程 (推荐)

从质量检测到轨迹分析的完整流程：

```bash
# 方法1: 分步执行 (便于调试)
python md_workflow.py /path/to/md/simulations --quality-only
python md_workflow.py /path/to/md/simulations --slurm-only  
python md_workflow.py /path/to/processed/results --analysis-only

# 方法2: 一次性执行 (质量检测 + SLURM脚本生成)
python md_workflow.py /path/to/md/simulations
```

### 2. 仅轨迹分析

如果已有PBC处理后的轨迹文件：

```bash
# 基本分析 (RMSD + 回转半径 + 距离)
python run_trajectory_analysis.py /path/to/processed/simulations

# 指定分析类型
python run_trajectory_analysis.py /path/to/processed/simulations \
    --analysis rmsd rg distances

# 并行分析 (4个进程)
python run_trajectory_analysis.py /path/to/processed/simulations \
    --max-workers 4

# 分析单个轨迹
python run_trajectory_analysis.py /path/to/processed/simulations \
    --single md_product_1_processed
```

### 3. 预览模式

查看将要分析的轨迹：

```bash
python run_trajectory_analysis.py /path/to/processed/simulations --dry-run
```

## 分析类型

### 支持的分析类型

| 分析类型 | 关键字 | 描述 | 默认选择 |
|---------|--------|------|----------|
| RMSD | `rmsd` | 均方根偏差 | ✅ 默认 |
| 回转半径 | `rg`, `radius_gyration` | 蛋白质紧实度 | ✅ 默认 |
| 距离分析 | `distances`, `distance` | 原子/基团距离 | ✅ 默认 |
| 径向分布函数 | `rdf` | RDF分析 | ❌ 可选 |
| 氢键分析 | `hbonds`, `hydrogen_bonds` | 氢键统计 | ❌ 可选 |

### RMSD分析

自动分析多种原子选择：
- **α碳原子**: `protein and name CA`
- **疏水残基**: `protein and name CA and resname ALA VAL LEU ILE PHE TRP MET`
- **极性残基**: `protein and name CA and resname SER THR ASN GLN`
- **蛋白质骨架**: `protein and backbone`

### 回转半径分析

分析蛋白质的紧实度：
- **整个蛋白质**: `protein`
- **α碳原子**: `protein and name CA`

### 距离分析

默认分析：
- **蛋白质质心距离**: 蛋白质几何中心的自身距离
- **蛋白质末端距离**: N端到C端的距离

## 输出结果

### 目录结构

```
analysis/
├── rmsd/                    # RMSD原始数据
├── radius_gyration/         # 回转半径原始数据
├── distances/               # 距离分析原始数据
├── plots/                   # 所有图表
│   ├── rmsd_*.png          # RMSD时间序列图
│   ├── radius_gyration_*.png # 回转半径时间序列图
│   └── distance_*.png      # 距离时间序列图
├── data/                    # XVG格式数据文件
│   ├── rmsd_*.xvg
│   ├── radius_gyration_*.xvg
│   └── distance_*.xvg
└── reports/                 # 分析报告
    ├── analysis_report.md  # Markdown格式报告
    └── analysis_results.json # JSON格式详细结果
```

### 生成的图表

- **RMSD时间序列图**: 显示结构稳定性变化
- **回转半径图**: 显示蛋白质紧实度变化  
- **距离变化图**: 显示特定距离随时间的变化
- **统计信息**: 平均值、标准差、最大/最小值

### 分析报告

每个轨迹生成详细的Markdown报告，包含：
- 输入文件信息
- 分析摘要和耗时
- 各分析结果的统计数据
- 生成文件列表

## 高级用法

### 自定义分析参数

通过修改分析脚本可以自定义：

```python
# 自定义RMSD选择
rmsd_selections = [
    "protein and name CA",
    "protein and resname ALA",  # 只分析丙氨酸
    "resid 1-50 and name CA"    # 只分析前50个残基
]

# 自定义距离对
distance_pairs = [
    {
        "name": "active_site_distance",
        "sel1": "resid 25 and name CA", 
        "sel2": "resid 100 and name CA",
        "type": "minimum_distance"
    }
]
```

### 并行处理优化

```bash
# 对于大量轨迹，建议使用并行处理
python run_trajectory_analysis.py /path/to/processed/simulations \
    --max-workers 8 \
    --log-level INFO

# 对于长轨迹，可以只运行快速分析
python run_trajectory_analysis.py /path/to/processed/simulations \
    --analysis rmsd rg \
    --max-workers 4
```

## 故障排除

### 常见问题

1. **轨迹文件未找到**
   ```
   ❌ 未找到任何已处理的轨迹文件
   💡 确保目录结构为: {MD_Product}_processed/{MD_Product}_processed.xtc
   ```
   - 检查目录结构是否正确
   - 确认PBC预处理已完成

2. **参考结构文件缺失**
   ```
   ⚠️  跳过 md_product_1_processed: 未找到参考结构文件
   ```
   - 确保存在 `md.gro` 或其他结构文件
   - 检查文件权限

3. **MDAnalysis加载失败**
   ```
   ❌ Failed to load trajectory for RMSD: ...
   ```
   - 检查轨迹文件完整性
   - 确认文件格式正确 (.xtc/.trr)

4. **内存不足**
   ```
   ❌ Memory error during analysis
   ```
   - 减少并行进程数 (`--max-workers 1`)
   - 考虑分批处理轨迹

### 调试模式

```bash
# 开启详细日志
python run_trajectory_analysis.py /path/to/processed/simulations \
    --log-level DEBUG \
    --verbose

# 预览模式检查输入
python run_trajectory_analysis.py /path/to/processed/simulations --dry-run
```

## 性能参考

| 轨迹长度 | 文件大小 | 分析时间 | 推荐配置 |
|---------|----------|----------|----------|
| 1ns | ~10MB | 1-2分钟 | 单进程 |
| 10ns | ~100MB | 5-10分钟 | 2-4进程 |
| 100ns | ~1GB | 20-40分钟 | 4-8进程 |
| 1μs | ~10GB | 2-4小时 | 并行+批处理 |

## 下一步

轨迹分析完成后，可以：

1. **查看分析报告**: `*/analysis/reports/analysis_report.md`
2. **检查生成图表**: `*/analysis/plots/`
3. **进一步定制分析**: 修改分析脚本添加特定分析
4. **批量统计**: 汇总多个轨迹的分析结果
5. **结果发布**: 将图表和数据用于论文和报告

## 技术支持

如有问题或建议，请：
- 查看日志文件: `trajectory_analysis.log`
- 检查错误信息和建议的解决方案
- 提交问题到项目仓库
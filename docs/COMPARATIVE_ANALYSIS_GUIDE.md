# MD Products对比分析指南

本指南介绍如何使用AfterMD的对比分析功能，对多个MD Product的分析结果进行统计对比和可视化。

## 对比分析概述

```
单个MD分析 → 结果收集 → 统计对比 → 可视化呈现
     ↓          ↓         ↓         ↓
  RMSD/Rg等   → 数据汇总  → 差异分析 → 图表报告
```

## 前置条件

对比分析需要以下输入：

1. **已完成的轨迹分析**: 每个MD Product都已运行过轨迹分析
2. **标准文件结构**: 分析结果按标准格式存储
3. **分析结果文件**: 每个MD Product包含 `analysis/reports/analysis_results.json`

```
processed_simulations/
├── md_product_1_processed/
│   └── analysis/
│       └── reports/
│           └── analysis_results.json  ✅ 必需
├── md_product_2_processed/
│   └── analysis/
│       └── reports/
│           └── analysis_results.json  ✅ 必需
└── md_product_N_processed/
    └── analysis/
        └── reports/
            └── analysis_results.json  ✅ 必需
```

## 使用方法

### 1. 基本对比分析

```bash
# 分析所有MD Products的结果
python run_comparative_analysis.py /path/to/processed/simulations

# 指定输出目录
python run_comparative_analysis.py /path/to/processed/simulations \
    --output ./comparative_results
```

### 2. 集成工作流程

```bash
# 完整流程：质量检测 → SLURM处理 → 轨迹分析 → 对比分析
python md_workflow.py /path/to/md/simulations --quality-only
python md_workflow.py /path/to/md/simulations --slurm-only
# ... 运行SLURM作业 ...
python md_workflow.py /path/to/processed/results --analysis-only
python md_workflow.py /path/to/processed/results --comparative-only
```

### 3. 预览模式

查看将要分析的MD Products：

```bash
python run_comparative_analysis.py /path/to/processed/simulations --dry-run
```

### 4. 详细日志

```bash
python run_comparative_analysis.py /path/to/processed/simulations \
    --log-level DEBUG --verbose
```

## 生成的对比图表

### 1. RMSD对比图表

**按原子选择类型分组的对比**:
- `rmsd_comparison_protein_and_name_CA.png` - α碳原子RMSD对比
- `rmsd_comparison_protein_and_backbone.png` - 骨架原子RMSD对比
- `rmsd_comparison_selection_X.png` - 其他选择类型对比

**RMSD热力图**:
- `rmsd_heatmap_all.png` - 所有MD Products和选择类型的热力图

**图表特点**:
- 条形图显示平均值和标准差
- 不同颜色区分各MD Product
- 网格线便于数值读取

### 2. 回转半径对比图表

- `rg_comparison_protein.png` - 整个蛋白质回转半径对比
- `rg_comparison_protein_and_name_CA.png` - α碳原子回转半径对比

### 3. 距离分析对比图表

- `distance_comparison_protein_center_of_mass.png` - 蛋白质质心距离对比
- `distance_comparison_protein_ends.png` - 蛋白质末端距离对比

### 4. 质量统计图表

**轨迹质量对比**:
- `trajectory_size_comparison.png` - 轨迹文件大小对比
- `analysis_time_comparison.png` - 分析耗时对比
- `analysis_success_rate.png` - 分析成功率对比

### 5. 综合仪表板

- `comparative_dashboard.png` - 四象限综合视图
  - RMSD汇总对比
  - 回转半径汇总对比
  - 文件大小 vs 分析时间散点图
  - 分析类型分布饼图

## 生成的数据文件

### CSV格式数据

存储在 `data/` 目录：

1. **`rmsd_comparison_data.csv`**:
   ```csv
   MD_Product,Selection,Mean_RMSD,Std_RMSD,Max_RMSD,Min_RMSD,Selection_Description
   md_product_1,protein_and_name_CA,1.25,0.15,2.10,0.95,"protein and name CA"
   md_product_2,protein_and_name_CA,1.45,0.22,2.35,1.05,"protein and name CA"
   ```

2. **`radius_gyration_comparison_data.csv`**:
   ```csv
   MD_Product,Selection,Mean_Rg,Std_Rg,Max_Rg,Min_Rg,Selection_Description
   md_product_1,protein,18.5,0.8,20.2,17.1,"protein"
   md_product_2,protein,17.9,0.6,19.8,16.8,"protein"
   ```

3. **`distance_comparison_data.csv`**:
   ```csv
   MD_Product,Distance_Type,Mean_Distance,Std_Distance,Max_Distance,Min_Distance
   md_product_1,protein_ends,45.2,3.1,52.8,38.9
   md_product_2,protein_ends,43.8,2.8,50.1,37.5
   ```

## 分析报告

### Markdown报告

`reports/comparative_analysis_report.md` 包含：

1. **基本统计**:
   - 分析的MD Products总数
   - 总轨迹大小和分析时间
   - 平均分析时间

2. **RMSD分析摘要**:
   - 按选择类型的统计数据
   - RMSD范围和平均值

3. **回转半径分析摘要**:
   - 各选择类型的Rg统计
   - 紧实度变化范围

4. **距离分析摘要**:
   - 各距离类型的统计
   - 距离变化范围

5. **生成文件列表**:
   - 所有图表文件
   - 所有数据文件

## 应用场景

### 1. 蛋白质稳定性比较

比较不同条件下蛋白质的结构稳定性：

```python
# 示例：比较不同pH条件下的MD结果
# md_product_pH7_processed/
# md_product_pH8_processed/
# md_product_pH9_processed/

# RMSD对比显示pH对结构稳定性的影响
# 回转半径对比显示紧实度变化
```

### 2. 突变体效应分析

比较野生型与突变体蛋白质：

```python
# 示例：比较突变体效应
# md_product_wildtype_processed/
# md_product_mutant_A_processed/
# md_product_mutant_B_processed/

# 距离分析显示突变对蛋白质构象的影响
# RMSD分析显示突变对稳定性的影响
```

### 3. 药物结合效应研究

比较有无配体结合的蛋白质：

```python
# 示例：配体结合效应
# md_product_apo_processed/        # 无配体
# md_product_ligand1_processed/    # 配体1结合
# md_product_ligand2_processed/    # 配体2结合

# 回转半径分析显示配体对蛋白质紧实度的影响
# 距离分析显示活性位点变化
```

### 4. 时间尺度效应研究

比较不同模拟时长的结果：

```python
# 示例：时间尺度效应
# md_product_10ns_processed/
# md_product_100ns_processed/
# md_product_1us_processed/

# RMSD分析显示平衡时间
# 质量统计显示计算资源消耗
```

## 高级分析

### 1. 自定义数据提取

```python
from aftermd.analysis.comparative_analysis import ComparativeAnalyzer

analyzer = ComparativeAnalyzer("./custom_output")
analyzer.collect_analysis_results("/path/to/processed")

# 提取特定数据
rmsd_df = analyzer.extract_rmsd_data()
filtered_data = rmsd_df[rmsd_df['Selection'].str.contains('CA')]

# 自定义图表
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.bar(filtered_data['MD_Product'], filtered_data['Mean_RMSD'])
plt.title('Custom RMSD Comparison')
plt.show()
```

### 2. 统计显著性检验

```python
import scipy.stats as stats

# 提取两组数据进行t检验
group1_rmsd = rmsd_df[rmsd_df['MD_Product'].str.contains('condition1')]['Mean_RMSD']
group2_rmsd = rmsd_df[rmsd_df['MD_Product'].str.contains('condition2')]['Mean_RMSD']

t_stat, p_value = stats.ttest_ind(group1_rmsd, group2_rmsd)
print(f"T-test: t={t_stat:.3f}, p={p_value:.3f}")
```

### 3. 相关性分析

```python
# 分析RMSD与回转半径的相关性
rmsd_df = analyzer.extract_rmsd_data()
rg_df = analyzer.extract_radius_gyration_data()

# 合并数据
merged_df = rmsd_df.merge(rg_df, on='MD_Product', suffixes=('_rmsd', '_rg'))

# 计算相关性
correlation = merged_df['Mean_RMSD'].corr(merged_df['Mean_Rg'])
print(f"RMSD vs Rg correlation: {correlation:.3f}")
```

## 性能考虑

### 数据量估算

| MD Products数量 | 生成图表数量 | 处理时间 | 内存使用 |
|----------------|------------|----------|----------|
| 5-10个 | 10-15个 | 10-30秒 | < 100MB |
| 20-50个 | 15-25个 | 30-60秒 | 100-500MB |
| 100+个 | 20-30个 | 1-3分钟 | 500MB-1GB |

### 优化建议

1. **大数据集处理**:
   - 分批处理MD Products
   - 使用子集进行初步分析

2. **图表优化**:
   - 大量MD Products时使用热力图
   - 限制显示的MD Product数量

3. **内存管理**:
   - 及时释放不需要的数据
   - 使用流式处理大文件

## 故障排除

### 常见问题

1. **未找到分析结果**:
   ```
   ❌ 未找到任何分析结果文件
   💡 确保已运行轨迹分析: python run_trajectory_analysis.py
   ```

2. **数据不一致**:
   ```
   ⚠️  部分MD Product缺少某些分析类型
   ```
   - 检查轨迹分析是否完整
   - 确认分析参数一致

3. **图表生成失败**:
   ```
   ❌ 图表生成过程出错
   ```
   - 检查matplotlib和seaborn版本
   - 确认有足够的内存

## 下一步

对比分析完成后，可以：

1. **深入分析**: 基于对比结果进行假设验证
2. **发表展示**: 使用生成的图表制作presentation
3. **进一步计算**: 基于CSV数据进行统计检验
4. **方法改进**: 根据对比结果优化MD参数
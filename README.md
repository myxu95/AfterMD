# AfterMD - GROMACS MD Analysis Toolkit

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![GROMACS](https://img.shields.io/badge/GROMACS-2023+-green.svg)](https://gromacs.org)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**AfterMD** 是一个强大的 GROMACS 分子动力学 (MD) 模拟分析工具包，专为大规模批量处理和 HPC 集群部署而设计。

## 🚀 核心特性

- **🔄 智能 PBC 处理**: 自动最短链检测 + 三步 PBC 校正流程
- **📦 批量处理**: 一键处理数百个 MD 轨迹，支持并行执行
- **🖥️ SLURM 集群支持**: 自动生成 SLURM 脚本，无缝集群部署
- **🧠 智能文件发现**: 自动识别和验证 MD 文件完整性
- **📊 丰富分析工具**: RMSD、RDF、距离、氢键等分析模块

## ⚡ 快速开始

### 安装

#### Conda 环境安装 (推荐)
```bash
git clone https://github.com/your-username/AfterMD.git
cd AfterMD
conda env create -f environment.yml
conda activate aftermd
pip install -e .
```

#### Pip 直接安装
```bash
git clone https://github.com/your-username/AfterMD.git
cd AfterMD
pip install -r requirements.txt
pip install -e .
```

### 基本使用

#### 1. 批量 PBC 处理
```python
from aftermd import process_md_tasks

# 一行代码处理所有 MD 轨迹
results = process_md_tasks("/path/to/md_simulations")
print(f"处理了 {results['successful']}/{results['total_tasks']} 个任务")
```

#### 2. 生成 SLURM 集群脚本
```python
from aftermd import generate_slurm_scripts_for_md_tasks

# 为 20 个任务生成 SLURM 脚本，每批 10 个
results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/md_simulations",
    tasks_per_batch=10
)
# 结果：生成 2 个 SLURM 脚本 + 批量提交脚本
```

#### 3. 命令行使用
```bash
# 批量处理
python -m aftermd.batch_process /data/simulations

# 生成 SLURM 脚本
python -m aftermd.utils.slurm_generator /data/simulations --tasks-per-batch 10

# 提交所有作业
bash slurm_scripts/submit_all_batches.sh
```

## 📁 支持的文件结构

AfterMD 自动发现以下文件结构：

```
your_simulations/
├── task1/
│   ├── md.xtc          # 直接在任务目录
│   └── md.tpr
├── task2/
│   └── prod/
│       ├── md.xtc      # 在 prod 子文件夹
│       └── md.tpr
└── task3/
    ├── md.xtc
    └── md.tpr
```

## 🔬 技术亮点

### 智能最短链检测
```python
# 自动检测最短肽链用于最优 PBC 中心化
# 1. 使用 gmx make_ndx splitch 分析所有链
# 2. 过滤离子和小分子
# 3. 选择最短有效肽链
# 4. 生成精确的 index 文件
```

### 三步 PBC 处理流程
```bash
# Step 1: 使用最短链进行中心化
gmx trjconv -center -pbc atom -n shortest_chain.ndx

# Step 2: 保持分子完整性
gmx trjconv -pbc whole

# Step 3: 使用 backbone 进行结构对齐
gmx trjconv -fit rot+trans
```

### 自动结构文件复制
```python
# AfterMD 自动复制 md.gro 文件到输出目录
# 优先级: md.gro > prod.gro > production.gro > {trajectory_name}.gro
# 便于后续轨迹可视化和分析
```

### SLURM 作业优化
```bash
# 生成的作业名称格式: amd_dataset_XofY
amd_antibody_sim_1of4    # 第1批，共4批
amd_antibody_sim_2of4    # 第2批，共4批
amd_antibody_sim_3of4    # 第3批，共4批
amd_antibody_sim_4of4    # 第4批，共4批

# 特点：
# ✓ 长度 ≤24 字符，squeue 显示友好
# ✓ 清晰的进度指示
# ✓ 数据集易于识别
```

## 📊 高级分析

### RMSD 分析
```python
from aftermd.analysis import RMSDCalculator

rmsd_calc = RMSDCalculator(trajectory, topology)
rmsd_data = rmsd_calc.calculate_backbone_rmsd()
```

### 径向分布函数
```python
from aftermd.analysis import RDFCalculator

rdf_calc = RDFCalculator(trajectory, topology)
rdf_data = rdf_calc.calculate_protein_water_rdf()
```

### 自定义分析流程
```python
from aftermd.preprocessing import PBCProcessor
from aftermd.analysis import RMSDCalculator, RDFCalculator

# 1. PBC 处理
pbc = PBCProcessor()
processed_traj = pbc.remove_pbc(trajectory, topology, output)

# 2. 多种分析
rmsd = RMSDCalculator(processed_traj, topology)
rdf = RDFCalculator(processed_traj, topology)

# 3. 生成结果
rmsd_data = rmsd.calculate_backbone_rmsd()
rdf_data = rdf.calculate_protein_water_rdf()
```

## 🖥️ 集群使用示例

### 基本集群处理
```bash
# 1. 生成 SLURM 脚本
python -m aftermd.utils.slurm_generator /data/simulations \
  --tasks-per-batch 10 \
  --partition gpu \
  --time 24:00:00

# 2. 提交作业
bash slurm_scripts/submit_all_batches.sh

# 3. 监控作业
squeue -u $USER
```

### 自定义集群配置
```python
slurm_params = {
    "partition": "gpu",
    "time": "24:00:00",
    "cpus_per_task": 16,
    "gres": "gpu:2",
    "memory": "64G",
    "conda_env": "aftermd"
}

results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/simulations",
    tasks_per_batch=5,
    slurm_params=slurm_params
)
```

## 📚 文档

- [批量处理指南](docs/batch_processing_guide.md)
- [SLURM 集群部署](docs/slurm_cluster_guide.md)
- [项目概览](PROJECT_OVERVIEW.md)
- [功能总结](FEATURES_SUMMARY.md)

## 🎯 使用场景

- **🧬 蛋白质动力学研究**: 蛋白质折叠、构象变化分析
- **💊 药物设计**: 药物-蛋白结合、药效评估
- **🧪 膜蛋白研究**: 膜蛋白稳定性、传输机制
- **⚗️ 酶学研究**: 酶催化机制、活性位点动态
- **🔬 分子相互作用**: 蛋白质-蛋白质、蛋白质-核酸相互作用

## 🌟 性能优势

| 特性 | 传统方法 | AfterMD |
|------|---------|---------|
| 批量处理 | 手动逐个处理 | 全自动批量处理 |
| 处理时间 | 数天到数周 | 数小时到数天 |
| 错误处理 | 停止整个流程 | 继续处理其他任务 |
| 集群部署 | 手写脚本 | 自动生成 |
| 文件管理 | 手动组织 | 智能发现 |

## 🛠️ 系统要求

- **Python**: 3.8+
- **GROMACS**: 2020+
- **依赖包**: numpy, pandas, matplotlib, pathlib
- **可选**: plotly (可视化), MDAnalysis (高级分析)

## 🔧 开发者指南

### 添加新的分析模块
```python
from aftermd.analysis import BaseAnalyzer

class MyAnalyzer(BaseAnalyzer):
    def __init__(self, trajectory, topology):
        super().__init__(trajectory, topology)
    
    def calculate_my_property(self):
        # 实现你的分析逻辑
        return results
```

### 扩展批量处理功能
```python
from aftermd.utils import BatchProcessor

def my_custom_processor(trajectory, topology, output_dir):
    # 实现自定义处理逻辑
    return results

# 集成到批量处理
batch = BatchProcessor()
results = batch.process_files(file_list, my_custom_processor)
```

## 🤝 贡献

欢迎贡献代码、报告 bug 或提出功能建议！

1. Fork 项目
2. 创建功能分支 (`git checkout -b feature/amazing-feature`)
3. 提交更改 (`git commit -m 'Add amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 开启 Pull Request

## 📄 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

## 🙏 致谢

- [GROMACS](https://gromacs.org) - 强大的 MD 模拟软件
- [MDAnalysis](https://mdanalysis.org) - 轨迹分析库
- 所有使用和改进 AfterMD 的研究者们

## 📞 联系

- 项目主页: [https://github.com/your-username/AfterMD](https://github.com/your-username/AfterMD)
- 问题报告: [Issues](https://github.com/your-username/AfterMD/issues)
- 文档: [Wiki](https://github.com/your-username/AfterMD/wiki)

---

**AfterMD - 让 MD 分析变得简单高效！** 🚀
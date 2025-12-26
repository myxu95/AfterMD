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

**简单用法 - 一行代码处理所有轨迹**
```python
from aftermd import process_md_tasks

# 最简单的使用方式
results = process_md_tasks("/path/to/md_simulations")
print(f"处理了 {results['successful']}/{results['total_tasks']} 个任务")
```

**完整参数使用**
```python
# 带完整参数的批量处理
results = process_md_tasks(
    simulations_path="/data/md_simulations",    # 输入目录
    output_dir="/data/processed_results",       # 输出目录
    dt=10.0,                                    # 轨迹下采样：每10ps取一帧
    max_workers=4,                              # 并行worker数量（同时处理4个任务）
    gmx_executable="gmx"                        # GROMACS可执行文件
)

# 查看详细结果
print(f"总任务数: {results['total_tasks']}")
print(f"成功处理: {results['successful']}")
print(f"失败任务: {results['failed']}")
```

**单个轨迹 PBC 处理**
```python
from aftermd import PBCProcessor

# 处理单个轨迹
pbc_processor = PBCProcessor()

# 基本PBC处理（三步骤：居中->整体化->拟合）
processed_traj = pbc_processor.remove_pbc(
    trajectory="md.xtc",
    topology="md.tpr",
    output="processed_trajectory.xtc",
    dt=10.0  # 可选：轨迹下采样
)

# 综合PBC处理（自动创建输出目录并复制结构文件）
results = pbc_processor.comprehensive_pbc_process(
    trajectory="md.xtc",
    topology="md.tpr",
    output_dir="processed_output"
)
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

**批量PBC处理**
```bash
# 最简单用法
python -m aftermd.batch_process /data/simulations

# 带参数的完整用法
python -m aftermd.batch_process /data/simulations \
  --output /data/processed \
  --dt 10.0 \
  --workers 4 \
  --verbose

# 仅检查任务状态（不处理）
python -m aftermd.batch_process /data/simulations --check-only

# 查看帮助
python -m aftermd.batch_process --help
```

**任务发现和状态检查**
```python
from aftermd import discover_md_tasks, check_task_status

# 发现所有有效的MD任务
tasks = discover_md_tasks("/data/simulations")
for task_name, (traj, topo) in tasks.items():
    print(f"{task_name}: {traj}, {topo}")

# 检查任务状态
status = check_task_status("/data/simulations")
for task, info in status.items():
    print(f"{task}: {info['status']} - {info['reason']}")
```

**SLURM 集群脚本生成**
```bash
# 基础脚本生成
python scripts/generate_slurm.py /data/simulations --batch-size 10

# 或使用模块方式
python -m aftermd.utils.slurm_generator /data/simulations --tasks-per-batch 10

# 带质量检查的完整流程（推荐）
python scripts/md_quality_check.py /data/simulations
python scripts/generate_slurm.py /data/simulations \
  --qualified-list ./quality_check_results/qualified_mds.txt \
  --batch-size 10 \
  --partition gpu \
  --time 24:00:00

# 提交所有作业
bash slurm_scripts/submit_all_batches.sh

# 监控作业状态
squeue -u $USER | grep amd_
```

## 📁 支持的文件结构

AfterMD 智能发现以下目录结构，**无需手动指定文件路径**：

### 基本结构（推荐）
```
your_simulations/
├── task1/
│   ├── md.xtc          # 轨迹文件：直接在任务目录
│   ├── md.tpr          # 拓扑文件：直接在任务目录
│   └── md.gro          # 结构文件：自动复制到输出目录
├── task2/
│   ├── md.xtc
│   ├── md.tpr
│   └── prod.gro        # 备选结构文件
└── task3/
    ├── md.xtc
    └── md.tpr
```

### Production子目录结构
```
your_simulations/
├── task1/
│   └── prod/           # production子目录
│       ├── md.xtc      # 在prod子文件夹中
│       ├── md.tpr
│       └── md.gro
├── task2/
│   └── prod/
│       ├── md.xtc
│       └── md.tpr
```

### 文件发现优先级

**MD输入文件搜索规则（必须成对存在）:**
1. **最高优先级**: 任务目录根目录中的 `md.xtc` + `md.tpr`
2. **次要优先级**: `prod/` 子目录中的 `md.xtc` + `md.tpr`

**重要说明:**
- ✅ 必须同时找到 `md.xtc` 和 `md.tpr` 两个文件
- ✅ 文件名必须完全匹配（区分大小写）
- ❌ 不支持其他命名如 `traj.xtc`, `production.xtc`
- ❌ 不支持 `.trr` 格式轨迹文件
- 🔍 任务名称 = 目录名称

**结构文件发现优先级（可选，用于可视化）:**
3. **结构文件优先级**: `md.gro` > `prod.gro` > `production.gro` > `{trajectory_name}.gro`

### 输出结构
```
processed_results/
├── task1/
│   ├── task1_processed.xtc    # 处理后的轨迹
│   ├── md.gro                 # 复制的结构文件
│   └── processing_log.txt     # 处理日志
├── task2/
│   ├── task2_processed.xtc
│   └── md.gro
```

## 🔬 技术亮点

### 严格最短链检测 (Strict Shortest Chain Detection)
```python
# 简化的严格最短链检测流程 - 基于GROMACS标准工具
# 1. 执行 'gmx make_ndx -f topology.tpr' 加载标准组 (0-17)
# 2. 执行 'splitch 1' 命令按链拆分Protein组
# 3. 解析输出找到新生成的链组 (18+) 中最短的一个
# 4. 选择最短链(peptide)作为center group，严格无fallback策略
```

**严格检测原理:**
- 🎯 **复合物系统必须**: 必须检测到≥2条蛋白质链
- 🔍 **peptide识别**: 最短链长度10-300原子 (理论peptide≤20AAs)，符合peptide特征
- ⚡ **gmx splitch工具**: 使用GROMACS内置splitch命令可靠拆分链
- 🚫 **无备选方案**: 检测失败必须停止处理，保证科学正确性

**检测输出示例:**
```bash
# gmx make_ndx splitch 1 输出：
# 18 Protein_chain_A     :  3456 atoms
# 19 Protein_chain_B     :  3456 atoms
# 20 Protein_chain_C     :  3456 atoms
# 21 Protein_chain_D     :  3456 atoms
# 22 Protein_chain_E     :   456 atoms  ← 选择最短链(peptide)

# 返回结果: group_id="22", index_file="chains.ndx"
```

### 三步 PBC 处理流程
```bash
# Step 1: 严格最短链居中 - 使用检测到的peptide组进行精确中心化
gmx trjconv -f trajectory.xtc -s topology.tpr -o temp_centered.xtc \
  -center -pbc mol -n chains.ndx
# 选择组22 (检测到的最短链peptide)

# Step 2: 分子完整性 - 确保所有分子保持完整
gmx trjconv -f temp_centered.xtc -s topology.tpr -o temp_whole.xtc \
  -pbc whole

# Step 3: 结构对齐 - 消除旋转和平移运动
gmx trjconv -f temp_whole.xtc -s topology.tpr -o processed.xtc \
  -fit rot+trans
```

**为什么采用严格检测？**
- ✅ **科学正确性**: 使用错误的center group会导致分析结果无效
- ✅ **peptide专用**: 复合物系统中peptide是最佳的居中选择
- ✅ **工具可靠性**: gmx splitch是GROMACS官方链拆分工具
- ✅ **处理安全性**: 检测失败立即停止，避免错误处理
- ❌ **无备选方案**: 不允许任何fallback策略，确保处理质量

### 自动结构文件复制
```python
# AfterMD 自动复制 md.gro 文件到输出目录
# 优先级: md.gro > prod.gro > production.gro > {trajectory_name}.gro
# 便于后续轨迹可视化和分析
```

### SLURM 作业优化和最佳实践

**智能作业命名**
```bash
# 生成的作业名称格式: amd_dataset_XofY
amd_Human_ClassI_1of4    # 第1批，共4批
amd_Human_ClassI_2of4    # 第2批，共4批
amd_Human_ClassI_3of4    # 第3批，共4批
amd_Human_ClassI_4of4    # 第4批，共4批

# 特点：
# ✓ 长度 ≤24 字符，squeue 显示友好
# ✓ 清晰的进度指示
# ✓ 数据集易于识别
```

**批次大小选择建议**
```bash
# 根据任务规模选择批次大小：
# 小规模 (<20个任务)    : --batch-size 5
# 中等规模 (20-100个)   : --batch-size 10-15
# 大规模 (100-500个)    : --batch-size 20-30
# 超大规模 (>500个)     : --batch-size 50+

# 考虑因素：
# - 单个任务处理时间（轨迹大小）
# - 集群排队情况
# - 失败重试代价
```

**监控和故障排除**
```bash
# 检查作业详细信息
scontrol show job JOBID

# 查看作业历史
sacct -j JOBID --format=JobID,JobName,State,ExitCode,Start,End

# 检查失败原因
tail -n 50 slurm-JOBID.out

# 重新提交失败的作业
sbatch slurm_scripts/aftermd_batch_X.sh

# 统计完成情况
find /output/dir -name "md_processed.xtc" | wc -l
```

**集群资源优化**
```bash
# CPU密集型任务（PBC处理）
--cpus 8-16    # 适中的CPU分配
--memory 16G   # 根据轨迹大小调整

# GPU加速（如果GROMACS支持）
--gres gpu:1   # 单GPU通常足够

# 时间估算（每个任务）
# 小轨迹 (<1GB): 10-30分钟
# 中等轨迹 (1-5GB): 30-90分钟
# 大轨迹 (>5GB): 1-3小时
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

### 批量PBC处理工作流

**第一步：检查任务状态**
```bash
# 检查有多少个有效的MD任务
python -m aftermd.batch_process /data/simulations --check-only

# 输出示例：
# ✅ task1: Complete MD files found in task directory
# ✅ task2: Complete MD files found in prod subfolder
# ❌ task3: Missing files: md.tpr
# Task Status Summary: 2/3 valid tasks
```

**第二步：本地批量处理（小规模）**
```bash
# 处理少量任务（<50个）- 本地并行处理
python -m aftermd.batch_process /data/simulations \
  --output /data/processed \
  --dt 10.0 \
  --workers 4 \
  --verbose

# 实时查看处理进度
tail -f processing.log
```

### 大规模SLURM集群处理（推荐）

**为什么使用SLURM？**
- ❌ **避免占用登录节点资源** - 直接运行会占用头节点CPU
- ✅ **专用计算节点分配** - 队列系统分配计算资源
- ✅ **并行处理能力** - 可同时运行多个批次
- ✅ **作业持久性** - 登出后继续运行
- ✅ **集群最佳实践** - 符合HPC环境规范

**第三步：生成SLURM脚本**
```bash
# 推荐完整流程（带质量检查）
python scripts/md_quality_check.py /data/simulations
python scripts/generate_slurm.py /data/simulations \
  --qualified-list ./quality_check_results/qualified_mds.txt \
  --batch-size 10 \
  --partition gpu \
  --time 24:00:00 \
  --cpus 11 \
  --dt 10.0

# 快速生成（跳过质量检查）
python scripts/generate_slurm.py /data/simulations \
  --skip-quality-check \
  --batch-size 10 \
  --partition gpu \
  --time 12:00:00

# 生成的文件：
# slurm_scripts/aftermd_batch_1.sh    (第1批任务)
# slurm_scripts/aftermd_batch_2.sh    (第2批任务)
# slurm_scripts/submit_all_batches.sh (批量提交脚本)
```

**第四步：提交和监控作业**
```bash
# 一键提交所有批次
bash slurm_scripts/submit_all_batches.sh

# 或手动提交单个批次
sbatch slurm_scripts/aftermd_batch_1.sh
sbatch slurm_scripts/aftermd_batch_2.sh

# 监控作业状态
squeue -u $USER | grep amd_

# 查看实时日志
tail -f slurm-*.out

# 检查完成情况
ls /data/simulations_processed/*/md_processed.xtc | wc -l
```

**SLURM脚本高级配置**
```bash
# 自定义集群资源
python scripts/generate_slurm.py /data/simulations \
  --batch-size 15 \
  --partition gpu \
  --time 48:00:00 \
  --cpus 16 \
  --memory 64G \
  --gpu gpu:2 \
  --dt 5.0

# 使用自定义模板
python scripts/generate_slurm.py /data/simulations \
  --template ./my_slurm_template.sh \
  --batch-size 20

# 查看配置（不生成脚本）
python scripts/generate_slurm.py /data/simulations --dry-run
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

### 用户指南
- [批量处理指南](docs/batch_processing_guide.md) - 本地和集群批量处理
- [SLURM 集群部署](docs/slurm_cluster_guide.md) - HPC集群使用详解
- [轨迹分析指南](docs/TRAJECTORY_ANALYSIS_GUIDE.md) - MD轨迹分析方法
- [质量控制指南](docs/quality_control_guide.md) - MD质量检查流程

### 开发者文档
- [项目概览](PROJECT_OVERVIEW.md) - 架构和设计理念
- [功能总结](FEATURES_SUMMARY.md) - 完整功能列表
- [API参考](docs/api_reference.md) - 编程接口说明

### 快速参考
```bash
# 常用命令速查
python scripts/md_quality_check.py /data/sims          # 质量检查
python scripts/generate_slurm.py /data/sims --help     # SLURM脚本生成
python -m aftermd.batch_process /data/sims --help      # 本地批量处理
```

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
| 并行处理 | 单任务串行 | 多Worker并行 |
| 错误处理 | 停止整个流程 | 继续处理其他任务 |
| 集群部署 | 手写脚本 | 自动生成 |
| 文件管理 | 手动组织 | 智能发现 |

### Worker并行处理说明

**Worker = 同时并行处理的任务数量**

```python
# 示例：10个MD任务，4个workers
max_workers=4

# 处理顺序：
# 第1轮：Worker1处理task1，Worker2处理task2，Worker3处理task3，Worker4处理task4
# 第2轮：Worker1处理task5，Worker2处理task6，Worker3处理task7，Worker4处理task8
# 第3轮：Worker1处理task9，Worker2处理task10
```

**Worker数量选择建议：**
- **CPU密集型**: `max_workers = CPU核心数 - 1`
- **内存限制**: 大轨迹文件时减少worker数量
- **默认推荐**: 4-8个workers适合大多数场景

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

## 📞 联系

- 项目主页: [https://github.com/myxu95/AfterMD](https://github.com/myxu95/AfterMD)
- 邮箱：myuxu@zju.edu.cn

---

**AfterMD - 让 MD 分析变得简单高效！** 🚀

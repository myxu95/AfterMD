# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Status
项目运行在服务器上的，但是修改都在本地进行，我会把我们修改的内容更新在服务器上来测试我们的项目不需要你进行本地测试
This appears to be a new or empty repository under `/Users/xumingyu/Documents/work/research/AfrerMD`. The project name suggests it may be related to "After Markdown" or markdown processing research.

## 语言
1.脚本中避免出现中文和emoji
2.对话采取中文对话

## 文件操作规则
1.及时清理测试脚本，保持项目整洁性
2.markdown只保留对用户必要的，对于某个步骤修改的解释可以放在一个临时文件夹里或者直接以对话的形式告诉我
3.创建文件时的命名要兼顾在整体程序中的功能和职责，也要保持精简

## 文件命名规范
### 处理脚本命名规则
对于PBC处理相关的脚本，按以下格式命名：
- **单任务处理**: `pbc_process.py` - 处理单个MD轨迹的PBC校正
- **批量处理**: `batch_pbc.py` - 本地批量处理多个MD任务
- **集群处理**: `batch_pbc_slurm.py` - 生成SLURM集群批处理脚本

### 命名原则
1. **功能优先**: 体现核心处理步骤(如pbc_process)
2. **规模区分**: single < batch < batch_slurm
3. **执行环境**: 本地处理 vs 集群调度(slurm)
4. **简洁明确**: 避免冗余词汇，保持功能性描述

## log板块
对于文件运行的反馈要抓住重点，只在关键的地方去做反馈


## Getting Started

Since this is an empty repository, you'll need to:

1. Initialize the project structure based on the intended technology stack
2. Set up package management files (package.json, requirements.txt, Cargo.toml, etc.)
3. Configure build and development tools
4. Establish testing framework

## Development Commands

### Installation and Setup
```bash
# Install in development mode
pip install -e .

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage
```bash
# Run example script
python examples/basic_usage.py

# Install and test import
python -c "import aftermd; print('AfterMD imported successfully')"
```

### Testing
```bash
# Run tests (when implemented)
pytest tests/

# Run specific test module
pytest tests/test_trajectory.py
```

### Code Quality
```bash
# Format code
black aftermd/
isort aftermd/

# Lint code
flake8 aftermd/
pylint aftermd/
```

## Architecture

### Project Structure
```
aftermd/
├── __init__.py                 # Main package imports
├── utils/                      # Utility modules (tools)
│   ├── __init__.py
│   ├── batch_processor.py     # BatchProcessor class
│   ├── plotting.py            # PlotManager class  
│   └── path_manager.py        # PathManager class
├── analysis/                   # Analysis modules (core functionality)
│   ├── __init__.py
│   ├── quality/               # Quality control modules
│   │   ├── __init__.py
│   │   ├── md_completeness.py    # MDCompletenessChecker
│   │   ├── structure_validator.py # StructureValidator  
│   │   ├── batch_tracker.py      # BatchTracker
│   │   └── quality_reporter.py   # QualityReporter
│   ├── trajectory/            # Trajectory analysis modules
│   │   ├── __init__.py
│   │   ├── rmsd.py           # RMSDCalculator
│   │   ├── rdf.py            # RDFCalculator  
│   │   ├── radius_gyration.py # RadiusGyrationCalculator
│   │   ├── distance.py       # DistanceCalculator
│   │   └── hydrogen_bonds.py # HydrogenBondAnalyzer
│   └── structure/             # Structure analysis modules
│       ├── __init__.py
│       ├── bfactor.py        # BFactorAnalyzer
│       ├── contact_map.py    # ContactMapCalculator
│       ├── geometry.py       # GeometryAnalyzer
│       └── atom_info.py      # AtomInfoExtractor
└── preprocessing/              # Preprocessing modules
    ├── __init__.py
    └── pbc_processor.py       # PBCProcessor class
```

### Module Organization

**Utils (���具模块)**: 
- `BatchProcessor`: Parallel processing of multiple MD files
- `PathManager`: File organization and directory structure management
- `PlotManager`: Publication-ready plotting from XVG and analysis data

**Analysis/Quality (质量控制模块)**:
- `MDCompletenessChecker`: MD simulation completeness verification
- `StructureValidator`: PDB structure quality validation and chain analysis
- `BatchTracker`: Batch processing tracking and duplicate detection
- `QualityReporter`: Comprehensive quality analysis reporting

**Analysis/Trajectory (轨迹分析模块)**:
- `RMSDCalculator`: RMSD calculation (MDAnalysis + GROMACS)
- `RDFCalculator`: Radial distribution function analysis
- `RadiusGyrationCalculator`: Radius of gyration and components
- `DistanceCalculator`: Distance measurements (COM, minimum, atom-atom)
- `HydrogenBondAnalyzer`: Hydrogen bond analysis

**Analysis/Structure (结构分析模块)**:
- `BFactorAnalyzer`: B-factor extraction and analysis by residue/atom
- `ContactMapCalculator`: Contact map calculation and analysis
- `GeometryAnalyzer`: Geometric properties and shape analysis
- `AtomInfoExtractor`: Comprehensive atom information extraction

**Preprocessing (预处理模块)**:
- `PBCProcessor`: PBC artifact removal and trajectory preparation

### Data Flow
1. Raw MD data → Quality modules (completeness and structure validation)
2. Validated data → PBCProcessor (PBC removal and fitting)
3. Processed trajectories → Individual trajectory analysis modules
4. PDB structures → Individual structure analysis modules
5. Analysis results → PlotManager (visualization)
6. Quality reports → QualityReporter (issue tracking and statistics)
7. BatchProcessor coordinates parallel execution across all modules

## Notes

-脚本中避免出现中文和emoji
-对于功能模块我们要综合应用的考虑来设计怎样的绘图可以呈现出直观的结果

## 主要语言
Python 

## 简介
本项目的目的为对gromacs MD production的结果做分析，主要利用MDAnalysis和gromacs内置的程序来进行处理分析

## 初步框架的构建
项目的功能会逐渐的增加，功能主要分为静态的结构的分析和动态的轨迹的分析
## 设计的模块：
### 批量处理模块
我们每个功能都有可能需要对一批文件进行处理，这里需要自动对批量文件处理的能力

### 路径管理模块
我们中间的路径错综复杂，需要集成一个路径管理的模块来专门负责这个部分的内容，比如在处理文件后，结果存放的位置等等

### 绘图模块
一般结果产生的都是xvg之类的文件，我们需要针对这些文件进行可视化的绘图，我们可以集成一些专业的包来处理绘图的问题

## 核心功能部分：
### MD Production质量分析模块

MD Production质量分析模块是AfterMD的核心质量控制功能，用于自动检测MD模拟的完整性、正确性和数据质量。

#### 主要功能

1. **MD完整性检测**
   - 检测MD模拟是否完整运行完成
   - 验证关键输出文件的存在性（md.gro, md.xtc/trr, md.log等）
   - 分析模拟时长和预期时间的匹配度

2. **结构质量验证** (仅限PDB文件)
   - PDB复合物结构分析和链数检测
   - 异常结构识别（链数异常、原子坐标异常等）
   - 体系完整性验证
   - 注意：GRO文件不包含链信息，无法进行链分析

3. **批次追踪和统计**
   - PDB处理次数统计和重复检测
   - 批次间数据对比分析
   - 遗漏数据识别和报告

4. **质量报告生成**
   - 生成详细的质量分析报告
   - 批次统计汇总和可视化
   - 异常数据标记和建议

#### 模块架构

```
analysis/
└── quality/
    ├── __init__.py
    ├── md_completeness.py      # MDCompletenessChecker
    ├── structure_validator.py  # StructureValidator  
    ├── batch_tracker.py       # BatchTracker
    └── quality_reporter.py    # QualityReporter
```

#### 核心类设计

**MDCompletenessChecker**: MD完整性检测器
- 检测MD模拟完成状态
- 验证输出文件完整性
- 分析模拟时长和质量指标

**StructureValidator**: 结构质量验证器 (仅限PDB文件)
- PDB文件结构分析和链数检测
- 蛋白质链异常识别 (GRO文件不包含链信息)
- 复合物完整性验证

**BatchTracker**: 批次追踪器
- PDB处理记录管理
- 重复和遗漏检测
- 批次间数据统计

**QualityReporter**: 质量报告生成器
- 综合质量分析报告
- 统计数据可视化
- 问题诊断和建议

#### 质量检查标准

1. **MD完整性标准**
   - md.gro文件存在且非空
   - 轨迹文件大小合理（> 1MB）
   - 日志文件显示正常结束
   - 模拟时长达到预期

2. **结构质量标准**
   - 标准体系：5条蛋白质链
   - 异常体系：链数 != 5
   - 原子坐标范围合理
   - 缺失残基比例 < 5%

3. **数据完整性标准**
   - 关键文件完整存在
   - 文件大小在合理范围
   - 时间戳逻辑正确

#### 报告输出格式

**质量检查报告**:
```json
{
  "summary": {
    "total_jobs": 150,
    "completed": 145,
    "failed": 3,
    "incomplete": 2,
    "success_rate": 96.7
  },
  "pdb_statistics": {
    "unique_pdbs": 75,
    "multiple_runs": 12,
    "missing_pdbs": 1
  },
  "issues": [
    {
      "pdb_id": "1ABC",
      "issue_type": "incomplete_md",
      "description": "Missing md.gro file",
      "severity": "high"
    }
  ]
}
```

#### 使用示例

```python
from aftermd.analysis.quality import (
    MDCompletenessChecker, StructureValidator, 
    BatchTracker, QualityReporter
)

# MD完整性检查
completeness_checker = MDCompletenessChecker(
    min_simulation_time_ps=5000.0
)
completeness_results = completeness_checker.batch_check("/path/to/md/results")

# PDB结构验证 (仅限PDB文件)
structure_validator = StructureValidator(expected_chain_count=5)
pdb_files = ["/path/to/structure1.pdb", "/path/to/structure2.pdb"]
validation_results = structure_validator.batch_validate(pdb_files)

# 批次追踪
batch_tracker = BatchTracker()
batch_data = batch_tracker.analyze_batch("/path/to/batch")

# 生成综合报告
reporter = QualityReporter()
comprehensive_report = reporter.generate_comprehensive_report(
    completeness_results, validation_results, batch_data
)
```

#### 可能面临的问题及解决方案

1. **PDB结构包含多个复合物**
   - **问题**: 复合物链数异常，正常体系为5条链
   - **检测**: 通过MDAnalysis分析PDB文件链数
   - **处理**: 标记异常结构，生成警告报告

2. **MD模拟未完整运行**  
   - **问题**: MD因异常中断，输出不完整
   - **检测**: 检查md.gro存在性和日志文件状态
   - **处理**: 标记为失败任务，建议重新运行

3. **批次数据重复和遗漏**
   - **问题**: 同一PDB多次处理或某些PDB被遗漏
   - **检测**: 基于PDB ID（文件名前4字符）统计
   - **处理**: 生成重复/遗漏报告，建议数据清理

### 预处理步骤
对MD的原始轨迹去除pbd，具体的流程可以是参考：[pbc_process.sh](./pbc_process.sh)，这是一个比较全面的处理轨迹的方式，可以解决大部分的pbc的问题

### 轨迹的分析
1.对于轨迹RMSD的计算并且绘图

### 静态结构分析
1.PBD文件原子的B-facter的信息获取并且绘图的展示

## 新增模块 (2025-12)

### Analysis/Angles (角度分析模块)

专门用于分子动力学轨迹中的角度相关分析，特别针对TCR-pMHC复合物对接姿态研究。

**核心类**:
- `PrincipalAxesCalculator`: 主轴计算基础类，基于惯性张量对角化
- `DockingAngleAnalyzer`: TCR-pMHC对接角度分析 (Twist/Tilt/Swing三种角度)
- `DihedralCalculator`: 通用二面角计算器
- `VectorAngleCalculator`: 向量夹角计算工具类

**模块架构**:
```
analysis/angles/
├── __init__.py
├── principal_axes.py      # 主轴计算基础类 (~150行)
├── docking_angles.py      # TCR-pMHC对接角度 (~250行)
├── dihedral.py            # 通用二面角 (~100行)
└── vector_angles.py       # 向量夹角 (~100行)
```

**主要功能**:

1. **主轴计算** (principal_axes.py)
   - 基于惯性张量的主轴计算 (复用GeometryAnalyzer方法)
   - 支持单帧和整条轨迹的主轴演化
   - 返回特征值和特征向量 (按惯性矩降序排列)

2. **TCR-pMHC对接角度** (docking_angles.py)
   - **Twist角**: TCR二硫键连线在MHC平面上投影与MHC主轴的夹角
   - **Tilt角**: TCR主轴在MHC平面上投影与MHC主轴的夹角
   - **Swing角**: TCR连线相对peptide的侧向偏移角度
   - 源自VMD脚本，完全Python实现

3. **二面角计算** (dihedral.py)
   - 支持任意4个原子定义的二面角
   - 主链φ/ψ角和侧链旋转角
   - 返回-180到180度范围

4. **向量夹角** (vector_angles.py)
   - 主轴-主轴夹角
   - 质心连线-参考方向夹角
   - 灵活的向量角度计算接口

**使用示例**:
```python
from aftermd.analysis.angles import DockingAngleAnalyzer

analyzer = DockingAngleAnalyzer("md.tpr", "md_pbc.xtc")
times, twist, tilt, swing = analyzer.calculate_docking_angles_trajectory(
    mhc_selection="segname PROA and resid 50:86 140:176",
    tcr_alpha_cys="segname PROD and resid 89:94 20:25 and resname CYS and name CA",
    tcr_beta_cys="segname PROE and resid 89:94 20:25 and resname CYS and name CA",
    tcr_v_selection="segname PROD PROE and resid 3:115",
    peptide_selection="segname PROC",
    output_file="docking_angles.csv"
)
```

**技术特点**:
- 纯Python实现，基于MDAnalysis和NumPy
- 复用现有GeometryAnalyzer的惯性张量计算 (geometry.py:55-68)
- 向量投影和几何计算优化
- 输出XVG/CSV格式，兼容GROMACS工具链

---

### Analysis/Interface (界面分析模块)

用于分析蛋白质-蛋白质相互作用界面的物理化学性质。

**核心类**:
- `BuriedSurfaceCalculator`: 埋藏表面积计算
- `ContactResidueAnalyzer`: 接触残基对统计
- `ContactAtomAnalyzer`: 接触原子详细信息
- `ContactHeatmapGenerator`: 界面接触热力图

**模块架构**:
```
analysis/interface/
├── __init__.py
├── buried_surface.py      # 埋藏表面积 (~120行)
├── contact_residues.py    # 接触残基统计 (~150行)
├── contact_atoms.py       # 接触原子详情 (~130行)
└── contact_heatmap.py     # 界面热力图 (~120行)
```

**主要功能**:

1. **埋藏表面积 (BSA)** (buried_surface.py)
   - 公式: BSA = (SASA_A + SASA_B - SASA_AB) / 2
   - 使用MDAnalysis.analysis.sasa模块
   - 探针半径1.4 Å (水分子)
   - 支持轨迹时间演化

2. **接触残基对统计** (contact_residues.py)
   - 基于距离截断 (默认4.0 Å) 识别接触
   - 统计接触频率 (接触帧数/总帧数)
   - 输出DataFrame: ResA_ID, ResA_Name, ResB_ID, ResB_Name, Frequency

3. **接触原子详细信息** (contact_atoms.py)
   - 逐帧输出所有接触原子对
   - 包含原子ID、名称、残基信息、距离
   - 使用距离矩阵加速计算
   - 支持stride参数减少计算量

4. **界面热力图** (contact_heatmap.py)
   - 生成残基-残基接触频率矩阵
   - 使用seaborn绘制热力图
   - 输出300 dpi高分辨率图像
   - 集成PlotManager样式

**使用示例**:
```python
from aftermd.analysis.interface import (
    BuriedSurfaceCalculator,
    ContactResidueAnalyzer,
    ContactHeatmapGenerator
)

# BSA分析
bsa_calc = BuriedSurfaceCalculator("md.tpr", "md.xtc")
times, bsa_values = bsa_calc.calculate_bsa_trajectory(
    selection_a="segname PROA PROB PROC",  # pMHC
    selection_b="segname PROD PROE",       # TCR
    output_file="bsa_trajectory.csv"
)

# 接触频率统计
contact_analyzer = ContactResidueAnalyzer("md.tpr", "md.xtc")
contact_df = contact_analyzer.calculate_contact_frequency(
    selection_a="segname PROD",  # TCR alpha
    selection_b="segname PROC",  # Peptide
    cutoff=4.0
)

# 热力图生成
heatmap_gen = ContactHeatmapGenerator("md.tpr", "md.xtc")
contact_matrix = heatmap_gen.generate_and_plot(
    selection_a="segname PROD",
    selection_b="segname PROC",
    output_file="contact_heatmap.png"
)
```

**性能优化**:
- 距离矩阵向量化计算 (MDAnalysis.lib.distances.distance_array)
- 轨迹分块处理 (每次1000帧)
- 流式写入大文件
- 并行计算多任务

---

### Utils/BatchAnalyzer (统一批量分析工具)

集成角度分析和界面分析的统一批量处理接口。

**文件**: `aftermd/utils/batch_analyzer.py` (~400行)

**核心类**:
```python
class BatchAnalyzer:
    def batch_docking_angles()        # 批量对接角度分析
    def batch_bsa_analysis()          # 批量BSA分析
    def batch_contact_analysis()      # 批量接触分析
    def comprehensive_analysis()      # 综合批量分析
```

**主要功能**:

1. **批量角度分析**
   - 并行处理多个MD任务的对接角度
   - 自动创建任务输出目录
   - 错误隔离：单任务失败不影响其他

2. **批量界面分析**
   - 支持frequency/atoms/heatmap三种模式
   - 统一的任务配置格式
   - 生成综合分析报告

3. **综合分析**
   - 灵活的模块组合 (docking_angles, bsa, contact_frequency, contact_heatmap)
   - 基于ProcessPoolExecutor的多进程并行
   - 详细的处理日志和错误报告

**使用示例**:
```python
from aftermd.utils import BatchAnalyzer

tasks = [
    {
        'name': '1AO7',
        'topology': 'data/1AO7/md.tpr',
        'trajectory': 'data/1AO7/md_pbc.xtc',
        'mhc_selection': '...',
        'tcr_alpha_cys': '...',
        # ... 其他选择参数
    },
    # 更多任务...
]

batch_analyzer = BatchAnalyzer()
results = batch_analyzer.comprehensive_analysis(
    task_list=tasks,
    analysis_modules=['docking_angles', 'bsa', 'contact_frequency', 'contact_heatmap'],
    max_workers=4,
    output_dir="./results"
)
```

**设计特点**:
- 基于现有BatchProcessor模式 (batch_processor.py)
- 完善的错误处理和日志记录
- 支持YAML/JSON配置文件
- 生成汇总统计报告

**输出目录结构**:
```
results/
└── task_name/
    ├── angles/
    │   ├── docking_angles.csv
    │   ├── twist_angle.xvg
    │   ├── tilt_angle.xvg
    │   └── swing_angle.xvg
    ├── interface/
    │   ├── bsa_trajectory.csv
    │   ├── contact_frequency.csv
    │   ├── contact_atoms.csv
    │   └── contact_heatmap.png
    └── plots/
        ├── docking_angles.png
        └── bsa_evolution.png
```

---

### 模块依赖关系

```
principal_axes.py (基础)
    ├── docking_angles.py (继承)
    ├── vector_angles.py (继承)
    └── 独立使用

dihedral.py (独立)

buried_surface.py (独立)
contact_residues.py (独立)
contact_atoms.py (独立)
contact_heatmap.py (依赖contact_residues算法)

batch_analyzer.py (调用所有上述模块)
```

### 数据流扩展

```
Processed trajectories
    ├── Angle Analysis
    │   ├── Principal axes calculation
    │   ├── Docking angles (Twist/Tilt/Swing)
    │   ├── Dihedral angles
    │   └── Vector angles
    └── Interface Analysis
        ├── Buried surface area (BSA)
        ├── Contact residue statistics
        ├── Contact atom details
        └── Contact heatmap visualization
            ↓
    PlotManager (统一可视化)
            ↓
    Batch analysis reports
```

### 参考来源

这些模块的算法源自：
- VMD脚本: `development/angle_analysis/Analysis-scripts/`
  - mkvmd_twistangle.sh
  - mkvmd_tiltangle.sh
  - mkvmd_swingangle.sh
  - mkvmd_buriedarea.sh
  - mkvmd_PairContact.sh

完全用Python/MDAnalysis重新实现，无需VMD依赖。


# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Status

This appears to be a new or empty repository under `/Users/xumingyu/Documents/work/research/AfrerMD`. The project name suggests it may be related to "After Markdown" or markdown processing research.

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
1. Raw MD data → PBCProcessor (PBC removal and fitting)
2. Processed trajectories → Individual trajectory analysis modules
3. PDB structures → Individual structure analysis modules
4. Analysis results → PlotManager (visualization)
5. BatchProcessor coordinates parallel execution across all modules

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
### 预处理步骤
对MD的原始轨迹去除pbd，具体的流程可以是参考：[pbc_process.sh](./pbc_process.sh)，这是一个比较全面的处理轨迹的方式，可以解决大部分的pbc的问题

### 轨迹的分析
1.对于轨迹RMSD的计算并且绘图

### 静态结构分析
1.PBD文件原子的B-facter的信息获取并且绘图的展示


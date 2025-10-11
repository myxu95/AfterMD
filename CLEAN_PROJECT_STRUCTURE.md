# AfterMD 项目结构

## 清理后的完整项目结构

```
AfterMD/
├── README.md                   # 项目主要文档
├── CLAUDE.md                   # Claude代码助手指令
├── LICENSE                     # 许可证文件
├── setup.py                    # Python包安装配置
├── requirements.txt            # Python依赖
├── environment.yml             # Conda环境配置
├── .gitignore                 # Git忽略文件
│
├── aftermd/                    # 主包目录
│   ├── __init__.py            # 包初始化
│   ├── batch_process.py       # 批量处理主程序
│   │
│   ├── analysis/              # 分析模块
│   │   ├── __init__.py
│   │   ├── analysis_pipeline.py    # 轨迹分析流水线
│   │   ├── comparative_analysis.py # 对比分析模块
│   │   │
│   │   ├── quality/           # 质量检查模块
│   │   │   ├── __init__.py
│   │   │   ├── md_completeness.py    # MD完整性检查
│   │   │   ├── structure_validator.py # 结构验证
│   │   │   ├── batch_tracker.py      # 批次跟踪
│   │   │   └── quality_reporter.py   # 质量报告
│   │   │
│   │   ├── trajectory/        # 轨迹分析模块
│   │   │   ├── __init__.py
│   │   │   ├── rmsd.py        # RMSD计算
│   │   │   ├── radius_gyration.py # 回转半径计算
│   │   │   ├── distance.py    # 距离分析
│   │   │   ├── rdf.py         # 径向分布函数
│   │   │   └── hydrogen_bonds.py # 氢键分析
│   │   │
│   │   └── structure/         # 结构分析模块
│   │       ├── __init__.py
│   │       ├── bfactor.py     # B因子分析
│   │       ├── contact_map.py # 接触图计算
│   │       ├── geometry.py    # 几何分析
│   │       └── atom_info.py   # 原子信息提取
│   │
│   ├── preprocessing/         # 预处理模块
│   │   ├── __init__.py
│   │   └── pbc_processor.py   # PBC处理
│   │
│   └── utils/                 # 工具模块
│       ├── __init__.py
│       ├── batch_processor.py # 批处理器
│       ├── path_manager.py    # 路径管理
│       ├── plotting.py        # 绘图工具
│       ├── group_selector.py  # 组选择器
│       ├── slurm_generator.py # SLURM脚本生成
│       └── cleanup_manager.py # 清理管理
│
├── scripts/                   # 执行脚本
│   ├── md_workflow.py         # 主工作流脚本
│   ├── md_quality_check.py    # 质量检查脚本
│   ├── generate_slurm.py      # SLURM脚本生成
│   ├── run_trajectory_analysis.py     # 轨迹分析脚本
│   └── run_comparative_analysis.py   # 对比分析脚本
│
├── examples/                  # 使用示例
│   ├── basic_usage.py         # 基本用法示例
│   ├── complete_analysis_example.py # 完整分析示例
│   ├── modular_usage.py       # 模块化使用示例
│   ├── quality_analysis_usage.py    # 质量分析示例
│   ├── gro_copy_usage.py      # GRO文件处理示例
│   └── group_selector_usage.py      # 组选择器示例
│
├── docs/                      # 文档目录
│   ├── TRAJECTORY_ANALYSIS_GUIDE.md # 轨迹分析指南
│   └── COMPARATIVE_ANALYSIS_GUIDE.md # 对比分析指南
│
└── templates/                 # 模板文件
    └── slurm_template.sh      # SLURM作业模板
```

## 主要改进

### 1. 清理内容
- ✅ 删除所有 `__pycache__/` 目录和 `.pyc` 文件
- ✅ 删除临时脚本和修复文件
- ✅ 删除重复的质量分析脚本
- ✅ 删除测试和调试文件
- ✅ 删除空的质量报告目录

### 2. 结构重组
- ✅ 创建 `scripts/` 目录，集中管理执行脚本
- ✅ 创建 `docs/` 目录，集中管理文档
- ✅ 创建 `templates/` 目录，集中管理模板文件
- ✅ 保持 `aftermd/` 作为核心Python包

### 3. .gitignore更新
- ✅ 添加AfterMD特定的忽略规则
- ✅ 忽略MD文件类型（.xtc, .trr, .gro, .pdb等）
- ✅ 忽略分析结果目录和SLURM输出
- ✅ 忽略临时文件和备份文件

### 4. 模块化设计
- ✅ 清晰的包结构：`analysis/`, `preprocessing/`, `utils/`
- ✅ 功能性子模块：`quality/`, `trajectory/`, `structure/`
- ✅ 独立的脚本和示例目录

## 使用方式

### 核心脚本位置
- **主工作流**: `scripts/md_workflow.py`
- **质量检查**: `scripts/md_quality_check.py`
- **轨迹分析**: `scripts/run_trajectory_analysis.py`
- **对比分析**: `scripts/run_comparative_analysis.py`

### 文档位置
- **轨迹分析指南**: `docs/TRAJECTORY_ANALYSIS_GUIDE.md`
- **对比分析指南**: `docs/COMPARATIVE_ANALYSIS_GUIDE.md`
- **项目说明**: `README.md`

### 示例代码位置
- **基本使用**: `examples/basic_usage.py`
- **完整流程**: `examples/complete_analysis_example.py`

## 下一步开发

1. **功能扩展**: 基于现有模块化结构添加新分析类型
2. **性能优化**: 利用批处理和并行处理能力
3. **用户界面**: 基于脚本开发GUI或Web界面
4. **文档完善**: 补充API文档和更多使用示例
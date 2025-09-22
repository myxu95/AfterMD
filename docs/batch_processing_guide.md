# AfterMD 批量处理指南

## 概述

AfterMD 提供了强大的批量处理功能，可以自动发现和处理多个MD模拟任务。系统采用智能文件发现机制，支持不同的目录组织结构，并能够并行处理多个任务以提高效率。

## 文件发现规则

### 输入文件命名约定

批量处理系统专门寻找以下文件名：
- **轨迹文件**: `md.xtc`
- **拓扑文件**: `md.tpr`

### 搜索策略

系统按照以下优先级顺序搜索文件：

1. **第一优先级：任务目录直接搜索**
   ```
   task_directory/
   ├── md.xtc
   └── md.tpr
   ```

2. **第二优先级：prod子文件夹搜索**
   ```
   task_directory/
   └── prod/
       ├── md.xtc
       └── md.tpr
   ```

### 完整性验证

- 两个文件（`md.xtc` 和 `md.tpr`）必须同时存在
- 如果某个位置文件不完整，系统会自动检查下一个搜索位置
- 只有当找到完整的文件对时，该任务才被认为是有效的

## 使用示例

### 基本批量处理

```python
from aftermd.utils.batch_processor import BatchProcessor

# 初始化批量处理器
batch_processor = BatchProcessor(max_workers=4)

# 发现所有有效的MD任务
discovered_tasks = batch_processor.discover_batch_tasks("/path/to/simulations")

# 执行批量PBC处理
results = batch_processor.batch_pbc_processing(
    base_directory="/path/to/simulations",
    output_base_dir="/path/to/output",
    dt=10.0,  # 可选：轨迹降采样
    gmx_executable="gmx"
)

print(f"成功处理: {results['successful']}/{results['total_tasks']} 个任务")
```

### 单个任务文件发现

```python
# 测试单个任务目录
md_files = batch_processor.find_md_input_files("/path/to/single_task")

if md_files:
    trajectory, topology = md_files
    print(f"轨迹文件: {trajectory}")
    print(f"拓扑文件: {topology}")
else:
    print("未找到完整的MD文件集")
```

### 自定义处理函数

```python
def custom_analysis_function(trajectory, topology, output_dir, **kwargs):
    """自定义分析函数示例"""
    # 执行自定义分析
    # 返回分析结果
    return {"status": "success", "output_files": ["result1.xvg", "result2.pdb"]}

# 使用自定义函数进行批量处理
task_list = list(discovered_tasks.items())
results = batch_processor.process_files(
    files=task_list,
    process_func=custom_analysis_function,
    use_multiprocessing=True,
    output_base_dir="/path/to/output"
)
```

## 目录结构示例

### 支持的目录结构

#### 结构1：直接文件组织
```
simulations/
├── antibody_complex/
│   ├── md.xtc
│   ├── md.tpr
│   └── other_files...
├── enzyme_substrate/
│   ├── md.xtc
│   ├── md.tpr
│   └── analysis/
└── membrane_protein/
    ├── md.xtc
    ├── md.tpr
    └── logs/
```

#### 结构2：prod子文件夹组织
```
simulations/
├── system_A/
│   ├── setup/
│   ├── prod/
│   │   ├── md.xtc
│   │   └── md.tpr
│   └── analysis/
├── system_B/
│   ├── equilibration/
│   ├── prod/
│   │   ├── md.xtc
│   │   └── md.tpr
│   └── results/
└── system_C/
    ├── preparation/
    └── prod/
        ├── md.xtc
        └── md.tpr
```

#### 结构3：混合组织（优先级演示）
```
simulations/
├── mixed_task1/
│   ├── md.xtc          # ← 优先选择（直接目录）
│   ├── md.tpr          # ← 优先选择（直接目录）
│   └── prod/
│       ├── md.xtc      # 备用位置
│       └── md.tpr      # 备用位置
├── mixed_task2/
│   ├── md.tpr          # 不完整（缺少md.xtc）
│   └── prod/
│       ├── md.xtc      # ← 选择这个位置（完整）
│       └── md.tpr      # ← 选择这个位置（完整）
└── invalid_task/
    ├── trajectory.xtc  # 文件名不正确
    └── topology.tpr    # 文件名不正确
```

## 批量PBC处理工作流

### 完整处理流程

1. **文件发现阶段**
   - 扫描基础目录中的所有子目录
   - 应用文件发现规则查找`md.xtc`和`md.tpr`
   - 验证文件完整性和可访问性

2. **任务准备阶段**
   - 为每个有效任务创建输出目录
   - 初始化PBC处理器
   - 准备并行处理队列

3. **智能PBC处理**（每个任务）
   - **最短链检测**: 使用两步法检测最短肽链
   - **Index文件生成**: 创建包含最短链的自定义index文件
   - **精确组号映射**: 从生成的index文件解析准确的组号
   - **三步PBC处理**:
     - Step 1: 使用最短链进行中心化（`-center -pbc atom`）
     - Step 2: 应用PBC whole保持分子完整性（`-pbc whole`）
     - Step 3: 使用backbone进行旋转平移拟合（`-fit rot+trans`）

4. **结果汇总阶段**
   - 收集所有任务的处理结果
   - 生成成功/失败统计
   - 创建详细的处理报告

### 输出结构

```
output_directory/
├── task1_antibody_complex/
│   ├── task1_antibody_complex_processed.xtc
│   └── processing_metadata.json
├── task2_enzyme_substrate/
│   ├── task2_enzyme_substrate_processed.xtc
│   └── processing_metadata.json
├── task3_membrane_protein/
│   ├── task3_membrane_protein_processed.xtc
│   └── processing_metadata.json
└── batch_processing_summary.json
```

## 高级功能

### 并行处理配置

```python
# 配置并行处理参数
batch_processor = BatchProcessor(max_workers=8)  # 使用8个工作进程

# 或者根据系统自动配置
import multiprocessing
batch_processor = BatchProcessor(max_workers=multiprocessing.cpu_count())
```

### 轨迹降采样

```python
# 处理时应用轨迹降采样
results = batch_processor.batch_pbc_processing(
    base_directory="/path/to/simulations",
    output_base_dir="/path/to/output",
    dt=20.0,  # 每20 ps采样一帧
    gmx_executable="gmx"
)
```

### 错误处理和恢复

```python
# 检查处理结果
for result in results['results']:
    task_name = result['task_name']
    status = result['status']
    
    if status == 'success':
        print(f"✅ {task_name}: 处理成功")
        print(f"   输出文件: {result['processed']}")
    else:
        print(f"❌ {task_name}: 处理失败")
        print(f"   错误信息: {result['error']}")
```

### 自定义GROMACS可执行文件

```python
# 使用自定义GROMACS安装
results = batch_processor.batch_pbc_processing(
    base_directory="/path/to/simulations",
    output_base_dir="/path/to/output",
    gmx_executable="/usr/local/gromacs/bin/gmx"  # 自定义路径
)
```

## 最佳实践

### 目录组织建议

1. **保持一致的命名**: 始终使用`md.xtc`和`md.tpr`作为文件名
2. **逻辑分组**: 将相关的模拟任务放在同一基础目录下
3. **备份重要数据**: 在批量处理前备份原始轨迹文件
4. **使用prod子文件夹**: 对于复杂的模拟项目，将生产运行文件放在`prod`子文件夹中

### 性能优化

1. **合理设置工作进程数**: 通常设置为CPU核心数或略少
2. **考虑内存使用**: 处理大轨迹文件时注意内存限制
3. **使用轨迹降采样**: 对于大型轨迹文件，考虑使用`dt`参数降采样
4. **并行I/O考虑**: 确保存储系统能够处理并行读写操作

### 故障排除

#### 常见问题和解决方案

1. **未找到MD文件**
   - 检查文件名是否正确（必须是`md.xtc`和`md.tpr`）
   - 验证目录结构是否符合搜索规则
   - 确认文件权限允许读取

2. **GROMACS命令失败**
   - 验证GROMACS安装和路径设置
   - 检查拓扑文件完整性
   - 确认轨迹文件与拓扑文件匹配

3. **内存不足错误**
   - 减少并行工作进程数
   - 使用轨迹降采样减少内存使用
   - 分批处理大型数据集

## 命令行工具

### 批量处理演示脚本

```bash
# 运行批量处理演示
python examples/batch_pbc_processing_demo.py --base-dir /path/to/simulations

# 创建演示目录结构
python examples/batch_pbc_processing_demo.py --create-demo

# 仅进行文件发现测试
python examples/md_file_discovery_demo.py /path/to/simulations
```

### 实际使用示例

```bash
# 批量处理所有任务，使用降采样
python -c "
from aftermd.utils.batch_processor import BatchProcessor
bp = BatchProcessor(max_workers=4)
results = bp.batch_pbc_processing(
    '/data/md_simulations',
    '/data/processed_output',
    dt=10.0
)
print(f'处理完成: {results[\"successful\"]}/{results[\"total_tasks\"]} 任务')
"
```

## 集成到工作流

### 与其他AfterMD工具集成

```python
from aftermd.utils.batch_processor import BatchProcessor
from aftermd.analysis.rmsd_analyzer import RMSDAnalyzer
from aftermd.analysis.rdf_analyzer import RDFAnalyzer

# 第一步：批量PBC处理
batch_processor = BatchProcessor()
pbc_results = batch_processor.batch_pbc_processing(
    base_directory="/path/to/raw_simulations",
    output_base_dir="/path/to/pbc_processed"
)

# 第二步：批量RMSD分析
for result in pbc_results['results']:
    if result['status'] == 'success':
        processed_traj = result['processed']
        topology = result['topology']  # 从原始结果获取
        
        rmsd_analyzer = RMSDAnalyzer(processed_traj, topology)
        rmsd_analyzer.calculate_rmsd(output_dir=f"/path/to/analysis/{result['task_name']}")

# 第三步：批量RDF分析
# 类似的方式处理其他分析...
```

## 总结

AfterMD的批量处理功能提供了：

✅ **智能文件发现**: 自动识别MD文件的多种组织结构  
✅ **最短链检测**: 集成的智能PBC处理策略  
✅ **并行处理**: 高效的多任务并行执行  
✅ **错误处理**: 完善的错误报告和恢复机制  
✅ **灵活配置**: 支持各种自定义参数和工作流  

这使得处理大量MD模拟数据变得简单、高效且可靠。
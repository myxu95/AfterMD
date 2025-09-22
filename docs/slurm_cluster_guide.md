# AfterMD 集群批量处理指南

## 概述

AfterMD 提供了强大的 SLURM 脚本生成功能，可以将大量的 MD 模拟任务自动分组并生成相应的集群作业脚本。这使得在 HPC 集群上批量处理数百个 MD 轨迹变得简单高效。

## 核心功能

### 🔄 任务分组
- 自动将大量 MD 任务分组成合适的批次
- 例如：20 个任务 → 每批 10 个 → 生成 2 个 SLURM 脚本
- 支持不均匀分组（如 23 个任务 → 2×10 + 1×3）

### ⚙️ SLURM 参数定制
- 灵活配置分区、时间、CPU、GPU 等资源
- 支持 conda 环境激活
- 自动处理 GROMACS 模块加载
- 基于模板的脚本生成

### 🚀 自动化提交
- 生成批量提交脚本
- 可配置作业间提交延迟
- 集成作业状态监控命令

## 基本使用方法

### 1. 最简单的使用

```python
from aftermd.utils.slurm_generator import generate_slurm_scripts_for_md_tasks

# 为 20 个 MD 任务生成 SLURM 脚本，每批处理 10 个
results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/md_simulations",
    tasks_per_batch=10
)

print(f"生成了 {len(results['scripts'])} 个 SLURM 脚本")
```

### 2. 自定义集群参数

```python
# 定义集群参数
slurm_params = {
    "partition": "gpu",
    "time": "24:00:00", 
    "cpus_per_task": 16,
    "gres": "gpu:2",
    "memory": "64G",
    "conda_env": "aftermd"
}

results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/md_simulations",
    tasks_per_batch=5,
    slurm_params=slurm_params,
    dt=10.0  # 轨迹降采样，每 10 ps 一帧
)
```

### 3. 命令行使用

```bash
# 基本用法
python -m aftermd.utils.slurm_generator /data/simulations --tasks-per-batch 10

# 自定义参数
python -m aftermd.utils.slurm_generator /data/simulations \
  --tasks-per-batch 5 \
  --partition gpu \
  --time 24:00:00 \
  --cpus-per-task 16 \
  --conda-env aftermd

# 使用自定义模板
python -m aftermd.utils.slurm_generator /data/simulations \
  --template my_slurm_template.sh
```

## 任务分组策略

### 分组计算示例

| 总任务数 | 每批任务数 | 生成脚本数 | 分组结果 |
|---------|----------|----------|----------|
| 20      | 10       | 2        | 10 + 10 |
| 23      | 10       | 3        | 10 + 10 + 3 |
| 100     | 15       | 7        | 6×15 + 1×10 |
| 5       | 10       | 1        | 5 |

### 选择合适的批次大小

**考虑因素：**
- **集群资源**: CPU/GPU 核心数，内存大小
- **任务复杂度**: 轨迹大小，处理时间
- **队列限制**: 最大作业时间，并发作业数
- **存储 I/O**: 避免过多并发读写

**推荐策略：**
```python
# 小任务（< 1GB 轨迹）
tasks_per_batch = 20

# 中等任务（1-5GB 轨迹）
tasks_per_batch = 10

# 大任务（> 5GB 轨迹）
tasks_per_batch = 5

# 超大任务（> 20GB 轨迹）
tasks_per_batch = 1
```

## 生成的脚本结构

### 典型 SLURM 脚本内容

```bash
#!/bin/bash
#SBATCH --job-name=aftermd_batch_001
#SBATCH -p gpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:2
#SBATCH --mem=64G

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate aftermd
echo "Activated conda environment: aftermd"

# Environment setup
export LD_LIBRARY_PATH=/public/software/lib/:$LD_LIBRARY_PATH
source /public/software/profile.d/apps_gromacs_2023.2.sh

# Job information
echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST" 
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"
echo "Processing 5 MD tasks in this batch"
echo ""

# AfterMD batch processing
python -m aftermd.batch_process "/data/md_simulations" \
  --output "/data/processed" \
  --dt 10.0 \
  --verbose

echo ""
echo "End time: $(date)"
echo "Batch processing completed"
```

### 提交脚本 (submit_all_batches.sh)

```bash
#!/bin/bash
# AfterMD Batch Job Submission Script

echo "Submitting 4 AfterMD batch jobs..."

# Submit batch 1/4
echo "Submitting aftermd_batch_001.sh..."
sbatch "aftermd_batch_001.sh"
sleep 5

# Submit batch 2/4
echo "Submitting aftermd_batch_002.sh..."
sbatch "aftermd_batch_002.sh"
sleep 5

# ... (继续其他批次)

echo "All 4 batch jobs submitted!"
echo "Use 'squeue -u $USER' to check job status"
```

## 实际使用工作流

### 步骤 1: 准备数据

```bash
# 确保 MD 数据组织正确
ls /data/md_simulations/
# task_001/md.xtc, md.tpr
# task_002/prod/md.xtc, md.tpr
# task_003/md.xtc, md.tpr
# ...
```

### 步骤 2: 生成 SLURM 脚本

```python
from aftermd.utils.slurm_generator import generate_slurm_scripts_for_md_tasks

results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/md_simulations",
    tasks_per_batch=10,
    output_script_dir="/home/user/slurm_scripts",
    slurm_params={
        "partition": "gpu",
        "time": "12:00:00",
        "conda_env": "aftermd"
    }
)
```

### 步骤 3: 检查生成的脚本

```bash
cd /home/user/slurm_scripts
ls -la
# aftermd_batch_001.sh
# aftermd_batch_002.sh  
# submit_all_batches.sh

# 检查脚本内容
head -20 aftermd_batch_001.sh
```

### 步骤 4: 提交作业

```bash
# 提交所有批次
bash submit_all_batches.sh

# 或者逐个提交
sbatch aftermd_batch_001.sh
sbatch aftermd_batch_002.sh
```

### 步骤 5: 监控作业

```bash
# 查看作业状态
squeue -u $USER

# 查看作业详情
scontrol show job <job_id>

# 查看输出日志
tail -f slurm-<job_id>.out
```

### 步骤 6: 收集结果

```bash
# 检查处理结果
ls /data/processed/
# task_001/task_001_processed.xtc
# task_002/task_002_processed.xtc
# ...

# 检查处理日志
grep -r "successfully" /data/processed/*/processing_metadata.json
```

## 高级配置

### 自定义 SLURM 模板

创建自定义模板文件 `my_template.sh`:

```bash
#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
{gres_line}
{memory_line}
#SBATCH --output=aftermd_%j.out
#SBATCH --error=aftermd_%j.err

{conda_activation}

# 自定义环境设置
module load gromacs/2023.2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 处理命令
{processing_commands}
```

使用自定义模板：

```python
results = generate_slurm_scripts_for_md_tasks(
    simulations_path="/data/simulations",
    template_file="my_template.sh",
    tasks_per_batch=10
)
```

### 不同集群配置示例

#### GPU 集群配置
```python
gpu_params = {
    "partition": "gpu",
    "time": "24:00:00",
    "cpus_per_task": 8,
    "gres": "gpu:1",
    "memory": "32G"
}
```

#### CPU 密集型配置
```python
cpu_params = {
    "partition": "cpu",
    "time": "48:00:00", 
    "cpus_per_task": 32,
    "memory": "128G"
}
```

#### 快速队列配置
```python
fast_params = {
    "partition": "fast",
    "time": "4:00:00",
    "cpus_per_task": 16
}
```

## 故障排除

### 常见问题

1. **SLURM 参数不兼容**
   ```bash
   # 检查集群分区信息
   sinfo
   
   # 检查资源限制
   scontrol show partition <partition_name>
   ```

2. **GROMACS 模块加载失败**
   ```bash
   # 检查可用模块
   module avail gromacs
   
   # 测试模块加载
   module load gromacs/2023.2
   gmx --version
   ```

3. **Conda 环境问题**
   ```bash
   # 验证环境存在
   conda env list
   
   # 测试环境激活
   conda activate aftermd
   python -c "import aftermd; print('OK')"
   ```

4. **文件路径访问问题**
   ```bash
   # 检查路径权限
   ls -la /data/md_simulations/
   
   # 测试文件发现
   python -c "from aftermd import discover_md_tasks; print(discover_md_tasks('/data/simulations'))"
   ```

### 调试技巧

1. **先测试小批次**
   ```python
   # 先生成 1 个任务的脚本进行测试
   results = generate_slurm_scripts_for_md_tasks(
       simulations_path="/data/simulations",
       tasks_per_batch=1
   )
   ```

2. **启用详细输出**
   ```bash
   python -m aftermd.utils.slurm_generator /data/simulations --verbose
   ```

3. **检查生成的脚本**
   ```bash
   # 手动运行生成的命令进行测试
   bash -x aftermd_batch_001.sh
   ```

## 性能优化建议

### 资源配置优化

1. **CPU 核心数**: 根据 GROMACS 性能测试确定最优核心数
2. **内存分配**: 通常 2-4GB 每个 CPU 核心
3. **GPU 配置**: 对于大系统使用 GPU 加速
4. **I/O 优化**: 避免过多并发读写同一存储

### 批次大小优化

```python
# 根据系统大小调整批次大小
def calculate_optimal_batch_size(total_tasks, avg_traj_size_gb):
    if avg_traj_size_gb < 1:
        return min(20, total_tasks)
    elif avg_traj_size_gb < 5:
        return min(10, total_tasks)
    else:
        return min(5, total_tasks)

optimal_batch = calculate_optimal_batch_size(100, 2.5)  # 100 个任务，平均 2.5GB
```

## 总结

AfterMD 的 SLURM 脚本生成器提供了：

✅ **自动化任务分组**: 智能分配大量 MD 任务到合适的批次  
✅ **灵活的集群配置**: 支持各种 HPC 集群环境  
✅ **完整的工作流集成**: 从数据发现到结果处理的端到端自动化  
✅ **强大的定制能力**: 模板系统支持各种特殊需求  
✅ **简单的使用接口**: 几行代码即可生成完整的作业脚本  

这使得在集群上批量处理数百个 MD 轨迹变得简单、高效且可靠。
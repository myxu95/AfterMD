#!/usr/bin/env python3
"""
MD轨迹处理完整工作流

集成质量检查和SLURM脚本生成的一站式工具。
先在前台进行质量检查，然后为合格的MD生成SLURM批处理脚本。

使用方法:
    python md_workflow.py /path/to/md/simulations
    python md_workflow.py /path/to/md/simulations --batch-size 15 --partition gpu
"""

import sys
import subprocess
import argparse
import time
from pathlib import Path
from datetime import datetime

def run_command(cmd, description):
    """运行命令并显示结果"""
    print(f"\n{'='*60}")
    print(f"🔄 {description}")
    print(f"💻 命令: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    start_time = time.time()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False, text=True)
        end_time = time.time()
        print(f"\n✅ {description} 完成 (耗时: {end_time - start_time:.1f}秒)")
        return True
    except subprocess.CalledProcessError as e:
        end_time = time.time()
        print(f"\n❌ {description} 失败 (耗时: {end_time - start_time:.1f}秒)")
        print(f"错误代码: {e.returncode}")
        return False
    except Exception as e:
        end_time = time.time()
        print(f"\n❌ {description} 出错 (耗时: {end_time - start_time:.1f}秒)")
        print(f"错误信息: {e}")
        return False

def main():
    """主工作流"""
    parser = argparse.ArgumentParser(
        description="MD轨迹处理完整工作流 - 质量检查 + SLURM脚本生成",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
工作流程:
  1. 前台快速质量检查所有MD
  2. 识别质量合格的MD (链数正确、模拟完整)
  3. 只为合格MD生成SLURM批处理脚本
  4. 提供作业提交指引

使用示例:
  基本用法:
    python md_workflow.py /data/md_simulations
    
  自定义参数:
    python md_workflow.py /data/md_simulations \\
        --expected-chains 3 \\
        --batch-size 15 \\
        --partition gpu \\
        --time 24:00:00
    
  仅质量检查:
    python md_workflow.py /data/md_simulations --quality-only
    
  仅生成SLURM脚本 (需要已有质量检查结果):
    python md_workflow.py /data/md_simulations --slurm-only
        """
    )
    
    # 必需参数
    parser.add_argument(
        "md_root",
        help="MD模拟根目录路径"
    )
    
    # 工作流控制
    workflow_group = parser.add_mutually_exclusive_group()
    workflow_group.add_argument(
        "--quality-only",
        action="store_true",
        help="仅执行质量检查，不生成SLURM脚本"
    )
    workflow_group.add_argument(
        "--slurm-only", 
        action="store_true",
        help="仅生成SLURM脚本 (需要已有质量检查结果)"
    )
    
    # 质量检查参数
    quality_group = parser.add_argument_group("质量检查参数")
    quality_group.add_argument(
        "--expected-chains",
        type=int,
        default=5,
        help="期望的蛋白质链数 (默认: 5)"
    )
    quality_group.add_argument(
        "--min-traj-size",
        type=float,
        default=1.0,
        help="最小轨迹文件大小(MB) (默认: 1.0)"
    )
    quality_group.add_argument(
        "--min-sim-time",
        type=float,
        default=5000.0,
        help="最小模拟时间(ps) (默认: 5000.0)"
    )
    
    # SLURM脚本参数
    slurm_group = parser.add_argument_group("SLURM脚本参数")
    slurm_group.add_argument(
        "--batch-size",
        type=int,
        default=10,
        help="每个SLURM作业的MD数量 (默认: 10)"
    )
    slurm_group.add_argument(
        "--partition",
        default="quick",
        help="SLURM分区 (默认: quick)"
    )
    slurm_group.add_argument(
        "--time",
        default="12:00:00", 
        help="作业时间限制 (默认: 12:00:00)"
    )
    slurm_group.add_argument(
        "--cpus",
        type=int,
        default=11,
        help="每个任务的CPU核数 (默认: 11)"
    )
    slurm_group.add_argument(
        "--gpu",
        help="GPU需求 (例如: gpu:1)"
    )
    
    # 输出参数
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="./md_workflow_results",
        help="工作流结果输出目录 (默认: ./md_workflow_results)"
    )
    
    # 其他选项
    parser.add_argument(
        "--force",
        action="store_true",
        help="强制覆盖已存在的结果"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="显示详细信息"
    )
    
    args = parser.parse_args()
    
    # 验证输入目录
    md_root = Path(args.md_root)
    if not md_root.exists():
        print(f"❌ 错误: MD目录不存在 - {md_root}")
        return 1
    
    if not md_root.is_dir():
        print(f"❌ 错误: 路径不是目录 - {md_root}")
        return 1
    
    # 创建输出目录
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("🚀 MD轨迹处理工作流")
    print("=" * 60)
    print(f"📂 MD根目录: {md_root}")
    print(f"📁 输出目录: {output_dir}")
    print(f"⏰ 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"🔧 期望链数: {args.expected_chains}")
    print(f"📦 批次大小: {args.batch_size}")
    print(f"🖥️  SLURM分区: {args.partition}")
    print("=" * 60)
    
    # 定义路径
    quality_results_dir = output_dir / "quality_check_results"
    qualified_list_file = quality_results_dir / "qualified_mds.txt"
    slurm_scripts_dir = output_dir / "slurm_scripts"
    
    workflow_success = True
    
    # 步骤1: 质量检查
    if not args.slurm_only:
        quality_cmd = [
            sys.executable, "md_quality_check.py",
            str(md_root),
            "--expected-chains", str(args.expected_chains),
            "--min-traj-size", str(args.min_traj_size),
            "--min-sim-time", str(args.min_sim_time),
            "--output", str(quality_results_dir)
        ]
        
        if not run_command(quality_cmd, "MD质量检查"):
            workflow_success = False
            print("\n❌ 质量检查失败，工作流中止")
            return 1
        
        # 检查是否有合格的MD
        if not qualified_list_file.exists():
            print(f"\n❌ 未找到合格MD列表文件: {qualified_list_file}")
            return 1
        
        # 读取合格MD数量
        try:
            with open(qualified_list_file, 'r') as f:
                qualified_count = len([line.strip() for line in f if line.strip()])
            
            if qualified_count == 0:
                print(f"\n❌ 没有发现合格的MD，无法继续生成SLURM脚本")
                print(f"📋 请检查质量报告: {quality_results_dir}")
                return 1
            
            print(f"\n✅ 发现 {qualified_count} 个合格MD，可以继续生成SLURM脚本")
            
        except Exception as e:
            print(f"\n❌ 读取合格MD列表时出错: {e}")
            return 1
    
    # 仅质量检查模式
    if args.quality_only:
        print(f"\n🎯 质量检查完成")
        print(f"📋 查看结果: {quality_results_dir}")
        return 0
    
    # 步骤2: 生成SLURM脚本
    if not args.slurm_only and not qualified_list_file.exists():
        print(f"\n❌ 未找到质量检查结果，请先运行质量检查")
        print(f"建议: python md_workflow.py {md_root} --quality-only")
        return 1
    
    slurm_cmd = [
        sys.executable, "generate_slurm.py",
        str(md_root),
        "--qualified-list", str(qualified_list_file),
        "--batch-size", str(args.batch_size),
        "--output-scripts", str(slurm_scripts_dir),
        "--partition", args.partition,
        "--time", args.time,
        "--cpus", str(args.cpus)
    ]
    
    # 添加可选参数
    if args.gpu:
        slurm_cmd.extend(["--gpu", args.gpu])
    
    if args.verbose:
        slurm_cmd.append("--verbose")
    
    if not run_command(slurm_cmd, "SLURM脚本生成"):
        workflow_success = False
        print("\n❌ SLURM脚本生成失败")
        return 1
    
    # 工作流完成
    end_time = datetime.now()
    
    print(f"\n{'='*60}")
    print("🎉 MD轨迹处理工作流完成!")
    print(f"{'='*60}")
    print(f"⏰ 完成时间: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📁 结果目录: {output_dir}")
    print(f"\n📋 生成的文件:")
    
    if quality_results_dir.exists():
        print(f"   🔍 质量检查: {quality_results_dir}")
        for result_file in ["quality_summary.txt", "qualified_mds.txt", "failed_mds.txt"]:
            file_path = quality_results_dir / result_file
            if file_path.exists():
                print(f"      • {result_file}")
    
    if slurm_scripts_dir.exists():
        print(f"   🖥️  SLURM脚本: {slurm_scripts_dir}")
        script_files = list(slurm_scripts_dir.glob("*.sh"))
        for script_file in sorted(script_files)[:5]:  # 只显示前5个
            print(f"      • {script_file.name}")
        if len(script_files) > 5:
            print(f"      • ... 以及其他 {len(script_files) - 5} 个脚本")
    
    print(f"\n🚀 下一步操作:")
    
    # 查找提交脚本
    submit_script = slurm_scripts_dir / "submit_all_batches.sh"
    if submit_script.exists():
        print(f"   1. 提交所有作业: bash {submit_script}")
        print(f"   2. 监控作业状态: squeue -u $USER")
        print(f"   3. 查看作业日志: ls {slurm_scripts_dir}/*.out")
    else:
        print(f"   1. 检查SLURM脚本: ls {slurm_scripts_dir}")
        print(f"   2. 手动提交作业: sbatch <script_name>")
    
    print(f"\n💡 提示:")
    print(f"   • 质量报告: {quality_results_dir / 'quality_summary.txt'}")
    print(f"   • 合格MD列表: {qualified_list_file}")
    print(f"   • 只有质量合格的MD会被处理，节省计算资源")
    
    return 0 if workflow_success else 1


if __name__ == "__main__":
    exit(main())
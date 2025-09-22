#!/usr/bin/env python3
"""
AfterMD SLURM脚本生成器
为MD模拟批量处理生成SLURM作业脚本

使用方法:
    python generate_slurm.py /path/to/simulations
    python generate_slurm.py /path/to/simulations --batch-size 5 --partition gpu
"""

import sys
import argparse
from pathlib import Path
from aftermd import generate_slurm_scripts_for_md_tasks

def create_parser():
    """创建命令行参数解析器"""
    parser = argparse.ArgumentParser(
        description="AfterMD SLURM脚本生成器 - 为MD模拟批量处理生成SLURM作业脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  带质量检查的推荐流程:
    # 1. 先进行质量检查
    python md_quality_check.py /data/simulations
    
    # 2. 只为合格MD生成SLURM脚本
    python generate_slurm.py /data/simulations \\
        --qualified-list ./quality_check_results/qualified_mds.txt
    
  跳过质量检查 (处理所有MD):
    python generate_slurm.py /data/simulations --skip-quality-check
    
  自定义参数:
    python generate_slurm.py /data/simulations \\
        --qualified-list ./quality_check_results/qualified_mds.txt \\
        --batch-size 8 \\
        --partition gpu \\
        --time 24:00:00
        
输出文件:
  - slurm_scripts/aftermd_batch_*.sh  # 批处理脚本
  - slurm_scripts/submit_all_batches.sh  # 批量提交脚本
        """
    )
    
    # 必需参数
    parser.add_argument(
        "simulations_path",
        help="包含MD模拟任务的目录路径"
    )
    
    # 质量过滤参数
    parser.add_argument(
        "--qualified-list",
        help="质量检查合格的MD列表文件 (来自md_quality_check.py的输出)"
    )
    
    parser.add_argument(
        "--skip-quality-check",
        action="store_true", 
        help="跳过质量检查，处理所有找到的MD目录"
    )
    
    # 批处理参数
    parser.add_argument(
        "--batch-size", "-b",
        type=int,
        default=10,
        help="每个SLURM作业处理的任务数量 (默认: 10)"
    )
    
    # 输出参数
    parser.add_argument(
        "--output-scripts", "-o",
        default="./slurm_scripts",
        help="SLURM脚本输出目录 (默认: ./slurm_scripts)"
    )
    
    parser.add_argument(
        "--output-data",
        help="数据处理结果输出目录 (默认: 自动生成)"
    )
    
    parser.add_argument(
        "--output-suffix",
        default="processed",
        help="自动生成输出目录的后缀 (默认: 'processed')"
    )
    
    # SLURM集群参数
    cluster_group = parser.add_argument_group("SLURM集群参数")
    
    cluster_group.add_argument(
        "--partition", "-p",
        default="quick",
        help="SLURM分区名称 (默认: quick)"
    )
    
    cluster_group.add_argument(
        "--time", "-t",
        default="12:00:00",
        help="作业时间限制 (默认: 12:00:00)"
    )
    
    cluster_group.add_argument(
        "--cpus",
        type=int,
        default=11,
        help="每个任务的CPU核数 (默认: 11)"
    )
    
    cluster_group.add_argument(
        "--memory",
        help="内存需求 (例如: 32G)"
    )
    
    cluster_group.add_argument(
        "--gpu",
        help="GPU需求 (例如: gpu:1)"
    )
    
    # MD处理参数
    md_group = parser.add_argument_group("MD处理参数")
    
    md_group.add_argument(
        "--dt",
        type=float,
        help="时间采样间隔 (ps)"
    )
    
    md_group.add_argument(
        "--template",
        help="自定义SLURM模板文件路径"
    )
    
    # 其他选项
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="只显示配置信息，不生成实际脚本"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="显示详细信息"
    )
    
    return parser

def validate_args(args):
    """验证命令行参数"""
    # 检查输入目录
    sim_path = Path(args.simulations_path)
    if not sim_path.exists():
        print(f"❌ 错误: 输入目录不存在: {args.simulations_path}")
        return False
    
    if not sim_path.is_dir():
        print(f"❌ 错误: 输入路径不是目录: {args.simulations_path}")
        return False
    
    # 检查质量列表文件
    if args.qualified_list:
        qualified_file = Path(args.qualified_list)
        if not qualified_file.exists():
            print(f"❌ 错误: 质量列表文件不存在: {args.qualified_list}")
            return False
        if not qualified_file.is_file():
            print(f"❌ 错误: 质量列表路径不是文件: {args.qualified_list}")
            return False
    
    # 检查批次大小
    if args.batch_size <= 0:
        print(f"❌ 错误: 批次大小必须大于0，当前值: {args.batch_size}")
        return False
    
    # 检查模板文件
    if args.template and not Path(args.template).exists():
        print(f"❌ 错误: 模板文件不存在: {args.template}")
        return False
    
    return True

def load_qualified_mds(qualified_list_file):
    """从质量检查结果文件加载合格的MD列表"""
    qualified_mds = []
    
    try:
        with open(qualified_list_file, 'r') as f:
            for line in f:
                md_path = line.strip()
                if md_path and Path(md_path).exists():
                    qualified_mds.append(Path(md_path))
                elif md_path:
                    print(f"⚠️  警告: 合格MD路径不存在，跳过: {md_path}")
        
        print(f"📋 从质量检查结果加载了 {len(qualified_mds)} 个合格MD")
        return qualified_mds
        
    except Exception as e:
        print(f"❌ 错误: 无法读取质量列表文件: {e}")
        return []

def build_slurm_params(args):
    """构建SLURM参数字典"""
    slurm_params = {
        "partition": args.partition,
        "time": args.time,
        "cpus_per_task": args.cpus
    }
    
    # 可选参数
    if args.memory:
        slurm_params["memory"] = args.memory
    
    if args.gpu:
        slurm_params["gres"] = args.gpu
    
    return slurm_params

def print_configuration(args, slurm_params, qualified_count=None):
    """打印配置信息"""
    print("🚀 AfterMD SLURM脚本生成器")
    print("=" * 50)
    print(f"📂 输入目录: {args.simulations_path}")
    
    # 质量过滤信息
    if args.qualified_list:
        print(f"🔍 质量过滤: 启用")
        print(f"📋 质量列表: {args.qualified_list}")
        if qualified_count is not None:
            print(f"✅ 合格MD数: {qualified_count} 个")
    elif args.skip_quality_check:
        print(f"🔍 质量过滤: 跳过 (处理所有MD)")
    else:
        print(f"🔍 质量过滤: 未启用 (建议先运行 md_quality_check.py)")
    
    print(f"📦 批次大小: {args.batch_size} 任务/作业")
    print(f"📁 脚本输出: {args.output_scripts}")
    
    if args.output_data:
        print(f"💾 数据输出: {args.output_data}")
    else:
        suffix = args.output_suffix
        sim_name = Path(args.simulations_path).name
        auto_output = f"{Path(args.simulations_path).parent}/{sim_name}_{suffix}"
        print(f"💾 数据输出: {auto_output} (自动生成)")
    
    print("\n🖥️  SLURM配置:")
    for key, value in slurm_params.items():
        print(f"   {key}: {value}")
    
    if args.dt:
        print(f"\n⏱️  采样间隔: {args.dt} ps")
    
    if args.template:
        print(f"\n📄 自定义模板: {args.template}")
    
    print("=" * 50)

def main():
    """主函数"""
    parser = create_parser()
    args = parser.parse_args()
    
    # 验证参数
    if not validate_args(args):
        return 1
    
    # 处理质量过滤
    qualified_mds = None
    if args.qualified_list:
        print("\n🔍 加载质量检查结果...")
        qualified_mds = load_qualified_mds(args.qualified_list)
        if not qualified_mds:
            print("❌ 没有找到合格的MD，无法生成SLURM脚本")
            return 1
    elif not args.skip_quality_check:
        print("\n⚠️  警告: 未启用质量过滤")
        print("   建议先运行: python md_quality_check.py /path/to/md/simulations")
        print("   然后使用: python generate_slurm.py /path/to/md/simulations --qualified-list ./quality_check_results/qualified_mds.txt")
        print("   或者使用 --skip-quality-check 跳过质量检查")
        
        user_input = input("\n是否继续处理所有MD? (y/N): ").strip().lower()
        if user_input not in ['y', 'yes']:
            print("操作已取消")
            return 0
    
    # 构建SLURM参数
    slurm_params = build_slurm_params(args)
    
    # 打印配置信息
    qualified_count = len(qualified_mds) if qualified_mds else None
    print_configuration(args, slurm_params, qualified_count)
    
    # 干运行模式
    if args.dry_run:
        print("\n🔍 干运行模式 - 仅显示配置，不生成脚本")
        return 0
    
    try:
        print("\n📋 开始生成SLURM脚本...")
        
        # 生成SLURM脚本 - 传递合格MD列表
        results = generate_slurm_scripts_for_md_tasks(
            simulations_path=args.simulations_path,
            tasks_per_batch=args.batch_size,
            output_script_dir=args.output_scripts,
            base_output_dir=args.output_data,
            output_suffix=args.output_suffix,
            slurm_params=slurm_params,
            dt=args.dt,
            template_file=args.template,
            qualified_mds=qualified_mds  # 新增参数
        )
        
        if results['success']:
            print("\n✅ SLURM脚本生成成功!")
            print("=" * 50)
            print(f"📊 发现MD任务: {results['total_tasks']} 个")
            print(f"📦 生成批次: {results['num_batches']} 个")
            print(f"📄 生成脚本: {len(results['scripts'])} 个")
            print(f"📁 脚本目录: {results['output_script_dir']}")
            print(f"🚀 提交脚本: {results['submission_script']}")
            print(f"💾 数据输出: {results['base_output_dir']}")
            
            print("\n🎯 下一步操作:")
            print(f"1. 检查生成的脚本: ls {results['output_script_dir']}")
            print(f"2. 提交所有作业: bash {results['submission_script']}")
            print("3. 监控作业状态: squeue -u $USER")
            
            # 显示生成的脚本列表
            if args.verbose:
                print(f"\n📋 生成的SLURM脚本:")
                for i, script in enumerate(results['scripts'], 1):
                    script_name = Path(script).name
                    print(f"   {i:2d}. {script_name}")
            
        else:
            print(f"\n❌ 脚本生成失败: {results['error']}")
            return 1
            
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1
    
    print("\n🎉 完成!")
    return 0

if __name__ == "__main__":
    exit(main())
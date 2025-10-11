#!/usr/bin/env python3
"""
MD轨迹分析运行脚本

基于质量检测和PBC预处理后的规范文件进行RMSD等后续分析。

使用方法:
    python run_trajectory_analysis.py /path/to/processed/simulations
    python run_trajectory_analysis.py /path/to/processed/simulations --analysis rmsd rg distances
    python run_trajectory_analysis.py /path/to/processed/simulations --single /path/to/specific/md
"""

import sys
import argparse
import logging
import time
from pathlib import Path
from typing import List, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from aftermd.analysis.analysis_pipeline import AnalysisPipeline, create_analysis_pipeline


def setup_logging(log_level: str = "INFO"):
    """设置日志配置"""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('trajectory_analysis.log')
        ]
    )


def find_processed_trajectories(root_dir: Path) -> List[Dict[str, str]]:
    """
    发现已处理的轨迹文件
    
    查找符合以下结构的目录:
    {MD_Product}_processed/
    ├── {MD_Product}_processed.xtc
    └── md.gro (或其他参考结构)
    
    Args:
        root_dir: 根目录路径
        
    Returns:
        轨迹文件信息列表
    """
    trajectories = []
    
    # 查找所有 *_processed 目录
    for processed_dir in root_dir.rglob("*_processed"):
        if not processed_dir.is_dir():
            continue
        
        # 查找轨迹文件
        trajectory_files = list(processed_dir.glob("*.xtc")) + list(processed_dir.glob("*.trr"))
        if not trajectory_files:
            continue
        
        # 选择最大的轨迹文件 (通常是主要的处理结果)
        trajectory_file = max(trajectory_files, key=lambda f: f.stat().st_size)
        
        # 查找参考结构文件
        reference_files = (
            list(processed_dir.glob("md.gro")) +
            list(processed_dir.glob("*.gro")) +
            list(processed_dir.glob("*.tpr")) +
            list(processed_dir.glob("*.pdb"))
        )
        
        if not reference_files:
            print(f"⚠️  跳过 {processed_dir.name}: 未找到参考结构文件")
            continue
        
        # 优先选择 md.gro
        reference_file = None
        for ref in reference_files:
            if ref.name == "md.gro":
                reference_file = ref
                break
        if not reference_file:
            reference_file = reference_files[0]
        
        trajectories.append({
            "name": processed_dir.name,
            "trajectory": str(trajectory_file),
            "reference": str(reference_file),
            "output_dir": str(processed_dir / "analysis")
        })
    
    return trajectories


def analyze_single_trajectory(trajectory_info: Dict[str, str], 
                            analysis_types: List[str],
                            **kwargs) -> Dict[str, str]:
    """
    分析单个轨迹
    
    Args:
        trajectory_info: 轨迹信息字典
        analysis_types: 要运行的分析类型列表
        **kwargs: 其他参数
        
    Returns:
        分析结果摘要
    """
    name = trajectory_info["name"]
    
    try:
        # 创建分析流水线
        pipeline = create_analysis_pipeline(
            processed_trajectory=trajectory_info["trajectory"],
            reference_structure=trajectory_info["reference"],
            output_dir=trajectory_info["output_dir"]
        )
        
        # 确定要运行的分析
        enable_rmsd = "rmsd" in analysis_types
        enable_rg = "rg" in analysis_types or "radius_gyration" in analysis_types
        enable_distances = "distances" in analysis_types or "distance" in analysis_types
        enable_rdf = "rdf" in analysis_types
        enable_hbonds = "hbonds" in analysis_types or "hydrogen_bonds" in analysis_types
        
        # 运行综合分析
        results = pipeline.run_comprehensive_analysis(
            enable_rmsd=enable_rmsd,
            enable_rg=enable_rg,
            enable_distances=enable_distances,
            enable_rdf=enable_rdf,
            enable_hbonds=enable_hbonds,
            **kwargs
        )
        
        # 保存结果
        pipeline.save_results()
        
        summary = results.get("analysis_summary", {})
        completed = len(summary.get("completed_analyses", []))
        total_time = summary.get("total_time_seconds", 0)
        
        return {
            "status": "success",
            "name": name,
            "completed_analyses": completed,
            "total_time": total_time,
            "output_dir": trajectory_info["output_dir"]
        }
        
    except Exception as e:
        return {
            "status": "failed",
            "name": name,
            "error": str(e)
        }


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="MD轨迹分析工具 - 基于预处理轨迹进行RMSD等分析",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  分析所有预处理轨迹:
    python run_trajectory_analysis.py /path/to/processed/simulations
    
  只进行RMSD和回转半径分析:
    python run_trajectory_analysis.py /path/to/processed/simulations --analysis rmsd rg
    
  分析单个轨迹:
    python run_trajectory_analysis.py /path/to/processed/simulations --single md_product_1_processed
    
  并行分析 (4个进程):
    python run_trajectory_analysis.py /path/to/processed/simulations --max-workers 4
    
  详细日志:
    python run_trajectory_analysis.py /path/to/processed/simulations --log-level DEBUG
        """
    )
    
    # 必需参数
    parser.add_argument(
        "processed_root",
        help="包含已处理轨迹的根目录"
    )
    
    # 分析选项
    parser.add_argument(
        "--analysis", "-a",
        nargs="+",
        choices=["rmsd", "rg", "radius_gyration", "distances", "distance", "rdf", "hbonds", "hydrogen_bonds"],
        default=["rmsd", "rg", "distances"],
        help="要运行的分析类型 (默认: rmsd rg distances)"
    )
    
    parser.add_argument(
        "--single", "-s",
        help="只分析指定的单个轨迹目录名"
    )
    
    # 并行处理
    parser.add_argument(
        "--max-workers",
        type=int,
        default=1,
        help="最大并行进程数 (默认: 1)"
    )
    
    # 日志选项
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="日志级别 (默认: INFO)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="显示详细信息"
    )
    
    # 其他选项
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="干运行模式，只显示将要分析的轨迹"
    )
    
    args = parser.parse_args()
    
    # 设置日志
    log_level = "DEBUG" if args.verbose else args.log_level
    setup_logging(log_level)
    logger = logging.getLogger(__name__)
    
    # 验证输入目录
    processed_root = Path(args.processed_root)
    if not processed_root.exists():
        print(f"❌ 错误: 目录不存在 - {processed_root}")
        return 1
    
    if not processed_root.is_dir():
        print(f"❌ 错误: 路径不是目录 - {processed_root}")
        return 1
    
    # 发现轨迹文件
    print("🔍 搜索已处理的轨迹文件...")
    trajectories = find_processed_trajectories(processed_root)
    
    if not trajectories:
        print("❌ 未找到任何已处理的轨迹文件")
        print("💡 确保目录结构为: {MD_Product}_processed/{MD_Product}_processed.xtc")
        return 1
    
    print(f"📁 发现 {len(trajectories)} 个已处理的轨迹")
    
    # 过滤单个轨迹
    if args.single:
        trajectories = [t for t in trajectories if args.single in t["name"]]
        if not trajectories:
            print(f"❌ 未找到指定的轨迹: {args.single}")
            return 1
        print(f"🎯 仅分析指定轨迹: {trajectories[0]['name']}")
    
    # 干运行模式
    if args.dry_run:
        print("\n🔍 干运行模式 - 将要分析的轨迹:")
        for i, traj in enumerate(trajectories, 1):
            print(f"  {i:2d}. {traj['name']}")
            print(f"      轨迹: {Path(traj['trajectory']).name}")
            print(f"      参考: {Path(traj['reference']).name}")
            print(f"      输出: {traj['output_dir']}")
        print(f"\n📊 分析类型: {', '.join(args.analysis)}")
        print(f"🔧 并行进程: {args.max_workers}")
        return 0
    
    # 显示分析计划
    print(f"\n🚀 开始轨迹分析")
    print(f"📊 分析类型: {', '.join(args.analysis)}")
    print(f"🔧 并行进程: {args.max_workers}")
    print(f"📁 输出目录: 各轨迹的 analysis/ 子目录")
    print("=" * 60)
    
    # 开始分析
    start_time = time.time()
    results = []
    
    if args.max_workers == 1:
        # 串行处理
        for i, traj_info in enumerate(trajectories, 1):
            print(f"\n[{i:3d}/{len(trajectories)}] 分析 {traj_info['name']}")
            result = analyze_single_trajectory(traj_info, args.analysis)
            results.append(result)
            
            if result["status"] == "success":
                print(f"   ✅ 成功 (耗时: {result['total_time']:.1f}s, 分析: {result['completed_analyses']})")
            else:
                print(f"   ❌ 失败: {result['error']}")
    
    else:
        # 并行处理
        print(f"\n🔄 使用 {args.max_workers} 个进程并行分析...")
        
        with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
            # 提交所有任务
            future_to_traj = {
                executor.submit(analyze_single_trajectory, traj_info, args.analysis): traj_info
                for traj_info in trajectories
            }
            
            # 收集结果
            for i, future in enumerate(as_completed(future_to_traj), 1):
                traj_info = future_to_traj[future]
                try:
                    result = future.result()
                    results.append(result)
                    
                    if result["status"] == "success":
                        print(f"[{i:3d}/{len(trajectories)}] ✅ {result['name']} (耗时: {result['total_time']:.1f}s)")
                    else:
                        print(f"[{i:3d}/{len(trajectories)}] ❌ {result['name']}: {result['error']}")
                        
                except Exception as e:
                    print(f"[{i:3d}/{len(trajectories)}] ❌ {traj_info['name']}: 进程异常 - {e}")
                    results.append({
                        "status": "failed",
                        "name": traj_info["name"],
                        "error": f"进程异常: {e}"
                    })
    
    # 统计结果
    end_time = time.time()
    total_time = end_time - start_time
    
    successful = [r for r in results if r["status"] == "success"]
    failed = [r for r in results if r["status"] == "failed"]
    
    print(f"\n{'='*60}")
    print("📊 轨迹分析完成!")
    print(f"{'='*60}")
    print(f"⏱️  总耗时: {total_time:.1f} 秒")
    print(f"✅ 成功: {len(successful)} 个轨迹")
    print(f"❌ 失败: {len(failed)} 个轨迹")
    
    if successful:
        total_analysis_time = sum(r["total_time"] for r in successful)
        avg_time = total_analysis_time / len(successful)
        print(f"📈 平均分析时间: {avg_time:.1f} 秒/轨迹")
        
        print(f"\n📁 分析结果位置:")
        for result in successful[:5]:  # 显示前5个
            print(f"   • {result['name']}: {result['output_dir']}")
        if len(successful) > 5:
            print(f"   ... 以及其他 {len(successful) - 5} 个")
    
    if failed:
        print(f"\n❌ 失败的轨迹:")
        for result in failed:
            print(f"   • {result['name']}: {result['error']}")
    
    print(f"\n💡 查看各轨迹的分析报告: */analysis/reports/analysis_report.md")
    print(f"💡 查看绘图结果: */analysis/plots/")
    
    return 0 if not failed else 1


if __name__ == "__main__":
    exit(main())
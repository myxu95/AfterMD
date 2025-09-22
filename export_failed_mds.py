#!/usr/bin/env python3
"""
异常MD Product导出工具

从质量检查结果中导出异常的MD，复制整理到指定目录，方便修正后重新处理。

使用方法:
    python export_failed_mds.py quality_check_results/failed_mds.txt
    python export_failed_mds.py quality_check_results/failed_mds.txt --output ./failed_mds_backup
"""

import sys
import argparse
import shutil
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple


class FailedMDExporter:
    """异常MD导出器"""
    
    def __init__(self, output_dir: Path):
        """
        初始化导出器
        
        Args:
            output_dir: 导出目录
        """
        self.output_dir = output_dir
        self.export_log = []
        
    def load_failed_mds(self, failed_list_file: Path) -> List[Dict]:
        """从失败列表文件加载异常MD信息"""
        failed_mds = []
        
        try:
            with open(failed_list_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    # 解析格式: MD路径\t失败原因
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        md_path = parts[0].strip()
                        reasons = parts[1].strip()
                    else:
                        md_path = parts[0].strip()
                        reasons = "未知原因"
                    
                    failed_mds.append({
                        'path': md_path,
                        'reasons': reasons,
                        'line_num': line_num
                    })
            
            print(f"📋 从失败列表加载了 {len(failed_mds)} 个异常MD")
            return failed_mds
            
        except Exception as e:
            print(f"❌ 读取失败列表文件出错: {e}")
            return []
    
    def categorize_failures(self, failed_mds: List[Dict]) -> Dict[str, List[Dict]]:
        """按失败原因分类异常MD"""
        categories = {}
        
        for md in failed_mds:
            reasons = md['reasons']
            
            # 识别主要失败类型
            if '链数异常' in reasons or '链数' in reasons:
                category = 'chain_count_error'
                category_name = '链数异常'
            elif 'MD未完成' in reasons or '未完成' in reasons:
                category = 'incomplete_md'
                category_name = 'MD未完成'
            elif '缺失文件' in reasons or '文件不存在' in reasons:
                category = 'missing_files'
                category_name = '文件缺失'
            elif '文件太小' in reasons or '大小' in reasons:
                category = 'file_size_error'
                category_name = '文件大小异常'
            elif '时间不足' in reasons or '模拟时间' in reasons:
                category = 'simulation_time_error'
                category_name = '模拟时间不足'
            else:
                category = 'other_errors'
                category_name = '其他错误'
            
            # 添加分类信息
            md_info = md.copy()
            md_info['category'] = category
            md_info['category_name'] = category_name
            
            if category not in categories:
                categories[category] = []
            categories[category].append(md_info)
        
        return categories
    
    def export_md_directory(self, md_path: str, target_dir: Path, dry_run: bool = False) -> bool:
        """导出单个MD目录"""
        source_path = Path(md_path)
        
        if not source_path.exists():
            self.export_log.append({
                'md_path': md_path,
                'status': 'error',
                'message': '源目录不存在'
            })
            return False
        
        if not source_path.is_dir():
            self.export_log.append({
                'md_path': md_path,
                'status': 'error', 
                'message': '源路径不是目录'
            })
            return False
        
        # 目标路径
        md_name = source_path.name
        target_path = target_dir / md_name
        
        if dry_run:
            self.export_log.append({
                'md_path': md_path,
                'target_path': str(target_path),
                'status': 'dry_run',
                'message': '模拟导出'
            })
            return True
        
        try:
            # 复制整个MD目录
            if target_path.exists():
                print(f"   ⚠️  目标已存在，跳过: {md_name}")
                self.export_log.append({
                    'md_path': md_path,
                    'target_path': str(target_path),
                    'status': 'skipped',
                    'message': '目标已存在'
                })
                return False
            
            shutil.copytree(source_path, target_path)
            
            # 计算目录大小
            total_size = sum(f.stat().st_size for f in target_path.rglob('*') if f.is_file())
            size_mb = total_size / (1024 * 1024)
            
            self.export_log.append({
                'md_path': md_path,
                'target_path': str(target_path),
                'status': 'success',
                'size_mb': size_mb,
                'message': f'导出成功，大小: {size_mb:.1f} MB'
            })
            
            print(f"   ✅ 导出: {md_name} ({size_mb:.1f} MB)")
            return True
            
        except Exception as e:
            self.export_log.append({
                'md_path': md_path,
                'target_path': str(target_path),
                'status': 'error',
                'message': f'导出失败: {str(e)}'
            })
            print(f"   ❌ 导出失败: {md_name} - {e}")
            return False
    
    def export_by_category(self, categories: Dict[str, List[Dict]], dry_run: bool = False) -> Dict:
        """按分类导出异常MD"""
        export_stats = {
            'total_attempted': 0,
            'successful': 0,
            'failed': 0,
            'skipped': 0,
            'categories': {}
        }
        
        for category, mds in categories.items():
            category_name = mds[0]['category_name'] if mds else category
            category_dir = self.output_dir / f"{category}_{category_name}"
            
            print(f"\n📁 导出分类: {category_name} ({len(mds)} 个MD)")
            print(f"   目标目录: {category_dir}")
            
            if not dry_run:
                category_dir.mkdir(parents=True, exist_ok=True)
            
            category_stats = {
                'attempted': len(mds),
                'successful': 0,
                'failed': 0,
                'skipped': 0
            }
            
            for md in mds:
                export_stats['total_attempted'] += 1
                category_stats['attempted'] += 1
                
                if self.export_md_directory(md['path'], category_dir, dry_run):
                    export_stats['successful'] += 1
                    category_stats['successful'] += 1
                else:
                    # 检查是跳过还是失败
                    last_log = self.export_log[-1] if self.export_log else {}
                    if last_log.get('status') == 'skipped':
                        export_stats['skipped'] += 1
                        category_stats['skipped'] += 1
                    else:
                        export_stats['failed'] += 1
                        category_stats['failed'] += 1
            
            export_stats['categories'][category] = category_stats
            
            # 创建分类说明文件
            if not dry_run and mds:
                readme_file = category_dir / "README.txt"
                with open(readme_file, 'w', encoding='utf-8') as f:
                    f.write(f"异常MD分类: {category_name}\n")
                    f.write(f"导出时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"MD数量: {len(mds)}\n\n")
                    f.write("失败原因:\n")
                    for md in mds:
                        f.write(f"- {Path(md['path']).name}: {md['reasons']}\n")
        
        return export_stats
    
    def save_export_report(self) -> Path:
        """保存导出报告"""
        report_file = self.output_dir / "export_report.json"
        
        report = {
            'export_time': datetime.now().isoformat(),
            'output_directory': str(self.output_dir),
            'export_log': self.export_log
        }
        
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        return report_file
    
    def print_summary(self, export_stats: Dict):
        """打印导出摘要"""
        print(f"\n{'='*60}")
        print("📊 异常MD导出摘要")
        print(f"{'='*60}")
        print(f"总计尝试: {export_stats['total_attempted']} 个MD")
        print(f"成功导出: {export_stats['successful']} 个")
        print(f"导出失败: {export_stats['failed']} 个")
        print(f"跳过重复: {export_stats['skipped']} 个")
        
        print(f"\n📁 分类统计:")
        for category, stats in export_stats['categories'].items():
            print(f"   {category}: {stats['successful']}/{stats['attempted']} 成功")
        
        total_size = sum(log.get('size_mb', 0) for log in self.export_log 
                        if log.get('status') == 'success')
        print(f"\n💾 总导出大小: {total_size:.1f} MB")
        print(f"📁 导出目录: {self.output_dir}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="异常MD Product导出工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  基本导出:
    python export_failed_mds.py quality_check_results/failed_mds.txt
    
  指定输出目录:
    python export_failed_mds.py quality_check_results/failed_mds.txt \\
        --output ./failed_mds_backup
    
  预览模式:
    python export_failed_mds.py quality_check_results/failed_mds.txt --dry-run
    
输出结构:
  failed_mds_backup/
  ├── chain_count_error_链数异常/     # 链数异常的MD
  ├── incomplete_md_MD未完成/         # 未完成的MD
  ├── missing_files_文件缺失/         # 文件缺失的MD
  └── export_report.json             # 导出报告
        """
    )
    
    # 必需参数
    parser.add_argument(
        "failed_list_file",
        help="质量检查产生的失败MD列表文件 (failed_mds.txt)"
    )
    
    # 输出参数
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default="./failed_mds_backup",
        help="导出目录 (默认: ./failed_mds_backup)"
    )
    
    # 操作选项
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="预览模式，只显示将要执行的操作，不实际复制文件"
    )
    
    parser.add_argument(
        "--force",
        action="store_true", 
        help="强制覆盖已存在的目标目录"
    )
    
    args = parser.parse_args()
    
    # 验证输入文件
    failed_list_file = Path(args.failed_list_file)
    if not failed_list_file.exists():
        print(f"❌ 错误: 失败列表文件不存在 - {failed_list_file}")
        return 1
    
    # 创建输出目录
    output_dir = args.output
    if output_dir.exists() and not args.force and not args.dry_run:
        print(f"❌ 错误: 输出目录已存在 - {output_dir}")
        print("使用 --force 强制覆盖，或选择其他目录")
        return 1
    
    print("🚀 异常MD Product导出工具")
    print("=" * 60)
    print(f"📂 失败列表: {failed_list_file}")
    print(f"📁 导出目录: {output_dir}")
    if args.dry_run:
        print("🔍 模式: 预览模式 (不实际复制文件)")
    print("=" * 60)
    
    try:
        # 创建导出器
        exporter = FailedMDExporter(output_dir)
        
        # 加载失败MD列表
        failed_mds = exporter.load_failed_mds(failed_list_file)
        if not failed_mds:
            print("❌ 没有找到异常MD需要导出")
            return 1
        
        # 按失败原因分类
        categories = exporter.categorize_failures(failed_mds)
        
        print(f"\n🔍 异常MD分类:")
        for category, mds in categories.items():
            category_name = mds[0]['category_name'] if mds else category
            print(f"   {category_name}: {len(mds)} 个")
        
        # 确认操作
        if not args.dry_run:
            print(f"\n将导出 {len(failed_mds)} 个异常MD到 {output_dir}")
            user_input = input("确认继续? (y/N): ").strip().lower()
            if user_input not in ['y', 'yes']:
                print("操作已取消")
                return 0
        
        # 创建输出目录
        if not args.dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)
        
        # 按分类导出
        export_stats = exporter.export_by_category(categories, args.dry_run)
        
        # 保存导出报告
        if not args.dry_run:
            report_file = exporter.save_export_report()
            print(f"\n📋 导出报告已保存: {report_file}")
        
        # 显示摘要
        exporter.print_summary(export_stats)
        
        if args.dry_run:
            print(f"\n💡 这是预览模式。要实际导出，请运行:")
            print(f"   python {sys.argv[0]} {args.failed_list_file} --output {args.output}")
        else:
            print(f"\n🎉 异常MD导出完成!")
            print(f"   可以修正这些MD后重新运行质量检查")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ 导出过程出错: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
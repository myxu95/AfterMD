#!/usr/bin/env python3
"""
测试简化的GRO链检测方案
"""

import logging
from pathlib import Path
from aftermd.utils.simple_gro_detector import create_shortest_chain_index

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def test_simple_gro_detection():
    """测试简化的GRO链检测"""

    md_gro_path = "/Users/xumingyu/Documents/work/research/AfrerMD/md.gro"
    topology_path = "/Users/xumingyu/Documents/work/research/AfrerMD/md.gro"
    output_dir = "/tmp/claude/test_output"

    print("=" * 50)
    print("测试简化GRO链检测")
    print("=" * 50)

    if not Path(md_gro_path).exists():
        print(f"错误: md.gro文件不存在: {md_gro_path}")
        return

    print(f"md.gro文件: {md_gro_path}")
    print(f"输出目录: {output_dir}")
    print()

    # 创建最短链index
    index_file = create_shortest_chain_index(
        md_gro_path=md_gro_path,
        topology_path=topology_path,
        output_dir=output_dir
    )

    if index_file:
        print(f"✓ 成功生成index文件: {index_file}")

        # 检查文件
        if Path(index_file).exists():
            file_size = Path(index_file).stat().st_size
            print(f"✓ 文件大小: {file_size} bytes")

            # 显示文件末尾（最短链部分）
            try:
                with open(index_file, 'r') as f:
                    lines = f.readlines()

                print(f"✓ 文件总行数: {len(lines)}")
                print("\n最后10行内容:")
                for line in lines[-10:]:
                    print(f"  {line.rstrip()}")

                # 查找Shortest_Chain组
                shortest_found = False
                for i, line in enumerate(lines):
                    if "Shortest_Chain" in line:
                        print(f"\n✓ 找到Shortest_Chain组 (第{i+1}行)")
                        # 显示该组的几行内容
                        for j in range(min(5, len(lines) - i)):
                            if i + j < len(lines):
                                print(f"  {lines[i+j].rstrip()}")
                        shortest_found = True
                        break

                if not shortest_found:
                    print("✗ 未找到Shortest_Chain组")

            except Exception as e:
                print(f"✗ 读取文件失败: {e}")
        else:
            print(f"✗ 文件不存在: {index_file}")
    else:
        print("✗ 生成index文件失败")

    print("\n" + "=" * 50)
    print("测试完成")

if __name__ == "__main__":
    test_simple_gro_detection()
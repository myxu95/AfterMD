import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class SimpleGroDetector:
    """简化的GRO文件链检测器，专注于从md.gro生成最短链index"""

    def __init__(self, gro_file: str, topology_file: str, gmx_executable: str = "gmx"):
        self.gro_file = gro_file
        self.topology_file = topology_file
        self.gmx = gmx_executable

        # 蛋白质残基类型
        self.protein_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'ACE', 'NME', 'HIE', 'HID', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN'
        }

    def generate_shortest_chain_index(self, output_dir: str) -> Optional[str]:
        """
        直接从md.gro生成包含最短链的index文件

        Args:
            output_dir: 输出目录

        Returns:
            生成的index文件路径
        """
        try:
            # 1. 检测链
            chains = self._detect_chains_from_gro()
            if not chains:
                logger.warning("未检测到有效链")
                return None

            # 2. 找最短链
            shortest_chain_atoms = self._find_shortest_chain_atoms(chains)
            if not shortest_chain_atoms:
                logger.warning("未找到最短链")
                return None

            # 3. 生成index文件
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            index_file = output_path / "shortest_chain.ndx"
            self._create_index_file_with_shortest_chain(str(index_file), shortest_chain_atoms)

            logger.info(f"生成最短链index文件: {index_file}")
            logger.info(f"最短链包含 {len(shortest_chain_atoms)} 个原子")

            return str(index_file)

        except Exception as e:
            logger.error(f"生成index文件失败: {e}")
            return None

    def _detect_chains_from_gro(self) -> Dict[int, List[int]]:
        """
        从GRO文件检测链，返回每条链的原子ID列表

        Returns:
            {chain_id: [atom_ids]}
        """
        try:
            with open(self.gro_file, 'r') as f:
                lines = f.readlines()

            # 跳过标题和原子数行，最后一行是盒子
            atom_lines = lines[2:-1]

            # 解析所有蛋白质原子
            protein_atoms = []
            for line in atom_lines:
                if len(line.strip()) < 44:
                    continue

                try:
                    # GRO格式: resnum+resname(5位) + atomname(5位) + atomnum(5位) + xyz
                    resid_resname = line[0:5].strip()  # 前5个字符包含resnum+resname
                    atom_name = line[5:10].strip()     # 接下来5个字符是atom name
                    atom_id = int(line[10:15].strip()) # 再5个字符是atom number

                    # 分离residue ID和name - 更灵活的匹配
                    match = re.match(r'(\d+)([A-Z]+)', resid_resname)
                    if not match:
                        # 尝试另一种格式：可能数字和字母之间有空格或格式不同
                        match = re.search(r'(\d+)\s*([A-Z]+)', resid_resname)
                        if not match:
                            continue

                    resid = int(match.group(1))
                    resname = match.group(2)

                    if resname in self.protein_residues:
                        protein_atoms.append({
                            'resid': resid,
                            'atom_id': atom_id
                        })

                except (ValueError, IndexError):
                    continue

            # 按residue ID排序
            protein_atoms.sort(key=lambda x: x['resid'])

            # 根据residue ID连续性分链
            chains = {}
            current_chain = 0
            current_atoms = []
            last_resid = None

            for atom in protein_atoms:
                if last_resid is not None and atom['resid'] - last_resid > 1:
                    # 新链开始
                    if len(current_atoms) >= 20:  # 至少20个原子才算有效链
                        chains[current_chain] = current_atoms
                        current_chain += 1
                    current_atoms = []

                current_atoms.append(atom['atom_id'])
                last_resid = atom['resid']

            # 添加最后一条链
            if len(current_atoms) >= 20:
                chains[current_chain] = current_atoms

            logger.info(f"检测到 {len(chains)} 条蛋白质链")
            for chain_id, atoms in chains.items():
                logger.info(f"  链 {chain_id}: {len(atoms)} 个原子")

            return chains

        except Exception as e:
            logger.error(f"从GRO文件检测链失败: {e}")
            return {}

    def _find_shortest_chain_atoms(self, chains: Dict[int, List[int]]) -> Optional[List[int]]:
        """找到最短链的原子ID列表"""
        if not chains:
            return None

        # 按原子数排序，选择最短的
        shortest_chain_id = min(chains.keys(), key=lambda x: len(chains[x]))
        shortest_atoms = chains[shortest_chain_id]

        logger.info(f"最短链: 链 {shortest_chain_id} ({len(shortest_atoms)} 原子)")
        return shortest_atoms

    def _create_index_file_with_shortest_chain(self, output_file: str, shortest_atoms: List[int]):
        """创建包含默认组和最短链的index文件"""

        # 1. 先尝试从topology文件生成基础index文件
        temp_base = output_file + ".temp"
        base_content = ""

        try:
            cmd = [self.gmx, "make_ndx", "-f", self.topology_file, "-o", temp_base]
            result = subprocess.run(
                cmd,
                input="q\n",
                text=True,
                capture_output=True,
                check=True,
                timeout=30
            )

            # 读取生成的基础内容
            if Path(temp_base).exists():
                with open(temp_base, 'r') as f:
                    base_content = f.read()
                Path(temp_base).unlink(missing_ok=True)

        except subprocess.CalledProcessError as e:
            logger.warning(f"从topology生成基础index失败: {e}")
            logger.info("尝试从GRO文件生成基础index")

            # 备选方案：从GRO文件生成基础index
            try:
                cmd = [self.gmx, "make_ndx", "-f", self.gro_file, "-o", temp_base]
                result = subprocess.run(
                    cmd,
                    input="q\n",
                    text=True,
                    capture_output=True,
                    check=True,
                    timeout=30
                )

                if Path(temp_base).exists():
                    with open(temp_base, 'r') as f:
                        base_content = f.read()
                    Path(temp_base).unlink(missing_ok=True)

            except subprocess.CalledProcessError as e2:
                logger.warning(f"从GRO生成基础index也失败: {e2}")
                # 最后备选：使用默认的index内容
                base_content = self._create_minimal_index_content()

        # 2. 创建最终文件
        try:
            # 确保输出目录存在
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # 格式化最短链组
            shortest_content = self._format_shortest_chain_group(shortest_atoms)

            # 写入最终文件
            with open(output_file, 'w') as f:
                if base_content:
                    f.write(base_content)
                    if not base_content.endswith('\n'):
                        f.write("\n")
                f.write(shortest_content)

            logger.info(f"成功创建包含最短链的index文件: {output_file}")

        except Exception as e:
            logger.error(f"创建最终index文件失败: {e}")
            # 清理可能的临时文件
            Path(temp_base).unlink(missing_ok=True)
            raise

    def _create_minimal_index_content(self) -> str:
        """创建最小的index文件内容（备选方案）"""
        return """[ System ]
   1

[ Protein ]
   1

[ Protein-H ]
   1

[ C-alpha ]
   1

[ Backbone ]
   1

[ MainChain ]
   1

[ MainChain+Cb ]
   1

[ MainChain+H ]
   1

[ SideChain ]
   1

[ SideChain-H ]
   1

"""

    def _format_shortest_chain_group(self, atom_ids: List[int]) -> str:
        """格式化最短链组"""
        lines = ["[ Shortest_Chain ]"]

        # 每行12个原子，右对齐5位
        for i in range(0, len(atom_ids), 12):
            line_atoms = atom_ids[i:i + 12]
            formatted_line = "".join(f"{atom_id:>5}" for atom_id in line_atoms)
            lines.append(formatted_line)

        return "\n".join(lines) + "\n"


def create_shortest_chain_index(md_gro_path: str, topology_path: str,
                               output_dir: str, gmx_executable: str = "gmx") -> Optional[str]:
    """
    便捷函数：从md.gro创建最短链index文件

    Args:
        md_gro_path: md.gro文件路径
        topology_path: 拓扑文件路径
        output_dir: 输出目录
        gmx_executable: gmx命令

    Returns:
        生成的index文件路径
    """
    if not Path(md_gro_path).exists():
        logger.error(f"md.gro文件不存在: {md_gro_path}")
        return None

    detector = SimpleGroDetector(md_gro_path, topology_path, gmx_executable)
    return detector.generate_shortest_chain_index(output_dir)
"""
AfterMD Core Module - 核心处理功能

包含PBC处理、批量处理等核心功能
"""

from .pbc_processor import PBCProcessor
from .batch_process import process_md_tasks, discover_md_tasks

# 简化的高层API名称
process_pbc = process_md_tasks
discover_tasks = discover_md_tasks

__all__ = [
    'PBCProcessor',
    'process_md_tasks',
    'discover_md_tasks',
    'process_pbc',
    'discover_tasks',
]
